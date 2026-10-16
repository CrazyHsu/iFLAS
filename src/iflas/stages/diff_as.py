"""Differential AS (legacy ``diff_as.py``) backed by rMATS-turbo.

``merge_isoforms`` builds the cross-sample isoform GTF rMATS annotates with, ``run_rmats`` runs every pairwise
condition comparison (BAMs from a 7-column table or from the samples' own NGS alignments), ``write_sig_diff_as``
keeps FDR <= 0.05 events per AS type, and ``diff_as_project`` is the project-level driver (optional GO enrichment).
"""

from __future__ import annotations

import copy
import csv
import os
import shutil
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import pandas as pd

from ..io.records import BedFile
from ..pu.cpc2 import _find_utils_dir
from ..runtime import ToolRunner, getCurrentTime, quote as q

RMATS_FDR = 0.05
MATS_TABLES = {"IR": "RI.MATS.JC.txt", "SE": "SE.MATS.JC.txt", "A3SS": "A3SS.MATS.JC.txt", "A5SS": "A5SS.MATS.JC.txt"}

_COMP_COLUMNS = {"project", "sampleName", "condition", "repeat", "bamFile", "paired", "readsLength"}


@dataclass(frozen=True)
class BamRecord:
    project: str
    sample_name: str
    condition: str
    repeat: str
    bam_file: str
    paired: str
    reads_length: int

    @property
    def comparison_name(self) -> str:
        return f"{self.project}_{self.condition}"


def _split_csv(value: str | Iterable[str] | None) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        items = value.split(",")
    else:
        items = [str(item) for item in value]
    return [item.strip() for item in items if item and item.strip()]


def _resolve_file(path: str | Path, *, base_dir: Path, label: str) -> str:
    resolved = Path(path).expanduser()
    if not resolved.is_absolute():
        resolved = base_dir / resolved
    resolved = resolved.resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    return str(resolved)


def load_bam_comparison_table(path: str | Path) -> list[BamRecord]:
    """Load the legacy 7-column comparison table used by ``diff_as``."""
    table = Path(path).expanduser().resolve()
    if not table.is_file():
        raise FileNotFoundError(f"comparison table does not exist: {path}")
    rows: list[BamRecord] = []
    with table.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = sorted(_COMP_COLUMNS - set(reader.fieldnames or []))
        if missing:
            raise ValueError("comparison table missing column(s): " + ", ".join(missing))
        for line_no, raw in enumerate(reader, start=2):
            try:
                reads_length = int(str(raw["readsLength"]).strip())
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid readsLength at {table}:{line_no}") from exc
            paired = str(raw["paired"]).strip()
            if paired not in {"paired", "single"}:
                raise ValueError(f"invalid paired value at {table}:{line_no}: {paired!r}")
            rows.append(
                BamRecord(
                    project=str(raw["project"]).strip(),
                    sample_name=str(raw["sampleName"]).strip(),
                    condition=str(raw["condition"]).strip(),
                    repeat=str(raw["repeat"]).strip(),
                    bam_file=_resolve_file(raw["bamFile"], base_dir=table.parent, label="BAM file"),
                    paired=paired,
                    reads_length=reads_length,
                )
            )
    if not rows:
        raise ValueError(f"comparison table is empty: {path}")
    return rows


def _conditions_in_order(records: list[BamRecord], requested: str | Iterable[str] | None) -> list[str]:
    available = []
    for record in records:
        if record.condition not in available:
            available.append(record.condition)
    requested_items = _split_csv(requested)
    if not requested_items:
        return available
    missing = [item for item in requested_items if item not in available]
    if missing:
        raise ValueError("unknown condition(s): " + ", ".join(missing))
    return requested_items


def _records_by_condition(records: list[BamRecord]) -> dict[str, list[BamRecord]]:
    grouped: dict[str, list[BamRecord]] = {}
    for record in records:
        grouped.setdefault(record.condition, []).append(record)
    return grouped


def _homogeneous(records: list[BamRecord], attr: str, condition: str):
    values = {getattr(record, attr) for record in records}
    if len(values) != 1:
        raise ValueError(f"condition {condition!r} has mixed {attr}: {sorted(values)}")
    return next(iter(values))


def bam_records_from_samples(config, conditions=None) -> list[BamRecord]:
    """BAM records for the samples' own NGS alignments (legacy ``getSample2Bam`` without a table):
    ``<base>/mapping/rna-seq/alignment/repeat<i>.sorted.bam`` per repeat; ``conditions`` restricts the samples.
    Comparison names are ``<project>_<condition>`` (legacy used ``<sample>_<repeat>`` and picked one at random)."""
    from ..orchestrator import parse_ngs_repeats, sample_base_dir

    records = []
    for sample in config.samples:
        if not sample.condition or (conditions and sample.condition not in conditions):
            continue
        align_dir = Path(sample_base_dir(config, sample)) / "mapping" / "rna-seq" / "alignment"
        for i, _ in enumerate(parse_ngs_repeats(sample)):
            records.append(BamRecord(sample.project_name, sample.sample_name, sample.condition, f"repeat{i}",
                                     str(align_dir / f"repeat{i}.sorted.bam"), sample.ngs_reads_paired or "paired",
                                     int(sample.ngs_reads_length or 0)))
    return records


def as_isoform_names(ordinary_as_dir) -> set[str]:
    """Isoforms taking part in confident AS events (legacy: cut -f 8,10 of A3SS/A5SS/IR.confident.bed6+ and
    -f 16,18 of SE.confident.bed12+, comma lists)."""
    lr = Path(ordinary_as_dir) / "LR"
    names = set()
    for fn, cols in (("A3SS.confident.bed6+", (7, 9)), ("A5SS.confident.bed6+", (7, 9)), ("IR.confident.bed6+", (7, 9)),
                     ("SE.confident.bed12+", (15, 17))):
        path = lr / fn
        if not path.is_file():
            raise FileNotFoundError(str(path))
        for line in open(path):
            c = line.rstrip("\n").split("\t")
            for col in cols:
                if len(c) > col:
                    names.update(x for x in c[col].split(",") if x)
    return names


def merge_isoforms(sample_bases, out_dir, *, pu_filter=False, runner=None) -> str:
    """Legacy ``mergeIsoforms``: per sample keep the AS isoforms (or, with ``pu_filter``, the PU-filtered HQ
    isoforms) of ``refine/isoformGrouped.bed12+``, prefix their names with the sample label, merge identical
    junction chains across samples (name = ``+``-joined, longest span representative) ->
    ``all_sample_merged_iso.bed`` -> ``bed2gpe.pl -g 13`` -> ``genePredToGtf`` -> returns the GTF path."""
    runner = runner or ToolRunner()
    out_dir = Path(out_dir)
    merged = {}
    for label, base in sample_bases:
        base = Path(base)
        grouped = base / "refine" / "isoformGrouped.bed12+"
        if pu_filter:
            hq = base / "hqIsoforms" / "hq.collapsed.bed12+"
            if not hq.is_file():
                raise FileNotFoundError(f"{hq} (run pu_filter first or drop --pu-filter)")
            keep = {l.split("\t")[3] for l in open(hq) if l.strip()}
        else:
            keep = as_isoform_names(base / "as_events" / "ordinary_as")
        for iso in BedFile(str(grouped), type="bed12+").reads.values():
            if iso.name not in keep:
                continue
            iso.name = f"{label}_{iso.name}"
            key = f"{iso.chrom}_{iso.juncChain}" if iso.juncChain else f"{iso.chrom}:{iso.chromStart}-{iso.chromEnd}"   # mono-exon: by span
            merged.setdefault(key, []).append(iso)
    merged_bed = out_dir / "all_sample_merged_iso.bed"
    with open(merged_bed, "w") as out:
        for isos in merged.values():
            rep = copy.copy(sorted(isos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)[0])
            rep.name = "+".join(x.name for x in isos)
            print(str(rep), file=out)
    runner.run("perl {u}/bed2gpe.pl -g 13 all_sample_merged_iso.bed > all_sample_merged_iso.gpe".format(u=q(_find_utils_dir())), cwd=str(out_dir))
    runner.run(["genePredToGtf", "file", "all_sample_merged_iso.gpe", "all_sample_merged_iso.gtf"], cwd=str(out_dir))
    return str(out_dir / "all_sample_merged_iso.gtf")


def write_sig_diff_as(comp_out_dir, fdr=RMATS_FDR) -> dict[str, str]:
    """``<comp>.sigDiffAS/<type>.sig.txt``: rMATS JC tables filtered to FDR <= ``fdr`` (legacy l.207-225)."""
    comp_out_dir = Path(comp_out_dir)
    sig_dir = comp_out_dir.parent / f"{comp_out_dir.name}.sigDiffAS"
    sig_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    for as_type, table in MATS_TABLES.items():
        df = pd.read_csv(comp_out_dir / table, sep="\t")
        sig = sig_dir / f"{as_type}.sig.txt"
        df.loc[df.FDR <= fdr].to_csv(sig, sep="\t", header=True, index=False)
        out[as_type] = str(sig)
    return out


def run_rmats(records, gtf_file, das_dir, *, conditions=None, threads=1, runner=None, sig=True) -> list[dict]:
    """rMATS-turbo for every pairwise comparison of ``conditions`` (default: table order) with the BAMs of
    ``records``; writes ``<comp>.rmats.log`` and, with ``sig``, the ``*.sigDiffAS`` tables."""
    selected = _conditions_in_order(records, conditions)
    if len(selected) < 2:
        raise ValueError("diff_as requires at least two conditions")
    gtf = _resolve_file(gtf_file, base_dir=Path.cwd(), label="GTF file")
    das_dir = Path(das_dir)
    das_dir.mkdir(parents=True, exist_ok=True)
    grouped = _records_by_condition(records)
    runner = runner or ToolRunner({"rmats.py": "iflas-diffsplice"})
    threads = max(1, int(threads))
    outputs = []
    for cond1, cond2 in combinations(selected, 2):
        r1, r2 = grouped[cond1], grouped[cond2]
        paired1, paired2 = _homogeneous(r1, "paired", cond1), _homogeneous(r2, "paired", cond2)
        length1, length2 = _homogeneous(r1, "reads_length", cond1), _homogeneous(r2, "reads_length", cond2)
        if paired1 != paired2 or length1 != length2:
            raise ValueError(f"rMATS requires matching paired/single mode and read length for {cond1!r} vs {cond2!r}")
        if int(length1) <= 0:
            raise ValueError(f"rMATS needs the NGS read length (ngs_reads_length / readsLength) for {cond1!r} vs {cond2!r}")
        comp1_name, comp2_name = r1[0].comparison_name, r2[0].comparison_name
        b1_file, b2_file = das_dir / f"{comp1_name}.txt", das_dir / f"{comp2_name}.txt"
        b1_file.write_text(",".join(r.bam_file for r in r1) + "\n")
        b2_file.write_text(",".join(r.bam_file for r in r2) + "\n")
        comp_name = f"{comp1_name}_vs_{comp2_name}"
        # rMATS >= 4.1 requires --tmp; --variable-read-length keeps trimmed reads (rMATS >= 4.0.2 drops
        # reads whose length differs from --readLength otherwise)
        cmd = ["rmats.py", "--b1", b1_file.name, "--b2", b2_file.name, "--gtf", gtf, "--od", comp_name,
               "--tmp", f"{comp_name}.tmp", "-t", paired1, "--readLength", str(length1), "--variable-read-length",
               "--tstat", str(threads), "--nthread", str(threads)]
        rmats_log = das_dir / f"{comp_name}.rmats.log"
        # rMATS' post step reads every *.rmats file in --tmp: a previous run's files for the same BAMs make it fail (BUG-174)
        shutil.rmtree(das_dir / f"{comp_name}.tmp", ignore_errors=True)
        with rmats_log.open("w", encoding="utf-8") as log_handle:
            runner.run(cmd, cwd=str(das_dir), stdout=log_handle, stderr=log_handle)
        comp = {"name": comp_name, "b1_file": str(b1_file), "b2_file": str(b2_file), "out_dir": str(das_dir / comp_name),
                "log_file": str(rmats_log), "command": cmd}
        if sig:
            comp["sig"] = write_sig_diff_as(das_dir / comp_name)
        outputs.append(comp)
    return outputs


def run_rmats_from_bam_table(comp_cond_file, gtf_file, *, out_dir=".", conditions=None, threads=1, runner=None,
                             sig=False) -> dict:
    """rMATS for all pairwise condition comparisons listed in a 7-column BAM table."""
    records = load_bam_comparison_table(Path(comp_cond_file).expanduser().resolve())
    das_dir = Path(out_dir).expanduser().resolve() / "das"
    outputs = run_rmats(records, gtf_file, das_dir, conditions=conditions, threads=threads, runner=runner, sig=sig)
    return {"das_dir": str(das_dir), "comparisons": outputs}


def diff_as_project(config, comp_cond=None, *, pu_filter=False, go=False, gene2go=None, threads=1, runner=None) -> dict:
    """Legacy ``diff_as`` driver over a loaded config: isoforms of every finished sample dir (raw samples and
    ``<strain>_all`` merges) are merged into one GTF; BAMs come from a 7-column table (``comp_cond`` file) or from
    the samples' own NGS alignments (``comp_cond`` = comma-separated conditions, default all); rMATS per condition
    pair -> ``das/<c1>_vs_<c2>[.sigDiffAS]``; ``go`` runs GO enrichment on the significant genes (``gene2go`` or
    ``options.gene2go``)."""
    from ..orchestrator import _merged_samples, sample_base_dir

    runner = runner or ToolRunner({"rmats.py": "iflas-diffsplice"})
    das_dir = Path(config.dirs.out_dir).expanduser().resolve() / "das"
    das_dir.mkdir(parents=True, exist_ok=True)
    print(getCurrentTime() + " Identify differential alternative spliced genes...")
    if comp_cond and os.path.isfile(str(comp_cond)):
        records = load_bam_comparison_table(Path(comp_cond).resolve())
        conditions = None
    else:
        conditions = [c for c in str(comp_cond).split(",") if c] if comp_cond else None
        records = bam_records_from_samples(config, conditions)
    if not records:
        raise ValueError("diff_as: no NGS BAMs found for the requested conditions (run mapping first)")

    sources = []
    for sample in list(config.samples) + _merged_samples(config):
        base = Path(sample_base_dir(config, sample))
        if conditions and getattr(sample, "condition", None) in conditions or not conditions or sample.sample_name.endswith("_all"):
            if (base / "refine" / "isoformGrouped.bed12+").is_file() and (base / "as_events" / "ordinary_as" / "LR").is_dir():
                sources.append((sample.sample_name, base))
    if not sources:
        raise ValueError("diff_as: no sample dir with refine/isoformGrouped.bed12+ and as_events/ordinary_as found")
    gtf = merge_isoforms(sources, das_dir, pu_filter=pu_filter, runner=runner)
    comparisons = run_rmats(records, gtf, das_dir, conditions=conditions, threads=threads, runner=runner, sig=True)

    if go:
        from .go import run_go_enrichment
        gene2go = gene2go or config.options.gene2go
        if not gene2go:
            print("You don't provide gene2go file, the GO enrichment will not be carried out!")
        else:
            for comp in comparisons:
                sig_dir = Path(comp["sig"]["SE"]).parent
                genes = set()
                for path in comp["sig"].values():
                    df = pd.read_csv(path, sep="\t")
                    genes.update(str(g) for g in df.get("GeneID", pd.Series(dtype=str)).dropna())
                dasg = sig_dir / "dasg.lst"
                dasg.write_text("".join(g + "\n" for g in sorted(genes)))
                comp["go"] = run_go_enrichment([str(dasg)], [comp["name"]], gene2go, out_name="sigDiff", out_dir=str(sig_dir), runner=runner)
    print(getCurrentTime() + " Identify differential alternative spliced genes done!")
    return {"das_dir": str(das_dir), "gtf": gtf, "isoform_sources": [n for n, _ in sources], "comparisons": comparisons}
