"""visual_as stage - gene-structure plots with pyGenomeTracks.

Legacy ``visual_as.py`` drew SpliceGrapher plots (Python-2 only); the migrated stage keeps its products and
layout with pyGenomeTracks: per gene (``<base>/visual_as/<gene>/``) the reference gene model, the iFLAS HQ
isoforms and, when the sample has NGS, per-repeat read coverage (``samtools depth`` -> bedgraph) plus the NGS
junctions of the region as arcs. ``visual_as_merged`` stacks several samples under one reference model
(legacy ``visual_as_merge`` -> ``isoViewer_sample_merged/<gene>/``).
"""

from __future__ import annotations

import glob
import os
from bisect import bisect_left
from itertools import groupby
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing import get_context
from pathlib import Path

from ..io.genepred import GenePredObj
from ..io.intervals import Bed12Plus
from ..runtime import SINGLE_THREAD_ENV, ToolRunner, get_logger


def _split_genes(values=None) -> list[str] | None:
    if values is None:
        return None
    items = []
    if isinstance(values, str):
        values = [values]
    for value in values:
        items.extend(part.strip() for part in str(value).split(",") if part.strip())
    return items or None


def _read_gene_file(path: str | Path | None) -> list[str]:
    if not path:
        return []
    genes = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genes.append(line.split()[0])
    return genes


def _load_isoforms_by_gene(isoform_bed: str | Path) -> dict[str, list[Bed12Plus]]:
    """Group isoforms by gene name: column 14 of ``hqIsoforms/hq.collapsed.bed12+`` (13 = cupcake gene id,
    15 = group key), column 13 of ``refine/isoformGrouped.bed12+``."""
    by_gene: dict[str, list[Bed12Plus]] = {}
    with open(isoform_bed) as f:
        for line in f:
            if not line.strip():
                continue
            iso = Bed12Plus(line.rstrip("\n"))
            gene = iso.record[13] if len(iso.record) >= 14 else iso.record[12]
            by_gene.setdefault(gene, []).append(iso)
    return by_gene


def _bed12(iso: Bed12Plus) -> str:
    """Plain bed12 line (pyGenomeTracks skips bed lines with extra columns)."""
    return "\t".join(map(str, iso.record[:12]))


def _target_genes(genes, gene_file, isoforms_by_gene) -> list[str]:
    requested = (_split_genes(genes) or []) + _read_gene_file(gene_file)
    targets = requested or sorted(isoforms_by_gene)
    if not targets:
        raise ValueError("no genes found to visualize")
    return targets


def _ngs_inputs(base: Path) -> tuple[list[Path], Path | None]:
    """Per-repeat NGS BAMs and the NGS junction bed of a sample dir, if mapping produced them."""
    bams = [Path(p) for p in sorted(glob.glob(str(base / "mapping" / "rna-seq" / "alignment" / "repeat*.sorted.bam")))]
    junctions = base / "mapping" / "rna-seq" / "reassembly" / "junctions.bed"
    return bams, (junctions if junctions.is_file() else None)


def _track_file(path) -> str:
    """The ``file =`` value of a track that sits next to its tracks.ini: pyGenomeTracks splits the value on spaces
    (several files per track), so a path with whitespace is written by name and resolved against tracks.ini."""
    path = str(path)
    return os.path.basename(path) if any(c.isspace() for c in path) else path


def _bed_section(name: str, path: Path, title: str, color: str, display: str) -> str:
    return (f"[{name}]\nfile = {_track_file(path)}\ntitle = {title}\nheight = {{height}}\ncolor = {color}\n"
            f"display = {display}\nlabels = true\nfile_type = bed\n")


def junction_index(junctions_bed) -> dict:
    """NGS ``junctions.bed`` read once: chrom -> (sorted intron starts, rows (start, end, file order, score), widest
    intron). Genes then pick their junctions by bisection instead of re-reading the file (quadratic in genes)."""
    by_chrom = {}
    with open(junctions_bed) as f:
        for order, line in enumerate(f):
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            try:
                sizes, starts = c[10].split(","), c[11].split(",")
                s, e = int(c[1]) + int(sizes[0]), int(c[1]) + int(starts[1])
            except (IndexError, ValueError):
                continue
            by_chrom.setdefault(c[0], []).append((s, e, order, c[4]))
    index = {}
    for chrom, rows in by_chrom.items():
        rows.sort()
        index[chrom] = ([r[0] for r in rows], rows, max(r[1] - r[0] for r in rows))
    return index


def write_ngs_tracks(gene_dir: Path, label: str, chrom: str, start: int, end: int, bams, junctions_bed, runner) -> list[str]:
    """Coverage bedgraph per repeat (``samtools depth``) and the NGS junctions of the region as a links file;
    returns the tracks.ini sections (``label`` distinguishes samples in merged plots). ``junctions_bed`` is a
    junctions.bed path or its :func:`junction_index`."""
    sections = []
    for i, bam in enumerate(bams):
        depth = runner.run(["samtools", "depth", "-a", "-r", f"{chrom}:{start + 1}-{end}", str(bam)]).stdout or ""
        bedgraph = gene_dir / f"repeat{i}.{label}.bedgraph"
        with open(bedgraph, "w") as out:
            # one interval per run of equal depth (-a lists every position): a 750-kb maize gene wrote 750k lines,
            # mostly intronic zeros, and took pyGenomeTracks over 40 CPU minutes to draw
            for d, run in groupby((line.split("\t")[:3] for line in depth.splitlines()), key=lambda c: c[2]):
                run = list(run)
                out.write(f"{run[0][0]}\t{int(run[0][1]) - 1}\t{run[-1][1]}\t{d}\n")
        sections.append(f"[coverage_repeat{i}_{label}]\nfile = {_track_file(bedgraph)}\ntitle = NGS coverage repeat{i} ({label})\n"
                        "height = 2\ncolor = #54a24b\nmin_value = 0\nfile_type = bedgraph\n")
    if junctions_bed:
        index = junctions_bed if isinstance(junctions_bed, dict) else junction_index(junctions_bed)
        hits = []
        if chrom in index:
            starts, rows, widest = index[chrom]
            # an intron ending after `start` must begin at or after `start - widest`; output keeps the file order
            window = rows[bisect_left(starts, start - widest):bisect_left(starts, end)]
            hits = sorted((r for r in window if r[1] > start), key=lambda r: r[2])
        links = gene_dir / f"{label}.junctions.links"
        with open(links, "w") as out:
            out.writelines(f"{chrom}\t{s - 1}\t{s}\t{chrom}\t{e}\t{e + 1}\t{score}\n" for s, e, _, score in hits)
        if hits:
            sections.append(f"[junctions_{label}]\nfile = {_track_file(links)}\ntitle = NGS junctions ({label})\nheight = 2\n"
                            "links_type = arcs\ncolor = #e45756\nline_width = 1\nfile_type = links\n")
    return sections


def _write_tracks(tracks: Path, sections: list[str]) -> None:
    tracks.write_text("[x-axis]\nwhere = top\n\n" + "\n".join(sections))


def _render_worker_init() -> None:
    os.environ.update(SINGLE_THREAD_ENV)


def _render(job) -> str:
    """One plot through pyGenomeTracks' own entry point, in a long-lived worker (imports paid once per worker)."""
    from pygenometracks import plotTracks

    tracks, region, pdf = job
    try:
        plotTracks.main(["--tracks", str(tracks), "--region", region, "--outFileName", str(pdf)])
    except SystemExit as exc:                                   # argparse / pyGenomeTracks errors exit
        raise RuntimeError(f"pyGenomeTracks failed for {tracks} ({region}): exit {exc.code}") from None
    return str(pdf)


def render_all(jobs, threads: int) -> None:
    """Render (tracks, region, pdf) jobs in ``threads`` spawned worker processes."""
    jobs = list(jobs)
    if not jobs:
        return
    with ProcessPoolExecutor(max_workers=max(1, min(int(threads or 1), len(jobs))), mp_context=get_context("spawn"),
                             initializer=_render_worker_init) as pool:
        list(pool.map(_render, jobs))      # one job at a time: a long gene must not hold up the genes queued behind it


def _plot(runner, tracks: Path, region: str, pdf: Path) -> None:
    runner.run(["pyGenomeTracks", "--tracks", str(tracks), "--region", region, "--outFileName", str(pdf)], env=SINGLE_THREAD_ENV)


def _reference_track(gene_dir: Path, gene: str, ref_models) -> tuple[Path | None, str]:
    """Write the reference gene model bed (empty for genes/novel clusters absent from the annotation)."""
    ref_bed = gene_dir / f"{gene}.reference.bed"
    ref_bed.write_text("".join(repr(model.go_to_bed()) + "\n" for model in ref_models))
    if not ref_models:
        return None, ""
    return ref_bed, _bed_section("reference", ref_bed, "reference gene model", "#4c78a8", "collapsed").format(height=3)


def build_visual_as_tracks(base_dir, ref_gpe, *, genes=None, gene_file=None, out_dir=None, isoform_bed=None,
                           runner=None, run_plot=False, ngs_bams=None, ngs_junctions=None, ref_index=None,
                           isoform_index=None) -> dict[str, str]:
    """Build (and optionally render) the pyGenomeTracks plot of one target gene for one sample. ``ref_index`` (the
    parsed reference, ``GenePredObj(...).geneName2gpeObj``) and ``isoform_index`` (``_load_isoforms_by_gene``) are
    passed by callers that plot many genes, so the files are parsed once rather than once per gene."""
    base = Path(base_dir)
    isoform_bed = Path(isoform_bed) if isoform_bed else base / "hqIsoforms" / "hq.collapsed.bed12+"
    requested = (_split_genes(genes) or []) + _read_gene_file(gene_file)
    if len(requested) != 1:
        raise ValueError("build_visual_as_tracks expects exactly one target gene")
    gene = requested[0]
    if ref_index is None:
        ref_index = GenePredObj(str(ref_gpe), False).geneName2gpeObj
    ref_models = ref_index.get(gene, [])
    isoforms = (isoform_index if isoform_index is not None else _load_isoforms_by_gene(isoform_bed)).get(gene, [])
    if not ref_models and not isoforms:
        raise KeyError(f"gene not found in ref_gpe or in the isoform bed: {gene}")
    gene_dir = (Path(out_dir) if out_dir else base / "visual_as") / gene
    gene_dir.mkdir(parents=True, exist_ok=True)
    runner = runner or ToolRunner()

    ref_bed, ref_section = _reference_track(gene_dir, gene, ref_models)
    iflas_bed = gene_dir / f"{gene}.iflas_isoforms.bed"
    iflas_bed.write_text("".join(_bed12(isoform) + "\n" for isoform in isoforms))
    starts = [m.txStart for m in ref_models] + [iso.chromStart for iso in isoforms]
    ends = [m.txEnd for m in ref_models] + [iso.chromEnd for iso in isoforms]
    chrom = ref_models[0].chrom if ref_models else isoforms[0].chrom
    start, end = min(starts), max(ends)
    region = f"{chrom}:{start}-{end}"

    sections = ([ref_section] if ref_section else []) + [
        _bed_section("iflas_isoforms", iflas_bed, "iFLAS HQ isoforms", "#f58518", "stacked").format(height=5)]
    ngs_sections = write_ngs_tracks(gene_dir, base.name, chrom, start, end, ngs_bams or [], ngs_junctions, runner)
    tracks, pdf = gene_dir / "tracks.ini", gene_dir / f"{gene}.pdf"
    _write_tracks(tracks, sections + ngs_sections)
    if run_plot:
        _plot(runner, tracks, region, pdf)
    return {"gene": gene, "gene_dir": str(gene_dir), "tracks": str(tracks), "reference_bed": str(ref_bed or ""),
            "isoform_bed": str(iflas_bed), "pdf": str(pdf), "region": region, "ngs_tracks": len(ngs_sections)}


def visual_as_sample(base_dir, ref_gpe, *, genes=None, gene_file=None, out_dir=None, isoform_bed=None, runner=None,
                     run_plot=True, ngs=None, threads=1) -> list[dict[str, str]]:
    """Plot every target gene (default: all genes with HQ isoforms, like legacy) of one sample; NGS coverage
    and junction tracks are added when ``ngs`` is None/True and the mapping stage produced them."""
    base = Path(base_dir)
    isoform_bed = Path(isoform_bed) if isoform_bed else base / "hqIsoforms" / "hq.collapsed.bed12+"
    isoforms_by_gene = _load_isoforms_by_gene(isoform_bed)
    targets = _target_genes(genes, gene_file, isoforms_by_gene)
    bams, junctions = _ngs_inputs(base) if ngs in (None, True) else ([], None)
    in_process = run_plot and runner is None          # a caller-supplied runner keeps one pyGenomeTracks command per gene
    runner = runner or ToolRunner()
    ref_index = GenePredObj(str(ref_gpe), False).geneName2gpeObj
    junctions = junction_index(junctions) if junctions else None

    def one(gene):
        return build_visual_as_tracks(base, ref_gpe, genes=[gene], out_dir=out_dir, isoform_bed=isoform_bed, runner=runner,
                                      run_plot=run_plot and not in_process, ngs_bams=bams, ngs_junctions=junctions,
                                      ref_index=ref_index, isoform_index=isoforms_by_gene)

    with ThreadPoolExecutor(max_workers=max(1, int(threads or 1))) as pool:
        outputs = list(pool.map(one, targets))
    if in_process:
        render_all([(o["tracks"], o["region"], o["pdf"]) for o in outputs], threads)
    return outputs


def visual_as_merged(base_dirs, ref_gpe, *, genes=None, gene_file=None, out_dir=None, runner=None, run_plot=True,
                     ngs=None, threads=1) -> list[dict[str, str]]:
    """Legacy ``visual_as_merge``: one plot per target gene with every sample's HQ isoforms (and NGS tracks)
    under the reference model, in ``<out_dir>/<gene>/`` (default ``isoViewer_sample_merged`` next to the samples'
    project dir). Samples are labelled by their base-dir name; samples without HQ isoforms for a gene are skipped."""
    bases = [Path(b) for b in base_dirs]
    out_root = Path(out_dir) if out_dir else bases[0].parent.parent / "isoViewer_sample_merged"
    ref = GenePredObj(str(ref_gpe), False).geneName2gpeObj
    per_sample = {b: _load_isoforms_by_gene(b / "hqIsoforms" / "hq.collapsed.bed12+") for b in bases}
    all_genes = {g for isos in per_sample.values() for g in isos}
    targets = _target_genes(genes, gene_file, {g: None for g in all_genes})
    in_process = run_plot and runner is None
    runner = runner or ToolRunner()
    get_logger().info("visual_as merged: %d gene(s) across %d sample(s)", len(targets), len(bases))
    ngs_inputs = {}
    if ngs in (None, True):
        for b in bases:
            bams, junctions = _ngs_inputs(b)
            ngs_inputs[b] = (bams, junction_index(junctions) if junctions else None)
    render_jobs = []

    outputs = []
    for gene in targets:
        ref_models = ref.get(gene, [])
        with_gene = [b for b in bases if per_sample[b].get(gene)]
        if not ref_models and not with_gene:
            raise KeyError(f"gene not found in ref_gpe or in any sample: {gene}")
        gene_dir = out_root / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        ref_bed, ref_section = _reference_track(gene_dir, gene, ref_models)
        sections = [ref_section] if ref_section else []
        starts = [m.txStart for m in ref_models]
        ends = [m.txEnd for m in ref_models]
        chrom = ref_models[0].chrom if ref_models else per_sample[with_gene[0]][gene][0].chrom
        for b in with_gene:
            isos = per_sample[b][gene]
            bed = gene_dir / f"{gene}.{b.name}.iflas_isoforms.bed"
            bed.write_text("".join(_bed12(iso) + "\n" for iso in isos))
            starts += [iso.chromStart for iso in isos]
            ends += [iso.chromEnd for iso in isos]
            sections.append(_bed_section(f"isoforms_{b.name}", bed, f"iFLAS HQ isoforms ({b.name})", "#f58518", "stacked").format(height=4))
        start, end = min(starts), max(ends)
        for b in with_gene if ngs_inputs else ():
            bams, junctions = ngs_inputs[b]
            sections += write_ngs_tracks(gene_dir, b.name, chrom, start, end, bams, junctions, runner)
        tracks, pdf, region = gene_dir / "tracks.ini", gene_dir / f"{gene}.pdf", f"{chrom}:{start}-{end}"
        _write_tracks(tracks, sections)
        if in_process:
            render_jobs.append((tracks, region, pdf))
        elif run_plot:
            _plot(runner, tracks, region, pdf)
        outputs.append({"gene": gene, "gene_dir": str(gene_dir), "tracks": str(tracks), "pdf": str(pdf), "region": region,
                        "samples": [b.name for b in with_gene]})
    render_all(render_jobs, threads)
    return outputs
