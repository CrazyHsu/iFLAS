"""Real-reference check (2026-10-16): reads simulated from a real genome and annotation by scripts/make_realref_data.py.

The generator test runs offline on synthetic v2 standing in for the reference. The checks of a whole-chain run on the
real maize reference (MaizeGDB B73 v5) run only when IFLAS_REALREF names the run directory (scripts/run_realref.sh).
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import os
import shutil
import sys

import pysam
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
V2 = os.path.join(ROOT, "test_data", "synthetic_v2")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import make_realref_data  # noqa: E402

from iflas.config import load_config  # noqa: E402
from iflas.refprep import init_ref  # noqa: E402


def _tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="needs bedtools (core iflas env)")
def test_generator_on_a_stand_in_reference(tmp_path):
    """Deterministic for a seed; the truth tables agree with the reads written; the configs load."""
    ref = init_ref(os.path.join(V2, "genome.fa"), gtf=os.path.join(V2, "genes.gtf"), outdir=str(tmp_path / "ref"), build_minimap2=False)
    base = ref["ref_gpe"][:-4]
    runs = []
    for name in ("a", "b"):
        out = tmp_path / name
        make_realref_data.main([base, str(out), "--genome", ref["ref_genome"], "--gtf", ref["ref_gtf"], "--threads", "2"])
        runs.append(out)
    a, b = runs
    for rel in ("pacbio/hifi_reads.bam", "ont/rna002/reads.fastq", "ont/rna004/dorado.bam", "ngs/condB_rep2_R2.fastq.gz",
                "truth/isoforms.tsv", "truth/designs.tsv", "truth/expected_as.tsv", "samples.tsv"):
        assert _sha(a / rel) == _sha(b / rel), rel

    isoforms = _tsv(a / "truth" / "isoforms.tsv")
    pacbio = _tsv(a / "pacbio" / "read_truth.tsv")
    with pysam.AlignmentFile(str(a / "pacbio" / "hifi_reads.bam"), check_sq=False) as bam:
        assert sum(1 for _ in bam) == len(pacbio) == sum(int(i["n_flnc"]) for i in isoforms)
        assert "READTYPE=CCS" in bam.header.to_dict()["RG"][0]["DS"]
    ont = {chem: _tsv(a / "ont" / chem / "read_truth.tsv") for chem in ("rna002", "rna004")}
    recovered = sum(int(r["flnc_recovered"]) for r in pacbio)
    assert sum(int(i["long_reads"]) for i in isoforms) == recovered + len(ont["rna002"]) + len(ont["rna004"])
    assert open(a / "ont" / "rna002" / "reads.fastq").read().count("\n") == 4 * len(ont["rna002"])
    assert os.path.isfile(a / "ont" / "rna002" / "signal" / "fast5" / "batch0.fast5")
    with gzip.open(a / "ngs" / "condA_rep1_R1.fastq.gz", "rt") as f:
        pairs = sum(1 for _ in f) // 4                    # transcripts shorter than a pair + 20 nt yield none
    assert 0.9 * sum(int(i["ngs_condA_rep1"]) for i in isoforms) <= pairs <= sum(int(i["ngs_condA_rep1"]) for i in isoforms)

    designs = _tsv(a / "truth" / "designs.tsv")
    assert {d["design"] for d in designs} == set(make_realref_data.DESIGNS)
    assert all(d["asas_expected"] == "1" for d in designs if d["design"] == "asas")
    expected = _tsv(a / "truth" / "expected_as.tsv")
    assert {e["type"] for e in expected} >= {"SE", "IR", "A5SS", "A3SS"}
    for config, samples in (("config.yaml", "samples.tsv"), ("config_conditions.yaml", "samples_conditions.tsv")):
        cfg = load_config(str(a / config), str(a / samples))
        assert {s.condition for s in cfg.samples} == {"condA", "condB"}


# ----------------------------------------------------------------------------- a run of scripts/run_realref.sh
REALREF = os.environ.get("IFLAS_REALREF")
needs_run = pytest.mark.skipif(not REALREF, reason="set IFLAS_REALREF to a dataset run by scripts/run_realref.sh")
AS_FILES = {"SE": "SE.bed12+", "IR": "IR.bed6+", "A5SS": "A5SS.bed6+", "A3SS": "A3SS.bed6+"}
# measured on MaizeGDB B73 v5 (2026-10-16): isoform recovery 92.9 % (long genes 92.6 %),
# event recall 92.0 % with 2.0 % extra events; the floors leave room for the PU filter's run-to-run noise
RECOVERY, EVENT_RECALL, EVENT_EXTRA = 0.88, 0.88, 0.05


def _rows(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f if line.strip()]


def _base():
    return os.path.join(REALREF, "iflas_out", "realref", "B73_all")


def _intron_chain(bed):
    s, sizes, starts = int(bed[1]), [int(x) for x in bed[10].strip(",").split(",")], [int(x) for x in bed[11].strip(",").split(",")]
    return ";".join(f"{s + a + n}-{s + b}" for a, n, b in zip(starts, sizes, starts[1:]))


def _hq_by_chain():
    return {(c[0], c[5], _intron_chain(c)): c for c in _rows(os.path.join(_base(), "hqIsoforms", "hq.collapsed.bed12+")) if int(c[9]) > 1}


def _truth(name):
    return _tsv(os.path.join(REALREF, "truth", name))


@needs_run
def test_realref_every_step_ran():
    """Every step exits 0 without a traceback; the merged run runs all eleven per-sample stages (diff_as runs on its own
    per-library run)."""
    import re
    for step in ("all", "report", "cond_mapping", "cond_lr", "diff_as", "report_diff"):
        assert "Exit status: 0" in open(os.path.join(REALREF, "logs", f"{step}.time")).read(), step
        assert "Traceback" not in open(os.path.join(REALREF, "logs", f"{step}.log")).read(), step
    lines = {m.group(1): m.group(0) for m in re.finditer(r"(?m)^  ([a-z_]+): .*$", open(os.path.join(REALREF, "logs", "all.log")).read())}
    for stage in ("preproc", "mapping", "collapse", "refine", "pu_filter", "find_as", "find_pa", "characterize", "visual_as", "asas", "ase", "palen_as"):
        assert "skipped" not in lines[stage], lines[stage]


@needs_run
def test_realref_keeps_the_reference_names():
    """Chromosome, scaffold and gene IDs come out as the reference has them."""
    import re
    ref = load_config(os.path.join(REALREF, "config.yaml")).reference[0]
    names = {line.split("\t")[0] for line in open(ref.ref_size)}
    genes = {line.split("\t")[11] for line in open(ref.ref_gpe)}            # reads may be assigned to any overlapping gene
    hq = _rows(os.path.join(_base(), "hqIsoforms", "hq.collapsed.bed12+"))
    assert {c[0] for c in hq} <= names
    assert any(c[0].startswith("scaf") for c in hq) or not any(n.startswith("scaf") for n in names)
    novel_cluster = re.compile(r"^(.+):[+-]:\d+$")                         # refine's id for reads of no reference gene
    assert all(c[13] in genes or (novel_cluster.match(c[13]) and novel_cluster.match(c[13]).group(1) in names) for c in hq)


@needs_run
def test_realref_recovers_the_simulated_isoforms():
    """Truth multi-exon isoforms with >= 10 long reads come back as HQ isoforms with their exact intron chain, the
    annotated ones labelled annotated; so do those of the long genes BUG-193 used to drop. The misses traced on this
    reference are aligner placements (equivalent or nearby splice sites, terminal exons of 2-30 nt, local duplications),
    collapse merging annotated sites <= 5 bp apart, and introns over minimap2's 50-kb limit."""
    hq = _hq_by_chain()
    genes = {g["gene"]: g for g in _truth("genes.tsv")}
    want = [i for i in _truth("isoforms.tsv") if i["junc_chain"] and int(i["long_reads"]) >= 10]
    got = [i for i in want if (i["chrom"], i["strand"], i["junc_chain"]) in hq]
    assert len(got) >= RECOVERY * len(want), (len(got), len(want))
    assert all(hq[(i["chrom"], i["strand"], i["junc_chain"])][18] == "annotated" for i in got if i["annotated"] == "1")
    long = [i for i in want if "long_gene" in genes[i["gene"]]["flags"]]
    assert long and sum((i["chrom"], i["strand"], i["junc_chain"]) in hq for i in long) >= RECOVERY * len(long)


@needs_run
def test_realref_finds_the_expected_as_events():
    """The AS events find_as reports on the simulated isoforms (expected_as.tsv) are found from the reads."""
    expected = {(e["type"], e["chrom"], int(e["start"]), int(e["end"])) for e in _truth("expected_as.tsv")}
    found = {(t, c[0], int(c[1]), int(c[2])) for t, f in AS_FILES.items()
             for c in _rows(os.path.join(_base(), "as_events", "ordinary_as", "LR", f))}
    assert len(expected & found) >= EVENT_RECALL * len(expected), (len(expected & found), len(expected))
    assert len(found - expected) <= EVENT_EXTRA * len(found), (len(found - expected), len(found))


def _designs(kind):
    return [d for d in _truth("designs.tsv") if d["design"] == kind]


@needs_run
def test_realref_calls_the_designed_allele_effects():
    """Allele-specific splicing and expression designed on two-isoform genes: most ASAS designs are called (the others
    lose an isoform before HQ, have a single-exon isoform (R-2) or miss the published >= 5 / <= 3 rule), every ASE design
    shows its imbalance in the FL reads, and few ASAS calls fall on genes without a variant (BUG-201)."""
    called = {c[2] for c in _rows(os.path.join(_base(), "asas", "partialAsRelatedHaplotype.txt"))}
    designed = {d["gene"] for d in _designs("asas")}
    assert len(called & designed) >= 0.6 * len(designed) and len(called - designed) <= 0.1 * len(called), (called & designed, called - designed)
    rows = {}
    for c in _rows(os.path.join(_base(), "ase", "repeat0.mergedASEfile.txt")):
        rows.setdefault(c[1], []).append(c)
    imbalanced = [d["gene"] for d in _designs("ase") if any(float(c[2]) + float(c[3]) and abs(float(c[2]) / (float(c[2]) + float(c[3])) - 0.5) >= 0.25
                                                            for c in rows.get(d["gene"], []))]
    assert len(imbalanced) >= 0.9 * len(_designs("ase"))


@needs_run
def test_realref_calls_the_designed_tail_length_shifts():
    """Tails differ by design between the isoforms of 50 genes (RNA002 via nanopolish on the simulated signal, RNA004 via
    dorado's pt:i): 39 are called; the rest lose an isoform before HQ or have a single-exon isoform (R-2). Before
    BUG-200 the strain's dorado table made palen_as skip nanopolish, and only 33 were."""
    import glob
    called = {c[13] for f in glob.glob(os.path.join(_base(), "palenAS", "mergeByJunc", "*.palenAndAS.sig.bed12+")) for c in _rows(f)}
    designed = {d["gene"] for d in _designs("palen")}
    assert len(called & designed) >= 0.7 * len(designed), sorted(designed - called)
    assert len(called - designed) <= 2, sorted(called - designed)


@needs_run
def test_realref_finds_the_designed_pa_sites_and_differential_splicing():
    """The second PA site 250 bp downstream of the canonical one is found; rMATS calls every designed differential event
    it can define from the merged HQ GTF (R-14: same flanking exons) and none of the equal-PSI controls."""
    import glob
    pa = {}
    for c in _rows(os.path.join(_base(), "as_events", "pa", "PA.bed6+")):
        pa.setdefault((c[0], c[5]), []).append((int(c[1]), int(c[2])))
    isoforms = {i["isoform"]: i for i in _truth("isoforms.tsv")}
    found = 0
    for d in _designs("apa"):
        full = isoforms[d["isoforms"].split(";")[0]]
        exons = [tuple(map(int, e.split("-"))) for e in full["exons"].split(";")]
        site = exons[-1][1] + 250 if full["strand"] == "+" else exons[0][0] - 250
        found += any(a - 30 <= site <= b + 30 for a, b in pa.get((full["chrom"], full["strand"]), []))
    assert found >= 0.9 * len(_designs("apa"))
    das = os.path.join(REALREF, "iflas_out_conditions", "das")
    comp = next(p for p in glob.glob(os.path.join(das, "*_vs_*")) if os.path.isdir(p) and not p.endswith((".sigDiffAS", ".tmp")))
    defined = {r["GeneID"].strip('"') for t in ("SE", "RI", "A3SS", "A5SS", "MXE") if os.path.isfile(os.path.join(comp, f"fromGTF.{t}.txt"))
               for r in _tsv(os.path.join(comp, f"fromGTF.{t}.txt"))}
    called = {r["GeneID"].strip('"') for f in glob.glob(os.path.join(das, "*.sigDiffAS", "*.sig.txt")) for r in _tsv(f)}
    diff, null = {d["gene"] for d in _designs("diff")}, {d["gene"] for d in _designs("diff_null")}
    assert diff & defined and diff & defined <= called and not null & called
    assert os.path.isfile(os.path.join(comp + ".sigDiffAS", "GO", "sigDiff.goEnrichResults.txt"))


@needs_run
def test_realref_reports_every_plot():
    sample = os.path.join(REALREF, "iflas_out", "reports", "realref_B73_all")
    for folder in ("basicStatistics", "asPattern", "geneStrucPlots", "alleleAsPlots", "palenAsPlots", "palenApaPlots"):
        assert any(os.path.getsize(os.path.join(sample, folder, f)) for f in os.listdir(os.path.join(sample, folder)) if f.endswith(".pdf")), folder
    assert os.path.getsize(os.path.join(REALREF, "iflas_out", "reports", "html", "iflas_report.html"))
