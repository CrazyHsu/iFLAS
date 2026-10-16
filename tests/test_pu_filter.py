"""pu_filter stage — integration tests on synthetic data.

getInputBed test (gtfToGenePred + bedtools): runs refine on the synthetic collapsed
GFF against a SUBSET reference (canonical iso1 per gene only), so the AS-variant
isoforms become "novel" while the canonical ones are "annotated" — the class balance
the PU learner needs. Asserts getInputBed's annotation split + per-junction-chain output.

Full-chain test (cupcake + NGS + LIBSVM): minimap2 -> cupcake collapse -> refine ->
align NGS reassembly (tmp.bam + junctions.bed) -> filter_hq_isoforms, asserting the
hqIsoforms/hq.collapsed.bed12+ file is produced. Route B: verifies the chain RUNS;
numerical PU equivalence to the published result still needs the golden baseline.
"""

import os
import shutil
import subprocess
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_synthetic_data as syn  # noqa: E402

from iflas.io.genepred import GenePredObj
from iflas.pu.iso_pu import getInputBed
from iflas.stages.refine import refine_isoforms


_TOOLS_BASIC = shutil.which("gtfToGenePred") and shutil.which("bedtools")
MIN_PU_NOVEL_MULTI_EXON = 20


def _write_collapsed_gff(path, only=None):
    """All designed isoforms as a cupcake-style collapsed GFF (1-based inclusive)."""
    with open(path, "w") as f:
        for gene, isos in syn.GENES.items():
            for name, chrom, strand, exons in isos:
                if only and name not in only:
                    continue
                attrs = f'gene_id "{gene}"; transcript_id "{name}";'
                f.write(f"{chrom}\tPacBio\ttranscript\t{exons[0][0] + 1}\t{exons[-1][1]}\t.\t{strand}\t.\t{attrs}\n")
                for (st, en) in exons:
                    f.write(f"{chrom}\tPacBio\texon\t{st + 1}\t{en}\t.\t{strand}\t.\t{attrs}\n")


def _write_subset_gtf(path):
    """Reference = canonical iso1 per gene only -> AS variants are 'novel'."""
    with open(path, "w") as f:
        for gene, isos in syn.GENES.items():
            for name, chrom, strand, exons in isos:
                if not name.endswith(".iso1"):
                    continue
                for (st, en) in exons:
                    f.write(f'{chrom}\tsyn\texon\t{st + 1}\t{en}\t.\t{strand}\t.\t'
                            f'gene_id "{gene}"; transcript_id "{name}";\n')


def _build_refs(gtf, ref_gpe, ref_bed):
    subprocess.run(["gtfToGenePred", "-genePredExt", gtf, ref_gpe], check=True)
    with open(ref_bed, "w") as f:
        for line in GenePredObj(ref_gpe, bincolumn=False).toBed(gene=True):
            f.write(str(line) + "\n")


def _write_group(path):
    """collapse group.txt: isoform -> comma-separated supporting read ids (5 per isoform)."""
    with open(path, "w") as f:
        for gene, isos in syn.GENES.items():
            for name, chrom, strand, exons in isos:
                reads = ",".join(f"{name}/ccs/{k}" for k in range(5))
                f.write(f"{name}\t{reads}\n")


def _rows(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [l.rstrip("\n").split("\t") for l in f if l.strip()]


def _expected_getinputbed_counts():
    annotated = 0
    novel = 0
    for isos in syn.GENES.values():
        for name, chrom, strand, exons in isos:
            if len(exons) < 2:
                continue
            if name.endswith(".iso1"):
                annotated += 1
            else:
                novel += 1
    return annotated, novel


@pytest.mark.skipif(not _TOOLS_BASIC, reason="needs gtfToGenePred + bedtools (core iflas env)")
def test_getinputbed_annotation_split(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    gff = os.path.join(d, "collapsed.gff")
    _write_collapsed_gff(gff)

    # refine uses the FULL annotation (gene assignment); annotation-ref for the PU
    # anno/novel split is a SUBSET (canonical iso1 only) so AS variants read as novel.
    full_gpe = os.path.join(d, "ref.full.gpe")
    full_bed = os.path.join(d, "ref.full.bed12+")
    _build_refs(os.path.join(d, "genes.gtf"), full_gpe, full_bed)
    subset_gtf = os.path.join(d, "ref.subset.gtf")
    _write_subset_gtf(subset_gtf)
    subset_bed = os.path.join(d, "ref.subset.bed12+")
    _build_refs(subset_gtf, os.path.join(d, "ref.subset.gpe"), subset_bed)

    # refine writes <d>/refine/tofu.collapsed.assigned.unambi.bed12+ (getInputBed's input)
    refine_isoforms(gff, f"{d}/genome.fa", full_gpe, full_bed, d)
    os.makedirs(os.path.join(d, "collapse"), exist_ok=True)
    _write_group(os.path.join(d, "collapse", "tofu.collapsed.group.txt"))

    prev = os.getcwd()
    os.chdir(d)
    try:
        input_bed = getInputBed(baseDir=d, refBed=subset_bed)
    finally:
        os.chdir(prev)

    rows = _rows(input_bed)
    # one row per multi-exon junction chain; the single-exon geneD.iso2_IR is excluded
    by_name = {r[3]: r for r in rows}
    assert "geneA.iso1" in by_name
    # last column is the annotation label
    anno = {r[3]: r[-1] for r in rows}
    # canonical iso1 chains are in the subset reference -> annotated
    assert anno["geneB.iso1"] == "annotated"
    # AS-variant chains are NOT in the subset reference -> novel
    assert anno["geneB.iso2_SE"] == "novel"
    assert anno["geneC.iso2_A5SS"] == "novel"
    assert anno["geneE.iso2_A3SS"] == "novel"
    # both classes present (PU learning needs positives + unlabeled)
    labels = set(anno.values())
    assert labels == {"annotated", "novel"}
    expected_annotated, expected_novel = _expected_getinputbed_counts()
    assert sum(1 for label in anno.values() if label == "annotated") == expected_annotated
    assert sum(1 for label in anno.values() if label == "novel") == expected_novel
    # The synthetic set must stay large enough for bagging PU OOB scoring.
    assert expected_novel >= MIN_PU_NOVEL_MULTI_EXON


def test_filter_hq_isoforms_uses_runner_for_final_filter(tmp_path, monkeypatch):
    """The stage driver should route the final filter.pl pipeline through ToolRunner."""
    from iflas.stages import pu_filter

    base = tmp_path / "base"
    base.mkdir()
    genome = tmp_path / "genome.fa"
    ref_bed = tmp_path / "ref.bed12+"
    input_bed = tmp_path / "input.bed12+"
    genome.write_text(">chr1\nACGT\n")
    ref_bed.write_text("")
    input_bed.write_text("")

    def fake_get_input_bed(baseDir, refBed):
        assert baseDir == str(base)
        assert refBed == str(ref_bed)
        return str(input_bed)

    feature_data = pd.DataFrame(
        {
            "isoLength": [1000, 1200],
            "flCount": [3, 2],
            "ratioIsoToGene": [0.8, 0.4],
            "exonNum": [2, 2],
            "GC": [50.0, 48.0],
            "orfLength": [600, 480],
            "orfIntegrity": [1, 1],
            "pepLength": [200, 160],
            "FickettScore": [0.5, 0.4],
            "pI": [6.0, 6.2],
            "codingP": [0.9, 0.3],
            "canJuncRatio": [1.0, 1.0],
            "sdJuncCov": [0.1, 0.2],
            "minJuncRPKM": [1.0, 0.8],
            "nIndelsAroundJunc": [0, 0],
            "ratioMinJuncCovToAllCov": [0.7, 0.6],
            "nJuncsWithIndels": [0, 0],
            "indelNearJunc": [False, False],
            "annotation": ["annotated", "novel"],
        },
        index=["anno1", "novel1"],
    )

    def fake_process_input_file(inputBed, genomeFa, refBed, junctionFile=None, samFile=None, runner=None):
        assert inputBed == str(input_bed)
        assert genomeFa == str(genome)
        assert refBed == str(ref_bed)
        assert runner is recording_runner
        return feature_data

    class DummySchemeModelEval:
        def __init__(self, trainingData, model, scheme, random_state=0):
            self.trainingData = trainingData
            self.model = model
            self.scheme = scheme
            self.predResults = pd.DataFrame(
                {"label": [1, 0], "pu_score": [1.0, 0.75]},
                index=["anno1", "novel1"],
            )

        def eval(self):
            assert self.model == "GB"
            assert self.scheme == "bagging"

        def filterIsoformsByScore(self, filterScore, outFile, lqAnnoIso):
            assert filterScore == 0.5
            with open(outFile, "w") as f:
                f.write("isoform\tlabel\tpu_score\nnovel1\t0\t0.75\n")

    class RecordingRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs, os.getcwd()))
            with open("hq.collapsed.bed12+", "w") as f:
                f.write("hq\n")

    recording_runner = RecordingRunner()
    monkeypatch.setattr(pu_filter, "getInputBed", fake_get_input_bed)
    monkeypatch.setattr(pu_filter, "processInputFile", fake_process_input_file)
    monkeypatch.setattr(pu_filter, "SchemeModelEval", DummySchemeModelEval)
    monkeypatch.setattr(pu_filter, "_find_utils_dir", lambda: "/repo/utils")

    out = pu_filter.filter_hq_isoforms(
        str(base), str(genome), str(ref_bed), runner=recording_runner
    )

    assert out["hq_bed"] == str(base / "hqIsoforms" / "hq.collapsed.bed12+")
    assert os.path.isfile(out["hq_bed"])
    assert len(recording_runner.calls) == 1
    cmd, kwargs, cwd = recording_runner.calls[0]
    assert cwd == str(base / "hqIsoforms")
    assert kwargs == {}
    assert "/repo/utils/filter.pl" in cmd
    assert str(input_bed) in cmd
    assert cmd.endswith("> hq.collapsed.bed12+")


def test_filter_hq_isoforms_reads_precomputed_feature_file(tmp_path, monkeypatch):
    """A configured feature_file should be read as isoFeatures data, not recomputed."""
    from iflas.stages import pu_filter

    base = tmp_path / "base"
    base.mkdir()
    genome = tmp_path / "genome.fa"
    ref_bed = tmp_path / "ref.bed12+"
    input_bed = tmp_path / "input.bed12+"
    feature_file = tmp_path / "precomputed_features.tsv"
    genome.write_text(">chr1\nACGT\n")
    ref_bed.write_text("")
    input_bed.write_text("")

    feature_data = pd.DataFrame(
        {
            "isoLength": [1000, 1200],
            "flCount": [3, 2],
            "ratioIsoToGene": [0.8, 0.4],
            "exonNum": [2, 2],
            "annotation": ["annotated", "novel"],
            "GC": [50.0, 48.0],
            "orfLength": [600, 480],
            "orfIntegrity": [1, 1],
            "pepLength": [200, 160],
            "FickettScore": [0.5, 0.4],
            "pI": [6.0, 6.2],
            "codingP": [0.9, 0.3],
            "codingLabel": ["coding", "noncoding"],
            "canJuncRatio": [1.0, 1.0],
            "sdJuncCov": [0.1, 0.2],
            "withNovelJunc": [False, True],
            "minNovelJuncRPKM": [0.0, 0.7],
            "bite": [False, False],
            "minJuncRPKM": [1.0, 0.8],
            "nIndelsAroundJunc": [0, 0],
            "ratioMinJuncCovToAllCov": [0.7, 0.6],
            "nJuncsWithIndels": [0, 0],
            "indelNearJunc": [False, False],
        },
        index=["anno1", "novel1"],
    )
    feature_data.to_csv(feature_file, sep="\t", index=True, index_label="isoform")

    def fake_get_input_bed(baseDir, refBed):
        assert baseDir == str(base)
        assert refBed == str(ref_bed)
        return str(input_bed)

    def fail_process_input_file(*args, **kwargs):
        raise AssertionError("feature_file path must not recompute PU features")

    class DummySchemeModelEval:
        def __init__(self, trainingData, model, scheme, random_state=0):
            self.trainingData = trainingData
            self.model = model
            self.scheme = scheme
            self.predResults = pd.DataFrame(
                {"label": [1, 0], "pu_score": [1.0, 0.75]},
                index=["anno1", "novel1"],
            )

        def eval(self):
            assert self.model == "DT"
            assert self.scheme == "bagging"
            assert set(self.trainingData.index) == {"anno1", "novel1"}
            assert "annotation" not in self.trainingData.columns

        def filterIsoformsByScore(self, filterScore, outFile, lqAnnoIso):
            assert filterScore == 0.25
            with open(outFile, "w") as f:
                f.write("isoform\tlabel\tpu_score\nnovel1\t0\t0.75\n")

    class RecordingRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs, os.getcwd()))
            assert os.path.isfile("isoFeatures.txt")
            with open("hq.collapsed.bed12+", "w") as f:
                f.write("hq\n")

    recording_runner = RecordingRunner()
    monkeypatch.setattr(pu_filter, "getInputBed", fake_get_input_bed)
    monkeypatch.setattr(pu_filter, "processInputFile", fail_process_input_file)
    monkeypatch.setattr(pu_filter, "SchemeModelEval", DummySchemeModelEval)
    monkeypatch.setattr(pu_filter, "_find_utils_dir", lambda: "/repo/utils")

    out = pu_filter.filter_hq_isoforms(
        str(base), str(genome), str(ref_bed), feature_file=str(feature_file),
        filter_score=0.25, model="DT", runner=recording_runner,
    )

    assert out["features"] == str(base / "hqIsoforms" / "isoFeatures.txt")
    assert pd.read_csv(out["features"], sep="\t", index_col=0).equals(feature_data)
    assert len(recording_runner.calls) == 1
    assert str(input_bed) in recording_runner.calls[0][0]


def _cupcake_available():
    if shutil.which("conda") is None:
        return False
    r = subprocess.run(["conda", "run", "-n", "iflas-cupcake",
                        "collapse_isoforms_by_sam.py", "--help"], capture_output=True)
    return r.returncode == 0


_TOOLS_FULL = all(shutil.which(t) for t in
                  ("minimap2", "hisat2", "hisat2-build", "stringtie", "regtools",
                   "samtools", "seqkit", "perl"))


@pytest.mark.skipif(not (_TOOLS_BASIC and _TOOLS_FULL and _cupcake_available()),
                    reason="needs cupcake env + minimap2/hisat2/stringtie/regtools/seqkit/perl")
def test_full_chain_filter_hq_isoforms_outputs_hq_bed(tmp_path):
    """Real cupcake collapse -> refine -> NGS reassembly -> filter_hq_isoforms.

    Verifies the tool-heavy feature pipeline end-to-end on real data: getInputBed (anno/
    novel split against the subset ref) + getFeatures (CPC2 ORF, canonical-junction ratio,
    GC) + predCodingP (CPC2 LIBSVM coding potential) + addIsoSupport + addJuncIndelInfo
    (pysam over the NGS tmp.bam + junction RPKM from junctions.bed) + PU scoring +
    utils/filter.pl -> hqIsoforms/hq.collapsed.bed12+.
    """
    from iflas.stages.align import minimap2_tgs, ngs_reassembly_junctions
    from iflas.stages.collapse import collapse_isoforms
    from iflas.stages.find_as import find_as_sample
    from iflas.stages.pu_filter import filter_hq_isoforms

    d = str(tmp_path)
    syn.generate(d)
    base = os.path.join(d, "sample")
    os.makedirs(base, exist_ok=True)

    # refine uses the FULL annotation; the PU anno/novel split uses a SUBSET ref
    full_gpe = os.path.join(d, "ref.full.gpe")
    full_bed = os.path.join(d, "ref.full.bed12+")
    _build_refs(os.path.join(d, "genes.gtf"), full_gpe, full_bed)
    subset_gtf = os.path.join(d, "ref.subset.gtf")
    _write_subset_gtf(subset_gtf)
    subset_bed = os.path.join(d, "ref.subset.bed12+")
    _build_refs(subset_gtf, os.path.join(d, "ref.subset.gpe"), subset_bed)

    # TGS -> cupcake collapse -> refine (into base/refine + base/collapse group)
    bam = os.path.join(d, "flnc.bam")
    minimap2_tgs(f"{d}/genome.fa", f"{d}/flnc.fa", bam, platform="pacbio", threads=2)
    tofu = os.path.join(base, "collapse", "tofu")
    os.makedirs(os.path.dirname(tofu), exist_ok=True)
    collapse_isoforms(f"{d}/flnc.fa", bam, tofu, fl_coverage=2)
    refine_isoforms(tofu + ".collapsed.good.gff", f"{d}/genome.fa", full_gpe, full_bed, base)

    # NGS reassembly -> base/mapping/rna-seq/reassembly/{tmp.bam,junctions.bed}
    # (ngs_reassembly_junctions appends its own reassembly/ subdir, so pass .../rna-seq)
    idx = os.path.join(d, "hisat2_idx")
    subprocess.run(["hisat2-build", "-q", f"{d}/genome.fa", idx], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ngs_reassembly_junctions(
        idx, os.path.join(base, "mapping", "rna-seq"),
        repeats=[{"left": f"{d}/ngs_R1.fastq.gz", "right": f"{d}/ngs_R2.fastq.gz"}],
        threads=2, filter_count=0)
    reasm = os.path.join(base, "mapping", "rna-seq", "reassembly")

    out = filter_hq_isoforms(
        base, f"{d}/genome.fa", subset_bed,
        junction_file=os.path.join(reasm, "junctions.bed"),
        sam_file=os.path.join(reasm, "tmp.bam"),
        filter_score=0.0)

    expected_annotated, expected_novel = _expected_getinputbed_counts()
    fd = pd.read_csv(out["features"], sep="\t", index_col=0)
    assert (fd.annotation == "annotated").sum() == expected_annotated
    assert (fd.annotation == "novel").sum() == expected_novel
    # NGS junction features wired through (RPKM from junctions.bed)
    assert (fd.minJuncRPKM > 0).all()
    # CPC2 features present (ORF finder + LIBSVM coding potential)
    for col in ("orfLength", "FickettScore", "codingP", "pepLength", "GC"):
        assert col in fd.columns
    score_rows = _rows(out["pu_score"])
    valid_rows = _rows(out["valid"])
    hq_rows = _rows(out["hq_bed"])
    with open(os.path.join(base, "collapse", "tofu.collapsed.group.txt")) as f:
        group_keys = {line.split("\t")[0] for line in f if line.strip()}
    assert score_rows[0] == ["", "label", "pu_score"]
    assert valid_rows[0] == ["isoform", "label", "pu_score"]
    assert len(score_rows) == expected_annotated + expected_novel + 1
    assert len(hq_rows) == expected_annotated + expected_novel
    assert all(len(r) == 19 for r in hq_rows)
    assert {r[3] for r in hq_rows} == set(fd.index)
    # find_as_lr reads Bed12Plus.otherList[2] (1-based col15) as the group key.
    assert all(set(r[14].split(",")).issubset(group_keys) for r in hq_rows)

    as_out = find_as_sample(base, full_gpe)
    assert "geneB:2200@2500-2650@3000" in [r[3] for r in _rows(as_out["SE"])]
    assert "geneC:4000-4050" in [r[3] for r in _rows(as_out["A5SS"])]
    assert "geneE:4950-5000" in [r[3] for r in _rows(as_out["A3SS"])]


def test_filter_hq_isoforms_accepts_relative_paths(tmp_path, monkeypatch):
    """BUG-147: the stage chdirs into hqIsoforms/, so relative CLI paths must be absolutized first."""
    from iflas.stages import pu_filter

    base = tmp_path / "base"
    base.mkdir()
    (tmp_path / "genome.fa").write_text(">chr1\nACGT\n")
    (tmp_path / "ref.bed12+").write_text("")
    (tmp_path / "input.bed12+").write_text("")
    pd.DataFrame(
        {"isoLength": [1000, 1200], "flCount": [3, 2], "ratioIsoToGene": [0.8, 0.4], "exonNum": [2, 2],
         "annotation": ["annotated", "novel"], "GC": [50.0, 48.0], "orfLength": [600, 480], "orfIntegrity": [1, 1],
         "pepLength": [200, 160], "FickettScore": [0.5, 0.4], "pI": [6.0, 6.2], "codingP": [0.9, 0.3],
         "codingLabel": ["coding", "noncoding"], "canJuncRatio": [1.0, 1.0], "sdJuncCov": [0.1, 0.2],
         "withNovelJunc": [False, True], "minNovelJuncRPKM": [0.0, 0.7], "bite": [False, False],
         "minJuncRPKM": [1.0, 0.8], "nIndelsAroundJunc": [0, 0], "ratioMinJuncCovToAllCov": [0.7, 0.6],
         "nJuncsWithIndels": [0, 0], "indelNearJunc": [False, False]},
        index=["anno1", "novel1"],
    ).to_csv(tmp_path / "features.tsv", sep="\t", index=True, index_label="isoform")
    seen = {}

    def fake_get_input_bed(baseDir, refBed):
        seen["base"], seen["ref_bed"] = baseDir, refBed
        return str(tmp_path / "input.bed12+")

    class DummySchemeModelEval:
        def __init__(self, trainingData, model, scheme, random_state=0):
            self.predResults = pd.DataFrame({"label": [1, 0], "pu_score": [1.0, 0.75]}, index=["anno1", "novel1"])

        def eval(self):
            pass

        def filterIsoformsByScore(self, filterScore, outFile, lqAnnoIso):
            with open(outFile, "w") as f:
                f.write("isoform\tlabel\tpu_score\nnovel1\t0\t0.75\n")

    class Runner:
        def run(self, cmd, **kwargs):
            seen["cmd"] = cmd
            with open("hq.collapsed.bed12+", "w") as f:
                f.write("hq\n")

    monkeypatch.setattr(pu_filter, "getInputBed", fake_get_input_bed)
    monkeypatch.setattr(pu_filter, "SchemeModelEval", DummySchemeModelEval)
    monkeypatch.chdir(tmp_path)
    out = pu_filter.filter_hq_isoforms("base", "genome.fa", "ref.bed12+", feature_file="features.tsv",
                                       sam_file="tmp.bam", junction_file="junctions.bed", runner=Runner())
    # every path handed to the core / the runner is absolute although the caller passed relative ones
    assert os.path.isabs(seen["base"]) and os.path.isabs(seen["ref_bed"])
    assert os.path.isabs(out["hq_bed"]) and out["hq_bed"].startswith(str(base))
    assert os.getcwd() == str(tmp_path)


def test_filter_hq_isoforms_can_select_the_best_model_and_threshold(tmp_path, monkeypatch):
    """F-10: select_best_model / auto_filter_score reach the stage (legacy iso_pu1 honoured them, the modern stage
    ignored the config fields): the chosen classifier is trained and the threshold is the mean of the per-fold
    95th-percentile negatives' scores."""
    from iflas.stages import pu_filter

    base = tmp_path / "base"
    base.mkdir()
    genome, ref_bed, input_bed = tmp_path / "genome.fa", tmp_path / "ref.bed12+", tmp_path / "input.bed12+"
    genome.write_text(">chr1\nACGT\n")
    ref_bed.write_text("")
    input_bed.write_text("")
    feature_data = pd.DataFrame(
        {"isoLength": [1000, 1200], "flCount": [3, 2], "ratioIsoToGene": [0.8, 0.01], "exonNum": [2, 2], "GC": [50.0, 48.0],
         "orfLength": [600, 480], "orfIntegrity": [1, 1], "pepLength": [200, 160], "FickettScore": [0.5, 0.4], "pI": [6.0, 6.2],
         "codingP": [0.9, 0.3], "canJuncRatio": [1.0, 1.0], "sdJuncCov": [0.1, 0.2], "minJuncRPKM": [1.0, 0.8],
         "nIndelsAroundJunc": [0, 0], "ratioMinJuncCovToAllCov": [0.7, 0.6], "nJuncsWithIndels": [0, 0],
         "indelNearJunc": [False, False], "annotation": ["annotated", "novel"]}, index=["anno1", "novel1"])
    seen = {}

    class DummySchemeModelEval:
        def __init__(self, trainingData, model, scheme, random_state=0):
            seen["model"] = model
            self.predResults = pd.DataFrame({"label": [1, 0], "pu_score": [1.0, 0.75]}, index=["anno1", "novel1"])

        def eval(self):
            pass

        def filterIsoformsByScore(self, filterScore, outFile, lqAnnoIso):
            seen["filter_score"] = filterScore
            with open(outFile, "w") as f:
                f.write("isoform\tlabel\tpu_score\nnovel1\t0\t0.75\n")

    def fake_select(featureData, drawAUC=False, selectPUscore=False, random_state=0):
        seen["draw_auc"], seen["select_pu_score"], seen["seed"] = drawAUC, selectPUscore, random_state
        return "RF", {"RF": [0.6, 0.7], "GB": [0.9]}

    class Runner:
        def run(self, cmd, **kwargs):
            open("hq.collapsed.bed12+", "w").write("hq\n")

    monkeypatch.setattr(pu_filter, "getInputBed", lambda baseDir, refBed: str(input_bed))
    monkeypatch.setattr(pu_filter, "processInputFile", lambda *a, **k: feature_data)
    monkeypatch.setattr(pu_filter, "SchemeModelEval", DummySchemeModelEval)
    monkeypatch.setattr(pu_filter, "selectBestModel", fake_select)
    monkeypatch.setattr(pu_filter, "_find_utils_dir", lambda: "/repo/utils")
    monkeypatch.setattr(pu_filter, "MIN_SELECTION_SET", 1)         # the two-row fixture is enough here

    out = pu_filter.filter_hq_isoforms(str(base), str(genome), str(ref_bed), runner=Runner(), select_best_model=True,
                                       auto_filter_score=True, draw_auc=True)

    assert seen == {"model": "RF", "filter_score": pytest.approx(0.65), "draw_auc": True, "select_pu_score": True, "seed": 0}
    assert out["model"] == "RF" and out["filter_score"] == pytest.approx(0.65)

    seen.clear()
    out = pu_filter.filter_hq_isoforms(str(base), str(genome), str(ref_bed), runner=Runner())
    assert seen == {"model": "GB", "filter_score": 0.5} and out["model"] == "GB"      # default path unchanged

    seen.clear()
    monkeypatch.setattr(pu_filter, "MIN_SELECTION_SET", 5)
    out = pu_filter.filter_hq_isoforms(str(base), str(genome), str(ref_bed), runner=Runner(), select_best_model=True)
    assert seen == {"model": "GB", "filter_score": 0.5} and out["model"] == "GB"      # too few negatives: fallback, no selection


def test_write_isoform_score_derives_the_novel_hq_as_report_input(tmp_path):
    """F-9: the six columns `report --novel-hq-as` plots come from the PU feature table."""
    from iflas.stages.pu_filter import write_isoform_score
    feats = pd.DataFrame({"flCount": [8, 2, 5], "ratioIsoToGene": [0.8, 0.2, 1.0], "annotation": ["annotated", "novel", "annotated"]},
                         index=["PB.1.1", "PB.1.2", "PB.2.1"])
    out = write_isoform_score(feats, str(tmp_path / "isoIdentity" / "isoformScore.txt"))
    assert [l.split("\t") for l in open(out).read().splitlines()] == [
        ["PB.1", "PB.1.1", "8", "10", "0.8", "Annotated"], ["PB.1", "PB.1.2", "2", "10", "0.2", "Novel"], ["PB.2", "PB.2.1", "5", "5", "1.0", "Annotated"]]


def test_pu_filter_without_ngs_says_what_it_needs(tmp_path):
    """Long-read-only samples: the PU features include NGS junction support (published design); feature extraction stops
    with what to provide instead of a KeyError on the missing feature columns (BUG-191)."""
    import pytest

    from iflas.pu.iso_pu import processInputFile

    with pytest.raises(ValueError, match="needs the sample's NGS alignments"):
        processInputFile(str(tmp_path / "in.bed"), str(tmp_path / "genome.fa"), str(tmp_path / "ref.bed"))
