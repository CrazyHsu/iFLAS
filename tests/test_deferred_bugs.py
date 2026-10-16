"""Regression tests of the fixes that change results compared with iFLAS 1 (listed in CHANGELOG.md).

Each test reproduces an iFLAS 1 defect on a minimal input and asserts the corrected behaviour.

Fixtures are tiny hand-written BED12+/genePredExt/FASTA/SAM files. bedtools (pybedtools),
gtfToGenePred and samtools (pysam) live in the core ``iflas`` env.
"""
from __future__ import annotations

import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest


def deferred(bug):
    return pytest.mark.xfail(
        strict=True, raises=AssertionError,
        reason=f"{bug} deferred (class b): awaiting user sign-off after the synthetic-v2 impact report")


needs_bedtools = pytest.mark.skipif(shutil.which("bedtools") is None,
                                    reason="needs bedtools (core iflas env)")
needs_refine_tools = pytest.mark.skipif(
    shutil.which("bedtools") is None or shutil.which("gtfToGenePred") is None,
    reason="needs bedtools + gtfToGenePred (core iflas env)")


# --------------------------------------------------------------------------- fixtures
def _bed12plus(chrom, name, strand, exons, extra):
    """0-based half-open exons -> 12 BED columns + extra columns (one line)."""
    start, end = exons[0][0], exons[-1][1]
    sizes = ",".join(str(e - s) for s, e in exons)
    starts = ",".join(str(s - start) for s, e in exons)
    return "\t".join(map(str, [chrom, start, end, name, 0, strand, start, end, "0,0,0",
                               len(exons), sizes, starts] + list(extra))) + "\n"


def _gpe(name, chrom, strand, exons, gene):
    """genePredExt line without bin column (what init-ref / gtfToGenePred -genePredExt write)."""
    starts = ",".join(str(s) for s, e in exons) + ","
    ends = ",".join(str(e) for s, e in exons) + ","
    frames = ",".join("0" for _ in exons) + ","
    return "\t".join(map(str, [name, chrom, strand, exons[0][0], exons[-1][1], exons[0][0],
                               exons[-1][1], len(exons), starts, ends, 0, gene, "cmpl", "cmpl",
                               frames])) + "\n"


def _junction_bed12(chrom, donor, acceptor, strand, score):
    """regtools-style BED12 junction whose single intron is (donor, acceptor), 0-based half-open."""
    cs, ce = donor - 1, acceptor + 1
    return "\t".join(map(str, [chrom, cs, ce, f"junc_{donor}_{acceptor}", score, strand, cs, ce,
                               "0,0,0", 2, "1,1", f"0,{acceptor - cs}"])) + "\n"


def _fasta(path, chrom, seq, width=60):
    with open(path, "w") as f:
        f.write(f">{chrom}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


def _rows(path):
    with open(path) as f:
        return [l.rstrip("\n").split("\t") for l in f if l.strip()]


# --------------------------------------------------------------------------- BUG-133
def test_bug133_canonical_splice_site_check_is_case_insensitive():
    """Soft-masked genomes (e.g. Ensembl *dna_sm*) hold repeats in lowercase; 'gt' is a GT donor."""
    from iflas.stages.refine import isCanonicalSite
    assert isCanonicalSite("+", "donor", "gt") == isCanonicalSite("+", "donor", "GT") == 1
    assert isCanonicalSite("-", "accept", "ct") == isCanonicalSite("-", "accept", "CT") == 1


@needs_bedtools
def test_bug133_canonical_junction_fraction_on_soft_masked_genome(tmp_path):
    from iflas.io.records import BedFile
    from iflas.pu.iso_pu import GetFeatures
    seq = list("a" * 3000)
    seq[1000:1002], seq[1998:2000] = "gt", "ag"          # canonical intron 1000-2000, lowercase
    genome = tmp_path / "genome.fa"
    _fasta(str(genome), "chr1", "".join(seq))
    bed = tmp_path / "iso.bed12+"
    bed.write_text(_bed12plus("chr1", "iso1", "+", [(500, 1000), (2000, 2500)], ["PB.1", "PB.1", "iso1"]))
    chains = BedFile(str(bed), type="bed12+").getJuncChainDict()
    freq = GetFeatures(inputBed=str(bed), genomeFa=str(genome)).getCanonicalJuncFreq(chains, str(genome))
    assert freq["chr1:1000-2000"] == 1.0            # 0.0 today: 'gt-ag' != 'GT-AG'


# --------------------------------------------------------------------------- BUG-134
def test_bug134_consensus_intron_tolerance_is_symmetric():
    """readsAssign's getConsensusIntronN (offset=10) must accept a reference intron shifted either way."""
    from iflas.intervalops import getConsensusIntronN
    read_starts, read_ends = [0, 2000], [1000, 3000]            # read intron 1000-2000
    assert getConsensusIntronN(read_starts, read_ends, [0, 1997], [997, 3000], 10) == 1   # ref shifted -3
    assert getConsensusIntronN(read_starts, read_ends, [0, 2003], [1003, 3000], 10) == 1   # ref shifted +3 (0 today)


# --------------------------------------------------------------------------- BUG-135
@needs_bedtools
def test_bug135_junc_chain_flank_match_respects_junction_boundaries(tmp_path, monkeypatch):
    from iflas.stages.refine import readsAssign
    trans, reads = tmp_path / "trans.bed12+", tmp_path / "reads.bed12+"
    trans.write_text(_bed12plus("chr1", "tx1", "+", [(0, 100), (200, 2300), (2400, 2500)], ["geneX"]))  # 100-200;2300-2400
    reads.write_text(_bed12plus("chr1", "readA", "+", [(0, 300), (2400, 2500)], ["PB.1"]))              # 300-2400
    monkeypatch.chdir(tmp_path)
    readsAssign(str(trans), str(reads), readsColNum=13, outPrefix="ra", group=True)
    row = {r[3]: r for r in _rows("ra.bed12+")}["readA"]
    # no shared intron, exon coverage 300/400 < 0.9, and "300-2400" is NOT a flank of
    # "100-200;2300-2400" -> readType must be "I" (today "E": the regex "300-2400$" matches ";2300-2400")
    assert row[13] == "I"


# --------------------------------------------------------------------------- BUG-136
@needs_bedtools
def test_bug136_novel_cluster_ids_do_not_leak_across_chromosomes(tmp_path, monkeypatch):
    from iflas.stages.refine import readsAssign
    trans, reads = tmp_path / "trans.bed12+", tmp_path / "reads.bed12+"
    trans.write_text(_bed12plus("chr9", "tx1", "+", [(100000, 101000)], ["geneX"]))      # overlaps nothing
    reads.write_text(_bed12plus("chr1", "a", "+", [(100, 500)], ["PB.1"])
                     + _bed12plus("chr2", "b", "+", [(200, 9000)], ["PB.2"])
                     + _bed12plus("chr1", "c", "+", [(600, 700)], ["PB.3"]))
    monkeypatch.chdir(tmp_path)
    readsAssign(str(trans), str(reads), readsColNum=13, outPrefix="ra", group=True)
    ids = {r[3]: r[13] for r in _rows("ra.unambi.bed12+")}
    assert ids["a"].startswith("chr1:+:") and ids["c"].startswith("chr1:+:") and ids["b"].startswith("chr2:+:")
    # a (100-500) and c (600-700) do not overlap -> different novel clusters
    # (today both are chr1:+:1 because chr2's 200-9000 read inflated clusterEnd)
    assert ids["a"] != ids["c"]


# --------------------------------------------------------------------------- BUG-138
@needs_bedtools
def test_bug138_single_exon_strand_evidence_survives_a_low_coverage_transcript(tmp_path, monkeypatch):
    from iflas.stages.refine import strandAdjust
    seq = list("A" * 4000)
    seq[2200:2202], seq[2398:2400] = "GT", "AG"             # canonical intron for the multi-exon read
    genome = tmp_path / "genome.fa"
    _fasta(str(genome), "chr1", "".join(seq))
    gpe = tmp_path / "ref.gpe"
    # txA fully covers r1 (coverage 1.0, '+'); txB barely touches it (coverage ~0.05). bedtools reports
    # hits in reference order, so txB is processed after txA.
    gpe.write_text(_gpe("txA", "chr1", "+", [(100, 1000)], "geneA")
                   + _gpe("txB", "chr1", "+", [(950, 2900)], "geneB"))
    bed = tmp_path / "reads.bed12+"
    bed.write_text(_bed12plus("chr1", "r1", "-", [(100, 1000)], ["PB.1"])                 # mis-stranded mono-exon
                   + _bed12plus("chr1", "r2", "+", [(2000, 2200), (2400, 2600)], ["PB.2"]))
    monkeypatch.chdir(tmp_path)
    strandAdjust(str(genome), str(gpe), str(bed), 0.5, 0.5,
                 strandAdjust="adj.bed12+", strandConfirmed="conf.bed12+")
    row = {r[3]: r for r in _rows("adj.bed12+")}["r1"]
    # txA's '+' evidence must survive; today the txB else-branch overwrites it with r1's own '-'
    assert row[5] == "+"


# --------------------------------------------------------------------------- BUG-139
@needs_bedtools
def test_bug139_refine_with_junc_orders_covering_junctions_numerically(tmp_path, monkeypatch):
    from iflas.stages.refine import refineWithJunc
    seq = list("A" * 3000)
    for pos, dinuc in ((900, "GT"), (998, "AG"), (1050, "GT"), (1498, "AG")):
        seq[pos:pos + 2] = dinuc
    genome = tmp_path / "genome.fa"
    _fasta(str(genome), "chr1", "".join(seq))
    bed = tmp_path / "reads.bed12+"
    bed.write_text(_bed12plus("chr1", "R", "+", [(0, 900), (1500, 2000)], ["geneX"]))   # juncChain 900-1500
    junc = tmp_path / "junctions.bed"
    # NGS says 900-1500 is really two introns flanking a 1000-1050 micro-exon (50 bp: inside the +/-5 % exon-length
    # band refineWithJunc accepts for a refined read)
    junc.write_text(_junction_bed12("chr1", 900, 1000, "+", 10) + _junction_bed12("chr1", 1050, 1500, "+", 10))
    monkeypatch.chdir(tmp_path)
    refineWithJunc(bedFile=str(bed), junctionFile=str(junc), outFile="refined.bed12+", genomeFasta=str(genome))
    rows = {r[3]: r for r in _rows("refined.bed12+")}
    # sorting "1050-1500" before "900-1000" *as strings* built an invalid chain and DROPPED the read (BUG-139)
    assert "R" in rows
    assert rows["R"][9] == "3"


# --------------------------------------------------------------------------- BUG-140
def test_bug140_indel_support_ratio_uses_true_division(tmp_path):
    import pysam
    from iflas.io.records import BedFile
    from iflas.pu.iso_pu import GetFeatures
    iso = tmp_path / "input.bed12+"
    iso.write_text(_bed12plus("chr1", "iso1", "+", [(800, 1000), (1200, 1400)],
                              ["PB.1", "PB.1", "iso1", "3", "3", "1.0", "novel"]))
    ref = tmp_path / "ref.bed12+"
    ref.write_text(_bed12plus("chr1", "tx1", "+", [(800, 1000), (1200, 1400)], ["geneX"]))
    junc = tmp_path / "junctions.bed"
    junc.write_text(_junction_bed12("chr1", 1000, 1200, "+", 4))              # junction coverage 4
    sam = tmp_path / "reads.sam"
    # three NGS reads, each with a 2-bp deletion right before the splice donor (1-based 999-1000)
    body = "".join(f"r{i}\t0\tchr1\t985\t60\t14M2D200N20M\t*\t0\t0\t{'A' * 34}\t{'I' * 34}\n" for i in range(3))
    sam.write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:5000\n" + body)
    bam = str(tmp_path / "reads.bam")
    pysam.sort("-o", bam, str(sam))
    pysam.index(bam)
    gf = GetFeatures(inputBed=str(iso))
    gf.inputBedObj = BedFile(str(iso), type="bed12+")
    gf.features = pd.DataFrame(index=["iso1"])
    feats = gf.addJuncIndelInfo(samFile=bam, junctionFile=str(junc), refBed=str(ref))
    # 3 indel-bearing reads / 4 junction reads = 0.75 >= 0.5 -> the indel is "valid" and flags the isoform.
    # Today int//int floors 3//4 to 0, so the 0.5 threshold silently means "indel reads >= junction reads".
    assert int(feats.loc["iso1", "indelNearJunc"]) == 1


# --------------------------------------------------------------------------- BUG-141
def test_bug141_pu_scores_are_reproducible_across_runs():
    from iflas.pu.iso_pu import SchemeModelEval
    rng = np.random.RandomState(0)
    X = np.vstack([rng.normal(1, 1, (12, 5)), rng.normal(0, 1, (12, 5))])
    df = pd.DataFrame(X, columns=list("abcde"))
    df["label"] = [1] * 12 + [0] * 12
    runs = []
    for _ in range(2):
        sme = SchemeModelEval(trainingData=df, model="GB", scheme="bagging", bagging_n_jobs=1)
        sme.eval()
        runs.append(sme.predResults.loc[sme.predResults.label == 0, "pu_score"].to_numpy())
    # identical input must give identical PU scores (today BaggingClassifierPU has no random_state)
    assert np.allclose(runs[0], runs[1])


# --------------------------------------------------------------------------- BUG-142
@needs_refine_tools
def test_bug142_refine_driver_passes_strand_adjust_thresholds_in_signature_order(tmp_path):
    """refine_isoforms calls strandAdjust(..., 0.8, 2) but the signature is (juncDiffScore, minCoverage):
    minCoverage=2 can never be reached by a coverage fraction, so reference-based strand correction of
    mono-exon isoforms is dead code (and juncDiffScore=0.8 lets one canonical site decide a strand)."""
    from iflas.io.genepred import GenePredObj
    from iflas.stages.refine import refine_isoforms
    d = tmp_path
    seq = list("A" * 3000)
    seq[1700:1702], seq[1998:2000] = "GT", "AG"
    _fasta(str(d / "genome.fa"), "chr1", "".join(seq))
    (d / "ref.gtf").write_text('chr1\tsyn\texon\t101\t1000\t.\t+\t.\tgene_id "geneX"; transcript_id "tx1";\n')
    subprocess.run(["gtfToGenePred", "-genePredExt", str(d / "ref.gtf"), str(d / "ref.gpe")], check=True)
    with open(d / "ref.bed12+", "w") as f:
        for line in GenePredObj(str(d / "ref.gpe"), bincolumn=False).toBed(gene=True):
            f.write(str(line) + "\n")
    a = 'gene_id "PB.1"; transcript_id "PB.1.1";'
    b = 'gene_id "PB.2"; transcript_id "PB.2.1";'
    (d / "collapsed.gff").write_text(
        f"chr1\tPacBio\ttranscript\t101\t1000\t.\t-\t.\t{a}\nchr1\tPacBio\texon\t101\t1000\t.\t-\t.\t{a}\n"     # mis-stranded mono-exon
        f"chr1\tPacBio\ttranscript\t1501\t2200\t.\t+\t.\t{b}\nchr1\tPacBio\texon\t1501\t1700\t.\t+\t.\t{b}\n"
        f"chr1\tPacBio\texon\t2001\t2200\t.\t+\t.\t{b}\n")
    refine_isoforms(str(d / "collapsed.gff"), str(d / "genome.fa"), str(d / "ref.gpe"), str(d / "ref.bed12+"), str(d))
    row = {r[3]: r for r in _rows(str(d / "refine" / "tofu.strandAdjusted.bed12+"))}["PB.1.1"]
    # PB.1.1 is fully covered (1.0) by a '+' reference transcript -> must be re-stranded to '+'
    assert row[5] == "+"


# --------------------------------------------------------------------------- BUG-111 / 128 / 137 / 150 (CP3, no ledger xfail before)
@needs_bedtools
def test_bug111_mono_exon_annotation_uses_the_real_overlap(tmp_path, monkeypatch):
    """A mono-exon read is 'annotated' when it covers >= 50 % of a known mono-exon transcript; legacy computed the
    overlap from a hard-coded (10,50)/(20,60) pair (always 30 bp)."""
    from iflas.stages.refine import identifyNovelIsoformsByJunctions
    gpe = tmp_path / "ref.gpe"
    gpe.write_text(_gpe("txA", "chr1", "+", [(100, 200)], "geneA") + _gpe("txB", "chr2", "+", [(100, 2000)], "geneB"))
    bed = tmp_path / "reads.bed12+"
    bed.write_text(_bed12plus("chr1", "rA", "+", [(100, 1000)], ["PB.1"])        # covers txA fully (100/100)
                   + _bed12plus("chr2", "rB", "+", [(100, 1000)], ["PB.2"]))     # covers 900/1900 of txB
    monkeypatch.chdir(tmp_path)
    identifyNovelIsoformsByJunctions(str(gpe), str(bed), anno="anno.bed", novel="novel.bed")
    assert [r[3] for r in _rows("anno.bed")] == ["rA"] and [r[3] for r in _rows("novel.bed")] == ["rB"]


@needs_bedtools
def test_bug128_ambiguous_multi_exon_read_takes_the_strand_of_the_overlapping_reference(tmp_path, monkeypatch):
    """One forward-canonical and one reverse-canonical junction leave the motif vote tied; the overlapping '+'
    reference must then decide (legacy compared a set to the integer 1, which is never true)."""
    from iflas.stages.refine import strandAdjust
    seq = list("A" * 4000)
    seq[1000:1002], seq[1198:1200] = "GT", "AG"              # junction 1: canonical on '+'
    seq[2000:2002], seq[2198:2200] = "CT", "AC"              # junction 2: canonical on '-'
    genome = tmp_path / "genome.fa"
    _fasta(str(genome), "chr1", "".join(seq))
    gpe = tmp_path / "ref.gpe"
    gpe.write_text(_gpe("tx1", "chr1", "+", [(500, 1000), (1200, 2000), (2200, 2500)], "geneX"))
    bed = tmp_path / "reads.bed12+"
    bed.write_text(_bed12plus("chr1", "r1", "-", [(500, 1000), (1200, 2000), (2200, 2500)], ["PB.1"]))
    monkeypatch.chdir(tmp_path)
    strandAdjust(str(genome), str(gpe), str(bed), 2, 0.8, strandAdjust="adj.bed12+", strandConfirmed="conf.bed12+")
    assert {r[3]: r[5] for r in _rows("adj.bed12+")}["r1"] == "+"


@needs_bedtools
def test_bug137_reads_are_assigned_only_to_same_strand_transcripts(tmp_path, monkeypatch):
    from iflas.stages.refine import readsAssign
    trans, reads = tmp_path / "trans.bed12+", tmp_path / "reads.bed12+"
    trans.write_text(_bed12plus("chr1", "tx1", "+", [(0, 1000)], ["geneX"]))
    reads.write_text(_bed12plus("chr1", "readA", "-", [(100, 900)], ["PB.1"]))     # antisense, fully inside tx1
    monkeypatch.chdir(tmp_path)
    readsAssign(str(trans), str(reads), readsColNum=13, outPrefix="ra", group=True)
    assigned = {r[3]: r[13] for r in _rows("ra.unambi.bed12+")}
    assert assigned["readA"].startswith("chr1:-:")                                  # a novel cluster, not geneX


def test_bug150_mono_exon_isoforms_are_merged_by_span_not_by_empty_junction_chain(tmp_path):
    from types import SimpleNamespace
    from iflas.stages.diff_as import merge_isoforms
    base = tmp_path / "s1"
    (base / "refine").mkdir(parents=True)
    lr = base / "as_events" / "ordinary_as" / "LR"
    lr.mkdir(parents=True)
    (base / "refine" / "isoformGrouped.bed12+").write_text(
        _bed12plus("chr1", "PB.1.1", "+", [(100, 900)], ["geneA", "annotated", "PB.1.1"])
        + _bed12plus("chr1", "PB.2.1", "+", [(5000, 5900)], ["geneB", "annotated", "PB.2.1"]))
    for fn, cols in (("A3SS.confident.bed6+", 10), ("A5SS.confident.bed6+", 10), ("SE.confident.bed12+", 18)):
        (lr / fn).write_text("")
    (lr / "IR.confident.bed6+").write_text("\t".join(["chr1", "100", "900", "e", "0", "+", "x", "PB.1.1", "y", "PB.2.1"]) + "\n")
    runner = SimpleNamespace(run=lambda *a, **k: None)
    (tmp_path / "das").mkdir()
    merge_isoforms([("s1", base)], tmp_path / "das", runner=runner)
    assert len(_rows(str(tmp_path / "das" / "all_sample_merged_iso.bed"))) == 2       # legacy: one record (key "chr1_")


# --------------------------------------------------------------------------- BUG-163 (adopted 2026-10-15)

def test_bug163_ngs_alignment_order_does_not_depend_on_thread_scheduling(tmp_path):
    """`hisat2 -p` wrote reads in completion order, so ties in the sorted BAM -- and with them StringTie's assembly and
    regtools' junction names -- differed between runs (legacy too). `--reorder` keeps input order at no measurable cost."""
    from iflas.stages.align import hisat2_ngs

    class Recorder:
        cmds = []

        def run(self, cmd, **kwargs):
            self.cmds.append(cmd if isinstance(cmd, str) else " ".join(cmd))

    hisat2_ngs("idx", str(tmp_path / "x.bam"), left="r1.fq", right="r2.fq", threads=4, runner=Recorder())
    assert "--reorder" in Recorder.cmds[0]


# --------------------------------------------------------------------------- BUG-178 (adopted at CP-ONT)

def test_bug178_ont_reads_are_tail_trimmed_by_default():
    """ONT direct-RNA reads keep part of their basecalled poly(A), which collapse's -c 0.9 holds against them (a third of
    the synthetic ONT reads dropped). options.ont_trim_polya is on by default; false is the published pipeline."""
    from iflas.config import Options

    assert Options().ont_trim_polya is True


# --------------------------------------------------------------------------- BUG-182
@needs_bedtools
def test_bug182_junction_replacement_leaves_other_junctions_alone(tmp_path, monkeypatch):
    """Refining junction 100-200 to the NGS junction 100-195 must not touch 1100-2005 of the same read, whose text
    contains "100-200" (the regex replace turned it into 1100-1955, inside the +/-5 % band, so silently kept)."""
    from iflas.stages.refine import refineWithJunc
    seq = list("A" * 3000)
    seq[100:102], seq[193:195] = "GT", "AG"
    genome = tmp_path / "genome.fa"
    _fasta(str(genome), "chr1", "".join(seq))
    bed = tmp_path / "reads.bed12+"
    bed.write_text(_bed12plus("chr1", "R", "+", [(0, 100), (200, 1100), (2005, 2500)], ["g1"]))
    junc = tmp_path / "junctions.bed"
    junc.write_text(_junction_bed12("chr1", 100, 195, "+", 10))
    monkeypatch.chdir(tmp_path)
    refineWithJunc(bedFile=str(bed), junctionFile=str(junc), outFile="refined.bed12+", genomeFasta=str(genome))
    row = _rows("refined.bed12+")[0]
    assert (row[10], row[11]) == ("100,905,495", "0,195,2005")       # exons (0,100) (195,1100) (2005,2500)


# --------------------------------------------------------------------------- BUG-181 (adopted 2026-10-15)
def test_bug181_single_intron_isoforms_are_grouped_by_their_junction(tmp_path):
    from iflas.stages.as_pairs import as_event_pairs
    iso = tmp_path / "hq.bed12+"
    iso.write_text(_bed12plus("chr1", "PB.1.1", "+", [(100, 300), (500, 900)], ["g", ".", "PB.1.1"])
                   + _bed12plus("chr1", "PB.1.2", "+", [(100, 350), (500, 900)], ["g", ".", "PB.1.2"])
                   + _bed12plus("chr1", "PB.1.3", "+", [(100, 250), (500, 900)], ["g", ".", "PB.1.3"]))
    group = tmp_path / "group.txt"
    group.write_text("PB.1.1\tr1,r2\nPB.1.2\tr3\nPB.1.3\tr4,r5\n")
    events = tmp_path / "A5SS.confident.bed6+"
    events.write_text("chr1\t250\t350\tPB.1:250-350\t500\t+\tPB.1\tPB.1.1,PB.1.2\t2\tPB.1.3\t1\tinc:300-500\texc:250-500\n")
    pairs = [(inc, exc) for _, _, inc, exc in as_event_pairs(str(events), str(group), str(iso), "A5SS", 0, True)]
    # PB.1.1 and PB.1.2 have different introns (300-500, 350-500): two inclusion groups, not one
    assert pairs == [(["PB.1.1"], ["PB.1.3"]), (["PB.1.2"], ["PB.1.3"])]


# ----------------------------------------------------------------------------- BUG-196 (real-reference check)
def _b196_bed(name, strand, exons, gene):
    s, e = exons[0][0], exons[-1][1]
    return "\t".join(map(str, ["chr1", s, e, name, 0, strand, s, e, "0,0,0", len(exons), ",".join(str(b - a) for a, b in exons),
                               ",".join(str(a - s) for a, _ in exons), gene])) + "\n"


def _b196_gpe(name, strand, exons, gene):
    starts, ends = "".join(f"{a}," for a, _ in exons), "".join(f"{b}," for _, b in exons)
    return "\t".join(map(str, [name, "chr1", strand, exons[0][0], exons[-1][1], exons[0][0], exons[-1][1], len(exons), starts, ends,
                               0, gene, "cmpl", "cmpl", "0," * len(exons)])) + "\n"


@needs_bedtools
def test_bug196_the_same_intron_chain_on_the_other_strand_is_another_isoform(tmp_path, monkeypatch):
    """getInputBed grouped isoforms by chromosome + intron chain only: a 2-read antisense artefact carrying the chain of a
    52-read isoform on the other strand absorbed it (the longer isoform represents the group) and the group was labelled
    by the annotation of either strand; refine's anno/novel split ignored the strand as well."""
    from iflas.pu.iso_pu import getInputBed
    from iflas.stages.refine import identifyNovelIsoformsByJunctions

    sense, antisense = [(100, 200), (300, 400)], [(90, 200), (300, 410)]
    (tmp_path / "refine").mkdir()
    (tmp_path / "collapse").mkdir()
    (tmp_path / "ref.bed").write_text(_b196_bed("tP", "+", sense, "gP") + _b196_bed("tM", "-", [(80, 420)], "gM"))
    (tmp_path / "ref.gpe").write_text(_b196_gpe("tP", "+", sense, "gP") + _b196_gpe("tM", "-", [(80, 420)], "gM"))
    isoforms = _b196_bed("PB.1.1", "+", sense, "gP") + _b196_bed("PB.2.1", "-", antisense, "gM")
    (tmp_path / "refine" / "tofu.collapsed.assigned.unambi.bed12+").write_text(isoforms)
    (tmp_path / "isoforms.bed12+").write_text(isoforms)
    (tmp_path / "collapse" / "tofu.collapsed.group.txt").write_text(
        "PB.1.1\t" + ",".join(f"r{i}" for i in range(52)) + "\nPB.2.1\ts1,s2\n")
    monkeypatch.chdir(tmp_path)
    rows = [line.rstrip("\n").split("\t") for line in open(getInputBed(baseDir=str(tmp_path), refBed=str(tmp_path / "ref.bed")))]
    assert sorted((r[3], r[5], r[-5], r[-4], r[-1]) for r in rows) == [("PB.1.1", "+", "PB.1.1", "52", "annotated"),
                                                                         ("PB.2.1", "-", "PB.2.1", "2", "novel")]
    identifyNovelIsoformsByJunctions(str(tmp_path / "ref.gpe"), str(tmp_path / "isoforms.bed12+"), anno="anno.bed", novel="novel.bed")
    assert [line.split("\t")[3] for line in open("anno.bed")] == ["PB.1.1"]
    assert [line.split("\t")[3] for line in open("novel.bed")] == ["PB.2.1"]
