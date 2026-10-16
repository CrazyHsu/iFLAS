"""Review Round 3 (2026-10-14): robustness guards found while reading the Phase D drivers. Each test reproduces the
legacy-inherited crash on degenerate inputs and asserts the stage now finishes with a notice instead."""

from __future__ import annotations

import os
from types import SimpleNamespace

import pytest

from iflas.stages.asas import ase_sample
from iflas.stages.characterize import splicesite2seq
from iflas.stages.diff_as import BamRecord, run_rmats
from iflas.stages.find_pa import find_pa_sample
from iflas.stages.palen_as import getPalenAS


class _FakeRunner:
    def __init__(self):
        self.calls = []

    def run(self, cmd, **kwargs):
        self.calls.append(cmd)
        return SimpleNamespace(returncode=0, stdout="", stderr="")


def _bed12plus(chrom, start, end, name, strand, blocks, extra):
    sizes = ",".join(str(e - s) for s, e in blocks)
    starts = ",".join(str(s - start) for s, _ in blocks)
    return "\t".join(map(str, [chrom, start, end, name, 0, strand, start, end, "0,0,0", len(blocks), sizes, starts] + extra))


def test_find_pa_without_any_pa_cluster_skips_the_motif_analysis(tmp_path):
    """One read per gene can never reach the 5-read PA cluster support: PA.bed6+ is empty and the motif tables
    would be empty files that pandas refused to read (legacy crashed here too)."""
    base = tmp_path / "s1"
    (base / "collapse").mkdir(parents=True)
    (base / "refine").mkdir()
    reads = [_bed12plus("chr1", 100 + i * 1000, 600 + i * 1000, f"read{i}", "+", [(100 + i * 1000, 600 + i * 1000)],
                        [f"PB.{i}.1", "annotated", f"gene{i}"]) for i in range(3)]
    (base / "refine" / "reads.assigned.unambi.bed12+").write_text("".join(r + "\n" for r in reads))
    (base / "refine" / "isoformGrouped.bed12+").write_text("".join(r + "\n" for r in reads))
    (base / "collapse" / "tofu.collapsed.group.txt").write_text("".join(f"PB.{i}.1\tread{i}\n" for i in range(3)))
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "A" * 5000 + "\n")
    (tmp_path / "genome.len.txt").write_text("chr1\t5000\n")
    runner = _FakeRunner()

    out = find_pa_sample(str(base), str(genome), str(tmp_path / "genome.len.txt"), filter_pa_by_count=5, runner=runner)

    pa_dir = base / "as_events" / "pa"
    assert out["pa_dir"] == str(pa_dir)
    assert (pa_dir / "PA.bed6+").stat().st_size == 0
    assert not (pa_dir / "motif").exists()                      # skipped, not crashed
    assert not any("lines.R" in str(c) for c in runner.calls)


def test_splicesite2seq_counts_only_junctions_with_both_sides_in_the_genome(tmp_path):
    """A junction whose acceptor window lies beyond the chromosome end made getfasta warn on stderr, which
    pybedtools raises as BEDToolsError (and the donor/acceptor pairing would then KeyError); such windows are
    dropped and the junction is left out of the count."""
    seq = list("A" * 100)
    seq[20:22] = "GT"
    seq[48:50] = "AG"
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "".join(seq) + "\n")
    sites = tmp_path / "IR.splicesite"
    sites.write_text("chr1\t20\t50\t+\nchr1\t30\t101\t+\n")
    out = tmp_path / "IR.tmp"

    splicesite2seq(str(genome), str(sites), outFile=str(out))

    assert out.read_text().splitlines() == ["GT-AG\t1"]


def test_run_rmats_rejects_a_missing_ngs_read_length(tmp_path):
    gtf = tmp_path / "merged.gtf"
    gtf.write_text("")
    records = [BamRecord("p", "a", "condA", "repeat0", "/x/a.bam", "paired", 0),
               BamRecord("p", "b", "condB", "repeat0", "/x/b.bam", "paired", 0)]
    runner = _FakeRunner()
    with pytest.raises(ValueError, match="read length"):
        run_rmats(records, str(gtf), tmp_path / "das", runner=runner)
    assert runner.calls == []


def test_ase_sample_without_a_representative_isoform_per_haplotype_skips_quantification(tmp_path):
    """One phased SNP but no tagged NGS read and no haplotype FASTA: salmon index on an empty FASTA used to fail."""
    base = tmp_path / "s1"
    locus = base / "asas" / "by_loci" / "PB.1_size5"
    locus.mkdir(parents=True)
    (base / "refine").mkdir()
    (base / "refine" / "isoformGrouped.bed12+").write_text(
        _bed12plus("chr1", 100, 900, "PB.1.1", "+", [(100, 300), (500, 900)], ["gene1", "annotated", "PB.1.1"]) + "\n")
    (locus / "phased.partial.cleaned.human_readable.txt").write_text("# IsoPhase\nhaplotype\tPB.1.1\nA\t5\nC\t0\n")
    (locus / "phased.partial.cleaned.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t150\t.\tA\tC\t.\tPASS\t.\n")
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "A" * 1000 + "\n")
    runner = _FakeRunner()

    out = ase_sample(str(base), str(genome), ngs_left="r1.fq", ngs_right="r2.fq", runner=runner)

    assert out["merged"] == []
    assert os.path.getsize(out["snp_positions"]) > 0            # the SNP was collected ...
    assert not any("salmon" in str(c) for c in runner.calls)    # ... but nothing was quantified


def test_getPalenAS_skips_events_whose_isoforms_are_absent_from_the_isoform_bed(tmp_path, monkeypatch):
    """The representative isoform is picked from the isoform bed by the member ids of column 15; when none of them
    is a record there the empty list was indexed (IndexError)."""
    iso_bed = tmp_path / "hq.bed12+"
    iso_bed.write_text(_bed12plus("chr1", 100, 900, "PB.1.1", "+", [(100, 300), (500, 900)], ["gene1", "annotated", "PB.1.1a"]) + "\n"
                       + _bed12plus("chr1", 100, 900, "PB.1.2", "+", [(100, 900)], ["gene1", "annotated", "PB.1.2a"]) + "\n")
    reads_bed = tmp_path / "reads.bed12+"
    reads_bed.write_text(_bed12plus("chr1", 100, 900, "r0", "+", [(100, 900)], ["x", "y", "z"]) + "\n")
    inc, exc = [f"r{i}" for i in range(10)], [f"s{i}" for i in range(10)]
    palen = {**{r: 20.0 + i for i, r in enumerate(inc)}, **{s: 80.0 + i for i, s in enumerate(exc)}}
    as_pairs = {"SE": {"PB.1": [[["PB.1.1"], ["PB.1.2"], "PB.1:chr1:300-500"]]}}
    monkeypatch.chdir(tmp_path)

    getPalenAS(palen, str(iso_bed), str(reads_bed), collapsedTrans2reads={"PB.1.1a": inc, "PB.1.2a": exc},
               asPairs=as_pairs, filterByCount=10, mergeByJunc=True)

    assert (tmp_path / "mergeByJunc" / "SE.palenAndAS.sig.bed12+").stat().st_size == 0
    assert os.getcwd() == str(tmp_path)


def test_ase_haplotype_without_ngs_reads_gets_zero_counts(tmp_path):
    """BUG-161 (found on the 1 % scale data): no NGS read of a replicate carried the ALT allele, the ALT read file
    was empty and salmon refused it. The helpers mark the haplotype unexpressed instead, and the merged ASE table
    reports 0 NGS support for it."""
    import gzip
    from iflas.stages.asas import _is_empty_fastq, mergeASEfromSalmonOut, write_zero_quant
    empty, full = tmp_path / "alt.fq.gz", tmp_path / "ref.fq.gz"
    with gzip.open(empty, "wb"):
        pass
    with gzip.open(full, "wb") as f:
        f.write(b"@r1\nACGT\n+\nIIII\n")
    assert _is_empty_fastq(str(empty)) and not _is_empty_fastq(str(full)) and not _is_empty_fastq(str(tmp_path / "missing.fq.gz"))
    (tmp_path / "alt.fa").write_text(">PB.1_size5_ALT\nACGTACGT\n")
    alt_quant = write_zero_quant(str(tmp_path / "alt.fa"), str(tmp_path / "repeat0.alt_salmon_quant"))
    (tmp_path / "ref_quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nPB.1_size5_REF\t8\t8\t1000000.0\t12.0\n")
    out = tmp_path / "merged.txt"
    mergeASEfromSalmonOut(str(tmp_path / "ref_quant.sf"), alt_quant, {"PB.1_size5_REF": [9, "g1"], "PB.1_size5_ALT": [3, "g1"]}, str(out))
    assert out.read_text().split("\t") == ["PB.1", "g1", "9.0", "3.0", "1000000.0", "12.0", "0.0", "0.0\n"]
