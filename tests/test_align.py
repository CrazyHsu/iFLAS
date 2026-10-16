"""Tests for iflas.stages.align (TGS minimap2 + NGS hisat2) on the synthetic dataset."""

import os
import sys
from types import SimpleNamespace

import pysam

from iflas.runtime import run
from iflas.stages.align import (
    filter_by_junction, hisat2_ngs, mapping_filter_and_add_tags, minimap2_tgs,
)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_synthetic_data as syn  # noqa: E402


def _introns(path):
    out = set()
    for r in pysam.AlignmentFile(path):          # SAM or BAM (autodetected)
        if r.is_unmapped:
            continue
        blocks = r.get_blocks()
        for i in range(len(blocks) - 1):
            out.add((r.reference_name, blocks[i][1], blocks[i + 1][0]))
    return out


def _bed12(chrom, starts, ends, name, strand="+"):
    chrom_start, chrom_end = starts[0], ends[-1]
    sizes = [e - s for s, e in zip(starts, ends)]
    rel_starts = [s - chrom_start for s in starts]
    return "\t".join(map(str, [
        chrom, chrom_start, chrom_end, name, 0, strand, chrom_start, chrom_end,
        "0,0,0", len(starts), ",".join(map(str, sizes)), ",".join(map(str, rel_starts)),
    ]))


def _junction_bed12(chrom, donor, acceptor, name, strand="+", score=12):
    start = donor - 10
    end = acceptor + 10
    return "\t".join(map(str, [
        chrom, start, end, name, score, strand, start, end,
        "0,0,0", 2, "10,10", "0,{}".format(acceptor - start),
    ]))


def _write_lines(path, lines):
    with open(path, "w") as f:
        for line in lines:
            f.write(line + "\n")


def _read_names(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


class _FakeRunner:
    def __init__(self, tmp_bed_lines):
        self.tmp_bed_lines = tmp_bed_lines
        self.calls = []

    def run(self, cmd, cwd=None, **kwargs):
        self.calls.append((cmd, cwd, kwargs))
        printable = cmd[-1] if isinstance(cmd, list) else cmd
        if "bamToBed" in printable:
            _write_lines(os.path.join(cwd, "tmp.bed12"), self.tmp_bed_lines)
        elif "> tmp.sam" in printable:
            open(os.path.join(cwd, "tmp.sam"), "w").close()
        elif "--output-fmt SAM" in printable:
            open(os.path.join(cwd, "flnc.mm2.sam"), "w").close()
        elif "samAddTag.pl" in printable:
            open(os.path.join(cwd, "flnc.addCVandID.bed12+"), "w").close()
        elif "samtools view -bS" in printable:
            open(os.path.join(cwd, "flnc.mm2.sorted.bam"), "w").close()
        return SimpleNamespace(returncode=0)


def test_minimap2_tgs_recovers_designed_junctions(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    sam = os.path.join(d, "flnc.mm2.sam")
    minimap2_tgs(f"{d}/genome.fa", f"{d}/flnc.fa", sam, platform="pacbio", threads=2)
    assert os.path.isfile(sam)
    introns = _introns(sam)
    for chrom, iset in syn.expected_introns().items():
        for s, e in iset:
            assert (chrom, s, e) in introns, f"missing {chrom}:{s}-{e}"
    assert ("chr1", 2200, 3000) in introns   # SE
    assert ("chr1", 4050, 4500) in introns    # A5SS
    assert ("chr1", 4600, 4950) in introns    # A3SS


def test_hisat2_ngs_maps(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    run(["hisat2-build", "-q", f"{d}/genome.fa", f"{d}/gidx"])
    bam = os.path.join(d, "ngs.bam")
    hisat2_ngs(f"{d}/gidx", bam, left=f"{d}/ngs_R1.fastq.gz", right=f"{d}/ngs_R2.fastq.gz", threads=2)
    assert os.path.isfile(bam + ".bai")
    mapped = sum(1 for r in pysam.AlignmentFile(bam, "rb") if not r.is_unmapped)
    assert mapped > 0


def test_filter_by_junction_preserves_legacy_support_rules(tmp_path):
    ref_bed = tmp_path / "ref.bed12+"
    ngs_bed = tmp_path / "junctions.bed"
    flnc_bed = tmp_path / "flnc.bed12"

    _write_lines(ref_bed, [
        _bed12("chr1", [100, 300], [200, 400], "ref_supported"),
        _bed12("chr1", [1200, 1400], [1300, 1500], "low_comb_ref"),
        _bed12("chr1", [3000, 5000], [3100, 5100], "long_ref"),
    ])
    _write_lines(ngs_bed, [
        _junction_bed12("chr1", 450, 500, "ngs_supported"),
    ])

    lines = [
        _bed12("chr1", [10], [80], "single_exon"),
        _bed12("chr1", [100, 300], [200, 400], "ref_a"),
        _bed12("chr1", [100, 300], [200, 400], "ref_b"),
        _bed12("chr1", [430, 500], [450, 550], "ngs_a"),
        _bed12("chr1", [430, 500], [450, 550], "ngs_b"),
        _bed12("chr1", [680, 800], [700, 860], "flnc_self_a"),
        _bed12("chr1", [680, 800], [700, 860], "flnc_self_b"),
        _bed12("chr1", [680, 800], [700, 860], "flnc_self_c"),
        _bed12("chr1", [880, 1000], [900, 1060], "unsupported_a"),
        _bed12("chr1", [880, 1000], [900, 1060], "unsupported_b"),
        _bed12("chr1", [1200, 1400, 1600], [1300, 1500, 1700], "low_comb"),
        _bed12("chr1", [3000, 5000], [3100, 5100], "too_long_a"),
        _bed12("chr1", [3000, 5000], [3100, 5100], "too_long_b"),
    ]
    for i in range(2):
        lines.append(_bed12(
            "chr1",
            [2000, 2200, 2400, 2600, 2800],
            [2050, 2250, 2450, 2650, 2850],
            "many_introns_{}".format(i),
        ))
    _write_lines(flnc_bed, lines)

    out = filter_by_junction(
        str(flnc_bed),
        "flnc",
        str(ref_bed),
        junc_bed_file=str(ngs_bed),
        max_intron_len=700,
        junc_comb_support_n=2,
        work_dir=str(tmp_path),
    )

    assert _read_names(out) == {
        "single_exon",
        "ref_a", "ref_b",
        "ngs_a", "ngs_b",
        "flnc_self_a", "flnc_self_b", "flnc_self_c",
        "many_introns_0", "many_introns_1",
    }


def test_filter_by_junction_limits_each_intron_not_their_sum(tmp_path):
    """A long gene whose introns are each within max_intron_len keeps its full-length reads (the published filter summed
    the introns and dropped them: 162 maize B73 transcripts; BUG-193)."""
    ref_bed, flnc_bed = tmp_path / "ref.bed12+", tmp_path / "flnc.bed12"
    starts, ends = [6000, 6300, 6600, 6900, 7200], [6100, 6400, 6700, 7000, 7300]          # 4 introns of 200 bp
    _write_lines(ref_bed, [_bed12("chr1", starts, ends, "long_gene")])
    _write_lines(flnc_bed, [_bed12("chr1", starts, ends, "long_a"), _bed12("chr1", starts, ends, "long_b"),
                            _bed12("chr1", [100, 900], [200, 1000], "one_long_intron_a"),
                            _bed12("chr1", [100, 900], [200, 1000], "one_long_intron_b")])
    out = filter_by_junction(str(flnc_bed), "flnc", str(ref_bed), max_intron_len=700, junc_comb_support_n=2,
                             work_dir=str(tmp_path))
    assert _read_names(out) == {"long_a", "long_b"}                  # 800 bp of introns in total, none over 700


def test_mapping_filter_and_add_tags_routes_legacy_commands(tmp_path):
    ref_bed = tmp_path / "ref.bed12+"
    sam = tmp_path / "flnc.input.sam"
    sam.write_text("")
    _write_lines(ref_bed, [
        _bed12("chr1", [100, 300], [200, 400], "ref_supported"),
    ])
    runner = _FakeRunner([
        _bed12("chr1", [100], [180], "single_exon"),
        _bed12("chr1", [100, 300], [200, 400], "ref_a"),
        _bed12("chr1", [100, 300], [200, 400], "ref_b"),
        _bed12("chr1", [600, 800], [700, 900], "unsupported_a"),
        _bed12("chr1", [600, 800], [700, 900], "unsupported_b"),
    ])

    out = mapping_filter_and_add_tags(
        sam_file=str(sam),
        out_prefix="flnc",
        ref_bed_file=str(ref_bed),
        threads=3,
        junc_comb_sup=2,
        work_dir=str(tmp_path),
        runner=runner,
    )

    assert _read_names(out["filtered_reads"]) == {"single_exon", "ref_a", "ref_b"}
    assert out["sam"] == os.path.join(str(tmp_path), "flnc.mm2.sam")
    assert out["bed12plus"] == os.path.join(str(tmp_path), "flnc.addCVandID.bed12+")
    assert out["bam"] == os.path.join(str(tmp_path), "flnc.mm2.sorted.bam")
    assert not os.path.exists(tmp_path / "tmp.bed12")
    assert not os.path.exists(tmp_path / "tmp.sam")

    commands = [c[0] for c in runner.calls]
    assert commands[0][:2] == ["bash", "-c"]
    assert "bamToBed" in commands[0][2]
    assert any("filter.pl" in c for c in commands if isinstance(c, str))
    assert any("samAddTag.pl" in c and "sam2bed.pl" in c for c in commands if isinstance(c, str))
    assert all(cwd == str(tmp_path) for _, cwd, _ in runner.calls)


# --------------------------------------------------------------------------- large genomes (BUG-189/190)
def test_index_bam_writes_a_csi_for_references_beyond_bai_reach(tmp_path):
    """A .bai cannot address positions past 2^29 (wheat, barley, ... chromosomes); such BAMs get a .csi (BUG-189)."""
    from iflas.runtime import ToolRunner
    from iflas.stages.align import index_bam

    for name, length, pos in (("long", 600_000_000, 550_000_000), ("short", 2_000_000, 1_000_000)):
        bam = str(tmp_path / f"{name}.bam")
        with pysam.AlignmentFile(bam, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"},
                                                    "SQ": [{"SN": "chrL", "LN": length}]}) as out:
            read = pysam.AlignedSegment(out.header)
            read.query_name, read.query_sequence, read.flag = "r1", "ACGT" * 25, 0
            read.reference_id, read.reference_start, read.cigarstring, read.mapping_quality = 0, pos, "100M", 60
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            out.write(read)
        index_bam(bam, runner=ToolRunner())
        assert os.path.isfile(bam + (".csi" if name == "long" else ".bai"))
        with pysam.AlignmentFile(bam, "rb") as f:
            assert [r.query_name for r in f.fetch("chrL", pos, pos + 50)] == ["r1"]


def test_minimap2_index_batch_keeps_large_genome_indexes_in_one_part():
    from iflas.stages.align import minimap2_index_batch

    assert minimap2_index_batch(2_300_000_000) is None                  # maize: command unchanged
    assert minimap2_index_batch(4_200_000_000) == "6G"                  # barley
    assert minimap2_index_batch(14_500_000_000) == "16G"                # wheat


def test_init_ref_builds_a_single_part_minimap2_index_for_large_genomes(tmp_path, monkeypatch):
    """Above the threshold init-ref passes -I covering the genome (the default would split the index; BUG-190)."""
    from iflas import refprep
    from iflas.runtime import ToolRunner
    from iflas.stages import align

    monkeypatch.setattr(align, "MINIMAP2_SPLIT_BASES", 100)               # the toy genome counts as large
    calls = []

    class Recorder(ToolRunner):
        def run(self, cmd, **kw):
            calls.append(cmd)
            return super().run(cmd, **kw)

    fa = tmp_path / "genome.fa"
    fa.write_text(">chr1\n" + "ACGT" * 50 + "\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('chr1\tiflas\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
    refprep.init_ref(str(fa), gtf=str(gtf), outdir=str(tmp_path), runner=Recorder())
    assert [c for c in calls if c[0] == "minimap2"] == [["minimap2", "-I", "2G", "-d", str(tmp_path / "genome.mm2"), str(fa)]]


def test_mapping_stops_on_a_multi_part_minimap2_index(tmp_path):
    """A multi-part index writes SAM without @SQ lines and each read once per part; mapping stops with what to do."""
    import pytest

    from iflas.stages.align import mapping_filter_and_add_tags, minimap2_tgs

    v2 = os.path.join(os.path.dirname(__file__), "..", "test_data", "synthetic_v2")
    index = str(tmp_path / "multi.mmi")
    run(["minimap2", "-I", "60k", "-d", index, os.path.join(v2, "genome.fa")])
    sam = minimap2_tgs(index, os.path.join(v2, "flnc.fa"), str(tmp_path / "flnc.mm2.sam"), threads=2)
    with pytest.raises(ValueError, match="split into parts"):
        mapping_filter_and_add_tags(sam_file=sam, ref_bed_file=None, work_dir=str(tmp_path / "work"))
    single = minimap2_tgs(os.path.join(v2, "genome.fa"), os.path.join(v2, "flnc.fa"), str(tmp_path / "ok.sam"),
                          threads=2, index_batch="1G")
    assert sum(line.startswith("@SQ") for line in open(single)) == 4
