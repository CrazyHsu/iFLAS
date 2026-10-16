"""Tests for the find_pa PAC clustering core, the cleavage-site / motif helpers and the CLI modes.

The full sample driver (``find_pa_sample``) is exercised on synthetic v2 in tests/test_synthetic_v2.py.
"""

import os
import shutil
from io import StringIO

import numpy as np
import pytest

from iflas.cli import main
from iflas.io.intervals import Bed12Plus
from iflas.stages import find_pa


def _bed12plus(name, chrom, start, end, strand, gene="gene1"):
    fields = [
        chrom,
        start,
        end,
        name,
        0,
        strand,
        start,
        end,
        "0,0,0",
        1,
        end - start,
        0,
        "isoform",
        "assigned",
        gene,
    ]
    return Bed12Plus("\t".join(map(str, fields)))


def test_get_pac_detects_and_ranks_separate_clusters():
    depth = np.array([0, 1, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 3, 5, 1, 0])

    peaks = find_pa.getPAC(depth, distance=3, windowSize=1, paSup=3)

    assert peaks == [
        (13, 9, [12, 13, 14]),
        (2, 7, [1, 2, 3]),
    ]


def test_get_pac_applies_rpkm_filter_to_cluster_count():
    depth = np.array([0, 1, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 3, 5, 1, 0])

    peaks = find_pa.getPAC(
        depth,
        distance=3,
        windowSize=1,
        paSup=3,
        rpkmPAC=800,
        allReadCount=1000,
        effLen=10,
    )

    assert peaks == [(13, 9, [12, 13, 14])]


def test_get_pac_includes_overlapping_right_boundary_peak_support():
    depth = np.array([0, 5, 0, 5, 0])

    peaks = find_pa.getPAC(depth, distance=3, windowSize=0, paSup=1)

    assert peaks == [(1, 10, [1, 3])]


def test_get_pac_returns_empty_when_support_is_too_low():
    depth = np.array([0, 1, 0, 1, 0])

    assert find_pa.getPAC(depth, windowSize=1, paSup=3) == []


def test_pa_cluster_writes_plus_strand_bed8plus_record():
    reads = [
        _bed12plus("read1", "chr1", 90, 105, "+"),
        _bed12plus("read2", "chr1", 91, 106, "+"),
        _bed12plus("read3", "chr1", 92, 107, "+"),
    ]
    out = StringIO()

    find_pa.paCluster(
        reads,
        distance=3,
        windowSize=1,
        outHandle=out,
        allReadCount=3,
        paSup=2,
    )

    assert out.getvalue().strip().split("\t") == [
        "chr1",
        "105",
        "107",
        "read1,read2,read3",
        "3",
        "+",
        "104",
        "105",
        "gene1",
        "Unknown",
    ]


def test_pa_cluster_marks_known_peak_from_confident_pa_dict():
    reads = [
        _bed12plus("read1", "chr1", 90, 105, "+"),
        _bed12plus("read2", "chr1", 91, 106, "+"),
        _bed12plus("read3", "chr1", 92, 107, "+"),
    ]
    out = StringIO()

    find_pa.paCluster(
        reads,
        distance=3,
        windowSize=1,
        outHandle=out,
        allReadCount=3,
        paSup=2,
        confidentPaDict={"chr1_105": ""},
    )

    assert out.getvalue().strip().split("\t")[-1] == "Known"


def test_get_pa_cluster_reads_bed12plus_and_group_file(tmp_path):
    reads_bed = tmp_path / "reads.assigned.unambi.bed12+"
    reads_bed.write_text(
        "\n".join(
            str(_bed12plus(name, "chr1", start, end, "+", gene="gene1"))
            for name, start, end in [
                ("read1", 90, 105),
                ("read2", 91, 106),
                ("read3", 92, 107),
            ]
        )
        + "\n"
    )
    group = tmp_path / "tofu.collapsed.group.txt"
    group.write_text("gene1.iso1\tread1,read2,read3\n")
    out = tmp_path / "paCluster.bed8+"

    find_pa.getPaCluster(
        readsBed=str(reads_bed),
        tofuGroup=str(group),
        filterByCount=2,
        paClusterOut=str(out),
        paDist=3,
        windowSize=1,
    )

    assert out.read_text().strip().split("\t") == [
        "chr1",
        "105",
        "107",
        "read1,read2,read3",
        "3",
        "+",
        "104",
        "105",
        "gene1",
        "Unknown",
    ]


def test_iflas_find_pa_cli_runs_ported_cluster_driver(tmp_path):
    reads_bed = tmp_path / "reads.assigned.unambi.bed12+"
    reads_bed.write_text(
        "\n".join(
            str(_bed12plus(name, "chr1", start, end, "+", gene="gene1"))
            for name, start, end in [
                ("read1", 90, 105),
                ("read2", 91, 106),
                ("read3", 92, 107),
            ]
        )
        + "\n"
    )
    group = tmp_path / "tofu.collapsed.group.txt"
    group.write_text("gene1.iso1\tread1,read2,read3\n")
    out = tmp_path / "paCluster.bed8+"

    rc = main([
        "find_pa",
        "--reads-bed",
        str(reads_bed),
        "--group",
        str(group),
        "-o",
        str(out),
        "--pa-support",
        "2",
        "--pa-dist",
        "3",
        "--window-size",
        "1",
    ])

    assert rc == 0
    assert out.read_text().strip().split("\t")[:6] == [
        "chr1",
        "105",
        "107",
        "read1,read2,read3",
        "3",
        "+",
    ]


def test_write_cleavage_sites_dedups_per_strand(tmp_path):
    reads = tmp_path / "reads.bed12+"
    reads.write_text("\n".join([
        str(_bed12plus("r1", "chr1", 100, 200, "+")),
        str(_bed12plus("r2", "chr1", 120, 200, "+")),      # same 3' end as r1
        str(_bed12plus("r3", "chr1", 100, 250, "-")),      # minus strand: 3' end = start
        str(_bed12plus("r4", "chr1", 100, 300, "-")),
    ]) + "\n")
    find_pa.write_cleavage_sites(str(reads), str(tmp_path / "cleavage.bed6"))
    assert (tmp_path / "cleavage.bed6").read_text() == (
        "chr1\t199\t200\tClevageSite1\t0\t+\n"
        "chr1\t100\t101\tClevageSite2\t0\t-\n"
    )


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="needs bedtools (pybedtools getfasta)")
def test_motif_around_pa_writes_position_profiles(tmp_path, monkeypatch):
    genome = tmp_path / "genome.fa"
    seq = ("ACGT" * 60) + "AATAAA" + ("CGTA" * 60)          # canonical PAS 6 bp before position 246
    genome.write_text(f">chr1\n{seq}\n")
    (tmp_path / "genome.len.txt").write_text(f"chr1\t{len(seq)}\n")
    pa = tmp_path / "PA.bed6+"
    pa.write_text("chr1\t245\t246\tr1,r2\t2\t+\t245\t246\tgene1\tUnknown\n")
    monkeypatch.chdir(tmp_path)
    find_pa.motifAroundPA(str(pa), up1=20, down1=20, up2=10, down2=0, refFasta=str(genome), chrLenFile=str(tmp_path / "genome.len.txt"))
    nucleotide = [l.split("\t") for l in (tmp_path / "A.nucleotide").read_text().splitlines()]
    assert [int(x[0]) for x in nucleotide] == list(range(-20, 0)) + list(range(1, 21))
    pas = {int(a): float(b) for a, b in (l.split("\t") for l in (tmp_path / "AATAAA.PAS").read_text().splitlines())}
    assert list(pas) == list(range(-10, 0))                  # up2 windows only (down2=0)
    assert pas[-5] == 1.0 and sum(pas.values()) == 1.0       # hexamer start (240) relative to the PA position (245)
    assert all((tmp_path / f"{m}.PAS").is_file() for m in find_pa.PAS_HEXAMERS)


def test_find_pa_cli_base_dir_mode_requires_reference_files(tmp_path, capsys):
    rc = main(["find_pa", "-b", str(tmp_path)])
    assert rc == 2
    assert "--ref-genome" in capsys.readouterr().err
