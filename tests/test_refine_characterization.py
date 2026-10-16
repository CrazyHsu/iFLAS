"""Characterization of refine's NGS-junction refinement and annotated/novel split, recorded from the code before the
G5 refactor (2026-10-15). The gate datasets never reach the replacement loop of ``refineWithJunc`` (no LR-only junction
is covered by NGS junctions there), so every outcome is built here: refined, two reads sharing the junction, kept for
a non-canonical motif, kept outside the +/-5 % exon-length band, dropped for an invalid chain, untouched when the
junction is itself an NGS junction or overlaps another LR junction."""
import json
import shutil

import pytest

from iflas.stages.refine import identifyNovelIsoformsByJunctions, refineWithJunc

pytestmark = pytest.mark.skipif(shutil.which("bedtools") is None, reason="needs bedtools (core iflas env)")

EXPECTED = json.loads(r'''{
 "anno": [
  "chr1\t10\t1990\ta_same\t0\t+\t10\t1990\t0,0,0\t2\t890,490\t0,1490\tPB.1",
  "chr1\t6400\t7200\ta_mono\t0\t+\t6400\t7200\t0,0,0\t1\t800\t0\tPB.1",
  "chr2\t200\t900\ta_mono_any_strand\t0\t+\t200\t900\t0,0,0\t1\t700\t0\tPB.1",
  "chr1\t3000\t5000\ta_full\t0\t+\t3000\t5000\t0,0,0\t3\t400,300,800\t0,600,1200\tPB.1"
 ],
 "filtered_junctions": [
  "chr1\t899\t1001\tjunc_900_1000\t10\t+\t899\t1001\t0,0,0\t2\t1,1\t0,101",
  "chr1\t1049\t1501\tjunc_1050_1500\t10\t+\t1049\t1501\t0,0,0\t2\t1,1\t0,451",
  "chr1\t3899\t4001\tjunc_3900_4000\t10\t+\t3899\t4001\t0,0,0\t2\t1,1\t0,101",
  "chr1\t4049\t4501\tjunc_4050_4500\t10\t+\t4049\t4501\t0,0,0\t2\t1,1\t0,451",
  "chr1\t6899\t7001\tjunc_6900_7000\t10\t+\t6899\t7001\t0,0,0\t2\t1,1\t0,101",
  "chr1\t7299\t7501\tjunc_7300_7500\t10\t+\t7299\t7501\t0,0,0\t2\t1,1\t0,201",
  "chr1\t9049\t10101\tjunc_9050_10100\t10\t+\t9049\t10101\t0,0,0\t2\t1,1\t0,1051",
  "chr1\t12399\t13001\tjunc_12400_13000\t10\t+\t12399\t13001\t0,0,0\t2\t1,1\t0,601",
  "chr1\t15399\t15951\tjunc_15400_15950\t10\t+\t15399\t15951\t0,0,0\t2\t1,1\t0,551"
 ],
 "novel": [
  "chr1\t3100\t3900\ta_sub\t0\t+\t3100\t3900\t0,0,0\t2\t300,300\t0,500\tPB.1",
  "chr1\t10\t1990\tn_shift\t0\t+\t10\t1990\t0,0,0\t2\t895,490\t0,1490\tPB.1",
  "chr1\t6900\t7300\tn_mono_small\t0\t+\t6900\t7300\t0,0,0\t1\t400\t0\tPB.1",
  "chr3\t100\t900\tn_mono_other_chrom\t0\t+\t100\t900\t0,0,0\t1\t800\t0\tPB.1"
 ],
 "refined": [
  "chr1\t0\t2000\tR1\t0\t+\t0\t2000\t0,0,0\t3\t900,50,500\t0,1000,1500\tg1",
  "chr1\t100\t1900\tR1b\t0\t+\t100\t1900\t0,0,0\t3\t800,50,400\t0,900,1400\tg1",
  "chr1\t3000\t5000\tR2\t0\t+\t3000\t5000\t0,0,0\t2\t900,500\t0,1500\tg2",
  "chr1\t6000\t8000\tR3\t0\t+\t6000\t8000\t0,0,0\t2\t900,500\t0,1500\tg3",
  "chr1\t12000\t14000\tR5\t0\t+\t12000\t14000\t0,0,0\t2\t400,1000\t0,1000\tg5",
  "chr1\t15000\t17000\tR6a\t0\t+\t15000\t17000\t0,0,0\t2\t400,1000\t0,1000\tg6",
  "chr1\t15000\t17000\tR6b\t0\t+\t15000\t17000\t0,0,0\t2\t500,1100\t0,900\tg6",
  "chr1\t17200\t17600\tR7\t0\t+\t17200\t17600\t0,0,0\t1\t400\t0\tg7"
 ]
}''')


def bed12plus(chrom, name, strand, exons, extra):
    s, e = exons[0][0], exons[-1][1]
    return "\t".join(map(str, [chrom, s, e, name, 0, strand, s, e, "0,0,0", len(exons), ",".join(str(b - a) for a, b in exons),
                               ",".join(str(a - s) for a, _ in exons)] + list(extra))) + "\n"


def gpe(name, chrom, strand, exons, gene):
    starts, ends = ",".join(str(a) for a, _ in exons) + ",", ",".join(str(b) for _, b in exons) + ","
    return "\t".join(map(str, [name, chrom, strand, exons[0][0], exons[-1][1], exons[0][0], exons[-1][1], len(exons),
                               starts, ends, 0, gene, "cmpl", "cmpl", ",".join("0" for _ in exons) + ","])) + "\n"


def junction(chrom, donor, acceptor, strand="+", score=10):
    cs, ce = donor - 1, acceptor + 1
    return "\t".join(map(str, [chrom, cs, ce, f"junc_{donor}_{acceptor}", score, strand, cs, ce, "0,0,0", 2, "1,1",
                               f"0,{acceptor - cs}"])) + "\n"


def snapshot(tmp_path):
    seq = list("A" * 18000)
    for pos, dinuc in ((900, "GT"), (998, "AG"), (1050, "GT"), (1498, "AG"),            # g1: refined
                       (6900, "GT"), (6998, "AG"), (7300, "GT"), (7498, "AG"),          # g3: outside the length band
                       (9050, "GT"), (10098, "AG"), (12400, "GT"), (12998, "AG")):      # g4: invalid chain, g5: NGS junction
        seq[pos:pos + 2] = dinuc
    (tmp_path / "genome.fa").write_text(">chr1\n" + "\n".join("".join(seq[i:i + 60]) for i in range(0, len(seq), 60)) + "\n")
    reads = [("R1", [(0, 900), (1500, 2000)], "g1"), ("R1b", [(100, 900), (1500, 1900)], "g1"),
             ("R2", [(3000, 3900), (4500, 5000)], "g2"), ("R3", [(6000, 6900), (7500, 8000)], "g3"),
             ("R4", [(9000, 9100), (9500, 10000)], "g4"), ("R5", [(12000, 12400), (13000, 14000)], "g5"),
             ("R6a", [(15000, 15400), (16000, 17000)], "g6"), ("R6b", [(15000, 15500), (15900, 17000)], "g6"),
             ("R7", [(17200, 17600)], "g7")]
    (tmp_path / "reads.bed12+").write_text("".join(bed12plus("chr1", n, "+", ex, [g]) for n, ex, g in reads))
    (tmp_path / "junctions.bed").write_text("".join(junction("chr1", a, b) for a, b in [
        (900, 1000), (1050, 1500), (3900, 4000), (4050, 4500), (6900, 7000), (7300, 7500), (9050, 10100),
        (12400, 13000), (15400, 15950)]))
    refineWithJunc(bedFile=str(tmp_path / "reads.bed12+"), junctionFile=str(tmp_path / "junctions.bed"), outFile="refined.bed12+",
                   genomeFasta=str(tmp_path / "genome.fa"))
    snap = {"refined": (tmp_path / "refined.bed12+").read_text().splitlines(),
            "filtered_junctions": (tmp_path / "filtered.junction.bed").read_text().splitlines()}

    (tmp_path / "ref2.gpe").write_text(
        gpe("tA", "chr1", "+", [(0, 900), (1500, 2000)], "gA") + gpe("tB", "chr1", "+", [(3000, 3400), (3600, 3900), (4200, 5000)], "gB")
        + gpe("tM", "chr1", "+", [(6000, 7000)], "gM") + gpe("tM2", "chr2", "-", [(100, 1100)], "gM2"))
    iso = [("chr1", "a_same", "+", [(10, 900), (1500, 1990)]), ("chr1", "a_sub", "+", [(3100, 3400), (3600, 3900)]),
           ("chr1", "n_shift", "+", [(10, 905), (1500, 1990)]), ("chr1", "a_mono", "+", [(6400, 7200)]),
           ("chr1", "n_mono_small", "+", [(6900, 7300)]), ("chr2", "a_mono_any_strand", "+", [(200, 900)]),
           ("chr3", "n_mono_other_chrom", "+", [(100, 900)]), ("chr1", "a_full", "+", [(3000, 3400), (3600, 3900), (4200, 5000)])]
    (tmp_path / "iso.bed12+").write_text("".join(bed12plus(c, n, s, ex, ["PB.1"]) for c, n, s, ex in iso))
    identifyNovelIsoformsByJunctions(str(tmp_path / "ref2.gpe"), str(tmp_path / "iso.bed12+"), anno="anno.bed", novel="novel.bed")
    snap["anno"] = (tmp_path / "anno.bed").read_text().splitlines()
    snap["novel"] = (tmp_path / "novel.bed").read_text().splitlines()
    return snap


def test_refine_matches_the_pre_refactor_record(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert snapshot(tmp_path) == EXPECTED
