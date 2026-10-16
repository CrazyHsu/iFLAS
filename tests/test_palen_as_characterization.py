"""Characterization of palen_as's two tests, recorded from the code before the G5 refactor (2026-10-15). getPalenAPA:
three sites of one gene, too few reads, one site, a 50-fold read imbalance, reads without tails, identical tails, and
the confident-site filter (PAC column 7 or reference transcript ends, +/-30 nt) the gate never runs. getPalenAS: SE
and IR pairs (the IR inclusion reads must cover the retained intron), a pair whose isoforms are not in the HQ file,
non-significant and too-small pairs, both output directories."""
import json
from types import SimpleNamespace

from iflas.stages.palen_as import getPalenAPA, getPalenAS

EXPECTED = json.loads(r'''{
 "apa_all": [
  "gA\t+\tchr1_1000;chr1_2000;chr1_3000\t12;12;12\t50.0,51.0,52.0,53.0,54.0,50.0,51.0,52.0,53.0,54.0,50.0,51.0;120.0,121.0,122.0,123.0,124.0,120.0,121.0,122.0,123.0,124.0,120.0,121.0;51.0,52.0,53.0,54.0,55.0,51.0,52.0,53.0,54.0,55.0,51.0,52.0"
 ],
 "apa_confident": [
  "gA\t+\tchr1_1000;chr1_2000\t12;12\t50.0,51.0,52.0,53.0,54.0,50.0,51.0,52.0,53.0,54.0,50.0,51.0;120.0,121.0,122.0,123.0,124.0,120.0,121.0,122.0,123.0,124.0,120.0,121.0"
 ],
 "mergeByJunc": {
  "IR.palenAndAS.sig.bed12+": [
   "chr2\t800\t1400\tPB.2.1\t0\t+\t800\t1400\t0,0,0\t1\t600\t0\tg\t.\tPB.2.1\tPB.2.1+PB.2.2\tinclusion\t30.0,31.0,32.0,33.0,34.0,30.0,31.0,32.0,33.0,34.0,30.0,31.0\t12",
   "chr2\t800\t1400\tPB.2.2\t0\t+\t800\t1400\t0,0,0\t2\t200,200\t0,400\tg\t.\tPB.2.2\tPB.2.1+PB.2.2\texclusion\t90.0,91.0,92.0,93.0,94.0,90.0,91.0,92.0,93.0,94.0,90.0,91.0\t12"
  ],
  "IR_gI_1.txt": [
   "gI\tIR\tPB.2.1\ti_ret_0\t30",
   "gI\tIR\tPB.2.1\ti_ret_1\t31",
   "gI\tIR\tPB.2.1\ti_ret_2\t32",
   "gI\tIR\tPB.2.1\ti_ret_3\t33",
   "gI\tIR\tPB.2.1\ti_ret_4\t34",
   "gI\tIR\tPB.2.1\ti_ret_5\t30",
   "gI\tIR\tPB.2.1\ti_ret_6\t31",
   "gI\tIR\tPB.2.1\ti_ret_7\t32",
   "gI\tIR\tPB.2.1\ti_ret_8\t33",
   "gI\tIR\tPB.2.1\ti_ret_9\t34",
   "gI\tIR\tPB.2.1\ti_ret_10\t30",
   "gI\tIR\tPB.2.1\ti_ret_11\t31",
   "gI\tIR\tPB.2.2\ti_exc_0\t90",
   "gI\tIR\tPB.2.2\ti_exc_1\t91",
   "gI\tIR\tPB.2.2\ti_exc_2\t92",
   "gI\tIR\tPB.2.2\ti_exc_3\t93",
   "gI\tIR\tPB.2.2\ti_exc_4\t94",
   "gI\tIR\tPB.2.2\ti_exc_5\t90",
   "gI\tIR\tPB.2.2\ti_exc_6\t91",
   "gI\tIR\tPB.2.2\ti_exc_7\t92",
   "gI\tIR\tPB.2.2\ti_exc_8\t93",
   "gI\tIR\tPB.2.2\ti_exc_9\t94",
   "gI\tIR\tPB.2.2\ti_exc_10\t90",
   "gI\tIR\tPB.2.2\ti_exc_11\t91"
  ],
  "SE.palenAndAS.sig.bed12+": [
   "chr1\t100\t1100\tPB.1.1\t0\t+\t100\t1100\t0,0,0\t3\t200,200,200\t0,400,800\tg\t.\tPB.1.1\tPB.1.1+PB.1.2\tinclusion\t40.0,41.0,42.0,43.0,44.0,40.0,41.0,42.0,43.0,44.0,40.0,41.0\t12",
   "chr1\t100\t1100\tPB.1.2\t0\t+\t100\t1100\t0,0,0\t2\t200,200\t0,800\tg\t.\tPB.1.2\tPB.1.1+PB.1.2\texclusion\t100.0,101.0,102.0,103.0,104.0,100.0,101.0,102.0,103.0,104.0,100.0,101.0\t12",
   "chr4\t50\t800\tPB.4.1_PB.4.3\t0\t+\t50\t800\t0,0,0\t2\t250,300\t0,450\tg\t.\tPB.4.3\tPB.4.1_PB.4.3+PB.4.2\tinclusion\t20.0,21.0,22.0,23.0,24.0,20.0,22.0,23.0,24.0,25.0,26.0,22.0\t12",
   "chr4\t100\t700\tPB.4.2\t0\t+\t100\t700\t0,0,0\t2\t200,100\t0,500\tg\t.\tPB.4.2\tPB.4.1_PB.4.3+PB.4.2\texclusion\t80.0,81.0,82.0,83.0,84.0,80.0,81.0,82.0,83.0,84.0,80.0,81.0\t12"
  ],
  "SE_gM_3.txt": [
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_0\t20",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_1\t21",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_2\t22",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_3\t23",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_4\t24",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_5\t20",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_0\t22",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_1\t23",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_2\t24",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_3\t25",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_4\t26",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_5\t22",
   "gM\tSE\tPB.4.2\tm_exc_0\t80",
   "gM\tSE\tPB.4.2\tm_exc_1\t81",
   "gM\tSE\tPB.4.2\tm_exc_2\t82",
   "gM\tSE\tPB.4.2\tm_exc_3\t83",
   "gM\tSE\tPB.4.2\tm_exc_4\t84",
   "gM\tSE\tPB.4.2\tm_exc_5\t80",
   "gM\tSE\tPB.4.2\tm_exc_6\t81",
   "gM\tSE\tPB.4.2\tm_exc_7\t82",
   "gM\tSE\tPB.4.2\tm_exc_8\t83",
   "gM\tSE\tPB.4.2\tm_exc_9\t84",
   "gM\tSE\tPB.4.2\tm_exc_10\t80",
   "gM\tSE\tPB.4.2\tm_exc_11\t81"
  ],
  "SE_gS_1.txt": [
   "gS\tSE\tPB.1.1\ts_inc_0\t40",
   "gS\tSE\tPB.1.1\ts_inc_1\t41",
   "gS\tSE\tPB.1.1\ts_inc_2\t42",
   "gS\tSE\tPB.1.1\ts_inc_3\t43",
   "gS\tSE\tPB.1.1\ts_inc_4\t44",
   "gS\tSE\tPB.1.1\ts_inc_5\t40",
   "gS\tSE\tPB.1.1\ts_inc_6\t41",
   "gS\tSE\tPB.1.1\ts_inc_7\t42",
   "gS\tSE\tPB.1.1\ts_inc_8\t43",
   "gS\tSE\tPB.1.1\ts_inc_9\t44",
   "gS\tSE\tPB.1.1\ts_inc_10\t40",
   "gS\tSE\tPB.1.1\ts_inc_11\t41",
   "gS\tSE\tPB.1.2\ts_exc_0\t100",
   "gS\tSE\tPB.1.2\ts_exc_1\t101",
   "gS\tSE\tPB.1.2\ts_exc_2\t102",
   "gS\tSE\tPB.1.2\ts_exc_3\t103",
   "gS\tSE\tPB.1.2\ts_exc_4\t104",
   "gS\tSE\tPB.1.2\ts_exc_5\t100",
   "gS\tSE\tPB.1.2\ts_exc_6\t101",
   "gS\tSE\tPB.1.2\ts_exc_7\t102",
   "gS\tSE\tPB.1.2\ts_exc_8\t103",
   "gS\tSE\tPB.1.2\ts_exc_9\t104",
   "gS\tSE\tPB.1.2\ts_exc_10\t100",
   "gS\tSE\tPB.1.2\ts_exc_11\t101"
  ],
  "SE_gS_2.txt": [
   "gS\tSE\tPB.9.9\ts_far_0\t40",
   "gS\tSE\tPB.9.9\ts_far_1\t41",
   "gS\tSE\tPB.9.9\ts_far_2\t42",
   "gS\tSE\tPB.9.9\ts_far_3\t43",
   "gS\tSE\tPB.9.9\ts_far_4\t44",
   "gS\tSE\tPB.9.9\ts_far_5\t40",
   "gS\tSE\tPB.9.9\ts_far_6\t41",
   "gS\tSE\tPB.9.9\ts_far_7\t42",
   "gS\tSE\tPB.9.9\ts_far_8\t43",
   "gS\tSE\tPB.9.9\ts_far_9\t44",
   "gS\tSE\tPB.9.9\ts_far_10\t40",
   "gS\tSE\tPB.9.9\ts_far_11\t41",
   "gS\tSE\tPB.1.2\ts_exc_0\t100",
   "gS\tSE\tPB.1.2\ts_exc_1\t101",
   "gS\tSE\tPB.1.2\ts_exc_2\t102",
   "gS\tSE\tPB.1.2\ts_exc_3\t103",
   "gS\tSE\tPB.1.2\ts_exc_4\t104",
   "gS\tSE\tPB.1.2\ts_exc_5\t100",
   "gS\tSE\tPB.1.2\ts_exc_6\t101",
   "gS\tSE\tPB.1.2\ts_exc_7\t102",
   "gS\tSE\tPB.1.2\ts_exc_8\t103",
   "gS\tSE\tPB.1.2\ts_exc_9\t104",
   "gS\tSE\tPB.1.2\ts_exc_10\t100",
   "gS\tSE\tPB.1.2\ts_exc_11\t101"
  ]
 },
 "no_mergeByJunc": {
  "IR.palenAndAS.sig.bed12+": [
   "chr2\t800\t1400\tPB.2.1\t0\t+\t800\t1400\t0,0,0\t1\t600\t0\tg\t.\tPB.2.1\tPB.2.1+PB.2.2\tinclusion\t30.0,31.0,32.0,33.0,34.0,30.0,31.0,32.0,33.0,34.0,30.0,31.0\t12",
   "chr2\t800\t1400\tPB.2.2\t0\t+\t800\t1400\t0,0,0\t2\t200,200\t0,400\tg\t.\tPB.2.2\tPB.2.1+PB.2.2\texclusion\t90.0,91.0,92.0,93.0,94.0,90.0,91.0,92.0,93.0,94.0,90.0,91.0\t12"
  ],
  "IR_gI_1.txt": [
   "gI\tIR\tPB.2.1\ti_ret_0\t30",
   "gI\tIR\tPB.2.1\ti_ret_1\t31",
   "gI\tIR\tPB.2.1\ti_ret_2\t32",
   "gI\tIR\tPB.2.1\ti_ret_3\t33",
   "gI\tIR\tPB.2.1\ti_ret_4\t34",
   "gI\tIR\tPB.2.1\ti_ret_5\t30",
   "gI\tIR\tPB.2.1\ti_ret_6\t31",
   "gI\tIR\tPB.2.1\ti_ret_7\t32",
   "gI\tIR\tPB.2.1\ti_ret_8\t33",
   "gI\tIR\tPB.2.1\ti_ret_9\t34",
   "gI\tIR\tPB.2.1\ti_ret_10\t30",
   "gI\tIR\tPB.2.1\ti_ret_11\t31",
   "gI\tIR\tPB.2.2\ti_exc_0\t90",
   "gI\tIR\tPB.2.2\ti_exc_1\t91",
   "gI\tIR\tPB.2.2\ti_exc_2\t92",
   "gI\tIR\tPB.2.2\ti_exc_3\t93",
   "gI\tIR\tPB.2.2\ti_exc_4\t94",
   "gI\tIR\tPB.2.2\ti_exc_5\t90",
   "gI\tIR\tPB.2.2\ti_exc_6\t91",
   "gI\tIR\tPB.2.2\ti_exc_7\t92",
   "gI\tIR\tPB.2.2\ti_exc_8\t93",
   "gI\tIR\tPB.2.2\ti_exc_9\t94",
   "gI\tIR\tPB.2.2\ti_exc_10\t90",
   "gI\tIR\tPB.2.2\ti_exc_11\t91"
  ],
  "SE.palenAndAS.sig.bed12+": [
   "chr1\t100\t1100\tPB.1.1\t0\t+\t100\t1100\t0,0,0\t3\t200,200,200\t0,400,800\tg\t.\tPB.1.1\tPB.1.1+PB.1.2\tinclusion\t40.0,41.0,42.0,43.0,44.0,40.0,41.0,42.0,43.0,44.0,40.0,41.0\t12",
   "chr1\t100\t1100\tPB.1.2\t0\t+\t100\t1100\t0,0,0\t2\t200,200\t0,800\tg\t.\tPB.1.2\tPB.1.1+PB.1.2\texclusion\t100.0,101.0,102.0,103.0,104.0,100.0,101.0,102.0,103.0,104.0,100.0,101.0\t12",
   "chr4\t50\t800\tPB.4.1_PB.4.3\t0\t+\t50\t800\t0,0,0\t2\t250,300\t0,450\tg\t.\tPB.4.3\tPB.4.1_PB.4.3+PB.4.2\tinclusion\t20.0,21.0,22.0,23.0,24.0,20.0,22.0,23.0,24.0,25.0,26.0,22.0\t12",
   "chr4\t100\t700\tPB.4.2\t0\t+\t100\t700\t0,0,0\t2\t200,100\t0,500\tg\t.\tPB.4.2\tPB.4.1_PB.4.3+PB.4.2\texclusion\t80.0,81.0,82.0,83.0,84.0,80.0,81.0,82.0,83.0,84.0,80.0,81.0\t12"
  ],
  "SE_gM_3.txt": [
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_0\t20",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_1\t21",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_2\t22",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_3\t23",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_4\t24",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc_5\t20",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_0\t22",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_1\t23",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_2\t24",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_3\t25",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_4\t26",
   "gM\tSE\tPB.4.1_PB.4.3\tm_inc3_5\t22",
   "gM\tSE\tPB.4.2\tm_exc_0\t80",
   "gM\tSE\tPB.4.2\tm_exc_1\t81",
   "gM\tSE\tPB.4.2\tm_exc_2\t82",
   "gM\tSE\tPB.4.2\tm_exc_3\t83",
   "gM\tSE\tPB.4.2\tm_exc_4\t84",
   "gM\tSE\tPB.4.2\tm_exc_5\t80",
   "gM\tSE\tPB.4.2\tm_exc_6\t81",
   "gM\tSE\tPB.4.2\tm_exc_7\t82",
   "gM\tSE\tPB.4.2\tm_exc_8\t83",
   "gM\tSE\tPB.4.2\tm_exc_9\t84",
   "gM\tSE\tPB.4.2\tm_exc_10\t80",
   "gM\tSE\tPB.4.2\tm_exc_11\t81"
  ],
  "SE_gS_1.txt": [
   "gS\tSE\tPB.1.1\ts_inc_0\t40",
   "gS\tSE\tPB.1.1\ts_inc_1\t41",
   "gS\tSE\tPB.1.1\ts_inc_2\t42",
   "gS\tSE\tPB.1.1\ts_inc_3\t43",
   "gS\tSE\tPB.1.1\ts_inc_4\t44",
   "gS\tSE\tPB.1.1\ts_inc_5\t40",
   "gS\tSE\tPB.1.1\ts_inc_6\t41",
   "gS\tSE\tPB.1.1\ts_inc_7\t42",
   "gS\tSE\tPB.1.1\ts_inc_8\t43",
   "gS\tSE\tPB.1.1\ts_inc_9\t44",
   "gS\tSE\tPB.1.1\ts_inc_10\t40",
   "gS\tSE\tPB.1.1\ts_inc_11\t41",
   "gS\tSE\tPB.1.2\ts_exc_0\t100",
   "gS\tSE\tPB.1.2\ts_exc_1\t101",
   "gS\tSE\tPB.1.2\ts_exc_2\t102",
   "gS\tSE\tPB.1.2\ts_exc_3\t103",
   "gS\tSE\tPB.1.2\ts_exc_4\t104",
   "gS\tSE\tPB.1.2\ts_exc_5\t100",
   "gS\tSE\tPB.1.2\ts_exc_6\t101",
   "gS\tSE\tPB.1.2\ts_exc_7\t102",
   "gS\tSE\tPB.1.2\ts_exc_8\t103",
   "gS\tSE\tPB.1.2\ts_exc_9\t104",
   "gS\tSE\tPB.1.2\ts_exc_10\t100",
   "gS\tSE\tPB.1.2\ts_exc_11\t101"
  ],
  "SE_gS_2.txt": [
   "gS\tSE\tPB.9.9\ts_far_0\t40",
   "gS\tSE\tPB.9.9\ts_far_1\t41",
   "gS\tSE\tPB.9.9\ts_far_2\t42",
   "gS\tSE\tPB.9.9\ts_far_3\t43",
   "gS\tSE\tPB.9.9\ts_far_4\t44",
   "gS\tSE\tPB.9.9\ts_far_5\t40",
   "gS\tSE\tPB.9.9\ts_far_6\t41",
   "gS\tSE\tPB.9.9\ts_far_7\t42",
   "gS\tSE\tPB.9.9\ts_far_8\t43",
   "gS\tSE\tPB.9.9\ts_far_9\t44",
   "gS\tSE\tPB.9.9\ts_far_10\t40",
   "gS\tSE\tPB.9.9\ts_far_11\t41",
   "gS\tSE\tPB.1.2\ts_exc_0\t100",
   "gS\tSE\tPB.1.2\ts_exc_1\t101",
   "gS\tSE\tPB.1.2\ts_exc_2\t102",
   "gS\tSE\tPB.1.2\ts_exc_3\t103",
   "gS\tSE\tPB.1.2\ts_exc_4\t104",
   "gS\tSE\tPB.1.2\ts_exc_5\t100",
   "gS\tSE\tPB.1.2\ts_exc_6\t101",
   "gS\tSE\tPB.1.2\ts_exc_7\t102",
   "gS\tSE\tPB.1.2\ts_exc_8\t103",
   "gS\tSE\tPB.1.2\ts_exc_9\t104",
   "gS\tSE\tPB.1.2\ts_exc_10\t100",
   "gS\tSE\tPB.1.2\ts_exc_11\t101"
  ]
 }
}''')


def names(prefix, n):
    return [f"{prefix}_{i}" for i in range(n)]


def spread(base, reads):
    return {r: base + (i % 5) for i, r in enumerate(reads)}


def bed12plus(name, chrom, strand, exons, extra):
    s, e = exons[0][0], exons[-1][1]
    return "\t".join(map(str, [chrom, s, e, name, 0, strand, s, e, "0,0,0", len(exons), ",".join(str(b - a) for a, b in exons),
                               ",".join(str(a - s) for a, _ in exons)] + list(extra))) + "\n"


def apa_snapshot(tmp_path):
    sites = [("chr1", 1000, "+", "gA", names("a1", 12), 50), ("chr1", 2000, "+", "gA", names("a2", 12), 120),
             ("chr1", 3000, "+", "gA", names("a3", 12), 51), ("chr1", 4000, "-", "gB", names("b1", 5), 50),
             ("chr1", 4100, "-", "gB", names("b2", 5), 120), ("chr2", 100, "+", "gC", names("c1", 12), 50),
             ("chr2", 5000, "+", "gD", names("d1", 500), 50), ("chr2", 6000, "+", "gD", names("d2", 10), 120),
             ("chr3", 100, "+", "gE", names("e1", 12), 50), ("chr3", 200, "+", "gE", names("e2", 12), 120),
             ("chr3", 300, "+", "gF", names("f1", 12), 60), ("chr3", 400, "+", "gF", names("f2", 12), 60)]
    tails = {}
    for chrom, pos, strand, gene, reads, base in sites:
        if gene == "gE" and pos == 100:
            reads = reads[:3]                                   # only three reads of this site have a tail
        tails.update({r: base for r in reads} if gene == "gF" else spread(base, reads))
    pac = tmp_path / "paCluster.bed8+"
    pac.write_text("".join(f"{c}\t{p - 1}\t{p}\t{','.join(r)}\t{len(r)}\t{s}\t.\t{p}\t{g}\n" for c, p, s, g, r, _ in sites))
    (tmp_path / "confident.pac").write_text("chr1\t990\t1000\tpac1\t0\t+\t1000\n")
    (tmp_path / "ref.bed").write_text("chr1\t1500\t2010\ttx1\t0\t+\nchr2\t5990\t7000\ttx2\t0\t-\n")
    snap = {}
    for tag, kw in (("all", {}), ("confident", {"confidentPacFiles": [str(tmp_path / "confident.pac")],
                                                "refParams": SimpleNamespace(ref_bed=str(tmp_path / "ref.bed"))})):
        out = tmp_path / f"apa_{tag}.txt"
        getPalenAPA(str(pac), tails, 10, str(out), **kw)
        snap["apa_" + tag] = out.read_text().splitlines()
    return snap


def as_snapshot(tmp_path):
    iso = [("PB.1.1", "chr1", [(100, 300), (500, 700), (900, 1100)], "PB.1.1"), ("PB.1.2", "chr1", [(100, 300), (900, 1100)], "PB.1.2"),
           ("PB.1.3", "chr1", [(100, 300), (500, 700), (900, 1000)], "PB.9.9"),
           ("PB.2.1", "chr2", [(800, 1400)], "PB.2.1"), ("PB.2.2", "chr2", [(800, 1000), (1200, 1400)], "PB.2.2"),
           ("PB.3.1", "chr3", [(100, 300), (500, 700)], "PB.3.1"), ("PB.3.2", "chr3", [(100, 300), (550, 700)], "PB.3.2"),
           ("PB.4.1", "chr4", [(100, 300), (500, 700)], "PB.4.1,PB.4.3"), ("PB.4.2", "chr4", [(100, 300), (600, 700)], "PB.4.2"),
           ("PB.4.3", "chr4", [(50, 300), (500, 800)], "PB.4.3"), ("PB.5.1", "chr5", [(100, 300)], "PB.5.1"),
           ("PB.5.2", "chr5", [(100, 200)], "PB.5.2")]
    (tmp_path / "hq.bed12+").write_text("".join(bed12plus(n, c, "+", ex, ["g", ".", cup]) for n, c, ex, cup in iso))
    groups = {"PB.1.1": names("s_inc", 12), "PB.1.2": names("s_exc", 12), "PB.9.9": names("s_far", 12),
              "PB.2.1": names("i_ret", 12) + names("i_spl", 3), "PB.2.2": names("i_exc", 12),
              "PB.3.1": names("t_inc", 12), "PB.3.2": names("t_exc", 12), "PB.4.1": names("m_inc", 6),
              "PB.4.3": names("m_inc3", 6), "PB.4.2": names("m_exc", 12), "PB.5.1": names("u_inc", 5), "PB.5.2": names("u_exc", 5)}
    reads = [bed12plus(r, "chr2", "+", [(800, 1400)], []) for r in names("i_ret", 12)]
    reads += [bed12plus(r, "chr2", "+", [(800, 1000), (1200, 1400)], []) for r in names("i_spl", 3)]
    (tmp_path / "reads.bed12+").write_text("".join(reads))
    tails = {}
    for prefix, base in (("s_inc", 40), ("s_exc", 100), ("s_far", 40), ("i_ret", 30), ("i_spl", 30), ("i_exc", 90),
                         ("t_inc", 60), ("t_exc", 61), ("m_inc", 20), ("m_inc3", 22), ("m_exc", 80), ("u_inc", 10), ("u_exc", 90)):
        tails.update(spread(base, groups[{"s_inc": "PB.1.1", "s_exc": "PB.1.2", "s_far": "PB.9.9", "i_ret": "PB.2.1", "i_spl": "PB.2.1",
                                          "i_exc": "PB.2.2", "t_inc": "PB.3.1", "t_exc": "PB.3.2", "m_inc": "PB.4.1", "m_inc3": "PB.4.3",
                                          "m_exc": "PB.4.2", "u_inc": "PB.5.1", "u_exc": "PB.5.2"}[prefix]]))
    asPairs = {"SE": {"gS": [[["PB.1.1"], ["PB.1.2"], "gS:300@500-700@900"], [["PB.1.3"], ["PB.1.2"], "gS:300@500-700@900"]],
                      "gT": [[["PB.3.1"], ["PB.3.2"], "gT:300@500-700@1000"]],
                      "gM": [[["PB.4.1"], ["PB.4.2"], "gM:300@500-700@1000"]],
                      "gU": [[["PB.5.1"], ["PB.5.2"], "gU:300@500-700@1000"]]},
               "IR": {"gI": [[["PB.2.1"], ["PB.2.2"], "gI:chr2:1000-1200"]]}}
    snap = {}
    for merge in (False, True):
        getPalenAS(tails, str(tmp_path / "hq.bed12+"), str(tmp_path / "reads.bed12+"), collapsedTrans2reads=groups, asPairs=asPairs,
                   filterByCount=10, mergeByJunc=merge)
        mode = tmp_path / ("mergeByJunc" if merge else "no_mergeByJunc")
        snap[mode.name] = {p.name: p.read_text().splitlines() for p in sorted(mode.iterdir())}
    return snap


def snapshot(tmp_path):
    return {**apa_snapshot(tmp_path), **as_snapshot(tmp_path)}


def test_palen_as_matches_the_pre_refactor_record(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert snapshot(tmp_path) == EXPECTED
