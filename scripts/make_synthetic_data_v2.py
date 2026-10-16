#!/usr/bin/env python
"""Synthetic iFLAS dataset v2: diverse genome/gene structures with machine-readable truth tables.

Everything is derived from the declarative gene catalog in
``CATALOG`` below; coordinates are laid out by a cursor per chromosome, so adding a gene never requires
hand-computing positions.  Formats of the core files (genome.fa / genes.gtf / flnc.fa / ngs/*.fastq.gz)
are identical to v1 so iFLAS 1 (Python 2) can read them unchanged.

Usage:
  python scripts/make_synthetic_data_v2.py test_data/synthetic_v2 [--seed 20261014] [--scale 1.0]
                                            [--error-rate 0.003] [--stranded]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import random
import uuid
import warnings
from dataclasses import dataclass, field

RC = str.maketrans("ACGTacgt", "TGCAtgca")
CONDITIONS = ("condA", "condB")
REPLICATES = ("rep1", "rep2")
PRIMER_5P = "AAGCAGTGGTATCAACGCAGAGTACATGGGG"   # official Iso-Seq Clontech_5p: the extra ATGGGG keeps it distinct from revcomp(3p),
PRIMER_3P = "GTACTCTGCGTTGATACCACTGCTT"         # official NEB_Clontech_3p (lima rejects primer pairs that are reverse complements)
MOVIE = "m64000_000000_000000"                 # PacBio movie name: RG PU and the read-name prefix lima/isoseq keep
ISOSEQ_MIN_POLYA = 20                           # isoseq refine --require-polya default --min-polya-length
MOTIFS = {"GT-AG": ("GT", "AG"), "GC-AG": ("GC", "AG"), "AT-AC": ("AT", "AC"), "GA-AG": ("GA", "AG")}
DEFAULT_PALEN = (60.0, 15.0)
LR_MIN_READS = 2          # find_as drops isoforms supported by < 2 reads
MICROEXON_MIN = 12        # minimap2 -ax splice smears a 9-bp internal exon into indels (unrecoverable); 15 bp is still recovered via NGS-refined junctions
AS_OFFSET = 5             # find_as A5SS/A3SS minimum site difference
JUNC_COMB_SUPPORT = 2     # legacy `mapping -jcs`: FLNC reads sharing a junction chain needed to keep them
NGS_JUNC_MIN_READS = 10   # legacy getJuncFromRegtools keeps NGS junctions with > 10 reads (regtools -a 5)
NGS_ANCHOR = 5            # designed NGS support = reads with >= 5 bp on both exons (what regtools -a 5 counts)
NGS_BAND = (7, 18)        # a junction whose filter consensus hinges on NGS must not have designed support in this
                          # open band: the aligner may add/lose a read or two around the > 10 decision


def revcomp(s: str) -> str:
    return s.translate(RC)[::-1]


# ----------------------------------------------------------------------------- catalog
@dataclass
class Gene:
    """One locus. ``exons``/``introns`` are lengths; ``variants`` are (suffix, op, n_flnc, annotated)."""
    name: str
    chrom: str
    strand: str
    exons: list
    introns: list
    variants: list = field(default_factory=list)
    n_full: int = 8                 # FLNC reads of the full-structure isoform (0 = no full isoform)
    full_annotated: bool = True
    category: str = "as"
    motif: str = "GT-AG"            # motif of every intron unless ``motif_by_intron`` overrides
    motif_by_intron: dict = field(default_factory=dict)
    softmask: bool = False
    gap: int = 1200                 # distance from previous gene end (negative = overlap it)
    pa_sites: tuple = ((0, 1.0),)   # (offset downstream of the canonical TES, fraction of reads)
    snps: int = 0                   # heterozygous exonic SNPs (diploid / ASAS genes)
    hap_bias: float = 0.5           # P(hapA | inclusion isoform); exclusion gets 1 - bias (allele-specific splicing)
    expr_bias: float = None         # P(hapA) for every isoform (allele-specific expression, no splicing bias)
    ngs_psi: dict = field(default_factory=dict)   # condition -> designed NGS PSI of the gene's event
    ngs_weight: dict = field(default_factory=dict)  # isoform suffix -> NGS depth multiplier (0 = LR-only isoform)
    palen: dict = field(default_factory=dict)     # isoform suffix -> (mean, sd) of poly(A) length
    note: str = ""


def _pu_background(n=30):
    """30 canonical/novel pairs whose read counts overlap: every 4th novel isoform dominates its gene
    (incomplete annotation), two novel isoforms have a single read (find_as min-read filter)."""
    genes = []
    for i in range(1, n + 1):
        full_n = 2 + (i * 7) % 19                      # 2..20 gradient
        novel_n = 1 + (i * 5) % 12                     # 1..12 gradient
        if i % 4 == 0:
            full_n, novel_n = 2 + i % 4, 14 + (i * 3) % 7   # novel-dominant genes
        op = ("skip", 1) if i % 3 else ("end", 0, 20)  # every third pair: novel junction instead of SE
        genes.append(Gene(f"genePU{i:02d}", "chr3", "+" if i % 2 else "-", [120, 120, 130], [180, 230],
                          [("novel", op, novel_n, False)], n_full=full_n, category="pu_background",
                          softmask=(i % 7 == 0), gap=400))
    return genes


CATALOG = [
    # ---- chr1: AS core, strands alternate -------------------------------------------------------
    Gene("seA", "chr1", "+", [300, 150, 200, 250], [400, 350, 500], [("skip", ("skip", 1), 6, False)], 10,
         ngs_psi={"condA": 0.8, "condB": 0.2}, palen={"full": DEFAULT_PALEN, "skip": DEFAULT_PALEN}),
    Gene("seB", "chr1", "-", [250, 120, 180, 220], [600, 300, 450], [("skip", ("skip", 1), 8, True)], 8),
    Gene("seC", "chr1", "+", [200, 100, 120, 150, 300], [350, 300, 400, 500], [("skip", ("skip", 1, 2), 5, False)], 7,
         ngs_psi={"condA": 0.7, "condB": 0.3}),
    Gene("seD", "chr1", "-", [220, 130, 140, 260], [300, 800, 350],
         [("skip1", ("skip", 1), 4, True), ("skip2", ("skip", 2), 3, False)], 9),
    # poly(A) tail-length genes: the published palen_as rule needs t-test AND Kruskal p <= 0.001, which small groups
    # cannot reach (n = 5 vs 5: Kruskal minimum p = 0.0079; R-16) -> 15 reads per isoform
    Gene("irA", "chr1", "+", [300, 200, 250], [400, 350], [("ir", ("ir", 0), 15, False)], 15,
         palen={"ir": (80.0, 15.0), "full": (40.0, 10.0)}),
    Gene("irB", "chr1", "-", [280, 150, 300], [500, 420], [("ir", ("ir", 1), 6, True)], 8, ngs_psi={"condA": 0.6, "condB": 0.4}),
    Gene("irC", "chr1", "+", [200, 180, 160, 240], [300, 300, 300], [("ir", ("ir", 1), 4, False)], 6, softmask=True),
    Gene("irD", "chr1", "-", [210, 170, 150, 230], [320, 310, 330], [("ir", ("ir", 0), 15, False)], 15,
         palen={"ir": (85.0, 15.0), "full": (45.0, 10.0)}),
    Gene("a5A", "chr1", "+", [300, 200, 250], [500, 400], [("alt", ("end", 0, 30), 6, False)], 9),
    Gene("a5B", "chr1", "-", [260, 220, 240], [450, 380], [("alt", ("start", 1, -60), 5, True)], 8, ngs_psi={"condA": 0.8, "condB": 0.2}),
    Gene("a5C", "chr1", "+", [280, 210, 260], [420, 410], [("alt", ("end", 0, 3), 6, False)], 6),
    Gene("a3A", "chr1", "+", [300, 200, 250], [500, 400], [("alt", ("start", 1, -30), 4, False)], 10, ngs_psi={"condA": 0.7, "condB": 0.3}),
    Gene("a3B", "chr1", "-", [270, 190, 230], [480, 390], [("alt", ("end", 0, 60), 15, True)], 15,
         palen={"alt": (90.0, 15.0), "full": (40.0, 10.0)}),
    Gene("a3C", "chr1", "+", [250, 200, 300], [400, 400], [("alt", ("start", 1, -3), 6, False)], 6),
    Gene("mxeA", "chr1", "+", [250, 150, 150, 300], [400, 300, 400],
         [("mxe1", ("skip", 2), 8, True), ("mxe2", ("skip", 1), 8, False)], 0),
    Gene("mxeB", "chr1", "-", [260, 140, 160, 280], [380, 320, 420],
         [("mxe1", ("skip", 2), 7, True), ("mxe2", ("skip", 1), 7, False)], 0),
    Gene("cpxA", "chr1", "+", [300, 250, 200], [400, 500], [("cpx", ("complex", 0, 40, 40), 6, False)], 8, ngs_psi={"condA": 0.6, "condB": 0.4}),
    Gene("cpxB", "chr1", "-", [280, 240, 220], [420, 480], [("cpx", ("complex", 0, 35, 45), 6, False)], 8),
    Gene("irseA", "chr1", "+", [250, 150, 200, 220], [350, 300, 450], [("ir", ("ir", 2), 5, False), ("skip", ("skip", 1), 4, False)], 10),
    Gene("tssA", "chr1", "+", [300, 200, 300], [400, 400], [("tss", ("tss", -150), 5, False)], 8, pa_sites=((0, 0.6), (120, 0.4))),
    Gene("tssB", "chr1", "-", [320, 180, 280], [380, 420], [("tss", ("tss", -130), 4, False)], 8, pa_sites=((0, 0.5), (90, 0.3), (200, 0.2))),
    Gene("longA", "chr1", "+", [300, 200, 250], [16000, 400], [("skip", ("skip", 1), 4, False)], 6, note="long_intron"),
    Gene("longB", "chr1", "-", [260, 180, 300], [17000, 450], [("alt", ("end", 0, 30), 4, False)], 6, note="long_intron"),
    Gene("wobbleA", "chr1", "+", [300, 200, 250], [500, 400], [("wob", ("wobble", 0, "end", 3), 3, False)], 10),
    Gene("wobbleB", "chr1", "-", [280, 220, 260], [460, 440], [("wob", ("wobble", 1, "start", -7), 3, False)], 9),
    Gene("wobbleC", "chr1", "+", [290, 210, 240], [430, 470], [("wob", ("wobble", 1, "start", 10), 3, False)], 9),
    Gene("lowA", "chr1", "+", [300, 200, 250], [400, 400], [("skip", ("skip", 1), 1, False)], 8),
    Gene("lowB", "chr1", "-", [280, 190, 240], [410, 390], [("skip", ("skip", 1), 2, False)], 5, ngs_weight={"skip": 0.0},
         note="2-read LR-only novel isoform: fails the legacy junction-combination consensus in every NGS configuration"),
    Gene("nullA", "chr1", "+", [300, 160, 260], [450, 350], [("skip", ("skip", 1), 8, True)], 8, ngs_psi={"condA": 0.5, "condB": 0.5}),
    Gene("nullB", "chr1", "-", [290, 170, 250], [440, 360], [("skip", ("skip", 1), 7, True)], 7, ngs_psi={"condA": 0.5, "condB": 0.5}),
    # ---- chr2: structural edge cases ------------------------------------------------------------
    Gene("monoA", "chr2", "+", [800], [], [], 8, category="mono_exon"),
    Gene("monoB", "chr2", "-", [600], [], [], 6, category="mono_exon", softmask=True),
    Gene("monoC", "chr2", "+", [500], [], [], 3, full_annotated=False, category="mono_exon", note="novel_mono"),
    Gene("diA", "chr2", "+", [300, 400], [500], [("ir", ("ir", 0), 5, False)], 8, category="di_exon"),
    Gene("diB", "chr2", "-", [350, 300], [600], [], 6, category="di_exon"),
    Gene("microA", "chr2", "+", [250, 30, 300], [400, 350], [("skip", ("skip", 1), 4, False)], 8, note="small_exon_30bp"),
    Gene("microB", "chr2", "-", [260, 9, 290], [380, 360], [("skip", ("skip", 1), 4, False)], 8, note="micro_exon_9bp"),
    Gene("microC", "chr2", "+", [240, 15, 280], [390, 370], [("skip", ("skip", 1), 4, False)], 8, note="micro_exon_15bp"),
    Gene("antiA", "chr2", "+", [300, 200, 300], [500, 400], [("skip", ("skip", 1), 4, False)], 8, category="antisense_pair"),
    Gene("antiB", "chr2", "-", [250, 250], [700], [], 7, gap=-1400, category="antisense_pair", note="overlaps antiA"),
    Gene("antiC", "chr2", "+", [300, 220, 280], [600, 450], [], 8, category="antisense_pair"),
    Gene("antiD", "chr2", "-", [200, 240], [600], [], 7, gap=-1520, category="antisense_pair", note="shares antiC intron 0 coordinates"),
    Gene("nearA", "chr2", "+", [280, 200, 260], [400, 380], [("skip", ("skip", 1), 4, False)], 8, category="adjacent"),
    Gene("nearB", "chr2", "+", [270, 210, 250], [390, 410], [], 8, gap=200, category="adjacent"),
    Gene("noncanA", "chr2", "+", [300, 200, 250], [400, 400], [], 7, motif_by_intron={0: "GA-AG"}, note="noncanonical_intron0"),
    Gene("gcagA", "chr2", "-", [280, 200, 240], [420, 430], [("skip", ("skip", 1), 3, False)], 6, motif="GC-AG"),
    Gene("atacA", "chr2", "+", [290, 190, 230], [410, 440], [], 6, motif="AT-AC"),
    # ASAS: the published test (asas.relationshipBetweenAlleleSpecifiAndAS) needs, per isoform, >= 5 reads on one
    # haplotype and <= 3 on the other, plus chi-square p <= 0.001 -> near-exclusive allelic isoforms (R-15). Haplotypes
    # are assigned by exact quota (see simulate_flnc), so 0.9 x 20-24 reads leaves 2 minority reads per isoform.
    Gene("asasA", "chr2", "+", [300, 180, 260], [450, 400], [("skip", ("skip", 1), 20, False)], 24, snps=3, hap_bias=0.9, category="asas", ngs_psi={"condA": 0.55, "condB": 0.55}),
    Gene("asasB", "chr2", "-", [280, 170, 250], [430, 410], [("skip", ("skip", 1), 20, True)], 24, snps=3, hap_bias=0.9, category="asas"),
    Gene("asasC", "chr2", "+", [310, 160, 270], [440, 420], [("skip", ("skip", 1), 20, False)], 24, snps=4, hap_bias=0.95, category="asas"),
    Gene("asasD", "chr2", "-", [290, 150, 240], [460, 430], [("skip", ("skip", 1), 20, False)], 24, snps=3, hap_bias=0.9, category="asas"),
    Gene("asasN1", "chr2", "+", [300, 170, 250], [450, 400], [("skip", ("skip", 1), 20, False)], 24, snps=3, hap_bias=0.5, category="asas_null"),
    Gene("asasN2", "chr2", "-", [280, 160, 260], [430, 410], [("skip", ("skip", 1), 20, True)], 24, snps=3, expr_bias=0.85, category="ase",
         note="allele-specific expression without splicing bias: both isoforms 85% hapA (ASAS negative, ASE positive)"),
    Gene("teA", "chr2", "+", [300, 200, 250], [1400, 400], [("skip", ("skip", 1), 3, False)], 6, note="te_repeat_in_intron0"),
    # ---- chr4: short chromosome -----------------------------------------------------------------
    Gene("endA", "chr4", "+", [300, 200, 250], [400, 400], [("skip", ("skip", 1), 4, False)], 8, note="near_chrom_end"),
    Gene("monoD", "chr4", "-", [450], [], [], 5, category="mono_exon", gap=300),
] + _pu_background()

PA_EXTRA = {"irA": ((0, 0.7), (80, 0.3)), "seB": ((0, 0.6), (70, 0.3), (160, 0.1)), "a3A": ((0, 0.7), (100, 0.3))}
for _g in CATALOG:
    if _g.name in PA_EXTRA:
        _g.pa_sites = PA_EXTRA[_g.name]


# ----------------------------------------------------------------------------- structures
@dataclass
class Isoform:
    name: str
    gene: Gene
    exons: list                      # absolute 0-based half-open
    n_flnc: int
    annotated: bool
    op: tuple
    wobble_shift: int = 0
    note: str = ""

    @property
    def chrom(self):
        return self.gene.chrom

    @property
    def strand(self):
        return self.gene.strand

    @property
    def introns(self):
        return [(self.exons[i][1], self.exons[i + 1][0]) for i in range(len(self.exons) - 1)]

    @property
    def junc_chain(self):
        return ";".join(f"{s}-{e}" for s, e in self.introns)


def apply_op(exons, op):
    """Return a new exon list for ``op`` applied to the full structure ``exons``."""
    ex = [list(e) for e in exons]
    kind = op[0]
    if kind == "full":
        pass
    elif kind == "skip":
        lo, hi = op[1], op[-1]
        ex = ex[:lo] + ex[hi + 1:]
    elif kind == "ir":
        i = op[1]
        ex[i][1] = ex[i + 1][1]
        del ex[i + 1]
    elif kind == "end":
        ex[op[1]][1] += op[2]
    elif kind == "start":
        ex[op[1]][0] += op[2]
    elif kind == "complex":            # exon i end +d1 and exon i+1 start -d2
        ex[op[1]][1] += op[2]
        ex[op[1] + 1][0] -= op[3]
    elif kind == "tss":
        ex[0][0] += op[1]
    elif kind == "wobble":             # ("wobble", i, "end"|"start", delta)
        if op[2] == "end":
            ex[op[1]][1] += op[3]
        else:
            ex[op[1]][0] += op[3]
    else:
        raise ValueError(op)
    return [tuple(e) for e in ex]


def layout(catalog):
    """Place every gene on its chromosome; return isoforms, gene spans and chromosome lengths."""
    cursor = {}
    isoforms, spans = [], {}
    for g in catalog:
        start = cursor.get(g.chrom, 1500) + g.gap
        exons, pos = [], start
        for k, length in enumerate(g.exons):
            exons.append((pos, pos + length))
            pos += length + (g.introns[k] if k < len(g.introns) else 0)
        spans[g.name] = (exons[0][0], exons[-1][1])
        cursor[g.chrom] = max(cursor.get(g.chrom, 0), exons[-1][1] + max(o for o, _ in g.pa_sites))
        if g.n_full:
            isoforms.append(Isoform(f"{g.name}.full", g, exons, g.n_full, g.full_annotated, ("full",)))
        for suffix, op, n, annotated in g.variants:
            iso = Isoform(f"{g.name}.{suffix}", g, apply_op(exons, op), n, annotated, op)
            if op[0] == "wobble":
                iso.wobble_shift = op[3]
                iso.note = "junction shifted; NGS supports only the canonical site"
            isoforms.append(iso)
    chrom_len = {c: end + 2000 for c, end in cursor.items()}
    return isoforms, spans, chrom_len


def as_events(isoforms, scale, jfilter):
    """Designed AS events with the detection expectation of the current find_as; events whose inclusion or
    exclusion isoform does not survive the legacy mapping junction-combination filter become ``no:junc-filter``."""
    by_gene = {}
    for iso in isoforms:
        by_gene.setdefault(iso.gene.name, []).append(iso)
    rows = []
    for gname, isos in by_gene.items():
        full = next((i for i in isos if i.op[0] == "full"), None)
        strand = isos[0].strand
        for iso in isos:
            op = iso.op
            if op[0] in ("full", "tss"):
                continue
            if op[0] == "wobble":                                                   # ("wobble", i, side, delta)
                i, side, d = op[1], op[2], op[3]
                site = full.exons[i][1] if side == "end" else full.exons[i][0]
                exp = _min_reads(iso, full, scale, "yes" if abs(d) > AS_OFFSET else "no:offset")
                rows.append(_event(_alt_type(side, strand), iso, full, (min(site, site + d), max(site, site + d)), scale, exp,
                                   f"alignment-wobble isoform ({d:+d} bp); LR-only, NGS confident filter should drop it",
                                   confident="no:ngs-unsupported" if exp == "yes" else exp))
                continue
            g = iso.gene
            if g.n_full == 0 and op[0] == "skip":                                   # MXE pair
                if iso.name.endswith("mxe2"):
                    other = next(i for i in isos if i.name.endswith("mxe1"))
                    rows.append(_event("MXE", other, iso, (other.exons[1][0], iso.exons[1][1]), scale, "no:MXE-unsupported", "mutually exclusive exons"))
                continue
            if op[0] == "skip":
                lo, hi = op[1], op[-1]
                coords = (full.exons[lo][0], full.exons[hi][1])
                micro = any(e - s_ < MICROEXON_MIN for s_, e in full.exons[lo:hi + 1])
                rows.append(_event("SE", full, iso, coords, scale, _min_reads(iso, full, scale, "no:aligner-microexon" if micro else "yes"),
                                   f"skips exon(s) {lo}-{hi}" if hi != lo else f"skips exon {lo}"))
            elif op[0] == "ir":
                i = op[1]
                coords = (full.exons[i][1], full.exons[i + 1][0])
                exp = _min_reads(iso, full, scale, "no:R-2" if len(full.exons) == 2 else "yes")
                rows.append(_event("IR", iso, full, coords, scale, exp, "retention isoform is mono-exon" if exp == "no:R-2" else ""))
            elif op[0] in ("end", "start"):
                i, d = op[1], op[2]
                site = full.exons[i][1] if op[0] == "end" else full.exons[i][0]
                coords = (min(site, site + d), max(site, site + d))
                exp = _min_reads(iso, full, scale, "yes" if abs(d) > AS_OFFSET else "no:offset")
                rows.append(_event(_alt_type(op[0], strand), iso, full, coords, scale, exp, f"site shift {d:+d} bp"))
            elif op[0] == "complex":
                i, d1, d2 = op[1], op[2], op[3]
                e = full.exons[i][1]
                s = full.exons[i + 1][0]
                rows.append(_event(_alt_type("end", strand), iso, full, (e, e + d1), scale, "both:A5SS+A3SS", "complex event: donor and acceptor both differ"))
                rows.append(_event(_alt_type("start", strand), iso, full, (s - d2, s), scale, "both:A5SS+A3SS", "complex event: donor and acceptor both differ"))
    for ev in rows:
        failing = [n for n in (ev["inc_isoforms"], ev["exc_isoforms"]) if jfilter.get(n, "pass") != "pass"]
        if failing:
            ev["expected_detection"] = ev["expected_confident"] = "no:junc-filter"
            ev["note"] = (ev["note"] + "; " if ev["note"] else "") + ", ".join(f"{n} {jfilter[n]}" for n in failing) + \
                         " (legacy mapping junction-combination filter removes its reads)"
    return rows


def _n_reads(iso, scale):
    return max(1, round(iso.n_flnc * scale)) if iso.n_flnc else 0


def _junction_evidence(isoforms, scale):
    """What legacy ``filterByJunc`` sees besides NGS: annotated junctions, FLNC reads per junction and per chain."""
    chain_reads, flnc_junc, annotated = {}, {}, set()
    for iso in isoforms:
        n = _n_reads(iso, scale)
        chain_reads[(iso.chrom, iso.junc_chain)] = chain_reads.get((iso.chrom, iso.junc_chain), 0) + n
        for s_, e in iso.introns:
            flnc_junc[(iso.chrom, s_, e)] = flnc_junc.get((iso.chrom, s_, e), 0) + n
            if iso.annotated:
                annotated.add((iso.chrom, s_, e))
    return annotated, flnc_junc, chain_reads


def ngs_reads(junction_reads, key, libs=None):
    """Designed NGS reads over ``libs`` (default: every library — the committed config merges all four into one
    hisat2 run, so that is what the pipeline's junction filter sees)."""
    counts = junction_reads.get(key, {})
    return sum(counts.get(l, 0) for l in (libs or counts))


def junction_filter(isoforms, junction_reads, scale, libs=None):
    """Legacy ``mapping.filterByJunc`` outcome per isoform. Mono-exon reads always pass; multi-exon reads need
    >= JUNC_COMB_SUPPORT reads with the same junction chain and either >= 4 introns or >= 1/3 consensus junctions
    (annotated, carried by >= 3 FLNC reads, or > NGS_JUNC_MIN_READS NGS reads). Returns (status per isoform,
    annotated junctions, FLNC reads per junction)."""
    annotated, flnc_junc, chain_reads = _junction_evidence(isoforms, scale)

    def consensus(key):
        return key in annotated or flnc_junc[key] >= 3 or ngs_reads(junction_reads, key, libs) > NGS_JUNC_MIN_READS

    status = {}
    for iso in isoforms:
        keys = [(iso.chrom, s_, e) for s_, e in iso.introns]
        if not keys or not _n_reads(iso, scale):
            status[iso.name] = "pass"
        elif chain_reads[(iso.chrom, iso.junc_chain)] < JUNC_COMB_SUPPORT:
            status[iso.name] = "fail:chain-support"
        elif len(keys) >= 4 or sum(map(consensus, keys)) / len(keys) >= 1 / 3:
            status[iso.name] = "pass"
        else:
            status[iso.name] = "fail:consensus"
    return status, annotated, flnc_junc


def palen_detectable(inc_lengths, exc_lengths, min_reads=5):
    """The published palen_as rule on the simulated PASS tail lengths of an AS isoform pair: >= ``min_reads``
    reads each side (the acceptance test runs with -pa_sup 5), Student t-test AND Kruskal-Wallis both p <= 0.001
    (palen_as.getPalenAS). Small groups cannot reach 0.001 with Kruskal (n=5 vs 5: minimum p = 0.0079), R-16."""
    from scipy import stats
    if len(inc_lengths) < min_reads or len(exc_lengths) < min_reads:
        return 0
    return int(stats.ttest_ind(inc_lengths, exc_lengths)[1] <= 0.001 and stats.kruskal(inc_lengths, exc_lengths)[1] <= 0.001)


def asas_significant(inc, exc):
    """The published asas criterion on designed per-haplotype read counts ([hapA, hapB] of the inclusion and
    exclusion isoforms): each isoform >= 5 reads on one haplotype and <= 3 on the other, opposite majority
    haplotypes, chi-square p <= 0.001 (asas.relationshipBetweenAlleleSpecifiAndAS)."""
    from scipy.stats import chi2_contingency
    shape = all(max(c) >= 5 and min(c) <= 3 for c in (inc, exc)) and inc.index(max(inc)) != exc.index(max(exc))
    return int(shape and chi2_contingency([[inc[0], exc[0]], [inc[1], exc[1]]])[1] <= 0.001)


def _min_reads(a, b, scale, expected):
    """find_as ignores isoforms with < LR_MIN_READS supporting reads, so such events cannot be called."""
    return "no:min-reads" if min(a.n_flnc, b.n_flnc) * scale < LR_MIN_READS else expected


def _alt_type(side, strand):
    """A change at an exon's genomic END is a donor change on '+' (A5SS) and an acceptor change on '-' (A3SS)."""
    return ("A5SS" if side == "end" else "A3SS") if strand == "+" else ("A3SS" if side == "end" else "A5SS")


def _event(typ, inc, exc, coords, scale, expected, note, confident=None):
    n_inc, n_exc = inc.n_flnc * scale, exc.n_flnc * scale
    return {"event_id": f"{inc.gene.name}:{typ}:{coords[0]}-{coords[1]}", "type": typ, "gene": inc.gene.name,
            "chrom": inc.chrom, "strand": inc.strand, "start": coords[0], "end": coords[1],
            "inc_isoforms": inc.name, "exc_isoforms": exc.name,
            "lr_psi": round(n_inc / (n_inc + n_exc), 3) if n_inc + n_exc else "",
            "expected_detection": expected,
            "expected_confident": confident or ("yes" if expected in ("yes", "both:A5SS+A3SS") else expected),
            "note": note}


# ----------------------------------------------------------------------------- genome
def build_genome(rng, isoforms, spans, chrom_len, te_unit):
    """Random sequence with splice motifs planted for every isoform intron; returns (seq dict, features)."""
    seq = {c: [rng.choice("ACGT") if rng.random() > 0.05 else rng.choice("AT") for _ in range(n)] for c, n in chrom_len.items()}
    features = []
    for iso in isoforms:
        g = iso.gene
        for k, (s, e) in enumerate(iso.introns):
            motif = g.motif_by_intron.get(k, g.motif) if iso.op[0] == "full" else g.motif
            _plant(seq[iso.chrom], s, e, iso.strand, motif)
            if motif == "GA-AG":
                features.append((iso.chrom, s, e, "noncanonical_splice", f"{iso.name} intron {k} {motif}"))
    # poly(A) signal upstream of every PA site makes nanopolish-style/PAS motif checks meaningful
    for iso in isoforms:
        if iso.op[0] != "full":
            continue
        for off, _ in iso.gene.pa_sites:
            site = iso.exons[-1][1] + off if iso.strand == "+" else iso.exons[0][0] - off
            pos = site - 22 if iso.strand == "+" else site + 16
            if 0 <= pos < len(seq[iso.chrom]) - 6:
                seq[iso.chrom][pos:pos + 6] = list("AATAAA" if iso.strand == "+" else "TTTATT")
    # TE-like repeat: 3 intergenic copies on chr2 + 1 on chr1 + one inside teA intron 0
    te_pos = [("chr2", p) for p in (200, 500 + len(te_unit) + 50, 900 + 2 * len(te_unit))] + [("chr1", 300)]
    tea = next(i for i in isoforms if i.name == "teA.full")
    te_pos.append(("chr2", tea.exons[0][1] + 200))
    for chrom, p in te_pos:
        seq[chrom][p:p + len(te_unit)] = list(te_unit.lower())
        features.append((chrom, p, p + len(te_unit), "te_repeat", "soft-masked repeat copy"))
    return seq, features


def _plant(seq, start, end, strand, motif):
    donor, acceptor = MOTIFS[motif]
    if strand == "+":
        seq[start:start + 2], seq[end - 2:end] = list(donor), list(acceptor)
    else:
        seq[start:start + 2], seq[end - 2:end] = list(revcomp(acceptor)), list(revcomp(donor))


def soft_mask(rng, seq, catalog, spans, features, target=0.30):
    """Lowercase ``softmask`` gene spans, TE copies and random intergenic windows (up to ``target``)."""
    masked = {c: [False] * len(s) for c, s in seq.items()}
    gene_bp = {c: [False] * len(s) for c, s in seq.items()}
    for chrom, s, e, feature, _ in features:
        if feature == "te_repeat":
            masked[chrom][s:e] = [True] * (e - s)
    for g in catalog:
        s, e = spans[g.name]
        gene_bp[g.chrom][s:e] = [True] * (e - s)
        if g.softmask:
            masked[g.chrom][s:e] = [True] * (e - s)
            features.append((g.chrom, s, e, "softmasked", f"{g.name} entire span incl. splice sites"))
    for c, s in seq.items():
        n = len(s)
        tries = 0
        while sum(masked[c]) < target * n and tries < 3000:
            tries += 1
            w = rng.randint(200, 600)
            p = rng.randint(0, n - w)
            if any(gene_bp[c][p:p + w]) or any(masked[c][p:p + w]):
                continue
            masked[c][p:p + w] = [True] * w
            features.append((c, p, p + w, "softmasked", "intergenic"))
    for c in seq:
        seq[c] = "".join(ch.lower() if m else ch.upper() for ch, m in zip(seq[c], masked[c]))
    return seq


def plant_snps(rng, seq_upper, isoforms):
    """Heterozygous exonic SNPs for diploid genes -> (hapB base map, vcf rows)."""
    hapb, vcf = {}, []
    for iso in isoforms:
        g = iso.gene
        if not g.snps or iso.op[0] != "full":
            continue
        candidates = [p for s, e in iso.exons for p in range(s + 15, e - 15)]
        chosen = []
        for p in rng.sample(candidates, len(candidates)):
            if all(abs(p - q) >= 100 for q in chosen):
                chosen.append(p)
            if len(chosen) == g.snps:
                break
        for p in sorted(chosen):
            ref = seq_upper[iso.chrom][p]
            alt = rng.choice([b for b in "ACGT" if b != ref])
            hapb[(iso.chrom, p)] = alt
            vcf.append((iso.chrom, p + 1, f"{g.name}_snp{len(vcf) + 1}", ref, alt, g.name))
    return hapb, vcf


# ----------------------------------------------------------------------------- reads
def transcript(seq_upper, iso, exons, hap, hapb):
    parts = []
    for s, e in exons:
        block = list(seq_upper[iso.chrom][s:e])
        if hap == "hapB":
            for p in range(s, e):
                if (iso.chrom, p) in hapb:
                    block[p - s] = hapb[(iso.chrom, p)]
        parts.append("".join(block))
    tx = "".join(parts)
    return tx if iso.strand == "+" else revcomp(tx)


def mutate(rng, s, rate):
    if rate <= 0:
        return s, 0
    out, n = list(s), 0
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
            n += 1
    return "".join(out), n


def simulate_flnc(rng, seq_upper, isoforms, hapb, scale, error_rate):
    """FLNC reads (poly(A)-trimmed, as `isoseq refine --require-polya` emits them) with TES jitter,
    5' truncation, haplotype origin and PA-site choice; the designed tail length is kept in the truth and
    only materialised in the CCS-like reads (ccs.bam / ccs.fastq)."""
    reads = []
    for iso in isoforms:
        g = iso.gene
        n = max(1, round(iso.n_flnc * scale)) if iso.n_flnc else 0
        mean, sd = g.palen.get(iso.name.split(".")[-1], DEFAULT_PALEN)
        haps = ["hapA"] * n
        if g.snps:                                   # exact haplotype quota: the designed bias is the truth, not a draw
            n_a = round(n * _p_hap_a(g, iso))
            haps = ["hapA"] * n_a + ["hapB"] * (n - n_a)
            rng.shuffle(haps)
        for k in range(n):
            site_off = rng.choices([o for o, _ in g.pa_sites], [f for _, f in g.pa_sites])[0]
            jitter = int(round(rng.gauss(0, 3)))
            trunc = rng.randint(0, 50) if rng.random() < 0.2 else 0
            exons = [list(e) for e in iso.exons]
            if iso.strand == "+":
                exons[-1][1] += site_off + jitter
                exons[0][0] += trunc
            else:
                exons[0][0] -= site_off + jitter
                exons[-1][1] -= trunc
            hap = haps[k]
            tx = transcript(seq_upper, iso, [tuple(e) for e in exons], hap, hapb)
            tx, nerr = (mutate(rng, tx, error_rate) if rng.random() < 0.5 else (tx, 0))
            palen = max(10, int(round(rng.gauss(mean, sd))))
            reads.append({"read_id": f"{iso.name}/ccs/{k}", "isoform": iso.name, "gene": g.name, "haplotype": hap,
                          "pa_site_offset": site_off, "polya_len": palen, "tes_jitter": jitter, "five_prime_trunc": trunc,
                          "wobble_shift": iso.wobble_shift, "n_errors": nerr, "seq": tx})
    rng.shuffle(reads)
    return reads


def simulate_ngs(rng, seq_upper, isoforms, hapb, scale, stranded, read_len=150):
    """Paired-end libraries per condition x replicate; designed PSI differences via ``Gene.ngs_psi``.

    Also returns the designed read support of every junction (reads with >= NGS_ANCHOR bp on both exons, the
    quantity ``regtools junctions extract -a 5`` reports). Junctions whose legacy junction-filter consensus hinges
    on NGS (unannotated, < 3 FLNC reads) should keep their total support out of NGS_BAND, so the "> 10 reads"
    decision never depends on a read the aligner may add or lose (warned here, asserted on the committed data)."""
    annotated, flnc_junc, chain_reads = _junction_evidence(isoforms, scale)
    by_gene = {}
    for iso in isoforms:
        by_gene.setdefault(iso.gene.name, []).append(iso)
    support = {iso.name: {} for iso in isoforms}
    junction_reads, libs = {}, {}
    for cond in CONDITIONS:
        for rep in REPLICATES:
            lib = f"{cond}_{rep}"
            pairs = []
            for gname, isos in by_gene.items():
                g = isos[0].gene
                weights = {iso.name: iso.n_flnc * 3.0 * scale for iso in isos}
                if g.ngs_psi and len(isos) == 2:                       # one designed event per diff gene
                    full = next(i for i in isos if i.op[0] == "full")
                    alt = next(i for i in isos if i is not full)
                    inc, exc = (alt, full) if alt.op[0] in ("ir", "end", "start", "complex") else (full, alt)
                    total = weights[inc.name] + weights[exc.name]
                    weights[inc.name], weights[exc.name] = total * g.ngs_psi[cond], total * (1 - g.ngs_psi[cond])
                for iso in isos:
                    if iso.wobble_shift:
                        weights[iso.name] = 0.0                        # NGS never supports shifted junctions
                    weights[iso.name] *= g.ngs_weight.get(iso.name.split(".")[-1], 1.0)
                    n = round(weights[iso.name] * rng.lognormvariate(0, 0.15))
                    support[iso.name][lib] = n
                    pairs.extend(_fragments(rng, seq_upper, iso, hapb, n, read_len, stranded, lib, junction_reads))
            rng.shuffle(pairs)
            libs[lib] = pairs
    # junctions whose NGS support decides the filter: unannotated, < 3 FLNC reads, on a chain that passes the
    # chain-support test (a 1-read chain is dropped before consensus is even looked at)
    deciding = {(iso.chrom, s_, e) for iso in isoforms if chain_reads[(iso.chrom, iso.junc_chain)] >= JUNC_COMB_SUPPORT
                for s_, e in iso.introns if (iso.chrom, s_, e) not in annotated and flnc_junc[(iso.chrom, s_, e)] < 3}
    edgy = [k for k in sorted(deciding) if NGS_BAND[0] < ngs_reads(junction_reads, k) < NGS_BAND[1]]
    if edgy:
        warnings.warn(f"designed NGS support of filter-deciding junctions inside {NGS_BAND}: {edgy}")
    return libs, support, junction_reads


def _junction_positions(iso):
    """Transcript-coordinate position of each intron of ``iso`` (same order as ``iso.introns``)."""
    length = sum(e - s_ for s_, e in iso.exons)
    out, cum = [], 0
    for s_, e in iso.exons[:-1]:
        cum += e - s_
        out.append(cum if iso.strand == "+" else length - cum)
    return out


def _count_spanning(iso, start, end, read_len, lib, junction_reads):
    """Count the mates of fragment [start, end) that cover each junction of ``iso`` with >= NGS_ANCHOR bp."""
    for (s_, e), p in zip(iso.introns, _junction_positions(iso)):
        for rs, re_ in ((start, start + read_len), (end - read_len, end)):
            if rs + NGS_ANCHOR <= p <= re_ - NGS_ANCHOR:
                reads = junction_reads.setdefault((iso.chrom, s_, e), {})
                reads[lib] = reads.get(lib, 0) + 1


def _p_hap_a(g, iso):
    """P(read of ``iso`` comes from haplotype A): expression bias (all isoforms), splicing bias (inclusion vs
    exclusion isoform) or 0.5 for the null genes."""
    if g.expr_bias is not None:
        return g.expr_bias
    if g.category == "asas_null":
        return 0.5
    return g.hap_bias if iso.op[0] != "full" else 1 - g.hap_bias


def _haplotype(rng, iso):
    g = iso.gene
    if not g.snps:
        return "hapA"
    return "hapA" if rng.random() < _p_hap_a(g, iso) else "hapB"


def _pair(rng, tx, start, end, read_len, stranded, name):
    frag = tx[start:end]
    r1, r2 = frag[:read_len], revcomp(frag)[:read_len]
    if stranded or rng.random() < 0.5:          # fr-firststrand: R1 from the reverse strand
        r1, r2 = r2, r1
    return (name, r1, r2)


def _fragments(rng, seq_upper, iso, hapb, n, read_len, stranded, lib, junction_reads):
    out = []
    for k in range(n):
        tx = transcript(seq_upper, iso, iso.exons, _haplotype(rng, iso), hapb)
        if len(tx) < read_len + 20:
            continue
        frag_len = min(len(tx), rng.randint(250, 450))
        start = rng.randint(0, len(tx) - frag_len)
        _count_spanning(iso, start, start + frag_len, read_len, lib, junction_reads)
        out.append(_pair(rng, tx, start, start + frag_len, read_len, stranded, f"{lib}_{iso.name}_{k}"))
    return out


# ----------------------------------------------------------------------------- writers
def write_fasta(path, seq, width=60):
    with open(path, "w") as f:
        for c, s in seq.items():
            f.write(f">{c}\n")
            for i in range(0, len(s), width):
                f.write(s[i:i + width] + "\n")


def write_gtf(path, isoforms, annotated_only):
    with open(path, "w") as f:
        for iso in isoforms:
            if annotated_only and not iso.annotated:
                continue
            for s, e in iso.exons:
                f.write(f'{iso.chrom}\tsyn\texon\t{s + 1}\t{e}\t.\t{iso.strand}\t.\tgene_id "{iso.gene.name}"; transcript_id "{iso.name}";\n')


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in columns})


def write_gz(path, text):
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
        gz.write(text.encode())
    with open(path, "wb") as f:
        f.write(buf.getvalue())


def write_ccs(rng, reads, bam_path, fastq_path, truth):
    import pysam
    header = {"HD": {"VN": "1.6", "SO": "unknown"},
              "RG": [{"ID": "synv2", "PL": "PACBIO", "PM": "SEQUEL", "PU": MOVIE, "SM": "synthetic", "DS": "READTYPE=CCS"}]}
    fq = []
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for zmw, (r, t) in enumerate(zip(reads, truth), start=1):
            p5 = p3 = True
            if rng.random() < 0.05:
                p5 = rng.random() < 0.5
                p3 = not p5
            # Iso-Seq layout of a sense-strand CCS read: 5p primer, insert, poly(A), then the 3p primer as listed in
            # primer.fa (lima scores a 5p/3p pair; appending its reverse complement would make a 5p--5p read)
            s = (PRIMER_5P if p5 else "") + r["seq"] + "A" * r["polya_len"] + (PRIMER_3P if p3 else "")
            if rng.random() < 0.5:
                s = revcomp(s)
            t["has_primer_5p"], t["has_primer_3p"] = int(p5), int(p3)
            t["zmw"], t["flnc_recovered"] = zmw, int(p5 and p3 and r["polya_len"] >= ISOSEQ_MIN_POLYA)
            name = f"{MOVIE}/{zmw}/ccs"
            a = pysam.AlignedSegment(bam.header)
            a.query_name, a.flag, a.query_sequence = name, 4, s
            a.query_qualities = pysam.qualitystring_to_array("~" * len(s))
            a.set_tags([("np", 12, "i"), ("rq", 0.995, "f"), ("zm", zmw, "i"), ("RG", "synv2", "Z")])
            bam.write(a)
            fq.append(f"@{name}\n{s}\n+\n{'~' * len(s)}\n")
    with open(fastq_path, "w") as f:
        f.writelines(fq)



# ONT direct RNA (2026-10-15): the FLNC molecules read again the way a nanopore basecaller reports them
ONT_CHEMISTRY = {
    # per-base error rates (RNA002 ~9 %, RNA004 ~3 %), share of the poly(A) tail the basecaller emits, quality char
    "rna002": {"sub": 0.03, "ins": 0.02, "del": 0.04, "tail": (0.5, 1.0), "qual": "+"},
    "rna004": {"sub": 0.01, "ins": 0.005, "del": 0.015, "tail": (0.5, 1.0), "qual": "5"},
}


def mutate_ont(rng, s, sub, ins, dele):
    """Substitutions, single-base deletions and insertions at per-base rates; returns (read, number of errors)."""
    out, n = [], 0
    for b in s:
        r = rng.random()
        if r < dele:
            n += 1
            continue
        if r < dele + sub:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
            n += 1
        else:
            out.append(b)
        if rng.random() < ins:
            out.append(rng.choice("ACGT"))
            n += 1
    return "".join(out), n


def simulate_ont(rng, seq_upper, isoforms, hapb, scale, chem):
    """Direct-RNA reads of the same design as the FLNC reads: sense orientation, the 5' end short by the motor-protein
    stall (8-20 nt), the chemistry's error profile with indels, part of the poly(A) tail basecalled, UUID read ids."""
    p = ONT_CHEMISTRY[chem]
    reads = []
    for t in simulate_flnc(rng, seq_upper, isoforms, hapb, scale, 0.0):
        stall = rng.randint(8, 20)
        seq, nerr = mutate_ont(rng, t["seq"][stall:], p["sub"], p["ins"], p["del"])
        emitted = round(t["polya_len"] * rng.uniform(*p["tail"]))
        reads.append({"read_id": str(uuid.UUID(int=rng.getrandbits(128))), "isoform": t["isoform"], "gene": t["gene"],
                      "haplotype": t["haplotype"], "pa_site_offset": t["pa_site_offset"], "polya_len": t["polya_len"],
                      "five_prime_trunc": t["five_prime_trunc"] + stall, "n_errors": nerr, "emitted_tail": emitted,
                      "template": t["seq"], "seq": seq + "A" * emitted})
    return reads


def write_ont(rng, reads, out_dir, chem):
    """RNA002: guppy-style FASTQ (U for T). RNA004: dorado output -- a gzipped FASTQ and an unaligned BAM (DNA alphabet)
    whose ``pt:i`` tag is the estimated tail length (-1 when no tail was found)."""
    import pysam
    os.makedirs(out_dir, exist_ok=True)
    q = ONT_CHEMISTRY[chem]["qual"]
    fastq = "".join(f"@{r['read_id']}\n{r['seq'].replace('T', 'U')}\n+\n{q * len(r['seq'])}\n" for r in reads)
    if chem == "rna002":
        with open(os.path.join(out_dir, "reads.fastq"), "w") as f:
            f.write(fastq)
    else:
        write_gz(os.path.join(out_dir, "reads.fastq.gz"), fastq)
        header = {"HD": {"VN": "1.6", "SO": "unknown"},
                  "RG": [{"ID": "synv2_rna004", "PL": "ONT", "DS": "basecall_model=rna004_130bps_sup@v5.1.0 runid=synv2"}]}
        with pysam.AlignmentFile(os.path.join(out_dir, "dorado.bam"), "wb", header=header) as bam:
            for r in reads:
                r["dorado_pt"] = max(1, round(r["polya_len"] + rng.gauss(0, 3))) if rng.random() < 0.95 else -1
                a = pysam.AlignedSegment(bam.header)
                a.query_name, a.flag, a.query_sequence = r["read_id"], 4, r["seq"]
                a.query_qualities = pysam.qualitystring_to_array(q * len(r["seq"]))
                a.set_tags([("pt", r["dorado_pt"], "i"), ("RG", "synv2_rna004", "Z")])
                bam.write(a)
    columns = ["read_id", "isoform", "gene", "haplotype", "pa_site_offset", "polya_len", "five_prime_trunc", "n_errors", "emitted_tail"]
    write_tsv(os.path.join(out_dir, "read_truth.tsv"), reads, columns + (["dorado_pt"] if chem == "rna004" else []))

def write_configs(out, scale):
    yaml_text = f"""# iFLAS synthetic v2 config (paths relative to this file). Reference aux files are derived by
# `iflas init-ref -f genome.fa --gtf genes.gtf -o ref/`; tests/test_synthetic_v2.py fills them in.
reference:
  - ref_strain: SYN
    ref_genome: genome.fa
    ref_gtf: genes.gtf
ccs:
  min_pass: 2
  min_predicted_accuracy: 0.9
minimap2:
  max_intron_length: 50000
  junc_comb_support: 2
collapse:
  min_identity: 0.9
  min_coverage: 0.9
pu_filter:
  filter_score: 0.5
  pos_fl_coverage: 2
  pos_min_junc_rpkm: 0.05
options:
  threads: 4
  memory: 8000M
  merge_data_from_same_strain: true
dirs:
  out_dir: iflas_out
"""
    with open(os.path.join(out, "config.yaml"), "w") as f:
        f.write(yaml_text)
    cols = ["sample_name", "project_name", "ref_strain", "strain", "condition", "replicate", "tgs_plat", "strategy",
            "data_location", "data_processed_location", "primer", "ngs_left_reads", "ngs_right_reads", "ngs_reads_paired",
            "ngs_reads_length", "use_fmlrc2"]
    rows = []
    for cond in CONDITIONS:
        for rep in REPLICATES:
            tgs = cond == "condA" and rep == "rep1"
            rows.append({"sample_name": f"{cond}_{rep}", "project_name": "synv2", "ref_strain": "SYN", "strain": "SYN",
                         "condition": cond, "replicate": rep, "tgs_plat": "pacbio" if tgs else "None",
                         "strategy": "sequel" if tgs else "None", "data_location": "None",
                         "data_processed_location": "flnc.fa" if tgs else "None", "primer": "primer.fa" if tgs else "None",
                         "ngs_left_reads": f"ngs/{cond}_{rep}_R1.fastq.gz", "ngs_right_reads": f"ngs/{cond}_{rep}_R2.fastq.gz",
                         "ngs_reads_paired": "paired", "ngs_reads_length": "150", "use_fmlrc2": "false"})
    write_tsv(os.path.join(out, "samples.tsv"), rows, cols)
    # collapse/hqIsoFilter values below are the legacy Config.py DEFAULTS (identical to the modern pydantic defaults)
    legacy = """[refSection]
section_type = refSection
ref_strain = SYN
ref_genome = {V2}/genome.fa
ref_size = {V2}/ref/genome.len.txt
ref_gtf = {V2}/genes.gtf
ref_gff = None
ref_gpe = {V2}/ref/genome.gpe
ref_bed = {V2}/ref/genome.bed
ref_mm2_index = {V2}/ref/genome.mm2
hisat2_index = {V2}/ref/genome_hisat2

[condA_rep1]
section_type = dataCfg
project_name = synv2
sample_name = condA_rep1
ref_strain = SYN
strain = SYN
condition = condA
tgs_plat = pacbio
strategy = sequel
data_location = None
primer = {V2}/primer.fa
data_processed_location = {V2}/flnc.fa
polya_location = None
ngs_left_reads = {V2}/ngs/condA_rep1_R1.fastq.gz
ngs_right_reads = {V2}/ngs/condA_rep1_R2.fastq.gz
ngs_reads_paired = paired
ngs_reads_length = 150
ngs_junctions = None
ngs_library_strand = 0
use_fmlrc2 = False
single_run_threads = 4

[ccsCalling]
section_type = ccsCalling
min_read_score = 0.75
min_subread_length = 50
min_ccs_length = 50
min_predicted_accuracy = 0.9
min_pass = 2

[minimap2]
section_type = minimap2
mm2_index = None
max_intron_length = 50000

[collapse]
section_type = collapse
min_identity = 0.9
min_coverage = 0.9
max_5_diff = 1000
max_3_diff = 500
fl_coverage = 2
max_fuzzy_junction = 5
dun_merge_5_shorter = True

[pu_filter]
section_type = pu_filter
filter_score = 0.5
draw_auc = False
pos_fl_coverage = 2
pos_min_junc_rpkm = 0.05

[optionTools]
section_type = optionTools
threads = 4
memory = 8G
merge_data_from_same_strain = False
gene2go = None

[dirs]
section_type = dirs
out_dir = {V2}/legacy_out
tmp_dir = {V2}/legacy_tmp
"""
    with open(os.path.join(out, "legacy.cfg.template"), "w") as f:
        f.write(legacy)


# ----------------------------------------------------------------------------- main
def generate(out, seed=20261014, scale=1.0, error_rate=0.003, stranded=False):
    rng = random.Random(seed)
    os.makedirs(os.path.join(out, "truth"), exist_ok=True)
    os.makedirs(os.path.join(out, "ngs"), exist_ok=True)
    isoforms, spans, chrom_len = layout(CATALOG)
    te_unit = "".join(rng.choice("ACGT") for _ in range(300))
    seq_list, features = build_genome(rng, isoforms, spans, chrom_len, te_unit)
    seq_upper = {c: "".join(s).upper() for c, s in seq_list.items()}
    hapb, vcf = plant_snps(rng, seq_upper, isoforms)
    masked = soft_mask(rng, {c: list(s) for c, s in seq_upper.items()}, CATALOG, spans, features)
    write_fasta(os.path.join(out, "genome.fa"), masked)
    write_fasta(os.path.join(out, "genome.unmasked.fa"), seq_upper)
    write_gtf(os.path.join(out, "genes.gtf"), isoforms, annotated_only=True)
    write_gtf(os.path.join(out, "truth", "genes.full.gtf"), isoforms, annotated_only=False)
    with open(os.path.join(out, "primer.fa"), "w") as f:
        f.write(f">5p\n{PRIMER_5P}\n>3p\n{PRIMER_3P}\n")

    reads = simulate_flnc(rng, seq_upper, isoforms, hapb, scale, error_rate)
    with open(os.path.join(out, "flnc.fa"), "w") as fa:
        for r in reads:
            fa.write(f">{r['read_id']}\n{r['seq']}\n")
    truth_reads = [{k: v for k, v in r.items() if k != "seq"} for r in reads]
    write_ccs(rng, reads, os.path.join(out, "ccs.bam"), os.path.join(out, "ccs.fastq"), truth_reads)
    with open(os.path.join(out, "polya_results.tsv"), "w") as f:
        f.write("readname\tcontig\tposition\tleader_start\tadapter_start\tpolya_start\ttranscript_start\tread_rate\tpolya_length\tqc_tag\n")
        for r, t in zip(reads, truth_reads):
            tag = "PASS" if rng.random() < 0.9 else rng.choice(["NOPASS", "ADAPTER", "SUFFCLIP"])
            t["polya_qc"] = tag
            f.write(f"{r['read_id']}\t{t['gene']}\t0\t0\t10\t20\t{20 + r['polya_len']}\t80.0\t{r['polya_len']:.2f}\t{tag}\n")
    libs, support, junction_reads = simulate_ngs(rng, seq_upper, isoforms, hapb, scale, stranded)
    jfilter, annotated_junc, flnc_junc = junction_filter(isoforms, junction_reads, scale)
    for lib, pairs in libs.items():
        write_gz(os.path.join(out, "ngs", f"{lib}_R1.fastq.gz"), "".join(f"@{n}/1\n{a}\n+\n{'I' * len(a)}\n" for n, a, _ in pairs))
        write_gz(os.path.join(out, "ngs", f"{lib}_R2.fastq.gz"), "".join(f"@{n}/2\n{b}\n+\n{'I' * len(b)}\n" for n, _, b in pairs))

    # ---- truth tables ----
    T = lambda name: os.path.join(out, "truth", name)
    write_tsv(T("genome_features.tsv"), [dict(zip(("chrom", "start", "end", "feature", "note"), f)) for f in sorted(features)],
              ["chrom", "start", "end", "feature", "note"])
    write_tsv(T("genes.tsv"), [{"gene": g.name, "chrom": g.chrom, "strand": g.strand, "start": spans[g.name][0], "end": spans[g.name][1],
                                "n_isoforms": sum(1 for i in isoforms if i.gene is g), "category": g.category,
                                "softmasked": int(g.softmask), "motif": g.motif, "note": g.note} for g in CATALOG],
              ["gene", "chrom", "strand", "start", "end", "n_isoforms", "category", "softmasked", "motif", "note"])
    write_tsv(T("isoforms.tsv"), [{"isoform": i.name, "gene": i.gene.name, "chrom": i.chrom, "strand": i.strand,
                                   "exons": ";".join(f"{s}-{e}" for s, e in i.exons), "junc_chain": i.junc_chain,
                                   "n_flnc": max(1, round(i.n_flnc * scale)) if i.n_flnc else 0, "annotated": int(i.annotated),
                                   "junction_filter": jfilter[i.name],
                                   "category": i.gene.category, "op": ":".join(map(str, i.op)), "wobble_shift": i.wobble_shift,
                                   "palen_mean": i.gene.palen.get(i.name.split(".")[-1], DEFAULT_PALEN)[0], "note": i.note} for i in isoforms],
              ["isoform", "gene", "chrom", "strand", "exons", "junc_chain", "n_flnc", "annotated", "junction_filter", "category", "op",
               "wobble_shift", "palen_mean", "note"])
    events = as_events(isoforms, scale, jfilter)
    write_tsv(T("as_events.tsv"), events, ["event_id", "type", "gene", "chrom", "strand", "start", "end", "inc_isoforms", "exc_isoforms",
                                           "lr_psi", "expected_detection", "expected_confident", "note"])
    libnames = [f"{c}_{r}" for c in CONDITIONS for r in REPLICATES]
    write_tsv(T("isoform_support.tsv"), [{"isoform": i.name, "n_flnc": max(1, round(i.n_flnc * scale)) if i.n_flnc else 0,
                                          **{f"ngs_{l}": support[i.name].get(l, 0) for l in libnames}} for i in isoforms],
              ["isoform", "n_flnc"] + [f"ngs_{l}" for l in libnames])
    carriers = {}
    for i in isoforms:
        for s_, e in i.introns:
            carriers.setdefault((i.chrom, s_, e), []).append(i)
    write_tsv(T("junction_support.tsv"),
              [{"chrom": c, "start": s_, "end": e, "strand": isos[0].strand, "isoforms": ";".join(i.name for i in isos),
                "annotated": int((c, s_, e) in annotated_junc), "flnc_reads": flnc_junc[(c, s_, e)],
                **{f"ngs_reads_{l}": junction_reads.get((c, s_, e), {}).get(l, 0) for l in libnames},
                "ngs_reads_total": ngs_reads(junction_reads, (c, s_, e)),
                "ngs_detectable": int(min(x[1] - x[0] for i in isos for x in i.exons if x[1] == s_ or x[0] == e) >= MICROEXON_MIN),
                "consensus": int((c, s_, e) in annotated_junc or flnc_junc[(c, s_, e)] >= 3
                                 or ngs_reads(junction_reads, (c, s_, e)) > NGS_JUNC_MIN_READS)}
               for (c, s_, e), isos in sorted(carriers.items())],
              ["chrom", "start", "end", "strand", "isoforms", "annotated", "flnc_reads"] + [f"ngs_reads_{l}" for l in libnames]
              + ["ngs_reads_total", "ngs_detectable", "consensus"])
    write_tsv(T("read_truth.tsv"), truth_reads, ["read_id", "isoform", "gene", "haplotype", "pa_site_offset", "polya_len", "polya_qc",
                                                 "tes_jitter", "five_prime_trunc", "wobble_shift", "n_errors", "has_primer_5p", "has_primer_3p",
                                                 "zmw", "flnc_recovered"])
    diff = []
    for ev in events:
        g = next(x for x in CATALOG if x.name == ev["gene"])
        if g.ngs_psi:
            a, b = g.ngs_psi["condA"], g.ngs_psi["condB"]
            delta = round(a - b, 3)                     # rounded first: 0.6 - 0.4 is 0.19999... in binary
            diff.append({"event_id": ev["event_id"], "gene": g.name, "condA_psi": a, "condB_psi": b, "delta": delta,
                         "expected_significant": int(abs(delta) >= 0.2)})
    write_tsv(T("diff_psi.tsv"), diff, ["event_id", "gene", "condA_psi", "condB_psi", "delta", "expected_significant"])
    hap_counts = {}
    for r in truth_reads:
        hap_counts[(r["isoform"], r["haplotype"])] = hap_counts.get((r["isoform"], r["haplotype"]), 0) + 1
    allele_rows = []
    for ev in events:
        g = next((x for x in CATALOG if x.name == ev["gene"] and x.snps), None)
        if g is None:
            continue
        inc = [hap_counts.get((ev["inc_isoforms"], h), 0) for h in ("hapA", "hapB")]
        exc = [hap_counts.get((ev["exc_isoforms"], h), 0) for h in ("hapA", "hapB")]
        allele_rows.append({"gene": g.name, "event_id": ev["event_id"], "hapA_inc_frac": g.hap_bias, "hapB_inc_frac": round(1 - g.hap_bias, 2),
                            "inc_hapA_reads": inc[0], "inc_hapB_reads": inc[1], "exc_hapA_reads": exc[0], "exc_hapB_reads": exc[1],
                            "expected_significant": asas_significant(inc, exc)})
    write_tsv(T("allele_bias.tsv"), allele_rows, ["gene", "event_id", "hapA_inc_frac", "hapB_inc_frac", "inc_hapA_reads", "inc_hapB_reads",
                                                   "exc_hapA_reads", "exc_hapB_reads", "expected_significant"])
    ase_rows = []
    for g in CATALOG:
        if not g.snps:
            continue
        a = sum(n for (iso, h), n in hap_counts.items() if h == "hapA" and iso.startswith(g.name + "."))
        b = sum(n for (iso, h), n in hap_counts.items() if h == "hapB" and iso.startswith(g.name + "."))
        ase_rows.append({"gene": g.name, "hapA_reads": a, "hapB_reads": b, "hapA_share": round(a / (a + b), 3),
                         "expected_imbalanced": int(abs(a / (a + b) - 0.5) >= 0.25)})
    write_tsv(T("ase.tsv"), ase_rows, ["gene", "hapA_reads", "hapB_reads", "hapA_share", "expected_imbalanced"])
    pa_rows = []
    for g in CATALOG:
        full = next((i for i in isoforms if i.gene is g and i.op[0] == "full"), None) or next(i for i in isoforms if i.gene is g)
        for off, frac in g.pa_sites:
            site = full.exons[-1][1] + off if g.strand == "+" else full.exons[0][0] - off
            pa_rows.append({"gene": g.name, "chrom": g.chrom, "strand": g.strand, "site": site, "offset": off, "frac": frac,
                            "n_reads_designed": round(frac * sum(max(1, round(i.n_flnc * scale)) for i in isoforms if i.gene is g and i.n_flnc))})
    write_tsv(T("pa_sites.tsv"), pa_rows, ["gene", "chrom", "strand", "site", "offset", "frac", "n_reads_designed"])
    pass_palen = {}
    for t in truth_reads:
        if t["polya_qc"] == "PASS":
            pass_palen.setdefault(t["isoform"], []).append(float(t["polya_len"]))
    write_tsv(T("palen.tsv"), [{"event_id": ev["event_id"], "gene": ev["gene"], "inc_isoforms": ev["inc_isoforms"], "exc_isoforms": ev["exc_isoforms"],
                                "inc_mean": g.palen.get(ev["inc_isoforms"].split(".")[-1], DEFAULT_PALEN)[0],
                                "exc_mean": g.palen.get(ev["exc_isoforms"].split(".")[-1], DEFAULT_PALEN)[0],
                                "inc_pass_reads": len(pass_palen.get(ev["inc_isoforms"], [])), "exc_pass_reads": len(pass_palen.get(ev["exc_isoforms"], [])),
                                "expected_significant": int(g.palen.get(ev["inc_isoforms"].split(".")[-1], DEFAULT_PALEN)[0]
                                                            != g.palen.get(ev["exc_isoforms"].split(".")[-1], DEFAULT_PALEN)[0]),
                                "detectable_published": palen_detectable(pass_palen.get(ev["inc_isoforms"], []), pass_palen.get(ev["exc_isoforms"], []))}
                               for ev in events for g in CATALOG if g.name == ev["gene"] and ev["type"] != "MXE"],
              ["event_id", "gene", "inc_isoforms", "exc_isoforms", "inc_mean", "exc_mean", "inc_pass_reads", "exc_pass_reads", "expected_significant",
               "detectable_published"])
    with open(T("snps.vcf"), "w") as f:
        f.write("##fileformat=VCFv4.2\n##INFO=<ID=GENE,Number=1,Type=String,Description=\"gene\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, vid, ref, alt, gene in sorted(vcf):
            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t60\tPASS\tGENE={gene}\n")
    ont_rng = random.Random(seed + 1)              # own stream: everything above stays byte-identical
    for chem in ONT_CHEMISTRY:
        write_ont(ont_rng, simulate_ont(ont_rng, seq_upper, isoforms, hapb, scale, chem), os.path.join(out, "ont", chem), chem)
    write_configs(out, scale)
    with open(T("README.md"), "w") as f:
        f.write(TRUTH_README)
    return out


TRUTH_README = """# synthetic_v2 truth tables (0-based, half-open coordinates unless stated)

| file | one row per | key columns |
|---|---|---|
| genes.tsv | gene | span, strand, category, softmasked, motif |
| isoforms.tsv | isoform (full set incl. novel) | exons, junc_chain, n_flnc, annotated (=in genes.gtf), junction_filter (legacy mapping filterByJunc outcome: pass / fail:chain-support / fail:consensus), op, wobble_shift, palen_mean |
| as_events.tsv | designed AS event | type, coords (find_as output convention), inc/exc isoforms, lr_psi, expected_detection in raw LR output (yes / no:junc-filter / no:R-2 / no:offset / no:MXE-unsupported / no:min-reads / no:aligner-microexon / both:A5SS+A3SS), expected_confident after the NGS confident filter (adds no:ngs-unsupported for alignment-wobble isoforms) |
| isoform_support.tsv | isoform | FLNC count and NGS fragment count per library |
| junction_support.tsv | designed junction | carrier isoforms, annotated, FLNC reads, designed NGS reads per library and in total (>= 5 bp anchors, as regtools -a 5 counts; the committed config merges all four libraries into one hisat2 run), ngs_detectable (0 when a flanking exon is shorter than the aligner can place, R-11), consensus for the legacy junction filter (annotated, or >= 3 FLNC reads, or > 10 NGS reads in total) |
| read_truth.tsv | FLNC read | isoform, haplotype, PA site, poly(A) length + nanopolish qc, jitter/truncation, wobble, errors, primers + zmw + `flnc_recovered` (what lima/isoseq refine keep from ccs.bam) |
| diff_psi.tsv | NGS differential event | designed PSI per condition, expected_significant |
| ase.tsv | gene with SNPs | designed FL reads per haplotype over all isoforms, hapA share, expected_imbalanced (|share - 0.5| >= 0.25: allele-specific expression) |
| allele_bias.tsv | ASAS event | P(hapA | inclusion), designed FL reads per haplotype of the inclusion/exclusion isoform (exact quota), expected_significant under the published criterion (>= 5 / <= 3 reads per haplotype per isoform, opposite majority haplotypes, chi-square p <= 0.001) |
| pa_sites.tsv | PA site | genomic site, fraction, designed reads |
| palen.tsv | AS event | designed poly(A) means of inc/exc isoforms, PASS reads per side, expected_significant (means differ by design), detectable_published (the published t-test + Kruskal <= 0.001 rule with >= 5 PASS reads per side can resolve it at the designed depth) |
| genome_features.tsv | interval | softmasked / noncanonical_splice / te_repeat |
| snps.vcf | heterozygous SNP | 1-based VCF, INFO GENE |
| genes.full.gtf | all isoforms | truth annotation (genes.gtf is the *reference subset* the pipeline sees) |
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", nargs="?", default="test_data/synthetic_v2")
    ap.add_argument("--seed", type=int, default=20261014)
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--error-rate", type=float, default=0.003)
    ap.add_argument("--stranded", action="store_true")
    a = ap.parse_args(argv)
    generate(a.out, a.seed, a.scale, a.error_rate, a.stranded)
    print(f"synthetic v2 written to {a.out}")


if __name__ == "__main__":
    main()
