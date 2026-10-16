#!/usr/bin/env python
"""Reads simulated from a real reference (real-reference check, 2026-10-16): the synthetic v2 read models -- PacBio
CCS reads with primers, ONT direct RNA (RNA002 FASTQ + simulated fast5, RNA004 dorado BAM), paired NGS in two
conditions x two replicates, heterozygous SNPs, PA sites, poly(A) tails -- applied to the annotated transcripts of a
seeded choice of real genes, so the whole chain runs on real chromosome and gene names, repeats, intron lengths and
overlapping genes.

  python scripts/make_realref_data.py <init-ref base> <out> [--genome G.fa] [--gtf G.gtf] [--go GOTerms.csv.gz]

<init-ref base> is the prefix `iflas init-ref` wrote (``<base>.gpe``, ``.len.txt``, ``.bed``, ``.mm2``, ``_hisat2``);
the genome and GTF default to ``<base>.fa`` and ``<base>.gtf`` (init-ref's copies of compressed or GFF3 inputs).
Truth tables go to <out>/truth; expected_as.tsv is what iFLAS's own find_as reports on the simulated multi-exon
isoforms with >= 2 long reads, i.e. the events a perfect reconstruction would give.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import random
import re
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_fake_fast5  # noqa: E402
import make_synthetic_data_v2 as v2  # noqa: E402

MAX_INTRON = 50_000          # minimap2.max_intron_length (minimap2 -G) of the iFLAS config
SCAFFOLD_LEN = 10_000_000    # sequences shorter than this count as scaffolds
QUOTA = {"long_gene": 100, "over_max_intron": 10, "many_isoforms": 100, "overlap": 300, "mono_exon": 150,
         "scaffold": 100, "multi_isoform": 1500}
DESIGNS = {"asas": 60, "ase": 40, "diff": 80, "diff_null": 40, "palen": 50, "apa": 80}
ONT_SCALE = 0.4              # ONT reads per chemistry relative to PacBio
NGS_SCALE = 2.8              # v2 draws 3 x n_flnc NGS pairs per isoform and library
MODEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "test_data", "synthetic_v2",
                     "r9.4_70bps.u_to_t_rna.5mer.template.model")
LIBS = [f"{c}_{r}" for c in v2.CONDITIONS for r in v2.REPLICATES]


class Genome:
    """The v2 simulators' ``seq_upper`` backed by an indexed FASTA: ``genome[chrom][s:e]`` and ``genome[chrom][p]``."""

    def __init__(self, path):
        import pysam
        self.fa, self.cache = pysam.FastaFile(path), {}

    def __getitem__(self, chrom):
        return _Chrom(self, chrom)

    def fetch(self, chrom, s, e):
        """Like slicing v2's in-memory genome: empty when a 5' truncation or TES jitter runs past a short terminal exon."""
        s = max(0, s)
        if e <= s:
            return ""
        if (chrom, s, e) not in self.cache:
            self.cache[(chrom, s, e)] = self.fa.fetch(chrom, s, e).upper()
        return self.cache[(chrom, s, e)]


class _Chrom:
    def __init__(self, genome, chrom):
        self.genome, self.chrom = genome, chrom

    def __getitem__(self, k):
        return self.genome.fetch(self.chrom, k.start, k.stop) if isinstance(k, slice) else self.genome.fetch(self.chrom, k, k + 1)


def read_transcripts(gtf):
    tx = {}
    with open(gtf) as f:
        for line in f:
            col = line.rstrip("\n").split("\t")
            if len(col) > 8 and col[2] == "exon":
                tid = re.search(r'transcript_id "([^"]+)"', col[8]).group(1)
                t = tx.setdefault(tid, {"gene": re.search(r'gene_id "([^"]+)"', col[8]).group(1), "chrom": col[0],
                                        "strand": col[6], "exons": []})
                t["exons"].append((int(col[3]) - 1, int(col[4])))
    for t in tx.values():
        t["exons"].sort()
    return tx


def gene_table(tx, chrom_len):
    """gene -> chrom, strand, span, transcripts and the structural flags the selection draws on."""
    genes = {}
    for tid in sorted(tx):
        t = tx[tid]
        genes.setdefault(t["gene"], {"chrom": t["chrom"], "strand": t["strand"], "tx": [], "partners": []})["tx"].append(tid)
    for g in genes.values():
        exons = [tx[t]["exons"] for t in g["tx"]]
        introns = [[b[0] - a[1] for a, b in zip(e, e[1:])] for e in exons]
        g["start"], g["end"] = min(e[0][0] for e in exons), max(e[-1][1] for e in exons)
        flags = {"many_isoforms" if len(exons) >= 10 else "multi_isoform" if len(exons) > 1 else "single_isoform"}
        if any(i and max(i) > MAX_INTRON for i in introns):
            flags.add("over_max_intron")
        elif any(sum(i) > MAX_INTRON for i in introns):
            flags.add("long_gene")                        # BUG-193: every intron fits, their sum does not
        if all(len(e) == 1 for e in exons):
            flags.add("mono_exon")
        if chrom_len[g["chrom"]] < SCAFFOLD_LEN:
            flags.add("scaffold")
        g["flags"] = flags
    by_chrom = {}
    for name, g in genes.items():
        by_chrom.setdefault(g["chrom"], []).append((g["start"], g["end"], name))
    for spans in by_chrom.values():
        reach, holder = -1, None
        for s, e, name in sorted(spans):
            if s < reach:                                 # nested, overlapping or antisense neighbours
                for a, b in ((name, holder), (holder, name)):
                    genes[a]["flags"].add("overlap")
                    genes[a]["partners"].append(b)
            if e > reach:
                reach, holder = e, name
    return genes


def choose_genes(rng, genes, n):
    chosen = {}
    for flag, quota in QUOTA.items():
        pool = sorted(name for name, g in genes.items() if flag in g["flags"] and name not in chosen)
        rng.shuffle(pool)
        for name in pool[:quota]:
            chosen.update(dict.fromkeys([name] + (genes[name]["partners"] if flag == "overlap" else [])))
    rest = sorted(name for name, g in genes.items() if name not in chosen and "single_isoform" in g["flags"])
    rng.shuffle(rest)
    chosen.update(dict.fromkeys(rest[:max(0, n - len(chosen))]))
    return list(chosen)


def build_isoforms(rng, tx, genes, chosen):
    """v2 Gene/Isoform objects on the annotated transcripts (the one with most exons is "full"); log-normal gene
    expression split at random between the isoforms; one gene in ten with >= 4 exons gets an unannotated exon skip."""
    isoforms = []
    for name in chosen:
        g = genes[name]
        tids = sorted(g["tx"], key=lambda t: (-len(tx[t]["exons"]), tx[t]["exons"][0][0] - tx[t]["exons"][-1][1], t))
        full = tx[tids[0]]["exons"]
        gene = v2.Gene(name, g["chrom"], g["strand"], [e - s for s, e in full], [b[0] - a[1] for a, b in zip(full, full[1:])],
                       n_full=0, category="")
        total = min(400, max(4, round(rng.lognormvariate(2.5, 0.8)))) * (5 if "many_isoforms" in g["flags"] else 1)
        weights = [rng.gammavariate(1, 1) for _ in tids]
        isos = [v2.Isoform(t, gene, list(tx[t]["exons"]), round(total * w / sum(weights)), True, ("full",) if i == 0 else ("alt",))
                for i, (t, w) in enumerate(zip(tids, weights))]
        if len(full) >= 4 and rng.random() < 0.1:
            k = rng.randint(1, len(full) - 2)
            novel = v2.Isoform(f"{name}_novel1", gene, full[:k] + full[k + 1:], max(3, total // 3), False, ("skip", k, k))
            if novel.junc_chain not in {i.junc_chain for i in isos}:
                isos.append(novel)
        isoforms.extend(isos)
    return isoforms


def assign_designs(rng, isoforms, eventful):
    """Allele-specific splicing / expression, differential PSI, poly(A) shifts and a second PA site on genes with exactly
    two annotated isoforms between which find_as reports an event (``eventful``; one design per gene); both isoforms
    get >= 30 PacBio reads."""
    by_gene = {}
    for iso in isoforms:
        by_gene.setdefault(iso.gene.name, []).append(iso)
    eligible = sorted(g for g, isos in by_gene.items() if g in eventful and len(isos) == 2 and all(i.annotated for i in isos))
    rng.shuffle(eligible)
    quota, labels = dict(DESIGNS), []
    while any(quota.values()):                           # round robin: every design gets genes when few are eligible
        for design in DESIGNS:
            if quota[design]:
                labels.append(design)
                quota[design] -= 1
    for name, design in zip(eligible, labels):
        isos = by_gene[name]
        g = isos[0].gene
        g.category = design
        for iso in isos:
            iso.n_flnc = max(iso.n_flnc, 30)
        if design == "asas":          # the published rule wants <= 3 minor-haplotype reads per isoform (R-15)
            g.snps, g.hap_bias = 5, 0.97
        elif design == "ase":
            g.snps, g.expr_bias = 5, 0.85
        elif design in ("diff", "diff_null"):
            g.ngs_psi = {"condA": 0.8, "condB": 0.2} if design == "diff" else {"condA": 0.5, "condB": 0.5}
        elif design == "palen":
            g.palen = {isos[0].name.split(".")[-1]: (40.0, 8.0), isos[1].name.split(".")[-1]: (100.0, 10.0)}
        else:
            g.pa_sites = ((0, 0.6), (250, 0.4))


def write_fastq_gz(path, records):
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=1, mtime=0) as gz:
        for k in range(0, len(records), 10000):
            gz.write("".join(f"@{n}\n{s}\n+\n{'I' * len(s)}\n" for n, s in records[k:k + 10000]).encode())


def expected_as(isoforms, long_reads, ref_gpe, work, multi_exon_only=True):
    """find_as on the isoforms with >= 2 long reads -> rows (type, chrom, start, end, name). Single-exon isoforms are
    left out by default: they never reach the HQ set (R-2)."""
    from iflas.stages.find_as import find_as_lr
    os.makedirs(work, exist_ok=True)
    bed, group = os.path.join(work, "isoforms.bed12+"), os.path.join(work, "isoforms.group.txt")
    with open(bed, "w") as b, open(group, "w") as gr:
        for iso in isoforms:
            names = long_reads.get(iso.name, [])
            if len(names) < 2 or (multi_exon_only and len(iso.exons) < 2):
                continue
            s, e = iso.exons[0][0], iso.exons[-1][1]
            b.write("\t".join(map(str, [iso.chrom, s, e, iso.name, 0, iso.strand, s, e, "0,0,0", len(iso.exons),
                                        ",".join(str(y - x) for x, y in iso.exons), ",".join(str(x - s) for x, _ in iso.exons),
                                        iso.gene.name, iso.gene.name, iso.name])) + "\n")
            gr.write(f"{iso.name}\t{','.join(names)}\n")
    out = find_as_lr(bed, group, ref_gpe, os.path.join(work, "find_as"))
    rows = []
    for typ in ("IR", "SE", "A5SS", "A3SS"):
        with open(out[typ]) as f:
            rows += [{"type": typ, "chrom": c[0], "start": c[1], "end": c[2], "strand": c[5], "event": c[3]}
                     for c in (line.rstrip("\n").split("\t") for line in f if line.strip())]
    return rows


def write_gene2go(go_csv, path):
    """MaizeGDB GOTerms table (gene_model, ..., obo_terms) -> iFLAS gene2go (gene<TAB>go)."""
    with gzip.open(go_csv, "rt", newline="") if go_csv.endswith(".gz") else open(go_csv, newline="") as f, open(path, "w") as out:
        out.write("gene\tgo\n")
        for row in csv.DictReader(f):
            for term in sorted(set(re.findall(r"GO:\d{7}", row.get("obo_terms") or ""))):
                out.write(f"{row['gene_model']}\t{term}\n")


def write_configs(out, ref, threads, gene2go):
    reference = "\n".join(f"    {k}: {v}" for k, v in ref.items())
    go_line = "\n  gene2go: gene2go.tsv" if gene2go else ""
    for name, merge, out_dir in (("config.yaml", "true", "iflas_out"), ("config_conditions.yaml", "false", "iflas_out_conditions")):
        with open(os.path.join(out, name), "w") as f:
            f.write(f"""# real-reference check (scripts/make_realref_data.py); values as in examples/config.yaml
reference:
  - ref_strain: B73
{reference}
ccs:
  min_pass: 2
  min_predicted_accuracy: 0.9
minimap2:
  max_intron_length: {MAX_INTRON}
collapse:
  min_identity: 0.9
  min_coverage: 0.9
pu_filter:
  filter_score: 0.5
  random_state: 0
  pos_fl_coverage: 2
  pos_min_junc_rpkm: 0.05
options:
  threads: {threads}
  memory: 8000M{go_line}
  merge_data_from_same_strain: {merge}
  ont_trim_polya: true
dirs:
  out_dir: {out_dir}
""")
    cols = ["sample_name", "project_name", "ref_strain", "strain", "condition", "replicate", "tgs_plat", "strategy",
            "data_location", "primer", "ngs_left_reads", "ngs_right_reads", "ngs_reads_paired", "ngs_reads_length",
            "use_fmlrc2", "fast5_dir"]
    tgs = {"pacbio": ("pacbio", "sequel", "pacbio/hifi_reads.bam", "primer.fa", "false", "None"),
           "ont002": ("nanopore", "None", "ont/rna002/reads.fastq", "None", "true", "ont/rna002/signal/fast5"),
           "ont004": ("nanopore", "None", "ont/rna004/dorado.bam", "None", "false", "None"),
           "none": ("None", "None", "None", "None", "false", "None")}

    def row(name, lib, kind):
        plat, strategy, data, primer, fmlrc, fast5 = tgs[kind]
        cond, rep = lib.split("_")
        return dict(zip(cols, [name, "realref", "B73", "B73", cond, rep, plat, strategy, data, primer, f"ngs/{lib}_R1.fastq.gz",
                               f"ngs/{lib}_R2.fastq.gz", "paired", "150", fmlrc, fast5]))
    # one strain sequenced on three platforms, each long-read run with its own NGS library (merged run) ...
    merged = [row("B73_pacbio", LIBS[0], "pacbio"), row("B73_ont002", LIBS[1], "ont002"), row("B73_ont004", LIBS[2], "ont004"),
              row("B73_ngs", LIBS[3], "none")]
    # ... and per library for diff_as (the PacBio run on condA_rep1)
    per_lib = [row(lib, lib, "pacbio" if i == 0 else "none") for i, lib in enumerate(LIBS)]
    v2.write_tsv(os.path.join(out, "samples.tsv"), merged, cols)
    v2.write_tsv(os.path.join(out, "samples_conditions.tsv"), per_lib, cols)


def generate(ref_base, out, *, genome=None, gtf=None, go=None, n_genes=3000, seed=20261016, threads=16):
    rng = random.Random(seed)
    genome, gtf = genome or ref_base + ".fa", gtf or ref_base + ".gtf"
    for d in ("truth", "pacbio", "ngs", "ont"):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    chrom_len = {c: int(n) for c, n in (line.split("\t")[:2] for line in open(ref_base + ".len.txt"))}
    tx = read_transcripts(gtf)
    genes = gene_table(tx, chrom_len)
    chosen = choose_genes(rng, genes, n_genes)
    isoforms = build_isoforms(rng, tx, genes, chosen)
    T = lambda name: os.path.join(out, "truth", name)  # noqa: E731
    os.makedirs(T("work"), exist_ok=True)
    structural = expected_as(isoforms, {i.name: ["a", "b"] for i in isoforms}, ref_base + ".gpe", T("work/structure"), False)
    assign_designs(rng, isoforms, {r["event"].split(":")[0] for r in structural})
    seq = Genome(genome)
    hapb, vcf = v2.plant_snps(rng, seq, isoforms)

    reads = v2.simulate_flnc(rng, seq, isoforms, hapb, 1.0, 0.003)
    truth = [{k: v for k, v in r.items() if k != "seq"} for r in reads]
    v2.write_ccs(rng, reads, os.path.join(out, "pacbio", "hifi_reads.bam"), os.devnull, truth)
    with open(os.path.join(out, "primer.fa"), "w") as f:
        f.write(f">5p\n{v2.PRIMER_5P}\n>3p\n{v2.PRIMER_3P}\n")
    v2.write_tsv(os.path.join(out, "pacbio", "read_truth.tsv"), truth,
                 ["read_id", "isoform", "gene", "haplotype", "pa_site_offset", "polya_len", "tes_jitter", "five_prime_trunc",
                  "n_errors", "has_primer_5p", "has_primer_3p", "zmw", "flnc_recovered"])
    long_reads = {}
    for t in truth:
        if t["flnc_recovered"]:
            long_reads.setdefault(t["isoform"], []).append(f"{v2.MOVIE}/{t['zmw']}/ccs")

    ont_rng = random.Random(seed + 1)
    ont = {}
    for chem in v2.ONT_CHEMISTRY:
        ont[chem] = v2.simulate_ont(ont_rng, seq, isoforms, hapb, ONT_SCALE, chem)
        v2.write_ont(ont_rng, ont[chem], os.path.join(out, "ont", chem), chem)
        for r in ont[chem]:
            long_reads.setdefault(r["isoform"], []).append(r["read_id"])
    palen_genes = {i.gene.name for i in isoforms if i.gene.category == "palen"}
    signal = [r for r in ont["rna002"] if r["gene"] in palen_genes] + ont_rng.sample(ont["rna002"], min(300, len(ont["rna002"])))
    records = list({r["read_id"]: (r["read_id"], r["seq"][:len(r["seq"]) - r["emitted_tail"]], r["polya_len"]) for r in signal}.values())
    make_fake_fast5.write_fast5(records, MODEL, os.path.join(out, "ont", "rna002", "signal"))

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "designed NGS support")      # v2's own band check is about its designed junctions
        libs, support, _ = v2.simulate_ngs(rng, seq, isoforms, hapb, NGS_SCALE, False)
    for lib, pairs in libs.items():
        write_fastq_gz(os.path.join(out, "ngs", f"{lib}_R1.fastq.gz"), [(f"{n}/1", a) for n, a, _ in pairs])
        write_fastq_gz(os.path.join(out, "ngs", f"{lib}_R2.fastq.gz"), [(f"{n}/2", b) for n, _, b in pairs])

    ref = {"ref_genome": os.path.abspath(genome), "ref_gtf": os.path.abspath(gtf)}
    ref.update({k: os.path.abspath(ref_base + s) for k, s in (("ref_size", ".len.txt"), ("ref_gpe", ".gpe"), ("ref_bed", ".bed"),
                                                              ("ref_mm2_index", ".mm2"), ("hisat2_index", "_hisat2"))})
    write_configs(out, ref, threads, bool(go))
    if go:
        write_gene2go(go, os.path.join(out, "gene2go.tsv"))

    v2.write_tsv(T("genes.tsv"), [{"gene": n, "chrom": genes[n]["chrom"], "strand": genes[n]["strand"], "start": genes[n]["start"],
                                   "end": genes[n]["end"], "n_annotated": len(genes[n]["tx"]), "flags": ",".join(sorted(genes[n]["flags"])),
                                   "design": next(i.gene.category for i in isoforms if i.gene.name == n)} for n in chosen],
                 ["gene", "chrom", "strand", "start", "end", "n_annotated", "flags", "design"])
    v2.write_tsv(T("isoforms.tsv"), [{"isoform": i.name, "gene": i.gene.name, "chrom": i.chrom, "strand": i.strand,
                                      "exons": ";".join(f"{s}-{e}" for s, e in i.exons), "junc_chain": i.junc_chain,
                                      "annotated": int(i.annotated), "n_flnc": i.n_flnc, "long_reads": len(long_reads.get(i.name, [])),
                                      "max_intron": max((e - s for s, e in i.introns), default=0),
                                      "sum_intron": sum(e - s for s, e in i.introns),
                                      **{f"ngs_{lib}": support[i.name].get(lib, 0) for lib in LIBS}} for i in isoforms],
                 ["isoform", "gene", "chrom", "strand", "exons", "junc_chain", "annotated", "n_flnc", "long_reads", "max_intron",
                  "sum_intron"] + [f"ngs_{lib}" for lib in LIBS])
    v2.write_tsv(T("expected_as.tsv"), expected_as(isoforms, long_reads, ref["ref_gpe"], T("work/reads")),
                 ["type", "chrom", "start", "end", "strand", "event"])
    hap = {}
    for r in truth + ont["rna002"] + ont["rna004"]:
        hap[(r["isoform"], r["haplotype"])] = hap.get((r["isoform"], r["haplotype"]), 0) + 1
    designed = {}
    for iso in isoforms:
        designed.setdefault(iso.gene.name, []).append(iso)
    rows = []
    for name, isos in designed.items():
        g = isos[0].gene
        if not g.category:
            continue
        counts = [[hap.get((i.name, h), 0) for h in ("hapA", "hapB")] for i in isos]
        rows.append({"gene": name, "design": g.category, "isoforms": ";".join(i.name for i in isos),
                     "hap_reads": ";".join(f"{a},{b}" for a, b in counts),
                     "asas_expected": v2.asas_significant(counts[1], counts[0]) if g.category == "asas" else "",
                     "psi_condA": g.ngs_psi.get("condA", ""), "psi_condB": g.ngs_psi.get("condB", ""),
                     "palen_means": ";".join(str(g.palen.get(i.name.split(".")[-1], v2.DEFAULT_PALEN)[0]) for i in isos),
                     "pa_sites": ";".join(str(i.exons[-1][1] + o if g.strand == "+" else i.exons[0][0] - o) for o, _ in g.pa_sites
                                          for i in isos[:1])})
    v2.write_tsv(T("designs.tsv"), rows, ["gene", "design", "isoforms", "hap_reads", "asas_expected", "psi_condA", "psi_condB",
                                          "palen_means", "pa_sites"])
    with open(T("snps.vcf"), "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, vid, r, a, gene in sorted(vcf):
            f.write(f"{chrom}\t{pos}\t{vid}\t{r}\t{a}\t60\tPASS\tGENE={gene}\n")
    return {"genes": len(chosen), "isoforms": len(isoforms), "pacbio": len(reads), "rna002": len(ont["rna002"]),
            "rna004": len(ont["rna004"]), "fast5": len(records), "ngs_pairs": sum(len(p) for p in libs.values())}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref_base", help="prefix written by `iflas init-ref` (<dir>/<genome stem>)")
    ap.add_argument("out")
    ap.add_argument("--genome", help="genome FASTA with .fai (default: <ref_base>.fa)")
    ap.add_argument("--gtf", help="annotation GTF (default: <ref_base>.gtf)")
    ap.add_argument("--go", help="MaizeGDB GOTerms CSV(.gz) -> <out>/gene2go.tsv")
    ap.add_argument("--genes", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=20261016)
    ap.add_argument("--threads", type=int, default=16, help="options.threads of the written configs")
    a = ap.parse_args(argv)
    counts = generate(a.ref_base, a.out, genome=a.genome, gtf=a.gtf, go=a.go, n_genes=a.genes, seed=a.seed, threads=a.threads)
    print(f"{a.out}: " + ", ".join(f"{k} {v:,}" for k, v in counts.items()))


if __name__ == "__main__":
    main()
