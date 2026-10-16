"""Characterization of asas, recorded from the code before the G5 refactor (2026-10-15).
relationshipBetweenAlleleSpecifiAndAS on partial and nopartial tables: significant on different haplotypes, counts
outside the tested range, one non-zero haplotype, the same peak haplotype, not significant, isoforms missing -- and the
published quirk that an event skipped on the partial table is not tested on the nopartial one. getHaploIsoformSeq's
representative isoform per haplotype (junction chain with most reads, longest isoform, isoforms missing from the
isoform file, a locus without a table). ase_sample's freebayes SNP filter and SNP table, the no-SNP exit and the input
errors (the gate never enables freebayes)."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from iflas.stages.asas import ase_sample, getHaploIsoformSeq, relationshipBetweenAlleleSpecifiAndAS

EXPECTED = json.loads(r'''{
 "ase_calls": [
  "['bash', '-c', 'freebayes-parallel split.region.txt 2 -f <tmp>/genome.fa <tmp>/s1/mapping/rna-seq/reassembly/tmp.bam -C 10 -F 0.2 --min-base-quality 20 > freebayes.SNP.vcf']",
  "['bamtools', 'split', '-in', '<tmp>/s1/mapping/rna-seq/reassembly/tmp.bam', '-reference']",
  "['bash', '-c', \"cut -f 1-12 representIso.bed | bedtools getfasta -fi <tmp>/genome.fa -bed - -name -split | seqkit replace -w 0 -p '(.*?):(.*)' -r '$1' | seqkit rmdup -n > representIso.fa\"]",
  "seqkit grep -p 'REF$' -r representIso.fa > representIso.REF.fa",
  "seqkit grep -p 'ALT$' -r representIso.fa > representIso.ALT.fa"
 ],
 "ase_files": {
  "alt_reads.lst": [],
  "freebayes.SNP.vcf": [
   "##fileformat=VCFv4.2",
   "##contig=<ID=chr1,length=1000>",
   "##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">",
   "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type\">",
   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
   "chr1\t150\t.\tA\tC\t.\tPASS\tAF=0.5;TYPE=snp",
   "chr1\t160\t.\tG\tT\t.\tPASS\tAF=1;TYPE=snp",
   "chr1\t190\t.\tG\tA\t.\tPASS\tAF=0.5;TYPE=snp",
   "chr1\t200\t.\tC\tT\t.\tPASS\tAF=0.5,0.5;TYPE=snp,del"
  ],
  "ref_reads.lst": [],
  "representIso.bed": [
   "chr1\t100\t900\tPB.1_size5_ALT\t0\t+\t100\t900\t0,0,0\t2\t100,600\t0,200\tgene1",
   "chr1\t100\t900\tPB.3_size4_REF\t0\t+\t100\t900\t0,0,0\t2\t100,600\t0,200\tgene1"
  ],
  "snp_position.txt": [
   "chr1\t150\tA\tREF",
   "chr1\t150\tC\tALT",
   "chr1\t190\tG\tREF",
   "chr1\t190\tA\tALT"
  ],
  "split.region.txt": [
   "chr1:0-1000"
  ]
 },
 "ase_no_snp": {
  "alt_fa": "None",
  "ase_dir": "<tmp>/s2/ase",
  "merged": "[]",
  "ref_fa": "None",
  "snp_positions": "<tmp>/s2/ase/snp_position.txt"
 },
 "ase_no_snp_calls": [],
 "ase_out": {
  "alt_fa": "<tmp>/s1/ase/representIso.ALT.fa",
  "ase_dir": "<tmp>/s1/ase",
  "merged": "[]",
  "ref_fa": "<tmp>/s1/ase/representIso.REF.fa",
  "snp_positions": "<tmp>/s1/ase/snp_position.txt"
 },
 "haplo_bed": [
  "chr1\t100\t900\tPB.1_size5_REF\t0\t+\t100\t900\t0,0,0\t3\t100,100,400\t0,200,400\tgene1"
 ],
 "haplo_calls": [
  "['bash', '-c', \"cut -f 1-12 representIso.bed | bedtools getfasta -fi <tmp>/genome.fa -bed - -name -split | seqkit replace -w 0 -p '(.*?):(.*)' -r '$1' | seqkit rmdup -n > representIso.fa\"]",
  "seqkit grep -p 'REF$' -r representIso.fa > representIso.REF.fa",
  "seqkit grep -p 'ALT$' -r representIso.fa > representIso.ALT.fa"
 ],
 "haplo_counts": {
  "PB.1_size5_REF": [
   9.0,
   "gene1"
  ]
 },
 "haplo_fa": [
  "<tmp>/representIso.REF.fa",
  "<tmp>/representIso.ALT.fa"
 ],
 "rel_False_True": {
  "nopartial": {
   "E1": {
    "hapA": [
     "PB.1.1"
    ],
    "hapB": [
     "PB.1.2"
    ]
   },
   "E2": {
    "hapA": [
     "PB.2.1"
    ],
    "hapB": [
     "PB.2.2"
    ]
   },
   "E4": {
    "hapA": [
     "PB.4.1"
    ],
    "hapB": [
     "PB.4.2"
    ]
   }
  },
  "partial": {}
 },
 "rel_True_False": {
  "nopartial": {},
  "partial": {
   "E1": {
    "hapA": "PB.1.1",
    "hapB": "PB.1.2"
   }
  }
 },
 "rel_True_True": {
  "nopartial": {
   "E1": {
    "hapA": [
     "PB.1.1"
    ],
    "hapB": [
     "PB.1.2"
    ]
   },
   "E2": {
    "hapA": [
     "PB.2.1"
    ],
    "hapB": [
     "PB.2.2"
    ]
   }
  },
  "partial": {
   "E1": {
    "hapA": "PB.1.1",
    "hapB": "PB.1.2"
   }
  }
 },
 "rel_none_no_events": {
  "nopartial": {},
  "partial": {}
 }
}''')

HAPS = ["hapA", "hapB", "hapC"]
PARTIAL = {"PB.1.1": [30, 0, 0], "PB.1.2": [0, 30, 1], "PB.2.1": [30, 0, 0], "PB.2.2": [1, 0, 0],
           "PB.3.1": [6, 4, 5], "PB.3.2": [6, 5, 4], "PB.4.1": [30, 0, 0], "PB.4.2": [6, 0, 0],
           "PB.5.1": [100, 0, 0], "PB.5.2": [60, 30, 3], "PB.6.1": [8, 3, 6], "PB.6.2": [7, 2, 9]}
NOPARTIAL = dict(PARTIAL, **{"PB.2.1": [30, 0, 1], "PB.2.2": [0, 30, 0], "PB.4.1": [30, 0, 1], "PB.4.2": [0, 30, 0]})
EVENTS = dict({f"E{i}": [[f"PB.{i}.1"], [f"PB.{i}.2"]] for i in range(1, 7)}, E7=[["PB.9.9"], ["PB.1.2"]])


def table(path, columns, haps=HAPS):
    rows = ["# IsoPhase", "haplotype\t" + "\t".join(columns)]
    rows += [h + "\t" + "\t".join(str(columns[c][i]) for c in columns) for i, h in enumerate(haps)]
    Path(path).write_text("\n".join(rows) + "\n")


def vcf(path, records, info=""):
    head = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000>",
            '##INFO=<ID=AF,Number=A,Type=Float,Description="af">', '##INFO=<ID=TYPE,Number=A,Type=String,Description="type">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    Path(path).write_text("\n".join(head + [f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{inf or info or '.'}" for pos, ref, alt, inf in records]) + "\n")


def bed12plus(chrom, name, exons, extra):
    s, e = exons[0][0], exons[-1][1]
    return "\t".join(map(str, [chrom, s, e, name, 0, "+", s, e, "0,0,0", len(exons), ",".join(str(b - a) for a, b in exons),
                               ",".join(str(a - s) for a, _ in exons)] + list(extra))) + "\n"


class Runner:
    def __init__(self, tmp):
        self.tmp, self.calls = str(tmp), []

    def run(self, cmd, **kwargs):
        self.calls.append(str(cmd).replace(self.tmp, "<tmp>"))
        if "freebayes-parallel" in str(cmd):
            vcf("freebayes.SNP.vcf", [(150, "A", "C", "AF=0.5;TYPE=snp"), (160, "G", "T", "AF=1;TYPE=snp"),
                                      (190, "G", "A", "AF=0.5;TYPE=snp"), (200, "C", "T", "AF=0.5,0.5;TYPE=snp,del")])
        return SimpleNamespace(returncode=0, stdout="", stderr="")


def rel_snapshot(tmp_path):
    table(tmp_path / "phased.partial.cleaned.human_readable.txt", PARTIAL)
    table(tmp_path / "phased.nopartial.cleaned.human_readable.txt", NOPARTIAL)
    snap = {f"rel_{p}_{n}": relationshipBetweenAlleleSpecifiAndAS(partial=p, nopartial=n, isoPairs=EVENTS)
            for p, n in ((True, True), (False, True), (True, False))}
    snap["rel_none_no_events"] = relationshipBetweenAlleleSpecifiAndAS(isoPairs={})
    with pytest.raises(Exception, match="partial or nopartial"):
        relationshipBetweenAlleleSpecifiAndAS(isoPairs=EVENTS)
    return snap


def haplo_snapshot(tmp_path):
    loci = [tmp_path / "loci" / "PB.1_size5", tmp_path / "loci" / "PB.2_size3"]
    for d in loci:
        d.mkdir(parents=True)
    table(loci[0] / "phased.partial.cleaned.human_readable.txt",
          {"PB.1.1": [5, 0], "PB.1.2": [4, 0], "PB.1.3": [6, 0], "PB.1.4": [10, 7], "PB.1.5": [0, 0]}, haps=["A", "C"])
    vcf(loci[0] / "phased.partial.cleaned.vcf", [(150, "A", "C", "")])
    (tmp_path / "iso.bed12+").write_text(
        bed12plus("chr1", "PB.1.1", [(100, 200), (300, 400), (500, 900)], ["gene1"])
        + bed12plus("chr1", "PB.1.2", [(150, 200), (300, 400), (500, 800)], ["gene1"])
        + bed12plus("chr1", "PB.1.3", [(100, 200), (300, 950)], ["gene1"])
        + bed12plus("chr1", "PB.1.5", [(100, 200), (300, 400), (500, 990)], ["gene1"]))
    runner = Runner(tmp_path)
    refFa, altFa, counts = getHaploIsoformSeq(loci, str(tmp_path / "iso.bed12+"), str(tmp_path / "genome.fa"), runner=runner)
    return {"haplo_bed": Path("representIso.bed").read_text().splitlines(), "haplo_counts": {k: [float(v[0]), v[1]] for k, v in counts.items()},
            "haplo_fa": [str(refFa).replace(str(tmp_path), "<tmp>"), str(altFa).replace(str(tmp_path), "<tmp>")], "haplo_calls": runner.calls}


def ase_snapshot(tmp_path):
    base = tmp_path / "s1"
    loci = {name: base / "asas" / "by_loci" / name for name in ("PB.1_size5", "PB.2_size3", "PB.3_size4")}
    for name, d in loci.items():
        d.mkdir(parents=True)
        if name != "PB.2_size3":
            table(d / "phased.partial.cleaned.human_readable.txt", {"PB.1.1": [5, 0]}, haps=["A", "C"])
    vcf(loci["PB.1_size5"] / "phased.partial.cleaned.vcf", [(150, "A", "C", ""), (160, "G", "T", ""), (170, "AT", "A", ""), (180, "A", "C,G", "")])
    vcf(loci["PB.2_size3"] / "phased.partial.cleaned.vcf", [(120, "A", "C", "")])
    vcf(loci["PB.3_size4"] / "phased.partial.cleaned.vcf", [(150, "A", "C", ""), (190, "G", "A", "")])
    (base / "refine").mkdir()
    (base / "refine" / "isoformGrouped.bed12+").write_text(bed12plus("chr1", "PB.1.1", [(100, 200), (300, 900)], ["gene1"]))
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "A" * 1000 + "\n")
    Path(str(genome) + ".fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    runner = Runner(tmp_path)
    out = ase_sample(str(base), str(genome), ngs_left="r1.fq", ngs_right="r2.fq", use_freebayes=True, threads=2, runner=runner)
    ase_dir = base / "ase"
    snap = {"ase_out": {k: str(v).replace(str(tmp_path), "<tmp>") for k, v in out.items()}, "ase_calls": runner.calls,
            "ase_files": {p.name: p.read_text().splitlines() for p in sorted(ase_dir.iterdir()) if p.is_file()}}
    # a phased VCF without records: nothing to quantify
    base2 = tmp_path / "s2"
    d = base2 / "asas" / "by_loci" / "PB.1_size5"
    d.mkdir(parents=True)
    table(d / "phased.partial.cleaned.human_readable.txt", {"PB.1.1": [5, 0]}, haps=["A", "C"])
    vcf(d / "phased.partial.cleaned.vcf", [])
    runner2 = Runner(tmp_path)
    out2 = ase_sample(str(base2), str(genome), ngs_left="r1.fq", ngs_right="r2.fq", runner=runner2)
    snap["ase_no_snp"] = {k: str(v).replace(str(tmp_path), "<tmp>") for k, v in out2.items()}
    snap["ase_no_snp_calls"] = runner2.calls
    for kw, msg in (({}, "needs the NGS reads"), ({"ngs_left": "a.fq;b.fq", "ngs_right": "a2.fq"}, "not match")):
        with pytest.raises(ValueError, match=msg):
            ase_sample(str(base2), str(genome), runner=runner2, **kw)
    return snap


def snapshot(tmp_path):
    return json.loads(json.dumps({**rel_snapshot(tmp_path), **haplo_snapshot(tmp_path), **ase_snapshot(tmp_path)}, sort_keys=True))


def test_asas_matches_the_pre_refactor_record(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert snapshot(tmp_path) == EXPECTED
