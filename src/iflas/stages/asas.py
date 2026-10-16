"""Allele-specific AS (legacy ``asas.py``): IsoPhase phasing of the FL reads per locus, then the association
between haplotypes and AS isoform pairs (chi-square, p <= 0.001) -> ``asas/partialAsRelatedHaplotype.txt``.

IsoPhase (``select_loci_to_phase.py`` / ``run_phaser.py``) lives in the ``iflas-cupcake`` env; minimap2 and
samtools mpileup run in the core env.
"""

from __future__ import annotations

import configparser
import copy
import glob
import gzip
import os
import shutil
from concurrent.futures import ThreadPoolExecutor
from difflib import SequenceMatcher
from math import ceil
from pathlib import Path

import pandas as pd
import pysam
from scipy.stats import chi2_contingency

from ..io.records import BedFile
from ..pu.cpc2 import _find_utils_dir
from ..runtime import SINGLE_THREAD_ENV, ToolRunner, get_logger, getCurrentTime, quote as q
from .as_pairs import as_event_pairs

CUPCAKE_ENV = "iflas-cupcake"
MIN_BQ = 13          # legacy asas.py globals
PLOIDY = 2

AS_FILE_NAMES = {
    "SE": "SE.confident.bed12+",
    "IR": "IR.confident.bed6+",
    "A5SS": "A5SS.confident.bed6+",
    "A3SS": "A3SS.confident.bed6+",
}


def getASpairedIsoforms(asFile, collapsedGroupFile, isoformFile, asType="SE", filterByCount=0, mergeByJunc=True):
    """``{cupcake gene: {asEvent: [inclusionIsos, exclusionIsos]}}`` of the confident events (the last pair of an event
    wins, as in the published code)."""
    asPairs = {}
    for gene, asEvent, inc, exc in as_event_pairs(asFile, collapsedGroupFile, isoformFile, asType, filterByCount, mergeByJunc):
        asPairs.setdefault(gene, {})[asEvent] = [inc, exc]
    return asPairs


HAPLOTYPE_TABLES = {"partial": "phased.partial.cleaned.human_readable.txt",
                    "nopartial": "phased.nopartial.cleaned.human_readable.txt"}
_SKIP = object()


def _allele_specific(table, isoPair):
    """One AS event on an IsoPhase table (FL reads per haplotype and isoform): the haplotypes on which the inclusion and
    the exclusion isoforms peak when their counts depend on the haplotype (chi-square p <= 0.001); ``None`` when the
    counts are outside the tested range or not significant; ``_SKIP`` when the event cannot be allele-specific (a side
    missing from the table, reads on one haplotype only, or both sides peaking on the same haplotype)."""
    incIsos = list(set(table.columns) & set(isoPair[0]))
    excIsos = list(set(table.columns) & set(isoPair[1]))
    if len(incIsos) == 0 or len(excIsos) == 0:
        return _SKIP
    incIsoDf, excIsoDf = table.loc[:, incIsos].sum(axis=1), table.loc[:, excIsos].sum(axis=1)
    mergedDf = pd.concat([incIsoDf, excIsoDf], axis=1)
    if not ((mergedDf.max() >= 5).all() and (mergedDf.min() <= 3).all()):
        return None
    observed = mergedDf[(mergedDf.T != 0).any()]
    if len(observed) == 1:
        return _SKIP
    # BUG-180: a haplotype without reads of either side stays out of the test (its zero row made chi2_contingency
    # raise "expected frequencies has a zero element" and abort the stage)
    if not chi2_contingency(observed)[1] <= 0.001:
        return None
    haploInc, haploExc = incIsoDf.idxmax(axis=0), excIsoDf.idxmax(axis=0)
    return _SKIP if haploInc == haploExc else (haploInc, haploExc)


def relationshipBetweenAlleleSpecifiAndAS(partial=False, nopartial=False, isoPairs=None):
    """Allele-specific AS events of the phased locus in the current directory, per haplotype table (see
    ``_allele_specific``). As published, an event skipped on the partial table is not tested on the nopartial one."""
    if isoPairs and not partial and not nopartial:
        raise Exception("You must specify either partial or nopartial file")
    resultDict, tables = {"partial": {}, "nopartial": {}}, {}
    kinds = [kind for kind, wanted in (("partial", partial), ("nopartial", nopartial)) if wanted]
    for asEvent, isoPair in isoPairs.items():
        for kind in kinds:
            if not os.path.exists(HAPLOTYPE_TABLES[kind]):
                continue
            if kind not in tables:
                tables[kind] = pd.read_csv(HAPLOTYPE_TABLES[kind], skiprows=[0], index_col=0, sep="\t")
            res = _allele_specific(tables[kind], isoPair)
            if res is _SKIP:
                break
            if res:
                inc, exc = "_".join(isoPair[0]), "_".join(isoPair[1])
                resultDict[kind][asEvent] = dict(zip(res, [inc, exc] if kind == "partial" else [[inc], [exc]]))
    return resultDict


def makeAbundanceFile(groupFile, outFile):
    """cupcake-style read_stat (id, length, is_fl, stat, pbid) from ``tofu.collapsed.group.txt`` (legacy port)."""
    with open(groupFile) as f, open(outFile, "w") as out:
        print("id\tlength\tis_fl\tstat\tpbid", file=out)
        for line in f:
            if not line.strip():
                continue
            pbid, members = line.strip().split("\t")
            for cid in members.split(","):
                print("\t".join([cid, "NA", "Y", "unique", pbid]), file=out)


def select_loci_to_phase(ref_genome, flnc, collapsed_gff, read_stat, work_dir, *, min_coverage=10, runner=None):
    """IsoPhase locus selection (cupcake ``select_loci_to_phase.py``, no shebang -> run through python in the
    cupcake env) -> ``<work_dir>/by_loci/<PB.X>_size<n>/``; returns the locus dirs."""
    runner = runner or ToolRunner()
    with open(flnc) as f:
        fq = f.read(1) == "@"
    runner.run("python $(which select_loci_to_phase.py) {g} {f} {gff} {stat} -c {c}{fq}".format(
        g=q(ref_genome), f=q(flnc), gff=q(collapsed_gff), stat=q(read_stat), c=min_coverage, fq=" --fq" if fq else ""),
        conda_env=CUPCAKE_ENV, cwd=str(work_dir))
    return sorted(p for p in Path(work_dir).glob("by_loci/*size*") if p.is_dir())


def run_phaser(locus_dir, *, runner=None, min_bq=MIN_BQ, ploidy=PLOIDY) -> bool:
    """Legacy ``runPhaser`` for one IsoPhase locus dir: minimap2 -> samtools mpileup -> ``run_phaser.py``
    (``--partial_ok``, BH-FDR 0.01). A phaser failure (e.g. a locus without usable SNPs) is logged, not fatal."""
    runner = runner or ToolRunner()
    locus_dir = str(locus_dir)
    cfg = configparser.RawConfigParser()
    with open(os.path.join(locus_dir, "config")) as f:
        cfg.read_string("[dummy_section]\n" + f.read())
    strand = cfg.get("dummy_section", "ref_strand")
    runner.run("minimap2 -ax splice fake.fasta ccs.fastq > ccs.sam 2> /dev/null", cwd=locus_dir)
    runner.run("samtools sort ccs.sam > ccs.sorted.bam 2> /dev/null && samtools mpileup --min-BQ {} -f fake.fasta -s ccs.sorted.bam "
               "> ccs.mpileup 2> /dev/null".format(min_bq), cwd=locus_dir)
    proc = runner.run("run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --partial_ok --strand {} "
                      "-o phased.partial -n {} --bhFDR 0.01 -p 0.01 > partial.stdout.txt 2>&1".format(strand, ploidy),
                      conda_env=CUPCAKE_ENV, cwd=locus_dir, check=False, env=SINGLE_THREAD_ENV)
    if proc.returncode:
        get_logger().warning("run_phaser.py failed in %s (rc=%s); see partial.stdout.txt", locus_dir, proc.returncode)
    return proc.returncode == 0


def phase_sample(base_dir, ref_genome, flnc, *, out_dir=None, collapsed_gff=None, group_file=None, min_coverage=10,
                 threads=1, runner=None):
    """IsoPhase for one sample (legacy ``asas()`` l.372-406): read_stat from the collapse groups, locus selection,
    then ``run_phaser`` per locus in parallel. Returns the locus dirs."""
    base = Path(base_dir).resolve()
    out_root = Path(out_dir).resolve() if out_dir else base / "asas"
    out_root.mkdir(parents=True, exist_ok=True)
    group_file = Path(group_file) if group_file else base / "collapse" / "tofu.collapsed.group.txt"
    collapsed_gff = Path(collapsed_gff) if collapsed_gff else base / "collapse" / "tofu.collapsed.good.gff"
    read_stat = out_root / "tofu.collapsed.read_stat.txt"
    makeAbundanceFile(str(group_file), str(read_stat))
    by_loci = out_root / "by_loci"
    if by_loci.exists():                                   # select_loci_to_phase.py refuses to overwrite it
        shutil.rmtree(by_loci)
    loci = select_loci_to_phase(os.path.abspath(ref_genome), os.path.abspath(flnc), str(collapsed_gff), str(read_stat), out_root,
                                min_coverage=min_coverage, runner=runner)
    with ThreadPoolExecutor(max_workers=max(1, int(threads or 1))) as pool:
        list(pool.map(lambda d: run_phaser(d, runner=runner), loci))
    return loci


def asas_sample(
    base_dir,
    *,
    out_dir=None,
    isoform_file=None,
    group_file=None,
    as_dir=None,
    loci_dir=None,
    filter_by_count=2,
    merge_by_junc=False,
    ref_genome=None,
    flnc=None,
    collapsed_gff=None,
    min_coverage=10,
    threads=1,
    runner=None,
    ngs_confirm=True,
):
    """Legacy ``asas()``: with ``ref_genome`` and ``flnc`` the FL reads are phased first (``phase_sample``);
    otherwise existing ``by_loci`` outputs are consumed. Then every confident AS isoform pair of each phased locus
    is tested for haplotype association -> ``partialAsRelatedHaplotype.txt``."""
    base = Path(base_dir).resolve()
    out_root = Path(out_dir).resolve() if out_dir else base / "asas"
    isoform_file = Path(isoform_file) if isoform_file else base / "hqIsoforms" / "hq.collapsed.bed12+"
    group_file = Path(group_file) if group_file else base / "collapse" / "tofu.collapsed.group.txt"
    as_dir = Path(as_dir) if as_dir else base / "as_events" / "ordinary_as" / "LR"
    loci_root = Path(loci_dir) if loci_dir else out_root / "by_loci"
    out_root.mkdir(parents=True, exist_ok=True)
    print(getCurrentTime() + " Identify allele-specific AS events...")
    if ref_genome and flnc:
        phase_sample(base, ref_genome, flnc, out_dir=out_root, collapsed_gff=collapsed_gff, group_file=group_file,
                     min_coverage=min_coverage, threads=threads, runner=runner)

    as_pairs = {}
    for as_type, file_name in AS_FILE_NAMES.items():
        as_file = as_dir / file_name
        if as_file.is_file():
            as_pairs[as_type] = getASpairedIsoforms(
                str(as_file),
                str(group_file),
                str(isoform_file),
                asType=as_type,
                filterByCount=filter_by_count,
                mergeByJunc=merge_by_junc,
            )

    resultDict = {}
    locus_dirs = sorted(p for p in loci_root.glob("*size*") if p.is_dir())
    if ngs_confirm and ref_genome:                       # BUG-201: phasing artefacts at splice sites are not SNPs
        locus_dirs = ngs_confirmed_loci(locus_dirs, ref_genome, base / "mapping" / "rna-seq" / "reassembly" / "tmp.bam",
                                        out_root, runner or ToolRunner())
    prevDir = os.getcwd()
    try:
        for locus in locus_dirs:
            os.chdir(locus)
            fakeGene = locus.name.split("_")[0]
            resultDict[fakeGene] = {}
            for asType in as_pairs:
                if fakeGene not in as_pairs[asType]:
                    continue
                isoPairs = as_pairs[asType][fakeGene]
                sigRelatedDict = relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=False, isoPairs=isoPairs)
                resultDict[fakeGene].update({asType: sigRelatedDict})
    finally:
        os.chdir(prevDir)

    partial_out = out_root / "partialAsRelatedHaplotype.txt"
    with open(partial_out, "w") as partialOut:
        for gene, byType in resultDict.items():
            for asType, kinds in byType.items():
                for asEvent, haplos in [item for events in kinds.values() for item in events.items()]:
                    (haplotype1, iso1), (haplotype2, iso2) = haplos.items()
                    refGene = "_".join(asEvent.split(":")[0:-1])
                    print("\t".join([gene, asType, refGene, asEvent, haplotype1, iso1, haplotype2, iso2]), file=partialOut)

    print(getCurrentTime() + " Identify allele-specific AS events done!")
    return {"asas_dir": out_root, "partial": partial_out, "loci": locus_dirs}


# ---------------------------------------------------------------------------------------------- ASE (legacy ase())
def generate_regions(fasta_index_file, size, outFile, chunks=False, chromosomes=None):
    """Genome regions of ``size`` bp (or ``size`` chunks per chromosome) for freebayes-parallel (legacy port)."""
    if not fasta_index_file.endswith(".fai"):
        fasta_index_file = fasta_index_file + ".fai"
    with open(fasta_index_file) as fasta_index, open(outFile, "w") as out:
        for line in fasta_index:
            chrom_name, chrom_length = line.strip().split("\t")[:2]
            chrom_length = int(chrom_length)
            if chromosomes is not None and chrom_name not in chromosomes:
                continue
            region_size = ceil(chrom_length / size) if chunks else size
            region_start = 0
            while region_start < chrom_length:
                region_end = min(region_start + region_size, chrom_length)
                print("{}:{}-{}".format(chrom_name, region_start, region_end), file=out)
                region_start = region_end


def vcfHaplo(cleanVCF):
    """{concatenated REF alleles: "REF", concatenated first ALT alleles: "ALT"} of a phased locus VCF (pysam)."""
    refStr, altStr = "", ""
    with pysam.VariantFile(cleanVCF) as vcf:
        for record in vcf:
            refStr += record.ref
            altStr += str(record.alts[0])
    return {refStr: "REF", altStr: "ALT"}


def assignHaplo(haplotypes, vcfHaploDict):
    """Label each IsoPhase haplotype string REF/ALT by closest match to the VCF allele strings (legacy port)."""
    haploAssigned = {}
    for haplo in haplotypes:
        best = max(vcfHaploDict, key=lambda v: SequenceMatcher(None, haplo, v).ratio())
        haploAssigned[haplo] = vcfHaploDict[best]
    return haploAssigned


def _representative(table, haplo, isoformBed):
    """The junction chain with most FL reads on a haplotype (single-exon isoforms pooled as "monoExon"; BUG-181: the
    published code pooled single-intron ones too) -> (a copy of its longest isoform, those reads); (None, 0) when no
    isoform of the haplotype is in the isoform file."""
    junc2iso = {}
    for j in table.columns[table.loc[haplo].to_numpy().nonzero()[0]]:
        if j in isoformBed:
            junc2iso.setdefault(isoformBed[j].juncChain if isoformBed[j].introns else "monoExon", []).append(j)
    if not junc2iso:
        return None, 0
    tmpJunc = max(junc2iso, key=lambda y: table.loc[haplo, junc2iso[y]].sum())
    longest = max((isoformBed[x] for x in junc2iso[tmpJunc]), key=lambda x: abs(x.chromStart - x.chromEnd))
    return copy.copy(longest), table.loc[haplo, junc2iso[tmpJunc]].sum()


def getHaploIsoformSeq(lociDir, isoformFile, refGenome, refFa=None, altFa=None, runner=None):
    """Per phased locus and haplotype, the junction chain with most FL reads is represented by its longest isoform
    (``<locus>_REF`` / ``<locus>_ALT``) -> ``representIso.bed`` -> genomic sequences -> ``representIso.{REF,ALT}.fa``
    (unless ``refFa``/``altFa`` are given). Returns (refFa, altFa, {name: [FL reads, gene]})."""
    runner = runner or ToolRunner()
    isoformBed = BedFile(isoformFile, type="bed12+").reads
    haplo2flncCount = {}
    with open("representIso.bed", "w") as representIsoOut:
        for i in lociDir:
            lociName = os.path.basename(str(i))
            partialHaplotype = os.path.join(str(i), "phased.partial.cleaned.human_readable.txt")
            cleanVCF = os.path.join(str(i), "phased.partial.cleaned.vcf")
            if not os.path.exists(partialHaplotype):
                continue
            partialDF = pd.read_csv(partialHaplotype, skiprows=[0], index_col=0, sep="\t")
            haploAssigned = assignHaplo(list(partialDF.index), vcfHaplo(cleanVCF))
            for haplo in partialDF.index:
                repIso, tmpCount = _representative(partialDF, haplo, isoformBed)
                if repIso is None:
                    continue
                repIso.name = "{}_{}".format(lociName, haploAssigned[haplo])
                haplo2flncCount[repIso.name] = [tmpCount, repIso.otherList[0]]
                print(str(repIso), file=representIsoOut)
    runner.run(["bash", "-c", "cut -f 1-12 representIso.bed | bedtools getfasta -fi {} -bed - -name -split | "
                              "seqkit replace -w 0 -p '(.*?):(.*)' -r '$1' | seqkit rmdup -n > representIso.fa".format(q(refGenome))])
    if not refFa or not os.path.isfile(refFa):
        runner.run("seqkit grep -p 'REF$' -r representIso.fa > representIso.REF.fa")
        refFa = os.path.join(os.getcwd(), "representIso.REF.fa")
    if not altFa or not os.path.isfile(altFa):
        runner.run("seqkit grep -p 'ALT$' -r representIso.fa > representIso.ALT.fa")
        altFa = os.path.join(os.getcwd(), "representIso.ALT.fa")
    return refFa, altFa, haplo2flncCount


def mergeASEfromSalmonOut(ref_quant_sf, alt_quant_sf, haplo2count, mergedASEfile):
    """Gene table of FL-read and salmon (TPM, NumReads) support per haplotype (legacy port; legacy iterated
    DataFrame columns and joined floats, which could not run — BUG-151)."""
    aseDict = {}
    for label, quant in (("REF", ref_quant_sf), ("ALT", alt_quant_sf)):
        for _, row in pd.read_csv(quant, sep="\t").iterrows():
            gene = str(row.Name).split("_")[0]
            d = aseDict.setdefault(gene, {"tgs": {"REF": 0, "ALT": 0}, "ngs": {"REF": [0, 0], "ALT": [0, 0]}, "refGene": ""})
            d["tgs"][label] = float(haplo2count[row.Name][0])
            d["ngs"][label] = [float(row.TPM), float(row.NumReads)]
            d["refGene"] = haplo2count[row.Name][1]
    with open(mergedASEfile, "w") as out:
        for gene, d in aseDict.items():
            tgsRef, tgsAlt = d["tgs"]["REF"], d["tgs"]["ALT"]
            (ngsRefTPM, ngsRefCount), (ngsAltTPM, ngsAltCount) = d["ngs"]["REF"], d["ngs"]["ALT"]
            print("\t".join(map(str, [gene, d["refGene"], tgsRef, tgsAlt, ngsRefTPM, ngsRefCount, ngsAltTPM, ngsAltCount])), file=out)


def _is_empty_fastq(path) -> bool:
    """True when a (gzipped) FASTQ exists and holds no read; a missing file is left for salmon to report."""
    if not os.path.isfile(path):
        return False
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as f:
        return not f.read(1)


def write_zero_quant(fasta, quant_dir) -> str:
    """A salmon-style quant.sf with zero TPM/NumReads for every sequence of ``fasta`` (an unexpressed haplotype)."""
    os.makedirs(quant_dir, exist_ok=True)
    rows, name, length = [], None, 0
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    rows.append((name, length))
                name, length = line[1:].split()[0], 0
            else:
                length += len(line.strip())
    if name is not None:
        rows.append((name, length))
    out = os.path.join(quant_dir, "quant.sf")
    with open(out, "w") as q:
        q.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        q.writelines(f"{n}\t{ln}\t{ln}\t0.0\t0.0\n" for n, ln in rows)
    return out


def _tagged_read_names(tagged_bam):
    """REF- and ALT-tagged read names, cut like the FASTQ ids they are looked up by (``seqkit grep --id-regexp
    '^([^/\\s]+)'``): HISAT2 keeps the ``/1`` of single-end names, so whole names never matched (BUG-173)."""
    ref_reads, alt_reads = [], []
    with pysam.AlignmentFile(tagged_bam, "rb") as bam:
        for record in bam:
            if not (record.has_tag("RG") and record.has_tag("NH")) or record.get_tag("NH") != 1:
                continue
            name = record.query_name.split("/", 1)[0]
            (ref_reads if record.get_tag("RG") == "REF" else alt_reads if record.get_tag("RG") == "ALT" else []).append(name)
    return ref_reads, alt_reads


def _freebayes_het_snps(ref_genome, hybrid_bam, threads, runner):
    """``chrom_pos`` of the heterozygous SNPs (AF 0.5, of type snp only) freebayes calls on the NGS alignments."""
    if not os.path.isfile(ref_genome + ".fai"):
        runner.run(["samtools", "faidx", ref_genome])
    generate_regions(ref_genome, 1000000, "split.region.txt")
    runner.run(["bash", "-c", "freebayes-parallel split.region.txt {} -f {} {} -C 10 -F 0.2 --min-base-quality 20 "
                              "> freebayes.SNP.vcf".format(threads, q(ref_genome), q(hybrid_bam))])
    snps = set()
    with pysam.VariantFile("freebayes.SNP.vcf") as vcf:
        for rec in vcf:
            af, typ = rec.info.get("AF"), rec.info.get("TYPE")
            af = af[0] if isinstance(af, tuple) else af
            typ = typ if isinstance(typ, tuple) else (typ,)
            if af == 0.5 and len(typ) == 1 and typ[0] == "snp":
                snps.add("{}_{}".format(rec.chrom, rec.pos))
    return snps


def _phased_records(loci):
    """(``chrom_pos``, record) of the phased VCFs of the loci with a haplotype table, first occurrence of a position only."""
    seen = set()
    for locus in loci:
        if not (locus / "phased.partial.cleaned.human_readable.txt").is_file():
            continue
        with pysam.VariantFile(str(locus / "phased.partial.cleaned.vcf")) as vcf:
            for rec in vcf:
                name = "{}_{}".format(rec.chrom, rec.pos)
                if name not in seen:
                    seen.add(name)
                    yield name, rec


def _locus_snps(locus):
    with pysam.VariantFile(str(locus / "phased.partial.cleaned.vcf")) as vcf:
        return {(rec.chrom, rec.pos) for rec in vcf}


def ngs_confirmed_loci(loci, ref_genome, ngs_bam, out_dir, runner):
    """The phased loci with at least one SNP that freebayes also calls on the NGS alignments (BUG-201). At an alternative
    splice site the reads of each isoform can align differently a few bases away; IsoPhase reads those bases as
    heterozygous SNPs, one "haplotype" per isoform, and the association test then calls the event allele-specific.
    freebayes runs over the phased SNPs' neighbourhoods with loose settings (-C 3 -F 0.05) so that the minor allele of a
    strongly allele-biased gene still counts. Without NGS alignments every locus is kept. The kept loci are listed in
    ``<out_dir>/ngs_confirmed_loci.txt`` for ase."""
    loci = [locus for locus in loci if (locus / "phased.partial.cleaned.human_readable.txt").is_file()]
    keep = loci
    if os.path.isfile(ngs_bam) and loci:
        regions = sorted({(c, max(0, p - 200), p + 200) for locus in loci for c, p in _locus_snps(locus)})
        bed, vcf = Path(out_dir) / "phased_snps.bed", Path(out_dir) / "phased_snps.ngs.vcf"
        bed.write_text("".join(f"{c}\t{a}\t{b}\n" for c, a, b in regions))
        runner.run(["bash", "-c", "bedtools merge -i {b} > {b}.merged && freebayes -f {g} -t {b}.merged -C 3 -F 0.05 "
                                  "--min-base-quality 20 {bam} > {v}".format(b=q(bed), g=q(ref_genome), bam=q(ngs_bam), v=q(vcf))])
        with pysam.VariantFile(str(vcf)) as calls:
            confirmed = {(rec.chrom, rec.pos) for rec in calls if "snp" in (rec.info.get("TYPE") or ())}
        keep = [locus for locus in loci if _locus_snps(locus) & confirmed]
    (Path(out_dir) / "ngs_confirmed_loci.txt").write_text("".join(locus.name + "\n" for locus in keep))
    return keep


def _write_snp_positions(loci, freebayes_snps):
    """``snp_position.txt`` (chrom, pos, allele, REF/ALT) of the phased single-base SNPs, confirmed by freebayes when it
    ran; returns whether the loci had any phased variant."""
    seen = False
    with open("snp_position.txt", "w") as snp_position:
        for name, rec in _phased_records(loci):
            seen = True
            if len(rec.ref) != 1 or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
                continue
            if freebayes_snps and name not in freebayes_snps:
                continue
            print("\t".join(map(str, [rec.chrom, rec.pos, rec.ref, "REF"])), file=snp_position)
            print("\t".join(map(str, [rec.chrom, rec.pos, rec.alts[0], "ALT"])), file=snp_position)
    return seen


def _tag_ngs_reads(hybrid_bam, util, threads, runner):
    """Split the NGS alignments per reference, tag the reads REF/ALT at the phased SNPs (biostar214299) and write the
    uniquely mapped REF/ALT read names to ``ref_reads.lst`` / ``alt_reads.lst``."""
    runner.run(["bamtools", "split", "-in", hybrid_bam, "-reference"])
    split_bams = sorted(glob.glob(os.path.join(os.path.dirname(hybrid_bam), "tmp.REF_*.bam")))

    def tag(bam):
        tagged = os.path.join(os.path.dirname(bam), "tagged." + os.path.basename(bam))
        runner.run(["java", "-jar", os.path.join(util, "biostar214299.jar"), "-p", "snp_position.txt", bam, "-o", tagged])
        return tagged

    with ThreadPoolExecutor(max_workers=threads) as pool:
        tagged_bams = list(pool.map(tag, split_bams))
    with open("ref_reads.lst", "w") as ref_out, open("alt_reads.lst", "w") as alt_out:
        for tagged in tagged_bams:
            ref_reads, alt_reads = _tagged_read_names(tagged)
            ref_out.write("".join(r + "\n" for r in ref_reads))
            alt_out.write("".join(r + "\n" for r in alt_reads))
    for f in split_bams + tagged_bams:
        os.remove(f)


def _quantify_repeat(i, lefts, rights, paired, ref_fa, alt_fa, threads, runner):
    """salmon quant of NGS repeat ``i``'s REF- and ALT-tagged reads against each haplotype's representative isoforms."""
    for label, lst in (("ref", "ref_reads.lst"), ("alt", "alt_reads.lst")):
        sel = []
        for side, files in (("left", lefts), ("right", rights)):
            if i < len(files):
                fq = "repeat{}.{}.{}_reads.fq.gz".format(i, label, side)
                # --id-regexp: match the BAM query names also when the FASTQ ids carry a /1 or /2 mate suffix
                reads_in = " ".join(q(x) for x in files[i].split(","))
                runner.run("seqkit grep -f {} --id-regexp '^([^/\\s]+)' {} | gzip > {}".format(lst, reads_in, fq))
                sel.append(fq)
        if _is_empty_fastq(sel[0]):          # BUG-161: salmon refuses an empty read file
            get_logger().warning("ase: no NGS read of repeat%d carries the %s alleles; their counts are 0", i, label.upper())
            write_zero_quant(ref_fa if label == "ref" else alt_fa, "repeat{}.{}_salmon_quant".format(i, label))
            continue
        reads_arg = "-1 {} -2 {}".format(*sel) if paired and len(sel) == 2 else "-r {}".format(sel[0])
        runner.run("salmon quant -l A -i {}_salmon_index {} -p {} -o repeat{}.{}_salmon_quant > repeat{}.{}_salmon_quant.log 2>&1".format(
            label, reads_arg, threads, i, label, i, label))


def ase_sample(base_dir, ref_genome, *, ngs_left=None, ngs_right=None, paired=True, out_dir=None, loci_dir=None,
               isoform_file=None, hybrid_bam=None, ref_fa=None, alt_fa=None, use_freebayes=False, threads=1, runner=None):
    """Legacy ``ase()``: allele-specific expression of the phased loci. Heterozygous SNP positions (IsoPhase VCFs,
    optionally confirmed by freebayes on the NGS BAM) tag the NGS reads REF/ALT (``bamtools split`` +
    ``biostar214299.jar``), the per-haplotype representative isoforms are quantified with salmon from those reads,
    and FL-read vs NGS support per haplotype is written to ``repeat<i>.mergedASEfile.txt``."""
    runner = runner or ToolRunner()
    base = Path(base_dir).resolve()
    out_root = Path(out_dir).resolve() if out_dir else base / "ase"
    out_root.mkdir(parents=True, exist_ok=True)
    loci = sorted(p for p in (Path(loci_dir) if loci_dir else base / "asas" / "by_loci").glob("*size*") if p.is_dir())
    confirmed = base / "asas" / "ngs_confirmed_loci.txt"
    if confirmed.is_file():                              # BUG-201: the loci asas kept
        names = set(confirmed.read_text().split())
        loci = [locus for locus in loci if locus.name in names]
    isoform_file = str(isoform_file or base / "refine" / "isoformGrouped.bed12+")
    hybrid_bam = str(hybrid_bam or base / "mapping" / "rna-seq" / "reassembly" / "tmp.bam")
    ref_genome = os.path.abspath(ref_genome)
    util = _find_utils_dir()
    threads = max(1, int(threads or 1))
    lefts = [r.strip() for r in (ngs_left or "").split(";") if r.strip()]
    rights = [r.strip() for r in (ngs_right or "").split(";") if r.strip()]
    if not lefts and not rights:
        raise ValueError("ase needs the NGS reads (ngs_left / ngs_right)")
    if paired and len(lefts) != len(rights):
        raise ValueError("The repeats of your NGS data not match between your left reads and right reads")
    skipped = {"ase_dir": str(out_root), "merged": [], "snp_positions": str(out_root / "snp_position.txt")}
    print(getCurrentTime() + " Identify allelic-specific expression genes...")
    prev = os.getcwd()
    os.chdir(out_root)
    try:
        freebayes_snps = _freebayes_het_snps(ref_genome, hybrid_bam, threads, runner) if use_freebayes else set()
        if not _write_snp_positions(loci, freebayes_snps):
            print(getCurrentTime() + " No phased heterozygous SNP under asas/by_loci: nothing to quantify, ASE skipped!")
            return dict(skipped, ref_fa=None, alt_fa=None)
        _tag_ngs_reads(hybrid_bam, util, threads, runner)
        ref_fa, alt_fa, haplo2flncCount = getHaploIsoformSeq(loci, isoform_file, ref_genome, refFa=ref_fa, altFa=alt_fa, runner=runner)
        if not all(os.path.isfile(f) and os.path.getsize(f) for f in (ref_fa, alt_fa)):
            print(getCurrentTime() + " No representative isoform for one haplotype: nothing to quantify, ASE skipped!")
            return dict(skipped, ref_fa=ref_fa, alt_fa=alt_fa)
        runner.run("salmon index -t {} -i ref_salmon_index --keepDuplicates -p {} > ase.ref_index.log 2>&1".format(q(ref_fa), threads))
        runner.run("salmon index -t {} -i alt_salmon_index --keepDuplicates -p {} > ase.alt_index.log 2>&1".format(q(alt_fa), threads))
        merged = []
        for i in range(len(lefts or rights)):
            _quantify_repeat(i, lefts, rights, paired, ref_fa, alt_fa, threads, runner)
            out_file = "repeat{}.mergedASEfile.txt".format(i)
            mergeASEfromSalmonOut("repeat{}.ref_salmon_quant/quant.sf".format(i), "repeat{}.alt_salmon_quant/quant.sf".format(i),
                                  haplo2flncCount, out_file)
            merged.append(str(out_root / out_file))
    finally:
        os.chdir(prev)
    print(getCurrentTime() + " Identify allelic-specific expression genes done!")
    return dict(skipped, merged=merged, ref_fa=ref_fa, alt_fa=alt_fa)
