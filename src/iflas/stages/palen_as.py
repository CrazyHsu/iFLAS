"""Poly(A)-length-related AS helpers migrated from legacy ``palen_as.py``."""

from __future__ import annotations

import copy
import itertools
import os
import re
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import scipy.stats as stats

from ..intervalops import getBlockLength, getOverlapOfTuple
from ..io.files import getDictFromFile, resolveDir
from ..io.records import BedFile
from ..pu.cpc2 import _find_utils_dir
from ..rscript import run_r_script
from ..runtime import ToolRunner, getCurrentTime, quote as q
from .as_pairs import as_event_pairs

AS_FILE_NAMES = {
    "SE": "SE.confident.bed12+",
    "IR": "IR.confident.bed6+",
    "A5SS": "A5SS.confident.bed6+",
    "A3SS": "A3SS.confident.bed6+",
}


def read_polya_lengths(palen_files):
    """Read one or more nanopolish ``polya_results.tsv`` files into read->length."""
    if isinstance(palen_files, (str, os.PathLike)):
        palen_files = [palen_files]
    palen_files = [str(p) for p in palen_files if str(p)]
    if not palen_files:
        raise ValueError("at least one poly(A) length file is required")

    frames = []
    for palen_file in palen_files:
        frame = pd.read_csv(palen_file, sep="\t", dtype={"contig": "string"})
        missing = {"readname", "polya_length", "qc_tag"} - set(frame.columns)
        if missing:
            raise ValueError(f"{palen_file}: missing column(s): {', '.join(sorted(missing))}")
        frames.append(frame)
    palen = pd.concat(frames)
    palen_pass = palen.loc[palen.qc_tag == "PASS", ]
    return dict(list(zip(palen_pass.readname, palen_pass.polya_length)))


def check_fast5_files(path):
    """Legacy ``checkFast5Files``: "fastq" for a fastq file, "fast5" for a directory of fast5 files, else None."""
    if os.path.isfile(path):
        return "fastq" if path.endswith(("fq", "fastq")) else None
    if not os.path.isdir(path):
        return None
    names = [x for x in os.listdir(path) if not x.endswith(".txt")]
    if names and all(x.endswith("fast5") for x in names):
        return "fast5"
    return None


def polya_len_calling(fast5_dir, raw_fq, mapped_bam, ref_genome, out_dir, *, sequencing_summary=None, threads=1, runner=None):
    """Legacy ``polyaLenCalling``: ``nanopolish index`` (fast5 dir + basecalling sequencing summary) and
    ``nanopolish polya`` on the mapped FLNC reads -> ``<out_dir>/polya_results.tsv`` (+ ``polyaTailLength.pdf``)."""
    runner = runner or ToolRunner()
    fast5_dirs = [fast5_dir] if isinstance(fast5_dir, (str, os.PathLike)) else list(fast5_dir)
    for d in fast5_dirs:
        if check_fast5_files(d) != "fast5":
            raise ValueError("Please specify the correct directory contain fast5 files! ({}; basecall externally, e.g. dorado, "
                             "and give the raw signal directory as fast5_dir)".format(d))
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = "--sequencing-summary={} ".format(q(sequencing_summary)) if sequencing_summary else ""
    runner.run("nanopolish index {} {}{}".format(" ".join("--directory={}".format(q(d)) for d in fast5_dirs), summary, q(raw_fq)),
               cwd=str(out_dir))
    runner.run("nanopolish polya --threads={} --reads={} --bam={} --genome={} > polya_results.tsv".format(
        max(1, int(threads or 1)), q(raw_fq), q(mapped_bam), q(ref_genome)), cwd=str(out_dir))
    run_r_script("grep 'PASS' polya_results.tsv | cut -f 9 | Rscript {}/hist.R -x='poly(A) tail Length' -y=Density -b=1 -d -x1=0 -x2=100 "
                 "-w=10 -p=polyaTailLength.pdf".format(q(_find_utils_dir())), cwd=str(out_dir), runner=runner)
    return str(out_dir / "polya_results.tsv")


def palen_as_sample(
    base_dir,
    palen_files=None,
    *,
    fast5_dir=None,
    ref_genome=None,
    sequencing_summary=None,
    threads=1,
    runner=None,
    out_dir=None,
    isoform_file=None,
    group_file=None,
    reads_file=None,
    as_dir=None,
    pa_cluster=None,
    filter_by_count=10,
    merge_by_junc=True,
    confident_pac_files=None,
    ref_bed=None,
):
    """Legacy ``palen_as``: poly(A) tail lengths per FL read (``polya_results.tsv`` files, or nanopolish on the
    sample's fast5 directory via ``polya_len_calling``) are tested between the isoforms of every confident AS event
    (``palenAS/mergeByJunc``) and between the PA sites of every gene (``palenAS/palenAPA``)."""
    base = Path(base_dir).resolve()
    out_root = Path(out_dir).resolve() if out_dir else base / "palenAS"
    print(getCurrentTime() + " Identify functional poly(A) tail length related to AS...")
    palen_files = [str(p) for p in (palen_files or [])]
    if not palen_files:
        if not (fast5_dir and ref_genome):
            raise ValueError("palen_as needs polya_results.tsv file(s) or a fast5 directory plus the reference genome for nanopolish")
        raw_fq = None
        for cand in (base / "preprocess" / "nanopore" / "rawFlnc.fq", base / "preprocess" / "pacbio" / "rawFlnc.fq"):
            if cand.is_file():
                raw_fq = str(cand)
        if not raw_fq:
            raise FileNotFoundError("nanopolish needs the basecalled FLNC reads at preprocess/nanopore/rawFlnc.fq")
        # nanopolish needs the alignments of the very reads it indexes: after fmlrc2 the processed FLNC are corrected
        # copies, and mapping also aligned the raw reads (legacy passed the corrected alignments: BUG-175)
        raw_bam = base / "mapping" / "rawFlnc.mm2.sorted.bam"
        mapped_bam = raw_bam if raw_bam.is_file() else base / "mapping" / "flnc.mm2.sorted.bam"
        palen_files = [polya_len_calling(fast5_dir, raw_fq, str(mapped_bam), ref_genome, out_root,
                                         sequencing_summary=sequencing_summary, threads=threads, runner=runner)]
    isoform_file = Path(isoform_file) if isoform_file else base / "hqIsoforms" / "hq.collapsed.bed12+"
    group_file = Path(group_file) if group_file else base / "collapse" / "tofu.collapsed.group.txt"
    reads_file = Path(reads_file) if reads_file else base / "mapping" / "flnc.addCVandID.bed12+"
    as_dir = Path(as_dir) if as_dir else base / "as_events" / "ordinary_as" / "LR"
    pa_cluster = Path(pa_cluster) if pa_cluster else base / "as_events" / "pa" / "paCluster.bed8+"

    flnc_reads_to_palen = read_polya_lengths(palen_files)
    collapsed_trans_to_reads = getDictFromFile(str(group_file), sep="\t", inlineSep=",", valueCol=2)
    as_pairs = {}
    for as_type, file_name in AS_FILE_NAMES.items():
        as_file = as_dir / file_name
        if as_file.is_file():
            as_pairs[as_type] = getASpairedIsoforms(
                str(as_file),
                str(group_file),
                str(isoform_file),
                asType=as_type,
                filterByCount=0,
                mergeByJunc=merge_by_junc,
            )

    out_root.mkdir(parents=True, exist_ok=True)
    prev_dir = os.getcwd()
    try:
        os.chdir(out_root)
        getPalenAS(
            flnc_reads_to_palen,
            str(isoform_file),
            str(reads_file),
            collapsedTrans2reads=collapsed_trans_to_reads,
            asPairs=as_pairs,
            filterByCount=filter_by_count,
            mergeByJunc=merge_by_junc,
        )
    finally:
        os.chdir(prev_dir)

    mode_dir = out_root / ("mergeByJunc" if merge_by_junc else "no_mergeByJunc")
    outputs = {
        "palen_dir": out_root,
        "as_dir": mode_dir,
        "as_sig": {as_type: mode_dir / f"{as_type}.palenAndAS.sig.bed12+" for as_type in as_pairs},
        "apa": None,
    }

    if pa_cluster.is_file():
        apa_dir = out_root / "palenAPA"
        apa_dir.mkdir(parents=True, exist_ok=True)
        confident = [str(p) for p in confident_pac_files or []]
        ref_params = None
        if confident:
            if not ref_bed:
                raise ValueError("--ref-bed is required when --confident-pac is provided")
            ref_params = SimpleNamespace(ref_bed=str(ref_bed))
        apa_out = apa_dir / "apaRelatedPalen.txt"
        getPalenAPA(
            str(pa_cluster),
            flnc_reads_to_palen,
            filter_by_count,
            str(apa_out),
            confidentPacFiles=confident or None,
            refParams=ref_params,
        )
        outputs["apa"] = apa_out

    return outputs


def getASpairedIsoforms(asFile, collapsedGroupFile, isoformFile, asType="SE", filterByCount=0, mergeByJunc=True):
    """``{cupcake gene: [[inclusionIsos, exclusionIsos, asEvent], ...]}`` of the confident events (every pair kept)."""
    asPairs = {}
    for gene, asEvent, inc, exc in as_event_pairs(asFile, collapsedGroupFile, isoformFile, asType, filterByCount, mergeByJunc):
        asPairs.setdefault(gene, []).append([inc, exc, asEvent])
    return asPairs


def _confident_pa_keys(confidentPacFiles, refBed):
    """``chrom_pos`` keys within 30 nt of a confident PA site (column 7 of the PAC files) or of a reference transcript's
    3' end."""
    sites = []
    for pacFile in confidentPacFiles:
        with open(pacFile) as f:
            for line in f:
                fields = line.strip("\n").split("\t")
                sites.append((fields[0], int(fields[6])))
    with open(refBed) as f:
        for line in f:
            bed = line.strip("\n").split("\t")
            sites.append((bed[0], int(bed[2] if bed[5] == "+" else bed[1])))
    return {"{}_{}".format(chrom, i) for chrom, pa in sites for i in range(pa - 30, pa + 30)}


def _pa_reads_by_gene(paClusterBed, filterByCount, confident):
    """``{gene: {chrom_pos: reads}}`` and the strand of each gene, for PA clusters with enough reads (and confident)."""
    gene2paReads, gene2strand = {}, {}
    with open(paClusterBed) as f:
        for line in f:
            infoList = line.strip("\n").split("\t")
            reads = infoList[3].split(",")
            if len(reads) < filterByCount:
                continue
            paSite = "{}_{}".format(infoList[0], infoList[7])
            if confident and paSite not in confident:
                continue
            gene2strand.setdefault(infoList[8], infoList[5])
            gene2paReads.setdefault(infoList[8], {})[paSite] = reads
    return gene2paReads, gene2strand


def _tails_differ(aPalen, bPalen):
    """t-test and Kruskal-Wallis both p <= 0.001 (a NaN p-value never passes)."""
    p1, p2 = stats.ttest_ind(aPalen, bPalen)[1], stats.kruskal(aPalen, bPalen)[1]
    return float(p1) <= 0.001 and float(p2) <= 0.001


def _differing_pa_sites(paReads, flncReads2Palen, filterByCount):
    """PA sites of a gene whose reads' tail lengths differ from those of another of its sites."""
    valid = {}
    for a, b in itertools.combinations(list(paReads), 2):
        aPalen = [float(flncReads2Palen[x]) for x in paReads[a] if x in flncReads2Palen]
        bPalen = [float(flncReads2Palen[x]) for x in paReads[b] if x in flncReads2Palen]
        if len(aPalen) < filterByCount or len(bPalen) < filterByCount:
            continue
        if len(aPalen) / float(len(bPalen)) >= 50 or len(bPalen) / float(len(aPalen)) >= 50:
            continue
        if _tails_differ(aPalen, bPalen):
            valid.setdefault(a, aPalen)
            valid.setdefault(b, bPalen)
    return valid


def getPalenAPA(paClusterBed, flncReads2Palen, filterByCount, palenAPA, confidentPacFiles=None, refParams=None):
    confident = _confident_pa_keys(confidentPacFiles, refParams.ref_bed) if confidentPacFiles else set()
    gene2paReads, gene2strand = _pa_reads_by_gene(paClusterBed, filterByCount, confident)
    with open(palenAPA, "w") as out:
        for gene, paReads in gene2paReads.items():
            if len(paReads) < 2:
                continue
            valid = _differing_pa_sites(paReads, flncReads2Palen, filterByCount)
            if len(valid) < 2:
                continue
            print("\t".join([gene, gene2strand[gene], ";".join(valid), ";".join(str(len(v)) for v in valid.values()),
                             ";".join(",".join(map(str, v)) for v in valid.values())]), file=out)


def _event_reads(isoNames, isoBedObj, collapsedTrans2reads):
    """Collapsed isoforms behind the HQ isoforms of one side of an event, and their FL reads."""
    isos = [i for x in isoNames for i in isoBedObj.reads[x].otherList[2].split(",")]
    return isos, [r for i in isos for r in collapsedTrans2reads[i]]


def _covering_intron(reads, asEvent, readBedObj):
    """IR inclusion reads that cover the whole retained intron (the last two numbers of the event name)."""
    retention = [list(map(int, re.split("[:|-]", asEvent)[-2:]))]
    return [x for x in reads if getBlockLength(getOverlapOfTuple(readBedObj.reads[x].exons, retention)) == getBlockLength(retention)]


def _write_palen_as(sigOut, gene, asType, count, sides, isoBedObj):
    """One ``<type>_<gene>_<n>.txt`` with every read's tail and, when both sides have isoforms in the HQ file, the
    longest isoform of each side as a sig.bed12+ line."""
    with open("{}_{}_{}.txt".format(asType, gene, count), "w") as out:
        for isos, tails, _ in sides:
            for read, palen in tails:
                print("\t".join(map(str, [gene, asType, "_".join(isos), read, palen])), file=out)
    objs = [[isoBedObj.reads[x] for x in isos if x in isoBedObj.reads] for isos, _, _ in sides]
    if not all(objs):
        return
    pair = "+".join("_".join(isos) for isos, _, _ in sides)
    for (isos, _, palen), isoObjs, kind in zip(sides, objs, ("inclusion", "exclusion")):
        rep = copy.copy(max(isoObjs, key=lambda x: abs(x.chromStart - x.chromEnd)))
        rep.name = "_".join(isos)
        print(str(rep) + "\t" + "\t".join([pair, kind, ",".join(map(str, palen)), str(len(palen))]), file=sigOut)


def _palen_as_type(asType, pairsByGene, flncReads2Palen, isoBedObj, readBedObj, collapsedTrans2reads, filterByCount):
    with open("{}.palenAndAS.sig.bed12+".format(asType), "w") as sigOut:
        count = 0
        for gene, pairs in pairsByGene.items():
            for inc, exc, asEvent in pairs:
                inclusionIsos, inclusionReads = _event_reads(inc, isoBedObj, collapsedTrans2reads)
                exclusionIsos, exclusionReads = _event_reads(exc, isoBedObj, collapsedTrans2reads)
                if asType == "IR":
                    inclusionReads = _covering_intron(inclusionReads, asEvent, readBedObj)
                incTails = [[x, flncReads2Palen[x]] for x in inclusionReads if x in flncReads2Palen]
                excTails = [[x, flncReads2Palen[x]] for x in exclusionReads if x in flncReads2Palen]
                if len(incTails) < filterByCount or len(excTails) < filterByCount:
                    continue
                aPalen, bPalen = [float(x[1]) for x in incTails], [float(x[1]) for x in excTails]
                if _tails_differ(aPalen, bPalen):
                    count += 1
                    _write_palen_as(sigOut, gene, asType, count, [(inclusionIsos, incTails, aPalen),
                                                                  (exclusionIsos, excTails, bPalen)], isoBedObj)


def getPalenAS(flncReads2Palen, isoformFile, readsFile, collapsedTrans2reads=None, asPairs=None, filterByCount=10, mergeByJunc=False):
    prevDir = os.getcwd()
    resolveDir("mergeByJunc" if mergeByJunc else "no_mergeByJunc")
    isoBedObj = BedFile(isoformFile, type="bed12+")
    readBedObj = BedFile(readsFile, type="bed12+")
    for asType, pairsByGene in asPairs.items():
        _palen_as_type(asType, pairsByGene, flncReads2Palen, isoBedObj, readBedObj, collapsedTrans2reads, filterByCount)
    os.chdir(prevDir)
