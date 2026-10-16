"""AS characterization — ported from legacy ``charaterize_as.py``.

  * ``getAnnoASList`` / ``getASstatistics`` / ``getSpliceSite`` — file-level helpers
  * ``build_anno_as_lists`` — annotated / confident / Common / LR-specific event lists (driver l.257-302)
  * ``classify_reference_as`` — AS events implied by the reference annotation itself (l.304-334)
  * ``splicesite2seq`` / ``drawSSmotif`` / ``splicesite2figure`` — splice-site dinucleotide profiles + bar plots
  * ``getDist2TTS`` — distance of long-read PA sites to annotated TTS
  * ``characterize_sample`` — the full driver over ``as_events/ordinary_as`` (+ ``as_events/pa`` when present)

Perl helpers run in the core env, R plots in the R env via ``rscript.run_r_script`` (non-fatal). Two problems
found during the port are fixed (BUG-117, BUG-118; see CHANGELOG.md).
"""

from __future__ import annotations

import os
from collections import Counter

import pybedtools
import pysam

from ..io.files import filterFile, makeLink, resolveDir
from ..io.genepred import GenePredObj
from ..io.intervals import Bed12
from ..rscript import run_r_script
from ..runtime import ToolRunner, getCurrentTime, quote as q
from .find_as import _find_utils_dir, _run_redirect, find_single_as, getAnnotation


def getAnnoASList(inFile, outFile, PA=False, append=False, uniq=False):
    """Event keys of a known-event file: ``chrom:name:strand`` of each line without "Novel", or (``PA``)
    ``chrom:start:end:sites`` of each line with more than one PA site; ``uniq`` drops repeated keys.
    BUG-117: legacy opened the file with "w+" in append mode, truncating the known A3SS/A5SS events."""
    seen = set()
    with open(inFile) as f, open(outFile, "a" if append else "w") as out:
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            if PA:
                if len(lineInfo[3].split(",")) <= 1:
                    continue
                asEvent = ":".join(lineInfo[0:4])
            else:
                if "Novel" in line:
                    continue
                asEvent = "{}:{}:{}".format(lineInfo[0], lineInfo[3], lineInfo[5])
            if uniq and asEvent in seen:
                continue
            seen.add(asEvent)
            print(asEvent, file=out)


def getASstatistics(asType="IR", asFile=None, annoFile=None, novelFile=None, outFile=None):
    out = open(outFile, "w")
    if asType == "PA":
        print("#Chr\tStrand\tKnown or Novel\tGene\tPA Sites", file=out)
        with open(annoFile) as f1:
            for line in f1:
                lineInfo = line.strip("\n").split("\t")
                if len(lineInfo[3].split(",")) > 1:
                    print("\t".join([lineInfo[0], lineInfo[1], "Known", lineInfo[2], lineInfo[3]]), file=out)
        with open(novelFile) as f2:
            for line in f2:
                lineInfo = line.strip("\n").split("\t")
                if len(lineInfo[3].split(",")) > 1:
                    print("\t".join([lineInfo[0], lineInfo[1], "Novel", lineInfo[2], lineInfo[3]]), file=out)
    else:
        asDict = {}
        with open(asFile) as f:
            for line in f:
                lineInfo = line.strip("\n").split("\t")
                key = ":".join([lineInfo[0], lineInfo[3], lineInfo[5]])
                asDict[key] = "\t".join([key, line.strip("\n")])
        annoList = filterFile(originFile=annoFile, targetFile=asDict, returnFlag=True)
        novelList = filterFile(originFile=novelFile, targetFile=asDict, returnFlag=True)
        os.remove("filtered.txt")
        geneCol = 7
        if asType == "IR":
            print("##AS ID is composed of Gene:Retained intron start-Retained intron end", file=out)
        elif asType == "SE":
            geneCol = 13
            print("##AS ID is composed of Gene:Left flanking constitutive exon end@Alternative exons locus@Right flanking constitutive exon start", file=out)
            print("##Alternative exons locus is composed of Alternative exon1 start-Alternative exon1 end[;Alternative exon2 start-Alternative exon2 end[;Alternative exon3 start-Alternative exon3 end...]", file=out)
        elif asType == "A3SS":
            print("##AS ID is composed of Gene:Alternative 5' splicing region start-Alternative 5' splicing region end", file=out)
        elif asType == "A5SS":
            print("##AS ID is composed of Gene:Alternative 3' splicing region start-Alternative 3' splicing region end", file=out)
        print("#Chr\tStrand\tKnown or Novel\tAS ID\tGene", file=out)

        for i in annoList:
            tmpList = asDict[i].strip().split("\t")
            print("\t".join([tmpList[1], tmpList[6], "Known", tmpList[4], tmpList[geneCol]]), file=out)
        for i in novelList:
            tmpList = asDict[i].strip().split("\t")
            print("\t".join([tmpList[1], tmpList[6], "Novel", tmpList[4], tmpList[geneCol]]), file=out)

    out.close()


def getSpliceSite(asType=None, asFile=None, outFile=None, runner=None):
    if asType == "IR":
        out = open(outFile, "w")
        with open(asFile) as f:
            for line in f:
                lineInfo = line.strip("\n").split(":")
                posList = lineInfo[2].split("-")
                print("\t".join([lineInfo[0], posList[0], posList[1], lineInfo[3]]), file=out)
        out.close()
    elif asType == "SE":
        utilDir = _find_utils_dir()
        _run_redirect(["perl", os.path.join(utilDir, "seDecompose.pl"), asFile or "confident.SE.lst"],
                      "SE.inc.splicesite", "SE.exc.splicesite", runner=runner)
    elif asType == "A3SS":
        utilDir = _find_utils_dir()
        _run_redirect(["perl", os.path.join(utilDir, "anssDecompose.pl"), "-n", "3", asFile or "confident.A3SS.lst"],
                      "A3SS.inc.splicesite", "A3SS.exc.splicesite", runner=runner)
    elif asType == "A5SS":
        utilDir = _find_utils_dir()
        _run_redirect(["perl", os.path.join(utilDir, "anssDecompose.pl"), "-n", "5", asFile or "confident.A5SS.lst"],
                      "A5SS.inc.splicesite", "A5SS.exc.splicesite", runner=runner)


def build_anno_as_lists(ordinaryAsDir, outDir, pa_grouped=None, ngs=False, runner=None):
    """Assemble the annotated AS/APA event lists + ``confident.*`` lists.

    Deterministic front-half of the legacy ``charaterize_as`` driver
    (``charaterize_as.py`` lines 257-302), ported line-by-line.  ``ordinaryAsDir``
    is the ``ordinary_as`` directory produced by ``find_as`` (containing ``LR/``
    and optional ``NGS/`` sub-trees); results are written into ``outDir`` (the
    ``as_events/characterization`` directory in the full layout).

    Two file-level adaptations of the otherwise verbatim port, consistent with
    the rest of the "file-level core first" migration:

      * the poly(A) ``LR.APA.lst`` is only built when ``pa_grouped`` is given
        (legacy reads a fixed ``../pa/reads.paGrouped.tsv``, a deferred real-data
        upstream);
      * the NGS Venn / Common / LR-specific branch runs only when ``ngs=True``.
    """
    prevDir = os.getcwd()
    utilDir = _find_utils_dir()
    lrDir = os.path.join(ordinaryAsDir, "LR")
    ngsDir = os.path.join(ordinaryAsDir, "NGS")
    resolveDir(outDir)
    try:
        getAnnoASList(os.path.join(lrDir, "SE.bed12+"), "LR.SE.lst")
        getAnnoASList(os.path.join(lrDir, "A5SS.bed6+"), "LR.A5SS.lst")
        getAnnoASList(os.path.join(lrDir, "A3SS.bed6+"), "LR.A3SS.lst")
        getAnnoASList(os.path.join(lrDir, "IR.bed6+"), "LR.IR.lst")
        if pa_grouped:
            getAnnoASList(pa_grouped, "LR.APA.lst", PA=True)

        if ngs:
            getAnnoASList(os.path.join(ngsDir, "SE.bed12+"), "NGS.SE.lst")
            getAnnoASList(os.path.join(ngsDir, "A5SS.known.bed6+"), "NGS.A5SS.lst")
            getAnnoASList(os.path.join(ngsDir, "A5SS.LR.bed6+"), "NGS.A5SS.lst", append=True, uniq=True)
            getAnnoASList(os.path.join(ngsDir, "A3SS.known.bed6+"), "NGS.A3SS.lst")
            getAnnoASList(os.path.join(ngsDir, "A3SS.LR.bed6+"), "NGS.A3SS.lst", append=True, uniq=True)

            filterFile(originFile="NGS.SE.lst", targetFile="LR.SE.lst", outFile="Common.SE.lst")
            filterFile(originFile="NGS.A5SS.lst", targetFile="LR.A5SS.lst", outFile="Common.A5SS.lst")
            filterFile(originFile="NGS.A3SS.lst", targetFile="LR.A3SS.lst", outFile="Common.A3SS.lst")

            makeLink("Common.SE.lst", "confident.SE.lst")
            makeLink("Common.A5SS.lst", "confident.A5SS.lst")
            makeLink("Common.A3SS.lst", "confident.A3SS.lst")

            for asType in ("SE", "A5SS", "A3SS"):
                run_r_script("Rscript {u}/venn.R NGS.{t}.lst LR.{t}.lst -rD=90 -cP=180,0 -p={t}.pdf".format(u=q(utilDir), t=asType),
                             cwd=outDir, runner=runner)

            filterFile(originFile="NGS.SE.lst", targetFile="LR.SE.lst", outFile="LR.specific.SE.lst", mode="e")
            filterFile(originFile="NGS.A5SS.lst", targetFile="LR.A5SS.lst", outFile="LR.specific.A5SS.lst", mode="e")
            filterFile(originFile="NGS.A3SS.lst", targetFile="LR.A3SS.lst", outFile="LR.specific.A3SS.lst", mode="e")
        else:
            makeLink("LR.SE.lst", "confident.SE.lst")
            makeLink("LR.A5SS.lst", "confident.A5SS.lst")
            makeLink("LR.A3SS.lst", "confident.A3SS.lst")
    finally:
        os.chdir(prevDir)


def classify_reference_as(ref_gpe):
    """AS events present in the reference annotation itself (legacy driver l.304-334): the reference
    transcripts are intersected with themselves and fed through find_as's annotation/AS detectors ->
    ``<type>.reference.bed6+`` / ``<type>.reference.lst`` plus ``reference.bed12+`` (for ``paGroup.pl``)."""
    # legacy: awk '{print $0"\t"$12}' ref.gpe | gpe2bed.pl -p ; GenePredObj.toBed(gene=True) yields the same bed12+
    transBedList = GenePredObj(ref_gpe, bincolumn=False).toBed(gene=True)
    with open("reference.bed12+", "w") as f:
        f.write("".join(line + "\n" for line in transBedList))
    transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
    readsBedObj = pybedtools.BedTool("\n".join("\t".join(Bed12(line).record[:12]) for line in transBedList), from_string=True)
    annoBedRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
    novelBedRes = readsBedObj.intersect(transBedObj, v=True, s=True)
    annoDict, novelDict = getAnnotation(annoBedRes, novelBedRes, offset=0)
    for asType in ("IR", "SE", "A5SS", "A3SS"):
        out = "{}.reference.bed6+".format(asType)
        find_single_as(asType=asType, gene2ReadsDict=annoDict, anno=True, outFile=out, outAppend=False)
        find_single_as(asType=asType, gene2ReadsDict=novelDict, anno=False, outFile=out, outAppend=True)
        getAnnoASList(out, "{}.reference.lst".format(asType))


def splicesite2seq(refFasta, spliceSite, noFiveEnd=False, noThreeEnd=False, outFile=None):
    """Count the donor/acceptor dinucleotide pairs (``GT-AG`` ...) of the junctions in ``spliceSite``
    (chrom, left, right, strand); legacy port using pybedtools getfasta."""
    donorList, acceptorList = [], []
    exonDist, intronDist = 0, 2
    with open(spliceSite) as f:
        for count, line in enumerate(l for l in f.readlines() if l.strip()):
            chrom, leftSite, rightSite, strand = line.strip("\n").split("\t")[:4]
            juncName = "junction_" + str(count + 1)
            # anssDecompose.pl writes "dummy" on the side that is not a splice site: convert only what is used
            donorSite, acceptorSite = (leftSite, rightSite) if strand == "+" else (rightSite, leftSite)
            if not noFiveEnd:
                d = int(donorSite)
                donor = [chrom, d - exonDist, d + intronDist] if strand == "+" else [chrom, d - intronDist, d + exonDist]
                donorList.append("\t".join(map(str, donor + [juncName, ".", strand])))
            if not noThreeEnd:
                a = int(acceptorSite)
                acceptor = [chrom, a - intronDist, a + exonDist] if strand == "+" else [chrom, a - exonDist, a + intronDist]
                acceptorList.append("\t".join(map(str, acceptor + [juncName, ".", strand])))

    def fasta(bedList):
        # windows beyond a chromosome end make getfasta warn on stderr, which pybedtools raises as BEDToolsError
        with pysam.FastaFile(refFasta) as fa:
            sizes = dict(zip(fa.references, fa.lengths))
        kept = [b for b in bedList if int(b.split("\t")[1]) >= 0 and int(b.split("\t")[2]) <= sizes.get(b.split("\t")[0], 0)]
        if not kept:
            return {}
        res = pybedtools.BedTool("\n".join(kept), from_string=True).sequence(refFasta, name=True, tab=True, s=True)
        out = {}
        for line in open(res.seqfn).read().split("\n")[:-1]:
            name, seq = line.split("\t")
            out[name.split(":")[0]] = seq
        return out

    donors, acceptors = fasta(donorList) if not noFiveEnd else {}, fasta(acceptorList) if not noThreeEnd else {}
    if not noFiveEnd and not noThreeEnd:
        junc2seq = {j: donors[j] + "-" + acceptors[j] for j in donors if j in acceptors}
    else:
        junc2seq = donors or acceptors
    d_tmp = sorted(((v, k) for k, v in Counter(junc2seq.values()).items()), reverse=True)
    with open(outFile, "w") as out:
        for v, k in d_tmp:
            print("\t".join(map(str, [k, v])), file=out)


def drawSSmotif(asMotif=None, outPrefix=None, asType="IR", runner=None):
    """Top-3 splice-site motifs (+ ``Other``) as fractions -> ``<outPrefix>.tmp.txt`` and a ``bar.R`` plot."""
    with open(asMotif) as f:
        lineList = [l for l in f.readlines() if l.strip()]
    mySum = sum(int(i.strip("\n").split("\t")[1]) for i in lineList)

    def label(motif):
        return "-{}".format(motif) if asType == "A3SS" else "{}-".format(motif) if asType == "A5SS" else motif

    with open("{}.tmp.txt".format(outPrefix), "w") as tmp:
        for i in lineList[:3]:
            lineInfo = i.strip("\n").split("\t")
            print("\t".join([label(lineInfo[0]), str(float(lineInfo[1]) / mySum)]), file=tmp)
        if len(lineList) >= 4:
            otherSum = sum(int(i.strip("\n").split("\t")[1]) for i in lineList[3:])
            print("\t".join(["Other", str(float(otherSum) / mySum)]), file=tmp)
    if lineList:            # no event of this class: bar.R fails on the empty table and logs a misleading warning (BUG-195)
        run_r_script("cat {p}.tmp.txt | Rscript {u}/bar.R -fillV=V1 -fp -lgPos=top -w=10 -p={p}.ssMotif.pdf".format(
            p=outPrefix, u=q(_find_utils_dir())), cwd=os.getcwd(), runner=runner)


def splicesite2figure(refGenome, runner=None):
    """Splice-site dinucleotide profiles + plots for every decomposed event class (legacy l.230-247)."""
    for prefix, kwargs, asType in (("IR", {}, "IR"), ("IR.anno", {}, "IR"), ("IR.novel", {}, "IR"),
                                   ("SE.inc", {}, "SE"), ("SE.exc", {}, "SE"),
                                   ("A5SS.inc", {"noThreeEnd": True}, "A5SS"), ("A5SS.exc", {"noThreeEnd": True}, "A5SS"),
                                   ("A3SS.inc", {"noFiveEnd": True}, "A3SS"), ("A3SS.exc", {"noFiveEnd": True}, "A3SS")):
        splicesite2seq(refGenome, prefix + ".splicesite", outFile=prefix + ".tmp", **kwargs)
        drawSSmotif(asMotif=prefix + ".tmp", outPrefix=prefix, asType=asType, runner=runner)


def getDist2TTS(refGpe, paGroup, runner=None):
    """Distance from each long-read PA site (``reads.paGrouped.bed6``) to the nearest annotated TTS ->
    ``lrPA2TTS.tsv`` + ``box.R`` plot (legacy l.130-148)."""
    util = _find_utils_dir()
    with open(paGroup) as f, open("lrPA.bed6", "w") as out:
        for line in f:
            c = line.strip("\n").split("\t")
            if c[5] == "+":
                print("\t".join(map(str, [c[0], int(c[2]) - 1, c[2]] + c[3:])), file=out)
            else:
                print("\t".join(map(str, [c[0], c[1], int(c[1]) + 1] + c[3:])), file=out)
    runner.run(["bash", "-c", "bedtools closest -a <(sort -k1,1 -k2,2n lrPA.bed6) -b <(perl {u}/gpeFeature.pl --tts {g} | "
                              "sort -k1,1 -k2,2n) -s -D a | perl {u}/select.pl -i 13,4 | sort -u > lrPA2TTS.tsv".format(u=q(util), g=q(refGpe))])
    run_r_script("cut -f1 lrPA2TTS.tsv | Rscript {u}/box.R -w=10 -height=6 -ng -nJ -no -y='Distance to TTS' -p=lrPA2TTS.pdf".format(u=q(util)),
                 cwd=os.getcwd(), runner=runner)


def characterize_sample(base_dir, out_dir=None, pa_grouped=None, ngs=None, splice_site=True, runner=None, *,
                        ref_gpe=None, ref_genome=None):
    """Legacy ``charaterize_as`` driver over ``<base_dir>/as_events/ordinary_as`` (find_as) and, when present,
    ``<base_dir>/as_events/pa`` (find_pa), writing ``<base_dir>/as_events/characterization``.

    Front half (always): annotated / confident / Common / LR-specific event lists (NGS branch auto-detected
    from ``ordinary_as/NGS``), ``LR.APA.lst`` when PA groups exist. Back half (needs ``ref_gpe``): reference AS
    classification -> anno/novel split -> ``statistic.*.tsv``; ``paGroup.pl``/``paCmp.pl`` APA comparison and
    ``getDist2TTS`` when PA groups exist; splice-site tables and (with ``ref_genome``) dinucleotide profiles.
    """
    runner = runner or ToolRunner()
    base_dir = os.path.abspath(base_dir)
    ordinaryAsDir = os.path.join(base_dir, "as_events", "ordinary_as")
    if not os.path.isdir(os.path.join(ordinaryAsDir, "LR")):
        raise FileNotFoundError(
            "ordinary_as/LR not found under {} (run find_as first)".format(ordinaryAsDir))
    out_dir = os.path.abspath(out_dir) if out_dir else os.path.join(base_dir, "as_events", "characterization")
    pa_dir = os.path.join(base_dir, "as_events", "pa")
    if pa_grouped is None and os.path.isfile(os.path.join(pa_dir, "reads.paGrouped.tsv")):
        pa_grouped = os.path.join(pa_dir, "reads.paGrouped.tsv")
    pa_grouped = os.path.abspath(pa_grouped) if pa_grouped else None
    if ngs is None:
        ngs = os.path.isdir(os.path.join(ordinaryAsDir, "NGS"))
    ref_gpe = os.path.abspath(ref_gpe) if ref_gpe else None
    ref_genome = os.path.abspath(ref_genome) if ref_genome else None
    lr = os.path.join(ordinaryAsDir, "LR")
    util = _find_utils_dir()

    print(getCurrentTime() + " ASE Characterization...")
    build_anno_as_lists(ordinaryAsDir, out_dir, pa_grouped=pa_grouped, ngs=ngs, runner=runner)

    prevDir = os.getcwd()
    os.chdir(out_dir)
    try:
        if ref_gpe:
            classify_reference_as(ref_gpe)
            if pa_grouped:
                runner.run("perl {u}/paGroup.pl reference.bed12+ > paGroup.reference.tsv 2> paGroup.reference.bed6".format(u=q(util)))
                runner.run("perl {u}/paCmp.pl -r paGroup.reference.tsv -a annoPA.tsv -n novelPA.tsv {pa} "
                           "> paGroup.novel.tsv 2> paGroup.anno.tsv".format(u=q(util), pa=q(pa_grouped)))
            filterFile(originFile="IR.reference.lst", targetFile="LR.IR.lst", outFile="IR.anno.lst", mode="i")
            filterFile(originFile="IR.reference.lst", targetFile="LR.IR.lst", outFile="IR.novel.lst", mode="e")
            for asType in ("SE", "A3SS", "A5SS"):
                filterFile(originFile="{}.reference.lst".format(asType), targetFile="confident.{}.lst".format(asType),
                           outFile="{}.anno.lst".format(asType), mode="i")
                filterFile(originFile="{}.reference.lst".format(asType), targetFile="confident.{}.lst".format(asType),
                           outFile="{}.novel.lst".format(asType), mode="e")
            getASstatistics(asType="IR", asFile=os.path.join(lr, "IR.bed6+"), annoFile="IR.anno.lst", novelFile="IR.novel.lst",
                            outFile="statistic.IR.tsv")
            if pa_grouped:
                getASstatistics(asType="PA", annoFile="paGroup.anno.tsv", novelFile="paGroup.novel.tsv", outFile="statistic.APA.tsv")
            for asType, fn in (("SE", "SE.confident.bed12+"), ("A3SS", "A3SS.confident.bed6+"), ("A5SS", "A5SS.confident.bed6+")):
                getASstatistics(asType=asType, asFile=os.path.join(lr, fn), annoFile="{}.anno.lst".format(asType),
                                novelFile="{}.novel.lst".format(asType), outFile="statistic.{}.tsv".format(asType))

        spliceSites = {}
        if splice_site:
            getSpliceSite(asType="IR", asFile="LR.IR.lst", outFile="IR.splicesite")
            if ref_gpe:
                getSpliceSite(asType="IR", asFile="IR.anno.lst", outFile="IR.anno.splicesite")
                getSpliceSite(asType="IR", asFile="IR.novel.lst", outFile="IR.novel.splicesite")
            getSpliceSite(asType="SE", asFile="confident.SE.lst", runner=runner)
            getSpliceSite(asType="A3SS", asFile="confident.A3SS.lst", runner=runner)
            getSpliceSite(asType="A5SS", asFile="confident.A5SS.lst", runner=runner)
            spliceSites = {
                "IR": os.path.join(out_dir, "IR.splicesite"),
                "SE": (os.path.join(out_dir, "SE.inc.splicesite"), os.path.join(out_dir, "SE.exc.splicesite")),
                "A5SS": (os.path.join(out_dir, "A5SS.inc.splicesite"), os.path.join(out_dir, "A5SS.exc.splicesite")),
                "A3SS": (os.path.join(out_dir, "A3SS.inc.splicesite"), os.path.join(out_dir, "A3SS.exc.splicesite")),
            }
            if ref_gpe and ref_genome:
                splicesite2figure(ref_genome, runner=runner)
        if ref_gpe and pa_grouped:
            getDist2TTS(ref_gpe, os.path.join(os.path.dirname(pa_grouped), "reads.paGrouped.bed6"), runner=runner)
    finally:
        os.chdir(prevDir)
    print(getCurrentTime() + " ASE Characterization done!")

    return {
        "out_dir": out_dir,
        "confident": {asType: os.path.join(out_dir, "confident.{}.lst".format(asType)) for asType in ("SE", "A5SS", "A3SS")},
        "splice_site": spliceSites,
        "statistics": {asType: os.path.join(out_dir, "statistic.{}.tsv".format(asType))
                       for asType in ("IR", "SE", "A3SS", "A5SS") + (("APA",) if pa_grouped else ())} if ref_gpe else {},
    }
