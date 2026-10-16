"""find_pa stage - poly(A) site clustering, ported from legacy find_pa.py.

``getPAC`` / ``paCluster`` / ``getPaCluster`` are the in-memory PAC clustering core;
``motifAroundPA`` the nucleotide / PAS-hexamer profiles around PA sites; ``find_pa_sample`` the
sample-directory driver (legacy ``find_pa()``), with the Perl helpers run in the core env and the
R plots in ``iflas-diffsplice``.
"""

from __future__ import annotations

import copy
import itertools
import os
from collections import Counter

import numpy as np
import pandas as pd
import pybedtools

from ..intervalops import isOverlap
from ..io.files import makeLink, resolveDir, validateFile
from ..io.intervals import Bed6Plus
from ..io.records import BedFile
from ..pu.cpc2 import _find_utils_dir
from ..rscript import run_r_script
from ..runtime import ToolRunner, getCurrentTime, quote as q

PAS_HEXAMERS = ["AATAAA", "AAATAA", "ATAAAA", "ATTAAA", "ATAAAT", "TAATAA",
                "ATAAAG", "AAAATA", "CAATAA", "ATAAAC", "AAAAAA", "AAAAAG"]
OTHER_PAS = ["ATAAAC.PAS", "ATAAAG.PAS", "CAATAA.PAS", "AAAAAA.PAS", "AAAAAG.PAS", "AAAATA.PAS"]


def getPAC(depth, distance=20, windowSize=2, paSup=5, rpkmPAC=0, allReadCount=1000000, effLen=1):
    N = len(depth)
    currPeaks = np.zeros(len(depth))
    for i in range(N):
        if np.sum(depth[max(0, i - windowSize):min(N, i + windowSize + 1)]) >= paSup:
            currPeaks[i] = depth[i] * 2 + np.mean(depth[max(0, i - windowSize):min(N, i + windowSize + 1)])
    currPeaksBak = copy.copy(currPeaks)
    pacDict = {}
    count = 0
    while np.max(currPeaks):
        cp = np.argmax(currPeaks)
        overlapFlag = 0
        maxDiff = 1000000
        closestPAC = ""
        for pac in pacDict:
            if isOverlap((cp-distance, cp+distance+1), (min(pacDict[pac]), max(pacDict[pac]))):
                overlapFlag = 1
                diff2pac = min(abs(cp - min(pacDict[pac])), abs(cp - max(pacDict[pac])))
                if diff2pac < maxDiff:
                    maxDiff = diff2pac
                    closestPAC = pac
        if overlapFlag:
            # PAC ranges are consumed as half-open slices below; include the
            # closed peak position and its one-past edge when extending a PAC.
            pacDict[closestPAC].extend([cp, min(cp + 1, N)])
        else:
            count += 1
            pacName = "PAC_{}".format(count)
            pacDict[pacName] = [cp, max([0, cp-windowSize]), min([cp+windowSize+1, N])]
        currPeaks[cp] = 0

    if len(pacDict) == 0: return []
    newPacDict = {list(pacDict.keys())[0]: copy.copy(pacDict[list(pacDict.keys())[0]])}
    for pac1 in list(pacDict.keys())[1:]:
        pac1Start = min(pacDict[pac1])
        pac1End = max(pacDict[pac1])
        for pac2 in newPacDict:
            pac2Start = min(newPacDict[pac2])
            pac2End = max(newPacDict[pac2])
            if isOverlap((pac1Start-distance, pac1End+distance+1), (pac2Start, pac2End)):
                newPacDict[pac2].extend(pacDict[pac1])
                break
        else:
            newPacDict[pac1] = copy.copy(pacDict[pac1])

    finalPeak2Count = []
    for pac in newPacDict:
        pacStart = min(newPacDict[pac])
        pacEnd = max(newPacDict[pac])
        peakPosInPac = int(np.argmax(currPeaksBak[pacStart: pacEnd]))
        peak = list(range(pacStart, pacEnd))[peakPosInPac]
        count = sum(depth[pacStart:pacEnd])
        if rpkmPAC:
            rpkm = (count * 1000000) / float(allReadCount*effLen)
            if rpkm >= rpkmPAC:
                finalPeak2Count.append((peak, count, [x for x in range(pacStart, pacEnd) if depth[x] > 0]))
        else:
            finalPeak2Count.append((peak, count, [x for x in range(pacStart, pacEnd) if depth[x] > 0]))
    if finalPeak2Count:
        finalPeak2Count = sorted(finalPeak2Count, key=lambda pair: pair[1], reverse=True)
    return finalPeak2Count


def paCluster(readsBedList, distance=20, windowSize=3, outHandle=None, allReadCount=1000000, paSup=5,
              rpkmPAC=0, confidentPaDict=None):
    chroms = [i.chrom for i in readsBedList]
    strands = [i.strand for i in readsBedList]
    refGene = Counter([i.otherList[2] for i in readsBedList]).most_common()[0][0]
    if len(set(strands)) != 1: return
    chrom = list(set(chroms))[0]
    strand = list(set(strands))[0]
    if strand == "+":
        sortedPAs = sorted(readsBedList, key=lambda x: x.chromEnd)
        paRange = [sortedPAs[0].chromEnd, sortedPAs[-1].chromEnd]
        effectiveLength = abs(paRange[0] - paRange[1]) + 1
        depth = np.zeros(paRange[1] - paRange[0] + 1, int)
        depth2reads = {}
        offest = paRange[0]
        for read in sortedPAs:
            depth[read.chromEnd-offest] += 1
            if read.chromEnd-offest not in depth2reads:
                depth2reads[read.chromEnd-offest] = []
            depth2reads[read.chromEnd-offest].append(read.name)
        finalPeak2Count = getPAC(depth, distance=distance, windowSize=windowSize, paSup=paSup, rpkmPAC=rpkmPAC, allReadCount=allReadCount, effLen=effectiveLength)
    else:
        sortedPAs = sorted(readsBedList, key=lambda x: x.chromStart)
        paRange = [sortedPAs[0].chromStart, sortedPAs[-1].chromStart]
        effectiveLength = abs(paRange[0] - paRange[1]) + 1
        depth = np.zeros(paRange[1] - paRange[0] + 1, int)
        depth2reads = {}
        offest = paRange[0]
        for read in sortedPAs:
            depth[read.chromStart-offest] += 1
            if read.chromStart-offest not in depth2reads:
                depth2reads[read.chromStart-offest] = []
            depth2reads[read.chromStart-offest].append(read.name)
        finalPeak2Count = getPAC(depth, distance=distance, windowSize=windowSize, paSup=paSup, rpkmPAC=rpkmPAC, allReadCount=allReadCount, effLen=effectiveLength)
    for peak, count, pac in finalPeak2Count:
        currPaSite = "{}_{}".format(chrom, peak+offest)
        if confidentPaDict:
            if currPaSite in confidentPaDict:
                annotation = "Known"
            else:
                annotation = "Novel"
        else:
            annotation = "Unknown"
        readNames = list(itertools.chain.from_iterable([depth2reads[x] for x in range(min(pac), max(pac)+1) if x in depth2reads]))
        print("\t".join(map(str, [chrom, offest + min(pac), offest + max(pac), ",".join(readNames),
                                                len(readNames), strand, offest + peak-1, offest + peak, refGene,
                                                annotation])), file=outHandle)


def getPaCluster(readsBed=None, tofuGroup=None, filterByCount=0, threads=None, paClusterOut=None, paDist=20,
                 windowSize=5, rpkmPAC=0, confidentPa=None):
    gene2reads = {}
    readsDict = BedFile(readsBed, type="bed12+").reads
    allReadCount = len(readsDict)
    with open(tofuGroup) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            gene = ".".join(infoList[0].split(".")[:2])
            if gene not in gene2reads:
                gene2reads[gene] = {"allReads": []}
            gene2reads[gene]["allReads"].extend(infoList[1].split(","))
    confidentPaDict = {}
    if confidentPa and validateFile(confidentPa):
        with open(confidentPa) as f:
            for i in f.readlines():
                infoList = i.strip("\n").split("\t")
                for x in range(max(0, int(infoList[2])-paDist-1), int(infoList[2])+paDist):
                    paSite = "{}_{}".format(infoList[0], x)
                    confidentPaDict.update({paSite: ""})
    outHandle = open(paClusterOut, "w")
    for z in gene2reads:
        if filterByCount:
            if len(gene2reads[z]["allReads"]) < filterByCount:
                continue
        readsList = gene2reads[z]["allReads"]
        readsBedList = [readsDict[r] for r in readsList if r in readsDict]
        if len(readsBedList) == 0: continue
        chroms = [r.chrom for r in readsBedList]
        mostChrom = Counter(chroms).most_common(1)[0][0]
        readsBedList = [r for r in readsBedList if r.chrom == mostChrom]
        paCluster(readsBedList, distance=paDist, windowSize=windowSize, outHandle=outHandle,
                       allReadCount=allReadCount, paSup=filterByCount, rpkmPAC=rpkmPAC, confidentPaDict=confidentPaDict)
    outHandle.close()


def motifAroundPA(bed6plus=None, up1=100, down1=100, up2=100, down2=100, refFasta=None, chrLenFile=None):
    """Per-position nucleotide (``<base>.nucleotide``) and PAS hexamer (``<hexamer>.PAS``) frequencies
    around the PA sites of ``bed6plus`` (legacy port; writes into the current directory)."""
    chrLenDict = {}
    with open(chrLenFile) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            chrLenDict[infoList[0]] = int(infoList[1])
    with open(bed6plus) as f:
        singleNucleotideUpBedList = []
        singleNucleotideDownBedList = []
        sixNucleotideUpBedList = []
        sixNucleotideDownBedList = []
        for i in f.readlines():
            bedObj = Bed6Plus(i)
            chrom, start, end, strand = bedObj.chrom, bedObj.chromStart, bedObj.chromEnd, bedObj.strand
            name = bedObj.name
            if strand == "+":
                if start > up1 and end + down1 < chrLenDict[chrom]:
                    singleNucleotideUpBedList.append(" ".join(map(str, [chrom, start - up1, start, name, ".", strand])))
                    singleNucleotideDownBedList.append(" ".join(map(str, [chrom, end, end + down1, name, ".", strand])))
                if start > up2 and end + down2 + 5 < chrLenDict[chrom]:
                    sixNucleotideUpBedList.append(" ".join(map(str, [chrom, start - up2, start + 5, name, ".", strand])))
                    sixNucleotideDownBedList.append(" ".join(map(str, [chrom, end, end + down2 + 5, name, ".", strand])))
            else:
                if start > down1 and end + up1 < chrLenDict[chrom]:
                    singleNucleotideUpBedList.append(" ".join(map(str, [chrom, end, end + up1, name, ".", strand])))
                    downLeft = start - down1 - 1 if start - down1 - 1 > 0 else 0
                    singleNucleotideDownBedList.append(" ".join(map(str, [chrom, downLeft, start - 1, name, ".", strand])))
                if start > down2 and end + up2 + 5 < chrLenDict[chrom]:
                    sixNucleotideUpBedList.append(" ".join(map(str, [chrom, end, end + up2 + 5, name, ".", strand])))
                    downLeft = start - down2 - 5 if start - down2 - 5 > 0 else 0
                    sixNucleotideDownBedList.append(" ".join(map(str, [chrom, downLeft, start, name, ".", strand])))

        def seqs(bedList, k):
            """name -> list of k-mers (k=1: bases) of the strand-aware sequence of each interval"""
            res = pybedtools.BedTool("\n".join(bedList), from_string=True).sequence(refFasta, name=True, tab=True, s=True)
            out = {}
            for line in open(res.seqfn).read().split("\n")[:-1]:
                name, seq = line.split("\t")
                out[name] = [seq.upper()[i:i + k] for i in range(0, len(seq) - k + 1)]
            return pd.DataFrame.from_dict(out, orient="index")

        singleNucleotideUpDf, singleNucleotideDownDf = seqs(singleNucleotideUpBedList, 1), seqs(singleNucleotideDownBedList, 1)
        sixNucleotideUpDf, sixNucleotideDownDf = seqs(sixNucleotideUpBedList, 6), seqs(sixNucleotideDownBedList, 6)

        def profile(outFile, motif, upDf, downDf, up):
            with open(outFile, "w") as out:
                upPercent = (upDf == motif).sum() / float(len(upDf))
                for j in range(len(upPercent)):
                    print("\t".join(map(str, [j - up, upPercent.iloc[j]])), file=out)
                downPercent = (downDf == motif).sum() / float(len(downDf))
                for j in range(len(downPercent)):
                    print("\t".join(map(str, [j + 1, downPercent.iloc[j]])), file=out)

        for i in ["A", "T", "C", "G"]:
            profile(i + ".nucleotide", i, singleNucleotideUpDf, singleNucleotideDownDf, up1)
        for i in PAS_HEXAMERS:
            profile(i + ".PAS", i, sixNucleotideUpDf, sixNucleotideDownDf, up2)


def write_cleavage_sites(reads_bed, out_file):
    """Unique read 3' ends (cleavage sites) of ``reads_bed`` as bed6 (legacy find_pa l.271-290)."""
    seen, count = set(), 0
    with open(reads_bed) as f, open(out_file, "w") as out:
        for line in f:
            c = line.rstrip("\n").split("\t")
            chrom, start, end, strand = c[0], int(c[1]), int(c[2]), c[5]
            if strand == "+":
                site = (chrom, end - 1, end, strand)
            elif strand == "-":
                site = (chrom, start, start + 1, strand)
            else:
                continue
            if site not in seen:
                seen.add(site)
                count += 1
                print("\t".join(map(str, [site[0], site[1], site[2], "ClevageSite{}".format(count), 0, strand])), file=out)


def find_pa_sample(base_dir, ref_genome, ref_size, *, filter_pa_by_rpkm=0, filter_pa_by_count=5,
                   confident_pa=None, threads=1, runner=None):
    """Legacy ``find_pa()`` in ``<base_dir>/as_events/pa``: cleavage sites -> PAC clustering (``paCluster.bed8+``)
    -> ``PA.bed6+`` -> isoform / read PA groups (``paGroup.pl``) -> 3'-end revision -> PA motif tables and
    plots (``motif/``). Perl helpers run in the core env, R plots in the R env (``run_r_script``)."""
    runner = runner or ToolRunner()
    base_dir = os.path.abspath(base_dir)
    ref_genome, ref_size = os.path.abspath(ref_genome), os.path.abspath(ref_size)
    util = q(_find_utils_dir())        # only used inside shell commands
    pa_dir = os.path.join(base_dir, "as_events", "pa")
    motif_dir = os.path.join(pa_dir, "motif")
    lines_r = "Rscript {}/lines.R".format(util)
    print(getCurrentTime() + " Polyadenylation analysis...")
    prev_dir = os.getcwd()
    resolveDir(pa_dir)
    try:
        for sub, name in (("collapse", "tofu.collapsed.group.txt"), ("refine", "reads.assigned.unambi.bed12+"),
                          ("refine", "isoformGrouped.bed12+")):
            makeLink(os.path.join(base_dir, sub, name), name)
        write_cleavage_sites("reads.assigned.unambi.bed12+", "cleavage.bed6")
        runner.run("perl {}/bedClosest.pl cleavage.bed6 > adjacentCleavage.tsv".format(util), cwd=pa_dir)
        run_r_script("cut -f13 adjacentCleavage.tsv | Rscript {}/distrCurve.R -w=10 -xl=2 -d "
                      "-x='Binned Distance between Adjacent Cleavage Site' -p=adjacentCleavageDistr.pdf".format(util), cwd=pa_dir, runner=runner)
        getPaCluster(readsBed="reads.assigned.unambi.bed12+", tofuGroup="tofu.collapsed.group.txt", threads=threads,
                     paClusterOut="paCluster.bed8+", paDist=24, windowSize=3, rpkmPAC=filter_pa_by_rpkm,
                     filterByCount=filter_pa_by_count, confidentPa=confident_pa)
        with open("paCluster.bed8+") as src, open("PA.bed6+", "w") as dst:
            for line in src:
                c = line.rstrip("\n").split("\t")
                if int(c[4]) >= filter_pa_by_count:
                    dst.write("\t".join([c[0], c[6], c[7]] + c[3:6] + c[9:]) + "\n")
        runner.run("perl {u}/paGroup.pl isoformGrouped.bed12+ > isoform.paGrouped.tsv 2> isoform.paGrouped.bed6".format(u=util), cwd=pa_dir)
        runner.run(["bash", "-c", "perl {u}/3endRevise.pl -p PA.bed6+ <(cut -f 1-12,15 reads.assigned.unambi.bed12+) "
                                  "> reads.3endRevised.bed12+".format(u=util)], cwd=pa_dir)
        runner.run("cut -f 4 PA.bed6+ | tr ',' '\\n' | perl {u}/filter.pl -o - reads.3endRevised.bed12+ -2 4 -m i | "
                   "perl {u}/paGroup.pl > reads.paGrouped.tsv 2> reads.paGrouped.bed6".format(u=util), cwd=pa_dir)

        if os.path.getsize("PA.bed6+") == 0:
            print(getCurrentTime() + " No PA cluster with >= {} reads: PA motif analysis skipped!".format(filter_pa_by_count))
            return _find_pa_outputs(pa_dir, motif_dir)
        resolveDir("motif")
        motifAroundPA("../PA.bed6+", up1=100, down1=100, up2=50, down2=0, refFasta=ref_genome, chrLenFile=ref_size)
        run_r_script("{r} -w=10 -y=Frequency -x='Distance Relative to PA' -m='Distribution of Nucleotide Frequency around PA' "
                      "-p=nucleotide.pdf *.nucleotide".format(r=lines_r), cwd=motif_dir, runner=runner)
        run_r_script("{r} -p=PAS.1.pdf {{AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}}.PAS -x1=-50 -x2=0 -w=10".format(r=lines_r), cwd=motif_dir, runner=runner)
        run_r_script("{r} -p=PAS.2.pdf {{ATAAAG,AAAATA,CAATAA,ATAAAC,AAAAAA,AAAAAG}}.PAS -x1=-50 -x2=0 -w=10".format(r=lines_r), cwd=motif_dir, runner=runner)
        other = pd.concat([pd.read_csv(f, sep="\t", index_col=0, header=None) for f in OTHER_PAS], axis=1).sum(axis=1)
        other.to_frame().to_csv("Other.PAS", sep="\t", header=None)
        run_r_script("{r} -p=PAS.pdf {{Other,AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}}.PAS -x1=-50 -x2=0 -w=10".format(r=lines_r), cwd=motif_dir, runner=runner)
    finally:
        os.chdir(prev_dir)
    print(getCurrentTime() + " Polyadenylation analysis done!")
    return _find_pa_outputs(pa_dir, motif_dir)


def _find_pa_outputs(pa_dir, motif_dir):
    return {"pa_dir": pa_dir,
            "pa_cluster": os.path.join(pa_dir, "paCluster.bed8+"),
            "pa_bed6plus": os.path.join(pa_dir, "PA.bed6+"),
            "reads_pa_grouped": os.path.join(pa_dir, "reads.paGrouped.tsv"),
            "reads_pa_grouped_bed6": os.path.join(pa_dir, "reads.paGrouped.bed6"),
            "isoform_pa_grouped": os.path.join(pa_dir, "isoform.paGrouped.tsv"),
            "reads_3end_revised": os.path.join(pa_dir, "reads.3endRevised.bed12+"),
            "motif_dir": motif_dir}
