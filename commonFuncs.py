#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: commonFuncs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

import os, re, datetime, psutil, subprocess
import pysam, pybedtools
from Bio import SeqIO

AStypes = ["IR", "SE", "A3SS", "A5SS"]


def batchCreateDir(dirList):
    for i in dirList:
        if not os.path.exists(i):
            os.makedirs(i)


def checkFast5Files(fast5Dir):
    if os.path.isfile(fast5Dir):
        if fast5Dir.endswith("fq") or fast5Dir.endswith("fastq"):
            return "fastq"
    else:
        fileCount = len([x for x in os.listdir(fast5Dir) if not x.endswith(".txt")])
        fast5FileCount = 0
        fqFileCount = 0
        for i in os.listdir(fast5Dir):
            if i.endswith("fast5"):
                fast5FileCount += 1
            if i.endswith("fq") or i.endswith("fastq"):
                fqFileCount += 1
        if fast5FileCount == fileCount:
            return "fast5"
        if fqFileCount == fileCount:
            return "fastq"
        if fast5FileCount != fileCount or fqFileCount != fileCount:
            raise Exception(
                "Please check the {} directory which is not a pure fast5 or fastq-containing directory".format(
                    fast5Dir))


def checkHisat2IndexExist(hisat2index):
    import glob
    hisat2indexFiles = glob.glob("{}*".format(hisat2index))

    if len(hisat2indexFiles) == 0:
        return False
    else:
        ht2Flag = 0
        for i in hisat2indexFiles:
            if i.endswith("exon") or i.endswith("ss") or i.endswith("log"):
                continue
            elif i.endswith("ht2"):
                ht2Flag = 1
                continue
            else:
                ht2Flag = 0
                break
        return True if ht2Flag == 1 else False


def checkFastpReadsExist(readStr):
    readsList = re.split("[,|;]", readStr)
    absReadsList = [os.path.join(os.path.dirname(x), "fastp." + os.path.basename(x)) for x in readsList]
    boolList = [validateFile(x) for x in absReadsList]
    if len(set(boolList)) == 1 and boolList[0] == True:
        return True
    else:
        return False


def filterFile(originFile=None, targetFile=None, originField=1, targetField=1, mode="i", outFile=None,
               returnFlag=False):
    originList, targetList = [], []
    if isinstance(originFile, list):
        originList = originFile
    elif os.path.isfile(originFile):
        with open(originFile) as f:
            originList = [line.strip("\n").split("\t")[originField - 1] for line in f.readlines()]

    returnList = []
    outFile = outFile if outFile else "filtered.txt"
    out = open(outFile, "w")
    if isinstance(targetFile, list):
        targetList = targetFile
        if mode == "i":
            returnList = list(set(originList) & set(targetList))
        else:
            returnList = list(set(targetList) - set(originList))

        for i in returnList:
            print >> out, i
    elif isinstance(targetFile, dict):
        if mode == "i":
            returnList = list(set(originList) & set(targetFile.keys()))
        else:
            returnList = list(set(targetFile.keys()) - set(originList))

        for i in returnList:
            print >> out, i
    elif os.path.isfile(targetFile):
        with open(targetFile) as f:
            targetDict = {}
            for line in f:
                lineInfo = line.strip("\n").split("\t")
                targetDict[lineInfo[targetField - 1]] = line.strip("\n")
            if mode == "i":
                returnList = list(set(originList) & set(targetDict.keys()))
            else:
                returnList = list(set(targetDict.keys()) - set(originList))
            for i in returnList:
                print >> out, targetDict[i]
    out.close()
    if returnFlag:
        return returnList


def makeLink(sourcePath, targetPath):
    if os.path.islink(targetPath) or os.path.exists(targetPath):
        os.remove(targetPath)
    os.symlink(sourcePath, targetPath)


def makeHisat2Index(refParams=None, dirSpec=None, threads=None):
    hisat2indexDir = os.path.join(dirSpec.out_dir, "hisat2indexDir")
    if not os.path.exists(hisat2indexDir):
        resolveDir(hisat2indexDir, chdir=False)
    refAnnoGTF = refParams.ref_gtf
    gtfPrefix = os.path.splitext(os.path.basename(refAnnoGTF))[0]
    hisat2index = os.path.join(hisat2indexDir, gtfPrefix)
    if not checkHisat2IndexExist(hisat2index):
        print getCurrentTime() + " Start Hisat2 indexing..."
        resolveDir(hisat2indexDir, chdir=False)

        refAnnoSS = "{}/{}.ss".format(hisat2indexDir, gtfPrefix)
        cmd = "hisat2_extract_splice_sites.py {} >{}".format(refAnnoGTF, refAnnoSS)
        subprocess.call(cmd, shell=True)
        refAnnoExon = "{}/{}.exon".format(hisat2indexDir, gtfPrefix)
        cmd = "hisat2_extract_exons.py {} >{}".format(refAnnoGTF, refAnnoExon)
        subprocess.call(cmd, shell=True)
        cmd = "hisat2-build -p {} --ss {} --exon {} {} {} 1>{}/hisat2build.log 2>&1"
        cmd = cmd.format(threads, refAnnoSS, refAnnoExon, refParams.ref_genome, hisat2index, hisat2indexDir)
        subprocess.call(cmd, shell=True)
        print getCurrentTime() + " End Hisat2 indexing!"
    else:
        print getCurrentTime() + " It seems like you have provided hisat2 index, so the index will be used for mapping."
    return hisat2index

def mergeSample(strain2data):
    sampleMergedToProcess = {}
    from commonObjs import MergedTgsSample
    for proj in strain2data:
        for ref_strain in strain2data[proj]:
            for strain in strain2data[proj][ref_strain]:
                dataObjs = strain2data[proj][ref_strain][strain]
                sampleMerged = MergedTgsSample()
                sampleMerged.project_name = proj
                sampleMerged.sample_name = "{}_all".format(strain)
                sampleMerged.ref_strain = ref_strain
                sampleMerged.strain = strain
                sampleMerged.condition = "{}_all".format(strain)
                sampleMerged.tgs_plat = dataObjs[0].tgs_plat
                sampleMerged.strategy = dataObjs[0].strategy
                sampleMerged.data_processed_location = [x.data_processed_location for x in dataObjs]
                polyaLocationDict = {}
                for x in dataObjs:
                    if x.polya_location != None and validateFile(x.polya_location):
                        polyaLocationDict.update({x.sample_name: x.polya_location})
                sampleMerged.polya_location = polyaLocationDict if len(polyaLocationDict) else None

                minLeftRepeats = [len(x.ngs_left_reads.split(";")) for x in dataObjs if x.ngs_left_reads != None]
                if minLeftRepeats:
                    minLeftRepeats = min(minLeftRepeats)
                    leftReadsRepeats = []
                    for x1 in range(minLeftRepeats):
                        leftReads = []
                        for x in dataObjs:
                            leftReads.append(x.ngs_left_reads.split(";")[x1])
                        leftReadsRepeats.append(",".join(leftReads))
                    sampleMerged.ngs_left_reads = ";".join(leftReadsRepeats)
                else:
                    sampleMerged.ngs_left_reads = None

                minRightRepeats = [len(x.ngs_right_reads.split(";")) for x in dataObjs if x.ngs_right_reads != None]
                if minRightRepeats:
                    minRightRepeats = min(minRightRepeats)
                    rightReadsRepeats = []
                    for x1 in range(minRightRepeats):
                        rightReads = []
                        for x in dataObjs:
                            rightReads.append(x.ngs_right_reads.split(";")[x1])
                        rightReadsRepeats.append(",".join(rightReads))
                    sampleMerged.ngs_right_reads = ";".join(rightReadsRepeats)
                else:
                    sampleMerged.ngs_right_reads = None

                if (not sampleMerged.ngs_left_reads and sampleMerged.ngs_right_reads) or (sampleMerged.ngs_left_reads and not sampleMerged.ngs_right_reads):
                    sampleMerged.ngs_reads_paired = "single"
                else:
                    sampleMerged.ngs_reads_paired = "paired"

                sampleMerged.use_fmlrc2 = dataObjs[0].use_fmlrc2

                if proj not in sampleMergedToProcess:
                    sampleMergedToProcess[proj] = {ref_strain: {strain: sampleMerged}}
                elif ref_strain not in sampleMergedToProcess[proj]:
                    sampleMergedToProcess[proj][ref_strain] = {strain: sampleMerged}
                else:
                    sampleMergedToProcess[proj][ref_strain].update({strain: sampleMerged})
    return sampleMergedToProcess


def checkAndMakeHisat2Index(refParams=None, threads=None, dirSpec=None):
    if not refParams.hisat2_index:
        refParams.hisat2_index = makeHisat2Index(refParams=refParams, dirSpec=dirSpec, threads=threads)
    else:
        if not checkHisat2IndexExist(refParams.hisat2_index):
            refParams.hisat2_index = makeHisat2Index(refParams=refParams, dirSpec=dirSpec, threads=threads)
        else:
            print getCurrentTime() + " It seems like you have provided hisat2 index, so the index will be used for mapping."


def sortTupleList(tupleList):
    return sorted(tupleList)


def mergeTupleList(tupleList):
    sortedTupleList = sortTupleList(tupleList)
    newTupleList = sortedTupleList[0:1]
    for myTuple in sortedTupleList[1:]:
        lastTuple = newTupleList.pop()
        if lastTuple[1] >= myTuple[0]:
            newTuple = (lastTuple[0], max(myTuple[1], lastTuple[1]))
            newTupleList.append(newTuple)
        else:
            newTupleList.extend([lastTuple, myTuple])
    return newTupleList


def isOverlap(listA, listB):
    sortedList = sorted([listA, listB])
    return True if sortedList[0][1] >= sortedList[1][0] else False


def getOverlapOfTuple(tupleListA, tupleListB):
    res = []
    i = j = 0
    while i < len(tupleListA) and j < len(tupleListB):
        lo = max(tupleListA[i][0], tupleListB[j][0])
        hi = min(tupleListA[i][1], tupleListB[j][1])
        if lo <= hi:
            res.append([lo, hi])
        if tupleListA[i][1] < tupleListB[j][1]:
            i += 1
        else:
            j += 1
    return res


def getIntrons(exonStarts, exonEnds):
    return exonEnds[:-1], exonStarts[1:]


def getSubSamByName(samFile, nameList, isBam=False, nameListIsFile=False, outPrefix=None, sort=True, threads=4):
    n = nameList
    if nameListIsFile:
        with open(nameList, 'r') as infile:
            n = infile.read().splitlines()
    if isBam:
        samHandle = pysam.AlignmentFile(samFile, 'rb', check_sq=False)
        suffix = ".bam"
    else:
        tmpBam = "tmp.bam"
        pysam.view("-o", tmpBam, "--output-fmt", "BAM", samFile, catch_stdout=False)
        samHandle = pysam.AlignmentFile(tmpBam, 'rb', check_sq=False)
        suffix = ".sam"
    name_indexed = pysam.IndexedReads(samHandle)
    name_indexed.build()
    header = samHandle.header.copy()
    if isBam:
        out = pysam.Samfile(outPrefix + suffix, 'wb', header=header)
    else:
        out = pysam.Samfile(outPrefix + suffix, 'w', header=header)
    for name in n:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                out.write(x)
    out.close()
    if not isBam:
        os.remove("tmp.bam")
    if sort:
        if isBam:
            pysam.sort("-o", outPrefix + ".sorted" + suffix, "-@", str(threads), outPrefix + suffix, catch_stdout=False)
        else:
            pysam.sort("-o", outPrefix + ".sorted" + suffix, "-@", str(threads), "--output-fmt", "SAM",
                       outPrefix + suffix, catch_stdout=False)
        os.remove(outPrefix + suffix)
        return outPrefix + ".sorted" + suffix
    else:
        return outPrefix + suffix


def getFxSequenceId(fxFile, isFa=False, isFq=False, outFile=None):
    recordId = []
    if isFa:
        for rec in SeqIO.parse(fxFile, "fasta"):
            recordId.append(rec.id)
        if outFile:
            with open(outFile, 'w') as f:
                for item in recordId:
                    print >> f, item
        else:
            return recordId
    if isFq:
        for rec in SeqIO.parse(fxFile, "fastq"):
            recordId.append(rec.id)
        if outFile:
            with open(outFile, 'w') as f:
                for item in recordId:
                    print >> f, item
        else:
            return recordId


def getBlockLength(blockList):
    return sum(map(lambda x: int(x[1]) - int(x[0]), blockList))


def getConsensusIntronN(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
    '''
    get consensus intron number, i.e. the intron info of reads is identical to the reference
    :return: the count of the consensus introns
    '''
    intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
    intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
    j, consensusN = 0, 0
    for i in range(len(intronStarts1)):
        for k in range(j, len(intronStarts2)):
            if intronStarts2[k] > intronStarts1[i]: break
            if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronEnds1[i] + offset and \
                    intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
                consensusN += 1
                j += 1
    return consensusN


def getSpliceMotifFromJuncDict(juncDict, genomeFasta, withStrand=False):
    dinucleotideBedList = []
    for juncName in juncDict:
        juncBed = juncDict[juncName]
        chrom = juncBed.split("\t")[0]
        juncStart = int(juncBed.split("\t")[1])
        juncEnd = int(juncBed.split("\t")[2])
        strand = juncBed.split("\t")[-1]
        leftDinucleotidePos = "\t".join(map(str, [chrom, juncStart, juncStart + 2, ":".join([juncName, "left"]), ".", strand]))
        rightDinucleotidePos = "\t".join(map(str, [chrom, juncEnd - 2, juncEnd, ":".join([juncName, "right"]), ".", strand]))
        dinucleotideBedList.extend([leftDinucleotidePos, rightDinucleotidePos])
    dinucleotideBedObj = pybedtools.BedTool("\n".join(dinucleotideBedList), from_string=True)
    dinucleotideBedRes = dinucleotideBedObj.sequence(genomeFasta, name=True, tab=True, s=True)

    juncMotifDict = {}
    for i in str(open(dinucleotideBedRes.seqfn).read()).split("\n")[:-1]:
        infoList = str(i).strip("\n").split("\t")
        strand = infoList[0][infoList[0].find("(") + 1:infoList[0].find(")")]
        if not withStrand:
            juncName = "{}:{}".format(":".join(infoList[0].split(":")[:2]), strand)
            dinucleotideType = infoList[0].split(":")[2]
        else:
            juncName = ":".join(infoList[0].split(":")[:3])
            dinucleotideType = infoList[0].split(":")[3]

        if juncName not in juncMotifDict:
            juncMotifDict.update({juncName: {dinucleotideType: infoList[1]}})
        else:
            juncMotifDict[juncName].update({dinucleotideType: infoList[1]})

    junc2motif = {}
    for juncName in juncMotifDict:
        strand = juncName.split(":")[-1]
        if strand == "+":
            spliceMotif = "{}-{}".format(juncMotifDict[juncName]["left"], juncMotifDict[juncName]["right"])
        else:
            spliceMotif = "{}-{}".format(juncMotifDict[juncName]["right"], juncMotifDict[juncName]["left"])
        if juncName not in junc2motif:
            junc2motif[juncName] = spliceMotif

    return junc2motif


def getCurrentTime():
    return str(datetime.datetime.now())


def getSizes(starts, ends):
    return map(lambda x: int(ends[x]) - int(starts[x]), range(len(starts)))


def getRelStarts(blockStarts):
    return map(lambda x: int(blockStarts[x]) - int(blockStarts[0]), range(len(blockStarts)))


def getDictFromFile(myFile, sep="\t", inlineSep=None, keyCol=1, valueCol=2):
    with open(myFile) as f:
        myDict = {}
        for line in f.readlines():
            infoList = line.strip("\n").split(sep)
            key = infoList[keyCol - 1]
            if valueCol:
                if inlineSep:
                    value = infoList[valueCol - 1].split(inlineSep)
                else:
                    value = infoList[valueCol - 1]
            else:
                value = infoList
            myDict[key] = value
        return myDict


def getColumnMapping(myFile, sep="\t", keyCol=1, valueCol=2, keySep=None, valueSep=None):
    with open(myFile) as f:
        mappingDict = {}
        for line in f.readlines():
            infoList = line.strip("\n").split(sep)
            if valueSep and valueSep:
                return

            if keySep:
                keys = infoList[keyCol - 1].split(keySep)
                for k in keys:
                    mappingDict[k] = infoList[valueCol - 1]

            if valueSep:
                values = infoList[valueCol - 1].split(valueSep)
                for v in values:
                    mappingDict[infoList[keyCol - 1]] = v

            if not keySep and not valueSep:
                mappingDict[infoList[keyCol - 1]] = infoList[valueCol - 1]
        return mappingDict


def initRefSetting(refParams=None, dirSpec=None):
    originDir = os.getcwd()
    if refParams.ref_genome == None:
        raise Exception("You must specify a genome file!")
    refParams.ref_genome = os.path.abspath(os.path.join(originDir, refParams.ref_genome))
    validateFile(refParams.ref_genome)

    if refParams.ref_gtf == None and refParams.ref_gff == None:
        raise Exception("You must specify a GTF or GFF3 file!")
    elif refParams.ref_gtf == None and refParams.ref_gff != None:
        refParams.ref_gff = os.path.abspath(os.path.join(originDir, refParams.ref_gff))
        validateFile(refParams.ref_gff)
        gffDir = os.path.dirname(refParams.ref_gff)
        ref_gtf = os.path.abspath(os.path.join(gffDir, "iFLAS.gtf"))
        cmd = "gffread -T -o {} {}".format(ref_gtf, refParams.ref_gff)
        subprocess.call(cmd, shell=True)
        refParams.ref_gtf = ref_gtf

    gtfDir = os.path.dirname(refParams.ref_gtf)
    if refParams.ref_size == None:
        refParams.ref_size = os.path.join(gtfDir, "iFLAS.len.txt")
        sizeOut = open(refParams.ref_size)
        for seqRecord in SeqIO.parse(refParams.ref_genome, "fasta"):
            print >> sizeOut, "\t".join([seqRecord.id, seqRecord.seq])
        sizeOut.close()
    else:
        refParams.ref_size = os.path.join(originDir, refParams.ref_size)

    if refParams.ref_gpe == None:
        refParams.ref_gpe = os.path.join(gtfDir, "iFLAS.gpe")
        cmd = "gtfToGenePred {} {} -genePredExt".format(refParams.ref_gtf, refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref_gpe = os.path.abspath(os.path.join(originDir, refParams.ref_gpe))

    if refParams.ref_bed == None:
        refParams.ref_bed = os.path.join(gtfDir, "iFLAS.annotation.bed")
        scriptDir = os.path.dirname(os.path.abspath(__file__))
        utilDir = os.path.join(scriptDir, "utils")
        cmd = "{}/gpe2bed.pl -g {} > {}".format(utilDir, refParams.ref_gpe, refParams.ref_bed)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref_bed = os.path.abspath(os.path.join(originDir, refParams.ref_bed))
    refParams.ref_mm2_index = os.path.abspath(os.path.join(originDir, refParams.ref_mm2_index))

    if dirSpec.tmp_dir == None:
        dirSpec.tmp_dir = "/tmp"
    else:
        if not os.path.exists(dirSpec.tmp_dir):
            batchCreateDir([dirSpec.tmp_dir])
    dirSpec.out_dir = os.path.abspath(os.path.join(originDir, dirSpec.out_dir))
    if not os.path.exists(dirSpec.out_dir):
        batchCreateDir([dirSpec.out_dir])


def initSysResourceSetting(optionTools=None):
    maxMem = psutil.virtual_memory().total >> 20
    if str(optionTools.memory).endswith("M"):
        if int(optionTools.memory[0:-1]) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    elif str(optionTools.memory).endswith("G"):
        optionTools.memory = str(int(optionTools.memory[0:-1]) * 1024) + "M"
        if int(optionTools.memory[0:-1]) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    elif optionTools.memory.isdigit():
        if int(optionTools.memory) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    else:
        raise Exception("Please check the format of memory you specified in optionTool section!")

    maxThreads = psutil.cpu_count()
    if str(optionTools.threads).isdigit():
        if int(optionTools.threads) > maxThreads * 0.75:
            optionTools.threads = int(maxThreads * 0.75)
        else:
            optionTools.threads = int(optionTools.threads)
    else:
        raise Exception("The number of threads should be digital, please check it!")


# def isFastaOrFastq(filename):
#     with open(filename, "r") as handle:
#         fasta = SeqIO.parse(handle, "fasta")
#         if any(fasta):
#             return "fasta"
#     with open(filename, "r") as handle:
#         fastq = SeqIO.parse(handle, "fastq")
#         if any(fastq):
#             return "fastq"
#     return None

def listStr2Int(myList):
    return [int(i) for i in myList]


def listInt2Str(myList):
    return [str(i) for i in myList]


def nestedDictValues(iterable, returned="key"):
    if isinstance(iterable, dict):
        for key, value in iterable.items():
            if returned == "key":
                yield key
            elif returned == "value":
                if not (isinstance(value, dict) or isinstance(value, list)):
                    yield value
            else:
                raise ValueError("'returned' keyword only accepts 'key' or 'value'.")
            for ret in nestedDictValues(value, returned=returned):
                yield ret
    elif isinstance(iterable, list):
        for el in iterable:
            for ret in nestedDictValues(el, returned=returned):
                yield ret


def resolveDir(dirName, chdir=True):
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    if chdir:
        os.chdir(dirName)


def removeFiles(myDir=None, fileList=None):
    for f in fileList:
        if myDir:
            os.remove(os.path.join(myDir, f.strip("\n")))
        else:
            os.remove(f.strip("\n"))


def renameNGSdata2fastp(dataObj=None):
    if "fastp" in dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
        return
    if dataObj.ngs_left_reads == None and "fastp" in dataObj.ngs_right_reads:
        return
    if "fastp" in dataObj.ngs_left_reads and "fastp" in dataObj.ngs_right_reads:
        return
    if dataObj.ngs_reads_paired == "paired":
        leftReadsRepeats = [i.strip() for i in dataObj.ngs_left_reads.split(";")]
        rightReadsRepeats = [i.strip() for i in dataObj.ngs_right_reads.split(";")]
        if len(leftReadsRepeats) != len(rightReadsRepeats):
            raise Exception("The repeats of your NGS data not match between your left reads and right reads")
        else:
            newLeftReadsRepeats = []
            newRightReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                rightReads = rightReadsRepeats[i].split(",")
                if len(leftReads) != len(rightReads):
                    raise Exception("please input the right paired NGS reads for the analysis")
                else:
                    newLeftReads = []
                    newRightReads = []
                    for j in range(len(leftReads)):
                        leftReadsDir = os.path.dirname(leftReads[j])
                        rightReadsDir = os.path.dirname(rightReads[j])
                        leftReadsBase = os.path.basename(leftReads[j])
                        rightReadsBase = os.path.basename(rightReads[j])
                        newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                        newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                        newLeftReads.append(newLeft)
                        newRightReads.append(newRight)
                    newLeftReadsRepeats.append(",".join(newLeftReads))
                    newRightReadsRepeats.append(",".join(newRightReads))
            dataObj.ngs_left_reads = ";".join(newLeftReadsRepeats)
            dataObj.ngs_right_reads = ";".join(newRightReadsRepeats)
    else:
        if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
            leftReadsRepeats = [i.strip() for i in dataObj.ngs_left_reads.split(";")]
            newLeftReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                newLeftReads = []
                for j in range(len(leftReads)):
                    leftReadsDir = os.path.dirname(leftReads[j])
                    leftReadsBase = os.path.basename(leftReads[j])
                    newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                    newLeftReads.append(newLeft)
                newLeftReadsRepeats.append(",".join(newLeftReads))
            dataObj.ngs_left_reads = ";".join(newLeftReadsRepeats)
        elif dataObj.ngs_right_reads and dataObj.ngs_left_reads == None:
            rightReadsRepeats = [i.strip() for i in dataObj.ngs_right_reads.split(";")]
            newRightReadsRepeats = []
            for i in range(len(rightReadsRepeats)):
                rightReads = rightReadsRepeats[i].split(",")
                newRightReads = []
                for j in range(len(rightReads)):
                    rightReadsDir = os.path.dirname(rightReads[j])
                    rightReadsBase = os.path.basename(rightReads[j])
                    newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                    newRightReads.append(newRight)
                newRightReadsRepeats.append(",".join(newRightReads))
            dataObj.ngs_right_reads = ";".join(newRightReadsRepeats)
        else:
            raise Exception("The NGS data seem not to be single, please check it")


def validateFaAndFqFile(myFile):
    cmd = "seqkit head -n 100 {} | seqkit stat".format(myFile)
    validateOutput = os.popen(cmd).read()
    if "FASTQ" in validateOutput:
        return "fastq"
    elif "FASTA" in validateOutput:
        return "fasta"
    else:
        return False


def validateFile(myFile):
    if isinstance(myFile, list):
        return False

    if not os.path.exists(myFile):
        return False

    if not os.path.isfile(myFile):
        return False

    return True


def validateDir(myDir):
    if not os.path.exists(myDir):
        return False

    if not os.path.isdir(myDir):
        return False

    return True