"""Filesystem and simple file-parsing helpers — faithfully ported from legacy commonFuncs.py."""

from __future__ import annotations

import os

from ..runtime import ToolError, ToolRunner



def makeLink(sourcePath, targetPath):
    if os.path.islink(targetPath) or os.path.exists(targetPath):
        os.remove(targetPath)
    os.symlink(sourcePath, targetPath)


def filterFile(originFile=None, targetFile=None, originField=1, targetField=1, mode="i", outFile=None,
               returnFlag=False):
    """Port of legacy commonFuncs.filterFile: the target entries (list items, dict keys, or file lines keyed by
    column ``targetField`` -- the last line per key) whose key is (mode "i") or is not (any other mode) one of the
    origin keys. They are written in target order: legacy wrote set order, which Python 3 changes from run to run
    (BUG-162)."""
    origin = set()
    if isinstance(originFile, list):
        origin = set(originFile)
    elif os.path.isfile(originFile):
        with open(originFile) as f:
            origin = {line.strip("\n").split("\t")[originField - 1] for line in f}
    entries = {}
    if isinstance(targetFile, (list, dict)):
        entries = {k: k for k in targetFile}
    elif os.path.isfile(targetFile):
        with open(targetFile) as f:
            for line in f:
                entries[line.strip("\n").split("\t")[targetField - 1]] = line.strip("\n")
    returnList = [k for k in entries if (k in origin) == (mode == "i")]
    with open(outFile or "filtered.txt", "w") as out:
        out.write("".join(entries[k] + "\n" for k in returnList))
    if returnFlag:
        return returnList


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


def validateFaAndFqFile(myFile, runner=None):
    "Detect fasta/fastq via seqkit (returns 'fastq'/'fasta'/False). Requires seqkit on PATH."
    runner = runner or ToolRunner()
    try:
        head_proc = runner.run(["seqkit", "head", "-n", "100", myFile], check=False)
    except ToolError:
        return False
    if head_proc.returncode != 0:
        return False

    try:
        stat_proc = runner.run(["seqkit", "stat"], input=head_proc.stdout or "", check=False)
    except ToolError:
        return False
    if stat_proc.returncode != 0:
        return False

    validateOutput = (stat_proc.stdout or "") + (stat_proc.stderr or "")
    if "FASTQ" in validateOutput:
        return "fastq"
    elif "FASTA" in validateOutput:
        return "fasta"
    else:
        return False


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


