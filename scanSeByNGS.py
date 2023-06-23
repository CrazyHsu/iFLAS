#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: scanSeByNGS.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

from itertools import combinations
from commonFuncs import *
from commonObjs import *

class Position(object):
    "chr:st-ed, zero based"

    def __init__(self, line):
        self.chrom, st_ed = line.strip().split(":")
        self.chromStart, self.chromEnd = (int(c) for c in st_ed.split("-"))
        assert self.chromStart < self.chromEnd

    def __repr__(self):
        return "%s:%d-%d" % (self.chrom, self.chromStart, self.chromEnd)


class SE_JUNC(Bed12):
    "SE events identified by RNA-seq and supported by reference or pacbio"

    def __init__(self, bed12_se, junc_left, junc_right):
        Bed12.__init__(self, repr(bed12_se))
        self.junc_left = junc_left
        self.junc_right = junc_right
        self.junction = Position("%s:%d-%d" % (self.chrom, junc_left, junc_right))
        donor = junc_left if self.strand == "+" else junc_right
        self.jpos = repr(self.junction)
        self.dnpos = "%s:%d:%s" % (self.chrom, donor, self.strand)
        self.posCode = "%d@%s@%d" % (
        self.junc_left, ";".join(["%d-%d" % (st, ed) for (st, ed) in zip(self.exonStarts, self.exonEnds)]),
        self.junc_right)
        self.supCount = 0
        self.totalCount = 0

    def __repr__(self):
        return "\t".join([Bed12.__repr__(self), repr(self.supCount),
                          repr(self.totalCount), self.jpos])


class SE(Bed12):
    def __init__(self, line):
        Bed12.__init__(self, line)
        self.geneName, self.se_info = self.name.split(":")
        self.supCount = int(self.record[12])
        self.totalCount = int(self.record[13])
        self.juncPoses = [Position(j) for j in self.record[14].split(",")]

    def __repr__(self):
        return "\t".join([Bed12.__repr__(self), repr(self.supCount),
                          repr(self.totalCount), ",".join([repr(j) for j in self.juncPoses])])


def check_ses(chrom, dn, acs, refDct, lrDct, jpos2count, dnpos2count, novel_count):
    """return a list of SE(Bed12+) for enumerated dn-ac(ac is element in acs)
    that can be supported by ref in refLst or lr in lrList.
    """
    # sort acs by position
    acs.sort()

    # check SE for dn with each 2 acs
    strand = "+" if dn < acs[0] else "-"
    seLst = []
    for ac1, ac2 in combinations(acs, 2):
        if strand == "+":
            ac_i, ac_o = ac1, ac2
        else:
            ac_i, ac_o = ac2, ac1
        seLst.extend(identify_ses(chrom, strand, dn, ac_i, ac_o, refDct, lrDct, jpos2count, dnpos2count, novel_count))
    return seLst


def find_refs(chrom, strand, dn, ac_i, ac_o, refDct, dnpos2count):
    "return a list of ref matching input situation"
    if chrom not in refDct:
        return []

    refLst = refDct[chrom]
    refs = []
    left = min(dn, ac_o)
    right = max(dn, ac_o)
    for ref in refLst:
        if ref.strand != strand:
            continue
        elif ref.exonStarts[0] >= right or ref.exonEnds[-1] <= left:
            continue
        else:
            if strand == "+":
                if (dn in ref.exonEnds and ac_i in ref.exonStarts and ac_o in ref.exonStarts):
                    ac_o_index = ref.exonStarts.index(ac_o)
                    dnPos = "%s:%d:%s" % (chrom, ref.exonEnds[ac_o_index - 1], strand)
                    if dnPos in dnpos2count:
                        refs.append(ref)
            else:
                if (dn in ref.exonStarts and ac_i in ref.exonEnds and ac_o in ref.exonEnds):
                    # refs.append(ref)
                    ac_o_index = ref.exonEnds.index(ac_o)
                    dnPos = "%s:%d:%s" % (chrom, ref.exonStarts[ac_o_index + 1], strand)
                    if dnPos in dnpos2count:
                        refs.append(ref)
    return refs


def find_se_refs(se, refLst):
    refs = []
    for ref in refLst:
        if ref.strand != se.strand:
            continue
        else:
            for junc in se.juncPoses:
                if (junc.chromStart in ref.exonEnds or junc.chromEnd in ref.exonStarts):
                    refs.append(ref)
                    break
    return refs


def se_structure(strand, dn, ac_i, ac_o, ref):
    "return SE events (in SE object) in ref"
    left = min([ac_i, ac_o])
    right = max([ac_i, ac_o])

    bed = Bed12()
    bed.chrom = ref.chrom
    bed.strand = strand
    exSts = []
    exEds = []
    for st, ed in zip(ref.exonStarts, ref.exonEnds):
        if st >= left and ed <= right:
            exSts.append(st)
            exEds.append(ed)
        elif st >= right:
            break
    bed.exonStarts = exSts
    bed.exonEnds = exEds
    bed.chromStart = bed.thickStart = exSts[0]
    bed.chromEnd = bed.thickEnd = exEds[-1]
    bed.blockStarts = getRelStarts(exSts)
    bed.blockSizes = getSizes(exSts, exEds)
    bed.blockCount = len(exSts)

    if strand == "+":
        junc_left = dn
        junc_right = ac_o
    else:
        junc_left = ac_o
        junc_right = dn
    se = SE_JUNC(bed, junc_left, junc_right)
    return se


def change_name(se, refDct, errorOut):
    "change name if necessary"
    if se.chrom not in refDct:
        return 1
    refs = find_se_refs(se, refDct[se.chrom])
    if refs:
        geneNames = set([ref.geneName for ref in refs])
        if len(geneNames) > 1:
            return 2
        else:
            gName = geneNames.pop()
            print >> errorOut, "%s\t%s" % (se.geneName, gName)
            se.geneName = gName
            se.name = ":".join([se.geneName, se.se_info])
            return 0
    else:
        return 1


def remove_redun(seLst):
    """return a sublist of seLst with redundent ones (same block and junction)
     removed"""
    newSes = []
    seCodeSet = set()    # se.posCode: dn@st-ed@ac
    for se in seLst:
        if se.posCode in seCodeSet:
            continue
        else:
            newSes.append(se)
            seCodeSet.add(se.posCode)
    return newSes


def union_ses(refSes, lrSes, novel_count):
    """return a list of se which is the union of the 2 ses without redundency,
    change name of the se if necessary"""
    uSes = []
    seCodeSet = set()
    for se in refSes:
        se.name = "%s:%s" % (se.geneName, se.posCode)
        seCodeSet.add(se.posCode)
        uSes.append(se)

    for se in lrSes:
        if se.posCode in seCodeSet:
            continue
        else:
            # global novel_count
            novel_count.count += 1
            se.name = "NovelSE_%d:%s" % (novel_count.count, se.posCode)
            seCodeSet.add(se.posCode)
            uSes.append(se)
    return uSes


def cal_psi(se, jpos2count, dnpos2count):
    "calculate supCount, totalCount and score for se"
    se.supCount = jpos2count[se.jpos]
    se.totalCount = sum([tup[1] for tup in dnpos2count[se.dnpos]])
    se.score = int(float(se.supCount)/se.totalCount * 1000)


def identify_ses(chrom, strand, dn, ac_i, ac_o, refDct, lrDct, jpos2count, dnpos2count, novel_count):
    """return a list of bed12+ of SEs, additional fields: jpos, dnpos.
    """
    # find refs or pacbio reads covering dn to ac_o
    has_ref = False
    has_lr = False

    refs = find_refs(chrom, strand, dn, ac_i, ac_o, refDct, dnpos2count)
    if refs:
        has_ref = True
    lrs = find_refs(chrom, strand, dn, ac_i, ac_o, lrDct, dnpos2count)
    if lrs:      # use pacbio reads to varify SE structure
        has_lr = True

    if not has_ref and not has_lr:
        return []

    # check se for each ref or lr
    refSes = []
    if has_ref:
        for ref in refs:
            se = se_structure(strand, dn, ac_i, ac_o, ref)
            se.geneName = ref.geneName
            refSes.append(se)
    lrSes = []
    if has_lr:
        for lr in lrs:
            se = se_structure(strand, dn, ac_i, ac_o, lr)
            lrSes.append(se)

    # remove redun in refSe and lrSe, then remove redun in between
    refSes = remove_redun(refSes)
    lrSes = remove_redun(lrSes)
    seLst = union_ses(refSes, lrSes, novel_count)

    newSeLst = []
    for se in seLst:
        cal_psi(se, jpos2count, dnpos2count)
        newSeLst.append(se)
        # print >>out, se
    return newSeLst


def scanSeByNGS(fref, has_bin, fjunc, flr, outFile):
    # load gene structures
    refDct = {}
    with open(fref) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ref = GenePredExtLine(line, bincolumn=has_bin)
            if ref.chrom in refDct:
                refDct[ref.chrom].append(ref)
            else:
                refDct[ref.chrom] = [ref]

    for chrom in refDct:
        refDct[chrom].sort(key=lambda r: r.txStart)

    # load long reads
    lrDct = {}
    with open(flr) as f:
        for line in f:
            if line.startswith("#"):
                continue
            read = Bed12(line)
            if read.chrom in lrDct:
                lrDct[read.chrom].append(read)
            else:
                lrDct[read.chrom] = [read]

    for chrom in lrDct:
        lrDct[chrom].sort(key=lambda r: r.chromStart)

    # load junctions into {chrom:{donor:[acceptors], ...}, ...}
    dn2ac = {}
    jpos2count = {}
    dnpos2count = {}      # {dnpos:[(ac,count), (ac,count), ...], ...}
    with open(fjunc) as f:
        for line in f:
            if line.startswith("#"):
                continue
            junc = Junction(line)

            jpos2count[junc.jPos] = junc.score
            if junc.dnPos in dnpos2count:
                dnpos2count[junc.dnPos].append((junc.ac, junc.score))
            else:
                dnpos2count[junc.dnPos] = [(junc.ac, junc.score)]

            if junc.strand == "+":
                dn = junc.donor
                ac = junc.acceptor
            else:
                dn = junc.acceptor
                ac = junc.donor

            if junc.chrom in dn2ac:
                if dn in dn2ac[junc.chrom]:
                    dn2ac[junc.chrom][dn].append(ac)
                else:
                    dn2ac[junc.chrom][dn] = [ac]
            else:
                dn2ac[junc.chrom] = {}
                dn2ac[junc.chrom][dn] = [ac]

    # check dn with multiple acs
    seLst_p = []
    seLst_n = []
    p_novel_count = Number()
    n_novel_count = Number()
    for chrom in dn2ac:
        for dn in dn2ac[chrom]:
            acs = list(set(dn2ac[chrom][dn]))
            if len(acs) > 1:
                pos_acs = [ac for ac in acs if ac > dn]
                neg_acs = [ac for ac in acs if ac < dn]
                if len(pos_acs) > 1:
                    seLst_p.extend(check_ses(chrom, dn, pos_acs, refDct, lrDct, jpos2count, dnpos2count, p_novel_count))
                if len(neg_acs) > 1:
                    seLst_n.extend(check_ses(chrom, dn, neg_acs, refDct, lrDct, jpos2count, dnpos2count, n_novel_count))
    out = open(outFile, "w")
    for i in seLst_p:
        print >>out, i
    for i in seLst_n:
        print >>out, i
    out.close()

def assignGeneNameToSE(fref, has_bin, fes, outFile, errorOutFile):
    # load gene structures
    refDct = {}
    with open(fref) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ref = GenePredExtLine(line, bincolumn=has_bin)
            if ref.chrom in refDct:
                refDct[ref.chrom].append(ref)
            else:
                refDct[ref.chrom] = [ref]
    for chrom in refDct:
        refDct[chrom].sort(key=lambda r: r.txStart)

    out = open(outFile, "w")
    errorOut = open(errorOutFile, "w")
    with open(fes) as f:
        for line in f:
            if line.startswith("#"):
                print >>out, line,
                continue
            se = SE(line)
            if se.geneName.startswith("Novel"):
                change_name(se, refDct, errorOut)
                print >>out, se
            else:
                print >>out, line,
    out.close()
    errorOut.close()

