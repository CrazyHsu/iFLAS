"""find_as stage — alternative-splicing event detection (ported from legacy find_as.py).

The four detectors (findIR/findSE/findA5SS/findA3SS), the annotation splitter (getAnnotation), the consensus-intron
helpers and find_single_as keep the published algorithm exactly: same events, PSI and output order
(tests/test_find_as_characterization.py).

The drivers read the post-PU bed12+ (``hqIsoforms/hq.collapsed.bed12+``) and, when NGS junctions exist, apply the
legacy NGS confident filters (perl utils).
"""

from __future__ import annotations

import itertools
import os
import re

from ..constants import AStypes
from ..intervalops import (getIntrons, getOverlapOfTuple, getRelStarts,
                           getSizes, isOverlap, listStr2Int)
from ..io.files import makeLink
from ..io.genepred import GenePredObj
from ..io.intervals import Bed12, Bed12Plus, Junction
from ..io.records import Gene2Reads, ReadLineStruc
from ..runtime import ToolRunner


def getCnsSite(intronStarts, transIntronStarts, offset):
    '''
    get consensus intron boundary site between reference and reads info
    :return: count of consensus site num
    '''
    s2, cnsN = 0, 0
    for intronStart1 in intronStarts:
        for i in range(s2, len(transIntronStarts)):
            intronStart2 = transIntronStarts[i]
            if intronStart2 > intronStart1 + offset: break
            if intronStart1 - offset <= intronStart2 and intronStart2 <= intronStart1 + offset:
                cnsN += 1
                s2 += 1
    return cnsN


def _cns_introns(intronStarts1, intronEnds1, intronStarts2, intronEnds2, offset):
    """Introns of structure 1 matched within ``offset`` at both ends by an intron of structure 2."""
    j, consensusN = 0, 0
    for i in range(len(intronStarts1)):
        for k in range(j, len(intronStarts2)):
            if intronStarts2[k] > intronEnds1[i]:
                break
            if abs(intronStarts2[k] - intronStarts1[i]) <= offset and abs(intronEnds2[k] - intronEnds1[i]) <= offset:
                consensusN += 1
                j += 1
    return consensusN


def getConsensusIntronNfromAS(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
    '''
    Intron evidence shared by a read and a reference transcript, counted once per AS type as published: whole introns
    for IR and SE, intron starts plus intron ends for A3SS and A5SS (callers only test ``>= 1``).
    '''
    intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
    intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
    introns = _cns_introns(intronStarts1, intronEnds1, intronStarts2, intronEnds2, offset)
    sites = getCnsSite(intronStarts1, intronStarts2, offset) + getCnsSite(intronEnds1, intronEnds2, offset)
    return sum(introns if asType in ("IR", "SE") else sites for asType in AStypes)


def filterIrByJunc(irBed, junctionBed, confidentIrBed):
    juncDict = {}
    with open(junctionBed) as f:
        for line in f.readlines():
            juncDict.update({Junction(line.strip("\n")).jPos: ""})

    out = open(confidentIrBed, "w")
    with open(irBed) as f:
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            juncPos = "{}:{}-{}".format(infoList[0], infoList[1], infoList[2])
            if juncPos in juncDict:
                print(line.strip("\n"), file=out)
    out.close()


def _span_inside(read, trans):
    """Fraction of the read's span that lies inside the transcript's span."""
    span, tspan = (read.chromStart, read.chromEnd), (trans.chromStart, trans.chromEnd)
    if not isOverlap(span, tspan):
        return 0
    return sum(e - s for s, e in getOverlapOfTuple([span], [tspan])) / float(read.chromEnd - read.chromStart)


def getAnnotation(annoBedRes, novelBedRes, offset=0):
    """Split reads into annotated genes and novel reads. A read joins the gene of an overlapping reference transcript
    when it shares an intron boundary with it, or, as a single-exon read, when at least half of its span lies inside
    the transcript's span; all other reads are novel, grouped by chromosome and strand."""
    annoDict, novelDict, firstSeen, annotated = {}, {}, {}, set()
    for line in annoBedRes:
        record = str(line).strip("\n").split("\t")
        read = Bed12("\t".join(record[:12]))
        trans = ReadLineStruc("\t".join(record[12:]))
        firstSeen.setdefault(read.name, read)
        if getConsensusIntronNfromAS(read.exonStarts, read.exonEnds, trans.exonStarts, trans.exonEnds, offset) >= 1 \
                or (read.blockCount == 1 and _span_inside(read, trans) >= 0.5):
            gene = annoDict.get(trans.geneName)
            if gene is None:
                gene = annoDict[trans.geneName] = Gene2Reads(trans.geneName)
            gene.update(read)
            gene.trans[trans.name] = trans
            annotated.add(read.name)
    novel = [Bed12(str(line)) for line in novelBedRes]
    novel += [read for name, read in firstSeen.items() if name not in annotated]
    for read in novel:
        novelDict.setdefault(read.chrom, {}).setdefault(read.strand, []).append(read)
    return annoDict, novelDict


def _psi(inc, exc, isoform2reads):
    """Inclusion level x1000 over the isoforms, or over their full-length reads when ``isoform2reads`` is given."""
    nInc, nExc = len(inc), len(exc)
    if isoform2reads:
        nInc, nExc = sum(len(isoform2reads[x]) for x in inc), sum(len(isoform2reads[x]) for x in exc)
    return float(nInc) / (nInc + nExc) * 1000


def _overlaps_transcript(start, end, trans):
    return any(start < t.chromEnd and end > t.chromStart for t in trans.values())


def _intron_reads(read1, i, others, offset):
    """Reads splicing intron ``i`` of read1 (the same two flanking boundaries) and reads retaining it (one exon spanning
    it, both ends more than ``offset`` away)."""
    intronStart, intronEnd, nextEnd = read1.exonEnds[i], read1.exonStarts[i + 1], read1.exonEnds[i + 1]
    ev = {"spliced": [read1]}
    for read2 in others:
        exons2 = read2.exons
        for k, (s2, e2) in enumerate(exons2):
            if s2 < intronStart and intronEnd < e2 and abs(s2 - intronStart) > offset and abs(intronEnd - e2) > offset:
                if "retention" not in ev:             # as published: only the first retaining read ends its scan here
                    ev["retention"] = [read2]
                    break
                ev["retention"].append(read2)
            if e2 == intronStart:
                if k + 1 >= len(exons2):
                    break
                if exons2[k + 1][0] == intronEnd:
                    ev["spliced"].append(read2)
            if s2 > nextEnd:
                break
    return ev


def findIR(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None):
    for gName, gene in gene2ReadsDict.items():
        reads = list(gene.reads.items())
        if len(reads) < 2:
            continue
        irDict = {}
        for name1, read1 in reads:
            others = [read2 for name2, read2 in reads if name2 != name1]
            for i in range(len(read1.exonStarts) - 1):
                ir = "%d-%d" % (read1.exonEnds[i], read1.exonStarts[i + 1])
                if ir not in irDict:
                    irDict[ir] = _intron_reads(read1, i, others, offset)
        events = {ir: ev for ir, ev in irDict.items() if "retention" in ev}
        if outAS:
            for ir, ev in events.items():
                irStart, irEnd = listStr2Int(ir.split("-"))
                if anno and not _overlaps_transcript(irStart, irEnd, gene.trans):
                    continue
                retention, spliced = [x.name for x in ev["retention"]], [x.name for x in ev["spliced"]]
                print("\t".join(map(str, [gene.chrom, irStart, irEnd, gName + ":" + ir, _psi(retention, spliced, isoform2reads),
                                          gene.strand, gName, ",".join(retention), len(retention), ",".join(spliced),
                                          len(spliced)])), file=out)
        gene.asDict["IR"] = events


def _inner_exon_runs(read):
    """Each run of inner exons of a read: (end of the exon before, the run as 'start-end;...', start of the exon after)."""
    exons = read.exonChain.split(";")
    for l in range(len(exons) - 2):
        for r in range(l + 2, len(exons)):
            yield int(exons[l].split("-")[1]), ";".join(exons[l + 1:r]), int(exons[r].split("-")[0])


def _skip_reads(read1, others, skipPat, keepPat):
    """Reads whose exon chain joins the two flanking boundaries directly (skip) or through the same inner exons (keep)."""
    ev = {"keep": [read1]}
    for read2 in others:
        if skipPat in read2.exonChain:                    # literal text: re.search == substring test
            ev.setdefault("skip", []).append(read2)
        elif keepPat in read2.exonChain:
            ev["keep"].append(read2)
    return ev


def _write_se_event(out, gName, gene, seChain, ev, anno, isoform2reads):
    lSpliceSite, SEs, rSpliceSite = seChain.split("@")
    exons = [listStr2Int(e.split("-")) for e in SEs.split(";")]
    starts, ends = [s for s, _ in exons], [e for _, e in exons]
    if anno and not _overlaps_transcript(starts[0], ends[-1], gene.trans):
        return
    keep, skip = [x.name for x in ev["keep"]], [x.name for x in ev["skip"]]
    print("\t".join(map(str, [gene.chrom, starts[0], ends[-1], "%s:%s" % (gName, seChain), _psi(keep, skip, isoform2reads),
                              gene.strand, starts[0], ends[-1], "0,0,0", len(starts), ",".join(map(str, getSizes(starts, ends))),
                              ",".join(map(str, getRelStarts(starts))), gName, lSpliceSite, rSpliceSite,
                              ",".join(keep), len(keep), ",".join(skip), len(skip)])), file=out)


def findSE(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None):
    for gName, gene in gene2ReadsDict.items():
        reads = list(gene.reads.items())
        if len(reads) < 2:
            continue
        seDict = {}
        for name1, read1 in reads:
            others = [read2 for name2, read2 in reads if name2 != name1]
            for lExonEnd, SEs, rExonStart in _inner_exon_runs(read1):
                seChain = "%d@%s@%d" % (lExonEnd, SEs, rExonStart)
                if seChain not in seDict:
                    seDict[seChain] = _skip_reads(read1, others, "-%d;%d-" % (lExonEnd, rExonStart),
                                                  "-%d;%s;%d-" % (lExonEnd, SEs, rExonStart))
        events = {se: ev for se, ev in seDict.items() if "skip" in ev}
        if outAS:
            for seChain, ev in events.items():
                _write_se_event(out, gName, gene, seChain, ev, anno, isoform2reads)
        gene.asDict["SE"] = events


def _record_pair(altDict, key, inc, exc, incJunc, excJunc):
    """One inclusion/exclusion read pair of an alternative splice site, keyed by its two boundaries."""
    entry = altDict.setdefault(key, {"inc": set(), "exc": set()})
    entry["inc"].add(inc)
    entry["exc"].add(exc)
    entry["incJunc"], entry["excJunc"] = incJunc, excJunc


def _write_alt_events(altDict, chrom, strand, gName, outAS, out, isoform2reads):
    """Order the read sets by name and, with ``outAS``, write one A5SS/A3SS event line per boundary pair."""
    for alt, ev in altDict.items():
        ev["inc"] = sorted(ev["inc"], key=lambda x: x.name)
        ev["exc"] = sorted(ev["exc"], key=lambda x: x.name)
        if not outAS:
            continue
        altStart, altEnd = alt.split("-")
        inc, exc = [x.name for x in ev["inc"]], [x.name for x in ev["exc"]]
        print("\t".join(map(str, [chrom, altStart, altEnd, "%s:%s" % (gName, alt), _psi(inc, exc, isoform2reads), strand, gName,
                                  ",".join(inc), len(inc), ",".join(exc), len(exc),
                                  "inc:{}-{}".format(*ev["incJunc"]), "exc:{}-{}".format(*ev["excJunc"])])), file=out)


def _alt_start(altDict, read1, read2, i, j, offset):
    """Exon starts that differ (the boundary after an intron): the read whose exon starts earlier is the inclusion read,
    provided its start lies after the other read's previous exon."""
    s1, s2 = read1.exons[i][0], read2.exons[j][0]
    prev1, prev2 = read1.exons[i - 1][1], read2.exons[j - 1][1]
    if s2 > s1 > prev2 and abs(s1 - prev2) > offset:
        _record_pair(altDict, "%d-%d" % (s1, s2), read1, read2, incJunc=(prev1, s1), excJunc=(prev2, s2))
    elif s1 > s2 > prev1 and abs(s2 - prev1) > offset:
        _record_pair(altDict, "%d-%d" % (s2, s1), read2, read1, incJunc=(prev2, s2), excJunc=(prev1, s1))


def _alt_end(altDict, read1, read2, i, j, offset):
    """Exon ends that differ (the boundary before an intron): the read whose exon ends later is the inclusion read,
    provided its end lies before the other read's next exon."""
    e1, e2 = read1.exons[i][1], read2.exons[j][1]
    next1, next2 = read1.exons[i + 1][0], read2.exons[j + 1][0]
    if e1 < e2 < next1 and abs(e2 - next1) > offset:
        _record_pair(altDict, "%d-%d" % (e1, e2), read2, read1, incJunc=(e2, next2), excJunc=(e1, next1))
    elif e2 < e1 < next2 and abs(e1 - next2) > offset:
        _record_pair(altDict, "%d-%d" % (e2, e1), read1, read2, incJunc=(e1, next1), excJunc=(e2, next2))


def _alt_sites_between(altDict, read1, read2, at_start, offset):
    """Compare each overlapping pair of the two reads' exons that border an intron on the compared side and whose
    compared boundaries differ by more than ``offset``."""
    compare, side, lo = (_alt_start, 0, 1) if at_start else (_alt_end, 1, 0)
    exons1, exons2 = read1.exons, read2.exons
    for i in range(lo, len(exons1) - 1 + lo):
        for j in range(lo, len(exons2) - 1 + lo):
            if exons1[i][0] < exons2[j][1] and exons1[i][1] > exons2[j][0] and abs(exons1[i][side] - exons2[j][side]) > offset:
                compare(altDict, read1, read2, i, j, offset)
            if exons2[j][0] > exons1[i][1]:
                break


def _find_alt_splice_sites(asType, gene2ReadsDict, outAS=None, out=None, offset=0, isoform2reads=None):
    """A3SS / A5SS: for every pair of reads of a gene, an internal exon overlap whose alternative boundary lies between
    the two reads' flanking exons (both differences > ``offset``) is one event, keyed ``start-end`` of the alternative
    region; the read with the longer exon is the inclusion read. On '+' A3SS and '-' A5SS the alternative boundary is
    the exon *start* (compared with the previous exon's end), otherwise the exon *end* (compared with the next
    exon's start); the published code spelled the two cases out per strand."""
    for gName, gene in gene2ReadsDict.items():
        if len(gene.reads) < 2:
            continue
        at_start = (asType == "A3SS") == (gene.strand == "+")
        altDict = {}
        for read1, read2 in itertools.combinations(gene.reads.values(), 2):
            _alt_sites_between(altDict, read1, read2, at_start, offset)
        _write_alt_events(altDict, gene.chrom, gene.strand, gName, outAS, out, isoform2reads)
        gene.asDict[asType] = altDict


def findA3SS(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None):
    _find_alt_splice_sites("A3SS", gene2ReadsDict, outAS, out, offset, isoform2reads)


def findA5SS(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None):
    _find_alt_splice_sites("A5SS", gene2ReadsDict, outAS, out, offset, isoform2reads)


def _cluster_novel(novelDict):
    """Novel reads as fake genes ``chrom:strand:Novel_<n>``: per chromosome and strand, reads chained by overlapping spans."""
    genes = {}
    for chrom, byStrand in novelDict.items():
        for strand, reads in byStrand.items():
            count, clusterEnd = 0, None
            for read in sorted(reads, key=lambda x: x.chromStart):
                if clusterEnd is None or read.chromStart >= clusterEnd:
                    count += 1
                    name = "%s:%s:Novel_%d" % (chrom, strand, count)
                    gene = genes[name] = Gene2Reads(name)
                    clusterEnd = read.chromEnd
                gene.update(read)
                clusterEnd = max(clusterEnd, read.chromEnd)
    return genes


def find_single_as(asType="IR", gene2ReadsDict=None, anno=False, outFile=None, outAppend=False, isoform2reads=None):
    """Write the events of one AS type of annotated genes, or (``anno=False``) of novel reads clustered into fake genes;
    with ``asType=None`` all four detectors only fill ``asDict``."""
    findASType = {"IR": findIR, "SE": findSE, "A3SS": findA3SS, "A5SS": findA5SS}
    if asType and asType.upper() not in findASType:
        raise Exception("You should input the right AS type!")
    genes = gene2ReadsDict if anno else _cluster_novel(gene2ReadsDict)
    with open(outFile, "a" if outAppend else "w") as out:
        if asType:
            findASType[asType.upper()](genes, asType, out, anno=anno, offset=5, isoform2reads=isoform2reads)
        else:
            for detect in findASType.values():
                detect(genes, offset=5, isoform2reads=isoform2reads)

# ---------------------------------------------------------------------------
# HQ BED + LR driver: ``find_as_sample`` resolves the legacy sample-directory paths, ``find_as_lr`` builds
# isoform2reads, splits isoforms into annotated/novel against the reference and runs the four detectors.
# ---------------------------------------------------------------------------


def collapsed_gff_to_bed12plus(gff, group_file, out_path):
    """Bridge cupcake ``collapse`` output to the bed12+ that ``find_as_lr`` expects.

    Reads the collapsed GFF (1-based, inclusive) + group file, and writes one
    bed12+ line per collapsed isoform present in the group file. The trailing
    columns are ``<gene placeholder>\t<gene placeholder>\t<group key>`` so that
    ``otherList[2]`` (the 15th column) carries the cupcake isoform id used to
    look reads up in the group file — matching the column ``find_all_as`` reads.
    """
    group_keys = set()
    with open(group_file) as f:
        for line in f:
            if line.strip():
                group_keys.add(line.split("\t")[0])

    tx = {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', cols[8])
            if not m:
                continue
            name = m.group(1)
            rec = tx.setdefault(name, {"chrom": cols[0], "strand": cols[6], "exons": []})
            rec["exons"].append((int(cols[3]) - 1, int(cols[4])))   # -> 0-based half-open

    out = open(out_path, "w")
    for name, rec in tx.items():
        if group_keys and name not in group_keys:
            continue
        exons = sorted(rec["exons"])
        chromStart, chromEnd = exons[0][0], exons[-1][1]
        blockSizes = ",".join(str(e - s) for s, e in exons)
        blockStarts = ",".join(str(s - chromStart) for s, e in exons)
        fields = [rec["chrom"], chromStart, chromEnd, name, 0, rec["strand"],
                  chromStart, chromEnd, "0,0,0", len(exons), blockSizes, blockStarts,
                  ".", ".", name]
        print("\t".join(map(str, fields)), file=out)
    out.close()
    return out_path


def _find_utils_dir():
    env_dir = os.environ.get("IFLAS_UTILS_DIR")
    here = os.path.abspath(os.path.dirname(__file__))
    repo_root = os.path.abspath(os.path.join(here, "..", "..", ".."))
    candidates = [
        env_dir,
        os.path.join(repo_root, "utils"),
        os.path.join(os.getcwd(), "utils"),
    ]
    for path in candidates:
        if path and os.path.isdir(path):
            return path
    raise FileNotFoundError("Cannot locate iFLAS utils directory (set IFLAS_UTILS_DIR)")


def _filter_file_by_field(originFile, targetFile, originField=1, targetField=1, outFile=None):
    originTags = set()
    if os.path.isfile(originFile):
        with open(originFile) as f:
            for line in f:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= originField:
                    originTags.add(fields[originField - 1])

    outFile = outFile if outFile else "filtered.txt"
    out = open(outFile, "w")
    if os.path.isfile(targetFile):
        with open(targetFile) as f:
            for line in f:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= targetField and fields[targetField - 1] in originTags:
                    print(line.rstrip("\n"), file=out)
    out.close()
    return outFile


def _run_redirect(cmd, stdout_path, stderr_path=None, input_text=None, runner=None):
    runner = runner or ToolRunner()
    stdout = open(stdout_path, "w")
    stderr = open(stderr_path, "w") if stderr_path else None
    try:
        proc = runner.run(cmd, input=input_text, stdout=stdout, stderr=stderr, check=False)
    finally:
        stdout.close()
        if stderr_path:
            stderr.close()
    if proc.returncode != 0:
        err = ""
        if stderr_path and os.path.isfile(stderr_path):
            err = open(stderr_path).read()
        elif proc.stderr:
            err = proc.stderr
        raise RuntimeError("command failed (exit {}): {}\n{}".format(
            proc.returncode, " ".join(map(str, cmd)), err.strip()))
    return stdout_path


def _link_confident(lrDir, source_name, target_name):
    target = os.path.join(lrDir, target_name)
    makeLink(source_name, target)
    return target


def _confident_from_lr(outFiles, out_dir):
    lrDir = os.path.join(out_dir, "LR")
    return {
        "IR": _link_confident(lrDir, "IR.bed6+", "IR.confident.bed6+"),
        "SE": _link_confident(lrDir, "SE.bed12+", "SE.confident.bed12+"),
        "A5SS": _link_confident(lrDir, "A5SS.bed6+", "A5SS.confident.bed6+"),
        "A3SS": _link_confident(lrDir, "A3SS.bed6+", "A3SS.confident.bed6+"),
    }


def _with_assigned_name(junction_assigned):
    lines = []
    with open(junction_assigned) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 5:
                fields[3] = "{}:{}".format(fields[3], fields[5])
            lines.append("\t".join(fields))
    return "\n".join(lines) + ("\n" if lines else "")


def _event_isoforms(path, columns):
    """Isoform names in the given comma-separated columns of an event file (none when the file does not exist)."""
    names = set()
    if os.path.isfile(path):
        with open(path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                names.update(x for idx in columns if len(fields) > idx for x in fields[idx].split(",") if x)
    return names


def _write_confident_isoforms(confident, isoform_bed12plus, out_dir):
    selected = set().union(*(_event_isoforms(confident[t], (7, 9)) for t in ("A3SS", "A5SS", "IR")),
                           _event_isoforms(confident["SE"], (15, 17)))
    out_path = os.path.join(out_dir, "isoformGrouped.AS.confident.bed12+")
    with open(out_path, "w") as out, open(isoform_bed12plus) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 3 and fields[3] in selected:
                print(line.rstrip("\n"), file=out)
    return out_path


def apply_ngs_confident_filters(outFiles, isoform_bed12plus, ref_gpe, ngs_junctions,
                                out_dir, runner=None):
    """Restore legacy find_all_as NGS confident filtering for LR AS events."""
    from .scan_se_by_ngs import assignGeneNameToSE, scanSeByNGS

    runner = runner or ToolRunner()
    utilDir = _find_utils_dir()
    lrDir = os.path.join(out_dir, "LR")
    ngsDir = os.path.join(out_dir, "NGS")
    os.makedirs(ngsDir, exist_ok=True)

    se_support = os.path.join(out_dir, "SEsupportByNGS.bed12+")
    ngs_se = os.path.join(ngsDir, "SE.bed12+")
    try:
        scanSeByNGS(fref=ref_gpe, has_bin=False, fjunc=ngs_junctions,
                    flr=isoform_bed12plus, outFile=se_support)
        assignGeneNameToSE(fref=ref_gpe, has_bin=False, fes=se_support,
                           outFile=ngs_se,
                           errorOutFile=os.path.join(ngsDir, "novelSE2gene.log"))
    except Exception as e:
        print(e)
        open(se_support, "w").close()
        open(ngs_se, "w").close()
        open(os.path.join(ngsDir, "novelSE2gene.log"), "w").close()

    _filter_file_by_field(ngs_se, outFiles["SE"], originField=4, targetField=4,
                          outFile=os.path.join(lrDir, "SE.NGS.bed12+"))
    _filter_file_by_field(outFiles["SE"], ngs_se, originField=4, targetField=4,
                          outFile=os.path.join(ngsDir, "SE.LR.bed12+"))

    junction_assigned = os.path.join(out_dir, "junction.assigned.bed12+")
    _run_redirect(["perl", os.path.join(utilDir, "juncAssign.pl"), "-g", ref_gpe, ngs_junctions],
                  junction_assigned, runner=runner)

    _run_redirect(["perl", os.path.join(utilDir, "AnSSconfirmByJunc.pl"), "-t", "5",
                   "-j", junction_assigned, outFiles["A5SS"]],
                  os.path.join(ngsDir, "A5SS.LR.bed6+"),
                  stderr_path=os.path.join(lrDir, "A5SS.NGS.bed6+"),
                  runner=runner)
    _run_redirect(["perl", os.path.join(utilDir, "AnSSconfirmByJunc.pl"), "-t", "3",
                   "-j", junction_assigned, outFiles["A3SS"]],
                  os.path.join(ngsDir, "A3SS.LR.bed6+"),
                  stderr_path=os.path.join(lrDir, "A3SS.NGS.bed6+"),
                  runner=runner)

    assigned_input = _with_assigned_name(junction_assigned)
    _run_redirect(["perl", os.path.join(utilDir, "getNgsA5SS.pl"), "-e", "0"],
                  os.path.join(ngsDir, "A5SS.known.bed6+"),
                  stderr_path=os.path.join(ngsDir, "A5SS.novel.bed6+"),
                  input_text=assigned_input, runner=runner)
    _run_redirect(["perl", os.path.join(utilDir, "getNgsA3SS.pl"), "-e", "0"],
                  os.path.join(ngsDir, "A3SS.known.bed6+"),
                  stderr_path=os.path.join(ngsDir, "A3SS.novel.bed6+"),
                  input_text=assigned_input, runner=runner)

    filterIrByJunc(outFiles["IR"], ngs_junctions, os.path.join(lrDir, "IR.NGS.bed6+"))

    confident = {
        "IR": _link_confident(lrDir, "IR.NGS.bed6+", "IR.confident.bed6+"),
        "SE": _link_confident(lrDir, "SE.NGS.bed12+", "SE.confident.bed12+"),
        "A5SS": _link_confident(lrDir, "A5SS.NGS.bed6+", "A5SS.confident.bed6+"),
        "A3SS": _link_confident(lrDir, "A3SS.NGS.bed6+", "A3SS.confident.bed6+"),
    }
    return confident


def find_as_sample(base_dir, ref_gpe, out_dir=None, isoform_bed12plus=None,
                   group_file=None, offset=0, ngs_junctions=None, runner=None):
    """Run LR AS detection using the legacy sample-directory contract.

    Defaults match published ``find_all_as`` paths:
    ``<base_dir>/hqIsoforms/hq.collapsed.bed12+`` for high-quality isoforms,
    ``<base_dir>/collapse/tofu.collapsed.group.txt`` for isoform->read support,
    and ``<base_dir>/as_events/ordinary_as/LR`` for outputs.
    """
    base_dir = os.path.abspath(base_dir)
    isoform_bed12plus = isoform_bed12plus or os.path.join(
        base_dir, "hqIsoforms", "hq.collapsed.bed12+")
    group_file = group_file or os.path.join(
        base_dir, "collapse", "tofu.collapsed.group.txt")
    out_dir = out_dir or os.path.join(base_dir, "as_events", "ordinary_as")
    if ngs_junctions is None:
        detected = os.path.join(base_dir, "mapping", "rna-seq", "reassembly", "junctions.bed")
        if os.path.isfile(detected):
            ngs_junctions = detected

    for label, path in (("isoform bed12+", isoform_bed12plus),
                        ("collapse group file", group_file),
                        ("reference genePredExt", ref_gpe)):
        if not os.path.isfile(path):
            raise FileNotFoundError("{} not found: {}".format(label, path))

    return find_as_lr(isoform_bed12plus, group_file, ref_gpe, out_dir, offset=offset,
                      ngs_junctions=ngs_junctions, runner=runner)


def find_as_lr(isoform_bed12plus, group_file, ref_gpe, out_dir, offset=0,
               ngs_junctions=None, runner=None):
    """Identify ordinary AS events from long-read isoforms (the LR core of find_all_as).

    Faithfully reproduces the published flow: build the reference transcript BED,
    map each isoform to its supporting reads (group file), intersect isoforms with
    the reference to split anno/novel, then run the four detectors. Writes
    ``<out_dir>/LR/{IR.bed6+,SE.bed12+,A5SS.bed6+,A3SS.bed6+}``, plus
    ``*.confident.*`` links/filters matching the legacy contract, and returns the
    path map.
    """
    import pybedtools  # lazy: only the driver needs bedtools

    transBedList = GenePredObj(ref_gpe, bincolumn=False).toBed(gene=True)
    transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)

    tmpReads = {}
    with open(group_file) as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            tmpReads[isoform] = reads.split(",")

    isoform2reads = {}
    readsBedList = []
    with open(isoform_bed12plus) as f:
        for line in f.readlines():
            readStruc = Bed12Plus(line)
            isoform2reads[readStruc.name] = list(itertools.chain.from_iterable(
                [tmpReads[i] for i in readStruc.otherList[2].split(",")]))
            if len(isoform2reads[readStruc.name]) < 2: continue
            readsBedList.append("\t".join(readStruc.record[:12]))

    readsBedObj = pybedtools.BedTool("\n".join(readsBedList), from_string=True)
    annoBedRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
    novelBedRes = readsBedObj.intersect(transBedObj, v=True, s=True)
    annoDict, novelDict = getAnnotation(annoBedRes, novelBedRes, offset=offset)

    lrDir = os.path.join(out_dir, "LR")
    os.makedirs(lrDir, exist_ok=True)
    outFiles = {
        "IR": os.path.join(lrDir, "IR.bed6+"),
        "SE": os.path.join(lrDir, "SE.bed12+"),
        "A5SS": os.path.join(lrDir, "A5SS.bed6+"),
        "A3SS": os.path.join(lrDir, "A3SS.bed6+"),
    }
    for asType in ("IR", "SE", "A5SS", "A3SS"):
        find_single_as(asType=asType, gene2ReadsDict=annoDict, anno=True,
                       outFile=outFiles[asType], outAppend=False, isoform2reads=isoform2reads)
        find_single_as(asType=asType, gene2ReadsDict=novelDict, anno=False,
                       outFile=outFiles[asType], outAppend=True, isoform2reads=isoform2reads)
    if ngs_junctions:
        if not os.path.isfile(ngs_junctions):
            raise FileNotFoundError("NGS junctions not found: {}".format(ngs_junctions))
        confident = apply_ngs_confident_filters(outFiles, isoform_bed12plus, ref_gpe,
                                                ngs_junctions, out_dir, runner=runner)
    else:
        confident = _confident_from_lr(outFiles, out_dir)
    outFiles["confident"] = confident
    outFiles["isoformGrouped_confident"] = _write_confident_isoforms(
        confident, isoform_bed12plus, out_dir)
    return outFiles
