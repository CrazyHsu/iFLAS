"""Inclusion/exclusion isoform pairs of the confident AS events, per cupcake gene (shared by asas and palen_as;
the published code carried one copy in each)."""

from __future__ import annotations

import itertools

from ..io.files import getDictFromFile
from ..io.records import BedFile


def as_event_pairs(asFile, collapsedGroupFile, isoformFile, asType="SE", filterByCount=0, mergeByJunc=True):
    """Yield ``(gene, asEvent, inclusionIsos, exclusionIsos)`` for every confident event whose isoforms, after the
    FL-read count filter, all belong to one cupcake gene (``PB.x``). With ``mergeByJunc`` the isoforms are first
    grouped by junction chain (single-exon isoforms together; BUG-181: the published code pooled single-intron ones
    too) and one pair is yielded per inclusion-group x exclusion-group product; the caller decides whether several
    pairs of one event are kept or the last one wins."""
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    isoBed = BedFile(isoformFile, type="bed12+").reads

    def fl_count(iso):
        return sum(len(collapsedTrans2reads[y]) for y in isoBed[iso].otherList[2].split(","))

    def chain_groups(isos):
        groups = {}
        for iso in isos:
            groups.setdefault(isoBed[iso].juncChain if isoBed[iso].introns else "monoExon", []).append(iso)
        return list(groups.values())

    incCol, excCol = (15, 17) if asType == "SE" else (7, 9)
    with open(asFile) as f:
        for line in f:
            if not line.strip():
                continue
            records = line.rstrip("\n").split("\t")
            inclusionIsos, exclusionIsos = records[incCol].split(","), records[excCol].split(",")
            candidates = (itertools.product(chain_groups(inclusionIsos), chain_groups(exclusionIsos)) if mergeByJunc
                          else [(inclusionIsos, exclusionIsos)])
            for inc, exc in candidates:
                inc = [x for x in inc if fl_count(x) >= filterByCount]
                exc = [x for x in exc if fl_count(x) >= filterByCount]
                if not inc or not exc:
                    continue
                genes = {".".join(x.split(".")[:2]) for x in inc + exc}
                if len(genes) == 1:
                    yield genes.pop(), records[3], inc, exc
