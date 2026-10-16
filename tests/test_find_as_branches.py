"""Characterization tests of the A3SS/A5SS detectors on the eight minimal branch cases (strand x which read carries
the longer exon x event type), recorded from the pre-refactor code (Phase G G3d). They pin the exact event key,
inclusion/exclusion assignment, junction coordinates and output line, independently of the synthetic-v2 coverage."""

from __future__ import annotations

import io
import itertools

import pytest

from iflas.io.intervals import Bed12
from iflas.io.records import Gene2Reads
from iflas.stages import find_as

EXPECTED = {
    "A3SS:+:long_first": {
        "alt": {
            "200-220": {
                "inc": [
                    "long"
                ],
                "exc": [
                    "short"
                ],
                "incJunc": [
                    100,
                    200
                ],
                "excJunc": [
                    100,
                    220
                ]
            }
        },
        "lines": [
            "chr1\t200\t220\tg:200-220\t500.0\t+\tg\tlong\t1\tshort\t1\tinc:100-200\texc:100-220"
        ]
    },
    "A3SS:+:short_first": {
        "alt": {
            "200-220": {
                "inc": [
                    "long"
                ],
                "exc": [
                    "short"
                ],
                "incJunc": [
                    100,
                    200
                ],
                "excJunc": [
                    100,
                    220
                ]
            }
        },
        "lines": [
            "chr1\t200\t220\tg:200-220\t500.0\t+\tg\tlong\t1\tshort\t1\tinc:100-200\texc:100-220"
        ]
    },
    "A3SS:-:long_first": {
        "alt": {
            "80-100": {
                "inc": [
                    "long"
                ],
                "exc": [
                    "short"
                ],
                "incJunc": [
                    100,
                    200
                ],
                "excJunc": [
                    80,
                    200
                ]
            }
        },
        "lines": [
            "chr1\t80\t100\tg:80-100\t500.0\t-\tg\tlong\t1\tshort\t1\tinc:100-200\texc:80-200"
        ]
    },
    "A3SS:-:short_first": {
        "alt": {
            "80-100": {
                "inc": [
                    "long"
                ],
                "exc": [
                    "short"
                ],
                "incJunc": [
                    100,
                    200
                ],
                "excJunc": [
                    80,
                    200
                ]
            }
        },
        "lines": [
            "chr1\t80\t100\tg:80-100\t500.0\t-\tg\tlong\t1\tshort\t1\tinc:100-200\texc:80-200"
        ]
    },
    "A5SS:+:long_first": {
        "alt": {
            "80-100": {
                "inc": [
                    "long"
                ],
                "exc": [
                    "short"
                ],
                "incJunc": [
                    100,
                    200
                ],
                "excJunc": [
                    80,
                    200
                ]
            }
        },
        "lines": [
            "chr1\t80\t100\tg:80-100\t500.0\t+\tg\tlong\t1\tshort\t1\tinc:100-200\texc:80-200"
        ]
    },
    "A5SS:+:short_first": {
        "alt": {
            "80-100": {
                "inc": [
                    "long"
                ],
                "exc": [
                    "short"
                ],
                "incJunc": [
                    100,
                    200
                ],
                "excJunc": [
                    80,
                    200
                ]
            }
        },
        "lines": [
            "chr1\t80\t100\tg:80-100\t500.0\t+\tg\tlong\t1\tshort\t1\tinc:100-200\texc:80-200"
        ]
    },
    "A5SS:-:long_first": {
        "alt": {
            "200-220": {
                "inc": [
                    "long"
                ],
                "exc": [
                    "short"
                ],
                "incJunc": [
                    100,
                    200
                ],
                "excJunc": [
                    100,
                    220
                ]
            }
        },
        "lines": [
            "chr1\t200\t220\tg:200-220\t500.0\t-\tg\tlong\t1\tshort\t1\tinc:100-200\texc:100-220"
        ]
    },
    "A5SS:-:short_first": {
        "alt": {
            "200-220": {
                "inc": [
                    "long"
                ],
                "exc": [
                    "short"
                ],
                "incJunc": [
                    100,
                    200
                ],
                "excJunc": [
                    100,
                    220
                ]
            }
        },
        "lines": [
            "chr1\t200\t220\tg:200-220\t500.0\t-\tg\tlong\t1\tshort\t1\tinc:100-200\texc:100-220"
        ]
    }
}


def _bed12(name, strand, exons):
    start, end = exons[0][0], exons[-1][1]
    sizes = ",".join(str(e - s) for s, e in exons)
    starts = ",".join(str(s - start) for s, e in exons)
    return Bed12("\t".join(map(str, ["chr1", start, end, name, 0, strand, start, end, "0,0,0", len(exons), sizes, starts])))


@pytest.mark.parametrize("as_type,strand,long_first", list(itertools.product(("A3SS", "A5SS"), ("+", "-"), (True, False))))
def test_alt_splice_site_branches_match_the_recorded_behaviour(as_type, strand, long_first):
    if (as_type == "A3SS") == (strand == "+"):          # the alternative site is the second exon's start
        long_exons, short_exons = [(0, 100), (200, 300)], [(0, 100), (220, 300)]
    else:                                              # ... or the first exon's end
        long_exons, short_exons = [(0, 100), (200, 300)], [(0, 80), (200, 300)]
    reads = [("long", long_exons), ("short", short_exons)]
    if not long_first:
        reads.reverse()
    gene = Gene2Reads("g")
    for name, exons in reads:
        gene.update(_bed12(name, strand, exons))
    out = io.StringIO()
    getattr(find_as, "find" + as_type)({"g": gene}, outAS=True, out=out, anno=True, offset=0)
    got = {k: {"inc": sorted(x.name for x in v["inc"]), "exc": sorted(x.name for x in v["exc"]),
               "incJunc": list(v["incJunc"]), "excJunc": list(v["excJunc"])} for k, v in gene.asDict[as_type].items()}
    expected = EXPECTED[f"{as_type}:{strand}:{'long_first' if long_first else 'short_first'}"]
    assert got == expected["alt"]
    assert out.getvalue().splitlines() == expected["lines"]
