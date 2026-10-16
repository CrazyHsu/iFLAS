"""Characterization of find_as on hand-built inputs, recorded from the code before the G5 refactor (2026-10-15):
read-to-gene annotation (shared intron, single-exon span overlap, novel grouping), novel clustering for all four event
types on both strands, the ``offset`` filter, several retention reads of one intron, the detectors run without output
(``asType=None``), PSI over the FL reads of the isoforms and the annotation-overlap filter of IR/SE -- the branches
the synthetic and scale runs of the refactor gate do not reach."""
import copy
import json
import os

from iflas.io.intervals import Bed12
from iflas.io.records import Gene2Reads, ReadLineStruc
from iflas.stages import find_as

EXPECTED = json.loads(r'''{
 "anno": {
  "gA": [
   [
    "r_mono_in",
    "r_mono_out",
    "r_shared"
   ],
   [
    "t1"
   ]
  ]
 },
 "anno_A3SS": [
  "chr5\t550\t600\tgC:550-600\t214.28571428571428\t+\tgC\talt3\t1\talt5,inc,near5,ret,ret2,skip\t6\tinc:200-550\texc:203-600",
  "chr6\t400\t450\tgD:400-450\t500.0\t-\tgD\tret2\t1\tinc,ret\t2\tinc:450-600\texc:400-600",
  "chr6\t200\t250\tgD:200-250\t357.14285714285717\t-\tgD\talt5\t1\talt3,inc,skip\t3\tinc:250-600\texc:200-550",
  "chr6\t200\t400\tgD:200-400\t272.7272727272727\t-\tgD\tret\t1\talt3,skip\t2\tinc:400-600\texc:200-550",
  "chr6\t200\t450\tgD:200-450\t333.3333333333333\t-\tgD\tret2\t1\talt3,skip\t2\tinc:450-600\texc:200-550",
  "chr6\t250\t400\tgD:250-400\t375.0\t-\tgD\tret\t1\talt5\t1\tinc:400-600\texc:250-600",
  "chr6\t203\t400\tgD:203-400\t300.0\t-\tgD\tret\t1\tnear5\t1\tinc:400-600\texc:203-600",
  "chr6\t250\t450\tgD:250-450\t444.4444444444444\t-\tgD\tret2\t1\talt5\t1\tinc:450-600\texc:250-600",
  "chr6\t203\t450\tgD:203-450\t363.6363636363636\t-\tgD\tret2\t1\tnear5\t1\tinc:450-600\texc:203-600",
  "chr6\t203\t250\tgD:203-250\t416.6666666666667\t-\tgD\talt5\t1\tnear5\t1\tinc:250-600\texc:203-600"
 ],
 "anno_A5SS": [
  "chr5\t400\t450\tgC:400-450\t500.0\t+\tgC\tret2\t1\tinc,ret\t2\tinc:450-600\texc:400-600",
  "chr5\t200\t250\tgC:200-250\t357.14285714285717\t+\tgC\talt5\t1\talt3,inc,skip\t3\tinc:250-600\texc:200-550",
  "chr5\t200\t400\tgC:200-400\t272.7272727272727\t+\tgC\tret\t1\talt3,skip\t2\tinc:400-600\texc:200-550",
  "chr5\t200\t450\tgC:200-450\t333.3333333333333\t+\tgC\tret2\t1\talt3,skip\t2\tinc:450-600\texc:200-550",
  "chr5\t250\t400\tgC:250-400\t375.0\t+\tgC\tret\t1\talt5\t1\tinc:400-600\texc:250-600",
  "chr5\t203\t400\tgC:203-400\t300.0\t+\tgC\tret\t1\tnear5\t1\tinc:400-600\texc:203-600",
  "chr5\t250\t450\tgC:250-450\t444.4444444444444\t+\tgC\tret2\t1\talt5\t1\tinc:450-600\texc:250-600",
  "chr5\t203\t450\tgC:203-450\t363.6363636363636\t+\tgC\tret2\t1\tnear5\t1\tinc:450-600\texc:203-600",
  "chr5\t203\t250\tgC:203-250\t416.6666666666667\t+\tgC\talt5\t1\tnear5\t1\tinc:250-600\texc:203-600",
  "chr6\t550\t600\tgD:550-600\t214.28571428571428\t-\tgD\talt3\t1\talt5,inc,near5,ret,ret2,skip\t6\tinc:200-550\texc:203-600"
 ],
 "anno_IR": [
  "chr5\t200\t300\tgC:200-300\t875.0\t+\tgC\tret,ret2\t2\tinc\t1"
 ],
 "anno_SE": [
  "chr5\t300\t400\tgC:200@300-400@600\t333.3333333333333\t+\t300\t400\t0,0,0\t1\t100\t0\tgC\t200\t600\tinc\t1\tskip\t1"
 ],
 "anno_asdict_A3SS": {
  "gC": {
   "A3SS": {
    "550-600": {
     "exc": [
      "alt5",
      "inc",
      "near5",
      "ret",
      "ret2",
      "skip"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "alt3"
     ],
     "incJunc": [
      200,
      550
     ]
    }
   },
   "A5SS": {},
   "IR": {},
   "SE": {}
  },
  "gD": {
   "A3SS": {
    "200-250": {
     "exc": [
      "alt3",
      "inc",
      "skip"
     ],
     "excJunc": [
      200,
      550
     ],
     "inc": [
      "alt5"
     ],
     "incJunc": [
      250,
      600
     ]
    },
    "200-400": {
     "exc": [
      "alt3",
      "skip"
     ],
     "excJunc": [
      200,
      550
     ],
     "inc": [
      "ret"
     ],
     "incJunc": [
      400,
      600
     ]
    },
    "200-450": {
     "exc": [
      "alt3",
      "skip"
     ],
     "excJunc": [
      200,
      550
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    },
    "203-250": {
     "exc": [
      "near5"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "alt5"
     ],
     "incJunc": [
      250,
      600
     ]
    },
    "203-400": {
     "exc": [
      "near5"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "ret"
     ],
     "incJunc": [
      400,
      600
     ]
    },
    "203-450": {
     "exc": [
      "near5"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    },
    "250-400": {
     "exc": [
      "alt5"
     ],
     "excJunc": [
      250,
      600
     ],
     "inc": [
      "ret"
     ],
     "incJunc": [
      400,
      600
     ]
    },
    "250-450": {
     "exc": [
      "alt5"
     ],
     "excJunc": [
      250,
      600
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    },
    "400-450": {
     "exc": [
      "inc",
      "ret"
     ],
     "excJunc": [
      400,
      600
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    }
   },
   "A5SS": {},
   "IR": {},
   "SE": {}
  }
 },
 "anno_asdict_A5SS": {
  "gC": {
   "A3SS": {},
   "A5SS": {
    "200-250": {
     "exc": [
      "alt3",
      "inc",
      "skip"
     ],
     "excJunc": [
      200,
      550
     ],
     "inc": [
      "alt5"
     ],
     "incJunc": [
      250,
      600
     ]
    },
    "200-400": {
     "exc": [
      "alt3",
      "skip"
     ],
     "excJunc": [
      200,
      550
     ],
     "inc": [
      "ret"
     ],
     "incJunc": [
      400,
      600
     ]
    },
    "200-450": {
     "exc": [
      "alt3",
      "skip"
     ],
     "excJunc": [
      200,
      550
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    },
    "203-250": {
     "exc": [
      "near5"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "alt5"
     ],
     "incJunc": [
      250,
      600
     ]
    },
    "203-400": {
     "exc": [
      "near5"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "ret"
     ],
     "incJunc": [
      400,
      600
     ]
    },
    "203-450": {
     "exc": [
      "near5"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    },
    "250-400": {
     "exc": [
      "alt5"
     ],
     "excJunc": [
      250,
      600
     ],
     "inc": [
      "ret"
     ],
     "incJunc": [
      400,
      600
     ]
    },
    "250-450": {
     "exc": [
      "alt5"
     ],
     "excJunc": [
      250,
      600
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    },
    "400-450": {
     "exc": [
      "inc",
      "ret"
     ],
     "excJunc": [
      400,
      600
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    }
   },
   "IR": {},
   "SE": {}
  },
  "gD": {
   "A3SS": {},
   "A5SS": {
    "550-600": {
     "exc": [
      "alt5",
      "inc",
      "near5",
      "ret",
      "ret2",
      "skip"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "alt3"
     ],
     "incJunc": [
      200,
      550
     ]
    }
   },
   "IR": {},
   "SE": {}
  }
 },
 "anno_asdict_IR": {
  "gC": {
   "A3SS": {},
   "A5SS": {},
   "IR": {
    "200-300": {
     "retention": [
      "ret",
      "ret2"
     ],
     "spliced": [
      "inc"
     ]
    }
   },
   "SE": {}
  },
  "gD": {
   "A3SS": {},
   "A5SS": {},
   "IR": {
    "200-300": {
     "retention": [
      "ret",
      "ret2"
     ],
     "spliced": [
      "inc"
     ]
    }
   },
   "SE": {}
  }
 },
 "anno_asdict_SE": {
  "gC": {
   "A3SS": {},
   "A5SS": {},
   "IR": {},
   "SE": {
    "200@300-400@600": {
     "keep": [
      "inc"
     ],
     "skip": [
      "skip"
     ]
    }
   }
  },
  "gD": {
   "A3SS": {},
   "A5SS": {},
   "IR": {},
   "SE": {
    "200@300-400@600": {
     "keep": [
      "inc"
     ],
     "skip": [
      "skip"
     ]
    }
   }
  }
 },
 "anno_none_asdict": {
  "gA": {
   "A3SS": {},
   "A5SS": {},
   "IR": {},
   "SE": {}
  },
  "gC": {
   "A3SS": {
    "550-600": {
     "exc": [
      "alt5",
      "inc",
      "near5",
      "ret",
      "ret2",
      "skip"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "alt3"
     ],
     "incJunc": [
      200,
      550
     ]
    }
   },
   "A5SS": {
    "200-250": {
     "exc": [
      "alt3",
      "inc",
      "skip"
     ],
     "excJunc": [
      200,
      550
     ],
     "inc": [
      "alt5"
     ],
     "incJunc": [
      250,
      600
     ]
    },
    "200-400": {
     "exc": [
      "alt3",
      "skip"
     ],
     "excJunc": [
      200,
      550
     ],
     "inc": [
      "ret"
     ],
     "incJunc": [
      400,
      600
     ]
    },
    "200-450": {
     "exc": [
      "alt3",
      "skip"
     ],
     "excJunc": [
      200,
      550
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    },
    "203-250": {
     "exc": [
      "near5"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "alt5"
     ],
     "incJunc": [
      250,
      600
     ]
    },
    "203-400": {
     "exc": [
      "near5"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "ret"
     ],
     "incJunc": [
      400,
      600
     ]
    },
    "203-450": {
     "exc": [
      "near5"
     ],
     "excJunc": [
      203,
      600
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    },
    "250-400": {
     "exc": [
      "alt5"
     ],
     "excJunc": [
      250,
      600
     ],
     "inc": [
      "ret"
     ],
     "incJunc": [
      400,
      600
     ]
    },
    "250-450": {
     "exc": [
      "alt5"
     ],
     "excJunc": [
      250,
      600
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    },
    "400-450": {
     "exc": [
      "inc",
      "ret"
     ],
     "excJunc": [
      400,
      600
     ],
     "inc": [
      "ret2"
     ],
     "incJunc": [
      450,
      600
     ]
    }
   },
   "IR": {
    "200-300": {
     "retention": [
      "ret",
      "ret2"
     ],
     "spliced": [
      "inc"
     ]
    }
   },
   "SE": {
    "200@300-400@600": {
     "keep": [
      "inc"
     ],
     "skip": [
      "skip"
     ]
    }
   }
  }
 },
 "anno_none_written": "",
 "novel": {
  "chr1": {
   "+": [
    "n_plus",
    "r_novel",
    "r_mono_edge"
   ],
   "-": [
    "n_minus"
   ]
  },
  "chr2": {
   "+": [
    "n_chr2"
   ]
  }
 },
 "novel_A3SS": [
  "chr3\t550\t600\tchr3:+:Novel_1:550-600\t142.85714285714286\t+\tchr3:+:Novel_1\talt3\t1\talt5,inc,near5,ret,ret2,skip\t6\tinc:200-550\texc:203-600",
  "chr3\t200\t250\tchr3:-:Novel_1:200-250\t333.3333333333333\t-\tchr3:-:Novel_1\tm3\t1\tm1,m4\t2\tinc:250-300\texc:200-350"
 ],
 "novel_A5SS": [
  "chr3\t400\t450\tchr3:+:Novel_1:400-450\t333.3333333333333\t+\tchr3:+:Novel_1\tret2\t1\tinc,ret\t2\tinc:450-600\texc:400-600",
  "chr3\t200\t450\tchr3:+:Novel_1:200-450\t333.3333333333333\t+\tchr3:+:Novel_1\tret2\t1\talt3,skip\t2\tinc:450-600\texc:200-550",
  "chr3\t250\t450\tchr3:+:Novel_1:250-450\t500.0\t+\tchr3:+:Novel_1\tret2\t1\talt5\t1\tinc:450-600\texc:250-600",
  "chr3\t203\t450\tchr3:+:Novel_1:203-450\t500.0\t+\tchr3:+:Novel_1\tret2\t1\tnear5\t1\tinc:450-600\texc:203-600",
  "chr3\t200\t250\tchr3:+:Novel_1:200-250\t250.0\t+\tchr3:+:Novel_1\talt5\t1\talt3,inc,skip\t3\tinc:250-600\texc:200-550",
  "chr3\t200\t400\tchr3:+:Novel_1:200-400\t333.3333333333333\t+\tchr3:+:Novel_1\tret\t1\talt3,skip\t2\tinc:400-600\texc:200-550",
  "chr3\t250\t400\tchr3:+:Novel_1:250-400\t500.0\t+\tchr3:+:Novel_1\tret\t1\talt5\t1\tinc:400-600\texc:250-600",
  "chr3\t203\t400\tchr3:+:Novel_1:203-400\t500.0\t+\tchr3:+:Novel_1\tret\t1\tnear5\t1\tinc:400-600\texc:203-600",
  "chr3\t203\t250\tchr3:+:Novel_1:203-250\t500.0\t+\tchr3:+:Novel_1\talt5\t1\tnear5\t1\tinc:250-600\texc:203-600",
  "chr3\t300\t350\tchr3:-:Novel_1:300-350\t666.6666666666666\t-\tchr3:-:Novel_1\tm1,m3\t2\tm4\t1\tinc:250-300\texc:200-350"
 ],
 "novel_IR": [
  "chr3\t200\t300\tchr3:+:Novel_1:200-300\t666.6666666666666\t+\tchr3:+:Novel_1\tret2,ret\t2\tinc\t1",
  "chr3\t200\t300\tchr3:-:Novel_1:200-300\t500.0\t-\tchr3:-:Novel_1\tm2\t1\tm1\t1",
  "chr3\t250\t300\tchr3:-:Novel_1:250-300\t500.0\t-\tchr3:-:Novel_1\tm2\t1\tm3\t1",
  "chr3\t200\t350\tchr3:-:Novel_1:200-350\t500.0\t-\tchr3:-:Novel_1\tm2\t1\tm4\t1"
 ],
 "novel_SE": [
  "chr3\t300\t400\tchr3:+:Novel_1:200@300-400@600\t500.0\t+\t300\t400\t0,0,0\t1\t100\t0\tchr3:+:Novel_1\t200\t600\tinc\t1\tskip\t1"
 ],
 "novel_append": [
  "header",
  "chr3\t200\t300\tchr3:+:Novel_1:200-300\t666.6666666666666\t+\tchr3:+:Novel_1\tret2,ret\t2\tinc\t1",
  "chr3\t200\t300\tchr3:-:Novel_1:200-300\t500.0\t-\tchr3:-:Novel_1\tm2\t1\tm1\t1",
  "chr3\t250\t300\tchr3:-:Novel_1:250-300\t500.0\t-\tchr3:-:Novel_1\tm2\t1\tm3\t1",
  "chr3\t200\t350\tchr3:-:Novel_1:200-350\t500.0\t-\tchr3:-:Novel_1\tm2\t1\tm4\t1"
 ],
 "one_read_asdict": {
  "g1": {
   "A3SS": {},
   "A5SS": {},
   "IR": {},
   "SE": {}
  }
 },
 "psi_IR_annoFalse": [
  "chr4\t200\t300\tgB:200-300\t600.0\t+\tgB\tc\t1\ta\t1"
 ],
 "psi_IR_annoTrue": [],
 "psi_SE_annoFalse": [
  "chr4\t300\t400\tgB:200@300-400@600\t666.6666666666666\t+\t300\t400\t0,0,0\t1\t100\t0\tgB\t200\t600\ta\t1\tb\t1"
 ],
 "psi_SE_annoTrue": []
}''')


def bed12(chrom, name, strand, exons):
    s, e = exons[0][0], exons[-1][1]
    return "\t".join(map(str, [chrom, s, e, name, 0, strand, s, e, "0,0,0", len(exons),
                               ",".join(str(b - a) for a, b in exons), ",".join(str(a - s) for a, _ in exons)]))


def gene(name, chrom, strand, isoforms, transcripts=()):
    g = Gene2Reads(name)
    for n, ex in isoforms:
        g.update(Bed12(bed12(chrom, n, strand, ex)))
    g.trans = {t: ReadLineStruc(bed12(chrom, t, strand, ex) + "\t" + name) for t, ex in transcripts}
    return g


def as_dict(genes):
    """asDict with reads as names and junction tuples as lists."""
    def plain(v):
        if isinstance(v, (list, tuple, set)):
            return [plain(x) for x in v]
        return getattr(v, "name", v)
    return {g: {t: {k: {c: plain(x) for c, x in ev.items()} for k, ev in events.items()} for t, events in v.asDict.items()}
            for g, v in genes.items()}


SE_IR_ALT = [("inc", [(100, 200), (300, 400), (600, 700)]), ("skip", [(100, 200), (600, 700)]),
             ("ret", [(100, 400), (600, 700)]), ("ret2", [(50, 450), (600, 700)]), ("alt5", [(100, 250), (600, 700)]),
             ("alt3", [(100, 200), (550, 700)]), ("near5", [(100, 203), (600, 700)])]          # 3 nt <= offset 5: no event


def snapshot(out_dir):
    T1 = [(1000, 1200), (1500, 1700), (2000, 2300)]                     # transcript t1 of gene gA (+)
    reads = {
        "r_shared": [(1050, 1200), (1500, 1700)],                       # shares intron 1200-1500 -> annotated
        "r_novel": [(1050, 1250), (1550, 1700)],                         # no shared intron -> novel
        "r_mono_in": [(1520, 1690)],                                     # single exon inside t1 exon 2 -> annotated
        "r_mono_out": [(1750, 1950)],                                    # single exon in t1's intron: the span counts -> annotated
        "r_mono_edge": [(2200, 2600)],                                   # single exon, 25 % inside t1's span -> novel
    }
    anno_rows = [bed12("chr1", n, "+", ex) + "\t" + bed12("chr1", "t1", "+", T1) + "\tgA" for n, ex in reads.items()]
    novel_rows = [bed12("chr1", "n_plus", "+", [(5000, 5100), (5300, 5400)]), bed12("chr1", "n_minus", "-", [(6000, 6100)]),
                  bed12("chr2", "n_chr2", "+", [(100, 300)])]
    annoDict, novelDict = find_as.getAnnotation(anno_rows, novel_rows, offset=0)
    snap = {"anno": {g: [sorted(v.reads), sorted(v.trans)] for g, v in annoDict.items()},
            "novel": {c: {s: [r.name for r in rs] for s, rs in d.items()} for c, d in novelDict.items()}}

    # find_single_as on novel clusters (anno=False): two clusters on '+', '-' strand alternative sites both ways
    nov = {"chr3": {"+": [Bed12(bed12("chr3", n, "+", ex)) for n, ex in SE_IR_ALT + [("far", [(5000, 5100), (5200, 5300)]),
                                                                                    ("far2", [(5050, 5100), (5200, 5300)])]],
                    "-": [Bed12(bed12("chr3", n, "-", ex)) for n, ex in [("m1", [(100, 200), (300, 400)]), ("m2", [(100, 400)]),
                                                                        ("m3", [(100, 250), (300, 400)]), ("m4", [(100, 200), (350, 400)])]]}}
    for t in ("SE", "IR", "A5SS", "A3SS"):
        path = os.path.join(out_dir, t)
        find_as.find_single_as(asType=t, gene2ReadsDict=copy.deepcopy(nov), anno=False, outFile=path)
        snap["novel_" + t] = open(path).read().splitlines()
    path = os.path.join(out_dir, "novel_append")
    open(path, "w").write("header\n")
    find_as.find_single_as(asType="IR", gene2ReadsDict=copy.deepcopy(nov), anno=False, outFile=path, outAppend=True)
    snap["novel_append"] = open(path).read().splitlines()

    # anno=True with no asType: every detector runs, nothing is written; asDict holds the events
    g = {"gA": copy.deepcopy(annoDict["gA"]), "gC": gene("gC", "chr5", "+", SE_IR_ALT, [("tc", [(100, 700)])])}
    path = os.path.join(out_dir, "none")
    find_as.find_single_as(asType=None, gene2ReadsDict=g, anno=True, outFile=path)
    snap["anno_none_written"] = open(path).read()
    snap["anno_none_asdict"] = as_dict(g)

    # annotated genes with one AS type: written only where an event overlaps a transcript of the gene
    for t in ("SE", "IR", "A5SS", "A3SS"):
        g = {"gC": gene("gC", "chr5", "+", SE_IR_ALT, [("tc", [(100, 700)])]),
             "gD": gene("gD", "chr6", "-", SE_IR_ALT, [("tfar", [(9000, 9100), (9200, 9300)])])}
        path = os.path.join(out_dir, "anno_" + t)
        find_as.find_single_as(asType=t, gene2ReadsDict=g, anno=True, outFile=path, isoform2reads={
            n: [n + str(i) for i in range(k)] for k, (n, _) in enumerate(SE_IR_ALT, 1)})
        snap["anno_" + t] = open(path).read().splitlines()
        snap["anno_asdict_" + t] = as_dict(g)

    # isoform2reads PSI and the annotation-overlap filter of IR/SE (transcripts far away -> no line written)
    gB = gene("gB", "chr4", "+", [("a", [(100, 200), (300, 400), (600, 700)]), ("b", [(100, 200), (600, 700)]),
                                  ("c", [(100, 400), (600, 700)])], [("tfar", [(9000, 9100), (9200, 9300)])])
    for anno in (True, False):
        for t, fn in (("IR", find_as.findIR), ("SE", find_as.findSE)):
            path = os.path.join(out_dir, f"{t}_{anno}")
            with open(path, "w") as out:
                fn({"gB": copy.deepcopy(gB)}, outAS=True, out=out, anno=anno, offset=0, isoform2reads={"a": ["x1", "x2"], "b": ["x3"], "c": ["x4", "x5", "x6"]})
            snap[f"psi_{t}_anno{anno}"] = open(path).read().splitlines()
    # a gene with one read: no event and asDict left as created
    one = {"g1": gene("g1", "chr7", "+", [("solo", [(100, 200), (300, 400)])])}
    for fn in (find_as.findIR, find_as.findSE, find_as.findA5SS, find_as.findA3SS):
        fn(one, offset=0)
    snap["one_read_asdict"] = as_dict(one)
    return json.loads(json.dumps(snap, sort_keys=True))


def test_find_as_matches_the_pre_refactor_record(tmp_path):
    assert snapshot(str(tmp_path)) == EXPECTED
