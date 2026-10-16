#!/usr/bin/env python
"""Fabricate a PacBio Sequel subreads BAM from FLNC sequences so the ccs step of `iflas preproc` can be exercised
without real data: every ZMW is one insert (transcript + poly(A) + primers) read in several passes with random
errors, written with the PacBio read-group / tag conventions (subreads.bam spec).

  python scripts/make_fake_pacbio_subreads.py --flnc test_data/synthetic_v2/flnc.fa --out subreads.bam --zmws 40

The reads only exercise the commands: ccs 6.4 accepts them, so the ccs -> lima -> isoseq chain runs end to end, but
they say nothing about how ccs behaves on real instrument data.
"""
from __future__ import annotations

import argparse
import random
from array import array

import pysam

MOVIE = "m64000_000000_000000"
ADAPTER = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"                       # SMRTbell hairpin
PRIMER_5P = "AAGCAGTGGTATCAACGCAGAGTACATGGGG"
PRIMER_3P = "GTACTCTGCGTTGATACCACTGCTT"
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(rng, seq, rate):
    out = []
    for base in seq:
        r = rng.random()
        if r < rate / 3:
            continue                                      # deletion
        if r < 2 * rate / 3:
            out.append(rng.choice("ACGT"))                # insertion
        if r < rate:
            out.append(rng.choice("ACGT".replace(base, "")))
        else:
            out.append(base)
    return "".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--flnc", required=True, help="FLNC FASTA (poly(A)-trimmed inserts)")
    ap.add_argument("--out", required=True, help="output subreads.bam")
    ap.add_argument("--zmws", type=int, default=40)
    ap.add_argument("--passes", type=int, nargs=2, default=(6, 10), help="min/max passes per ZMW")
    ap.add_argument("--error-rate", type=float, default=0.03)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    rng = random.Random(a.seed)
    seqs = []
    with open(a.flnc) as f:
        name, buf = None, []
        for line in f:
            if line.startswith(">"):
                if buf:
                    seqs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if buf:
            seqs.append((name, "".join(buf)))
    rng.shuffle(seqs)
    header = {"HD": {"VN": "1.6", "SO": "unknown", "pb": "5.0.0"},
              "RG": [{"ID": "7e3a77a1", "PL": "PACBIO", "PU": MOVIE, "PM": "SEQUEL",
                      "DS": "READTYPE=SUBREAD;BINDINGKIT=101-789-500;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000"}]}
    truth = []
    with pysam.AlignmentFile(a.out, "wb", header=header) as bam:
        for zmw, (name, insert) in enumerate(seqs[:a.zmws], start=1):
            template = PRIMER_5P + insert + "A" * rng.randint(20, 60) + PRIMER_3P
            pos = 0
            n_pass = rng.randint(*a.passes)
            for k in range(n_pass):
                pos += len(ADAPTER)                                           # the adapter is not part of a subread
                read = mutate(rng, template if k % 2 == 0 else revcomp(template), a.error_rate)
                qs, qe = pos, pos + len(read)
                rec = pysam.AlignedSegment(bam.header)
                rec.query_name = f"{MOVIE}/{zmw}/{qs}_{qe}"
                rec.flag = 4
                rec.query_sequence = read
                rec.query_qualities = pysam.qualitystring_to_array("~" * len(read))
                cx = (1 if k > 0 else 0) | (2 if k < n_pass - 1 else 0)   # adapter before / after
                rec.set_tags([("RG", "7e3a77a1", "Z"), ("qs", qs, "i"), ("qe", qe, "i"), ("zm", zmw, "i"), ("cx", cx, "i"),
                              ("rq", 0.8, "f"), ("sn", [10.0, 6.0, 7.0, 11.0], "f")])
                # ccs >= 6 refuses subreads without kinetics: inter-pulse duration / pulse width in frames per base
                rec.set_tag("ip", array("B", [rng.randint(2, 40) for _ in read]))
                rec.set_tag("pw", array("B", [rng.randint(2, 25) for _ in read]))
                bam.write(rec)
                pos = qe
            truth.append((zmw, name, n_pass))
    with open(a.out + ".truth.tsv", "w") as f:
        f.write("zmw\tflnc_read\tpasses\n")
        for zmw, name, n in truth:
            f.write(f"{zmw}\t{name}\t{n}\n")
    print(f"{a.out}: {len(truth)} ZMWs, {sum(n for _, _, n in truth)} subreads")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
