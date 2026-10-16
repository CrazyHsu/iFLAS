#!/usr/bin/env python
"""Fabricate ONT direct-RNA inputs for the nanopolish poly(A) step without an instrument: a multi-read fast5 file
(written with ont_fast5_api) whose raw signal is a simulated squiggle -- leader, adapter, a poly(A) stretch of the
designed length and the transcript body synthesised from nanopolish's own RNA 5-mer pore model (test_data/synthetic_v2/
r9.4_70bps.u_to_t_rna.5mer.template.model, time order 3' -> 5' as direct RNA translocates) -- plus the matching
basecalled FASTQ (read ids are the fast5 read ids), a sequencing_summary.txt and a truth table of designed tail lengths.

  python scripts/make_fake_fast5.py --flnc test_data/synthetic_v2/flnc.fa --model test_data/synthetic_v2/r9.4_70bps.u_to_t_rna.5mer.template.model \
      --out-dir fake_ont --reads 20
"""
from __future__ import annotations

import argparse
import os
import random
import uuid

import numpy as np


def write_fast5(records, model_path, out_dir, *, samples_per_base=40, seed=0):
    """``records`` = (read_id, transcript sequence 5'->3' without the tail, tail length) -> ``<out_dir>/fast5/batch0.fast5``
    (one multi-read file; signal: leader, adapter, poly(A), transcript 3'->5') and ``<out_dir>/sequencing_summary.txt``."""
    from ont_fast5_api.multi_fast5 import MultiFast5File

    rng, np_rng = random.Random(seed), np.random.RandomState(seed)
    model = {}
    with open(model_path) as f:
        for line in f:
            if line[0] in "#k" and not line.startswith(("A", "C", "G", "T")):
                continue
            kmer, level_mean, level_stdv = line.split("\t")[:3]
            model[kmer] = (float(level_mean), float(level_stdv))
    digitisation, offset, rng_pa, sampling_rate = 8192.0, 5.0, 1402.9, 3012.0
    scale = rng_pa / digitisation                                 # pA per raw unit

    def segment(mean_pa, stdv_pa, n):
        return np.clip(np.round(np_rng.normal(mean_pa, stdv_pa, n) / scale - offset), -32768, 32767).astype(np.int16)

    def transcript_signal(seq):
        seq = seq.upper().replace("U", "T")
        pieces = [segment(*model.get(seq[i:i + 5], (90.0, 3.0)), max(8, int(np_rng.poisson(samples_per_base))))
                  for i in reversed(range(len(seq) - 4))]           # 3' end enters the pore first
        return np.concatenate(pieces)

    os.makedirs(os.path.join(out_dir, "fast5"), exist_ok=True)
    run_id = uuid.UUID(int=rng.getrandbits(128)).hex
    polya_mean, polya_stdv = model["AAAAA"]
    with MultiFast5File(os.path.join(out_dir, "fast5", "batch0.fast5"), mode="w") as h5, \
            open(os.path.join(out_dir, "sequencing_summary.txt"), "w") as summary:
        summary.write("filename\tread_id\trun_id\tchannel\tstart_time\tduration\tnum_events\tsequence_length_template\tmean_qscore_template\n")
        for i, (read_id, seq, tail) in enumerate(records):
            signal = np.concatenate([segment(115.0, 2.0, 2000), segment(80.0, 3.0, 1500),
                                     segment(polya_mean, polya_stdv, tail * samples_per_base), transcript_signal(seq)])
            n_samples = len(signal)
            read = h5.create_empty_read(read_id, run_id)
            read.add_raw_data(signal, attrs={"duration": n_samples, "median_before": 250.0, "read_id": read_id,
                                             "read_number": i + 1, "start_mux": 1, "start_time": i * 100000})
            read.add_channel_info({"channel_number": str(i % 512 + 1), "digitisation": digitisation, "offset": offset,
                                   "range": rng_pa, "sampling_rate": sampling_rate})
            read.add_tracking_id({"run_id": run_id, "flow_cell_id": "FAK00000", "device_id": "MN00000", "exp_start_time": "2026-10-15T00:00:00Z"})
            read.add_context_tags({"sample_frequency": str(int(sampling_rate)), "experiment_type": "rna", "sequencing_kit": "sqk-rna002"})
            summary.write(f"batch0.fast5\t{read_id}\t{run_id}\t{i % 512 + 1}\t{i * 100000 / sampling_rate:.3f}\t{n_samples / sampling_rate:.3f}"
                          f"\t{n_samples // 10}\t{len(seq) + tail}\t12.0\n")
    return os.path.join(out_dir, "fast5")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--flnc", required=True, help="FLNC FASTA whose sequences become the basecalled reads")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--model", required=True, help="nanopolish RNA 5-mer pore model (kmer level_mean level_stdv ...)")
    ap.add_argument("--reads", type=int, default=20)
    ap.add_argument("--samples-per-base", type=int, default=40, help="~3012 Hz / 70 bases per s")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    rng = random.Random(a.seed)
    seqs, buf = [], []
    with open(a.flnc) as f:
        for line in f:
            if line.startswith(">"):
                if buf:
                    seqs.append("".join(buf))
                buf = []
            else:
                buf.append(line.strip())
    if buf:
        seqs.append("".join(buf))
    rng.shuffle(seqs)
    records = [(str(uuid.UUID(int=rng.getrandbits(128))), seq, rng.randint(30, 120)) for seq in seqs[:a.reads]]
    write_fast5(records, a.model, a.out_dir, samples_per_base=a.samples_per_base, seed=a.seed)
    with open(os.path.join(a.out_dir, "reads.fastq"), "w") as fq, open(os.path.join(a.out_dir, "tail_truth.tsv"), "w") as truth:
        truth.write("read_id\tpolya_len\tseq_len\n")
        for i, (read_id, seq, tail) in enumerate(records):
            basecalled = seq + "A" * tail                             # basecallers emit the tail (partly) as A's
            fq.write(f"@{read_id} read={i + 1} ch={i % 512 + 1}\n{basecalled}\n+\n{'~' * len(basecalled)}\n")
            truth.write(f"{read_id}\t{tail}\t{len(seq)}\n")
    print(f"{a.out_dir}: {len(records)} reads -> fast5/batch0.fast5, reads.fastq, sequencing_summary.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
