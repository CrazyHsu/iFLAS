"""Preprocessing stage (legacy ``preprocess.py`` high-level driver).

This ports the command orchestration for raw PacBio/Nanopore/NGS preprocessing to
the packaged CLI. External tools remain external and are run through ``ToolRunner``;
the science-sensitive command choices are kept aligned with the legacy code and
documented CLI delta notes.
"""

from __future__ import annotations

import gzip
import glob
from collections import Counter
import os
import re
import shutil
from pathlib import Path

from ..config import CcsParams, Sample
from ..io.files import makeLink, validateDir, validateFile
from ..runtime import ToolRunner, getCurrentTime, get_logger, quote as q

_PACBIO_ENV_MAP = {
    "bax2bam": "iflas-pacbio",
    "ccs": "iflas-pacbio",
    "dataset": "iflas-pacbio",
    "isoseq": "iflas-pacbio",
    "lima": "iflas-pacbio",
    "pbindex": "iflas-pacbio",
}


def _runner(runner=None, env_map=None):
    return runner if runner is not None else ToolRunner(env_map)


def _has_ngs(sample: Sample) -> bool:
    return bool(sample.ngs_left_reads or sample.ngs_right_reads)


def _split_repeats(value: str | None) -> list[list[str]]:
    if not value:
        return []
    repeats = []
    for repeat in str(value).split(";"):
        reads = [x.strip() for x in repeat.split(",") if x.strip()]
        if reads:
            repeats.append(reads)
    return repeats


def _join_repeats(repeats: list[list[str]]) -> str | None:
    if not repeats:
        return None
    return ";".join(",".join(reads) for reads in repeats)


def _flatten_reads(value: str | None) -> list[str]:
    if not value:
        return []
    return [x.strip() for x in re.split("[,;]", str(value)) if x.strip()]


def _fastp_dir(base_dir: str) -> str:
    return os.path.join(base_dir, "preprocess", "rna-seq")


def _fastp_targets(reads: list[str], out_dir: str) -> dict[str, str]:
    """fastp output path per input read: ``<out_dir>/fastp.<basename>``; inputs of one sample that share a basename
    get their parent directory in the name."""
    names = Counter(os.path.basename(r) for r in reads)
    out = {}
    for read in reads:
        base = os.path.basename(read)
        if names[base] > 1:
            base = os.path.basename(os.path.dirname(os.path.abspath(read))) + "." + base
        out[read] = os.path.join(out_dir, "fastp." + base)
    return out


def _existing_fastp(read: str, target: str, project_dir: str) -> str | None:
    """An already filtered copy of ``read``: this sample's target, the legacy copy beside the input, or a sibling
    sample's product whose ``.source`` sidecar names the same input (a library shared by samples is filtered once)."""
    legacy = os.path.join(os.path.dirname(read), "fastp." + os.path.basename(read))
    for cand in (target, legacy):
        if validateFile(cand):
            return cand
    source = os.path.abspath(read)
    for cand in sorted(glob.glob(os.path.join(project_dir, "*", "preprocess", "rna-seq", "fastp.*.source"))):
        product = cand[:-len(".source")]
        with open(cand) as f:
            if f.read().strip() == source and validateFile(product):
                return product
    return None


def _trimmed_length(sample: Sample, max_reads_length_trimmed: int) -> int:
    if sample.ngs_reads_length is None:
        raise ValueError("preproc fastp requires ngs_reads_length in samples.tsv")
    return int(sample.ngs_reads_length) - int(max_reads_length_trimmed)


def process_rnaseq(
    sample: Sample,
    base_dir: str,
    *,
    threads: int,
    runner=None,
    max_reads_length_trimmed: int = 1,
) -> dict:
    """Filter the short reads with fastp and point the sample at the filtered files.

    Legacy repeat syntax is preserved (semicolons separate repeats, commas files within a repeat). Products go to
    ``<base_dir>/preprocess/rna-seq/fastp.<name>`` (legacy wrote beside the raw reads, BUG-152) with a ``.source``
    sidecar; inputs already filtered (``fastp`` in the name, a legacy copy beside the input, or another sample's
    product of the same file) are reused.
    """
    runner = _runner(runner)
    log_dir = os.path.join(base_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    out_dir = _fastp_dir(base_dir)
    project_dir = os.path.dirname(os.path.abspath(base_dir))
    paired = (sample.ngs_reads_paired or "paired").lower() == "paired"
    left_repeats, right_repeats = _split_repeats(sample.ngs_left_reads), _split_repeats(sample.ngs_right_reads)
    paired = paired and bool(left_repeats) and bool(right_repeats)
    if paired and (len(left_repeats) != len(right_repeats) or any(len(l) != len(r) for l, r in zip(left_repeats, right_repeats))):
        raise ValueError("paired NGS repeats require equal left/right counts")
    reads = [r for repeat in left_repeats + right_repeats for r in repeat]
    if not reads or all("fastp" in os.path.basename(r) for r in reads):
        return {"ngs_left_reads": sample.ngs_left_reads, "ngs_right_reads": sample.ngs_right_reads}
    trim_len = _trimmed_length(sample, max_reads_length_trimmed)
    targets = _fastp_targets(reads, out_dir)
    resolved = {r: _existing_fastp(r, targets[r], project_dir) for r in reads}

    def options(target):
        # QC report and messages per filtered file (one shared fastp.json/html was overwritten by every repeat; BUG-188)
        report = os.path.join(log_dir, os.path.basename(target))
        return f"-w {threads} -q 20 -l {trim_len} -j {q(report + '.json')} -h {q(report + '.html')} 2> {q(report + '.log')}"

    def run(cmd, produced):
        os.makedirs(out_dir, exist_ok=True)
        runner.run(cmd)
        for read, path in produced:
            with open(path + ".source", "w") as f:
                f.write(os.path.abspath(read) + "\n")
            resolved[read] = path

    if paired:
        for lefts, rights in zip(left_repeats, right_repeats):
            for left, right in zip(lefts, rights):
                if not (resolved[left] and resolved[right]):
                    run(f"fastp -i {q(left)} -I {q(right)} -o {q(targets[left])} -O {q(targets[right])} {options(targets[left])}",
                        [(left, targets[left]), (right, targets[right])])
    else:
        for read in reads:
            if not resolved[read]:
                run(f"fastp -i {q(read)} -o {q(targets[read])} {options(targets[read])}", [(read, targets[read])])
    sample.ngs_left_reads = _join_repeats([[resolved[r] for r in repeat] for repeat in left_repeats])
    sample.ngs_right_reads = _join_repeats([[resolved[r] for r in repeat] for repeat in right_repeats])
    return {"ngs_left_reads": sample.ngs_left_reads, "ngs_right_reads": sample.ngs_right_reads}


def correct_with_fmlrc2(
    sample: Sample,
    base_dir: str,
    *,
    threads: int,
    runner=None,
    use_fmlrc2: bool | None = None,
) -> str | None:
    """Correct long reads with NGS reads via fmlrc2, preserving legacy commands."""
    use_fmlrc2 = sample.use_fmlrc2 if use_fmlrc2 is None else use_fmlrc2
    if not use_fmlrc2:
        return sample.data_processed_location
    if not _has_ngs(sample):
        raise ValueError("fmlrc2 correction requires NGS reads")

    runner = _runner(runner)
    work_dir = os.path.join(base_dir, "preprocess", "fmlrc")
    log_dir = os.path.join(base_dir, "log")
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    if sample.ngs_left_reads and sample.ngs_right_reads:
        left_reads = _flatten_reads(sample.ngs_left_reads)
        right_reads = _flatten_reads(sample.ngs_right_reads)
        runner.run("zcat -f {} > leftReads.fastq".format(" ".join(map(q, left_reads))), cwd=work_dir)
        runner.run("zcat -f {} > rightReads.fastq".format(" ".join(map(q, right_reads))), cwd=work_dir)
        seq_cmd = "cat leftReads.fastq rightReads.fastq | awk 'NR % 4 == 2' | tr NT TN | "
    elif sample.ngs_left_reads:
        left_reads = _flatten_reads(sample.ngs_left_reads)
        runner.run("zcat -f {} > leftReads.fastq".format(" ".join(map(q, left_reads))), cwd=work_dir)
        seq_cmd = "cat leftReads.fastq | awk 'NR % 4 == 2' | tr NT TN | "
    elif sample.ngs_right_reads:
        right_reads = _flatten_reads(sample.ngs_right_reads)
        runner.run("zcat -f {} > rightReads.fastq".format(" ".join(map(q, right_reads))), cwd=work_dir)
        seq_cmd = "cat rightReads.fastq | awk 'NR % 4 == 2' | tr NT TN | "
    else:
        raise ValueError("fmlrc2 correction requires NGS reads")

    runner.run(
        seq_cmd
        + "ropebwt2 -LR 2>{}/{}.ropebwt2.log | tr NT TN | "
        "fmlrc2-convert comp_msbwt.npy 1>{}/{}.fmlrc_convert.log 2>&1".format(
            q(log_dir), q(sample.sample_name), q(log_dir), q(sample.sample_name)
        ),
        cwd=work_dir,
    )

    raw_flnc = os.path.join(work_dir, "rawFlnc.fq")
    if sample.data_processed_location is None and sample.data_location:
        sample.data_processed_location = os.path.join(
            base_dir, "preprocess", (sample.tgs_plat or "tgs").lower(), "rawFlnc.fq"
        )
    if sample.data_processed_location:
        makeLink(sample.data_processed_location, raw_flnc)

    runner.run("seqkit seq --rna2dna rawFlnc.fq -w 0 > raw.dna.fastq", cwd=work_dir)
    runner.run(
        "fmlrc2 -t {} -k 25 59 -C 10 comp_msbwt.npy raw.dna.fastq "
        "fmlrc_corrected.fasta 1>{}/{}.fmlrc.log 2>&1".format(threads, q(log_dir), q(sample.sample_name)),
        cwd=work_dir,
    )
    sample.data_processed_location = os.path.join(work_dir, "fmlrc_corrected.fasta")
    return sample.data_processed_location


def retrieve_pacbio(
    sample: Sample,
    ccs_params: CcsParams,
    base_dir: str,
    *,
    threads: int,
    runner=None,
) -> str:
    """Retrieve PacBio FLNC reads (CCS/lima/isoseq refine)."""
    runner = _runner(runner, _PACBIO_ENV_MAP)
    work_dir = os.path.join(base_dir, "preprocess", "pacbio")
    os.makedirs(work_dir, exist_ok=True)
    strategy = (sample.strategy or "").lower()
    if strategy not in {"rsii", "sequel"}:
        raise ValueError("PacBio preproc requires strategy 'rsii' or 'sequel'")
    if not sample.primer:
        raise ValueError("PacBio preproc requires primer in samples.tsv")

    ccs_inputs = _collect_pacbio_inputs(sample.data_location) if strategy == "sequel" else []
    if ccs_inputs and all(is_ccs_bam(path) for path in ccs_inputs):
        # already-called CCS/HiFi reads (READTYPE=CCS): skip ccs; lima reads a plain BAM (no .pbi needed)
        ccs_bam = ccs_inputs[0] if len(ccs_inputs) == 1 else _merge_bams(ccs_inputs, os.path.join(work_dir, "CCS.bam"), runner)
    else:
        if strategy == "rsii":
            subreads = _retrieve_rsii_subreads(sample, work_dir, runner)
        else:
            subreads = _retrieve_sequel_subreads(sample, work_dir, runner)
        version = _ccs_version(runner)
        runner.run(_ccs_command(version, subreads, "CCS.bam", ccs_params, threads), cwd=work_dir)
        ccs_bam = "CCS.bam"

    runner.run(
        ["lima", ccs_bam, sample.primer, "fl.bam", "--isoseq", "--dump-clips", "-j", str(threads)],
        cwd=work_dir,
    )
    lima_bam = lima_output_bam(work_dir, sample.sample_name)     # BUG-101: named after the primer pair, not the sample
    runner.run(
        [
            "isoseq",
            "refine",
            lima_bam,
            sample.primer,
            "flnc.bam",
            "--require-polya",
            "-j",
            str(threads),
        ],
        cwd=work_dir,
    )
    raw_flnc = os.path.join(work_dir, "rawFlnc.fq")
    with open(raw_flnc, "w") as stdout:
        runner.run(["samtools", "fastq", "flnc.bam"], cwd=work_dir, stdout=stdout)
    sample.data_processed_location = raw_flnc
    return raw_flnc


def is_ccs_bam(path: str) -> bool:
    """True when a PacBio BAM declares READTYPE=CCS in its read-group description (ccs/HiFi output)."""
    if not path.endswith(".bam"):
        return False
    import pysam

    try:
        with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
            return any("READTYPE=CCS" in (rg.get("DS") or "") for rg in bam.header.to_dict().get("RG", []))
    except (OSError, ValueError):               # missing or not a BAM: let ccs report the problem
        return False


def _merge_bams(paths, out_bam, runner) -> str:
    runner.run(["samtools", "merge", "-f", out_bam, *paths])
    return out_bam


def lima_output_bam(work_dir: str, sample_name: str) -> str:
    """The demultiplexed BAM lima writes for ``fl.bam``: ``fl.<5p>--<3p>.bam`` named after the primer pair.
    Legacy assumed ``fl.5p--<sample>_3p.bam`` (BUG-101); with several primer pairs the one carrying the sample
    name wins, otherwise the choice must be unique."""
    candidates = sorted(os.path.basename(p) for p in glob.glob(os.path.join(work_dir, "fl.*--*.bam")))
    if not candidates:
        raise FileNotFoundError("lima produced no demultiplexed BAM (fl.<5p>--<3p>.bam) in {}".format(work_dir))
    if len(candidates) > 1:
        named = [c for c in candidates if sample_name in c]
        if len(named) != 1:
            raise ValueError("lima produced several demultiplexed BAMs {}; name the sample's 3' primer after the sample".format(candidates))
        return named[0]
    return candidates[0]


def _retrieve_rsii_subreads(sample: Sample, work_dir: str, runner) -> str:
    if not sample.data_location or not validateDir(sample.data_location):
        raise ValueError("RSII data_location must be a directory")
    fofn = "{}.fofn".format(sample.sample_name)
    runner.run("find {} -name '*.bax.h5' > {}".format(q(sample.data_location), q(fofn)), cwd=work_dir)
    runner.run(["bax2bam", "-f", fofn, "-o", sample.sample_name], cwd=work_dir)
    return "{}.subreads.bam".format(sample.sample_name)


def _retrieve_sequel_subreads(sample: Sample, work_dir: str, runner) -> str:
    files = _collect_pacbio_inputs(sample.data_location)
    for path in files:
        if path.endswith(".bam") and not os.path.exists(path + ".pbi"):
            runner.run(["pbindex", path])
    subreads = "{}.subreadset.xml".format(sample.sample_name)
    runner.run(
        [
            "dataset",
            "create",
            "--type",
            "SubreadSet",
            "--force",
            "--name",
            sample.project_name,
            subreads,
            *files,
        ],
        cwd=work_dir,
    )
    return subreads


def _collect_pacbio_inputs(data_location: str | None) -> list[str]:
    if not data_location:
        raise ValueError("PacBio preproc requires data_location")
    if validateDir(data_location):
        xml_files = sorted(str(p) for p in Path(data_location).rglob("*.subreadset.xml"))
        if xml_files:
            return xml_files
        # a Sequel run folder also holds *.scraps.bam (adapters, barcodes, low-quality regions), which are not subreads (BUG-185)
        bam_files = sorted(str(p) for p in Path(data_location).rglob("*.bam") if not p.name.endswith(".scraps.bam"))
        if bam_files:
            return bam_files
        raise ValueError("PacBio data directory must contain *.subreadset.xml or *.bam")
    files = [x.strip() for x in str(data_location).split(",") if x.strip()]
    if files and all(validateFile(x) for x in files):
        return files
    raise ValueError("PacBio data_location must be existing file(s) or directory")


def _ccs_version(runner) -> str:
    proc = runner.run(["ccs", "--version"], check=False)
    return ((getattr(proc, "stdout", "") or getattr(proc, "stderr", "")) or "").strip()


def _ccs_command(version: str, subreads: str, out_bam: str, params: CcsParams, threads: int) -> list[str]:
    if "ccs 3." in version:
        return [
            "ccs",
            subreads,
            out_bam,
            "-j",
            str(threads),
            "--noPolish",
            "--minPasses",
            str(params.min_pass),
            "--minLength",
            str(params.min_subread_length),
            "--minPredictedAccuracy",
            str(params.min_predicted_accuracy),
            "--minReadScore",
            str(params.min_read_score),
        ]
    # ccs >= 4 (BUG-102): --min-length is the minimum CCS draft length (config min_ccs_length), not a subread filter;
    # min_subread_length / min_read_score have no equivalent and are not passed.
    return [
        "ccs",
        subreads,
        out_bam,
        "-j",
        str(threads),
        "--min-passes",
        str(params.min_pass),
        "--min-length",
        str(params.min_ccs_length),
        "--min-rq",
        str(params.min_predicted_accuracy),
    ]


def polya_start(seq: str, *, min_tail: int = 10, min_a: float = 0.9, window: int = 400) -> int:
    """Start of the basecalled 3' poly(A) of a direct-RNA read: the longest 3' suffix (within the last ``window`` nt) of at
    least ``min_tail`` nt that is at least ``min_a`` A, moved to its first A; ``len(seq)`` when there is none."""
    import numpy as np

    tail = np.frombuffer(seq[-window:].upper().encode(), dtype=np.uint8)[::-1] == ord("A")
    k = np.arange(1, len(tail) + 1)
    ok = (k >= min_tail) & (np.cumsum(tail) >= min_a * k)
    if not ok.any():
        return len(seq)
    start = len(seq) - int(k[ok].max())
    while seq[start] not in "Aa":
        start += 1
    return start


def trim_polya_tails(fastq: str, out_fastq: str) -> int:
    """Copy of ``fastq`` without the basecalled 3' poly(A): direct-RNA reads keep part of their tail, which is
    soft-clipped and lowers read coverage below collapse's -c 0.9 (PacBio FLNC arrive tail-trimmed by isoseq).
    Reads that are poly(A) from end to end are left out (they would be empty records; BUG-186). Returns the number
    of reads trimmed."""
    trimmed = dropped = 0
    with open(fastq) as src, open(out_fastq, "w") as out:
        for head in src:
            seq, _, qual = src.readline().rstrip("\n"), src.readline(), src.readline().rstrip("\n")
            cut = polya_start(seq)
            trimmed += cut < len(seq)
            if cut == 0:
                dropped += 1
                continue
            out.write(f"{head}{seq[:cut]}\n+\n{qual[:cut]}\n")
    if dropped:
        get_logger().info("poly(A) trimming: %d read(s) were poly(A) only and are left out of %s", dropped, out_fastq)
    return trimmed


def retrieve_nanopore(
    sample: Sample,
    base_dir: str,
    *,
    threads: int,
    runner=None,
    trim_polya: bool = False,
) -> str:
    """Retrieve Nanopore reads, preserving legacy guppy/fastq concatenation path. With ``trim_polya`` the processed
    reads are a tail-trimmed copy; ``rawFlnc.fq`` keeps the basecalls for nanopolish."""
    runner = _runner(runner)
    if not sample.data_location:
        raise ValueError("Nanopore preproc requires data_location")
    work_dir = os.path.join(base_dir, "preprocess", "nanopore")
    os.makedirs(work_dir, exist_ok=True)
    raw_flnc = os.path.join(work_dir, "rawFlnc.fq")
    file_type = check_fast5_files(sample.data_location)
    if file_type == "bam":                       # dorado basecalls (F-11)
        runner.run("samtools fastq {} > rawFlnc.fq".format(q(sample.data_location)), cwd=work_dir)
        if dorado_polya_table(sample.data_location, os.path.join(work_dir, DORADO_POLYA)) and not sample.polya_location:
            sample.polya_location = os.path.join(work_dir, DORADO_POLYA)
    elif file_type == "fast5":
        # ONT's Guppy ships the guppy_basecaller binary (legacy called guppy_basecaller.py; BUG-187)
        if shutil.which("guppy_basecaller") is None:
            raise ValueError("a fast5 directory as data_location needs ONT's guppy_basecaller (GPU; Guppy is "
                             "discontinued). Basecall elsewhere with dorado or guppy and give the FASTQ, or dorado's "
                             "unaligned BAM, as data_location; keep the fast5 directory as fast5_dir for nanopolish")
        basecalling_out = "guppy_basecalling"
        runner.run(
            "guppy_basecaller --flowcell {} --kit {} -i {} --save_path {} "
            "--cpu_threads_per_caller {} --number_callers {}".format(
                sample.flowcell_type or "FLO-MIN106",
                sample.kit_type or "SQK-RNA002",
                q(sample.data_location),
                basecalling_out,
                threads,
                1,
            ),
            cwd=work_dir,
        )
        runner.run("cat {}/*.fastq > rawFlnc.fq".format(basecalling_out), cwd=work_dir)
    elif validateFile(sample.data_location):
        runner.run("zcat -f {} > rawFlnc.fq".format(q(sample.data_location)), cwd=work_dir)
    else:
        runner.run("zcat -f {}/* > rawFlnc.fq".format(q(sample.data_location)), cwd=work_dir)
    sample.data_processed_location = raw_flnc
    if trim_polya:
        sample.data_processed_location = os.path.join(work_dir, "flnc.polyA_trimmed.fq")
        trim_polya_tails(raw_flnc, sample.data_processed_location)
    return sample.data_processed_location


DORADO_POLYA = "dorado_polya.tsv"


def dorado_polya_table(bam: str, out_tsv: str) -> int:
    """dorado's poly(A) tail estimates (the ``pt:i`` tag of ``--estimate-poly-a``) as a table palen_as reads like
    nanopolish's ``polya_results.tsv`` (PASS when a tail was found). Returns the number of reads carrying the tag; with
    none (dorado ran without ``--estimate-poly-a``) no table is written."""
    import pysam

    rows = []
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as f:
        for r in f:
            if not (r.is_secondary or r.is_supplementary) and r.has_tag("pt"):
                pt = int(r.get_tag("pt"))
                rows.append(f"{r.query_name}\t.\t{pt}\t{'PASS' if pt > 0 else 'NOREGION'}\n")
    if rows:
        with open(out_tsv, "w") as out:
            out.write("readname\tcontig\tpolya_length\tqc_tag\n" + "".join(rows))
    return len(rows)


def check_fast5_files(path: str) -> str:
    if validateFile(path):
        if path.endswith(("fq", "fastq", "fq.gz", "fastq.gz")):
            return "fastq"
        if path.endswith(".bam"):
            return "bam"
        raise ValueError("Nanopore data file must be fastq/fq (optionally .gz) or a dorado BAM")
    if not validateDir(path):
        raise ValueError("Nanopore data_location must be a file or directory")
    entries = [x for x in os.listdir(path) if not x.endswith(".txt")]
    if not entries:
        raise ValueError("Nanopore data directory is empty")
    fast5 = [x for x in entries if x.endswith("fast5")]
    fastq = [x for x in entries if x.endswith(("fq", "fastq", "fq.gz", "fastq.gz"))]
    if len(fast5) == len(entries):
        return "fast5"
    if len(fastq) == len(entries):
        return "fastq"
    raise ValueError("Nanopore data directory must contain only fast5 or only fastq/fq files")


def _looks_like_fasta_or_fastq(path: str) -> bool:
    opener = gzip.open if path.endswith(".gz") else open
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    return line[0] in {">", "@"}
    except OSError:
        return False
    return False


def find_preprocessed_flnc(base_dir: str) -> str | None:
    candidates = [
        os.path.join(base_dir, "preprocess", "fmlrc", "fmlrc_corrected.fasta"),
        os.path.join(base_dir, "preprocess", "pacbio", "rawFlnc.fq"),
        os.path.join(base_dir, "preprocess", "nanopore", "flnc.polyA_trimmed.fq"),
        os.path.join(base_dir, "preprocess", "nanopore", "rawFlnc.fq"),
    ]
    for path in candidates:
        if validateFile(path):
            return path
    return None


def preprocess_sample(
    sample: Sample,
    ccs_params: CcsParams,
    base_dir: str,
    *,
    threads: int,
    runner=None,
    trim_polya: bool = False,
) -> dict:
    """Run migrated preproc for one sample and return updated sample paths."""
    print(
        getCurrentTime()
        + " Preprocess raw reads for project {} sample {}...".format(
            sample.project_name, sample.sample_name
        )
    )
    os.makedirs(base_dir, exist_ok=True)

    if sample.data_processed_location and validateFile(sample.data_processed_location):
        if not _looks_like_fasta_or_fastq(sample.data_processed_location):
            raise ValueError("data_processed_location must be FASTA/FASTQ")
        if _has_ngs(sample):
            process_rnaseq(sample, base_dir, threads=threads, runner=runner)
            if sample.use_fmlrc2:
                correct_with_fmlrc2(sample, base_dir, threads=threads, runner=runner)
    elif sample.data_location and (validateFile(sample.data_location) or validateDir(sample.data_location)):
        platform = (sample.tgs_plat or "").lower()
        ngs_done = False
        if platform == "pacbio":
            retrieve_pacbio(sample, ccs_params, base_dir, threads=threads, runner=runner)
        elif platform == "nanopore":
            retrieve_nanopore(sample, base_dir, threads=threads, runner=runner, trim_polya=trim_polya)
        elif _has_ngs(sample):
            process_rnaseq(sample, base_dir, threads=threads, runner=runner)
            ngs_done = True
        else:
            raise ValueError("preproc requires tgs_plat pacbio/nanopore or NGS reads")
        if _has_ngs(sample) and not ngs_done:
            process_rnaseq(sample, base_dir, threads=threads, runner=runner)
            if sample.data_processed_location and sample.use_fmlrc2:
                correct_with_fmlrc2(sample, base_dir, threads=threads, runner=runner)
    elif _has_ngs(sample):
        process_rnaseq(sample, base_dir, threads=threads, runner=runner)
    else:
        raise ValueError("preproc requires valid data_location, data_processed_location, or NGS reads")

    print(
        getCurrentTime()
        + " Preprocess raw reads for project {} sample {} done!".format(
            sample.project_name, sample.sample_name
        )
    )
    return {
        "flnc": sample.data_processed_location,
        "ngs_left_reads": sample.ngs_left_reads,
        "ngs_right_reads": sample.ngs_right_reads,
        "polya_location": sample.polya_location,
    }
