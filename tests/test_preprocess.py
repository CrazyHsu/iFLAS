"""Tests for the migrated preproc stage command contracts."""

from __future__ import annotations

import os

from types import SimpleNamespace

import pytest

from iflas.config import CcsParams, Sample
from iflas.stages import preprocess
from iflas.stages.preprocess import (
    _ccs_command,
    _collect_pacbio_inputs,
    correct_with_fmlrc2,
    is_ccs_bam,
    lima_output_bam,
    process_rnaseq,
    retrieve_nanopore,
    retrieve_pacbio,
    trim_polya_tails,
)


class _FakeRunner:
    def __init__(self, ccs_version="ccs 6.4.0"):
        self.calls = []
        self.ccs_version = ccs_version

    def run(self, cmd, **kwargs):
        self.calls.append((cmd, kwargs))
        if cmd == ["ccs", "--version"]:
            return SimpleNamespace(returncode=0, stdout=self.ccs_version + "\n", stderr="")
        if cmd[0] == "lima":                                      # lima names its output after the primer pair
            open(os.path.join(kwargs["cwd"], "fl.5p--3p.bam"), "w").close()
        return SimpleNamespace(returncode=0, stdout="", stderr="")


def _ngs_sample(left, right, **kw):
    return Sample(sample_name=kw.pop("sample_name", "s1"), project_name="demo", ref_strain="B73", ngs_left_reads=str(left),
                  ngs_right_reads=str(right), ngs_reads_paired="paired", ngs_reads_length=151, **kw)


def test_process_rnaseq_paired_fastp_contract(tmp_path):
    left, right = tmp_path / "reads_R1.fastq.gz", tmp_path / "reads_R2.fastq.gz"
    left.write_text("x")
    right.write_text("x")
    sample = _ngs_sample(left, right)
    runner = _FakeRunner()
    base = tmp_path / "out" / "demo" / "s1"

    out = process_rnaseq(sample, str(base), threads=3, runner=runner)

    # BUG-152: products under the sample dir (legacy wrote beside the raw reads), with a .source sidecar each
    expected_left = str(base / "preprocess" / "rna-seq" / "fastp.reads_R1.fastq.gz")
    expected_right = str(base / "preprocess" / "rna-seq" / "fastp.reads_R2.fastq.gz")
    assert out == {"ngs_left_reads": expected_left, "ngs_right_reads": expected_right}
    assert sample.ngs_left_reads == expected_left
    assert runner.calls == [(
        f"fastp -i {left} -I {right} -o {expected_left} -O {expected_right} "
        f"-w 3 -q 20 -l 150 -j {base / 'log'}/fastp.reads_R1.fastq.gz.json -h {base / 'log'}/fastp.reads_R1.fastq.gz.html "
        f"2> {base / 'log'}/fastp.reads_R1.fastq.gz.log", {})]
    assert (base / "preprocess" / "rna-seq" / "fastp.reads_R1.fastq.gz.source").read_text().strip() == str(left.resolve())


def test_process_rnaseq_passes_through_already_filtered_reads(tmp_path):
    sample = _ngs_sample(tmp_path / "fastp.reads_R1.fastq.gz", tmp_path / "fastp.reads_R2.fastq.gz")
    runner = _FakeRunner()
    out = process_rnaseq(sample, str(tmp_path / "out" / "demo" / "s1"), threads=3, runner=runner)
    assert out["ngs_left_reads"] == str(tmp_path / "fastp.reads_R1.fastq.gz") and runner.calls == []


def test_process_rnaseq_reuses_products_of_the_same_library(tmp_path):
    """A library already filtered beside the input (legacy layout) or by another sample of the project (same
    .source) is not filtered again; a sample whose inputs share a basename gets parent-directory-qualified names."""
    raw = tmp_path / "raw"
    raw.mkdir()
    left, right = raw / "reads_R1.fastq.gz", raw / "reads_R2.fastq.gz"
    left.write_text("x")
    right.write_text("x")
    (raw / "fastp.reads_R1.fastq.gz").write_text("filtered")                      # legacy copy beside the input
    sibling = tmp_path / "out" / "demo" / "s0" / "preprocess" / "rna-seq"
    sibling.mkdir(parents=True)
    (sibling / "fastp.reads_R2.fastq.gz").write_text("filtered")                  # another sample's product
    (sibling / "fastp.reads_R2.fastq.gz.source").write_text(str(right.resolve()) + "\n")
    sample = _ngs_sample(left, right)
    runner = _FakeRunner()

    out = process_rnaseq(sample, str(tmp_path / "out" / "demo" / "s1"), threads=2, runner=runner)

    assert runner.calls == []
    assert out == {"ngs_left_reads": str(raw / "fastp.reads_R1.fastq.gz"), "ngs_right_reads": str(sibling / "fastp.reads_R2.fastq.gz")}

    a, b = tmp_path / "lib1", tmp_path / "lib2"
    for d in (a, b):
        d.mkdir()
        (d / "R1.fq.gz").write_text("x")
        (d / "R2.fq.gz").write_text("x")
    sample = _ngs_sample(f"{a / 'R1.fq.gz'};{b / 'R1.fq.gz'}", f"{a / 'R2.fq.gz'};{b / 'R2.fq.gz'}", sample_name="s2")
    runner = _FakeRunner()
    out = process_rnaseq(sample, str(tmp_path / "out" / "demo" / "s2"), threads=2, runner=runner)
    rna = tmp_path / "out" / "demo" / "s2" / "preprocess" / "rna-seq"
    assert out["ngs_left_reads"] == f"{rna / 'fastp.lib1.R1.fq.gz'};{rna / 'fastp.lib2.R1.fq.gz'}"
    assert len(runner.calls) == 2


def test_correct_with_fmlrc2_contract(tmp_path):
    flnc = tmp_path / "flnc.fq"
    left = tmp_path / "fastp.R1.fq.gz"
    right = tmp_path / "fastp.R2.fq.gz"
    flnc.write_text("@r\nAC\n+\n!!\n")
    left.write_text("x")
    right.write_text("x")
    sample = Sample(
        sample_name="s1",
        project_name="demo",
        ref_strain="B73",
        tgs_plat="pacbio",
        data_processed_location=str(flnc),
        ngs_left_reads=str(left),
        ngs_right_reads=str(right),
    )
    runner = _FakeRunner()

    out = correct_with_fmlrc2(sample, str(tmp_path / "out" / "demo" / "s1"), threads=4, runner=runner)

    work_dir = tmp_path / "out" / "demo" / "s1" / "preprocess" / "fmlrc"
    log_dir = tmp_path / "out" / "demo" / "s1" / "log"
    assert out == str(work_dir / "fmlrc_corrected.fasta")
    assert sample.data_processed_location == out
    assert runner.calls == [
        (f"zcat -f {left} > leftReads.fastq", {"cwd": str(work_dir)}),       # -f: plain FASTQ passes through
        (f"zcat -f {right} > rightReads.fastq", {"cwd": str(work_dir)}),
        (
            "cat leftReads.fastq rightReads.fastq | awk 'NR % 4 == 2' | tr NT TN | "
            f"ropebwt2 -LR 2>{log_dir}/s1.ropebwt2.log | tr NT TN | "
            f"fmlrc2-convert comp_msbwt.npy 1>{log_dir}/s1.fmlrc_convert.log 2>&1",
            {"cwd": str(work_dir)},
        ),
        ("seqkit seq --rna2dna rawFlnc.fq -w 0 > raw.dna.fastq", {"cwd": str(work_dir)}),
        (
            f"fmlrc2 -t 4 -k 25 59 -C 10 comp_msbwt.npy raw.dna.fastq "
            f"fmlrc_corrected.fasta 1>{log_dir}/s1.fmlrc.log 2>&1",
            {"cwd": str(work_dir)},
        ),
    ]


def test_retrieve_pacbio_sequel_current_tool_contract(tmp_path):
    raw = tmp_path / "s1.subreads.bam"
    primer = tmp_path / "primer.fa"
    raw.write_text("bam")
    primer.write_text(">p\nAC\n")
    sample = Sample(
        sample_name="s1",
        project_name="demo",
        ref_strain="B73",
        tgs_plat="pacbio",
        strategy="sequel",
        data_location=str(raw),
        primer=str(primer),
    )
    runner = _FakeRunner()

    out = retrieve_pacbio(
        sample,
        CcsParams(min_pass=3, min_subread_length=75, min_predicted_accuracy=0.91, min_read_score=0.8),
        str(tmp_path / "out" / "demo" / "s1"),
        threads=5,
        runner=runner,
    )

    work_dir = tmp_path / "out" / "demo" / "s1" / "preprocess" / "pacbio"
    assert out == str(work_dir / "rawFlnc.fq")
    assert sample.data_processed_location == out
    assert runner.calls[:6] == [
        (["pbindex", str(raw)], {}),
        (
            [
                "dataset",
                "create",
                "--type",
                "SubreadSet",
                "--force",
                "--name",
                "demo",
                "s1.subreadset.xml",
                str(raw),
            ],
            {"cwd": str(work_dir)},
        ),
        (["ccs", "--version"], {"check": False}),
        (
            [
                "ccs",
                "s1.subreadset.xml",
                "CCS.bam",
                "-j",
                "5",
                "--min-passes",
                "3",
                "--min-length",
                "50",                     # ccs>=4: draft length = min_ccs_length default, not min_subread_length (BUG-102)
                "--min-rq",
                "0.91",
            ],
            {"cwd": str(work_dir)},
        ),
        (
            ["lima", "CCS.bam", str(primer), "fl.bam", "--isoseq", "--dump-clips", "-j", "5"],
            {"cwd": str(work_dir)},
        ),
        (
            [
                "isoseq",
                "refine",
                "fl.5p--3p.bam",
                str(primer),
                "flnc.bam",
                "--require-polya",
                "-j",
                "5",
            ],
            {"cwd": str(work_dir)},
        ),
    ]
    assert runner.calls[6][0] == ["samtools", "fastq", "flnc.bam"]
    assert runner.calls[6][1]["cwd"] == str(work_dir)


def test_ccs_command_maps_config_to_the_installed_ccs_generation():
    params = CcsParams(min_pass=2, min_subread_length=50, min_ccs_length=200, min_predicted_accuracy=0.9, min_read_score=0.75)
    # ccs >= 4: --min-length is the draft (CCS) length; subread length / read score have no equivalent (BUG-102)
    assert _ccs_command("ccs 6.4.0 (commit v6.4.0)", "in.xml", "CCS.bam", params, 4) == [
        "ccs", "in.xml", "CCS.bam", "-j", "4", "--min-passes", "2", "--min-length", "200", "--min-rq", "0.9"]
    assert _ccs_command("ccs 3.4.1", "in.xml", "CCS.bam", params, 4)[5:] == [
        "--noPolish", "--minPasses", "2", "--minLength", "50", "--minPredictedAccuracy", "0.9", "--minReadScore", "0.75"]


def test_lima_output_bam_is_named_after_the_primer_pair(tmp_path):
    (tmp_path / "fl.5p--3p.bam").write_text("")
    (tmp_path / "fl.lima.summary").write_text("")
    assert lima_output_bam(str(tmp_path), "s1") == "fl.5p--3p.bam"           # BUG-101: not fl.5p--s1_3p.bam
    (tmp_path / "fl.5p--s1_3p.bam").write_text("")
    assert lima_output_bam(str(tmp_path), "s1") == "fl.5p--s1_3p.bam"        # several pairs: the sample's own wins
    import pytest
    with pytest.raises(ValueError):
        lima_output_bam(str(tmp_path), "s2")
    with pytest.raises(FileNotFoundError):
        lima_output_bam(str(tmp_path / "empty"), "s1")


def test_is_ccs_bam_reads_the_read_group_description(tmp_path):
    import pysam
    for name, ds in (("ccs.bam", "READTYPE=CCS"), ("subreads.bam", "READTYPE=SUBREAD")):
        with pysam.AlignmentFile(str(tmp_path / name), "wb", header={"HD": {"VN": "1.6"}, "RG": [{"ID": "x", "DS": ds}]}):
            pass
    assert is_ccs_bam(str(tmp_path / "ccs.bam")) is True
    assert is_ccs_bam(str(tmp_path / "subreads.bam")) is False
    assert is_ccs_bam(str(tmp_path / "reads.fq")) is False



def test_fastp_reports_are_kept_per_repeat(tmp_path):
    """Each filtered file gets its own fastp QC report and message log (one shared fastp.json/html was overwritten
    by every repeat: BUG-188)."""
    files = {name: tmp_path / name for name in ("a_R1.fq.gz", "a_R2.fq.gz", "b_R1.fq.gz", "b_R2.fq.gz")}
    for f in files.values():
        f.write_text("x")
    sample = _ngs_sample(f"{files['a_R1.fq.gz']};{files['b_R1.fq.gz']}", f"{files['a_R2.fq.gz']};{files['b_R2.fq.gz']}")
    runner = _FakeRunner()
    process_rnaseq(sample, str(tmp_path / "out" / "demo" / "s1"), threads=1, runner=runner)
    reports = [c[0].split(" -j ")[1].split()[0] for c in runner.calls]
    assert [os.path.basename(r) for r in reports] == ["fastp.a_R1.fq.gz.json", "fastp.b_R1.fq.gz.json"]


def test_pacbio_directory_without_xml_leaves_out_scraps(tmp_path):
    """A Sequel run folder copied without its XML holds *.subreads.bam and *.scraps.bam; scraps are not subreads
    (legacy collected every *.bam: BUG-185)."""
    for name in ("m1.subreads.bam", "m1.scraps.bam", "m2.subreads.bam"):
        (tmp_path / name).write_text("x")
    assert [os.path.basename(p) for p in _collect_pacbio_inputs(str(tmp_path))] == ["m1.subreads.bam", "m2.subreads.bam"]


def test_trim_polya_tails_leaves_out_reads_that_are_only_polya(tmp_path):
    """Direct-RNA runs yield reads that are poly(A) end to end; trimming would make them empty records (BUG-186)."""
    fq, out = tmp_path / "in.fq", tmp_path / "out.fq"
    fq.write_text("@r1\nACGUACGUACGU" + "A" * 20 + "\n+\n" + "I" * 32 + "\n@r2\n" + "A" * 25 + "\n+\n" + "I" * 25 + "\n")
    assert trim_polya_tails(str(fq), str(out)) == 2
    assert out.read_text().splitlines() == ["@r1", "ACGUACGUACGU", "+", "I" * 12]


def test_fast5_input_needs_guppy_basecaller(tmp_path, monkeypatch):
    """ONT's Guppy ships guppy_basecaller (legacy called guppy_basecaller.py, which does not exist: BUG-187); without it
    the error says to basecall first."""
    fast5 = tmp_path / "fast5"
    fast5.mkdir()
    (fast5 / "a.fast5").write_text("x")
    sample = Sample(sample_name="s1", project_name="demo", ref_strain="B73", tgs_plat="nanopore", data_location=str(fast5))
    monkeypatch.setattr(preprocess.shutil, "which", lambda name: None)
    with pytest.raises(ValueError, match="dorado"):
        retrieve_nanopore(sample, str(tmp_path / "s1"), threads=1, runner=_FakeRunner())
    monkeypatch.setattr(preprocess.shutil, "which", lambda name: "/opt/ont-guppy/bin/" + name)
    runner = _FakeRunner()
    retrieve_nanopore(sample, str(tmp_path / "s1"), threads=1, runner=runner, trim_polya=False)
    assert runner.calls[0][0].startswith("guppy_basecaller --flowcell ")
