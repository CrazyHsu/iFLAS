"""Tests for the modern pipeline orchestrator (BUG-105 / iflas all)."""

from __future__ import annotations

import os

from iflas.config import (
    CollapseParams,
    IflasConfig,
    MergedTgsSample,
    Minimap2Params,
    Options,
    PuFilterParams,
    RefSection,
    Sample,
)
from iflas.orchestrator import (
    DEFAULT_STAGES,
    SampleJob,
    _stage_mapping,
    _stage_pu_filter,
    parse_ngs_repeats,
    run_pipeline,
    sample_base_dir,
)


class _Future:
    def __init__(self, value):
        self._value = value

    def result(self):
        return self._value


class _ImmediateExecutor:
    instances = []

    def __init__(self, max_workers=None):
        self.max_workers = max_workers
        self.submitted = []
        self.__class__.instances.append(self)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def submit(self, fn, *args, **kwargs):
        self.submitted.append((fn, args, kwargs))
        return _Future(fn(*args, **kwargs))


def _config(tmp_path):
    return IflasConfig(
        reference=[
            RefSection(
                ref_strain="B73",
                ref_genome="/ref/genome.fa",
                ref_gpe="/ref/genes.gpe",
                ref_bed="/ref/genes.bed12+",
                ref_mm2_index="/ref/genome.mm2",
                hisat2_index="/ref/hisat2/B73",
                ref_gtf="/ref/genes.gtf",
            )
        ],
        minimap2=Minimap2Params(max_intron_length=70000),
        collapse=CollapseParams(fl_coverage=3),
        pu_filter=PuFilterParams(filter_score=0.25),
        options=Options(threads=8),
        samples=[
            Sample(
                sample_name="s1",
                project_name="demo",
                ref_strain="B73",
                tgs_plat="pacbio",
                data_processed_location="/reads/s1.flnc.fa",
                ngs_left_reads="/reads/s1.R1.fq;/reads/s1b.R1.fq",
                ngs_right_reads="/reads/s1.R2.fq;/reads/s1b.R2.fq",
            ),
            Sample(
                sample_name="s2",
                project_name="demo",
                ref_strain="B73",
                tgs_plat="nanopore",
                data_processed_location="/reads/s2.flnc.fa",
            ),
        ],
        dirs={"out_dir": str(tmp_path / "out")},
    )


def test_sample_base_dir_uses_legacy_project_sample_layout(tmp_path):
    cfg = _config(tmp_path)
    assert sample_base_dir(cfg, cfg.samples[0]) == os.path.abspath(
        tmp_path / "out" / "demo" / "s1"
    )


def test_parse_ngs_repeats_supports_legacy_semicolon_replicates(tmp_path):
    sample = _config(tmp_path).samples[0]
    assert parse_ngs_repeats(sample) == [
        {"left": "/reads/s1.R1.fq", "right": "/reads/s1.R2.fq"},
        {"left": "/reads/s1b.R1.fq", "right": "/reads/s1b.R2.fq"},
    ]


def test_parse_ngs_repeats_rejects_single_end_with_both_read_sides(tmp_path):
    sample = Sample(
        sample_name="s1",
        project_name="demo",
        ref_strain="B73",
        ngs_reads_paired="single",
        ngs_left_reads="/reads/s1.R1.fq",
        ngs_right_reads="/reads/s1.R2.fq",
    )

    try:
        parse_ngs_repeats(sample)
    except ValueError as exc:
        assert "single-end NGS reads require exactly one read side" in str(exc)
    else:
        raise AssertionError("single-end NGS with both left/right reads should be rejected")


def test_run_pipeline_dispatches_stage_barriers_with_process_executor(tmp_path):
    cfg = _config(tmp_path)
    calls = []

    def record(job: SampleJob):
        calls.append((job.stage, job.sample.sample_name, job.threads, job.base_dir))
        return {"stage": job.stage, "sample": job.sample.sample_name}

    handlers = {"mapping": record, "collapse": record}
    results = run_pipeline(
        cfg,
        stages=("mapping", "collapse"),
        max_workers=2,
        handlers=handlers,
        executor_factory=_ImmediateExecutor,
    )

    assert [r.stage for r in results] == ["mapping", "mapping", "collapse", "collapse"]
    assert calls == [
        ("mapping", "s1", 4, os.path.abspath(tmp_path / "out" / "demo" / "s1")),
        ("mapping", "s2", 4, os.path.abspath(tmp_path / "out" / "demo" / "s2")),
        ("collapse", "s1", 4, os.path.abspath(tmp_path / "out" / "demo" / "s1")),
        ("collapse", "s2", 4, os.path.abspath(tmp_path / "out" / "demo" / "s2")),
    ]
    assert [ex.max_workers for ex in _ImmediateExecutor.instances[-2:]] == [2, 2]


def test_run_pipeline_preproc_outputs_feed_following_stage(tmp_path):
    cfg = IflasConfig(
        reference=[
            RefSection(
                ref_strain="B73",
                ref_genome="/ref/genome.fa",
                ref_gtf="/ref/genes.gtf",
            )
        ],
        options=Options(threads=4),
        samples=[
            Sample(
                sample_name="s1",
                project_name="demo",
                ref_strain="B73",
                data_location="/raw/s1.subreads.bam",
                ngs_left_reads="/reads/s1.R1.fq",
                ngs_right_reads="/reads/s1.R2.fq",
            )
        ],
        dirs={"out_dir": str(tmp_path / "out")},
    )
    seen = {}

    def fake_preproc(job: SampleJob):
        return {
            "flnc": "/processed/s1.flnc.fa",
            "ngs_left_reads": "/processed/fastp.R1.fq",
            "ngs_right_reads": "/processed/fastp.R2.fq",
        }

    def fake_mapping(job: SampleJob):
        seen["flnc"] = job.sample.data_processed_location
        seen["left"] = job.sample.ngs_left_reads
        seen["right"] = job.sample.ngs_right_reads
        return {}

    run_pipeline(
        cfg,
        stages=("preproc", "mapping"),
        max_workers=1,
        handlers={"preproc": fake_preproc, "mapping": fake_mapping},
        executor_factory=_ImmediateExecutor,
    )

    assert seen == {
        "flnc": "/processed/s1.flnc.fa",
        "left": "/processed/fastp.R1.fq",
        "right": "/processed/fastp.R2.fq",
    }


def test_run_pipeline_can_merge_samples_by_strain_for_migrated_stages(tmp_path):
    cfg = IflasConfig(
        reference=[
            RefSection(
                ref_strain="B73",
                ref_genome="/ref/genome.fa",
                ref_gtf="/ref/genes.gtf",
            )
        ],
        options=Options(threads=8, merge_data_from_same_strain=True),
        samples=[
            Sample(
                sample_name="s1",
                project_name="demo",
                ref_strain="B73",
                strain="B73",
                tgs_plat="pacbio",
                data_processed_location="/processed/s1.flnc.fa",
            ),
            Sample(
                sample_name="s2",
                project_name="demo",
                ref_strain="B73",
                strain="B73",
                tgs_plat="pacbio",
                data_processed_location="/processed/s2.flnc.fa",
            ),
        ],
        dirs={"out_dir": str(tmp_path / "out")},
    )
    calls = []

    def record(job: SampleJob):
        calls.append(
            (
                job.stage,
                job.sample.sample_name,
                job.sample.data_processed_location,
                job.threads,
                job.base_dir,
            )
        )
        return {}

    results = run_pipeline(
        cfg,
        stages=("mapping", "collapse"),
        max_workers=2,
        merge_samples=True,
        handlers={"mapping": record, "collapse": record},
        executor_factory=_ImmediateExecutor,
    )

    assert [r.sample_name for r in results] == ["B73_all", "B73_all"]
    assert calls == [
        (
            "mapping",
            "B73_all",
            ["/processed/s1.flnc.fa", "/processed/s2.flnc.fa"],
            8,
            os.path.abspath(tmp_path / "out" / "demo" / "B73_all"),
        ),
        (
            "collapse",
            "B73_all",
            ["/processed/s1.flnc.fa", "/processed/s2.flnc.fa"],
            8,
            os.path.abspath(tmp_path / "out" / "demo" / "B73_all"),
        ),
    ]


def test_run_pipeline_preproc_stays_unmerged_then_mapping_uses_merged_outputs(tmp_path):
    cfg = IflasConfig(
        reference=[
            RefSection(
                ref_strain="B73",
                ref_genome="/ref/genome.fa",
                ref_gtf="/ref/genes.gtf",
            )
        ],
        options=Options(threads=6, merge_data_from_same_strain=True),
        samples=[
            Sample(
                sample_name="s1",
                project_name="demo",
                ref_strain="B73",
                strain="B73",
                data_location="/raw/s1.bam",
            ),
            Sample(
                sample_name="s2",
                project_name="demo",
                ref_strain="B73",
                strain="B73",
                data_location="/raw/s2.bam",
            ),
        ],
        dirs={"out_dir": str(tmp_path / "out")},
    )
    seen = []

    def fake_preproc(job: SampleJob):
        seen.append((job.stage, job.sample.sample_name, job.threads))
        return {"flnc": f"/processed/{job.sample.sample_name}.flnc.fa"}

    def fake_mapping(job: SampleJob):
        seen.append(
            (
                job.stage,
                job.sample.sample_name,
                job.sample.data_processed_location,
                job.threads,
            )
        )
        return {}

    results = run_pipeline(
        cfg,
        stages=("preproc", "mapping"),
        max_workers=2,
        merge_samples=True,
        handlers={"preproc": fake_preproc, "mapping": fake_mapping},
        executor_factory=_ImmediateExecutor,
    )

    assert [r.sample_name for r in results] == ["s1", "s2", "B73_all"]
    assert seen == [
        ("preproc", "s1", 3),
        ("preproc", "s2", 3),
        ("mapping", "B73_all", ["/processed/s1.flnc.fa", "/processed/s2.flnc.fa"], 6),
    ]


def test_stage_mapping_concatenates_merged_flnc_inputs_like_legacy_mapping(tmp_path, monkeypatch):
    flnc1 = tmp_path / "s1.fq"
    flnc2 = tmp_path / "s2.fq"
    flnc1.write_text("@r1\nAC\n+\nII\n")
    flnc2.write_text("@r2\nGT\n+\nII\n")

    cfg = IflasConfig(
        reference=[
            RefSection(
                ref_strain="B73",
                ref_genome="/ref/genome.fa",
                ref_gtf="/ref/genes.gtf",
                ref_mm2_index="/ref/genome.mm2",
                ref_bed="/ref/genes.bed12+",
            )
        ],
        dirs={"out_dir": str(tmp_path / "out")},
    )
    sample = MergedTgsSample(
        project_name="demo",
        sample_name="B73_all",
        ref_strain="B73",
        strain="B73",
        tgs_plat="pacbio",
        data_processed_location=[str(flnc1), str(flnc2)],
    )
    ref = cfg.ref_for("B73")
    job = SampleJob(
        stage="mapping",
        config=cfg,
        sample=sample,
        ref=ref,
        base_dir=sample_base_dir(cfg, sample),
        threads=4,
    )
    seen = {}

    def fake_minimap2(target, flnc, out_sam, **kwargs):
        seen["target"] = target
        seen["flnc"] = flnc
        seen["out_sam"] = out_sam
        return out_sam

    def fake_filter(sam_file=None, out_prefix="flnc", work_dir=None, **kwargs):
        seen["filter"] = dict(kwargs, sam_file=sam_file, out_prefix=out_prefix)
        return {"bam": os.path.join(work_dir, f"{out_prefix}.mm2.sorted.bam"),
                "bed12plus": os.path.join(work_dir, f"{out_prefix}.addCVandID.bed12+")}

    import iflas.stages.align

    monkeypatch.setattr(iflas.stages.align, "minimap2_tgs", fake_minimap2)
    monkeypatch.setattr(iflas.stages.align, "mapping_filter_and_add_tags", fake_filter)
    monkeypatch.setattr(iflas.stages.align, "flnc_to_fasta", lambda flnc, out_fa, **kw: out_fa)

    out = _stage_mapping(job)

    merged_flnc = tmp_path / "out" / "demo" / "B73_all" / "preprocess" / "pacbio" / "rawFlnc.fq"
    assert merged_flnc.read_text() == "@r1\nAC\n+\nII\n@r2\nGT\n+\nII\n"
    assert seen["flnc"] == str(merged_flnc)
    assert seen["target"] == "/ref/genome.mm2"
    assert seen["out_sam"] == str(tmp_path / "out" / "demo" / "B73_all" / "mapping" / "flnc.mm2.sam")
    # BUG-143: the legacy junction-combination filter runs on the SAM with the reference bed (no NGS here)
    assert seen["filter"]["sam_file"] == seen["out_sam"]
    assert seen["filter"]["ref_bed_file"] == "/ref/genes.bed12+"
    assert seen["filter"]["junc_bed_file"] is None
    assert seen["filter"]["junc_comb_sup"] == 2
    assert out["flnc"] == str(merged_flnc)
    assert out["tgs_bam"] == str(tmp_path / "out" / "demo" / "B73_all" / "mapping" / "flnc.mm2.sorted.bam")


def test_stage_pu_filter_passes_configured_feature_file(tmp_path, monkeypatch):
    cfg = IflasConfig(
        reference=[
            RefSection(
                ref_strain="B73",
                ref_genome="/ref/genome.fa",
                ref_bed="/ref/genes.bed12+",
                ref_gtf="/ref/genes.gtf",
            )
        ],
        pu_filter=PuFilterParams(filter_score=0.25, feature_file="/features/isoFeatures.txt"),
        dirs={"out_dir": str(tmp_path / "out")},
    )
    sample = Sample(sample_name="s1", project_name="demo", ref_strain="B73")
    job = SampleJob(
        stage="pu_filter",
        config=cfg,
        sample=sample,
        ref=cfg.ref_for("B73"),
        base_dir=sample_base_dir(cfg, sample),
        threads=4,
    )
    seen = {}

    def fake_filter_hq_isoforms(base_dir, genome_fa, ref_bed, **kwargs):
        seen.update({"base_dir": base_dir, "genome_fa": genome_fa, "ref_bed": ref_bed, **kwargs})
        return {"hq_bed": "/out/hq.collapsed.bed12+"}

    import iflas.stages.pu_filter

    monkeypatch.setattr(iflas.stages.pu_filter, "filter_hq_isoforms", fake_filter_hq_isoforms)

    assert _stage_pu_filter(job) == {"hq_bed": "/out/hq.collapsed.bed12+"}
    assert seen["feature_file"] == "/features/isoFeatures.txt"
    assert seen["filter_score"] == 0.25
    assert seen["bagging_n_jobs"] == 4          # G4: the job's threads (bagging runs on threads, no loky processes: BUG-124)


def test_stage_mapping_ngs_only_sample_runs_short_read_mapping_only(tmp_path, monkeypatch):
    """Legacy mapping() ran srMapping alone for samples without long reads (diff_as conditions often are)."""
    cfg = IflasConfig(
        reference=[RefSection(ref_strain="B73", ref_genome="/ref/genome.fa", ref_gtf="/ref/genes.gtf",
                              hisat2_index="/ref/hisat2/B73")],
        dirs={"out_dir": str(tmp_path / "out")},
    )
    sample = Sample(sample_name="n1", ref_strain="B73", project_name="demo", condition="B",
                    ngs_left_reads="/ngs/R1.fq.gz", ngs_right_reads="/ngs/R2.fq.gz")
    job = SampleJob(stage="mapping", config=cfg, sample=sample, ref=cfg.ref_for("B73"),
                    base_dir=sample_base_dir(cfg, sample), threads=2)
    seen = {}

    def fake_ngs(index, out_dir, *, repeats, threads, max_intron):
        seen["ngs"] = (index, out_dir, repeats)
        return {"tmp_bam": f"{out_dir}/reassembly/tmp.bam", "junctions_bed": f"{out_dir}/reassembly/junctions.bed"}

    import iflas.stages.align

    monkeypatch.setattr(iflas.stages.align, "ngs_reassembly_junctions", fake_ngs)
    monkeypatch.setattr(iflas.stages.align, "minimap2_tgs", lambda *a, **k: (_ for _ in ()).throw(AssertionError("no long reads to map")))

    out = _stage_mapping(job)

    assert seen["ngs"][0] == "/ref/hisat2/B73" and seen["ngs"][2] == [{"left": "/ngs/R1.fq.gz", "right": "/ngs/R2.fq.gz"}]
    assert out["flnc"] is None and "tgs_bam" not in out and out["junctions_bed"].endswith("junctions.bed")


def test_default_chain_skips_stages_whose_inputs_are_absent(tmp_path):
    """`iflas all` without --stages runs the whole published chain, skipping (with a reason) ase without NGS reads,
    palen_as without a poly(A) table on a non-nanopore sample and diff_as without per-sample alignments of two
    conditions; an explicit stage list is strict and runs the stage for every sample."""
    cfg = _config(tmp_path)                 # s1: pacbio + NGS, no poly(A) table; s2: nanopore, no NGS
    fast5 = tmp_path / "fast5"
    fast5.mkdir()
    (fast5 / "batch0.fast5").write_text("")
    cfg.samples[1].fast5_dir = str(fast5)   # basecalled reads + raw signal directory (BUG-159 layout)
    calls = []

    def record(job: SampleJob):
        calls.append((job.stage, job.sample.sample_name))
        return {}

    handlers = {stage: record for stage in DEFAULT_STAGES if stage != "diff_as"}
    results = run_pipeline(cfg, stages=None, max_workers=1, handlers=handlers, executor_factory=_ImmediateExecutor)
    skipped = {(r.stage, r.sample_name): r.outputs["skipped"] for r in results if "skipped" in r.outputs}
    assert set(skipped) == {("ase", "s2"), ("palen_as", "s1"), ("diff_as", "project")}
    assert "no NGS reads" in skipped[("ase", "s2")]
    assert "poly(A)" in skipped[("palen_as", "s1")]
    assert "two conditions" in skipped[("diff_as", "project")]
    assert ("ase", "s1") in calls and ("palen_as", "s2") in calls      # nanopore with a fast5 directory: tails can be called

    cfg.samples[1].fast5_dir = None                                       # nanopore sample with basecalled reads only (BUG-159)
    results = run_pipeline(cfg, stages=None, max_workers=1, handlers=handlers, executor_factory=_ImmediateExecutor)
    skipped = {(r.stage, r.sample_name): r.outputs["skipped"] for r in results if "skipped" in r.outputs}
    assert "no fast5 directory" in skipped[("palen_as", "s2")]
    assert [r.stage for r in results if r.sample_name == "s1"] == [s for s in DEFAULT_STAGES if s != "diff_as"]

    calls.clear()
    run_pipeline(cfg, stages=("palen_as", "ase"), max_workers=1, handlers=handlers, executor_factory=_ImmediateExecutor)
    assert calls == [("palen_as", "s1"), ("palen_as", "s2"), ("ase", "s1"), ("ase", "s2")]


def test_rows_without_long_reads_skip_the_long_read_stages(tmp_path):
    """BUG-202: with merging off, an NGS-only row went into collapse and the run stopped; it now skips every
    long-read stage, also when the stages are listed explicitly, and still runs preproc and mapping."""
    cfg = _config(tmp_path)
    cfg.samples.append(Sample(sample_name="s3", project_name="demo", ref_strain="B73",
                              ngs_left_reads="/reads/s3.R1.fq", ngs_right_reads="/reads/s3.R2.fq"))
    calls = []

    def record(job: SampleJob):
        calls.append((job.stage, job.sample.sample_name))
        return {}

    handlers = {stage: record for stage in DEFAULT_STAGES if stage != "diff_as"}
    for stages in (None, ("mapping", "collapse", "refine", "pu_filter", "find_as")):
        calls.clear()
        results = run_pipeline(cfg, stages=stages, max_workers=1, handlers=handlers, executor_factory=_ImmediateExecutor)
        skipped = {r.stage for r in results if r.sample_name == "s3" and "skipped" in r.outputs}
        assert {stage for stage, name in calls if name == "s3"} == ({"preproc", "mapping"} if stages is None else {"mapping"})
        assert skipped == {s for s in (stages or DEFAULT_STAGES) if s not in ("preproc", "mapping", "diff_as")}
        assert ("collapse", "s1") in calls and ("collapse", "s2") in calls


def test_nanopore_reads_map_against_the_genome(tmp_path):
    """A prebuilt minimap2 index fixes the k/w it was built with, so nanopore samples map against the FASTA, where the
    command's -k14 and the splice preset's -w5 apply (BUG-192); PacBio keeps the index."""
    from types import SimpleNamespace

    from iflas.orchestrator import _mapping_target

    ref = SimpleNamespace(ref_genome="/ref/genome.fa", ref_mm2_index="/ref/genome.mm2")
    config = SimpleNamespace(minimap2=SimpleNamespace(mm2_index=None))
    job = lambda plat: SimpleNamespace(sample=SimpleNamespace(tgs_plat=plat), ref=ref, config=config)
    assert _mapping_target(job("nanopore")) == "/ref/genome.fa"
    assert _mapping_target(job("pacbio")) == "/ref/genome.mm2"


# ----------------------------------------------------------------------------- BUG-198 (real-reference check)
def _strain_config(tmp_path, samples):
    return IflasConfig(reference=[RefSection(ref_strain="B73", ref_genome="/ref/genome.fa", ref_gtf="/ref/genes.gtf")],
                       options=Options(threads=2, merge_data_from_same_strain=True), samples=samples,
                       dirs={"out_dir": str(tmp_path / "out")})


def test_runs_without_preproc_use_what_preproc_produced(tmp_path, monkeypatch):
    """BUG-198: preproc's results (processed FLNC, fastp-filtered NGS, dorado's tail table) reached the later stages
    only within the same run; a resumed `iflas all --stages ...` used the raw NGS reads, had no FLNC for samples
    given as raw data (asas stopped) and lost the dorado tails. preproc now writes them to
    preprocess/preproc_outputs.json and runs without preproc read them back."""
    import json

    import iflas.stages.preprocess as preprocess
    from iflas.orchestrator import _stage_preproc
    samples = [Sample(sample_name="pb", project_name="demo", ref_strain="B73", strain="B73", tgs_plat="pacbio",
                      data_location="/raw/pb.bam", ngs_left_reads="/raw/a_R1.fq.gz", ngs_right_reads="/raw/a_R2.fq.gz"),
               Sample(sample_name="ont", project_name="demo", ref_strain="B73", strain="B73", tgs_plat="nanopore",
                      data_location="/raw/dorado.bam")]
    cfg = _strain_config(tmp_path, samples)
    tails = tmp_path / "dorado_polya.tsv"
    tails.write_text("readname\tpolya_length\tqc_tag\n")
    produced = {"pb": {"flnc": "/p/pb.fq", "ngs_left_reads": "/p/fastp.a_R1.fq.gz", "ngs_right_reads": "/p/fastp.a_R2.fq.gz",
                       "polya_location": None},
                "ont": {"flnc": "/p/ont.fq", "ngs_left_reads": None, "ngs_right_reads": None, "polya_location": str(tails)}}
    monkeypatch.setattr(preprocess, "preprocess_sample", lambda sample, *a, **k: produced[sample.sample_name])
    for sample in samples:
        _stage_preproc(SampleJob("preproc", cfg, sample, cfg.reference[0], sample_base_dir(cfg, sample), 1))
        assert json.load(open(os.path.join(sample_base_dir(cfg, sample), "preprocess", "preproc_outputs.json"))) == produced[sample.sample_name]

    seen = []
    resumed = _strain_config(tmp_path, [s.model_copy() for s in samples])          # a new invocation: the sheet as written
    run_pipeline(resumed, stages=("asas",), merge_samples=True, executor_factory=_ImmediateExecutor,
                 handlers={"asas": lambda job: seen.append(job.sample) or {}})
    (merged,) = seen
    assert merged.data_processed_location == ["/p/pb.fq", "/p/ont.fq"]
    assert (merged.ngs_left_reads, merged.ngs_right_reads) == ("/p/fastp.a_R1.fq.gz", "/p/fastp.a_R2.fq.gz")
    assert merged.polya_location == {"ont": str(tails)}


def test_a_strain_missing_the_flnc_of_one_long_read_sample_stops(tmp_path):
    """BUG-198: the merged FLNC left out a long-read sample whose reads were not found (legacy did the same), so its reads
    were silently missing from every later stage; an NGS-only sample is not part of that list."""
    import pytest

    from iflas.orchestrator import _concatenate_flnc_inputs, _merged_samples
    flnc = tmp_path / "pb.fa"
    flnc.write_text(">r1\nACGT\n")
    samples = [Sample(sample_name="ngs", project_name="demo", ref_strain="B73", strain="B73", ngs_left_reads="/n_R1.fq.gz"),
               Sample(sample_name="pb", project_name="demo", ref_strain="B73", strain="B73", tgs_plat="pacbio",
                      data_processed_location=str(flnc)),
               Sample(sample_name="ont", project_name="demo", ref_strain="B73", strain="B73", tgs_plat="nanopore",
                      data_location="/raw/ont.fq")]
    cfg = _strain_config(tmp_path, samples)
    (merged,) = _merged_samples(cfg)
    assert merged.tgs_plat == "pacbio" and merged.data_processed_location == [str(flnc), None]   # legacy: None (NGS row first)
    job = SampleJob("mapping", cfg, merged, cfg.reference[0], sample_base_dir(cfg, merged), 1)
    with pytest.raises(ValueError, match="FLNC reads missing for 1 long-read sample"):
        _concatenate_flnc_inputs(job, merged.data_processed_location)


def test_a_merged_strain_reads_each_ont_runs_signal_against_its_own_reads(tmp_path, monkeypatch):
    """BUG-200: in a strain merging an RNA004 run (dorado's tail table) with an RNA002 run (fast5), the tail table made
    palen_as skip nanopolish, so the RNA002 tails were lost; without a table nanopolish read the signal against the
    strain's processed (fmlrc2-corrected, poly(A)-trimmed) FLNC. Each run's own basecalled reads are aligned for it now."""
    import iflas.runtime as runtime
    import iflas.stages.align as align
    import iflas.stages.palen_as as palen
    from iflas.orchestrator import _merged_samples, _stage_palen_as
    tails, fast5 = tmp_path / "dorado_polya.tsv", tmp_path / "fast5"
    tails.write_text("readname\tpolya_length\tqc_tag\n")
    fast5.mkdir()
    samples = [Sample(sample_name="ont004", project_name="demo", ref_strain="B73", strain="B73", tgs_plat="nanopore",
                      data_location="/raw/dorado.bam", polya_location=str(tails)),
               Sample(sample_name="ont002", project_name="demo", ref_strain="B73", strain="B73", tgs_plat="nanopore",
                      data_location="/raw/reads.fastq", fast5_dir=str(fast5))]
    cfg = _strain_config(tmp_path, samples)
    raw = tmp_path / "out" / "demo" / "ont002" / "preprocess" / "nanopore" / "rawFlnc.fq"
    raw.parent.mkdir(parents=True)
    raw.write_text("@r1\nACGU\n+\nIIII\n")
    calls = []

    def fake_minimap2(target, reads, out_sam, **kw):
        os.makedirs(os.path.dirname(out_sam), exist_ok=True)
        open(out_sam, "w").close()
        calls.append(("minimap2", reads, kw["platform"]))
        return out_sam
    monkeypatch.setattr(align, "minimap2_tgs", fake_minimap2)
    monkeypatch.setattr(align, "index_bam", lambda bam, runner: calls.append(("index", os.path.basename(bam))))
    monkeypatch.setattr(runtime.ToolRunner, "run", lambda self, cmd, **kw: calls.append(("run", cmd.split()[0])))
    monkeypatch.setattr(palen, "polya_len_calling", lambda f5, reads, bam, genome, out_dir, **kw:
                        calls.append(("nanopolish", f5, reads)) or str(tmp_path / "polya_results.tsv"))
    monkeypatch.setattr(palen, "palen_as_sample", lambda base, files, **kw:
                        calls.append(("palen_as", files, kw["fast5_dir"])) or {"palen_dir": "p", "as_dir": "a"})
    (merged,) = _merged_samples(cfg)
    _stage_palen_as(SampleJob("palen_as", cfg, merged, cfg.reference[0], sample_base_dir(cfg, merged), 2))
    assert calls[0] == ("minimap2", str(raw), "nanopore")
    assert ("nanopolish", str(fast5), str(raw)) in calls
    assert calls[-1] == ("palen_as", [str(tails), str(tmp_path / "polya_results.tsv")], None)
