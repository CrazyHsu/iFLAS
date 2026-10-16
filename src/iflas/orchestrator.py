"""Modern pipeline orchestrator.

This replaces the legacy ``iflas.py`` split dispatch / missing ``oneCommandRun`` path
with a small stage-barrier scheduler. The stage handlers only wire migrated drivers;
algorithmic work stays inside the stage modules.
"""

from __future__ import annotations

import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

from .config import IflasConfig, MergedTgsSample, RefSection, Sample, merge_samples_by_strain

# the published chain, in order; `iflas all` runs it when no --stages is given (legacy `all` never ran: BUG-105)
DEFAULT_STAGES = ("preproc", "mapping", "collapse", "refine", "pu_filter", "find_as", "find_pa", "characterize",
                  "visual_as", "diff_as", "asas", "ase", "palen_as")
PipelineSample = Sample | MergedTgsSample


@dataclass(frozen=True)
class SampleJob:
    stage: str
    config: IflasConfig
    sample: PipelineSample
    ref: RefSection
    base_dir: str
    threads: int


@dataclass(frozen=True)
class StageResult:
    stage: str
    sample_name: str
    base_dir: str
    outputs: dict


StageHandler = Callable[[SampleJob], dict]


def sample_base_dir(config: IflasConfig, sample: PipelineSample) -> str:
    """Legacy output layout: ``<out_dir>/<project_name>/<sample_name>``."""
    return os.path.abspath(os.path.join(config.dirs.out_dir, sample.project_name, sample.sample_name))


def _split_repeats(value: str | None) -> list[str]:
    if not value:
        return []
    repeats = []
    for part in str(value).split(";"):
        reads = ",".join(token.strip() for token in part.split(",") if token.strip())
        if reads:
            repeats.append(reads)
    return repeats


def parse_ngs_repeats(sample: Sample) -> list[dict[str, str]]:
    """Parse legacy NGS repeat syntax into ``ngs_reassembly_junctions`` inputs."""
    paired = (sample.ngs_reads_paired or "paired").lower()
    left = _split_repeats(sample.ngs_left_reads)
    right = _split_repeats(sample.ngs_right_reads)
    if paired == "paired":
        if not left and not right:
            return []
        if len(left) != len(right):
            raise ValueError(
                f"sample {sample.sample_name!r}: paired NGS repeats require equal left/right counts"
            )
        return [{"left": l, "right": r} for l, r in zip(left, right)]

    if left and right:
        raise ValueError(
            f"sample {sample.sample_name!r}: single-end NGS reads require exactly one read side"
        )
    single = left or right
    return [{"single": reads} for reads in single]


def _fast5_dirs(sample: PipelineSample) -> list[str]:
    """Raw ONT signal directories for nanopolish: the ``fast5_dir`` column(s), else a nanopore ``data_location``
    that is itself a directory of fast5 files (the published layout, where basecalling ran inside the pipeline)."""
    from .stages.palen_as import check_fast5_files

    dirs = sample.fast5_dir
    dirs = list(dirs.values()) if isinstance(dirs, dict) else ([dirs] if dirs else [])
    if not dirs and (sample.tgs_plat or "").lower() == "nanopore":
        locations = sample.data_location if isinstance(sample.data_location, list) else [sample.data_location]
        dirs = [d for d in locations if d and check_fast5_files(d) == "fast5"]
    return dirs


LONG_READ_STAGES = frozenset({"collapse", "refine", "pu_filter", "find_as", "find_pa", "characterize", "visual_as",
                              "asas", "ase", "palen_as"})


def no_long_reads(stage: str, sample: PipelineSample) -> str | None:
    """A row with NGS reads only has nothing for the long-read stages; it skips them whether or not they were listed
    explicitly (BUG-202: collapse stopped on it with "Input file None does not exist")."""
    flnc = sample.data_processed_location
    flnc = [f for f in flnc if f] if isinstance(flnc, list) else flnc
    if stage in LONG_READ_STAGES and not flnc and (sample.tgs_plat or "").lower() not in {"pacbio", "nanopore"}:
        return "no long reads (tgs_plat, data_location or data_processed_location)"
    return None


def _skip_reason(stage: str, sample: PipelineSample) -> str | None:
    """Why ``stage`` cannot run for ``sample``: an input it needs is absent. The default chain skips it and says so."""
    if stage == "ase" and not parse_ngs_repeats(sample):
        return "no NGS reads (ngs_left_reads / ngs_right_reads)"
    if stage == "palen_as" and not sample.polya_location and not _fast5_dirs(sample):
        return "no poly(A) length table (polya_location) and no fast5 directory (fast5_dir, or a fast5 data_location)"
    return None


def _project_skip_reason(config: IflasConfig, stage: str) -> str | None:
    if stage != "diff_as":
        return None
    from .stages.diff_as import bam_records_from_samples

    records = bam_records_from_samples(config)
    if len({r.condition for r in records}) < 2:
        return "needs NGS alignments from at least two conditions"
    if any(not os.path.isfile(r.bam_file) for r in records):
        return ("per-sample NGS alignments not found: merge_data_from_same_strain pools a strain's libraries into one "
                "alignment; run `iflas mapping` with merging off, then `iflas diff_as`")
    return None


def _skipped(stage: str, name: str, base_dir: str, reason: str) -> StageResult:
    print(f"Skip {stage} for {name}: {reason}")
    return StageResult(stage=stage, sample_name=name, base_dir=base_dir, outputs={"skipped": reason})


def _normalize_stages(stages: Sequence[str] | None) -> tuple[str, ...]:
    selected = tuple(stages or DEFAULT_STAGES)
    unknown = [s for s in selected if s not in DEFAULT_HANDLERS and s not in PROJECT_HANDLERS]
    if unknown:
        raise ValueError("unknown stage(s): " + ", ".join(unknown))
    return selected


def _workers(samples: Sequence[PipelineSample], max_workers: int | None) -> int:
    if not samples:
        return 0
    if max_workers is None:
        return len(samples)
    return max(1, min(int(max_workers), len(samples)))


def _threads_per_job(config: IflasConfig, workers: int) -> int:
    total = max(1, int(config.options.threads))
    return max(1, total // max(1, workers))


def _merged_samples(config: IflasConfig) -> list[MergedTgsSample]:
    merged = merge_samples_by_strain(config.samples)
    samples: list[MergedTgsSample] = []
    for by_ref in merged.values():
        for by_strain in by_ref.values():
            samples.extend(by_strain.values())
    return samples


def _samples_for_stage(
    config: IflasConfig,
    stage: str,
    *,
    merge_samples: bool,
) -> list[PipelineSample]:
    if merge_samples and stage != "preproc":
        return _merged_samples(config)
    return list(config.samples)


def _jobs_for_stage(
    config: IflasConfig,
    stage: str,
    workers: int,
    samples: Sequence[PipelineSample],
) -> list[SampleJob]:
    threads = _threads_per_job(config, workers)
    jobs = []
    for sample in samples:
        jobs.append(
            SampleJob(
                stage=stage,
                config=config,
                sample=sample,
                ref=config.ref_for(sample.ref_strain),
                base_dir=sample_base_dir(config, sample),
                threads=threads,
            )
        )
    return jobs


def run_pipeline(
    config: IflasConfig,
    *,
    stages: Sequence[str] | None = None,
    max_workers: int | None = None,
    merge_samples: bool = False,
    handlers: Mapping[str, StageHandler] | None = None,
    executor_factory=ProcessPoolExecutor,
    skip_unavailable: bool | None = None,
) -> list[StageResult]:
    """Run selected migrated stages, one stage barrier at a time.

    ``merge_samples=True`` restores the legacy project/ref_strain/strain grouping for
    migrated analysis stages. Explicit ``preproc`` still runs per original sample so its
    produced FLNC paths can be merged before mapping.

    ``skip_unavailable`` (default: only when ``stages`` is None, i.e. the whole chain) skips a stage for a sample
    that lacks its inputs (``ase`` without NGS, ``palen_as`` without a poly(A) table, ``diff_as`` without
    per-sample alignments of two conditions) and records ``{"skipped": reason}``; explicitly requested stages must run.
    Rows without long reads skip the long-read stages in either case (``no_long_reads``).
    """
    if config.dirs.tmp_dir:                        # BUG-184: temporary files go where the config says
        from .runtime import use_tmp_dir
        use_tmp_dir(config.dirs.tmp_dir)
    selected = _normalize_stages(stages)
    if skip_unavailable is None:
        skip_unavailable = stages is None
    if not config.samples:
        raise ValueError("pipeline config has no samples")

    if "preproc" not in selected:
        restore_preproc_outputs(config)
    stage_handlers = dict(DEFAULT_HANDLERS)
    if handlers:
        stage_handlers.update(handlers)

    results: list[StageResult] = []
    for stage in selected:
        if stage in PROJECT_HANDLERS:          # once per project, not per sample (diff_as)
            reason = _project_skip_reason(config, stage) if skip_unavailable else None
            if reason:
                results.append(_skipped(stage, "project", config.dirs.out_dir, reason))
                continue
            outputs = PROJECT_HANDLERS[stage](config)
            results.append(StageResult(stage=stage, sample_name="project", base_dir=outputs.get("das_dir", config.dirs.out_dir),
                                       outputs=outputs))
            continue
        runnable = []
        for sample in _samples_for_stage(config, stage, merge_samples=merge_samples):
            reason = no_long_reads(stage, sample) or (_skip_reason(stage, sample) if skip_unavailable else None)
            if reason:
                results.append(_skipped(stage, sample.sample_name, sample_base_dir(config, sample), reason))
            else:
                runnable.append(sample)
        stage_samples = runnable
        if not stage_samples:
            continue
        workers = _workers(stage_samples, max_workers)
        handler = stage_handlers[stage]
        jobs = _jobs_for_stage(config, stage, workers, stage_samples)
        stage_results: list[StageResult] = []
        with executor_factory(max_workers=workers) as executor:
            futures = [(job, executor.submit(handler, job)) for job in jobs]
            for job, future in futures:
                outputs = future.result()
                result = StageResult(
                    stage=stage,
                    sample_name=job.sample.sample_name,
                    base_dir=job.base_dir,
                    outputs=outputs or {},
                )
                results.append(result)
                stage_results.append(result)
        _apply_stage_outputs(stage_samples, stage, stage_results)
    return results


def _apply_stage_outputs(
    samples: Sequence[PipelineSample],
    stage: str,
    results: Sequence[StageResult],
) -> None:
    if stage != "preproc":
        return
    by_name = {sample.sample_name: sample for sample in samples}
    for result in results:
        if result.sample_name in by_name:
            _apply_preproc_outputs(by_name[result.sample_name], result.outputs)


def _apply_preproc_outputs(sample: PipelineSample, outputs: dict) -> None:
    """What preproc hands to the later stages: the processed FLNC, the fastp-filtered NGS reads, dorado's tail table."""
    if outputs.get("flnc"):
        sample.data_processed_location = outputs["flnc"]
    if outputs.get("ngs_left_reads"):
        sample.ngs_left_reads = outputs["ngs_left_reads"]
    if outputs.get("ngs_right_reads"):
        sample.ngs_right_reads = outputs["ngs_right_reads"]
    if outputs.get("polya_location") and not sample.polya_location:     # dorado's tail table (F-11)
        sample.polya_location = outputs["polya_location"]


PREPROC_OUTPUTS = os.path.join("preprocess", "preproc_outputs.json")


def restore_preproc_outputs(config: IflasConfig) -> None:
    """Apply each sample's last preproc results when a run does not preprocess (a resumed ``iflas all --stages ...``,
    ``iflas mapping`` after ``iflas preproc``). BUG-198: such runs used the raw NGS reads instead of the fastp-filtered
    ones, had no FLNC for samples given as raw data and lost dorado's tail table."""
    for sample in config.samples:
        path = os.path.join(sample_base_dir(config, sample), PREPROC_OUTPUTS)
        if os.path.isfile(path):
            with open(path) as f:
                _apply_preproc_outputs(sample, json.load(f))


def _mapping_target(job: SampleJob) -> str:
    # nanopore reads map against the genome itself: a prebuilt index fixes the k/w it was built with, so the command's
    # -k14 and the splice preset's -w5 never applied (BUG-192, adopted 2026-10-16; PacBio alignments are identical)
    if (job.sample.tgs_plat or "").lower() == "nanopore":
        return job.ref.ref_genome
    return job.config.minimap2.mm2_index or job.ref.ref_mm2_index or job.ref.ref_genome


def _ngs_junctions(job: SampleJob) -> str | None:
    if job.sample.ngs_junctions:
        return job.sample.ngs_junctions
    detected = os.path.join(job.base_dir, "mapping", "rna-seq", "reassembly", "junctions.bed")
    return detected if os.path.isfile(detected) else None


def _ngs_tmp_bam(job: SampleJob) -> str | None:
    detected = os.path.join(job.base_dir, "mapping", "rna-seq", "reassembly", "tmp.bam")
    return detected if os.path.isfile(detected) else None


def _collapsed_gff(job: SampleJob) -> str:
    for name in ("tofu.collapsed.good.gff", "tofu.collapsed.gff"):
        path = os.path.join(job.base_dir, "collapse", name)
        if os.path.isfile(path):
            return path
    return os.path.join(job.base_dir, "collapse", "tofu.collapsed.good.gff")


def _read_kind(path: str) -> bytes:
    """First byte of a sequence file: b">" FASTA, b"@" FASTQ, the gzip magic otherwise."""
    with open(path, "rb") as f:
        return f.read(1)


def _concatenate_flnc_inputs(job: SampleJob, paths: Sequence[str | None]) -> str | None:
    """The FLNC of a strain's rows as one file. Plain files of one format are concatenated as they are; a mix (PacBio
    FASTA with ONT FASTQ, or gzipped input) is written as FASTA by seqkit, since a byte concatenation of FASTA and
    FASTQ records is read as one format and loses the other's reads (BUG-177)."""
    from .runtime import ToolRunner, quote

    if not paths:                           # a strain sequenced with NGS only
        return None
    missing = [path for path in paths if not (path and os.path.isfile(path))]
    if missing:                             # legacy left them out silently, losing those samples' reads
        raise ValueError(f"sample {job.sample.sample_name!r}: FLNC reads missing for {len(missing)} long-read sample(s) of the "
                         f"strain ({', '.join(map(str, missing))}); run preproc first or set their data_processed_location")
    valid_files = list(paths)
    platform = (job.sample.tgs_plat or "pacbio").lower()
    out_flnc = os.path.join(job.base_dir, "preprocess", platform, "rawFlnc.fq")
    os.makedirs(os.path.dirname(out_flnc), exist_ok=True)
    kinds = {_read_kind(path) for path in valid_files}
    if len(kinds) == 1 and kinds <= {b">", b"@"}:
        with open(out_flnc, "wb") as out:
            for path in valid_files:
                with open(path, "rb") as src:
                    shutil.copyfileobj(src, out)
    else:
        ToolRunner().run("({}) > {}".format("; ".join(f"seqkit fq2fa -w 0 {quote(p)}" for p in valid_files), quote(out_flnc)))
    return out_flnc


def _prepared_flnc(job: SampleJob) -> str | None:
    flnc = job.sample.data_processed_location
    if isinstance(flnc, list):
        return _concatenate_flnc_inputs(job, flnc)
    return flnc or _preprocessed_flnc(job)


def _raw_flnc(job: SampleJob, flnc: str) -> str | None:
    """Uncorrected FLNC written by preproc; legacy maps it too (raw-vs-corrected mapping report)."""
    platform = (job.sample.tgs_plat or "pacbio").lower()
    raw = os.path.join(job.base_dir, "preprocess", platform, "rawFlnc.fq")
    if os.path.isfile(raw) and not os.path.samefile(raw, flnc):
        return raw
    return None


def _stage_mapping(job: SampleJob) -> dict:
    """Legacy ``mapping()`` order: NGS first (its junctions feed the FLNC junction-combination filter), then
    minimap2 -> SAM -> ``mapping_filter_and_add_tags`` -> ``flnc.mm2.sorted.bam`` + ``flnc.addCVandID.bed12+``
    (BUG-143), and ``flnc.processed.fa`` for collapse."""
    from .stages.align import (flnc_to_fasta, genome_bases, mapping_filter_and_add_tags, minimap2_index_batch,
                               minimap2_tgs, ngs_reassembly_junctions)

    flnc = _prepared_flnc(job)
    repeats = parse_ngs_repeats(job.sample)
    if not flnc and not repeats:
        raise ValueError(
            f"sample {job.sample.sample_name!r}: mapping requires data_processed_location, a prior preproc "
            "stage, or NGS reads"
        )
    if flnc and not job.ref.ref_bed:
        raise ValueError(
            f"sample {job.sample.sample_name!r}: mapping requires reference.ref_bed for the junction "
            "filter (run `iflas init-ref`)"
        )

    mapping_dir = os.path.join(job.base_dir, "mapping")
    os.makedirs(mapping_dir, exist_ok=True)
    outputs = {"flnc": flnc}
    if repeats:
        if not job.ref.hisat2_index:
            raise ValueError(
                f"sample {job.sample.sample_name!r}: NGS mapping requires reference.hisat2_index"
            )
        outputs.update(ngs_reassembly_junctions(
            job.ref.hisat2_index,
            os.path.join(mapping_dir, "rna-seq"),
            repeats=repeats,
            threads=job.threads,
            max_intron=job.config.minimap2.max_intron_length,
        ))

    if not flnc:                       # NGS-only sample (legacy mapping() ran srMapping alone)
        return outputs
    max_intron = job.config.minimap2.max_intron_length
    junctions = _ngs_junctions(job)
    for prefix, reads in (("flnc", flnc), ("rawFlnc", _raw_flnc(job, flnc))):
        if not reads:
            continue
        sam = minimap2_tgs(
            _mapping_target(job),
            reads,
            os.path.join(mapping_dir, f"{prefix}.mm2.sam"),
            platform=(job.sample.tgs_plat or "pacbio").lower(),
            max_intron=max_intron,
            threads=job.threads,
            index_batch=minimap2_index_batch(genome_bases(job.ref.ref_size or job.ref.ref_genome + ".fai")),
        )
        filtered = mapping_filter_and_add_tags(
            sam_file=sam,
            out_prefix=prefix,
            max_length=max_intron,
            ref_bed_file=job.ref.ref_bed,
            junc_bed_file=junctions,
            threads=job.threads,
            junc_comb_sup=job.config.minimap2.junc_comb_support,
            work_dir=mapping_dir,
        )
        outputs[f"{prefix}_bed12plus"] = filtered["bed12plus"]
    outputs["tgs_bam"] = os.path.join(mapping_dir, "flnc.mm2.sorted.bam")
    outputs["processed_fa"] = flnc_to_fasta(flnc, os.path.join(mapping_dir, "flnc.processed.fa"))
    return outputs


def _stage_collapse(job: SampleJob) -> dict:
    from .stages.collapse import collapse_isoforms

    processed = os.path.join(job.base_dir, "mapping", "flnc.processed.fa")   # legacy collapse input
    flnc = processed if os.path.isfile(processed) else _prepared_flnc(job)
    tgs_bam = os.path.join(job.base_dir, "mapping", "flnc.mm2.sorted.bam")
    out_prefix = os.path.join(job.base_dir, "collapse", "tofu")
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    gff = collapse_isoforms(
        flnc,
        tgs_bam,
        out_prefix,
        max_5_diff=job.config.collapse.max_5_diff,
        max_3_diff=job.config.collapse.max_3_diff,
        fl_coverage=job.config.collapse.fl_coverage,
        min_coverage=job.config.collapse.min_coverage,
        min_identity=job.config.collapse.min_identity,
        max_fuzzy_junction=job.config.collapse.max_fuzzy_junction,
        dun_merge_5_shorter=job.config.collapse.dun_merge_5_shorter,
    )
    return {"collapsed_gff": gff, "group": out_prefix + ".collapsed.group.txt"}


def _stage_refine(job: SampleJob) -> dict:
    from .stages.refine import refine_isoforms

    return refine_isoforms(
        _collapsed_gff(job),
        job.ref.ref_genome,
        job.ref.ref_gpe,
        job.ref.ref_bed,
        job.base_dir,
        ngs_junctions=_ngs_junctions(job),
    )


def _stage_pu_filter(job: SampleJob) -> dict:
    from .stages.pu_filter import filter_hq_isoforms

    return filter_hq_isoforms(
        job.base_dir,
        job.ref.ref_genome,
        job.ref.ref_bed,
        sam_file=_ngs_tmp_bam(job),
        junction_file=_ngs_junctions(job),
        pos_fl_coverage=job.config.pu_filter.pos_fl_coverage,
        pos_min_junc_rpkm=job.config.pu_filter.pos_min_junc_rpkm,
        filter_score=job.config.pu_filter.filter_score,
        feature_file=job.config.pu_filter.feature_file,
        bagging_n_jobs=job.threads,          # bagging runs its jobs as threads: no processes inside the pool worker (BUG-124)
        random_state=job.config.pu_filter.random_state,
        select_best_model=job.config.pu_filter.select_best_model,
        auto_filter_score=job.config.pu_filter.auto_filter_score,
        draw_auc=job.config.pu_filter.draw_auc,
    )


def _stage_find_as(job: SampleJob) -> dict:
    from .stages.find_as import find_as_sample

    return find_as_sample(job.base_dir, job.ref.ref_gpe, ngs_junctions=_ngs_junctions(job))


def _stage_find_pa(job: SampleJob) -> dict:
    from .stages.find_pa import find_pa_sample

    if not job.ref.ref_size:
        raise ValueError(
            f"sample {job.sample.sample_name!r}: find_pa requires reference.ref_size (genome.len.txt from `iflas init-ref`)"
        )
    return find_pa_sample(job.base_dir, job.ref.ref_genome, job.ref.ref_size, threads=job.threads)


def _stage_characterize(job: SampleJob) -> dict:
    from .stages.characterize import characterize_sample

    return characterize_sample(job.base_dir, ref_gpe=job.ref.ref_gpe, ref_genome=job.ref.ref_genome)


def _stage_visual_as(job: SampleJob) -> dict:
    from .stages.visual_as import visual_as_sample

    plots = visual_as_sample(job.base_dir, job.ref.ref_gpe, threads=job.threads)
    return {"visual_as_dir": os.path.join(job.base_dir, "visual_as"), "genes": len(plots)}


def _stage_asas(job: SampleJob) -> dict:
    from .stages.asas import asas_sample

    flnc = _prepared_flnc(job)
    if not flnc:
        raise ValueError(f"sample {job.sample.sample_name!r}: asas needs the FLNC reads (data_processed_location or preproc)")
    out = asas_sample(job.base_dir, ref_genome=job.ref.ref_genome, flnc=flnc, threads=job.threads,
                      ngs_confirm=job.config.options.asas_ngs_confirm)
    return {"asas_dir": str(out["asas_dir"]), "partial": str(out["partial"]), "loci": len(out["loci"])}


def _stage_ase(job: SampleJob) -> dict:
    from .stages.asas import ase_sample

    s = job.sample
    out = ase_sample(job.base_dir, job.ref.ref_genome, ngs_left=s.ngs_left_reads, ngs_right=s.ngs_right_reads,
                     paired=(s.ngs_reads_paired or "paired") == "paired", threads=job.threads)
    return {"ase_dir": out["ase_dir"], "merged": out["merged"]}


def _stage_palen_as(job: SampleJob) -> dict:
    from .stages.palen_as import palen_as_sample

    s = job.sample
    polya = s.polya_location
    palen_files = list(polya.values()) if isinstance(polya, dict) else ([polya] if polya else [])
    if isinstance(s.fast5_dir, dict):          # a merged strain: each ONT run's signal against its own reads (BUG-200)
        palen_files += _run_nanopolish_tables(job, skip=set(polya or {}))
        fast5 = None
    else:
        fast5 = None if palen_files else (_fast5_dirs(s) or None)
    summary = os.path.join(job.base_dir, "preprocess", "nanopore", "guppy_basecalling", "sequencing_summary.txt")
    out = palen_as_sample(job.base_dir, palen_files, fast5_dir=fast5, ref_genome=job.ref.ref_genome,
                          sequencing_summary=summary if os.path.isfile(summary) else None, threads=job.threads)
    return {"palen_dir": str(out["palen_dir"]), "as_dir": str(out["as_dir"]), "apa": str(out["apa"]) if out.get("apa") else None}


def _run_nanopolish_tables(job: SampleJob, skip: set) -> list[str]:
    """nanopolish tail tables of the ONT runs of a merged strain that brought raw signal and no tail table. Legacy and
    BUG-200: a strain with any tail table (e.g. dorado's) skipped nanopolish altogether, and without one it read the signal
    against the strain's processed FLNC (fmlrc2-corrected and poly(A)-trimmed copies; BUG-175). Each run's basecalled
    reads (``preprocess/nanopore/rawFlnc.fq``) are aligned on their own here, as the reads nanopolish indexes."""
    from .runtime import ToolRunner, quote
    from .stages.align import genome_bases, index_bam, minimap2_index_batch, minimap2_tgs
    from .stages.palen_as import polya_len_calling

    runner, tables = ToolRunner(), []
    for run, fast5 in sorted(job.sample.fast5_dir.items()):
        if run in skip:
            continue
        run_dir = os.path.abspath(os.path.join(job.config.dirs.out_dir, job.sample.project_name, run))
        raw = os.path.join(run_dir, "preprocess", "nanopore", "rawFlnc.fq")
        if not os.path.isfile(raw):
            raise ValueError(f"sample {run!r}: nanopolish needs its basecalled reads ({raw}); run preproc first")
        work = os.path.join(job.base_dir, "palenAS", "nanopolish", run)
        sam = minimap2_tgs(job.ref.ref_genome, raw, os.path.join(work, "rawFlnc.sam"), platform="nanopore",
                           max_intron=job.config.minimap2.max_intron_length, threads=job.threads,
                           index_batch=minimap2_index_batch(genome_bases(job.ref.ref_size or job.ref.ref_genome + ".fai")))
        bam = os.path.join(work, "rawFlnc.sorted.bam")
        runner.run(f"samtools sort -@ {job.threads} -o {quote(bam)} {quote(sam)}")
        index_bam(bam, runner=runner)
        os.remove(sam)
        summary = os.path.join(run_dir, "preprocess", "nanopore", "guppy_basecalling", "sequencing_summary.txt")
        tables.append(polya_len_calling(fast5, raw, bam, job.ref.ref_genome, work,
                                        sequencing_summary=summary if os.path.isfile(summary) else None, threads=job.threads))
    return tables


def _stage_preproc(job: SampleJob) -> dict:
    from .stages.preprocess import preprocess_sample

    outputs = preprocess_sample(
        job.sample,
        job.config.ccs,
        job.base_dir,
        threads=job.threads,
        trim_polya=job.config.options.ont_trim_polya,
    )
    path = os.path.join(job.base_dir, PREPROC_OUTPUTS)                  # read back by runs without preproc
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(outputs, f, indent=1)
    return outputs


def _preprocessed_flnc(job: SampleJob) -> str | None:
    from .stages.preprocess import find_preprocessed_flnc

    return find_preprocessed_flnc(job.base_dir)


def _project_diff_as(config: IflasConfig) -> dict:
    from .stages.diff_as import diff_as_project

    return diff_as_project(config, threads=config.options.threads)


PROJECT_HANDLERS: dict[str, Callable[[IflasConfig], dict]] = {"diff_as": _project_diff_as}

DEFAULT_HANDLERS: dict[str, StageHandler] = {
    "preproc": _stage_preproc,
    "mapping": _stage_mapping,
    "collapse": _stage_collapse,
    "refine": _stage_refine,
    "pu_filter": _stage_pu_filter,
    "find_as": _stage_find_as,
    "find_pa": _stage_find_pa,
    "characterize": _stage_characterize,
    "visual_as": _stage_visual_as,
    "asas": _stage_asas,
    "ase": _stage_ase,
    "palen_as": _stage_palen_as,
}
