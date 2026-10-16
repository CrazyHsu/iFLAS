"""File-level report plot helpers migrated from legacy ``report.py``."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Iterable

import pandas as pd

from ..config import IflasConfig, RefSection, Sample, merge_samples_by_strain
from ..io.files import getDictFromFile, validateFaAndFqFile
from ..io.records import BedFile
from ..rscript import r_literal, run_r_calls, run_r_function
from ..runtime import ToolRunner, get_logger, quote as q
from ..stages.align import BAI_MAX_LENGTH, index_bam, longest_reference

RSAMTOOLS_MAX_END = 2 ** 29          # Bioconductor's Rsamtools (Gviz alignment tracks) refuses regions ending later
from . import legacy_r

AS_TYPES = ("IR", "SE", "A3SS", "A5SS")
ASAS_HIGHLIGHT_COLORS = {"IR": "#0000FF", "SE": "#00EE00", "A3SS": "#FFA500", "A5SS": "#DD77FF"}

PLOT_READS_CORRECTED_EVAL_SOURCE = r'''
    library(ggplot2)
    library(scales)
    plotReadsCorrectedEval <- function(rawMappedBed, correctMappedBed, outPdf){
        rawMapped <- read.csv(rawMappedBed, sep="\t", header=FALSE)
        correctMapped <- read.csv(correctMappedBed, sep="\t", header=FALSE)
        rawMapped$accuracy <- (rawMapped$V13 * rawMapped$V14)/100
        correctMapped$accuracy <- (correctMapped$V13 * correctMapped$V14)/100
        rawMapped$type <- "Raw"
        correctMapped$type <- "Corrected"
        pdf(outPdf, height=6, width=10)
        p <- ggplot(rbind(rawMapped,correctMapped), aes(x=accuracy, fill=type)) +
        geom_histogram(alpha=0.6, position = 'identity', binwidth=1) +
        scale_fill_manual(values=c("#ff6666", "#C0C0C0")) +
        xlab("Accuracy (%)") + ylab("Count (x1e4)") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size=28), text = element_text(size=28),
              panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"), strip.background = element_blank(), axis.text = element_text(size=24),
              strip.text.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 24),
              legend.position=c(0.2, 0.8), legend.key.size = unit(1, "cm"), panel.spacing = unit(0, "lines")) +
        scale_y_continuous(labels = label_number(scale = 1 / 1e4)) +
        scale_x_continuous(limits=c(50,101))
        print(p)
        dev.off()
    }
'''

PLOT_NOVEL_HQ_AS_SOURCE = r'''
    library(ggplot2)
    plotNovelHqAS <- function(isoformScoreFile, outPdf){
        isoformScore <- read.csv(isoformScoreFile, sep="\t", header=FALSE)
        names(isoformScore) <- c("Gene", "Isos", "Count", "Total_count", "Freq", "Annotation")
        isoformScore <- isoformScore[order(-isoformScore$Freq),]
        isoformScoreNovel <- isoformScore[isoformScore$Annotation == "Novel", ]
        isoformScoreAnno <- isoformScore[isoformScore$Annotation == "Annotated", ]

        isoformScoreNovel$Isos <- factor(isoformScoreNovel$Isos, levels=as.vector(isoformScoreNovel$Isos))
        isoformScoreAnno$Isos <- factor(isoformScoreAnno$Isos, levels=as.vector(isoformScoreAnno$Isos))
        pdf(outPdf)
        p1 <- ggplot(data=isoformScoreNovel, aes(x=Isos, y=Freq, group=1)) +
            geom_point() +
            geom_line() +
            theme_bw() + ggtitle("Novel Isoforms Rank") +
            theme(plot.title = element_text(hjust = 0.5, size=20), axis.text.x=element_blank(),
                axis.title.x=element_blank(), text = element_text(size=12), panel.border = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"), strip.background = element_blank(),
                strip.text.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8),
                legend.key.size = unit(0.5, "cm"), panel.spacing = unit(0, "lines"))
        print(p1)
        p2 <- ggplot(data=isoformScoreAnno, aes(x=Isos, y=Freq, group=1)) +
            geom_point() +
            geom_line() +
            theme_bw() + ggtitle("Annotated Isoforms Rank") +
            theme(plot.title = element_text(hjust = 0.5, size=20), axis.text.x=element_blank(),
                axis.title.x=element_blank(), text = element_text(size=12), panel.border = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"), strip.background = element_blank(),
                strip.text.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8),
                legend.key.size = unit(0.5, "cm"), panel.spacing = unit(0, "lines"))
        print(p2)
        dev.off()
    }
'''

PLOT_DIFF_AS_SOURCE = r'''
    library(ggplot2)
    plotDiffAS <- function(diffAsFile, outPdf){
        diffAS <- read.csv(diffAsFile, sep="\t")
        pdf(outPdf, height=6 ,width=10)
        p <- ggplot(diffAS, aes(x=AS_type, y=Count, fill=AS_type)) + geom_bar(stat = 'identity') +
            theme_bw() + xlab("AS types") +
            scale_y_continuous(expand = c(0.02, 0)) +
            theme(plot.title = element_text(hjust = 0.5, size=20),
                panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.text = element_text(size=20), axis.title = element_text(size=20),
                strip.background = element_blank(), strip.text.x = element_blank(), legend.title = element_blank(),
                legend.text = element_text(size = 16), legend.key.size = unit(0.75, "cm"),
                panel.spacing = unit(0, "lines"))
        print(p)
        dev.off()
    }
'''


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _utils_dir() -> Path:
    env_dir = os.environ.get("IFLAS_UTILS_DIR")
    candidates = [
        Path(env_dir).expanduser() if env_dir else None,
        _repo_root() / "utils",
        Path.cwd() / "utils",
    ]
    for path in candidates:
        if path and path.is_dir():
            return path
    raise FileNotFoundError("Cannot locate iFLAS utils directory (set IFLAS_UTILS_DIR)")


def _sample_base_dir(config: IflasConfig, sample: Sample) -> Path:
    return Path(config.dirs.out_dir).expanduser().resolve() / sample.project_name / sample.sample_name


def _sample_report_dir(config: IflasConfig, sample: Sample) -> Path:
    report_dir = (
        Path(config.dirs.out_dir).expanduser().resolve()
        / "reports"
        / f"{sample.project_name}_{sample.sample_name}"
    )
    report_dir.mkdir(parents=True, exist_ok=True)
    return report_dir


def _run_bash(command: str, *, cwd: Path, runner: ToolRunner, conda_env: str | None = None):
    return runner.run(["bash", "-c", command], cwd=str(cwd), conda_env=conda_env)


def _reset_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def _write_sorted_unique_lines(inputs: Iterable[Path], out_file: Path) -> int:
    lines = set()
    for in_file in inputs:
        if not in_file.is_file():
            continue
        with in_file.open(encoding="utf-8") as handle:
            for line in handle:
                line = line.rstrip("\n")
                if line:
                    lines.add(line)
    with out_file.open("w", encoding="utf-8") as handle:
        for line in sorted(lines):
            handle.write(line + "\n")
    return len(lines)


def _collect_gene_structure_pdfs(base_dir: Path) -> list[Path]:
    pdfs: list[Path] = []
    for root_name in ("visual_as", "isoViewer"):
        root = base_dir / root_name
        if not root.is_dir():
            continue
        pdfs.extend(sorted(path for path in root.glob("*/*.pdf") if path.is_file()))
    return pdfs


def select_report_samples(config: IflasConfig, sample_names: Iterable[str] | None = None) -> list:
    """Samples whose outputs a report reads. As in `iflas all` and the published report, with
    ``merge_data_from_same_strain`` the rows of a strain are reported once, from ``<project>/<strain>_all``; when only
    one of the two layouts exists on disk that one is used (BUG-160: the report used to look only for per-row
    directories, so after a default merged run every report came out empty). ``sample_names`` may name sample-sheet
    rows or merged ``<strain>_all`` samples."""
    try:
        groups = merge_samples_by_strain(config.samples)
    except ValueError:                       # rows that cannot be merged (e.g. NGS reads on only some rows): per row
        groups = {}
    out, seen, known = [], set(), set()
    for sample in config.samples:
        merged = groups.get(sample.project_name, {}).get(sample.ref_strain, {}).get(sample.strain)
        names = {sample.sample_name} | ({merged.sample_name} if merged is not None else set())
        known |= names
        if sample_names and not names & set(sample_names):
            continue
        candidates = [merged, sample] if config.options.merge_data_from_same_strain else [sample, merged]
        candidates = [c for c in candidates if c is not None]
        target = next((c for c in candidates if _sample_base_dir(config, c).is_dir()), sample)
        key = (target.project_name, target.sample_name)
        if key not in seen:
            seen.add(key)
            out.append(target)
    missing = sorted(set(sample_names or ()) - known)
    if missing:
        raise ValueError("unknown sample(s): " + ", ".join(missing))
    if not out:
        raise ValueError("report requires at least one sample")
    return out


def _per_sample_report(plot_fn, config: IflasConfig, sample_names: Iterable[str] | None, **kwargs) -> dict:
    """Run a per-sample plotter over the selected samples -> {"plots": all, "by_sample": {sample: plots}}."""
    by_sample = {sample.sample_name: plot_fn(config, sample, **kwargs) for sample in select_report_samples(config, sample_names)}
    return {"plots": [p for plots in by_sample.values() for p in plots], "by_sample": by_sample}


def plot_reads_corrected_eval_for_sample(
    config: IflasConfig,
    sample: Sample,
    *,
    runner=None,
) -> list[str]:
    """Plot the legacy corrected-vs-raw read accuracy panel for one sample."""

    base_dir = _sample_base_dir(config, sample)
    mapping_dir = base_dir / "mapping"
    raw_mapped = mapping_dir / "rawFlnc.addCVandID.bed12+"
    correct_mapped = mapping_dir / "flnc.addCVandID.bed12+"
    if not (raw_mapped.is_file() and correct_mapped.is_file()):
        return []

    basic_dir = _sample_report_dir(config, sample) / "basicStatistics"
    basic_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = "readsCorrectResult.pdf"
    run_r_function(
        PLOT_READS_CORRECTED_EVAL_SOURCE,
        "plotReadsCorrectedEval",
        [str(raw_mapped.resolve()), str(correct_mapped.resolve()), out_pdf],
        cwd=str(basic_dir),
        runner=runner,
    )
    return [str((basic_dir / out_pdf).resolve())]


def _gc_percent(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc * 100.0 / len(seq)


def _write_gc_across_read(fx_file: Path, out_file: Path, *, interval: int = 20, runner: ToolRunner) -> None:
    from Bio import SeqIO

    file_type = validateFaAndFqFile(str(fx_file), runner=runner)
    if not file_type:
        raise ValueError(f"cannot detect FASTA/FASTQ format: {fx_file}")

    with out_file.open("w", encoding="utf-8") as handle:
        for record in SeqIO.parse(str(fx_file), file_type):
            seq = str(record.seq)
            chunk_size = max(1, len(seq) // interval)           # legacy Py2 integer division, kept (BUG-104)
            values = [record.name]
            for start in range(0, len(seq), chunk_size):
                values.append(_gc_percent(seq[start:start + chunk_size]))
            handle.write("\t".join(map(str, values)) + "\n")


def _sample_flnc_file(config: IflasConfig, sample: Sample) -> Path | None:
    location = sample.data_processed_location                 # a list on merged samples: use the concatenated FLNC below
    if isinstance(location, str) and location:
        flnc = Path(location).expanduser().resolve()
        if flnc.is_file():
            return flnc

    tgs_plat = (sample.tgs_plat or "tgs").lower()
    fallback = _sample_base_dir(config, sample) / "preprocess" / tgs_plat / "rawFlnc.fq"
    return fallback if fallback.is_file() else None


def plot_reads_content_eval_for_sample(
    config: IflasConfig,
    sample: Sample,
    ref: RefSection,
    *,
    runner=None,
) -> list[str]:
    """Generate legacy basic read-content plots for one sample."""

    flnc = _sample_flnc_file(config, sample)
    if flnc is None:
        return []
    if not ref.ref_gpe:
        raise ValueError(f"sample {sample.sample_name!r}: report -basic requires ref_gpe")

    runner = runner or ToolRunner()
    util_dir = _utils_dir()
    basic_dir = _sample_report_dir(config, sample) / "basicStatistics"
    basic_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[str] = []

    _run_bash(
        f"seqkit fx2tab -n -g {q(flnc)} | cut -f 1,2 > GC_of_raw_flnc_reads.log",
        cwd=basic_dir,
        runner=runner,
    )
    _run_bash(
        "cut -f2 GC_of_raw_flnc_reads.log | "
        f"{q(util_dir / 'distrCurve.R')} -d -x='Binned GC%' -y='Fraction of reads' "
        "-v=50 -w=10 -mainS=28 -p=GC_of_raw_flnc_reads.pdf",
        cwd=basic_dir,
        runner=runner,
        conda_env="iflas-diffsplice",
    )
    outputs.append(str((basic_dir / "GC_of_raw_flnc_reads.pdf").resolve()))

    _write_gc_across_read(flnc, basic_dir / "GC_across_raw_flnc_read.log", runner=runner)
    _run_bash(
        "cut -f2- GC_across_raw_flnc_read.log | "
        f"{q(util_dir / 'box.R')} -stack -nJ -ho=50 -x=Interval -y=GC% -oS=0.5 "
        "-w=10 -height=6 -mainS=28 -p=GC_across_raw_flnc_read.pdf",
        cwd=basic_dir,
        runner=runner,
        conda_env="iflas-diffsplice",
    )
    outputs.append(str((basic_dir / "GC_across_raw_flnc_read.pdf").resolve()))

    _run_bash(
        f"seqkit fx2tab -n -l {q(flnc)} | cut -f 2 > readsLength.lst",
        cwd=basic_dir,
        runner=runner,
    )
    _run_bash(
        f"{q(util_dir / 'gpe2bed.pl')} {q(ref.ref_gpe)} | "
        f"{q(util_dir / 'bedLength.pl')} | cut -f 13 > refGeneLength.lst",
        cwd=basic_dir,
        runner=runner,
    )
    _run_bash(
        f"{q(util_dir / 'distrCurves.R')} -lgPosX=0.8 -lgPosY=0.8 -x1=0 -x2=10000 "
        "-d -x='Binned length (limited in 0-10000)' -w=10 *.lst -b=200 -mainS=28 "
        "-p=LengthDistribution.curve.pdf",
        cwd=basic_dir,
        runner=runner,
        conda_env="iflas-diffsplice",
    )
    outputs.append(str((basic_dir / "LengthDistribution.curve.pdf").resolve()))
    _run_bash(
        f"{q(util_dir / 'boxes.R')} -ng -no -xlab='Category' "
        "-ylab='Length of reads or transcripts' *.lst -w=10 -p=LengthDistribution.box.pdf",
        cwd=basic_dir,
        runner=runner,
        conda_env="iflas-diffsplice",
    )
    outputs.append(str((basic_dir / "LengthDistribution.box.pdf").resolve()))

    base_dir = _sample_base_dir(config, sample)
    junctions = Path(sample.ngs_junctions).expanduser().resolve() if sample.ngs_junctions else (
        base_dir / "mapping" / "rna-seq" / "reassembly" / "junctions.bed"
    )
    isoform_bed = base_dir / "refine" / "isoformGrouped.bed12+"
    if junctions.is_file() and isoform_bed.is_file():
        _run_bash(
            f"awk '$10>1' {q(isoform_bed)} | {q(util_dir / 'bed2gpe.pl')} | "
            f"{q(util_dir / 'transSupportByJunction.pl')} -j {q(junctions)} "
            "> supportedByRNAseq.tsv 2> supportedByRNAseq.summary",
            cwd=basic_dir,
            runner=runner,
        )
        _run_bash(
            "awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$4/$3}' supportedByRNAseq.summary | "
            f"{q(util_dir / 'summary2d.R')} -binWidthX=0.9999999 -binWidthY=0.9999999 "
            "-x='Junction count of long reads' -y='Supported Junction Count' "
            "-fL=gray -fH=red -mainS=28 -w=10 -p=supportedByRNAseq.pdf",
            cwd=basic_dir,
            runner=runner,
            conda_env="iflas-diffsplice",
        )
        outputs.append(str((basic_dir / "supportedByRNAseq.pdf").resolve()))

    return outputs


def plot_basic_for_sample(
    config: IflasConfig,
    sample: Sample,
    *,
    runner=None,
) -> list[str]:
    """Generate migrated ``report -basic`` plots for one sample."""

    runner = runner or ToolRunner()
    plots: list[str] = []
    if sample.use_fmlrc2:
        plots.extend(plot_reads_corrected_eval_for_sample(config, sample, runner=runner))
    plots.extend(
        plot_reads_content_eval_for_sample(
            config,
            sample,
            config.ref_for(sample.ref_strain),
            runner=runner,
        )
    )
    return plots


def plot_basic_report(config: IflasConfig, *, sample_names: Iterable[str] | None = None, runner=None) -> dict:
    """Generate migrated basic report plots for selected samples."""
    return _per_sample_report(plot_basic_for_sample, config, sample_names, runner=runner or ToolRunner())


def plot_novel_hq_as_for_sample(
    config: IflasConfig,
    sample: Sample,
    *,
    runner=None,
) -> list[str]:
    """Plot legacy novel/high-quality isoform score ranking for one sample."""

    out_dir = Path(config.dirs.out_dir).expanduser().resolve()
    score_file = out_dir / sample.project_name / sample.sample_name / "isoIdentity" / "isoformScore.txt"
    if not score_file.is_file():
        return []

    report_dir = out_dir / "reports" / f"{sample.project_name}_{sample.sample_name}"
    report_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = "isoformScore.pdf"
    run_r_function(
        PLOT_NOVEL_HQ_AS_SOURCE,
        "plotNovelHqAS",
        [str(score_file.resolve()), out_pdf],
        cwd=str(report_dir),
        runner=runner,
    )
    return [str((report_dir / out_pdf).resolve())]


def plot_novel_hq_as_report(config: IflasConfig, *, sample_names: Iterable[str] | None = None, runner=None) -> dict:
    """Generate migrated novel-HQ-AS report plots for selected samples."""
    return _per_sample_report(plot_novel_hq_as_for_sample, config, sample_names, runner=runner)


def _line_count(
    path: Path,
    *,
    sep: str = "\t",
    grep_str: str | None = None,
    col: int = 0,
    header: bool = False,
) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        return 0
    try:
        if header:
            df = pd.read_csv(path, sep=sep)
        else:
            df = pd.read_csv(path, sep=sep, header=None)
    except pd.errors.EmptyDataError:
        return 0
    if grep_str:
        df = df.loc[df.iloc[:, col] == grep_str]
    return len(df)


def _read_splice_site_table(path: Path, as_type: str, category: str) -> pd.DataFrame:
    columns = ["Dinucleotide", "Frequency"]
    if not path.is_file() or path.stat().st_size == 0:
        table = pd.DataFrame(columns=columns)
    else:
        try:
            table = pd.read_csv(
                path,
                sep="\t",
                header=None,
                names=columns,
            )
        except pd.errors.EmptyDataError:
            table = pd.DataFrame(columns=columns)
    table["AS_type"] = as_type
    table["Category"] = category
    return table


def plot_as_pattern_for_sample(
    config: IflasConfig,
    sample: Sample,
    *,
    runner=None,
) -> list[str]:
    """Generate legacy AS-pattern count and splice-site plots for one sample."""

    base_dir = _sample_base_dir(config, sample)
    charac_as_dir = base_dir / "as_events" / "characterization"
    if not charac_as_dir.is_dir():
        return []

    as_pattern_dir = _sample_report_dir(config, sample) / "asPattern"
    as_pattern_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{sample.project_name}_{sample.sample_name}"
    as_list_inputs = [
        charac_as_dir / f"{as_type}.{annotation}.lst"
        for as_type in AS_TYPES
        for annotation in ("anno", "novel")
    ]
    if not any(path.is_file() for path in as_list_inputs):
        return []

    count_rows = [
        ("IR", "Anno", _line_count(charac_as_dir / "IR.anno.lst")),
        ("IR", "Novel", _line_count(charac_as_dir / "IR.novel.lst")),
        ("SE", "Anno", _line_count(charac_as_dir / "SE.anno.lst")),
        ("SE", "Novel", _line_count(charac_as_dir / "SE.novel.lst")),
        ("A3SS", "Anno", _line_count(charac_as_dir / "A3SS.anno.lst")),
        ("A3SS", "Novel", _line_count(charac_as_dir / "A3SS.novel.lst")),
        ("A5SS", "Anno", _line_count(charac_as_dir / "A5SS.anno.lst")),
        ("A5SS", "Novel", _line_count(charac_as_dir / "A5SS.novel.lst")),
        ("APA", "Anno", _line_count(charac_as_dir / "annoPA.tsv")),
        ("APA", "Novel", _line_count(charac_as_dir / "novelPA.tsv")),
        ("PA", "Anno", _line_count(charac_as_dir / "statistic.APA.tsv", grep_str="Known", col=2)),
        ("PA", "Novel", _line_count(charac_as_dir / "statistic.APA.tsv", grep_str="Novel", col=2)),
    ]
    as_anno_file = as_pattern_dir / f"{stem}.AS_annotation.txt"
    with as_anno_file.open("w", encoding="utf-8") as handle:
        handle.write("AS_type\tAnnotation\tCount\n")
        for as_type, annotation, count in count_rows:
            handle.write(f"{as_type}\t{annotation}\t{count}\n")

    splice_site = pd.concat(
        [
            _read_splice_site_table(charac_as_dir / "IR.tmp.txt", "IR", "Inc"),
            _read_splice_site_table(charac_as_dir / "SE.inc.tmp.txt", "SE", "Inc"),
            _read_splice_site_table(charac_as_dir / "SE.exc.tmp.txt", "SE", "Exc"),
            _read_splice_site_table(charac_as_dir / "A3SS.inc.tmp.txt", "A3SS", "Inc"),
            _read_splice_site_table(charac_as_dir / "A3SS.exc.tmp.txt", "A3SS", "Exc"),
            _read_splice_site_table(charac_as_dir / "A5SS.inc.tmp.txt", "A5SS", "Inc"),
            _read_splice_site_table(charac_as_dir / "A5SS.exc.tmp.txt", "A5SS", "Exc"),
        ]
    )
    as_splice_site_file = as_pattern_dir / f"{stem}.AS_spliceSite.txt"
    splice_site.to_csv(as_splice_site_file, sep="\t", header=True, index=False)

    as_anno_pdf = f"{stem}.AS_annotation.pdf"
    run_r_function(
        legacy_r.plotAsCountStatisticsStr,
        "plotAsCountStatistics",
        [as_anno_file.name, as_anno_pdf],
        cwd=str(as_pattern_dir),
        runner=runner,
    )
    plots = [str((as_pattern_dir / as_anno_pdf).resolve())]
    if not splice_site.empty:
        as_splice_site_pdf = f"{stem}.AS_spliceSite.pdf"
        run_r_function(
            legacy_r.plotAsDinucleotideStatisticsStr,
            "plotAsDinucleotideStatistics",
            [as_splice_site_file.name, as_splice_site_pdf],
            cwd=str(as_pattern_dir),
            runner=runner,
        )
        plots.append(str((as_pattern_dir / as_splice_site_pdf).resolve()))
    else:
        as_splice_site_file.unlink(missing_ok=True)
    return plots


def plot_as_pattern_report(config: IflasConfig, *, sample_names: Iterable[str] | None = None, runner=None) -> dict:
    """Generate migrated AS-pattern report plots for selected samples."""
    return _per_sample_report(plot_as_pattern_for_sample, config, sample_names, runner=runner)


def plot_target_gene_structure_for_sample(
    config: IflasConfig,
    sample: Sample,
) -> list[str]:
    """Merge existing visual_as/legacy isoViewer gene-structure PDFs for one sample."""

    base_dir = _sample_base_dir(config, sample)
    gene_pdfs = _collect_gene_structure_pdfs(base_dir)
    if not gene_pdfs:
        return []

    gene_struc_dir = _sample_report_dir(config, sample) / "geneStrucPlots"
    gene_struc_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = gene_struc_dir / "allTargetGeneStructure.pdf"
    _merge_pdfs(gene_pdfs, out_pdf)
    return [str(out_pdf.resolve())]


def plot_target_gene_structure_report(config: IflasConfig, *, sample_names: Iterable[str] | None = None) -> dict:
    """Generate migrated target-gene-structure PDF merges for selected samples."""
    return _per_sample_report(plot_target_gene_structure_for_sample, config, sample_names)


def plot_diff_as_for_comparison(
    sig_diff_as_dir: str | Path,
    *,
    runner=None,
) -> list[str]:
    """Plot the legacy differential-AS summary for one ``*.sigDiffAS`` directory."""

    comp_dir = Path(sig_diff_as_dir).expanduser().resolve()
    if not comp_dir.is_dir():
        raise FileNotFoundError(str(comp_dir))

    rows = []
    for as_type in AS_TYPES:
        sig_file = comp_dir / f"{as_type}.sig.txt"
        if not sig_file.is_file():
            raise FileNotFoundError(str(sig_file))
        rows.append((as_type, len(pd.read_csv(sig_file, sep="\t"))))

    summary_file = comp_dir / "sigDiff.AS.txt"
    with summary_file.open("w", encoding="utf-8") as handle:
        handle.write("AS_type\tCount\n")
        for as_type, count in rows:
            handle.write(f"{as_type}\t{count}\n")

    out_pdf = "sigDiff.AS_distribution.pdf"
    run_r_function(
        PLOT_DIFF_AS_SOURCE,
        "plotDiffAS",
        [summary_file.name, out_pdf],
        cwd=str(comp_dir),
        runner=runner,
    )
    return [str((comp_dir / out_pdf).resolve())]


def plot_diff_as_report(
    config: IflasConfig,
    *,
    runner=None,
) -> dict:
    """Generate migrated differential-AS summary plots from existing rMATS outputs."""

    das_dir = Path(config.dirs.out_dir).expanduser().resolve() / "das"
    comparisons = sorted(path for path in das_dir.glob("*.sigDiffAS") if path.is_dir())
    by_comparison = {}
    plots = []
    for comp_dir in comparisons:
        comp_plots = plot_diff_as_for_comparison(comp_dir, runner=runner)
        by_comparison[comp_dir.name] = comp_plots
        plots.extend(comp_plots)
    return {"plots": plots, "by_comparison": by_comparison}


def _plot_pa_tail_apa_source(ref_gtf: str | Path) -> str:
    make_txdb = '''
        options(ucscChromosomeNames=FALSE)
        refGtfDb <- txdbmaker::makeTxDbFromGFF(%s, format = "gtf")   # defunct in GenomicFeatures >= 1.61.1 (BUG-164)
        refGtfTrack <- GeneRegionTrack(refGtfDb, name="Gene model", transcriptAnnotation = "transcript", stackHeight = 0.5)
    ''' % (r_literal(str(Path(ref_gtf).expanduser().resolve())))
    return legacy_r.plotPaTailASStr + "\n" + make_txdb


def plot_palen_for_sample(
    config: IflasConfig,
    sample: Sample,
    *,
    runner=None,
) -> list[str]:
    """Generate migrated poly(A)-tail AS/APA report plots for one sample."""

    base_dir = _sample_base_dir(config, sample)
    palen_dir = base_dir / "palenAS"
    if not palen_dir.is_dir():
        return []

    tool_runner = runner or ToolRunner()
    report_dir = _sample_report_dir(config, sample)
    plots: list[str] = []

    palen_as_plots_dir = report_dir / "palenAsPlots"
    _reset_dir(palen_as_plots_dir)
    palen_as_inputs = [
        palen_dir / "mergeByJunc" / f"{as_type}.palenAndAS.sig.bed12+"
        for as_type in AS_TYPES
    ]
    merged_as = palen_as_plots_dir / "palenAndAS.sig.bed12+"
    if _write_sorted_unique_lines(palen_as_inputs, merged_as):
        run_r_function(
            legacy_r.plotPaTailASStr,
            "plotPaTailAsStructure",
            [merged_as.name],
            cwd=str(palen_as_plots_dir),
            runner=tool_runner,
        )
        plots.extend(str(path.resolve()) for path in sorted(palen_as_plots_dir.glob("*.palenAndAS.pdf")))

    palen_apa_file = palen_dir / "palenAPA" / "apaRelatedPalen.txt"
    read_bed = base_dir / "refine" / "reads.assigned.unambi.bed12+"
    align_bam = base_dir / "mapping" / "flnc.mm2.sorted.bam"
    if palen_apa_file.is_file() and read_bed.is_file():
        if not align_bam.is_file():
            raise FileNotFoundError(str(align_bam))
        index_bam(str(align_bam.resolve()), runner=tool_runner)
        needs_csi = longest_reference(align_bam) > BAI_MAX_LENGTH
        gene_pos = BedFile(str(read_bed), type="bed12+").getGenePos(bedType="bed12+", geneCol=15)
        palen_apa_plots_dir = report_dir / "palenApaPlots"
        _reset_dir(palen_apa_plots_dir)
        apa_source = _plot_pa_tail_apa_source(config.ref_for(sample.ref_strain).ref_gtf)
        apa_calls = []
        with palen_apa_file.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                info = line.rstrip("\n").split("\t")
                if len(info) < 5:
                    raise ValueError(f"malformed palen APA line in {palen_apa_file}: {line.rstrip()}")
                gene_name, _strand, pa_site, apa_support, pa_len_list = info[:5]
                if gene_name not in gene_pos:
                    raise ValueError(f"{read_bed}: cannot locate gene {gene_name!r} for palen APA plot")
                chrom, chrom_start, chrom_end = gene_pos[gene_name]
                if int(chrom_end) > RSAMTOOLS_MAX_END:
                    _skip_beyond_rsamtools(f"palen APA plot of {gene_name}", chrom, chrom_end)
                    continue
                gene_bam = align_bam.resolve()
                if needs_csi:            # Rsamtools only finds a .bai: hand R this gene's reads with one
                    gene_bam = palen_apa_plots_dir / f"{gene_name}.region.bam"
                    _run_bash(f"samtools view -b {q(align_bam.resolve())} {chrom}:{int(chrom_start) + 1}-{chrom_end} "
                              f"> {q(gene_bam)} && samtools index {q(gene_bam)}", cwd=palen_apa_plots_dir, runner=tool_runner)
                apa_calls.append([gene_name, pa_site, pa_len_list, apa_support, str(gene_bam), chrom, int(chrom_start), int(chrom_end)])
        if apa_calls:                   # one R session: the reference TxDb is built once, not once per gene
            run_r_calls(apa_source, "plotPaTailApaStructure", apa_calls, cwd=str(palen_apa_plots_dir), runner=tool_runner)
        plots.extend(str(path.resolve()) for path in sorted(palen_apa_plots_dir.glob("*.palenAndAPA.pdf")))

    return plots


def plot_palen_report(config: IflasConfig, *, sample_names: Iterable[str] | None = None, runner=None) -> dict:
    """Generate migrated poly(A)-tail AS/APA report plots for selected samples."""
    return _per_sample_report(plot_palen_for_sample, config, sample_names, runner=runner or ToolRunner())


def _asas_partial_file(base_dir: Path) -> Path | None:
    """Return the migrated ASAS result file, keeping the legacy alleleAS fallback."""

    for candidate in (
        base_dir / "asas" / "partialAsRelatedHaplotype.txt",
        base_dir / "alleleAS" / "partialAsRelatedHaplotype.txt",
    ):
        if candidate.is_file():
            return candidate
    return None


def _skip_beyond_rsamtools(what, chrom, end) -> None:
    get_logger().warning("%s skipped: it reaches %s:%s, and Bioconductor's Rsamtools cannot read alignments beyond "
                         "%d bp (chromosomes of wheat, barley, ...; BUG-189)", what, chrom, end, RSAMTOOLS_MAX_END)


def _parse_asas_highlight(as_type: str, as_event: str, iso_names: list[str], iso_bed: BedFile) -> tuple[int, int]:
    try:
        if as_type == "SE" and "@" in as_event:
            skipped = as_event.split("@")[1]
            parts = skipped.split("-")
            return int(parts[0]), int(parts[-1])
        coord = as_event.split(":")[-1]
        start, end = coord.split("-")[:2]
        return int(start), int(end)
    except (IndexError, ValueError):
        exons = []
        for iso_name in iso_names:
            exons.extend(iso_bed.reads[iso_name].exons)
        if not exons:
            raise ValueError(f"cannot infer ASAS highlight interval for {as_event!r}")
        return min(start for start, _end in exons), max(end for _start, end in exons)


def _write_asas_iso_bed(out_file: Path, iso_names: list[str], iso_bed: BedFile) -> None:
    with out_file.open("w", encoding="utf-8") as handle:
        for iso_name in iso_names:
            if iso_name not in iso_bed.reads:
                raise KeyError(f"{iso_bed.bedFile}: missing ASAS isoform {iso_name!r}")
            handle.write(str(iso_bed.reads[iso_name]) + "\n")


def _write_asas_read_names(out_file: Path, iso_names: list[str], isoform_to_reads: dict) -> None:
    reads = []
    for iso_name in iso_names:
        if iso_name not in isoform_to_reads:
            raise KeyError(f"collapse group file: missing reads for ASAS isoform {iso_name!r}")
        reads.extend(isoform_to_reads[iso_name])
    with out_file.open("w", encoding="utf-8") as handle:
        for read_name in dict.fromkeys(reads):
            handle.write(read_name + "\n")


def _merge_pdfs(pdf_files: Iterable[Path], out_pdf: Path) -> None:
    from pypdf import PdfReader, PdfWriter

    writer = PdfWriter()
    for pdf_file in pdf_files:
        with pdf_file.open("rb") as handle:
            reader = PdfReader(handle)
            for page in reader.pages:
                writer.add_page(page)
    with out_pdf.open("wb") as handle:
        writer.write(handle)


def plot_asas_for_sample(
    config: IflasConfig,
    sample: Sample,
    *,
    runner=None,
) -> list[str]:
    """Generate migrated allele-specific AS structure plots for one sample."""

    base_dir = _sample_base_dir(config, sample)
    partial_file = _asas_partial_file(base_dir)
    if partial_file is None:
        return []

    isoform_file = base_dir / "hqIsoforms" / "hq.collapsed.bed12+"
    group_file = base_dir / "collapse" / "tofu.collapsed.group.txt"
    flnc_bam = base_dir / "mapping" / "flnc.mm2.sorted.bam"
    for required in (isoform_file, group_file, flnc_bam):
        if not required.is_file():
            raise FileNotFoundError(str(required))

    tool_runner = runner or ToolRunner()
    util_dir = _utils_dir()
    ref_genome = Path(config.ref_for(sample.ref_strain).ref_genome).expanduser().resolve()
    iso_bed = BedFile(str(isoform_file), type="bed12+")
    isoform_to_reads = getDictFromFile(str(group_file), sep="\t", inlineSep=",", valueCol=2)
    report_dir = _sample_report_dir(config, sample)
    plots_dir = report_dir / "alleleAsPlots"
    _reset_dir(plots_dir)

    as_haplo = pd.read_csv(
        partial_file,
        header=None,
        sep="\t",
        names=["gene", "asType", "refGene", "asEvent", "haplo1", "haplo1isos", "haplo2", "haplo2isos"],
    )
    if as_haplo.empty:
        return []

    event_pdfs: list[Path] = []
    r_calls = []
    for (gene, ref_gene, as_type), group in as_haplo.groupby(["gene", "refGene", "asType"], sort=True):
        out_name = f"{ref_gene}.{gene}.{as_type}.allele_as"      # BUG-148: clusters sharing a refGene must not collide
        event_dir = plots_dir / out_name
        _reset_dir(event_dir)

        haplo1_iso_names: list[str] = []
        haplo2_iso_names: list[str] = []
        chrom_starts: list[int] = []
        chrom_ends: list[int] = []
        chroms = set()
        highlight_starts: list[int] = []
        highlight_ends: list[int] = []

        for row in group.itertuples(index=False):
            row_haplo1 = str(row.haplo1isos).split("_")
            row_haplo2 = str(row.haplo2isos).split("_")
            haplo1_iso_names.extend(row_haplo1)
            haplo2_iso_names.extend(row_haplo2)
            highlight_start, highlight_end = _parse_asas_highlight(
                as_type,
                str(row.asEvent),
                row_haplo1 + row_haplo2,
                iso_bed,
            )
            highlight_starts.append(highlight_start)
            highlight_ends.append(highlight_end)

        all_iso_names = haplo1_iso_names + haplo2_iso_names
        for iso_name in all_iso_names:
            if iso_name not in iso_bed.reads:
                raise KeyError(f"{isoform_file}: missing ASAS isoform {iso_name!r}")
            iso = iso_bed.reads[iso_name]
            chrom_starts.append(iso.chromStart)
            chrom_ends.append(iso.chromEnd)
            chroms.add(iso.chrom)
        chrom = sorted(chroms)[0]
        chrom_start, chrom_end = min(chrom_starts), max(chrom_ends)
        if chrom_end > RSAMTOOLS_MAX_END:
            _skip_beyond_rsamtools(f"allele-specific AS plot {out_name}", chrom, chrom_end)
            continue

        _write_asas_iso_bed(event_dir / "haplo1isosOut.bed", haplo1_iso_names, iso_bed)
        _write_asas_iso_bed(event_dir / "haplo2isosOut.bed", haplo2_iso_names, iso_bed)
        _write_asas_read_names(event_dir / "haplo1.reads.txt", haplo1_iso_names, isoform_to_reads)
        _write_asas_read_names(event_dir / "haplo2.reads.txt", haplo2_iso_names, isoform_to_reads)

        _run_bash(
            f"{q(util_dir / 'bed2gpe.pl')} -g 13 haplo1isosOut.bed > haplo1isosOut.gpe; "
            "genePredToGtf file haplo1isosOut.gpe haplo1isosOut.gtf; "
            f"{q(util_dir / 'bed2gpe.pl')} -g 13 haplo2isosOut.bed > haplo2isosOut.gpe; "
            "genePredToGtf file haplo2isosOut.gpe haplo2isosOut.gtf",
            cwd=event_dir,
            runner=tool_runner,
        )
        threads = max(1, int(getattr(sample, "single_run_threads", None) or config.options.threads or 1))
        _run_bash(
            f"samtools view -b -N haplo1.reads.txt {q(flnc_bam.resolve())} | "
            f"samtools sort -@ {threads} -o haplo1.flnc.sorted.bam -; "
            f"samtools view -b -N haplo2.reads.txt {q(flnc_bam.resolve())} | "
            f"samtools sort -@ {threads} -o haplo2.flnc.sorted.bam -; "
            "samtools index haplo1.flnc.sorted.bam; "
            "samtools index haplo2.flnc.sorted.bam",
            cwd=event_dir,
            runner=tool_runner,
        )
        _run_bash(
            f"samtools cat haplo1.flnc.sorted.bam haplo2.flnc.sorted.bam | "
            f"samtools sort -@ {threads} -o {out_name}.flnc.sorted.bam -; "
            f"samtools index {out_name}.flnc.sorted.bam",
            cwd=event_dir,
            runner=tool_runner,
        )

        target_ngs_bam = ""
        ngs_bam = base_dir / "mapping" / "rna-seq" / "reassembly" / "tmp.bam"
        if (sample.ngs_left_reads or sample.ngs_right_reads) and ngs_bam.is_file():
            _run_bash(
                f"samtools view -h {q(ngs_bam.resolve())} {chrom}:{chrom_start}-{chrom_end} | "
                "samtools sort -o ngsReads.sorted.bam -; "
                "samtools index ngsReads.sorted.bam",
                cwd=event_dir,
                runner=tool_runner,
            )
            target_ngs_bam = str((event_dir / "ngsReads.sorted.bam").resolve())

        out_base = plots_dir / out_name
        r_calls.append([
            str(ref_genome),
            f"{event_dir / 'haplo1isosOut.gtf'},{event_dir / 'haplo2isosOut.gtf'}",
            str(event_dir / f"{out_name}.flnc.sorted.bam"),
            f"{event_dir / 'haplo1.flnc.sorted.bam'},{event_dir / 'haplo2.flnc.sorted.bam'}",
            target_ngs_bam,
            chrom,
            chrom_start,
            chrom_end,
            ",".join(map(str, highlight_starts)),
            ",".join(map(str, highlight_ends)),
            ASAS_HIGHLIGHT_COLORS[as_type],
            str(out_base),
        ])
        event_pdfs.append(out_base.parent / f"{out_base.name}.pdf")

    if r_calls:                         # one R session for every event (Gviz loads once)
        run_r_calls(legacy_r.plotAlleleAsStructureStr, "plotAlleleAsStructure", r_calls, cwd=str(plots_dir), runner=tool_runner)
    aggregate_pdf = report_dir / "alleleAS.pdf"
    _merge_pdfs(event_pdfs, aggregate_pdf)
    return [str(path.resolve()) for path in event_pdfs] + [str(aggregate_pdf.resolve())]


def plot_asas_report(config: IflasConfig, *, sample_names: Iterable[str] | None = None, runner=None) -> dict:
    """Generate migrated allele-specific AS report plots for selected samples."""
    return _per_sample_report(plot_asas_for_sample, config, sample_names, runner=runner or ToolRunner())
