"""End-to-end acceptance test for the orchestrated core chain on synthetic data.

This covers the migrated core chain through the modern config/orchestrator path:
examples-style YAML/TSV -> init-ref -> mapping -> collapse -> refine -> pu_filter
-> find_as. The test is tool-gated because it exercises the real external
programs and the isolated cupcake environment.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys

import pytest
import yaml

from iflas.cli import main as cli_main
from iflas.config import load_config
from iflas.io.genepred import GenePredObj
from iflas.orchestrator import run_pipeline

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_synthetic_data as syn  # noqa: E402


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXAMPLES = os.path.join(ROOT, "examples")
SYNTHETIC = os.path.join(ROOT, "test_data", "synthetic")

_TOOLS_FULL = all(
    shutil.which(t)
    for t in (
        "gtfToGenePred",
        "bedtools",
        "minimap2",
        "hisat2",
        "hisat2-build",
        "stringtie",
        "regtools",
        "samtools",
        "seqkit",
        "perl",
    )
)


def _cupcake_available():
    if shutil.which("conda") is None:
        return False
    r = subprocess.run(
        ["conda", "run", "-n", "iflas-cupcake", "collapse_isoforms_by_sam.py", "--help"],
        capture_output=True,
    )
    return r.returncode == 0


pytestmark = pytest.mark.skipif(
    not (_TOOLS_FULL and _cupcake_available()),
    reason="needs full core toolchain + iflas-cupcake env",
)


class _Future:
    def __init__(self, value):
        self._value = value

    def result(self):
        return self._value


class _SerialExecutor:
    """Run submitted stage jobs inline; keeps this tool-heavy e2e deterministic."""

    def __init__(self, max_workers=None):
        self.max_workers = max_workers

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def submit(self, fn, *args, **kwargs):
        return _Future(fn(*args, **kwargs))


def _rows(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f if line.strip()]


def _write_subset_bed(ref_dir):
    """Canonical-isoform ref_bed for the PU annotated/novel split.

    ``init-ref`` still derives the main aux files from the full synthetic GTF.
    The PU learner needs both annotated and novel rows; the synthetic PU tests
    use this same canonical-only reference bed to keep that class contract.
    """
    subset_gtf = os.path.join(ref_dir, "genes.canonical.gtf")
    with open(subset_gtf, "w") as f:
        for gene, isos in syn.GENES.items():
            for name, chrom, strand, exons in isos:
                if not name.endswith(".iso1"):
                    continue
                for start, end in exons:
                    f.write(
                        f'{chrom}\tsyn\texon\t{start + 1}\t{end}\t.\t{strand}\t.\t'
                        f'gene_id "{gene}"; transcript_id "{name}";\n'
                    )

    subset_gpe = os.path.join(ref_dir, "genes.canonical.gpe")
    subset_bed = os.path.join(ref_dir, "genes.canonical.bed12+")
    subprocess.run(["gtfToGenePred", "-genePredExt", subset_gtf, subset_gpe], check=True)
    with open(subset_bed, "w") as out:
        for line in GenePredObj(subset_gpe, bincolumn=False).toBed(gene=True):
            out.write(str(line) + "\n")
    return subset_bed


def _write_examples_style_config(tmp_path, data_dir, ref_dir):
    config_path = tmp_path / "config.yaml"
    samples_path = tmp_path / "samples.tsv"

    assert cli_main(
        [
            "init-ref",
            "-f",
            str(data_dir / "genome.fa"),
            "--gtf",
            str(data_dir / "genes.gtf"),
            "-o",
            str(ref_dir),
            "--hisat2",
        ]
    ) == 0

    subset_bed = _write_subset_bed(str(ref_dir))
    with open(os.path.join(EXAMPLES, "config.yaml")) as f:
        cfg = yaml.safe_load(f)
    cfg["reference"][0].update(
        {
            "ref_genome": str(data_dir / "genome.fa"),
            "ref_gtf": str(data_dir / "genes.gtf"),
            "ref_size": str(ref_dir / "genome.len.txt"),
            "ref_gpe": str(ref_dir / "genome.gpe"),
            "ref_bed": subset_bed,
            "ref_mm2_index": str(ref_dir / "genome.mm2"),
            "hisat2_index": str(ref_dir / "genome_hisat2"),
        }
    )
    cfg["dirs"]["out_dir"] = str(tmp_path / "iflas_out")
    cfg["options"]["threads"] = 2
    cfg["pu_filter"]["filter_score"] = 0.0
    with open(config_path, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)

    with open(os.path.join(EXAMPLES, "samples.tsv"), newline="") as f:
        row = next(csv.DictReader(f, delimiter="\t"))
    row.update(
        {
            "sample_name": "synthetic",
            "project_name": "demo",
            "ref_strain": "B73",
            "strain": "B73",
            "condition": "synthetic",
            "replicate": "rep1",
            "tgs_plat": "pacbio",
            "strategy": "sequel",
            "data_location": "None",
            "data_processed_location": str(data_dir / "flnc.fa"),
            "primer": str(data_dir / "primer.fa"),
            "ngs_left_reads": str(data_dir / "ngs_R1.fastq.gz"),
            "ngs_right_reads": str(data_dir / "ngs_R2.fastq.gz"),
            "ngs_reads_paired": "paired",
            "ngs_reads_length": "100",
            "use_fmlrc2": "false",
        }
    )
    fieldnames = list(row)
    if "data_processed_location" not in fieldnames:
        fieldnames.append("data_processed_location")
    with open(samples_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    return config_path, samples_path


def _assert_core_chain_outputs(base):
    expected_files = [
        base / "mapping" / "flnc.mm2.sorted.bam",
        base / "mapping" / "flnc.mm2.sorted.bam.bai",
        base / "mapping" / "flnc.addCVandID.bed12+",
        base / "mapping" / "flnc.processed.fa",
        base / "mapping" / "rna-seq" / "reassembly" / "tmp.bam",
        base / "mapping" / "rna-seq" / "reassembly" / "junctions.bed",
        base / "collapse" / "tofu.collapsed.good.gff",
        base / "collapse" / "tofu.collapsed.group.txt",
        base / "refine" / "isoformGrouped.bed12+",
        base / "refine" / "tofu.collapsed.assigned.unambi.bed12+",
        base / "refine" / "reads.assigned.unambi.bed12+",
        base / "refine" / "reads.unambi.gtf",
        base / "hqIsoforms" / "isoFeatures.txt",
        base / "hqIsoforms" / "pu_score.txt",
        base / "hqIsoforms" / "validIsoforms.lst",
        base / "hqIsoforms" / "hq.collapsed.bed12+",
    ]
    for path in expected_files:
        assert path.is_file() and path.stat().st_size > 0, path

    hq_rows = _rows(base / "hqIsoforms" / "hq.collapsed.bed12+")
    assert hq_rows and all(len(row) == 19 for row in hq_rows)

    lr_dir = base / "as_events" / "ordinary_as" / "LR"
    for name in (
        "IR.confident.bed6+",
        "SE.confident.bed12+",
        "A5SS.confident.bed6+",
        "A3SS.confident.bed6+",
    ):
        path = lr_dir / name
        assert path.is_file() or path.is_symlink(), path

    assert "geneB:2200@2500-2650@3000" in [
        row[3] for row in _rows(lr_dir / "SE.bed12+")
    ]


def test_run_pipeline_synthetic_core_chain_outputs_acceptance_files(tmp_path):
    data_dir = tmp_path / "synthetic"
    ref_dir = tmp_path / "ref"
    shutil.copytree(SYNTHETIC, data_dir)
    ref_dir.mkdir()

    config_path, samples_path = _write_examples_style_config(tmp_path, data_dir, ref_dir)

    cfg = load_config(str(config_path), str(samples_path))
    results = run_pipeline(
        cfg,
        stages=("mapping", "collapse", "refine", "pu_filter", "find_as"),
        max_workers=1,
        executor_factory=_SerialExecutor,
    )
    assert [result.stage for result in results] == [
        "mapping",
        "collapse",
        "refine",
        "pu_filter",
        "find_as",
    ]

    base = tmp_path / "iflas_out" / "demo" / "synthetic"
    _assert_core_chain_outputs(base)


def test_cli_all_synthetic_core_chain_returns_with_process_executor(tmp_path):
    data_dir = tmp_path / "synthetic"
    ref_dir = tmp_path / "ref"
    shutil.copytree(SYNTHETIC, data_dir)
    ref_dir.mkdir()

    config_path, samples_path = _write_examples_style_config(tmp_path, data_dir, ref_dir)

    result = subprocess.run(
        [
            "iflas",
            "all",
            "-c",
            str(config_path),
            "-s",
            str(samples_path),
            "--stages",
            "mapping,collapse,refine,pu_filter,find_as",
            "-j",
            "1",
        ],
        capture_output=True,
        text=True,
        timeout=240,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Ran 5 stage/sample jobs" in result.stdout
    assert "pu_filter: B73_all" in result.stdout
    _assert_core_chain_outputs(tmp_path / "iflas_out" / "demo" / "B73_all")
