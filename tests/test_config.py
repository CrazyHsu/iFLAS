"""Unit tests for iflas.config (YAML + TSV sample sheet + pydantic validation)."""

import os

import pytest
from pydantic import ValidationError

from iflas.config import MergedTgsSample, Sample, load_config, merge_samples_by_strain, resolve_path

CONFIG_YAML = """\
reference:
  - ref_strain: B73
    ref_genome: genome.fa
    ref_gtf: genes.gtf
ccs:
  min_pass: 3
options:
  threads: 8
  gene2go: gene2go.txt
dirs:
  out_dir: out
"""

SAMPLES_TSV = (
    "sample_name\tref_strain\tcondition\treplicate\tngs_left_reads\tngs_reads_length\tuse_fmlrc2\n"
    "s1\tB73\troot\trep1\tR1.fastq.gz\t150\ttrue\n"
    "s2\tB73\tleaf\trep1\tNone\t\tfalse\n"
)


def _write(tmp_path, yaml_text=CONFIG_YAML, tsv_text=SAMPLES_TSV):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(yaml_text)
    tsv = tmp_path / "samples.tsv"
    tsv.write_text(tsv_text)
    return str(cfg), str(tsv)


def test_resolve_path():
    assert resolve_path(None, "/base") is None
    assert resolve_path("rel.fa", "/base") == "/base/rel.fa"
    assert resolve_path("/abs.fa", "/base") == "/abs.fa"
    os.environ["IFLAS_TEST_DIR"] = "/data/x"
    assert resolve_path("${IFLAS_TEST_DIR}/g.fa", "/base") == "/data/x/g.fa"


def test_load_defaults_and_overrides(tmp_path):
    cfg_path, tsv_path = _write(tmp_path)
    cfg = load_config(cfg_path, tsv_path)
    # reference + relative-path resolution
    assert cfg.reference[0].ref_strain == "B73"
    assert cfg.reference[0].ref_genome == str(tmp_path / "genome.fa")
    assert cfg.ref_for("B73").ref_gtf == str(tmp_path / "genes.gtf")
    # overrides + defaults
    assert cfg.ccs.min_pass == 3
    assert cfg.ccs.min_read_score == 0.75
    assert cfg.options.threads == 8
    assert cfg.options.gene2go == str(tmp_path / "gene2go.txt")
    assert cfg.collapse.min_identity == 0.9   # whole section defaulted


def test_pu_filter_feature_file_path_is_resolved(tmp_path):
    yaml_text = CONFIG_YAML + "pu_filter:\n  feature_file: features/isoFeatures.txt\n"
    cfg_path, tsv_path = _write(tmp_path, yaml_text=yaml_text)

    cfg = load_config(cfg_path, tsv_path)

    assert cfg.pu_filter.feature_file == str(tmp_path / "features" / "isoFeatures.txt")


def test_samples_parsing(tmp_path):
    cfg_path, tsv_path = _write(tmp_path)
    cfg = load_config(cfg_path, tsv_path)
    assert [s.sample_name for s in cfg.samples] == ["s1", "s2"]
    s1, s2 = cfg.samples
    assert s1.use_fmlrc2 is True and s1.ngs_reads_length == 150
    assert s1.ngs_left_reads == str(tmp_path / "R1.fastq.gz")
    assert s2.use_fmlrc2 is False
    assert s2.ngs_reads_length is None      # empty -> None
    assert s2.ngs_left_reads is None        # "None" token -> None


def test_sample_boolean_yes_no_tokens_are_parsed(tmp_path):
    tsv = (
        "sample_name\tref_strain\tuse_fmlrc2\n"
        "s1\tB73\tyes\n"
        "s2\tB73\tno\n"
    )
    cfg_path, tsv_path = _write(tmp_path, tsv_text=tsv)

    cfg = load_config(cfg_path, tsv_path)

    assert cfg.samples[0].use_fmlrc2 is True
    assert cfg.samples[1].use_fmlrc2 is False


def test_invalid_sample_boolean_is_rejected(tmp_path):
    tsv = "sample_name\tref_strain\tuse_fmlrc2\ns1\tB73\tmaybe\n"
    cfg_path, tsv_path = _write(tmp_path, tsv_text=tsv)

    with pytest.raises(ValidationError):
        load_config(cfg_path, tsv_path)


def test_gff_only_reference_points_to_init_ref(tmp_path):
    """Legacy configs could name only ref_gff (converted at start-up); the error now says how to derive the GTF."""
    gff_yaml = CONFIG_YAML.replace("    ref_gtf: genes.gtf\n", "    ref_gff: genes.gff3\n")
    cfg_path, tsv_path = _write(tmp_path, yaml_text=gff_yaml)
    with pytest.raises(ValidationError, match="iflas init-ref"):
        load_config(cfg_path, tsv_path)


def test_compressed_reference_points_to_init_ref(tmp_path):
    """BUG-194: a gzipped genome or GTF named directly in the config stops at load time (samtools faidx, hisat2-build
    and freebayes need plain files) and says that `iflas init-ref` writes the uncompressed copies."""
    import gzip
    with gzip.open(tmp_path / "genome.fa.gz", "wt") as f:
        f.write(">chr1\nACGT\n")
    cfg_path, tsv_path = _write(tmp_path, yaml_text=CONFIG_YAML.replace("genome.fa\n", "genome.fa.gz\n"))
    with pytest.raises(ValidationError, match="ref_genome is compressed.*iflas init-ref"):
        load_config(cfg_path, tsv_path)


def test_tmp_dir_receives_the_temporary_files(tmp_path, monkeypatch):
    """dirs.tmp_dir: once a run starts, pybedtools, Python and the external tools write their temporary files there
    (legacy set pybedtools' temp dir from it; it was read and ignored: BUG-184). Loading a config changes nothing."""
    import subprocess
    import tempfile
    from types import SimpleNamespace

    import pybedtools

    from iflas.cli import _load_config

    monkeypatch.setenv("TMPDIR", os.environ.get("TMPDIR", "/tmp"))      # restored after the test
    monkeypatch.setattr(tempfile, "tempdir", tempfile.tempdir)
    previous = pybedtools.get_tempdir()
    scratch = tmp_path / "scratch"
    try:
        cfg_path, tsv_path = _write(tmp_path, yaml_text=CONFIG_YAML.replace("  out_dir: out\n", f"  out_dir: out\n  tmp_dir: {scratch}\n"))
        assert load_config(cfg_path, tsv_path).dirs.tmp_dir == str(scratch) and not scratch.exists()
        _load_config(SimpleNamespace(config=cfg_path, samples=tsv_path))
        bed = pybedtools.BedTool("chr1\t1\t10\n", from_string=True)
        assert os.path.dirname(bed.fn) == str(scratch) and tempfile.gettempdir() == str(scratch)
        assert subprocess.run(["mktemp"], capture_output=True, text=True).stdout.startswith(str(scratch))
    finally:
        pybedtools.set_tempdir(previous)


def test_unknown_ref_strain_rejected(tmp_path):
    bad_tsv = SAMPLES_TSV.replace("s2\tB73", "s2\tMo17")
    cfg_path, tsv_path = _write(tmp_path, tsv_text=bad_tsv)
    with pytest.raises(ValidationError):
        load_config(cfg_path, tsv_path)


def test_extra_field_forbidden(tmp_path):
    bad_yaml = CONFIG_YAML + "  bogus_key: 1\n"
    cfg_path, tsv_path = _write(tmp_path, yaml_text=bad_yaml)
    with pytest.raises(ValidationError):
        load_config(cfg_path, tsv_path)


def test_merged_tgs_sample_string_matches_legacy_shape():
    sample = MergedTgsSample()
    sample.project_name = "demo"
    sample.sample_name = "B73_all"
    sample.ref_strain = "B73"
    sample.tgs_plat = "pacbio"
    sample.strategy = "iso-seq"
    sample.data_processed_location = ["s1.flnc.fa", "s2.flnc.fa"]
    sample.ngs_left_reads = "s1.R1,s2.R1"
    sample.ngs_right_reads = "s1.R2,s2.R2"
    sample.ngs_reads_length = 150
    sample.ngs_junctions = "junctions.bed"
    sample.single_run_threads = 3

    assert str(sample) == (
        "demo:B73_all:B73, pacbio:iso-seq:['s1.flnc.fa', 's2.flnc.fa'], "
        "s1.R1,s2.R1:s1.R2,s2.R2:paired:150:junctions.bed, 3"
    )


def test_merge_samples_by_strain_preserves_legacy_layout_and_reads(tmp_path):
    polya = tmp_path / "s1.polya.tsv"
    polya.write_text("read\tpolya\n")
    samples = [
        Sample(
            sample_name="s1",
            project_name="demo",
            ref_strain="B73",
            strain="B73",
            tgs_plat="pacbio",
            strategy="iso-seq",
            data_processed_location="/flnc/s1.fa",
            polya_location=str(polya),
            ngs_left_reads="s1.L0;s1.L1;s1.L2",
            ngs_right_reads="s1.R0;s1.R1",
            use_fmlrc2=False,
        ),
        Sample(
            sample_name="s2",
            project_name="demo",
            ref_strain="B73",
            strain="B73",
            tgs_plat="nanopore",
            strategy="direct",
            data_processed_location="/flnc/s2.fa",
            polya_location=str(tmp_path / "missing.polya.tsv"),
            ngs_left_reads="s2.L0;s2.L1",
            ngs_right_reads="s2.R0;s2.R1;s2.R2",
            use_fmlrc2=True,
        ),
    ]

    merged = merge_samples_by_strain(samples)
    sample = merged["demo"]["B73"]["B73"]

    assert sample.project_name == "demo"
    assert sample.sample_name == "B73_all"
    assert sample.ref_strain == "B73"
    assert sample.strain == "B73"
    assert sample.condition == "B73_all"
    assert sample.tgs_plat == "pacbio"
    assert sample.strategy == "iso-seq"
    assert sample.data_processed_location == ["/flnc/s1.fa", "/flnc/s2.fa"]
    assert sample.polya_location == {"s1": str(polya)}
    assert sample.ngs_left_reads == "s1.L0,s2.L0;s1.L1,s2.L1"
    assert sample.ngs_right_reads == "s1.R0,s2.R0;s1.R1,s2.R1"
    assert sample.ngs_reads_paired == "paired"
    assert sample.use_fmlrc2 is False


def test_minimap2_defaults_carry_the_legacy_junction_filter_support(tmp_path):
    cfg = load_config(*_write(tmp_path))
    assert cfg.minimap2.junc_comb_support == 2      # legacy `mapping -jcs` default


def test_sample_fast5_dir_is_resolved_relative_to_the_sheet(tmp_path):
    from iflas.config import load_config, merge_samples_by_strain, read_samples_tsv

    (tmp_path / "signal").mkdir()
    (tmp_path / "samples.tsv").write_text(
        "sample_name\tref_strain\ttgs_plat\tdata_processed_location\tfast5_dir\n"
        "s1\tB73\tnanopore\treads.fq\tsignal\n")
    rows = read_samples_tsv(str(tmp_path / "samples.tsv"))
    assert rows[0]["fast5_dir"] == str(tmp_path / "signal")
    (tmp_path / "config.yaml").write_text(
        "reference:\n  - ref_strain: B73\n    ref_genome: genome.fa\n    ref_gtf: genes.gtf\ndirs:\n  out_dir: out\n")
    cfg = load_config(str(tmp_path / "config.yaml"), str(tmp_path / "samples.tsv"))
    merged = next(iter(next(iter(merge_samples_by_strain(cfg.samples).values())).values())).values()
    assert next(iter(merged)).fast5_dir == {"s1": str(tmp_path / "signal")}


def test_output_and_temporary_dirs_are_relative_to_the_config_file(tmp_path, monkeypatch):
    """BUG-199: examples/config.yaml says relative paths are resolved against the config file's directory, but
    dirs.out_dir/tmp_dir were taken relative to the working directory, so `iflas all -c project/config.yaml` run
    elsewhere wrote its results there."""
    project = tmp_path / "project"
    project.mkdir()
    cfg_path, tsv_path = _write(project, yaml_text=CONFIG_YAML.replace("  out_dir: out\n", "  out_dir: out\n  tmp_dir: scratch\n"))
    monkeypatch.chdir(tmp_path)
    cfg = load_config(cfg_path, tsv_path)
    assert (cfg.dirs.out_dir, cfg.dirs.tmp_dir) == (str(project / "out"), str(project / "scratch"))
