"""Unit tests for iflas.cli (entry point + migrate-config command + examples)."""

import os
from pathlib import Path
import io
from types import SimpleNamespace

import pytest

import iflas.cli as cli
from iflas.cli import main
from iflas.config import load_config

LEGACY_INI = os.path.join(os.path.dirname(__file__), "..", "test_data", "config.cfg")
EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


class _FakeMappingRunner:
    def __init__(self):
        self.calls = []

    def run(self, cmd, **kwargs):
        self.calls.append((cmd, kwargs))
        if isinstance(cmd, list) and cmd[:3] == ["regtools", "junctions", "extract"]:
            raw_bed = cmd[-1]
            with open(raw_bed, "w") as fh:
                fh.write("chr1\t90\t210\tj_low\t10\t+\t90\t210\t0,0,0\t2\t10,10\t0,110\n")
                fh.write("chr1\t290\t410\tj_high\t11\t+\t290\t410\t0,0,0\t2\t10,10\t0,110\n")
        elif isinstance(cmd, list) and cmd[:2] == ["bash", "-c"] and "bamToBed" in cmd[2]:
            # the junction filter reads the bed12 of the minimap2 alignments: one single-exon read
            with open(os.path.join(kwargs["cwd"], "tmp.bed12"), "w") as fh:
                fh.write("chr1\t100\t180\tr_single\t0\t+\t100\t180\t0,0,0\t1\t80\t0\n")
        elif isinstance(cmd, str) and "> tmp.sam" in cmd:
            open(os.path.join(kwargs["cwd"], "tmp.sam"), "w").close()   # removed by the filter afterwards


class _FakeCollapseRunner:
    def __init__(self):
        self.calls = []

    def run(self, cmd, **kwargs):
        self.calls.append((cmd, kwargs))
        out_prefix = cmd[cmd.index("-o") + 1]
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
        with open(out_prefix + ".collapsed.good.gff", "w") as fh:
            fh.write(
                "chr1\tPacBio\ttranscript\t1\t100\t.\t+\t.\t"
                'gene_id "PB.1"; transcript_id "PB.1.1";\n'
            )
        return SimpleNamespace(returncode=1, stderr="empty bad.gff")


def test_cli_version():
    with pytest.raises(SystemExit) as ei:
        main(["--version"])
    assert ei.value.code == 0


def test_cli_no_command_prints_help():
    assert main([]) == 1


def test_cli_diff_as_missing_args_returns_2():
    assert main(["diff_as"]) == 2


def test_cli_count_nonempty_lines_closes_handle(monkeypatch):
    close_calls = []

    class TrackingHandle(io.StringIO):
        def close(self):
            close_calls.append(True)
            super().close()

    def fake_open(path):
        assert path == "events.tsv"
        return TrackingHandle("\nfirst\n \nsecond\n")

    monkeypatch.setattr(cli, "open", fake_open, raising=False)

    assert cli._count_nonempty_lines("events.tsv") == 2
    assert close_calls == [True]


def test_cli_accepts_safe_legacy_flag_aliases_for_migrated_commands():
    parser = cli.build_parser()

    diff_args = parser.parse_args(["diff_as", "-d", "comparisons.tsv", "--gtf", "merged.gtf"])
    assert diff_args.comp_cond == "comparisons.tsv"

    visual_args = parser.parse_args(["visual_as", "-b", "sample", "-r", "ref.gpe", "-g", "geneA,geneB"])
    assert visual_args.gene == ["geneA,geneB"]

    find_pa_args = parser.parse_args([
        "find_pa",
        "--reads-bed",
        "reads.bed12+",
        "-g",
        "tofu.group.txt",
        "-o",
        "paCluster.bed8+",
        "-pa_sup",
        "7",
        "-pa_rpkm",
        "1.5",
        "-conf_pa",
        "confident.pa.bed",
    ])
    assert find_pa_args.pa_support == 7
    assert find_pa_args.pa_rpkm == 1.5
    assert find_pa_args.confident_pa == "confident.pa.bed"

    pu_args = parser.parse_args([
        "pu_filter",
        "-b",
        "sample",
        "--ref-genome",
        "genome.fa",
        "--ref-bed",
        "ref.bed12+",
        "-feature_file",
        "isoFeatures.txt",
        "-filter_score",
        "0.7",
        "-pos_fl_cov",
        "3",
        "-pos_min_junc_rpkm",
        "0.2",
    ])
    assert pu_args.feature_file == "isoFeatures.txt"
    assert pu_args.filter_score == 0.7
    assert pu_args.pos_fl_coverage == 3
    assert pu_args.pos_min_junc_rpkm == 0.2


def test_cli_all_loads_config_and_runs_orchestrator(tmp_path, monkeypatch):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        """\
reference:
  - ref_strain: B73
    ref_genome: genome.fa
    ref_gtf: genes.gtf
    ref_gpe: genes.gpe
    ref_bed: genes.bed12+
options:
  threads: 6
dirs:
  out_dir: out
"""
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_name\tproject_name\tref_strain\tdata_processed_location\n"
        "s1\tdemo\tB73\tflnc.fa\n"
    )
    seen = {}

    def fake_run_pipeline(config, *, stages, max_workers, merge_samples=False):
        seen["sample_names"] = [s.sample_name for s in config.samples]
        seen["stages"] = stages
        seen["max_workers"] = max_workers
        seen["merge_samples"] = merge_samples
        return []

    import iflas.orchestrator

    monkeypatch.setattr(iflas.orchestrator, "run_pipeline", fake_run_pipeline)
    rc = main(["all", "-c", str(cfg), "-s", str(samples), "--stages", "mapping,find_as", "--jobs", "3"])
    assert rc == 0
    assert seen == {
        "sample_names": ["s1"],
        "stages": ("mapping", "find_as"),
        "max_workers": 3,
        "merge_samples": True,
    }
    rc = main(["all", "-c", str(cfg), "-s", str(samples)])         # no --stages: the whole published chain, gated
    assert rc == 0 and seen["stages"] is None and seen["max_workers"] is None


def test_cli_mapping_runs_stage_logic_with_toolrunner_contract(tmp_path, monkeypatch):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        f"""\
reference:
  - ref_strain: B73
    ref_genome: {tmp_path / "genome.fa"}
    ref_gtf: {tmp_path / "genes.gtf"}
    ref_gpe: {tmp_path / "genes.gpe"}
    ref_bed: {tmp_path / "genes.bed12+"}
    ref_mm2_index: {tmp_path / "genome.mm2"}
    hisat2_index: {tmp_path / "hisat2" / "B73"}
minimap2:
  max_intron_length: 70000
options:
  threads: 6
dirs:
  out_dir: {tmp_path / "out"}
"""
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_name\tproject_name\tref_strain\ttgs_plat\tdata_processed_location\t"
        "ngs_left_reads\tngs_right_reads\n"
        f"s1\tdemo\tB73\tnanopore\t{tmp_path / 'flnc.fa'}\t"
        f"{tmp_path / 'ngs_R1.fastq.gz'}\t{tmp_path / 'ngs_R2.fastq.gz'}\n"
    )
    (tmp_path / "genes.bed12+").write_text(
        "chr1\t100\t400\tref1\t0\t+\t100\t400\t0,0,0\t2\t100,100\t0,200\n"
    )
    runner = _FakeMappingRunner()

    import iflas.stages.align

    monkeypatch.setattr(iflas.stages.align, "ToolRunner", lambda: runner)
    rc = main(["mapping", "-c", str(cfg), "-s", str(samples), "--sample", "s1"])

    base = tmp_path / "out" / "demo" / "s1"
    assert rc == 0
    commands = [call[0] for call in runner.calls]
    # legacy mapping() order: NGS first (its junctions feed the FLNC junction-combination filter), then TGS
    assert commands[0] == (
        f"hisat2 -x {tmp_path / 'hisat2' / 'B73'} -1 {tmp_path / 'ngs_R1.fastq.gz'} "
        f"-2 {tmp_path / 'ngs_R2.fastq.gz'} --dta -p 6 --reorder --max-intronlen 70000 | "
        f"samtools sort -@ 6 -o {base / 'mapping' / 'rna-seq' / 'alignment' / 'repeat0.sorted.bam'} -"
    )
    assert commands[2] == [
        "stringtie",
        str(base / "mapping" / "rna-seq" / "alignment" / "repeat0.sorted.bam"),
        "-o",
        str(base / "mapping" / "rna-seq" / "alignment" / "repeat0.gtf"),
        "-p",
        "6",
    ]
    assert commands[6][:4] == ["regtools", "junctions", "extract", "-a"]
    sam = base / "mapping" / "flnc.mm2.sam"
    assert commands[7] == (          # nanopore maps against the genome, where -k14 applies (BUG-192)
        f"minimap2 -ax splice -G 70000 -uf -k14 --secondary=no --MD "
        f"{tmp_path / 'genome.fa'} {tmp_path / 'flnc.fa'} -t 6 > {sam}"
    )
    assert commands[8][:2] == ["bash", "-c"] and "bamToBed" in commands[8][2]
    shell = [c for c in commands[9:] if isinstance(c, str)]
    assert any("filter.pl" in c and str(sam) in c for c in shell)
    assert any("samAddTag.pl" in c and "flnc.addCVandID.bed12+" in c for c in shell)
    assert any("samtools view -bS -@ 6 flnc.mm2.sam > flnc.mm2.sorted.bam" in c for c in shell)
    assert commands[-2] == ["samtools", "index", str(base / "mapping" / "flnc.mm2.sorted.bam")]
    assert commands[-1] == (
        f"seqkit fq2fa -w 0 {tmp_path / 'flnc.fa'} > {base / 'mapping' / 'flnc.processed.fa'}"
    )
    assert (base / "mapping" / "flnc.filterReads.lst").read_text() == "r_single\n"
    junctions = base / "mapping" / "rna-seq" / "reassembly" / "junctions.bed"
    assert junctions.read_text() == (
        "chr1\t290\t410\tj_high\t11\t+\t290\t410\t0,0,0\t2\t10,10\t0,110\n"
    )


def test_cli_preproc_runs_stage_logic_with_toolrunner_contract(tmp_path, monkeypatch):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        f"""\
reference:
  - ref_strain: B73
    ref_genome: {tmp_path / "genome.fa"}
    ref_gtf: {tmp_path / "genes.gtf"}
ccs:
  min_pass: 3
  min_subread_length: 75
  min_predicted_accuracy: 0.91
dirs:
  out_dir: {tmp_path / "out"}
"""
    )
    raw = tmp_path / "s1.subreads.bam"
    primer = tmp_path / "primer.fa"
    raw.write_text("bam")
    primer.write_text(">p\nAC\n")
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_name\tproject_name\tref_strain\ttgs_plat\tstrategy\tdata_location\tprimer\n"
        f"s1\tdemo\tB73\tpacbio\tsequel\t{raw}\t{primer}\n"
    )

    class FakePreprocRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs))
            if cmd == ["ccs", "--version"]:
                return SimpleNamespace(returncode=0, stdout="ccs 6.4.0\n", stderr="")
            if cmd[0] == "lima":
                (Path(kwargs["cwd"]) / "fl.5p--3p.bam").write_text("")
            return SimpleNamespace(returncode=0, stdout="", stderr="")

    runner = FakePreprocRunner()

    import iflas.stages.preprocess

    monkeypatch.setattr(iflas.stages.preprocess, "ToolRunner", lambda env_map=None: runner)
    rc = main(["preproc", "-c", str(cfg), "-s", str(samples), "--sample", "s1", "--threads", "5"])

    base = tmp_path / "out" / "demo" / "s1"
    assert rc == 0
    commands = [call[0] for call in runner.calls]
    assert commands[:6] == [
        ["pbindex", str(raw)],
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
        ["ccs", "--version"],
        [
            "ccs",
            "s1.subreadset.xml",
            "CCS.bam",
            "-j",
            "5",
            "--min-passes",
            "3",
            "--min-length",
            "50",       # ccs>=4: draft (CCS) length = min_ccs_length default, not min_subread_length (BUG-102)
            "--min-rq",
            "0.91",
        ],
        ["lima", "CCS.bam", str(primer), "fl.bam", "--isoseq", "--dump-clips", "-j", "5"],
        ["isoseq", "refine", "fl.5p--3p.bam", str(primer), "flnc.bam", "--require-polya", "-j", "5"],
    ]
    assert commands[6] == ["samtools", "fastq", "flnc.bam"]
    assert (base / "preprocess" / "pacbio" / "rawFlnc.fq").is_file()


def test_cli_collapse_runs_stage_logic_with_toolrunner_contract(tmp_path, monkeypatch):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        f"""\
reference:
  - ref_strain: B73
    ref_genome: {tmp_path / "genome.fa"}
    ref_gtf: {tmp_path / "genes.gtf"}
collapse:
  min_identity: 0.92
  min_coverage: 0.85
  fl_coverage: 4
  max_fuzzy_junction: 9
  max_5_diff: 77
  max_3_diff: 33
  dun_merge_5_shorter: false
dirs:
  out_dir: {tmp_path / "out"}
"""
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_name\tproject_name\tref_strain\tdata_processed_location\n"
        f"s1\tdemo\tB73\t{tmp_path / 'flnc.fa'}\n"
    )
    runner = _FakeCollapseRunner()

    import iflas.stages.collapse

    monkeypatch.setattr(iflas.stages.collapse, "ToolRunner", lambda: runner)
    rc = main(["collapse", "-c", str(cfg), "-s", str(samples), "--sample", "s1"])

    base = tmp_path / "out" / "demo" / "s1"
    assert rc == 0
    assert runner.calls == [
        (
            [
                "collapse_isoforms_by_sam.py",
                "--input",
                str(tmp_path / "flnc.fa"),
                "-b",
                str(base / "mapping" / "flnc.mm2.sorted.bam"),
                "-o",
                str(base / "collapse" / "tofu"),
                "-c",
                "0.85",
                "-i",
                "0.92",
                "--max_fuzzy_junction",
                "9",
                "--max_5_diff",
                "77",
                "--max_3_diff",
                "33",
                "--flnc_coverage",
                "4",
            ],
            {"conda_env": "iflas-cupcake", "check": False},
        )
    ]
    assert (base / "collapse" / "tofu.collapsed.good.gff").is_file()


def test_cli_collapse_skips_rows_without_long_reads(tmp_path, monkeypatch, capsys):
    """BUG-202: `iflas collapse` over a sample table stopped at the first NGS-only row ("Input file None does not
    exist"); such rows are skipped with a notice, and selecting only such rows is an error."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text(f"""\
reference:
  - ref_strain: B73
    ref_genome: {tmp_path / "genome.fa"}
    ref_gtf: {tmp_path / "genes.gtf"}
dirs:
  out_dir: {tmp_path / "out"}
""")
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_name\tproject_name\tref_strain\tdata_processed_location\tngs_left_reads\n"
        f"s1\tdemo\tB73\t{tmp_path / 'flnc.fa'}\tNone\n"
        f"s2\tdemo\tB73\tNone\t{tmp_path / 's2_R1.fq.gz'}\n"
    )
    runner = _FakeCollapseRunner()

    import iflas.stages.collapse

    monkeypatch.setattr(iflas.stages.collapse, "ToolRunner", lambda: runner)
    assert main(["collapse", "-c", str(cfg), "-s", str(samples)]) == 0
    assert [call[0][2] for call in runner.calls] == [str(tmp_path / "flnc.fa")]
    assert "Skip collapse for s2: no long reads" in capsys.readouterr().out

    assert main(["collapse", "-c", str(cfg), "-s", str(samples), "--sample", "s2"]) == 2
    assert "none of the selected samples has long reads" in capsys.readouterr().err


def _bed12plus(name, exons, group, chrom="chr1", strand="+", gene="geneX"):
    start = exons[0][0]
    end = exons[-1][1]
    sizes = ",".join(str(e - s) for s, e in exons)
    rel_starts = ",".join(str(s - start) for s, e in exons)
    fields = [
        chrom,
        start,
        end,
        name,
        0,
        strand,
        start,
        end,
        "0,0,0",
        len(exons),
        sizes,
        rel_starts,
        gene,
        gene,
        group,
    ]
    return "\t".join(map(str, fields))


def test_cli_palen_as_runs_file_level_stage(tmp_path, capsys):
    base = tmp_path / "sample"
    hq_dir = base / "hqIsoforms"
    collapse_dir = base / "collapse"
    mapping_dir = base / "mapping"
    as_dir = base / "as_events" / "ordinary_as" / "LR"
    pa_dir = base / "as_events" / "pa"
    for d in (hq_dir, collapse_dir, mapping_dir, as_dir, pa_dir):
        d.mkdir(parents=True)

    inc = [f"inc{i}" for i in range(30)]
    exc = [f"exc{i}" for i in range(30)]
    (hq_dir / "hq.collapsed.bed12+").write_text(
        "\n".join(
            [
                _bed12plus("PB.1.1", [(100, 200), (300, 400)], "PB.1.1"),
                _bed12plus("PB.1.2", [(100, 200), (500, 600)], "PB.1.2"),
            ]
        )
        + "\n"
    )
    (collapse_dir / "tofu.collapsed.group.txt").write_text(
        "PB.1.1\t" + ",".join(inc) + "\n" + "PB.1.2\t" + ",".join(exc) + "\n"
    )
    (mapping_dir / "flnc.addCVandID.bed12+").write_text("")
    se = [
        "chr1",
        "200",
        "300",
        "PB.1:se",
        "0",
        "+",
        "200",
        "300",
        "0,0,0",
        "1",
        "100",
        "0",
        "geneX",
        "200",
        "300",
        "PB.1.1",
        "30",
        "PB.1.2",
        "30",
    ]
    (as_dir / "SE.confident.bed12+").write_text("\t".join(se) + "\n")
    for name in ("IR.confident.bed6+", "A3SS.confident.bed6+", "A5SS.confident.bed6+"):
        (as_dir / name).write_text("")
    (pa_dir / "paCluster.bed8+").write_text(
        "\n".join(
            [
                "\t".join(["chr1", "100", "120", ",".join(inc), "30", "+", "0", "110", "PB.1"]),
                "\t".join(["chr1", "200", "220", ",".join(exc), "30", "+", "0", "210", "PB.1"]),
            ]
        )
        + "\n"
    )
    palen = tmp_path / "polya_results.tsv"
    palen.write_text(
        "readname\tqc_tag\tpolya_length\n"
        + "\n".join([f"{read}\tPASS\t{10 + i}" for i, read in enumerate(inc)])
        + "\n"
        + "\n".join([f"{read}\tPASS\t{200 + i}" for i, read in enumerate(exc)])
        + "\n"
    )

    rc = main(["palen_as", "-b", str(base), "--palen", str(palen), "--pa-support", "5"])

    assert rc == 0
    assert (base / "palenAS" / "mergeByJunc" / "SE.palenAndAS.sig.bed12+").is_file()
    assert (base / "palenAS" / "palenAPA" / "apaRelatedPalen.txt").is_file()
    assert "Wrote palen_as outputs" in capsys.readouterr().out


def test_cli_migrate_config(tmp_path):
    out = tmp_path / "config.yaml"
    tsv = tmp_path / "samples.tsv"
    rc = main(["migrate-config", LEGACY_INI, "-o", str(out), "-s", str(tsv)])
    assert rc == 0 and out.exists() and tsv.exists()
    cfg = load_config(str(out), str(tsv))
    assert cfg.reference[0].ref_strain == "B73"


def test_example_templates_load():
    cfg = load_config(os.path.join(EXAMPLES, "config.yaml"),
                      os.path.join(EXAMPLES, "samples.tsv"))
    assert cfg.reference[0].ref_strain == "B73"
    assert [s.sample_name for s in cfg.samples] == \
        ["B73_root_rep1", "B73_root_rep2", "B73_leaf_rep1"]
    assert cfg.options.threads == 8
