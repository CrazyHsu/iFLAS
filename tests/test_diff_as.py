import subprocess
from pathlib import Path

import pandas as pd


def _write_comp_table(path: Path, bam_a: Path, bam_b: Path, bam_c: Path) -> None:
    path.write_text(
        "project\tsampleName\tcondition\trepeat\tbamFile\tpaired\treadsLength\n"
        f"demo\tA1\tA\trep1\t{bam_a}\tpaired\t150\n"
        f"demo\tA2\tA\trep2\t{bam_b}\tpaired\t150\n"
        f"demo\tB1\tB\trep1\t{bam_c}\tpaired\t150\n"
    )


class FakeRunner:
    def __init__(self):
        self.calls = []

    def run(self, cmd, **kwargs):
        if kwargs.get("stdout") is not None:
            kwargs["stdout"].write("rmats stdout\n")
        if kwargs.get("stderr") is not None:
            kwargs["stderr"].write("rmats stderr\n")
        self.calls.append((cmd, kwargs))
        if isinstance(cmd, list) and cmd[0] == "rmats.py":      # fake rMATS: empty JC tables in --od
            od = Path(kwargs["cwd"]) / cmd[cmd.index("--od") + 1]
            od.mkdir(parents=True, exist_ok=True)
            for table in ("RI.MATS.JC.txt", "SE.MATS.JC.txt", "A3SS.MATS.JC.txt", "A5SS.MATS.JC.txt"):
                (od / table).write_text("ID\tGeneID\tFDR\n")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")


def test_run_rmats_from_bam_table_writes_inputs_and_command_contract(tmp_path):
    from iflas.stages.diff_as import run_rmats_from_bam_table

    bam_a = tmp_path / "A1.bam"
    bam_b = tmp_path / "A2.bam"
    bam_c = tmp_path / "B1.bam"
    for bam in (bam_a, bam_b, bam_c):
        bam.write_text("bam")
    comp = tmp_path / "comp.tsv"
    _write_comp_table(comp, bam_a, bam_b, bam_c)
    gtf = tmp_path / "merged.gtf"
    gtf.write_text("gtf")

    runner = FakeRunner()
    out = run_rmats_from_bam_table(comp, gtf, out_dir=tmp_path / "out", threads=3, runner=runner)

    das_dir = tmp_path / "out" / "das"
    b1_file = das_dir / "demo_A.txt"
    b2_file = das_dir / "demo_B.txt"
    comp_out = das_dir / "demo_A_vs_demo_B"
    rmats_log = das_dir / "demo_A_vs_demo_B.rmats.log"
    assert b1_file.read_text() == f"{bam_a.resolve()},{bam_b.resolve()}\n"
    assert b2_file.read_text() == f"{bam_c.resolve()}\n"
    assert len(runner.calls) == 1
    cmd, kwargs = runner.calls[0]
    assert cmd == [
        "rmats.py",
        "--b1",
        "demo_A.txt",
        "--b2",
        "demo_B.txt",
        "--gtf",
        str(gtf.resolve()),
        "--od",
        "demo_A_vs_demo_B",
        "--tmp",
        "demo_A_vs_demo_B.tmp",
        "-t",
        "paired",
        "--readLength",
        "150",
        "--variable-read-length",
        "--tstat",
        "3",
        "--nthread",
        "3",
    ]
    assert kwargs["cwd"] == str(das_dir)
    assert kwargs["stdout"].name == str(rmats_log)
    assert kwargs["stderr"].name == str(rmats_log)
    assert rmats_log.read_text() == "rmats stdout\nrmats stderr\n"
    assert out["das_dir"] == str(das_dir)
    assert out["comparisons"][0]["out_dir"] == str(comp_out)
    assert out["comparisons"][0]["log_file"] == str(rmats_log)


def test_iflas_diff_as_cli_uses_diffsplice_toolrunner(tmp_path, monkeypatch, capsys):
    from iflas.cli import main
    import iflas.stages.diff_as as diff_stage

    bam_a = tmp_path / "A1.bam"
    bam_b = tmp_path / "A2.bam"
    bam_c = tmp_path / "B1.bam"
    for bam in (bam_a, bam_b, bam_c):
        bam.write_text("bam")
    comp = tmp_path / "comp.tsv"
    _write_comp_table(comp, bam_a, bam_b, bam_c)
    gtf = tmp_path / "merged.gtf"
    gtf.write_text("gtf")

    runner = FakeRunner()
    created = {}

    def fake_toolrunner(env_map=None):
        created["env_map"] = env_map
        return runner

    monkeypatch.setattr(diff_stage, "ToolRunner", fake_toolrunner)

    rc = main([
        "diff_as",
        "--comp-cond",
        str(comp),
        "--gtf",
        str(gtf),
        "-o",
        str(tmp_path / "out"),
        "--threads",
        "4",
    ])

    assert rc == 0
    assert created["env_map"] == {"rmats.py": "iflas-diffsplice"}
    assert runner.calls[0][0][-1] == "4"
    assert "Ran diff_as via rMATS-turbo" in capsys.readouterr().out


def _bed12plus(name, chrom, exons, gene, strand="+"):
    start, end = exons[0][0], exons[-1][1]
    sizes = ",".join(str(e - s) for s, e in exons)
    rel = ",".join(str(s - start) for s, e in exons)
    return "\t".join(map(str, [chrom, start, end, name, 0, strand, start, end, "0,0,0", len(exons), sizes, rel, gene]))


def _sample_dir(root, name, isoforms, as_isos):
    """isoformGrouped.bed12+ with `isoforms` and confident AS files naming `as_isos` (inc/exc pairs)."""
    base = root / "out" / "demo" / name
    (base / "refine").mkdir(parents=True)
    (base / "refine" / "isoformGrouped.bed12+").write_text("".join(_bed12plus(*iso) + "\n" for iso in isoforms))
    lr = base / "as_events" / "ordinary_as" / "LR"
    lr.mkdir(parents=True)
    for fn in ("A3SS.confident.bed6+", "A5SS.confident.bed6+", "IR.confident.bed6+"):
        (lr / fn).write_text("")
    inc, exc = as_isos
    (lr / "SE.confident.bed12+").write_text("\t".join(["chr1", "150", "200", "geneA:100@150-200@300", "500", "+", "150", "200", "0,0,0", "1", "50", "0",
                                                        "geneA", "100", "300", inc, "1", exc, "1"]) + "\n")
    return base


def test_merge_isoforms_merges_identical_junction_chains_across_samples(tmp_path):
    from iflas.stages.diff_as import merge_isoforms

    s1 = _sample_dir(tmp_path, "s1", [("PB.1.1", "chr1", [(100, 150), (200, 300)], "geneA"),
                                      ("PB.1.2", "chr1", [(100, 300)], "geneA"),
                                      ("PB.2.1", "chr1", [(1000, 1100)], "geneB")], ("PB.1.1", "PB.1.2"))
    s2 = _sample_dir(tmp_path, "s2", [("PB.7.1", "chr1", [(90, 150), (200, 320)], "geneA"),        # same chain as s1 PB.1.1, longer span
                                      ("PB.7.2", "chr1", [(90, 320)], "geneA")], ("PB.7.1", "PB.7.2"))
    runner = FakeRunner()
    das = tmp_path / "das"
    das.mkdir()

    gtf = merge_isoforms([("s1", s1), ("s2", s2)], das, runner=runner)

    rows = [l.rstrip("\n").split("\t") for l in open(das / "all_sample_merged_iso.bed")]
    by_name = {r[3]: r for r in rows}
    # PB.2.1 has no AS event -> dropped. The two mono-exon isoforms have different spans and stay apart
    # (BUG-150 fixed at CP3: mono-exon isoforms merge only on an identical span, not on the empty junction chain).
    assert set(by_name) == {"s1_PB.1.1+s2_PB.7.1", "s1_PB.1.2", "s2_PB.7.2"}
    assert by_name["s1_PB.1.1+s2_PB.7.1"][1:3] == ["90", "320"]                           # longest span is the representative
    assert by_name["s1_PB.1.1+s2_PB.7.1"][12] == "geneA"
    assert gtf == str(das / "all_sample_merged_iso.gtf")
    assert "bed2gpe.pl -g 13 all_sample_merged_iso.bed" in runner.calls[0][0]
    assert runner.calls[1][0] == ["genePredToGtf", "file", "all_sample_merged_iso.gpe", "all_sample_merged_iso.gtf"]


def test_merge_isoforms_pu_filter_requires_hq_isoforms(tmp_path):
    import pytest
    from iflas.stages.diff_as import merge_isoforms

    s1 = _sample_dir(tmp_path, "s1", [("PB.1.1", "chr1", [(100, 150), (200, 300)], "geneA")], ("PB.1.1", "PB.1.1"))
    with pytest.raises(FileNotFoundError):
        merge_isoforms([("s1", s1)], tmp_path, pu_filter=True, runner=FakeRunner())


def test_write_sig_diff_as_filters_by_fdr(tmp_path):
    from iflas.stages.diff_as import write_sig_diff_as

    comp = tmp_path / "das" / "demo_A_vs_demo_B"
    comp.mkdir(parents=True)
    for table in ("RI.MATS.JC.txt", "SE.MATS.JC.txt", "A3SS.MATS.JC.txt", "A5SS.MATS.JC.txt"):
        pd.DataFrame({"ID": [1, 2, 3], "GeneID": ["g1", "g2", "g3"], "FDR": [0.01, 0.5, 0.05]}).to_csv(comp / table, sep="\t", index=False)

    out = write_sig_diff_as(comp)

    sig_dir = tmp_path / "das" / "demo_A_vs_demo_B.sigDiffAS"
    assert set(out) == {"IR", "SE", "A3SS", "A5SS"} and out["IR"] == str(sig_dir / "IR.sig.txt")
    assert pd.read_csv(sig_dir / "SE.sig.txt", sep="\t")["GeneID"].tolist() == ["g1", "g3"]


def test_bam_records_from_samples_uses_sample_alignments(tmp_path):
    from iflas.config import load_config
    from iflas.stages.diff_as import bam_records_from_samples

    cfg = tmp_path / "config.yaml"
    cfg.write_text(f"reference:\n  - ref_strain: B73\n    ref_genome: genome.fa\n    ref_gtf: genes.gtf\ndirs:\n  out_dir: {tmp_path / 'out'}\n")
    tsv = tmp_path / "samples.tsv"
    tsv.write_text("sample_name\tproject_name\tref_strain\tcondition\tngs_left_reads\tngs_right_reads\tngs_reads_length\n"
                   "a1\tdemo\tB73\tA\ta1_R1.fq;a1b_R1.fq\ta1_R2.fq;a1b_R2.fq\t150\n"
                   "b1\tdemo\tB73\tB\tb1_R1.fq\tb1_R2.fq\t150\n"
                   "c1\tdemo\tB73\tC\tc1_R1.fq\tc1_R2.fq\t150\n")
    config = load_config(str(cfg), str(tsv))

    records = bam_records_from_samples(config, conditions=["A", "B"])

    assert [(r.sample_name, r.condition, r.repeat, r.comparison_name) for r in records] == [
        ("a1", "A", "repeat0", "demo_A"), ("a1", "A", "repeat1", "demo_A"), ("b1", "B", "repeat0", "demo_B")]
    assert records[1].bam_file == str(tmp_path / "out" / "demo" / "a1" / "mapping" / "rna-seq" / "alignment" / "repeat1.sorted.bam")
    assert records[0].paired == "paired" and records[0].reads_length == 150
