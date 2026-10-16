import subprocess
from pathlib import Path


def test_run_r_function_writes_script_and_routes_to_diffsplice(tmp_path):
    from iflas.rscript import run_r_function

    class DummyRunner:
        def __init__(self):
            self.cmd = None
            self.kwargs = None
            self.script = None

        def run(self, cmd, **kwargs):
            self.cmd = cmd
            self.kwargs = kwargs
            self.script = Path(cmd[1]).read_text()
            return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")

    runner = DummyRunner()
    proc = run_r_function(
        "plotThing <- function(path, cutoff, flag, nothing) { print(path) }",
        "plotThing",
        ["a b", 0.05, True, None],
        cwd=tmp_path,
        runner=runner,
    )

    assert proc.returncode == 0
    assert runner.cmd[0] == "Rscript"
    assert runner.kwargs["cwd"] == str(tmp_path)
    assert runner.kwargs["conda_env"] == "iflas-diffsplice"
    assert 'plotThing("a b", 0.05, TRUE, NULL)' in runner.script
    assert not Path(runner.cmd[1]).exists()


def test_iflas_go_cli_uses_packaged_rscript_helper(tmp_path, monkeypatch):
    import iflas.stages.go as go_stage
    from iflas.cli import main

    target = tmp_path / "target.txt"
    target.write_text("gene1\n")
    gene2go = tmp_path / "gene2go.tsv"
    gene2go.write_text("Gene\tGO\n")

    calls = []

    def fake_run_r_function(source, function_name, args, cwd=None, runner=None):
        calls.append((source, function_name, args, cwd, runner))
        return subprocess.CompletedProcess(["Rscript"], 0)

    monkeypatch.setattr(go_stage, "run_r_function", fake_run_r_function)

    rc = main([
        "go",
        "-tg",
        str(target),
        "-s",
        "S1",
        "-bg",
        str(gene2go),
        "-o",
        "S1",
        "--outdir",
        str(tmp_path / "out"),
        "-cutoff",
        "0.05",
        "-filterBy",
        "p.adjust",
        "-showCategory",
        "10",
    ])

    assert rc == 0
    assert len(calls) == 1
    source, function_name, r_args, cwd, runner = calls[0]
    assert "plotTargetGenesGoEnrichment <- function" in source
    assert function_name == "plotTargetGenesGoEnrichment"
    assert r_args == [str(target), "S1", str(gene2go), "S1", 0.05, "p.adjust", 10]
    assert cwd == str(tmp_path / "out" / "GO")
    assert runner is None


def test_report_cli_requires_an_output_selection(capsys):
    from iflas.cli import main

    assert main(["report"]) == 2
    assert "no output selected" in capsys.readouterr().err


def test_r_calls_define_the_function_once_and_call_it_per_argument_list():
    """Report plots per event share one Rscript session (packages and the reference TxDb load once)."""
    from iflas.rscript import build_r_calls, build_r_script

    script = build_r_calls("plotIt <- function(a, b) NULL", "plotIt", [["x", 1], ["y", 2]])
    assert script == 'plotIt <- function(a, b) NULL\n\nplotIt("x", 1)\nplotIt("y", 2)\n'
    assert build_r_script("f <- function(a) a", "f", ["z"]) == 'f <- function(a) a\n\nf("z")\n'
