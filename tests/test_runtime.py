"""Unit tests for iflas.runtime (logging + guarded command runner)."""

import logging
import multiprocessing
import os
import sys

import pytest

from iflas import runtime
from iflas.runtime import (
    MyPool, NoDaemonProcess, ToolError, ToolRunner, conda_env_prefix, get_logger, require_tools, run, tool_version,
    which,
)


def _double(value):
    return value * 2


def _spawn_child_pool(value):
    with multiprocessing.Pool(processes=1) as pool:
        return pool.map(_double, [value])[0]


def test_which_and_require():
    assert which("python") is not None
    assert which("definitely_missing_tool_xyz_123") is None
    assert require_tools(["python"]) is True
    with pytest.raises(ToolError):
        require_tools(["definitely_missing_tool_xyz_123"])


def test_run_capture_stdout():
    proc = run(["echo", "hello"])
    assert proc.returncode == 0
    assert proc.stdout.strip() == "hello"


def test_run_redirect_stdout_and_capture_stderr(tmp_path):
    out = tmp_path / "stdout.txt"
    with open(out, "w") as stdout:
        proc = run(
            [
                sys.executable,
                "-c",
                "import sys; data=sys.stdin.read(); print(data.upper()); print('err', file=sys.stderr)",
            ],
            input="abc",
            stdout=stdout,
        )

    assert out.read_text().strip() == "ABC"
    assert proc.stderr.strip() == "err"


def test_run_str_shell():
    proc = run("echo hi from shell")
    assert proc.returncode == 0
    assert "hi from shell" in proc.stdout


def test_run_nonzero_check_false():
    proc = run(["false"], check=False)
    assert proc.returncode != 0


def test_run_nonzero_raises():
    with pytest.raises(ToolError) as ei:
        run(["false"], check=True)
    assert ei.value.returncode != 0


def test_run_missing_tool_raises():
    with pytest.raises(ToolError) as ei:
        run(["definitely_missing_tool_xyz_123"])
    assert ei.value.returncode == 127


def test_tool_version():
    v = tool_version("python")
    assert v is not None and v.startswith("Python")
    assert tool_version("definitely_missing_tool_xyz_123") is None


def test_toolerror_attributes():
    err = ToolError("some cmd", 3, "boom")
    assert err.cmd == "some cmd" and err.returncode == 3 and "boom" in str(err)


def test_get_logger():
    assert isinstance(get_logger("iflas.test"), logging.Logger)


def test_toolrunner_basic_and_routing():
    tr = ToolRunner({"ccs": "iflas-pacbio"})
    # no env mapping for echo -> runs locally
    proc = tr.run(["echo", "ok"])
    assert proc.stdout.strip() == "ok"
    # routing lookup resolves the basename
    assert tr._tool_of(["/some/dir/ccs", "x"]) == "ccs"
    assert tr.env_map["ccs"] == "iflas-pacbio"


def test_no_daemon_process_ignores_daemon_assignment():
    proc = NoDaemonProcess()
    assert proc.daemon is False
    proc.daemon = True
    assert proc.daemon is False


def test_my_pool_workers_can_spawn_child_pools():
    with MyPool(processes=1) as pool:
        assert pool.map(_spawn_child_pool, [21]) == [42]


@pytest.mark.skipif(which("perl") is None, reason="perl not on PATH")
def test_perl_hash_order_is_fixed_unless_the_user_sets_a_seed(monkeypatch):
    """BUG-162: the legacy Perl helpers write hash-iteration order, which Perl otherwise changes on every run."""
    keys = 'my %h = map { $_ => 1 } ("a".."z"); print join(",", keys %h)'
    assert len({run(["perl", "-e", keys]).stdout for _ in range(3)}) == 1
    monkeypatch.setenv("PERL_HASH_SEED", "7")
    assert run(["perl", "-e", "print $ENV{PERL_HASH_SEED}"]).stdout == "7"


@pytest.mark.parametrize("via_conda_run", [False, True])
def test_conda_env_commands_run_inside_that_env(monkeypatch, via_conda_run):
    """Isolated tools run with the env's bin first on PATH, whether activated directly or through `conda run`."""
    prefix = conda_env_prefix("iflas-pacbio")
    if prefix is None:
        pytest.skip("conda env iflas-pacbio not installed")
    monkeypatch.setattr(runtime, "_CONDA_PREFIXES", {})
    if via_conda_run:
        monkeypatch.setenv("IFLAS_CONDA_RUN", "1")
    proc = run(["bash", "-c", 'echo "$CONDA_DEFAULT_ENV|$(command -v python)|$1"', "_", "two words"], conda_env="iflas-pacbio")
    assert proc.stdout.strip() == "iflas-pacbio|{}|two words".format(os.path.join(prefix, "bin", "python"))
    assert runtime._CONDA_PREFIXES == {"iflas-pacbio": None if via_conda_run else prefix}


def _raise_tool_error():
    raise ToolError("tool --flag", 3, "boom")


def test_tool_error_from_a_pool_worker_reaches_the_caller():
    """BUG-167: pickling kept only the message, so unpickling a worker's ToolError failed and the pool died with
    "terminated abruptly" instead of reporting the failed command."""
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(1, mp_context=multiprocessing.get_context("fork")) as pool:
        with pytest.raises(ToolError) as exc:
            pool.submit(_raise_tool_error).result()
    assert exc.value.cmd == "tool --flag" and exc.value.returncode == 3 and "boom" in str(exc.value)
