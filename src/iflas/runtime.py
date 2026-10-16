"""Runtime utilities: logging + a guarded external-command runner.

Replaces the legacy pattern of scattered ``subprocess.call(cmd, shell=True)`` (163
call sites) and ``print`` (210) with:
  * list-argument commands (no shell by default),
  * return-code checking with a clear ``ToolError``,
  * stdout/stderr captured to the logger,
  * optional per-conda-env execution (for the isolated iflas-pacbio / iflas-diffsplice envs),
  * a startup preflight (``require_tools``) and version probe (``tool_version``).
"""

from __future__ import annotations

import datetime
import logging
import multiprocessing
import multiprocessing.pool
import json
import os
import shlex
import shutil
import subprocess
import tempfile
from typing import Mapping, Sequence

_DEFAULT_FMT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"


def getCurrentTime() -> str:
    """Timestamp string for stage progress messages (ported from commonFuncs)."""
    return str(datetime.datetime.now())


def get_logger(name: str = "iflas") -> logging.Logger:
    """Return a module logger, ensuring at least a basic config exists."""
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format=_DEFAULT_FMT)
    return logging.getLogger(name)


class NoDaemonProcess(multiprocessing.Process):
    """Process whose daemon flag stays False, matching legacy commonObjs.py."""

    def _get_daemon(self) -> bool:
        return False

    def _set_daemon(self, value) -> None:
        pass

    daemon = property(_get_daemon, _set_daemon)


class _NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess


class MyPool(multiprocessing.pool.Pool):
    """multiprocessing Pool whose workers may create child processes."""

    def __init__(self, *args, context=None, **kwargs):
        super().__init__(*args, context=context or _NoDaemonContext(), **kwargs)


class ToolError(RuntimeError):
    """Raised when an external command fails (non-zero exit) or is missing."""

    def __init__(self, cmd: str, returncode: int, stderr: str = ""):
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr or ""
        msg = f"command failed (exit {returncode}): {cmd}"
        if self.stderr:
            msg += "\n" + self.stderr.strip()
        super().__init__(msg)

    def __reduce__(self):
        # rebuilt from the constructor arguments when a pool worker sends it back (BUG-167: the default pickling
        # passed only the message, so unpickling failed and the pool died with "terminated abruptly")
        return ToolError, (self.cmd, self.returncode, self.stderr)


def which(tool: str) -> str | None:
    """Absolute path of *tool* on PATH, or None."""
    return shutil.which(tool)


def require_tools(tools: Sequence[str], logger: logging.Logger | None = None) -> bool:
    """Preflight: ensure each tool is on PATH; raise ToolError listing any missing."""
    logger = logger or get_logger()
    missing = [t for t in tools if shutil.which(t) is None]
    for t in tools:
        if t not in missing:
            logger.debug("tool found: %s -> %s", t, shutil.which(t))
    if missing:
        raise ToolError(" ".join(missing), 127, "missing required tools: " + ", ".join(missing))
    return True


def quote(arg) -> str:
    """``arg`` (a path, usually) as one shell word, for the commands that need a pipe, redirect or glob."""
    return shlex.quote(str(arg))


def use_tmp_dir(path: str) -> str:
    """Send temporary files to ``path``: pybedtools and Python's tempfile here, and through TMPDIR the worker processes,
    R and the external tools (legacy iflas.py set pybedtools' temp dir from ``dirs.tmp_dir``; BUG-184)."""
    path = os.path.abspath(path)
    os.makedirs(path, exist_ok=True)
    os.environ["TMPDIR"] = path
    tempfile.tempdir = path
    import pybedtools  # lazy: the light subcommands do not need it

    pybedtools.set_tempdir(path)
    return path


def tool_version(tool: str, flag: str = "--version") -> str | None:
    """First line of ``tool flag`` output, or None if the tool is unavailable."""
    if shutil.which(tool) is None:
        return None
    try:
        proc = subprocess.run([tool, flag], capture_output=True, text=True)
        out = (proc.stdout or proc.stderr or "").strip().splitlines()
        return out[0] if out else ""
    except Exception:
        return None


SINGLE_THREAD_ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"}
"""For subprocesses that iFLAS itself runs many at a time (one per locus / gene): multithreaded math libraries
would otherwise start a thread per core in every one of them and oversubscribe the machine."""

DETERMINISTIC_ENV = {"PERL_HASH_SEED": "0", "PERL_PERTURB_KEYS": "0"}
"""Given to every external command (a value the user exported wins): the legacy Perl helpers write hash-iteration
order, which Perl otherwise changes from run to run (BUG-162)."""

_CONDA_PREFIXES: dict[str, str | None] = {}


def conda_env_prefix(name: str) -> str | None:
    """Install prefix of the named conda environment, resolved once per process; None when it cannot be resolved
    (or when ``IFLAS_CONDA_RUN=1`` asks for plain `conda run`)."""
    if name not in _CONDA_PREFIXES:
        prefix = None
        if os.environ.get("IFLAS_CONDA_RUN") != "1":
            try:
                out = subprocess.run(["conda", "env", "list", "--json"], capture_output=True, text=True, check=True).stdout
                prefix = next((p for p in json.loads(out).get("envs", []) if os.path.basename(p) == name), None)
            except (OSError, subprocess.CalledProcessError, ValueError):
                prefix = None
        _CONDA_PREFIXES[name] = prefix
    return _CONDA_PREFIXES[name]


def _activated(cmd, name: str, prefix: str, env: Mapping[str, str]):
    """``cmd`` run inside conda env ``name`` the way `conda activate` sets it up -- the current env's deactivate.d
    scripts, PATH and CONDA_PREFIX switched, the target's activate.d scripts -- without `conda run`'s ~1.5 s
    Python start-up per call (which dominated per-locus IsoPhase runs). ``env`` is the complete child environment."""
    run_env = dict(env)
    parent = run_env.get("CONDA_PREFIX")
    hooks = []
    if parent and os.path.realpath(parent) != os.path.realpath(prefix):
        hooks.append('for f in %s/etc/conda/deactivate.d/*.sh; do [ -r "$f" ] && . "$f"; done' % shlex.quote(parent))
        run_env["PATH"] = os.pathsep.join(p for p in run_env.get("PATH", "").split(os.pathsep) if p != os.path.join(parent, "bin"))
    run_env["PATH"] = os.path.join(prefix, "bin") + os.pathsep + run_env.get("PATH", "")
    hooks.append("export CONDA_PREFIX=%s CONDA_DEFAULT_ENV=%s" % (shlex.quote(prefix), shlex.quote(name)))
    hooks.append('for f in "$CONDA_PREFIX"/etc/conda/activate.d/*.sh; do [ -r "$f" ] && . "$f"; done')
    prologue = "; ".join(hooks) + "; "
    if isinstance(cmd, str):
        return ["bash", "-c", prologue + cmd], run_env
    return ["bash", "-c", prologue + 'exec "$@"', "bash", *map(str, cmd)], run_env


def run(cmd, *, cwd=None, env: Mapping[str, str] | None = None, check: bool = True,
        capture: bool = True, shell: bool | None = None, conda_env: str | None = None,
        logger: logging.Logger | None = None, input=None, stdout=None, stderr=None,
        text: bool = True) -> subprocess.CompletedProcess:
    """Run an external command with logging and return-code checking.

    ``cmd``: list[str] (preferred; no shell) or str (runs via shell).
    ``conda_env``: if set, execute inside that conda env (activated directly; ``conda run`` when the env's
    prefix cannot be resolved or ``IFLAS_CONDA_RUN=1``).
    Returns the CompletedProcess; raises :class:`ToolError` on failure when ``check``.
    """
    logger = logger or get_logger()
    if shell is None:
        shell = isinstance(cmd, str)

    run_env = {**DETERMINISTIC_ENV, **os.environ, **(env or {})}
    printable = cmd if isinstance(cmd, str) else " ".join(map(str, cmd))
    if conda_env:
        prefix = conda_env_prefix(conda_env)
        if prefix:
            printable = "[%s] %s" % (conda_env, printable)
            cmd, run_env = _activated(cmd, conda_env, prefix, run_env)
        else:
            cmd = (["conda", "run", "-n", conda_env, "bash", "-c", cmd] if isinstance(cmd, str)
                   else ["conda", "run", "-n", conda_env, *map(str, cmd)])
            printable = " ".join(cmd)
        shell = False

    logger.info("RUN: %s", printable)
    run_kwargs = {
        "cwd": cwd,
        "env": run_env,
        "shell": shell,
        "text": text,
    }
    if input is not None:
        run_kwargs["input"] = input
    if capture:
        stdout = subprocess.PIPE if stdout is None else stdout
        stderr = subprocess.PIPE if stderr is None else stderr
    if stdout is not None:
        run_kwargs["stdout"] = stdout
    if stderr is not None:
        run_kwargs["stderr"] = stderr

    try:
        proc = subprocess.run(cmd, **run_kwargs)
    except FileNotFoundError as e:
        raise ToolError(printable, 127, str(e))

    if capture:
        if isinstance(proc.stdout, str) and proc.stdout:
            logger.debug("stdout: %s", proc.stdout.strip())
        if isinstance(proc.stderr, str) and proc.stderr:
            logger.debug("stderr: %s", proc.stderr.strip())
    if check and proc.returncode != 0:
        raise ToolError(printable, proc.returncode, proc.stderr if isinstance(proc.stderr, str) else "")
    return proc


class ToolRunner:
    """Runs external tools, routing isolated tools to their conda env.

    Example::

        tr = ToolRunner({"ccs": "iflas-pacbio", "isoseq": "iflas-pacbio",
                         "rmats.py": "iflas-diffsplice"})
        tr.run(["ccs", "in.bam", "out.bam"])   # auto -> conda run -n iflas-pacbio ccs ...
    """

    def __init__(self, env_map: Mapping[str, str] | None = None, logger: logging.Logger | None = None):
        self.env_map = dict(env_map or {})
        self.logger = logger or get_logger()

    def _tool_of(self, cmd) -> str:
        first = cmd.split()[0] if isinstance(cmd, str) else str(cmd[0])
        return os.path.basename(first)

    def run(self, cmd, **kwargs) -> subprocess.CompletedProcess:
        conda_env = kwargs.pop("conda_env", None) or self.env_map.get(self._tool_of(cmd))
        kwargs.setdefault("logger", self.logger)
        return run(cmd, conda_env=conda_env, **kwargs)

    def require(self, tools: Sequence[str]) -> bool:
        return require_tools(tools, logger=self.logger)
