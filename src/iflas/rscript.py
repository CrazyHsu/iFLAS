"""Helpers for running legacy R plotting functions through ``Rscript``.

The published pipeline embedded R source strings through rpy2.  The migration keeps
those R functions intact but writes a short temporary ``.R`` file and executes it
with the isolated R/Bioconductor environment.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from numbers import Integral, Real
from pathlib import Path
from typing import Iterable

from .runtime import ToolRunner

DEFAULT_R_CONDA_ENV = "iflas-diffsplice"


def default_r_conda_env() -> str | None:
    """Return the conda env used for Rscript, allowing local execution override."""
    val = os.environ.get("IFLAS_R_CONDA_ENV")
    if val is None:
        return DEFAULT_R_CONDA_ENV
    val = val.strip()
    if val.lower() in {"", "none", "local"}:
        return None
    return val


def r_literal(value) -> str:
    """Render a Python value as a small R literal for function-call scripts."""
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, Integral):
        return str(value)
    if isinstance(value, Real):
        if not math.isfinite(float(value)):
            raise ValueError(f"non-finite numeric value cannot be passed to R: {value!r}")
        return repr(float(value))
    if isinstance(value, os.PathLike):
        value = os.fspath(value)
    return json.dumps(str(value), ensure_ascii=False)


def build_r_script(function_source: str, function_name: str, args: Iterable = ()) -> str:
    """Compose an R script that defines *function_source* and calls *function_name*."""
    return build_r_calls(function_source, function_name, [args])


def build_r_calls(function_source: str, function_name: str, calls: Iterable[Iterable]) -> str:
    """An R script that defines *function_source* once and calls *function_name* with each argument list."""
    body = "".join(f"{function_name}({', '.join(r_literal(arg) for arg in args)})\n" for args in calls)
    return f"{function_source.rstrip()}\n\n{body}"


def run_r_function(function_source: str, function_name: str, args: Iterable = (), **kwargs):
    """Run an R function by writing a temporary script and invoking ``Rscript``."""
    return run_r_calls(function_source, function_name, [args], **kwargs)


def run_r_calls(
    function_source: str,
    function_name: str,
    calls: Iterable[Iterable],
    *,
    cwd: str | os.PathLike | None = None,
    runner: ToolRunner | None = None,
    conda_env: str | None = None,
    keep_script: bool = False,
):
    """Call an R function once per argument list in a single ``Rscript`` session: packages and the source (e.g. a
    reference TxDb) load once, not once per plot (a report plot per event cost about 20 s of R start-up)."""
    workdir = Path(cwd or os.getcwd())
    script_text = build_r_calls(function_source, function_name, calls)
    fd, script_path = tempfile.mkstemp(
        prefix=f"iflas_{function_name}_",
        suffix=".R",
        dir=str(workdir),
        text=True,
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(script_text)
        runner = runner or ToolRunner()
        env_name = default_r_conda_env() if conda_env is None else conda_env
        return runner.run(["Rscript", script_path], cwd=str(workdir), conda_env=env_name)
    finally:
        if not keep_script:
            try:
                os.unlink(script_path)
            except FileNotFoundError:
                pass


def run_r_script(cmd: str, *, cwd=None, runner: ToolRunner | None = None, conda_env: str | None = None,
                 fatal: bool = False) -> bool:
    """Run a shell command that drives a ``utils/*.R`` script (``Rscript path/x.R args``) in the R env.

    Plots are not scientific products: a failure is logged as a warning (legacy ignored the exit code
    silently) unless ``fatal``. Returns True on success. Brace expansion / globs in ``cmd`` need bash, which
    ``conda run`` provides; without a conda env the command is wrapped in ``bash -c`` explicitly."""
    runner = runner or ToolRunner()
    env_name = default_r_conda_env() if conda_env is None else conda_env
    proc = runner.run(cmd if env_name else ["bash", "-c", cmd], conda_env=env_name,
                      cwd=str(cwd) if cwd else None, check=fatal)
    if proc.returncode:
        runner.logger.warning("R plot failed (rc=%s): %s", proc.returncode, cmd)
    return proc.returncode == 0
