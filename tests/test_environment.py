from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def _deps(name):
    env = yaml.safe_load((ROOT / name).read_text(encoding="utf-8"))
    return [dep for dep in env["dependencies"] if isinstance(dep, str)]


def test_core_environment_pins_pygenometracks_to_py3_compatible_release():
    assert any(dep.startswith("pygenometracks>=3.9") for dep in _deps("environment.yml"))


def test_environments_bound_the_libraries_whose_next_major_changes_behaviour():
    """A fresh install must not pick untested majors: pandas 3 (copy-on-write, no silent downcasting), numpy, scikit-learn
    (PU bagging), the R stack the report code was run with (ggplot2 4 / Bioconductor 3.22), StringTie 3 / regtools 1 CLIs."""
    pinned = {dep.split(">=")[0].split("<")[0].split("=")[0]: dep for name in ("environment.yml", "environment-diffsplice.yml")
              for dep in _deps(name)}
    for lib in ("numpy", "pandas", "scipy", "scikit-learn", "matplotlib-base", "pysam", "stringtie", "regtools", "r-base", "r-ggplot2",
                "bioconductor-clusterprofiler", "bioconductor-genomicfeatures", "rmats"):
        assert "<" in pinned[lib], pinned[lib]
    for name in ("iflas", "iflas-cupcake", "iflas-pacbio", "iflas-diffsplice"):
        lock = ROOT / "envs" / "lock" / f"{name}.linux-64.txt"
        assert "@EXPLICIT" in lock.read_text(encoding="utf-8"), lock
