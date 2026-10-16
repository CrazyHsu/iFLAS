"""The user documentation names the commands that exist and the tool versions that are installed."""
import re
from pathlib import Path

from iflas.cli import build_parser

ROOT = Path(__file__).resolve().parents[1]


def _subcommands():
    action = next(a for a in build_parser()._actions if a.dest == "command")
    return set(action.choices or ())


def test_readme_shows_every_command_and_the_cupcake_environment():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    for command in _subcommands():
        assert re.search(rf"^\s+{re.escape(command)}\b", text, re.MULTILINE), command
    for example in ("iflas all", "iflas mapping", "iflas collapse", "iflas refine", "iflas pu_filter",
                    "iflas find_as", "iflas report -all"):
        assert example in text
    assert "conda env create -f environment-cupcake.yml" in text


def _locked(lock, package):
    m = re.search(rf"/{re.escape(package)}-([0-9][^-/]*)-", (ROOT / "envs" / "lock" / lock).read_text())
    assert m, (lock, package)
    return m.group(1)


def test_external_tools_page_lists_the_locked_versions():
    text = (ROOT / "docs" / "external_tools.md").read_text(encoding="utf-8")
    for lock, package, name in [
        ("iflas.linux-64.txt", "minimap2", "minimap2"),
        ("iflas.linux-64.txt", "hisat2", "HISAT2"),
        ("iflas.linux-64.txt", "stringtie", "StringTie"),
        ("iflas.linux-64.txt", "regtools", "regtools"),
        ("iflas.linux-64.txt", "nanopolish", "nanopolish"),
        ("iflas.linux-64.txt", "pygenometracks", "pyGenomeTracks"),
        ("iflas-pacbio.linux-64.txt", "pbccs", "ccs"),
        ("iflas-pacbio.linux-64.txt", "lima", "lima"),
        ("iflas-pacbio.linux-64.txt", "isoseq", "isoseq"),
        ("iflas-diffsplice.linux-64.txt", "rmats", "rMATS"),
        ("iflas-diffsplice.linux-64.txt", "r-base", "R"),
        ("iflas-diffsplice.linux-64.txt", "bioconductor-clusterprofiler", "clusterProfiler"),
        ("iflas-diffsplice.linux-64.txt", "bioconductor-gviz", "Gviz"),
    ]:
        assert f"{name} {_locked(lock, package)}" in text, (name, _locked(lock, package))
