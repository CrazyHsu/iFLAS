"""Release metadata checks."""

from pathlib import Path
import tomllib

import iflas


ROOT = Path(__file__).resolve().parents[1]
RELEASE_VERSION = "2.0.0b1"
RELEASE_DATE = "2026-10-16"


def test_release_version_matches_package_metadata():
    with (ROOT / "pyproject.toml").open("rb") as fh:
        project = tomllib.load(fh)["project"]

    assert project["version"] == RELEASE_VERSION
    assert iflas.__version__ == RELEASE_VERSION


def test_changelog_has_current_release_section():
    changelog = (ROOT / "CHANGELOG.md").read_text()

    assert f"## [{RELEASE_VERSION}] - {RELEASE_DATE}" in changelog
