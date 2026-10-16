"""The Perl helpers in utils/ must not print import warnings: find_as uses some scripts' stderr as a data channel (R-8)."""
import os
import shutil
import subprocess

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
UTILS = os.path.join(ROOT, "utils")


@pytest.mark.skipif(shutil.which("perl") is None, reason="needs perl")
@pytest.mark.parametrize("script", ["filter.pl", "select.pl"])
def test_pm_common_users_emit_no_import_warning(tmp_path, script):
    """BUG-146: `use pm::common qw/.../` on a module without Exporter printed
    'Attempt to call missing import method ... via package "pm::common"' on every run."""
    origin = tmp_path / "ids.txt"
    target = tmp_path / "target.tsv"
    origin.write_text("a\n")
    target.write_text("a\t1\nb\t2\n")
    if script == "filter.pl":
        cmd = ["perl", os.path.join(UTILS, script), "-o", str(origin), str(target), "-1", "1", "-2", "1", "-m", "i"]
    else:
        cmd = ["perl", os.path.join(UTILS, script), "-i", "1", str(target)]
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=UTILS)
    assert "Attempt to call missing import method" not in proc.stderr
    assert "pm::common" not in proc.stderr
