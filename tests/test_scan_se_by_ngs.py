"""Regression tests for the NGS SE scanner helpers (scan_se_by_ngs.py)."""

from iflas.stages.scan_se_by_ngs import assignGeneNameToSE

_SE_ROW = ("chr1\t118340\t118500\tnullA:117890@118340-118500@118850\t527\t+\t118340\t118500\t0,0,0\t1\t160\t0"
           "\t19\t36\tchr1:117890-118850\n")


def test_assign_gene_name_to_se_copies_known_rows_verbatim(tmp_path):
    """BUG-149: Py2 `print >>out, line,` wrote each (newline-terminated) row as is; the 2to3 form
    `print(line, end=' ')` prefixed every following row with a space and left a dangling " " line, which
    broke the characterize NGS/LR comparisons (keys like " chr1:...") and getAnnoASList (IndexError)."""
    ref = tmp_path / "ref.gpe"
    ref.write_text("")                                   # no reference genes: nothing to rename
    fes = tmp_path / "SEsupportByNGS.bed12+"
    fes.write_text("#header\n" + _SE_ROW + _SE_ROW.replace("nullA", "nullB"))
    out, err = tmp_path / "SE.bed12+", tmp_path / "novelSE2gene.log"

    assignGeneNameToSE(str(ref), False, str(fes), str(out), str(err))

    assert out.read_text() == fes.read_text()
    assert err.read_text() == ""
