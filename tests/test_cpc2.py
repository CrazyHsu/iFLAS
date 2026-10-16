"""Verify CPC2 (coding-potential) works under Python 3 (prerequisite for pu_filter).

Two things are checked:
  1. FindCDS.longest_orf — pure-Python ORF finder; characterizes the already-working
     coding-region search on a sequence with a single, known ORF.
  2. calculate_potential — the full LIBSVM pipeline (Fickett + pI + ORF features →
     utils/svm-scale | utils/svm-predict → cpc2out.txt). This is the path iso_pu's
     predCodingP drives. It crashes under Py3 until the Py2 `file()` builtin is replaced.

The end-to-end test is gated on the LIBSVM binaries being present + executable.
NOTE (route B): this verifies CPC2 RUNS correctly under Py3 / modern Biopython; numerical
equivalence to the published Python-2 result still requires the real golden baseline.
"""

import os
import io
import subprocess

import pytest

import iflas.pu.cpc2 as cpc2
from iflas.pu.cpc2 import FindCDS, calculate_potential

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")
_UTILS = os.path.join(REPO_ROOT, "utils")
_HAS_SVM = all(
    os.access(os.path.join(_UTILS, tool), os.X_OK)
    for tool in ("svm-scale", "svm-predict")
) and all(
    os.path.isfile(os.path.join(_UTILS, f)) for f in ("cpc2.range", "cpc2.model")
)

# A sequence carrying exactly one ORF: ATG AAA TTT GGG TAA (15 nt), flanked by non-coding ends.
_KNOWN_ORF_SEQ = "TT" + "ATGAAATTTGGGTAA" + "CC"

# Two ~190/120 nt transcripts (valid ACGT) used only to exercise the LIBSVM pipeline
# structurally — the coding/noncoding call itself is not asserted.
_ISO1 = "GCGC" + "ATG" + "GCTGGCAAGCTGGAA" * 12 + "TAA" + "GCGC"
_ISO2 = "ACGT" * 20 + "ATGCCCGGGAAATTT" + "ACGTACGT"


def test_findcds_longest_orf_finds_known_orf():
    orf, start, strand, integrity = FindCDS(_KNOWN_ORF_SEQ).longest_orf("+")
    assert orf == "ATGAAATTTGGGTAA"
    assert start == 2          # 0-based index of the ATG
    assert strand == "+"
    assert integrity == 1      # stop codon reached → complete ORF


@pytest.mark.skipif(not _HAS_SVM, reason="LIBSVM binaries / model not available in utils/")
def test_calculate_potential_end_to_end(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">iso1\n{}\n>iso2\n{}\n".format(_ISO1, _ISO2))
    outfile = str(tmp_path / "cpc2out")

    rc = calculate_potential(str(fasta), "+", 0, outfile)
    assert rc == 0, "calculate_potential did not complete the LIBSVM prediction"

    result = tmp_path / "cpc2out.txt"
    assert result.exists(), "expected cpc2out.txt to be written"

    lines = result.read_text().rstrip("\n").split("\n")
    header, rows = lines[0], lines[1:]
    assert header.startswith("#ID")
    assert len(rows) == 2, "expected one row per input transcript"
    for row in rows:
        cols = row.split("\t")
        # 8 columns: #ID, transcript_length, peptide_length, Fickett, pI,
        #            ORF_integrity, coding_probability, label
        assert len(cols) == 8
        coding_prob = float(cols[6])
        assert 0.0 <= coding_prob <= 1.0
        assert cols[7] in ("coding", "noncoding")


def test_calculate_potential_routes_libsvm_through_toolrunner(tmp_path, monkeypatch):
    utils = tmp_path / "utils"
    utils.mkdir()
    for name in ("svm-scale", "svm-predict"):
        exe = utils / name
        exe.write_text("#!/bin/sh\nexit 0\n")
        exe.chmod(0o755)
    (utils / "cpc2.range").write_text("range\n")
    (utils / "cpc2.model").write_text("model\n")
    monkeypatch.setenv("IFLAS_UTILS_DIR", str(utils))

    class FakeRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs))
            tool = os.path.basename(cmd[0])
            if tool == "svm-scale":
                kwargs["stdout"].write("scaled\n")
            elif tool == "svm-predict":
                with open(cmd[-1], "w") as out:
                    out.write("labels 1 -1\n1 0.7 0.3\n-1 0.2 0.8\n")
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    fasta = tmp_path / "input.fa"
    fasta.write_text(">iso1\n{}\n>iso2\n{}\n".format(_ISO1, _ISO2))
    outfile = str(tmp_path / "cpc2out")
    runner = FakeRunner()

    assert calculate_potential(str(fasta), "+", 0, outfile, runner=runner) == 0

    assert [os.path.basename(call[0][0]) for call in runner.calls] == ["svm-scale", "svm-predict"]
    assert runner.calls[0][0][1:] == ["-r", str(utils / "cpc2.range"), outfile + ".tmp.1"]
    assert runner.calls[1][0][1:] == [
        "-b", "1", "-q", outfile + ".tmp.2", str(utils / "cpc2.model"), outfile + ".tmp.out",
    ]

    rows = (tmp_path / "cpc2out.txt").read_text().rstrip("\n").split("\n")
    assert rows[1].endswith("\t0.7\tcoding")
    assert rows[2].endswith("\t0.2\tnoncoding")
    for suffix in (".tmp.1", ".tmp.2", ".tmp.out", ".feat", ""):
        assert not os.path.exists(outfile + suffix)


def test_calculate_potential_closes_feature_handles_when_fasta_read_fails(monkeypatch):
    handles = []

    class TrackingHandle(io.StringIO):
        pass

    def fake_open(path, mode="r"):
        handle = TrackingHandle()
        handles.append(handle)
        return handle

    def broken_fasta_read(path):
        raise RuntimeError("broken FASTA")

    monkeypatch.setattr(cpc2, "open", fake_open, raising=False)
    monkeypatch.setattr(cpc2, "fasta_read", broken_fasta_read)

    with pytest.raises(RuntimeError, match="broken FASTA"):
        cpc2.calculate_potential("input.fa", "+", 0, "cpc2out")

    assert len(handles) == 3
    assert all(handle.closed for handle in handles)


def test_fasta_read_closes_handle_when_parser_fails(monkeypatch):
    handles = []

    class TrackingHandle:
        closed = False

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            self.close()

        def close(self):
            self.closed = True

    def fake_gz_file(path, mode):
        handle = TrackingHandle()
        handles.append(handle)
        return handle

    def broken_parse(handle, fmt):
        raise RuntimeError("parse failed")

    monkeypatch.setattr(cpc2, "gz_file", fake_gz_file)
    monkeypatch.setattr(cpc2.SeqIO, "parse", broken_parse)

    with pytest.raises(RuntimeError, match="parse failed"):
        list(cpc2.fasta_read("input.fa"))

    assert len(handles) == 1
    assert handles[0].closed


@pytest.mark.skipif(not _HAS_SVM, reason="LIBSVM binaries / model not available in utils/")
def test_predcoding_read_path(tmp_path):
    """Replicate iso_pu.predCodingP's consumer read, with the pandas-2 set_axis form."""
    import pandas as pd

    fasta = tmp_path / "input.fa"
    fasta.write_text(">iso1\n{}\n>iso2\n{}\n".format(_ISO1, _ISO2))
    outfile = str(tmp_path / "cpc2out")
    assert calculate_potential(str(fasta), "+", 0, outfile) == 0

    coding = pd.read_csv(
        tmp_path / "cpc2out.txt", sep="\t", skiprows=1, header=None,
        usecols=[0, 2, 3, 4, 5, 6, 7], index_col=0,
    )
    coding = coding.set_axis(
        ["pepLength", "FickettScore", "pI", "orfIntegrity", "codingP", "codingLabel"], axis=1
    )
    assert list(coding.columns) == ["pepLength", "FickettScore", "pI", "orfIntegrity", "codingP", "codingLabel"]
    assert list(coding.index) == ["iso1", "iso2"]
    assert coding["codingP"].between(0.0, 1.0).all()
