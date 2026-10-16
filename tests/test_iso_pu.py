"""Smoke tests for the ported PU-learning core (src/iflas/pu/iso_pu.py).

Verifies the package-internal port (iflas.pu.iso_pu → iflas.pu.bagging / iflas.pu.cpc2)
imports and that SchemeModelEval's bagging scheme fits + scores under sklearn 1.9.
NOTE (route B): verifies the pipeline RUNS; numerical equivalence to the published
Python-2 result still requires the real golden baseline. The tool-heavy feature
extraction (getFeatures/predCodingP/addJuncIndelInfo) and the stage driver (stages/pu_filter.py) are covered
by tests/test_pu_filter.py and the synthetic-v2 chain.
"""

import os

import numpy as np
import pandas as pd

import iflas.pu.iso_pu as iso_pu_module
from iflas.io.records import BedFile
from iflas.pu.iso_pu import (GetFeatures, SchemeModelEval, getInputBed, processInputFile, selectBestModel)


def _synthetic_features(n_pos=20, n_unl=20, seed=0):
    """A minimal PU training frame: numeric features + trailing binary `label`."""
    rng = np.random.RandomState(seed)
    pos = rng.normal(1.0, 0.3, size=(n_pos, 5))
    unl = rng.normal(0.0, 0.3, size=(n_unl, 5))
    X = np.vstack([pos, unl])
    df = pd.DataFrame(X, columns=["f{}".format(i) for i in range(5)],
                      index=["iso{}".format(i) for i in range(n_pos + n_unl)])
    df["label"] = [1] * n_pos + [0] * n_unl
    return df


def test_port_exposes_public_api():
    # the line-by-line port keeps the original entry points
    for fn in (SchemeModelEval, GetFeatures, getInputBed, processInputFile, selectBestModel):
        assert callable(fn)


def test_scheme_model_eval_bagging_runs():
    data = _synthetic_features()
    sme = SchemeModelEval(trainingData=data, model="DT", scheme="bagging")
    sme.eval()

    res = sme.predResults
    assert list(res.columns) == ["label", "pu_score"]
    assert len(res) == len(data)
    finite = res["pu_score"][np.isfinite(res["pu_score"])]
    assert finite.size > 0
    assert finite.min() >= 0.0 and finite.max() <= 1.0


def test_scheme_model_eval_can_limit_bagging_jobs():
    data = _synthetic_features()
    sme = SchemeModelEval(
        trainingData=data, model="DT", scheme="bagging", bagging_n_jobs=1
    )
    sme.eval()

    assert sme.finalEstimator.n_jobs == 1


def test_filter_isoforms_by_score(tmp_path):
    data = _synthetic_features(n_pos=12, n_unl=12)
    sme = SchemeModelEval(trainingData=data, model="DT", scheme="bagging")
    sme.eval()

    out = tmp_path / "validIsoforms.lst"
    # filterScore=0.0 keeps every annotated (label==1) + every novel (label==0, score>=0)
    sme.filterIsoformsByScore(filterScore=0.0, outFile=str(out))
    assert out.exists()

    written = pd.read_csv(out, sep="\t")
    assert "pu_score" in written.columns
    assert len(written) == len(data)


def test_getfeatures_input_fasta_pipeline_uses_runner(tmp_path):
    input_bed = tmp_path / "input.bed12+"
    genome = tmp_path / "genome.fa"
    input_bed.write_text("")
    genome.write_text(">chr1\nACGT\n")

    class RecordingRunner:
        def __init__(self):
            self.calls = []

        def run(self, cmd, **kwargs):
            self.calls.append((cmd, kwargs, os.getcwd()))
            (tmp_path / "tmpInput.fa").write_text(">iso1\nACGT\n")

    runner = RecordingRunner()
    prev = os.getcwd()
    os.chdir(tmp_path)
    try:
        gf = GetFeatures(inputBed=str(input_bed), genomeFa=str(genome), runner=runner)
        out = gf.getInputFa()
    finally:
        os.chdir(prev)

    assert out == str(tmp_path / "tmpInput.fa")
    assert (tmp_path / "tmpInput.fa").is_file()
    assert len(runner.calls) == 1
    cmd, kwargs, cwd = runner.calls[0]
    assert cwd == str(tmp_path)
    assert kwargs == {}
    assert cmd[:2] == ["bash", "-c"]
    pipeline = cmd[2]
    assert str(input_bed) in pipeline
    assert str(genome) in pipeline
    assert "bedtools getfasta" in pipeline
    assert "seqkit replace" in pipeline


def test_add_junc_indel_info_uses_true_division_for_the_indel_ratio(tmp_path, monkeypatch):
    """1 indel-bearing read of 2 junction reads = 0.5 >= 0.5 flags the junction (BUG-140 fixed at CP3; Py2 legacy
    floored 1 // 2 to 0 so the threshold was unreachable)."""
    input_bed = tmp_path / "input.bed12+"
    junction_bed = tmp_path / "junctions.bed"
    ref_bed = tmp_path / "ref.bed12+"
    bed12 = "chr1\t100\t300\tiso1\t2\t+\t100\t300\t0,0,0\t2\t50,50\t0,150"
    input_bed.write_text(f"{bed12}\tgene1\n")
    junction_bed.write_text(f"{bed12}\n")
    ref_bed.write_text(f"{bed12}\tgene1\n")

    class FakeRead:
        def __init__(self, cigar):
            self.is_unmapped = False
            self.cigar = cigar
            self.pos = 100
            self.reference_name = "chr1"
            self.reference_id = 0

    class FakeAlignment:
        def fetch(self):
            return iter([
                FakeRead([(0, 45), (1, 3), (0, 5), (3, 100), (0, 50)]),
                FakeRead([(0, 50), (3, 100), (0, 50)]),
            ])

        def close(self):
            pass

    monkeypatch.setattr(
        iso_pu_module.pysam, "AlignmentFile", lambda *args, **kwargs: FakeAlignment()
    )

    gf = GetFeatures(inputBed=str(input_bed))
    gf.inputBedObj = BedFile(str(input_bed), type="bed12+")

    features = gf.addJuncIndelInfo(
        samFile=str(tmp_path / "reads.sam"),
        junctionFile=str(junction_bed),
        refBed=str(ref_bed),
    )

    assert features.loc["iso1", "nIndelsAroundJunc"] == 1
    assert features.loc["iso1", "nJuncsWithIndels"] == 1
    assert features.loc["iso1", "indelNearJunc"] == 1
