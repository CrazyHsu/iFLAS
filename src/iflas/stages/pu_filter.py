"""pu_filter stage — high-quality novel-isoform identification via PU learning.

The PU core (SchemeModelEval / GetFeatures / getInputBed / processInputFile) is the
line-by-line port in ``iflas.pu.iso_pu`` (algorithm unchanged).
This module only adds a path-based DRIVER, ``filter_hq_isoforms``, that reproduces the
LR core of legacy ``iso_pu.iso_pu1`` with explicit paths (analog to refine_isoforms /
find_as_lr): it wires getInputBed -> processInputFile (features incl. CPC2 coding
potential + NGS junction/indel) -> SchemeModelEval(bagging) -> filterIsoformsByScore
-> ``utils/filter.pl`` and writes ``hqIsoforms/hq.collapsed.bed12+`` — the file the
find_as stage consumes via ``find_as_sample``.

Inputs expected under ``base_dir`` (produced by the upstream stages):
  refine/tofu.collapsed.assigned.unambi.bed12+   (refine)
  collapse/tofu.collapsed.group.txt              (collapse)
  mapping/rna-seq/reassembly/{tmp.bam,junctions.bed}  (align NGS reassembly; auto-detected)

NGS is effectively required: processInputFile's fixed column order includes the junction/
indel features, and the positive training set needs ``minJuncRPKM`` (from the NGS junctions).
"""

from __future__ import annotations

import os
import shutil

import numpy as np
import pandas as pd

from ..io.files import resolveDir
from ..pu.cpc2 import _find_utils_dir
from ..pu.iso_pu import (SchemeModelEval, getInputBed, processInputFile, selectBestModel)
from ..runtime import ToolRunner, getCurrentTime, quote as q

# feature columns fed to the PU classifier (identical to legacy iso_pu1 usedFeatures)
USED_FEATURES = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength",
                 "orfIntegrity", "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio",
                 "sdJuncCov", "minJuncRPKM", "nIndelsAroundJunc", "ratioMinJuncCovToAllCov",
                 "nJuncsWithIndels", "indelNearJunc", "label"]


MIN_SELECTION_SET = 5      # 5-fold cross-validation on the positives, 80/20 split of the reliable negatives


def _reliable_negatives(novelIsoData) -> int:
    """Novel isoforms selectBestModel treats as reliable negatives: < 5 % of their gene's reads, or an unsupported
    novel junction (RPKM < 0.05)."""
    weak = novelIsoData.ratioIsoToGene < 0.05
    if "withNovelJunc" in novelIsoData.columns:
        weak |= (novelIsoData.withNovelJunc == True) & (novelIsoData.minNovelJuncRPKM < 0.05)
    return int(weak.sum())


def write_isoform_score(featureData, out_file) -> str:
    """``isoIdentity/isoformScore.txt`` for ``iflas report --novel-hq-as`` (F-9: the published code plotted this file
    but nothing wrote it): per isoform its cupcake gene, FL read count, the gene's FL read total, the fraction
    (``ratioIsoToGene``) and Annotated/Novel, without header, in the isoform-feature order."""
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    genes = pd.Series([".".join(str(i).split(".")[:2]) for i in featureData.index], index=featureData.index)
    totals = featureData["flCount"].groupby(genes).transform("sum")
    with open(out_file, "w") as out:
        for iso, row in featureData.iterrows():
            out.write("\t".join(map(str, [genes[iso], iso, int(row["flCount"]), int(totals[iso]), row["ratioIsoToGene"],
                                          str(row["annotation"]).capitalize()])) + "\n")
    return out_file


def filter_hq_isoforms(base_dir, genome_fa, ref_bed, *, sam_file=None, junction_file=None,
                       pos_fl_coverage=2, pos_min_junc_rpkm=0.05, filter_score=0.5,
                       model="GB", out_dir=None, feature_file=None, bagging_n_jobs=None,
                       runner=None, random_state=0, select_best_model=False, auto_filter_score=False, draw_auc=False):
    """Reproduce iso_pu1's LR core with explicit paths; returns the hq.collapsed.bed12+ path.

    base_dir holds refine/ + collapse/ (+ mapping/ NGS). NGS ``tmp.bam`` / ``junctions.bed``
    are auto-detected under ``base_dir/mapping/rna-seq/reassembly`` when not given.
    ``select_best_model`` picks the classifier with the best cross-validated / held-out AUC (legacy
    ``-select_best_model``), ``auto_filter_score`` then takes the 95th percentile of the negatives' PU scores as
    the threshold (``-auto_filter_score``), ``draw_auc`` writes ``PU.bagging.ROC_AUC.pdf`` (``-draw_auc``).
    """
    print(getCurrentTime() + " Start finding high quality novel isoforms ...")
    prevDir = os.getcwd()
    # legacy cwd semantics: the stage chdirs into hqIsoforms/, so every input path must be absolute first (BUG-147)
    base_dir, genome_fa, ref_bed = os.path.abspath(base_dir), os.path.abspath(genome_fa), os.path.abspath(ref_bed)
    sam_file = os.path.abspath(sam_file) if sam_file else None
    junction_file = os.path.abspath(junction_file) if junction_file else None
    feature_file = os.path.abspath(feature_file) if feature_file else None
    hqIsoformDir = os.path.join(os.path.abspath(out_dir) if out_dir else base_dir, "hqIsoforms")
    utilDir = _find_utils_dir()
    resolveDir(hqIsoformDir)
    runner = runner or ToolRunner()

    try:
        reassembly = os.path.join(base_dir, "mapping", "rna-seq", "reassembly")
        if sam_file is None and os.path.exists(os.path.join(reassembly, "tmp.bam")):
            sam_file = os.path.join(reassembly, "tmp.bam")
        if junction_file is None and os.path.exists(os.path.join(reassembly, "junctions.bed")):
            junction_file = os.path.join(reassembly, "junctions.bed")

        inputBed = getInputBed(baseDir=base_dir, refBed=ref_bed)
        if feature_file:
            if not os.path.isfile(feature_file):
                raise FileNotFoundError(feature_file)
            local_features = os.path.join(hqIsoformDir, "isoFeatures.txt")
            if os.path.abspath(feature_file) != os.path.abspath(local_features):
                shutil.copyfile(feature_file, local_features)
            featureData = pd.read_csv(local_features, sep="\t", index_col=0)
        else:
            featureData = processInputFile(inputBed, genome_fa, ref_bed,
                                           junctionFile=junction_file, samFile=sam_file,
                                           runner=runner)

        # --- PU training set (verbatim from iso_pu1) ---
        annoIsoData = featureData.loc[featureData.annotation == "annotated", ]
        novelIsoData = featureData.loc[featureData.annotation == "novel", ]
        posIsoData = annoIsoData.loc[(annoIsoData.flCount >= int(pos_fl_coverage)) &
                                     (annoIsoData.minJuncRPKM >= float(pos_min_junc_rpkm)), ]
        lqAnnoIsoData = annoIsoData.loc[~annoIsoData.index.isin(posIsoData.index), ]

        posIsoDataCopy = posIsoData.copy()
        posIsoDataCopy["label"] = 1
        lqAnnoIsoDataCopy = lqAnnoIsoData.copy()
        lqAnnoIsoDataCopy["label"] = 1
        novelIsoDataCopy = novelIsoData.copy()
        novelIsoDataCopy["label"] = 0

        dataToClassify = pd.concat([posIsoDataCopy, novelIsoDataCopy])
        dataToClassify = dataToClassify.loc[:, USED_FEATURES]
        dataToClassify = dataToClassify.sample(frac=1, random_state=random_state)
        # booleans as 0/1 (object dtype when NaN is present) without pandas' deprecated replace-downcast
        dataToClassify = dataToClassify.apply(lambda col: col.map({True: 1, False: 0}) if col.dtype in (bool, object) else col)

        if select_best_model:
            negatives, positives = _reliable_negatives(novelIsoData), len(posIsoData)
            if negatives < MIN_SELECTION_SET or positives < MIN_SELECTION_SET:
                # selectBestModel splits the reliable negatives 80/20 and cross-validates the positives 5-fold
                print(getCurrentTime() + " WARNING: select_best_model needs >= {} reliable negatives and positives "
                      "(found {} / {}); using model {} and filter_score {}".format(MIN_SELECTION_SET, negatives, positives, model, filter_score))
            else:
                model, thresholds = selectBestModel(featureData, drawAUC=draw_auc, selectPUscore=auto_filter_score, random_state=random_state)
                if auto_filter_score:
                    filter_score = float(np.mean(thresholds[model]))
                print(getCurrentTime() + " PU model selected: {} (filter_score {:.3f})".format(model, float(filter_score)))
        sme_kwargs = {"trainingData": dataToClassify, "model": model, "scheme": "bagging", "random_state": random_state}
        if bagging_n_jobs is not None:
            sme_kwargs["bagging_n_jobs"] = bagging_n_jobs
        sme = SchemeModelEval(**sme_kwargs)
        sme.eval()
        sme.predResults.to_csv("pu_score.txt", sep="\t")
        sme.filterIsoformsByScore(filterScore=float(filter_score), outFile="validIsoforms.lst",
                                  lqAnnoIso=lqAnnoIsoDataCopy)

        cmd = ('''(grep 'novel' isoFeatures.txt | awk '{if($3>=15 && $4>=0.4){print $1}}'; '''
               '''cut -f 1 validIsoforms.lst) | %s/filter.pl -o - %s -2 4 -m i > hq.collapsed.bed12+'''
               % (q(utilDir), q(inputBed)))
        runner.run(cmd)

        hqBed = os.path.join(hqIsoformDir, "hq.collapsed.bed12+")
        isoform_score = write_isoform_score(featureData, os.path.join(os.path.dirname(hqIsoformDir), "isoIdentity", "isoformScore.txt"))
        print(getCurrentTime() + " Finding high quality novel isoforms done!")
        return {"hq_bed": hqBed, "pu_score": os.path.join(hqIsoformDir, "pu_score.txt"),
                "features": os.path.join(hqIsoformDir, "isoFeatures.txt"),
                "valid": os.path.join(hqIsoformDir, "validIsoforms.lst"), "model": model, "filter_score": float(filter_score),
                "isoform_score": isoform_score}
    finally:
        os.chdir(prevDir)
