"""Smoke test for the sklearn-1.9-adapted BaggingClassifierPU (BUG-107).

Confirms the PU-learning classifier imports, fits, and yields out-of-bag scores in [0, 1].
NOTE: this verifies it RUNS correctly under modern sklearn; numerical equivalence to the
published Python-2 result still requires the real golden baseline.
"""

import numpy as np
from sklearn.tree import DecisionTreeClassifier

from iflas.pu.bagging import BaggingClassifierPU


def test_bagging_pu_fits_and_scores():
    rng = np.random.RandomState(0)
    X = rng.rand(60, 5)
    y = np.array([1] * 30 + [0] * 30)
    cf = BaggingClassifierPU(DecisionTreeClassifier(max_depth=3), n_estimators=20,
                             max_samples=30, n_jobs=1, random_state=0)
    cf.fit(X, y)
    oob = cf.oob_decision_function_
    assert oob.shape == (60, 2)
    scores = oob[:, 1]
    finite = scores[np.isfinite(scores)]
    assert finite.size > 0
    assert finite.min() >= 0.0 and finite.max() <= 1.0


def test_set_oob_score_skips_estimators_with_empty_oob_mask():
    """BUG-132: with few unlabeled points a bootstrap can draw every unlabeled
    sample, leaving that estimator with an empty OOB set. Modern sklearn
    rejects a (0, n_features) predict_proba input, so _set_oob_score must skip
    such estimators. 12 pos + 12 unl, 100 trees, max_samples=12, seed 13
    deterministically produces an empty mask (tree #58 before the fix)."""
    rng = np.random.RandomState(0)
    X = np.vstack([rng.normal(1, 1, (12, 5)), rng.normal(0, 1, (12, 5))])
    y = np.array([1] * 12 + [0] * 12)
    cf = BaggingClassifierPU(DecisionTreeClassifier(), n_estimators=100,
                             max_samples=12, n_jobs=1, random_state=13)
    cf.fit(X, y)                                     # raised ValueError before the fix
    assert any(not mask.any() for mask in (~s for s in cf.estimators_samples_)), \
        "fixture no longer produces an empty OOB mask; pick another seed"
    unl = cf.oob_decision_function_[12:, 1]
    assert np.isfinite(unl).all() and (unl >= 0).all() and (unl <= 1).all()
