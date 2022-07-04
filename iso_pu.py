#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: iso_pu.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2022-04-23 14:22:38
Last modified: 2022-04-23 14:22:38
'''

import subprocess, os, pybedtools, pysam, itertools, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib.patches as patches

# from Bio import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
from collections import Counter
from scipy import interp
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split

from baggingPU import BaggingClassifierPU
from commonFuncs import *
from commonObjs import *

warnings.filterwarnings("ignore")

# import itertools
# import numbers
# # import numpy as np
# from warnings import warn
# from abc import ABCMeta, abstractmethod
#
# from sklearn.base import ClassifierMixin, RegressorMixin
# from sklearn.externals.joblib import Parallel, delayed
# from sklearn.externals.six import with_metaclass
# # from sklearn.externals.six.moves import zip
# from sklearn.metrics import r2_score, accuracy_score
# from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
# from sklearn.utils import check_random_state, check_X_y, check_array, column_or_1d
# from sklearn.utils.random import sample_without_replacement
# from sklearn.utils.validation import has_fit_parameter, check_is_fitted
# from sklearn.utils import indices_to_mask, check_consistent_length
# from sklearn.utils.metaestimators import if_delegate_has_method
# from sklearn.utils.multiclass import check_classification_targets
#
# from sklearn.ensemble.base import BaseEnsemble, _partition_estimators
#
# __all__ = ["BaggingClassifierPU"]
#
# MAX_INT = np.iinfo(np.int32).max


# def _generate_indices(random_state, bootstrap, n_population, n_samples):
#     """Draw randomly sampled indices."""
#     # Draw sample indices
#     if bootstrap:
#         indices = random_state.randint(0, n_population, n_samples)
#     else:
#         indices = sample_without_replacement(n_population, n_samples,
#                                              random_state=random_state)
#
#     return indices
#
#
# def _generate_bagging_indices(random_state, bootstrap_features,
#                               bootstrap_samples, n_features, n_samples,
#                               max_features, max_samples):
#     """Randomly draw feature and sample indices."""
#     # Get valid random state
#     random_state = check_random_state(random_state)
#
#     # Draw indices
#     feature_indices = _generate_indices(random_state, bootstrap_features,
#                                         n_features, max_features)
#     sample_indices = _generate_indices(random_state, bootstrap_samples,
#                                        n_samples, max_samples)
#
#     return feature_indices, sample_indices
#
#
# def _parallel_build_estimators(n_estimators, ensemble, X, y, sample_weight,
#                                seeds, total_n_estimators, verbose):
#     """Private function used to build a batch of estimators within a job."""
#     # Retrieve settings
#     n_samples, n_features = X.shape
#     max_features = ensemble._max_features
#     max_samples = ensemble._max_samples
#     bootstrap = ensemble.bootstrap
#     bootstrap_features = ensemble.bootstrap_features
#     support_sample_weight = has_fit_parameter(ensemble.base_estimator_,
#                                               "sample_weight")
#     if not support_sample_weight and sample_weight is not None:
#         raise ValueError("The base estimator doesn't support sample weight")
#
#     # Build estimators
#     estimators = []
#     estimators_features = []
#
#     for i in range(n_estimators):
#         if verbose > 1:
#             print("Building estimator %d of %d for this parallel run "
#                   "(total %d)..." % (i + 1, n_estimators, total_n_estimators))
#
#         random_state = np.random.RandomState(seeds[i])
#         estimator = ensemble._make_estimator(append=False,
#                                              random_state=random_state)
#
#         ################ MAIN MODIFICATION FOR PU LEARNING ##################
#         iP = [pair[0] for pair in enumerate(y) if pair[1] == 1]
#         iU = [pair[0] for pair in enumerate(y) if pair[1] < 1]
#         features, indices = _generate_bagging_indices(random_state,
#                                                       bootstrap_features,
#                                                       bootstrap, n_features,
#                                                       len(iU), max_features,
#                                                       max_samples)
#         indices = [iU[i] for i in indices] + iP
#         #####################################################################
#
#         # Draw samples, using sample weights, and then fit
#         if support_sample_weight:
#             if sample_weight is None:
#                 curr_sample_weight = np.ones((n_samples,))
#             else:
#                 curr_sample_weight = sample_weight.copy()
#
#             if bootstrap:
#                 sample_counts = np.bincount(indices, minlength=n_samples)
#                 curr_sample_weight *= sample_counts
#             else:
#                 not_indices_mask = ~indices_to_mask(indices, n_samples)
#                 curr_sample_weight[not_indices_mask] = 0
#
#             estimator.fit(X[:, features], y, sample_weight=curr_sample_weight)
#
#         # Draw samples, using a mask, and then fit
#         else:
#             estimator.fit((X[indices])[:, features], y[indices])
#
#         estimators.append(estimator)
#         estimators_features.append(features)
#
#     return estimators, estimators_features
#
#
# def _parallel_predict_proba(estimators, estimators_features, X, n_classes):
#     """Private function used to compute (proba-)predictions within a job."""
#     n_samples = X.shape[0]
#     proba = np.zeros((n_samples, n_classes))
#
#     for estimator, features in zip(estimators, estimators_features):
#         if hasattr(estimator, "predict_proba"):
#             proba_estimator = estimator.predict_proba(X[:, features])
#
#             if n_classes == len(estimator.classes_):
#                 proba += proba_estimator
#
#             else:
#                 proba[:, estimator.classes_] += \
#                     proba_estimator[:, range(len(estimator.classes_))]
#
#         else:
#             # Resort to voting
#             predictions = estimator.predict(X[:, features])
#
#             for i in range(n_samples):
#                 proba[i, predictions[i]] += 1
#
#     return proba
#
#
# def _parallel_predict_log_proba(estimators, estimators_features, X, n_classes):
#     """Private function used to compute log probabilities within a job."""
#     n_samples = X.shape[0]
#     log_proba = np.empty((n_samples, n_classes))
#     log_proba.fill(-np.inf)
#     all_classes = np.arange(n_classes, dtype=np.int)
#
#     for estimator, features in zip(estimators, estimators_features):
#         log_proba_estimator = estimator.predict_log_proba(X[:, features])
#
#         if n_classes == len(estimator.classes_):
#             log_proba = np.logaddexp(log_proba, log_proba_estimator)
#
#         else:
#             log_proba[:, estimator.classes_] = np.logaddexp(
#                 log_proba[:, estimator.classes_],
#                 log_proba_estimator[:, range(len(estimator.classes_))])
#
#             missing = np.setdiff1d(all_classes, estimator.classes_)
#             log_proba[:, missing] = np.logaddexp(log_proba[:, missing],
#                                                  -np.inf)
#
#     return log_proba
#
#
# def _parallel_decision_function(estimators, estimators_features, X):
#     """Private function used to compute decisions within a job."""
#     return sum(estimator.decision_function(X[:, features])
#                for estimator, features in zip(estimators,
#                                               estimators_features))


# class BaseBaggingPU(with_metaclass(ABCMeta, BaseEnsemble)):
#     """Base class for Bagging PU meta-estimator.
#     Warning: This class should not be used directly. Use derived classes
#     instead.
#     """
#
#     @abstractmethod
#     def __init__(self,
#                  base_estimator=None,
#                  n_estimators=10,
#                  max_samples=1.0,
#                  max_features=1.0,
#                  bootstrap=True,
#                  bootstrap_features=False,
#                  oob_score=True,
#                  warm_start=False,
#                  n_jobs=1,
#                  random_state=None,
#                  verbose=0):
#         super(BaseBaggingPU, self).__init__(
#             base_estimator=base_estimator,
#             n_estimators=n_estimators)
#
#         self.max_samples = max_samples
#         self.max_features = max_features
#         self.bootstrap = bootstrap
#         self.bootstrap_features = bootstrap_features
#         self.oob_score = oob_score
#         self.warm_start = warm_start
#         self.n_jobs = n_jobs
#         self.random_state = random_state
#         self.verbose = verbose
#
#     def fit(self, X, y, sample_weight=None):
#         """Build a Bagging ensemble of estimators from the training
#            set (X, y).
#         Parameters
#         ----------
#         X : {array-like, sparse matrix} of shape = [n_samples, n_features]
#             The training input samples. Sparse matrices are accepted only if
#             they are supported by the base estimator.
#         y : array-like, shape = [n_samples]
#             The target values (1 for positive, 0 for unlabeled).
#         sample_weight : array-like, shape = [n_samples] or None
#             Sample weights. If None, then samples are equally weighted.
#             Note that this is supported only if the base estimator supports
#             sample weighting.
#         Returns
#         -------
#         self : object
#             Returns self.
#         """
#         return self._fit(X, y, self.max_samples, sample_weight=sample_weight)
#
#     def _fit(self, X, y, max_samples=None, max_depth=None, sample_weight=None):
#         """Build a Bagging ensemble of estimators from the training
#            set (X, y).
#         Parameters
#         ----------
#         X : {array-like, sparse matrix} of shape = [n_samples, n_features]
#             The training input samples. Sparse matrices are accepted only if
#             they are supported by the base estimator.
#         y : array-like, shape = [n_samples]
#             The target values (1 for positive, 0 for unlabeled).
#         max_samples : int or float, optional (default=None)
#             Argument to use instead of self.max_samples.
#         max_depth : int, optional (default=None)
#             Override value used when constructing base estimator. Only
#             supported if the base estimator has a max_depth parameter.
#         sample_weight : array-like, shape = [n_samples] or None
#             Sample weights. If None, then samples are equally weighted.
#             Note that this is supported only if the base estimator supports
#             sample weighting.
#         Returns
#         -------
#         self : object
#             Returns self.
#         """
#         random_state = check_random_state(self.random_state)
#
#         self.y = y
#
#         # Convert data
#         X, y = check_X_y(X, y, ['csr', 'csc'])
#         if sample_weight is not None:
#             sample_weight = check_array(sample_weight, ensure_2d=False)
#             check_consistent_length(y, sample_weight)
#
#         # Remap output
#         n_samples, self.n_features_ = X.shape
#         self._n_samples = n_samples
#         y = self._validate_y(y)
#
#         # Check parameters
#         self._validate_estimator()
#
#         if max_depth is not None:
#             self.base_estimator_.max_depth = max_depth
#
#         # Validate max_samples
#         if max_samples is None:
#             max_samples = self.max_samples
#         elif not isinstance(max_samples, (numbers.Integral, np.integer)):
#             max_samples = int(max_samples * sum(y < 1))
#
#         if not (0 < max_samples <= sum(y < 1)):
#             raise ValueError("max_samples must be positive"
#                              " and no larger than the number of unlabeled points")
#
#         # Store validated integer row sampling value
#         self._max_samples = max_samples
#
#         # Validate max_features
#         if isinstance(self.max_features, (numbers.Integral, np.integer)):
#             max_features = self.max_features
#         else:  # float
#             max_features = int(self.max_features * self.n_features_)
#
#         if not (0 < max_features <= self.n_features_):
#             raise ValueError("max_features must be in (0, n_features]")
#
#         # Store validated integer feature sampling value
#         self._max_features = max_features
#
#         # Other checks
#         if not self.bootstrap and self.oob_score:
#             raise ValueError("Out of bag estimation only available"
#                              " if bootstrap=True")
#
#         if self.warm_start and self.oob_score:
#             raise ValueError("Out of bag estimate only available"
#                              " if warm_start=False")
#
#         if hasattr(self, "oob_score_") and self.warm_start:
#             del self.oob_score_
#
#         if not self.warm_start or not hasattr(self, 'estimators_'):
#             # Free allocated memory, if any
#             self.estimators_ = []
#             self.estimators_features_ = []
#
#         n_more_estimators = self.n_estimators - len(self.estimators_)
#
#         if n_more_estimators < 0:
#             raise ValueError('n_estimators=%d must be larger or equal to '
#                              'len(estimators_)=%d when warm_start==True'
#                              % (self.n_estimators, len(self.estimators_)))
#
#         elif n_more_estimators == 0:
#             warn("Warm-start fitting without increasing n_estimators does not "
#                  "fit new trees.")
#             return self
#
#         # Parallel loop
#         n_jobs, n_estimators, starts = _partition_estimators(n_more_estimators,
#                                                              self.n_jobs)
#         total_n_estimators = sum(n_estimators)
#
#         # Advance random state to state after training
#         # the first n_estimators
#         if self.warm_start and len(self.estimators_) > 0:
#             random_state.randint(MAX_INT, size=len(self.estimators_))
#
#         seeds = random_state.randint(MAX_INT, size=n_more_estimators)
#         self._seeds = seeds
#
#         all_results = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
#             delayed(_parallel_build_estimators)(
#                 n_estimators[i],
#                 self,
#                 X,
#                 y,
#                 sample_weight,
#                 seeds[starts[i]:starts[i + 1]],
#                 total_n_estimators,
#                 verbose=self.verbose)
#             for i in range(n_jobs))
#
#         # Reduce
#         self.estimators_ += list(itertools.chain.from_iterable(
#             t[0] for t in all_results))
#         self.estimators_features_ += list(itertools.chain.from_iterable(
#             t[1] for t in all_results))
#
#         if self.oob_score:
#             self._set_oob_score(X, y)
#
#         return self
#
#     @abstractmethod
#     def _set_oob_score(self, X, y):
#         """Calculate out of bag predictions and score."""
#
#     def _validate_y(self, y):
#         # Default implementation
#         return column_or_1d(y, warn=True)
#
#     def _get_estimators_indices(self):
#         # Get drawn indices along both sample and feature axes
#         for seed in self._seeds:
#             # Operations accessing random_state must be performed identically
#             # to those in `_parallel_build_estimators()`
#             random_state = np.random.RandomState(seed)
#
#             ############ MAIN MODIFICATION FOR PU LEARNING ###############
#             iP = [pair[0] for pair in enumerate(self.y) if pair[1] == 1]
#             iU = [pair[0] for pair in enumerate(self.y) if pair[1] < 1]
#
#             feature_indices, sample_indices = _generate_bagging_indices(
#                 random_state, self.bootstrap_features, self.bootstrap,
#                 self.n_features_, len(iU), self._max_features,
#                 self._max_samples)
#
#             sample_indices = [iU[i] for i in sample_indices] + iP
#             ###############################################################
#
#             yield feature_indices, sample_indices
#
#     @property
#     def estimators_samples_(self):
#         """The subset of drawn samples for each base estimator.
#         Returns a dynamically generated list of boolean masks identifying
#         the samples used for fitting each member of the ensemble, i.e.,
#         the in-bag samples.
#         Note: the list is re-created at each call to the property in order
#         to reduce the object memory footprint by not storing the sampling
#         data. Thus fetching the property may be slower than expected.
#         """
#         sample_masks = []
#         for _, sample_indices in self._get_estimators_indices():
#             mask = indices_to_mask(sample_indices, self._n_samples)
#             sample_masks.append(mask)
#
#         return sample_masks
#
#
# class BaggingClassifierPU(BaseBaggingPU, ClassifierMixin):
#     """A Bagging PU classifier.
#     Adapted from sklearn.ensemble.BaggingClassifier, based on
#     A bagging SVM to learn from positive and unlabeled examples (2013) by Mordelet and Vert
#     http://dx.doi.org/10.1016/j.patrec.2013.06.010
#     http://members.cbio.mines-paristech.fr/~jvert/svn/bibli/local/Mordelet2013bagging.pdf
#
#     Parameters
#     ----------
#     base_estimator : object or None, optional (default=None)
#         The base estimator to fit on random subsets of the dataset.
#         If None, then the base estimator is a decision tree.
#     n_estimators : int, optional (default=10)
#         The number of base estimators in the ensemble.
#     max_samples : int or float, optional (default=1.0)
#         The number of unlabeled samples to draw to train each base estimator.
#     max_features : int or float, optional (default=1.0)
#         The number of features to draw from X to train each base estimator.
#         - If int, then draw `max_features` features.
#         - If float, then draw `max_features * X.shape[1]` features.
#     bootstrap : boolean, optional (default=True)
#         Whether samples are drawn with replacement.
#     bootstrap_features : boolean, optional (default=False)
#         Whether features are drawn with replacement.
#     oob_score : bool, optional (default=True)
#         Whether to use out-of-bag samples to estimate
#         the generalization error.
#     warm_start : bool, optional (default=False)
#         When set to True, reuse the solution of the previous call to fit
#         and add more estimators to the ensemble, otherwise, just fit
#         a whole new ensemble.
#     n_jobs : int, optional (default=1)
#         The number of jobs to run in parallel for both `fit` and `predict`.
#         If -1, then the number of jobs is set to the number of cores.
#     random_state : int, RandomState instance or None, optional (default=None)
#         If int, random_state is the seed used by the random number generator;
#         If RandomState instance, random_state is the random number generator;
#         If None, the random number generator is the RandomState instance used
#         by `np.random`.
#     verbose : int, optional (default=0)
#         Controls the verbosity of the building process.
#     Attributes
#     ----------
#     base_estimator_ : estimator
#         The base estimator from which the ensemble is grown.
#     estimators_ : list of estimators
#         The collection of fitted base estimators.
#     estimators_samples_ : list of arrays
#         The subset of drawn samples (i.e., the in-bag samples) for each base
#         estimator. Each subset is defined by a boolean mask.
#     estimators_features_ : list of arrays
#         The subset of drawn features for each base estimator.
#     classes_ : array of shape = [n_classes]
#         The classes labels.
#     n_classes_ : int or list
#         The number of classes.
#     oob_score_ : float
#         Score of the training dataset obtained using an out-of-bag estimate.
#     oob_decision_function_ : array of shape = [n_samples, n_classes]
#         Decision function computed with out-of-bag estimate on the training
#         set. Positive data points, and perhaps some of the unlabeled,
#         are left out during the bootstrap. In these cases,
#         `oob_decision_function_` contains NaN.
#     """
#
#     def __init__(self,
#                  base_estimator=None,
#                  n_estimators=10,
#                  max_samples=1.0,
#                  max_features=1.0,
#                  bootstrap=True,
#                  bootstrap_features=False,
#                  oob_score=True,
#                  warm_start=False,
#                  n_jobs=1,
#                  random_state=None,
#                  verbose=0):
#
#         super(BaggingClassifierPU, self).__init__(
#             base_estimator,
#             n_estimators=n_estimators,
#             max_samples=max_samples,
#             max_features=max_features,
#             bootstrap=bootstrap,
#             bootstrap_features=bootstrap_features,
#             oob_score=oob_score,
#             warm_start=warm_start,
#             n_jobs=n_jobs,
#             random_state=random_state,
#             verbose=verbose)
#
#     def _validate_estimator(self):
#         """Check the estimator and set the base_estimator_ attribute."""
#         super(BaggingClassifierPU, self)._validate_estimator(
#             default=DecisionTreeClassifier())
#
#     def _set_oob_score(self, X, y):
#         n_samples = y.shape[0]
#         n_classes_ = self.n_classes_
#         classes_ = self.classes_
#
#         predictions = np.zeros((n_samples, n_classes_))
#
#         for estimator, samples, features in zip(self.estimators_,
#                                                 self.estimators_samples_,
#                                                 self.estimators_features_):
#             # Create mask for OOB samples
#             mask = ~samples
#
#             if hasattr(estimator, "predict_proba"):
#                 predictions[mask, :] += estimator.predict_proba(
#                     (X[mask, :])[:, features])
#
#             else:
#                 p = estimator.predict((X[mask, :])[:, features])
#                 j = 0
#
#                 for i in range(n_samples):
#                     if mask[i]:
#                         predictions[i, p[j]] += 1
#                         j += 1
#
#         # Modified: no warnings about non-OOB points (i.e. positives)
#         with np.errstate(invalid='ignore'):
#             oob_decision_function = (predictions /
#                                      predictions.sum(axis=1)[:, np.newaxis])
#             oob_score = accuracy_score(y, np.argmax(predictions, axis=1))
#
#         self.oob_decision_function_ = oob_decision_function
#         self.oob_score_ = oob_score
#
#     def _validate_y(self, y):
#         y = column_or_1d(y, warn=True)
#         check_classification_targets(y)
#         self.classes_, y = np.unique(y, return_inverse=True)
#         self.n_classes_ = len(self.classes_)
#
#         return y
#
#     def predict(self, X):
#         """Predict class for X.
#         The predicted class of an input sample is computed as the class with
#         the highest mean predicted probability. If base estimators do not
#         implement a ``predict_proba`` method, then it resorts to voting.
#         Parameters
#         ----------
#         X : {array-like, sparse matrix} of shape = [n_samples, n_features]
#             The training input samples. Sparse matrices are accepted only if
#             they are supported by the base estimator.
#         Returns
#         -------
#         y : array of shape = [n_samples]
#             The predicted classes.
#         """
#         predicted_probabilitiy = self.predict_proba(X)
#         return self.classes_.take((np.argmax(predicted_probabilitiy, axis=1)),
#                                   axis=0)
#
#     def predict_proba(self, X):
#         """Predict class probabilities for X.
#         The predicted class probabilities of an input sample is computed as
#         the mean predicted class probabilities of the base estimators in the
#         ensemble. If base estimators do not implement a ``predict_proba``
#         method, then it resorts to voting and the predicted class probabilities
#         of an input sample represents the proportion of estimators predicting
#         each class.
#         Parameters
#         ----------
#         X : {array-like, sparse matrix} of shape = [n_samples, n_features]
#             The training input samples. Sparse matrices are accepted only if
#             they are supported by the base estimator.
#         Returns
#         -------
#         p : array of shape = [n_samples, n_classes]
#             The class probabilities of the input samples. The order of the
#             classes corresponds to that in the attribute `classes_`.
#         """
#         check_is_fitted(self, "classes_")
#         # Check data
#         X = check_array(X, accept_sparse=['csr', 'csc'])
#
#         if self.n_features_ != X.shape[1]:
#             raise ValueError("Number of features of the model must "
#                              "match the input. Model n_features is {0} and "
#                              "input n_features is {1}."
#                              "".format(self.n_features_, X.shape[1]))
#
#         # Parallel loop
#         n_jobs, n_estimators, starts = _partition_estimators(self.n_estimators,
#                                                              self.n_jobs)
#
#         all_proba = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
#             delayed(_parallel_predict_proba)(
#                 self.estimators_[starts[i]:starts[i + 1]],
#                 self.estimators_features_[starts[i]:starts[i + 1]],
#                 X,
#                 self.n_classes_)
#             for i in range(n_jobs))
#
#         # Reduce
#         proba = sum(all_proba) / self.n_estimators
#
#         return proba
#
#     def predict_log_proba(self, X):
#         """Predict class log-probabilities for X.
#         The predicted class log-probabilities of an input sample is computed as
#         the log of the mean predicted class probabilities of the base
#         estimators in the ensemble.
#         Parameters
#         ----------
#         X : {array-like, sparse matrix} of shape = [n_samples, n_features]
#             The training input samples. Sparse matrices are accepted only if
#             they are supported by the base estimator.
#         Returns
#         -------
#         p : array of shape = [n_samples, n_classes]
#             The class log-probabilities of the input samples. The order of the
#             classes corresponds to that in the attribute `classes_`.
#         """
#         check_is_fitted(self, "classes_")
#         if hasattr(self.base_estimator_, "predict_log_proba"):
#             # Check data
#             X = check_array(X, accept_sparse=['csr', 'csc'])
#
#             if self.n_features_ != X.shape[1]:
#                 raise ValueError("Number of features of the model must "
#                                  "match the input. Model n_features is {0} "
#                                  "and input n_features is {1} "
#                                  "".format(self.n_features_, X.shape[1]))
#
#             # Parallel loop
#             n_jobs, n_estimators, starts = _partition_estimators(
#                 self.n_estimators, self.n_jobs)
#
#             all_log_proba = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
#                 delayed(_parallel_predict_log_proba)(
#                     self.estimators_[starts[i]:starts[i + 1]],
#                     self.estimators_features_[starts[i]:starts[i + 1]],
#                     X,
#                     self.n_classes_)
#                 for i in range(n_jobs))
#
#             # Reduce
#             log_proba = all_log_proba[0]
#
#             for j in range(1, len(all_log_proba)):
#                 log_proba = np.logaddexp(log_proba, all_log_proba[j])
#
#             log_proba -= np.log(self.n_estimators)
#
#             return log_proba
#
#         else:
#             return np.log(self.predict_proba(X))
#
#     @if_delegate_has_method(delegate='base_estimator')
#     def decision_function(self, X):
#         """Average of the decision functions of the base classifiers.
#         Parameters
#         ----------
#         X : {array-like, sparse matrix} of shape = [n_samples, n_features]
#             The training input samples. Sparse matrices are accepted only if
#             they are supported by the base estimator.
#         Returns
#         -------
#         score : array, shape = [n_samples, k]
#             The decision function of the input samples. The columns correspond
#             to the classes in sorted order, as they appear in the attribute
#             ``classes_``. Regression and binary classification are special
#             cases with ``k == 1``, otherwise ``k==n_classes``.
#         """
#         check_is_fitted(self, "classes_")
#
#         # Check data
#         X = check_array(X, accept_sparse=['csr', 'csc'])
#
#         if self.n_features_ != X.shape[1]:
#             raise ValueError("Number of features of the model must "
#                              "match the input. Model n_features is {0} and "
#                              "input n_features is {1} "
#                              "".format(self.n_features_, X.shape[1]))
#
#         # Parallel loop
#         n_jobs, n_estimators, starts = _partition_estimators(self.n_estimators,
#                                                              self.n_jobs)
#
#         all_decisions = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
#             delayed(_parallel_decision_function)(
#                 self.estimators_[starts[i]:starts[i + 1]],
#                 self.estimators_features_[starts[i]:starts[i + 1]],
#                 X)
#             for i in range(n_jobs))
#
#         # Reduce
#         decisions = sum(all_decisions) / self.n_estimators
#
#         return decisions


# class Bed6(object):
#     "BED6 format gene structure"
#
#     class BedError(Exception):
#         "Error in manipulating bed12 structures"
#
#         def __init__(self, value):
#             self.value = value
#
#         def __str__(self):
#             return repr(self.value)
#
#     def __init__(self, line=None):
#         if line:
#             self.record = line.strip().split("\t")
#             (self.chrom, self.chromStart, self.chromEnd, self.name,
#              self.score, self.strand) = self.record[:6]
#             self.chromStart = int(self.chromStart)
#             self.chromEnd = int(self.chromEnd)
#             try:
#                 self.score = int(float(self.score))
#             except ValueError:
#                 pass
#             self.pos = (str(self.chrom) + ":" + str(self.chromStart) +
#                         "-" + str(self.chromEnd))
#         else:
#             self.empty()
#
#     def empty(self):
#         (self.chrom, self.chromStart, self.chromEnd, self.name,
#          self.score, self.strand) = ("", 0, 0, "", 0, "")
#
#     def __repr__(self):
#         "return a line of bed6 format, without newline ending"
#         fields = [self.chrom, str(self.chromStart), str(self.chromEnd),
#                   self.name, str(self.score), self.strand]
#         return "\t".join(fields)
#
#
# class Bed12(object):
#     "BED12 format gene structure."
#
#     class BedError(Exception):
#         "Error in manipulating Bed12 structures"
#
#         def __init__(self, value):
#             self.value = value
#
#         def __str__(self):
#             return repr(self.value)
#
#     def __init__(self, line=None):
#         if line:
#             self.record = line.strip().split("\t")
#             (self.chrom, self.chromStart, self.chromEnd, self.name,
#              self.score, self.strand, self.thickStart, self.thickEnd,
#              self.itemRgb, self.blockCount, self.blockSizes,
#              self.blockStarts) = self.record[:12]
#
#             self.chromStart = int(self.chromStart)
#             self.chromEnd = int(self.chromEnd)
#             self.score = int(float(self.score))
#             self.thickStart = int(self.thickStart)
#             self.thickEnd = int(self.thickEnd)
#             self.blockCount = int(self.blockCount)
#
#             self.blockSizes = self.blockSizes.strip(",").split(",")
#             self.blockStarts = self.blockStarts.strip(",").split(",")
#
#             assert len(self.blockStarts) == len(self.blockSizes)
#             assert len(self.blockStarts) == self.blockCount
#
#             for i in range(self.blockCount):
#                 self.blockSizes[i] = int(self.blockSizes[i])
#                 self.blockStarts[i] = int(self.blockStarts[i])
#
#             self.exonStarts = []
#             self.exonEnds = []
#             for i in range(self.blockCount):
#                 self.exonStarts.append(self.chromStart + self.blockStarts[i])
#                 self.exonEnds.append(self.exonStarts[i] + self.blockSizes[i])
#
#             self.exons = self.parse_exon()
#             self.introns = self.parse_intron()
#             self.exonChain = ';'.join(map(lambda x: str(self.exons[x][0])+"-"+str(self.exons[x][1]), range(len(self.exons))))
#             self.juncChain = ";".join(map(lambda x: str(self.introns[x][0])+"-"+str(self.introns[x][1]), range(len(self.introns))))
#         else:
#             self.empty()
#
#     def empty(self):
#         "return an empty class with all values None, '', or []."
#         self.chrom = ""
#         self.chromStart = self.chromEnd = 0
#         self.name = ""
#         self.score = 0
#         self.strand = ""
#         self.thickStart = self.thickEnd = 0
#         self.itemRgb = "0,0,0"
#         self.blockCount = 0
#         self.blockSizes = []
#         self.blockStarts = []
#
#     def __len__(self):
#         "the length of transcript"
#         return sum(self.blockSizes)
#
#     def parse_exon(self):
#         "return a list of exon pos [(st, ed), (st, ed) , ... ]"
#         exons = []
#         for i in range(self.blockCount):
#             st = self.exonStarts[i]
#             ed = self.exonEnds[i]
#             exons.append((st, ed))
#         return exons
#
#     def parse_intron(self):
#         "return a list of intron pos [(st, ed], (st, ed], ... ]"
#         introns = []
#         for i in range(self.blockCount - 1):
#             st = self.exonEnds[i]
#             ed = self.exonStarts[i + 1]
#             introns.append((st, ed))
#         return introns
#
#     def __repr__(self):
#         "return a line of bed12 format, without newline ending"
#         fields = [self.chrom, str(self.chromStart), str(self.chromEnd),
#                   self.name, str(self.score), self.strand, str(self.thickStart),
#                   str(self.thickEnd), self.itemRgb, str(self.blockCount)]
#
#         blockSizesline = ','.join(repr(i) for i in self.blockSizes)
#         blockStartsline = ','.join(repr(i) for i in self.blockStarts)
#         fields.extend((blockSizesline, blockStartsline))
#         return "\t".join(fields)
#
#
# class Bed6Plus(Bed6):
#     def __init__(self, line=None):
#         if line:
#             Bed6.__init__(self, line)
#             self.otherList = self.record[6:]
#
#     def __str__(self):
#         return Bed6.__str__(self) + "\t" + "\t".join(self.otherList)
#
#
# class Bed12Plus(Bed12):
#     def __init__(self, line=None):
#         if line:
#             Bed12.__init__(self, line)
#             self.otherList = self.record[12:]
#
#     def __str__(self):
#         return Bed12.__str__(self) + "\t" + "\t".join(self.otherList)
#
#
# class BedFile(object):
#     def __init__(self, bedFile, type=None):
#         self.bedFile = bedFile
#         self.reads = self.getReadsInfo(type)
#
#     def getReadsInfo(self, bedType):
#         readsDict = {}
#         with open(self.bedFile) as f:
#             for line in f:
#                 if bedType == "bed12":
#                     b = Bed12(line)
#                     readsDict.__setitem__(b.name, b)
#                 elif bedType == "bed12+":
#                     b = Bed12Plus(line)
#                     readsDict.__setitem__(b.name, b)
#                 elif bedType == "bed6":
#                     b = Bed6(line)
#                     readsDict.__setitem__(b.name, b)
#                 elif bedType == "bed6+":
#                     b = Bed6Plus(line)
#                     readsDict.__setitem__(b.name, b)
#         return readsDict
#
#     def getJuncChainDict(self):
#         juncChainDict = {}
#         annoSingleExonList = []
#         for i in self.reads:
#             if len(self.reads[i].exons) > 1:
#                 juncChainInfo = "{}:{}".format(self.reads[i].chrom, self.reads[i].juncChain)
#                 if juncChainInfo not in juncChainDict:
#                     juncChainDict[juncChainInfo] = self.reads[i]
#             # else:
#             #     singleExon = "\t".join([self.reads[i].chrom, self.reads[i].chromStart, self.reads[i].chromEnd, self.reads[i].name, ".", self.reads[i].strand])
#             #     annoSingleExonList.append(singleExon)
#         # return juncDict, annoSingleExonList
#         return juncChainDict
#
#     def getAllJuncDict(self):
#         allJuncDict = {}
#         for i in self.reads:
#             if len(self.reads[i].exons) > 1:
#                 for junc in self.reads[i].introns:
#                     juncName = "{}:{}-{}".format(self.reads[i].chrom, junc[0]+1, junc[1])
#                     if juncName not in allJuncDict:
#                         allJuncDict[juncName] = "\t".join(map(str, [self.reads[i].chrom, junc[0], junc[1], juncName, self.reads[i].score, self.reads[i].strand]))
#         return allJuncDict
#
#     def getAllExonDict(self):
#         allExonDict = {}
#         for i in self.reads:
#             for exon in self.reads[i].exons:
#                 exonName = "{}:{}-{}".format(self.reads[i].chrom, exon[0], exon[1])
#                 if exonName not in allExonDict:
#                     allExonDict[exonName] = "\t".join(map(str, [self.reads[i].chrom, exon[0], exon[1], exonName, self.reads[i].score, self.reads[i].strand]))
#         return allExonDict


# class PUAdapter(object):
#     """
#     Adapts any probabilistic binary classifier to positive-unlabled learning using the PosOnly method proposed by
#     Elkan and Noto:
#
#     Elkan, Charles, and Keith Noto. \"Learning classifiers from only positive and unlabeled data.\"
#     Proceeding of the 14th ACM SIGKDD international conference on Knowledge discovery and data mining. ACM, 2008.
#     """
#
#     def __init__(self, estimator, hold_out_ratio=0.1, precomputed_kernel=False):
#         """
#         estimator -- An estimator of p(s=1|x) that must implement:
#                      * predict_proba(X): Takes X, which can be a list of feature vectors or a precomputed
#                                          kernel matrix and outputs p(s=1|x) for each example in X
#                      * fit(X,y): Takes X, which can be a list of feature vectors or a precomputed
#                                  kernel matrix and takes y, which are the labels associated to the
#                                  examples in X
#         hold_out_ratio -- The ratio of training examples that must be held out of the training set of examples
#                           to estimate p(s=1|y=1) after training the estimator
#         precomputed_kernel -- Specifies if the X matrix for predict_proba and fit is a precomputed kernel matrix
#         """
#         self.estimator = estimator
#         self.c = 1.0
#         self.hold_out_ratio = hold_out_ratio
#
#         if precomputed_kernel:
#             self.fit = self.__fit_precomputed_kernel
#         else:
#             self.fit = self.__fit_no_precomputed_kernel
#
#         self.estimator_fitted = False
#
#     def __str__(self):
#         return 'Estimator:' + str(self.estimator) + '\n' + 'p(s=1|y=1,x) ~= ' + str(self.c) + '\n' + \
#                'Fitted: ' + str(self.estimator_fitted)
#
#     def __fit_precomputed_kernel(self, X, y):
#         """
#         Fits an estimator of p(s=1|x) and estimates the value of p(s=1|y=1) using a subset of the training examples
#
#         X -- Precomputed kernel matrix
#         y -- Labels associated to each example in X (Positive label: 1.0, Negative label: -1.0)
#         """
#         positives = np.where(y == 1.)[0]
#         hold_out_size = int(np.ceil(len(positives) * self.hold_out_ratio))
#
#         if len(positives) <= hold_out_size:
#             raise ('Not enough positive examples to estimate p(s=1|y=1,x). Need at least ' + str(
#                 hold_out_size + 1) + '.')
#
#         np.random.shuffle(positives)
#         hold_out = positives[:hold_out_size]
#
#         # Hold out test kernel matrix
#         X_test_hold_out = X[hold_out]
#         keep = list(set(np.arange(len(y))) - set(hold_out))
#         X_test_hold_out = X_test_hold_out[:, keep]
#
#         # New training kernel matrix
#         X = X[:, keep]
#         X = X[keep]
#
#         y = np.delete(y, hold_out)
#
#         self.estimator.fit(X, y)
#
#         hold_out_predictions = self.estimator.predict_proba(X_test_hold_out)
#
#         try:
#             hold_out_predictions = hold_out_predictions[:, 1]
#         except:
#             pass
#
#         c = np.mean(hold_out_predictions)
#         self.c = c
#
#         self.estimator_fitted = True
#
#     def __fit_no_precomputed_kernel(self, X, y):
#         """
#         Fits an estimator of p(s=1|x) and estimates the value of p(s=1|y=1,x)
#
#         X -- List of feature vectors
#         y -- Labels associated to each feature vector in X (Positive label: 1.0, Negative label: -1.0)
#         """
#         positives = np.where(y == 1.)[0]
#         hold_out_size = int(np.ceil(len(positives) * self.hold_out_ratio))
#
#         if len(positives) <= hold_out_size:
#             raise ('Not enough positive examples to estimate p(s=1|y=1,x). Need at least ' + str(
#                 hold_out_size + 1) + '.')
#
#         np.random.shuffle(positives)
#         hold_out = positives[:hold_out_size]
#         X_hold_out = X[hold_out]
#         X = np.delete(X, hold_out, 0)
#         y = np.delete(y, hold_out)
#
#         self.estimator.fit(X, y)
#
#         hold_out_predictions = self.estimator.predict_proba(X_hold_out)
#
#         try:
#             hold_out_predictions = hold_out_predictions[:, 1]
#         except:
#             pass
#
#         c = np.mean(hold_out_predictions)
#         self.c = c
#
#         self.estimator_fitted = True
#
#     def predict_proba(self, X):
#         """
#         Predicts p(y=1|x) using the estimator and the value of p(s=1|y=1) estimated in fit(...)
#
#         X -- List of feature vectors or a precomputed kernel matrix
#         """
#         if not self.estimator_fitted:
#             raise Exception('The estimator must be fitted before calling predict_proba(...).')
#
#         probabilistic_predictions = self.estimator.predict_proba(X)
#
#         try:
#             probabilistic_predictions = probabilistic_predictions[:, 1]
#         except:
#             pass
#
#         return probabilistic_predictions / self.c
#
#     def predict(self, X, treshold=0.5):
#         """
#         Assign labels to feature vectors based on the estimator's predictions
#
#         X -- List of feature vectors or a precomputed kernel matrix
#         treshold -- The decision treshold between the positive and the negative class
#         """
#         if not self.estimator_fitted:
#             raise Exception('The estimator must be fitted before calling predict(...).')
#
#         return np.array([1. if p > treshold else -1. for p in self.predict_proba(X)])


class SchemeModelEval(object):
    def __init__(self, trainingData=None, model=None, scheme=None, featureData=None):
        self.trainingData = trainingData
        self.model = model
        self.scheme = scheme

        self.finalEstimator = None
        self.predResults = None
        self.f1_score = None
        self.precision = None
        self.recall = None

        self.featureData = featureData

    def eval(self):
        X_train = self.trainingData.iloc[:, :-1]
        y_train = self.trainingData.iloc[:, -1]

        classifier = self.getClassifier()
        if self.scheme == "bagging":
            X = X_train
            y = y_train
            y = pd.Series(y)

            labels = Counter(y)
            sorted_labels = sorted(labels.items(), key=lambda x: x[1], reverse=True)
            max_samples = sorted_labels[0][1] if sorted_labels[0][0] == 0 else sorted_labels[1][1]

            # from baggingPU import BaggingClassifierPU
            cf = BaggingClassifierPU(
                classifier,
                n_estimators=100,
                max_samples=max_samples,
                n_jobs=-1
            )
            cf.fit(X, y)

            pu_score = cf.oob_decision_function_[:, 1]
            self.predResults = pd.DataFrame({
                "label": y,
                "pu_score": pu_score
            }, columns=["label", "pu_score"])
            self.finalEstimator = cf
        elif self.scheme == "2step":
            X = X_train
            y = y_train
            y = pd.Series(y)

            cf = classifier
            cf.fit(X, y)

            ys = 2*y -1
            pred = cf.predict_proba(X)[:, 1]

            # Find the range of scores given to positive data points
            range_P = [min(pred * (ys > 0)), max(pred * (ys > 0))]

            # STEP 1
            # If any unlabeled point has a score above all known positives,
            # or below all known positives, label it accordingly
            iP_new = ys[(ys < 0) & (pred >= range_P[1])].index
            iN_new = ys[(ys < 0) & (pred <= range_P[0])].index
            ys.loc[iP_new] = 1
            ys.loc[iN_new] = 0

            # Classifier to be used for step 2
            cf2 = classifier

            # Limit to 10 iterations (this is arbitrary, but
            # otherwise this approach can take a very long time)
            for i in range(10):
                # If step 1 didn't find new labels, we're done
                if len(iP_new) + len(iN_new) == 0 and i > 0:
                    break

                # STEP 2
                # Retrain on new labels and get new scores
                cf2.fit(X, ys)
                pred = cf2.predict_proba(X)[:, -1]

                # Find the range of scores given to positive data points
                range_P = [min(pred * (ys > 0)), max(pred * (ys > 0))]

                # Repeat step 1
                iP_new = ys[(ys < 0) & (pred >= range_P[1])].index
                iN_new = ys[(ys < 0) & (pred <= range_P[0])].index
                ys.loc[iP_new] = 1
                ys.loc[iN_new] = 0

            # Lastly, get the scores assigned by this approach
            self.predResults = pd.DataFrame({
                "label": y,
                "pu_score": pred
            }, columns=["label", "pu_score"])
            self.finalEstimator = cf2
        else:
            X = X_train
            y = y_train
            y = pd.Series(y)

            cf = classifier
            cf.fit(X, y)

            self.predResults = pd.DataFrame({
                "label": y,
                "pu_score": cf.predict_proba(X)[:, 1]
            }, columns=["label", "pu_score"])
            self.finalEstimator = cf

        # self.cfTest(self.finalEstimator, X_test, y_test)
        # self.filterIsoformsByScore(filterScore=0.8, outFile="validIsoforms.lst")

    def getClassifier(self):
        if self.model == "SVM":
            from sklearn.svm import SVC
            return SVC(kernel='rbf', gamma='auto', random_state=0)
        elif self.model == "DT":
            from sklearn.tree import DecisionTreeClassifier
            return DecisionTreeClassifier()
        elif self.model == "RF":
            from sklearn.ensemble import RandomForestClassifier
            return RandomForestClassifier(n_estimators=10, random_state=0)
        elif self.model == "GB":
            from sklearn.ensemble import GradientBoostingClassifier
            return GradientBoostingClassifier(n_estimators=10, random_state=0)
        elif self.model == "NB":
            from sklearn.naive_bayes import GaussianNB
            return GaussianNB()

    def estimatorTest(self, testData=None):
        X_test = testData.iloc[:, :-1]
        y_test = testData.iloc[:, -1]

        from sklearn.metrics import precision_recall_fscore_support
        y_pred = self.finalEstimator.predict(X_test)
        y_pred_1 = y_pred.copy()
        y_pred_1[np.where(y_pred_1 == 0)[0]] = -1
        precision, recall, f1_score, _ = precision_recall_fscore_support(y_test, y_pred_1)

        self.f1_score = f1_score[1]
        self.precision = precision[1]
        self.recall = recall[1]
        print("F1 score:", f1_score[1])
        print("Precision:", precision[1])
        print("Recall:", recall[1])


    def filterIsoformsByScore(self, filterScore=0.0, outFile="validIsoform.lst"):
        hqNovelIsoRes = self.predResults.loc[(self.predResults.label==0) & (self.predResults.pu_score>=filterScore), ]
        annoIsoRes = self.predResults.loc[self.predResults.label==1,]
        outResults = pd.concat([annoIsoRes, hqNovelIsoRes])
        outResults.to_csv(outFile, sep="\t", index=True, index_label="isoform")


class GetFeatures(object):
    def __init__(self, inputBed=None, genomeFa=None, inputFa=None):
        self.inputBed = inputBed
        self.genomeFa = genomeFa
        self.inputFa = inputFa
        self.features = None
        self.inputBedObj = None

    def addIsoSupport(self):
        myDict = {}
        with open(self.inputBed) as f:
            for line in f.readlines():
                infoList = line.strip("\n").split("\t")
                isoName = infoList[3]
                flCount = int(infoList[-4])
                ratioIsoToGene = float(infoList[-2])
                myDict[isoName] = {"flCount": flCount, "ratioIsoToGene": ratioIsoToGene}

        newFeatures = pd.DataFrame.from_dict(myDict, orient="index")
        return pd.concat([self.features, newFeatures], axis=1)

    # def addIsoNgsSupp(self, isoNgsExp=None):
    #     self.features.update()
    #     return self.features

    def addJuncIndelInfo(self, samFile=None, junctionFile=None, refBed=None):
        juncDict = BedFile(junctionFile, type="bed12").getAllJuncDict()
        annoJuncDict = BedFile(refBed, type="bed12+").getAllJuncDict()
        sam = pysam.AlignmentFile(samFile, "r")

        indelsTotal = {}
        juncWithIndel = {}
        readCount = 0
        for read in sam.fetch():
            if read.is_unmapped: continue
            readCount += 1
            cigar = read.cigar
            pos = read.pos + 1
            spliceSites = []

            for cigarType, cigarLength in cigar:
                if cigarType not in [1, 4, 5, 7, 8]:
                    posEnd = pos + cigarLength - 1
                    if cigarType == 3:
                        # juncName = "{}:{}-{}".format(read.reference_id, pos, posEnd)
                        # if juncName in juncDict:
                        #     spliceSites.append([pos, posEnd])
                        spliceSites.append([pos, posEnd])
                    pos = posEnd + 1

            pos = read.pos + 1
            posEnd = pos
            for cigarType, cigarLength in cigar:
                if cigarType not in [1, 4, 5, 7, 8]:
                    posEnd = pos + cigarLength - 1

                if cigarType == 1:
                    indelStartPos = pos
                    indelEndPos = pos
                    indelLength = cigarLength
                    indelName = "{}:{}:{}+{}".format(read.reference_name, "insertion", indelStartPos, indelLength)

                    if indelName not in indelsTotal:
                        indelsTotal[indelName] = 1
                    else:
                        indelsTotal[indelName] += 1

                    for i in spliceSites:
                        juncName = "{}:{}-{}".format(read.reference_name, i[0], i[1])
                        if any(abs(indelStartPos-e) < 10 for e in i) or any(abs(indelEndPos-e) < 10 for e in i):
                            if juncName not in juncWithIndel:
                                juncWithIndel[juncName] = [indelName]
                            else:
                                juncWithIndel[juncName].append(indelName)

                if cigarType == 2:
                    indelStartPos = pos
                    indelEndPos = posEnd
                    indelLength = cigarLength
                    indelName = "{}:{}:{}-{}".format(read.reference_name, "deletion", indelStartPos, indelLength)

                    if indelName not in indelsTotal:
                        indelsTotal[indelName] = 1
                    else:
                        indelsTotal[indelName] += 1

                    for i in spliceSites:
                        juncName = "{}:{}-{}".format(read.reference_name, i[0], i[1])
                        if any(abs(indelStartPos-e) < 10 for e in i) or any(abs(indelEndPos-e) < 10 for e in i):
                            if juncName not in juncWithIndel:
                                juncWithIndel[juncName] = [indelName]
                            else:
                                juncWithIndel[juncName].append(indelName)

                if cigarType not in [1, 4, 5, 7, 8]:
                    pos = posEnd + 1
        sam.close()

        myDict = {}
        for isoName in self.inputBedObj.reads:
            if len(self.inputBedObj.reads[isoName].exons) > 1:
                myDict[isoName] = {
                    "nIndelsAroundJunc": 0,
                    "nJuncsWithIndels": 0,
                    "indelNearJunc": 0,
                    "ratioMinJuncCovToAllCov": 0,
                    "sdJuncCov": 0,
                    "minJuncRPKM": 0,
                    "withNovelJunc": False,
                    "minNovelJuncRPKM": 0
                }
                isoObj = self.inputBedObj.reads[isoName]
                isoLength = sum(isoObj.blockSizes)
                tmpDict = {}
                for junc in isoObj.introns:
                    juncName = "{}:{}-{}".format(isoObj.chrom, junc[0]+1, junc[1])
                    if juncName in juncDict:
                        juncCov = int(juncDict[juncName].split("\t")[4])
                        indelsNearJunc = list(set(juncWithIndel[juncName])) if juncName in juncWithIndel else []
                        validIndelsNearJunc = [indel for indel in indelsNearJunc if float(indelsTotal[indel]/juncCov) >= 0.5]
                        juncRPKM = juncCov / (isoLength/1000.0 * readCount/1000000.0)
                        juncAnno = True if juncName in annoJuncDict else False
                        tmpDict.update({juncName: {"juncCov": juncCov, "validIndels": validIndelsNearJunc,
                                                   "juncRPKM": juncRPKM, "juncAnno": juncAnno}})

                isoJuncCovs = [tmpDict[i]["juncCov"] for i in tmpDict]
                novelJuncRPKM = [tmpDict[i]["juncRPKM"] for i in tmpDict if tmpDict[i]["juncAnno"] == False]
                isoValidIndels = list(set(itertools.chain.from_iterable([tmpDict[i]["validIndels"] for i in tmpDict])))
                if isoValidIndels:
                    myDict[isoName]["indelNearJunc"] = 1
                    myDict[isoName]["nIndelsAroundJunc"] = len(isoValidIndels)
                    myDict[isoName]["nJuncsWithIndels"] = len([i for i in tmpDict if len(tmpDict[i]["validIndels"])])
                if isoJuncCovs:
                    myDict[isoName]["ratioMinJuncCovToAllCov"] = float(min(isoJuncCovs)) / sum(isoJuncCovs)
                    myDict[isoName]["sdJuncCov"] = np.std(isoJuncCovs)
                    myDict[isoName]["minJuncRPKM"] = min([tmpDict[i]["juncRPKM"] for i in tmpDict])
                if novelJuncRPKM:
                    myDict[isoName]["withNovelJunc"] = True
                    myDict[isoName]["minNovelJuncRPKM"] = min(novelJuncRPKM)


        newFeatures = pd.DataFrame.from_dict(myDict, orient="index")
        return pd.concat([self.features, newFeatures], axis=1)


    def getInputFa(self):
        tmpInputFa = os.path.join(os.getcwd(), "tmpInput.fa")
        cmd = '''cut -f 1-12 {} | 
            bedtools getfasta -fi {} -bed - -name -split -s | 
            seqkit replace -w 0 -p "(.*?):(.*)" -r '$1' > {}
        '''.format(self.inputBed, self.genomeFa, tmpInputFa)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        return tmpInputFa


    def getFeatures(self, refBed):
        from CPC2 import FindCDS
        self.inputFa = self.inputFa if self.inputFa else self.getInputFa()
        inputFaDict = SeqIO.to_dict(SeqIO.parse(self.inputFa, "fasta"))

        self.inputBedObj = BedFile(self.inputBed, type="bed12+")
        refBedObj = BedFile(refBed, type="bed12+")

        allAnnoJuncDict = refBedObj.getAllJuncDict()
        allAnnoExonDict = refBedObj.getAllExonDict()
        isoJuncChainDict = self.inputBedObj.getJuncChainDict()
        isos = self.inputBedObj.reads
        canonicalJuncFreq = self.getCanonicalJuncFreq(isoJuncChainDict, self.genomeFa)

        myDict = {}
        novelJuncList = []
        for isoName in isos:
            myDict[isoName] = {}
            juncChain = "{}:{}".format(isos[isoName].chrom, isos[isoName].juncChain)
            strand = isos[isoName].strand

            isoLength = sum(isos[isoName].blockSizes)
            exonNum = len(isos[isoName].blockSizes)
            annotation = isos[isoName].otherList[-1]
            orfSeq, startPos, orfStrand, orfFullness = FindCDS(inputFaDict[isoName].seq.upper()).longest_orf(strand)
            orfLen = len(orfSeq)
            canJuncRatio = canonicalJuncFreq[juncChain]
            gc = GC(inputFaDict[isoName].seq)
            myDict[isoName].update({
                "isoLength": isoLength,
                "exonNum": exonNum,
                "annotation": annotation,
                "orfLength": orfLen,
                "canJuncRatio": canJuncRatio,
                "GC": gc,
                "bite": False
            })

            for junc in isos[isoName].introns:
                juncStart, juncEnd = junc
                juncName = "{}:{}-{}".format(isos[isoName].chrom, juncStart+1, juncEnd)
                if juncName not in allAnnoJuncDict:
                    juncInfo = "\t".join(map(str, [isos[isoName].chrom, juncStart+1, juncEnd, isoName, ".", strand]))
                    novelJuncList.append(juncInfo)

        novelJuncObj = pybedtools.BedTool("\n".join(novelJuncList), from_string=True)
        annoExonObj = pybedtools.BedTool("\n".join(allAnnoExonDict.values()), from_string=True)
        intersectRes = novelJuncObj.intersect(annoExonObj, wa=True, wb=True, s=True)
        biteIsoforms = {}
        for i in intersectRes:
            infoList = str(i).strip("\n").split("\t")
            if infoList[3] not in biteIsoforms:
                biteIsoforms[infoList[3]] = ""

        for i in biteIsoforms:
            myDict[i]["bite"] = True

        return pd.DataFrame.from_dict(myDict, orient="index")


    def getCanonicalJuncFreq(self, isoJuncChainDict, genomeFasta):
        dinucleotideBedList = []
        for juncChain in isoJuncChainDict:
            for j in juncChain.split(":")[1].split(";"):
                chrom = isoJuncChainDict[juncChain].chrom
                strand = isoJuncChainDict[juncChain].strand
                juncStart, juncEnd = map(int, j.split("-"))
                juncName = "{}:{}-{}".format(chrom, juncStart+1, juncEnd)
                leftDinucleotidePos = "\t".join(map(str, [chrom, juncStart, juncStart + 2, ":".join([juncName, "left"]), ".", strand]))
                rightDinucleotidePos = "\t".join(map(str, [chrom, juncEnd - 2, juncEnd, ":".join([juncName, "right"]), ".", strand]))
                dinucleotideBedList.extend([leftDinucleotidePos, rightDinucleotidePos])

        dinucleotideBedObj = pybedtools.BedTool("\n".join(dinucleotideBedList), from_string=True)
        dinucleotideBedRes = dinucleotideBedObj.sequence(genomeFasta, name=True, tab=True, s=True)
        juncCanonicalDict = {}
        for i in str(open(dinucleotideBedRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            juncName, dinucleotideType = ":".join(infoList[0].split(":")[:2]), infoList[0].split(":")[2]

            if juncName not in juncCanonicalDict:
                juncCanonicalDict.update({juncName: {dinucleotideType: infoList[1]}})
            else:
                juncCanonicalDict[juncName].update({dinucleotideType: infoList[1]})

        canonicalJuncFreq = {}
        for juncChain in isoJuncChainDict:
            juncCount = len(juncChain.split(":")[1].split(";"))
            strand = isoJuncChainDict[juncChain].strand
            canonicalCount = 0
            for j in juncChain.split(":")[1].split(";"):
                juncStart = int(j.split("-")[0])+1
                juncEnd = int(j.split("-")[1])
                juncName = "{}:{}-{}".format(isoJuncChainDict[juncChain].chrom, juncStart, juncEnd)
                if strand == "+":
                    spliceMotif = "{}-{}".format(juncCanonicalDict[juncName]["left"],
                                                 juncCanonicalDict[juncName]["right"])
                else:
                    spliceMotif = "{}-{}".format(juncCanonicalDict[juncName]["right"],
                                                 juncCanonicalDict[juncName]["left"])
                if spliceMotif in ["GT-AG", "GC-AG", "AT-AC"]:
                    canonicalCount += 1
            canonicalJuncFreq[juncChain] = float(canonicalCount) / juncCount
        return canonicalJuncFreq


    def predCodingP(self, refBed):
        self.inputFa = self.inputFa if self.inputFa else self.getInputFa()
        self.features = self.features if not self.features.empty else self.getFeatures(refBed)

        from CPC2 import calculate_potential
        calculate_potential(self.inputFa, "+", 0, "cpc2out")
        codingPotential = pd.read_csv("cpc2out.txt", sep="\t", skiprows=1, header=None,
                                      usecols=[0, 2, 3, 4, 5, 6, 7], index_col=0)
        codingPotential.set_axis(["pepLength", "FickettScore", "pI", "orfIntegrity", "codingP", "codingLabel"], axis=1, inplace=True)

        return pd.concat([self.features, codingPotential], axis=1)


############################
# def getBlockLength(blockList):
#     return sum(map(lambda x: int(x[1]) - int(x[0]), blockList))
#
#
# def getDictFromFile(myFile, sep="\t", inlineSep=None, keyCol=1, valueCol=None):
#     with open(myFile) as f:
#         myDict = {}
#         for line in f.readlines():
#             infoList = line.strip("\n").split(sep)
#             key = infoList[keyCol - 1]
#             if valueCol:
#                 if inlineSep:
#                     value = infoList[valueCol - 1].split(inlineSep)
#                 else:
#                     value = infoList[valueCol - 1]
#             else:
#                 value = infoList
#             myDict[key] = value
#         return myDict


def getInvolvedIsos(asFile, asType="IR"):
    isoDict = {}
    if asType == "SE":
        incIndex = 15
        excIndex = 17
    else:
        incIndex = 7
        excIndex = 9

    with open(asFile) as f:
        for i in f.readlines():
            infoList = i.strip().split("\t")
            incIsos = infoList[incIndex]
            excIsos = infoList[excIndex]
            for j in incIsos.split(","):
                if j not in isoDict:
                    isoDict[j] = {}
            for j in excIsos.split(","):
                if j not in isoDict:
                    isoDict[j] = {}
    return isoDict


def getInputBed(baseDir=None, refBed=None, useAsFile=False):
    isoformFile = os.path.join(baseDir, "refine", "tofu.collapsed.assigned.unambi.bed12+")
    collapsedTrans2reads = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    isoformBed = BedFile(isoformFile, type="bed12+").reads

    isoDict = {}
    if useAsFile:
        irFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "IR.confident.bed6+")
        seFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "SE.confident.bed12+")
        a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "A3SS.confident.bed6+")
        a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "A5SS.confident.bed6+")

        isoDict.update(getInvolvedIsos(irFile, asType="IR"))
        isoDict.update(getInvolvedIsos(seFile, asType="SE"))
        isoDict.update(getInvolvedIsos(a5ssFile, asType="A5SS"))
        isoDict.update(getInvolvedIsos(a3ssFile, asType="A3SS"))
    else:
        isoDict = dict.fromkeys(isoformBed.keys(), "")

    # annoJuncDict, annoSingleExonList = BedFile(refBed, type="bed12+").getJuncChainDict()
    annoJuncDict = BedFile(refBed, type="bed12+").getJuncChainDict()
    isoform2reads = getDictFromFile(collapsedTrans2reads, sep="\t", inlineSep=",", valueCol=2)

    juncDict = {}
    gene2isoDict = {}
    isoSingleExonList = []
    for iso in isoDict:
        gene = isoformBed[iso].otherList[0]
        if gene not in gene2isoDict:
            gene2isoDict[gene] = {"isos": [iso], "count": len(isoform2reads[iso])}
        else:
            gene2isoDict[gene]["isos"].append(iso)
            gene2isoDict[gene]["count"] += len(isoform2reads[iso])

        if isoformBed[iso].juncChain == "":
            isoSingleExonList.append("{}\t{}\t{}\t{}\t{}\t{}".format(isoformBed[iso].chrom, isoformBed[iso].chromStart,
                                                                     isoformBed[iso].chromEnd, iso, ".",
                                                                     isoformBed[iso].strand))
            continue
        juncInfo = "{}:{}".format(isoformBed[iso].chrom, isoformBed[iso].juncChain)
        if juncInfo not in juncDict:
            isoLength = getBlockLength(isoformBed[iso].exons)
            juncDict[juncInfo] = {"iso": [iso], "gene": [isoformBed[iso].otherList[0]], "longest": [iso, isoLength]}
        else:
            juncDict[juncInfo]["iso"].append(iso)
            juncDict[juncInfo]["gene"].append(isoformBed[iso].otherList[0])
            if getBlockLength(isoformBed[iso].exons) > juncDict[juncInfo]["longest"][1]:
                juncDict[juncInfo]["longest"][0] = iso
                juncDict[juncInfo]["longest"][1] = getBlockLength(isoformBed[iso].exons)

    tmpOut = open(os.path.join(os.getcwd(), "inputBed.tmp"), "w")
    for junc in juncDict:
        longestIso = juncDict[junc]["longest"][0]
        gene = isoformBed[longestIso].otherList[0]
        isoSupport = sum([len(isoform2reads[x]) for x in juncDict[junc]["iso"]])
        geneSupport = gene2isoDict[gene]["count"]
        ratio = float(isoSupport) / geneSupport
        annotation = "annotated" if junc in annoJuncDict else "novel"
        print >> tmpOut, "\t".join(map(str, [str(isoformBed[longestIso]), ",".join(juncDict[junc]["iso"]), isoSupport, geneSupport, ratio, annotation]))
    tmpOut.close()
    return os.path.join(os.getcwd(), "inputBed.tmp")


# def processFeatureData(featureData, crossValid=False):
#     usedFeatures = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength", "orfIntegrity",
#                     "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio", "sdJuncCov", "minJuncRPKM",
#                     "nIndelsAroundJunc", "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc", "label"]
#
#     annoIsoData = featureData.loc[featureData.annotation == "annotated", ]
#     novelIsoData = featureData.loc[featureData.annotation == "novel", ]
#
#     if crossValid == False:
#         ########### construct test dataframe
#         topN = int(len(novelIsoData) * 0.1)
#         leastNonCodingIsos = novelIsoData.loc[novelIsoData.codingLabel == "noncoding", ].sort_values(by=["codingP"]).head(topN)
#         leastExpIsos = novelIsoData.loc[(novelIsoData.withNovelJunc == True) &
#                                         (novelIsoData.ratioIsoToGene < 0.1) &
#                                         (novelIsoData.flCount < 5), ]
#         unreliableJuncIsos = novelIsoData.loc[(novelIsoData.withNovelJunc == True) &
#                                               (novelIsoData.minNovelJuncRPKM < 0.1), ]
#         validNegIndex = set(leastNonCodingIsos.index) & (set(leastExpIsos.index) | set(unreliableJuncIsos.index))
#
#         posIsoData = annoIsoData.loc[(annoIsoData.codingLabel == "coding") & (annoIsoData.flCount >= 5) &
#                                      (annoIsoData.ratioIsoToGene >= 0.5) & (annoIsoData.minJuncRPKM >= 0.1), ]
#         validPosIndex = np.random.choice(posIsoData.index, replace=False, size=len(validNegIndex))
#
#         testNegData = featureData.loc[validNegIndex, :]
#         testNegData["label"] = -1
#         testPosData = featureData.loc[validPosIndex, :]
#         testPosData["label"] = +1
#         testData = pd.concat([testPosData, testNegData])
#         testData = testData.loc[:, usedFeatures]
#
#         ########### construct training dataframe
#         trainingPosData = posIsoData.loc[list(set(posIsoData.index) - set(validPosIndex)), ]
#         trainingUnlabelData = novelIsoData
#         trainingData = pd.concat([trainingPosData, trainingUnlabelData])
#         trainingData["label"] = 1
#         trainingData.loc[trainingData.annotation == "novel", "label"] = 0
#         trainingData = trainingData.loc[:, usedFeatures]
#
#         trainingData = trainingData.sample(frac=1)
#         testData = testData.sample(frac=1)
#
#         return trainingData, testData
#     else:
#         from sklearn.model_selection import KFold
#         posIsoData = annoIsoData.loc[(annoIsoData.codingLabel == "coding") & (annoIsoData.flCount >= 5) &
#                                      (annoIsoData.ratioIsoToGene >= 0.5) & (annoIsoData.minJuncRPKM >= 0.1),]
#         kf = KFold(n_splits=10, random_state=42, shuffle=True)
#
#         posTrainDataList = []
#         posTestDataList = []
#         unlabelTrainDataList = []
#         unlabelTestDataList = []
#         for train_index, test_index in kf.split(posIsoData):
#             posTrainDataList.append(posIsoData.iloc[train_index, ])
#             posTestDataList.append(posIsoData.iloc[test_index, ])
#
#         for train_index, test_index in kf.split(novelIsoData):
#             unlabelTrainDataList.append(novelIsoData.iloc[train_index, ])
#             unlabelTestDataList.append(novelIsoData.iloc[test_index, ])
#
#         trainingDataList = []
#         testDataList = []
#         for i in zip(posTrainDataList, unlabelTrainDataList):
#             tmpTrainingData = pd.concat(i)
#             tmpTrainingData["label"] = 1
#             tmpTrainingData.loc[tmpTrainingData.annotation == "novel", "label"] = 0
#             tmpTrainingData = tmpTrainingData.loc[:, usedFeatures]
#
#             tmpTrainingData = tmpTrainingData.sample(frac=1)
#             trainingDataList.append(tmpTrainingData)
#         for i in zip(posTestDataList, unlabelTestDataList):
#             tmpTestData = pd.concat(i)
#             tmpTestData["label"] = 1
#             tmpTestData.loc[tmpTestData.annotation == "novel", "label"] = 0
#             tmpTestData = tmpTestData.loc[:, usedFeatures]
#
#             tmpTestData = tmpTestData.sample(frac=1)
#             testDataList.append(tmpTestData)
#
#         return trainingDataList, testDataList
#
#
# def processFeatureData1(featureData, crossValid=False):
#     usedFeatures = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength", "orfIntegrity",
#                     "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio", "sdJuncCov", "minJuncRPKM",
#                     "nIndelsAroundJunc", "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc", "label"]
#
#     annoIsoData = featureData.loc[featureData.annotation == "annotated", ]
#     novelIsoData = featureData.loc[featureData.annotation == "novel", ]
#
#     if crossValid == False:
#         ########### construct test dataframe
#         topN = int(len(novelIsoData) * 0.1)
#         leastNonCodingIsos = novelIsoData.loc[novelIsoData.codingLabel == "noncoding", ].sort_values(by=["codingP"]).head(topN)
#         leastExpIsos = novelIsoData.loc[(novelIsoData.withNovelJunc == True) &
#                                         (novelIsoData.ratioIsoToGene < 0.1) &
#                                         (novelIsoData.flCount < 5), ]
#         unreliableJuncIsos = novelIsoData.loc[(novelIsoData.withNovelJunc == True) &
#                                               (novelIsoData.minNovelJuncRPKM < 0.1), ]
#         validNegIndex = set(leastNonCodingIsos.index) & (set(leastExpIsos.index) | set(unreliableJuncIsos.index))
#
#         posIsoData = annoIsoData.loc[(annoIsoData.codingLabel == "coding") & (annoIsoData.flCount >= 5) &
#                                      (annoIsoData.ratioIsoToGene >= 0.5) & (annoIsoData.minJuncRPKM >= 0.1), ]
#         validPosIndex = np.random.choice(posIsoData.index, replace=False, size=len(validNegIndex))
#
#         testNegData = featureData.loc[validNegIndex, :]
#         testNegData["label"] = -1
#         testPosData = featureData.loc[validPosIndex, :]
#         testPosData["label"] = +1
#         testData = pd.concat([testPosData, testNegData])
#         testData = testData.loc[:, usedFeatures]
#
#         ########### construct training dataframe
#         trainingPosData = posIsoData.loc[list(set(posIsoData.index) - set(validPosIndex)), ]
#         trainingUnlabelData = novelIsoData
#         trainingData = pd.concat([trainingPosData, trainingUnlabelData])
#         trainingData["label"] = 1
#         trainingData.loc[trainingData.annotation == "novel", "label"] = 0
#         trainingData = trainingData.loc[:, usedFeatures]
#
#         trainingData = trainingData.sample(frac=1)
#         testData = testData.sample(frac=1)
#
#         return trainingData, testData
#     else:
#         from sklearn.model_selection import KFold
#         posIsoData = annoIsoData.loc[(annoIsoData.codingLabel == "coding") & (annoIsoData.flCount >= 5) &
#                                      (annoIsoData.ratioIsoToGene >= 0.5) & (annoIsoData.minJuncRPKM >= 0.1),]
#         kf = KFold(n_splits=10, random_state=42, shuffle=True)
#
#         posTrainDataList = []
#         posTestDataList = []
#         unlabelTrainDataList = []
#         unlabelTestDataList = []
#         for train_index, test_index in kf.split(posIsoData):
#             posTrainDataList.append(posIsoData.iloc[train_index, ])
#             posTestDataList.append(posIsoData.iloc[test_index, ])
#
#         for train_index, test_index in kf.split(novelIsoData):
#             unlabelTrainDataList.append(novelIsoData.iloc[train_index, ])
#             unlabelTestDataList.append(novelIsoData.iloc[test_index, ])
#
#         trainingDataList = []
#         testDataList = []
#         for i in zip(posTrainDataList, unlabelTrainDataList):
#             tmpTrainingData = pd.concat(i)
#             tmpTrainingData["label"] = 1
#             tmpTrainingData.loc[tmpTrainingData.annotation == "novel", "label"] = 0
#             tmpTrainingData = tmpTrainingData.loc[:, usedFeatures]
#
#             tmpTrainingData = tmpTrainingData.sample(frac=1)
#             trainingDataList.append(tmpTrainingData)
#         for i in zip(posTestDataList, unlabelTestDataList):
#             tmpTestData = pd.concat(i)
#             tmpTestData["label"] = 1
#             tmpTestData.loc[tmpTestData.annotation == "novel", "label"] = 0
#             tmpTestData = tmpTestData.loc[:, usedFeatures]
#
#             tmpTestData = tmpTestData.sample(frac=1)
#             testDataList.append(tmpTestData)
#
#         return trainingDataList, testDataList


def processInputFile(inputBed, genomeFa, refBed, junctionFile=None, samFile=None):
    gf = GetFeatures(inputBed=inputBed, genomeFa=genomeFa)
    gf.features = gf.getFeatures(refBed)
    gf.features = gf.predCodingP(refBed)
    gf.features = gf.addIsoSupport()

    # if isoNgsExpFile:
    #     features = gf.addIsoNgsSupp(isoNgsExpFile)

    if samFile:
        gf.features = gf.addJuncIndelInfo(samFile=samFile, junctionFile=junctionFile, refBed=refBed)

    gf.features = gf.features.round(3)
    colOrder = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "annotation", "GC", "orfLength", "orfIntegrity",
                "pepLength", "FickettScore", "pI", "codingP", "codingLabel", "canJuncRatio", "sdJuncCov",
                "withNovelJunc", "minNovelJuncRPKM", "bite", "minJuncRPKM", "nIndelsAroundJunc",
                "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc"]
    gf.features = gf.features.loc[:, colOrder]
    gf.features.to_csv("isoFeatures.txt", sep="\t", index=True, index_label="isoform")

    return gf.features


def selectBestModel(featureData, drawAUC=False):
    from sklearn.model_selection import KFold

    usedFeatures = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength", "orfIntegrity",
                    "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio", "sdJuncCov", "minJuncRPKM",
                    "nIndelsAroundJunc", "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc", "label"]

    annoIsoData = featureData.loc[featureData.annotation == "annotated",]
    novelIsoData = featureData.loc[featureData.annotation == "novel",]

    leastExpIsos = novelIsoData.loc[(novelIsoData.ratioIsoToGene < 0.05),]
    unreliableJuncIsos = novelIsoData.loc[(novelIsoData.withNovelJunc == True) &
                                          (novelIsoData.minNovelJuncRPKM < 0.05),]
    validNegIndex = list(set(leastExpIsos.index) | set(unreliableJuncIsos.index))

    posIsoData = annoIsoData.loc[(annoIsoData.flCount >= 2) & (annoIsoData.minJuncRPKM >= 0.05),]
    posIsoDataInner, posIsoDataOuter = train_test_split(posIsoData, test_size=0.2, random_state=0)
    validNegIndexInner, validNegOuter = train_test_split(novelIsoData.loc[validNegIndex,], test_size=0.2,
                                                         random_state=0)

    outerPosData = posIsoDataInner.copy()
    outerPosData["label"] = 1
    posIsoDataOuterCopy = posIsoDataOuter.copy()
    posIsoDataOuterCopy["label"] = 1
    validNegOuterCopy = validNegOuter.copy()
    validNegOuterCopy["label"] = 0

    outerUnlabeledData = pd.concat([posIsoDataOuterCopy, validNegOuterCopy])
    outerUnlabeledData["label"] = 0

    outerPosData = outerPosData.loc[:, usedFeatures]
    outerPosData = outerPosData.sample(frac=1)
    outerUnlabeledData = outerUnlabeledData.loc[:, usedFeatures]
    outerUnlabeledData = outerUnlabeledData.sample(frac=1)

    outerPosData.replace({False: 0, True: 1}, inplace=True)
    outerUnlabeledData.replace({False: 0, True: 1}, inplace=True)
    outerTrainingData = pd.concat([outerPosData, outerUnlabeledData])

    kf = KFold(n_splits=5, random_state=42, shuffle=True)
    colors = ["#1f497d", "#f79646", "#9bbb59", "#7f7f7f", "#8064a2"]
    models = ["RF", "GB", "DT", "SVM", "NB"]
    names = ["RF", "GB", "DT", "SVM", "NB"]

    aucScore = dict.fromkeys(models, {})
    for scheme in ["bagging"]:
        if drawAUC:
            fig1 = plt.figure(figsize=[6, 6])
            ax1 = fig1.add_subplot(111, aspect='equal')
        for model, color, name in zip(models, colors, names):
            tprs = []
            # aucs = []
            mean_fpr = np.linspace(0, 1, 100)
            # i = 1
            for train_i, test_i in kf.split(posIsoDataInner):
                train_index = posIsoDataInner.index[train_i]
                test_index = posIsoDataInner.index[test_i]
                posIsoDataTraining = posIsoDataInner.loc[train_index,]
                posIsoDataTest = posIsoDataInner.loc[test_index,]
                negIsoDataTest = novelIsoData.loc[validNegIndexInner.index,]

                posIsoDataTraining["label"] = 1
                posIsoDataTest["label"] = 0
                negIsoDataTest["label"] = 0

                trainingData = pd.concat([posIsoDataTraining, posIsoDataTest, negIsoDataTest])
                trainingData = trainingData.loc[:, usedFeatures]
                trainingData = trainingData.sample(frac=1)

                posIsoDataTest["label"] = 1
                testData = pd.concat([posIsoDataTest, negIsoDataTest])
                testData = testData.loc[:, usedFeatures]
                testData = testData.sample(frac=1)

                trainingData.replace({False: 0, True: 1}, inplace=True)
                testData.replace({False: 0, True: 1}, inplace=True)

                sme = SchemeModelEval(trainingData=trainingData, model=model, scheme=scheme)
                sme.eval()
                prediction = sme.finalEstimator.predict_proba(testData.iloc[:, :-1])
                pu_score = pd.DataFrame({"pu_score": prediction[:, 1]}, index=testData.index)
                results = pd.DataFrame({
                    "true_label": testData.label,
                    "train_label": trainingData.loc[trainingData.label == 0,].label,
                    "pu_score": pu_score.pu_score
                }, columns=["true_label", "train_label", "pu_score"])

                fpr, tpr, t = roc_curve(results.true_label, results.pu_score)
                tprs.append(interp(mean_fpr, fpr, tpr))
                # roc_auc = auc(fpr, tpr)
                # aucs.append(roc_auc)
                # i = i + 1

            mean_tpr = np.mean(tprs, axis=0)
            mean_auc = auc(mean_fpr, mean_tpr)
            aucScore[model].update({"trainingAUC": mean_auc})
            if drawAUC:
                plt.plot(mean_fpr, mean_tpr, lw=2, linestyle='--', color=color,
                         label='%s (CV AUC = %0.3f )' % (name, mean_auc))

        ######################

        for model, color, name in zip(models, colors, names):
            sme = SchemeModelEval(trainingData=outerTrainingData, model=model, scheme=scheme)
            sme.eval()
            prediction = sme.finalEstimator.predict_proba(outerUnlabeledData.iloc[:, :-1])
            pu_score = pd.DataFrame({"pu_score": prediction[:, 1]}, index=outerUnlabeledData.index)
            results = pd.DataFrame({
                "true_label": pd.concat([posIsoDataOuterCopy, validNegOuterCopy]).label,
                "train_label": outerUnlabeledData.label,
                "pu_score": pu_score.pu_score
            }, columns=["true_label", "train_label", "pu_score"])

            fpr, tpr, t = roc_curve(results.true_label, results.pu_score)
            roc_auc = auc(fpr, tpr)
            aucScore[model].update({"testAUC": roc_auc})
            if drawAUC:
                plt.plot(fpr, tpr, lw=2, color=color, label='%s (Test AUC = %0.3f )' % (name, roc_auc))

        # plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='grey')
        if drawAUC:
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('ROC_AUC')
            plt.legend(loc="lower right")

            plt.savefig('PU.{}.ROC_AUC.pdf'.format(scheme))
            plt.close()

    tmpScore = 0
    tmpModel = None
    for i in aucScore:
        if np.mean(aucScore[i].values()) > tmpScore:
            tmpModel = i
            tmpScore = np.mean(aucScore[i].values())
    return tmpModel

# def iso_pu():
#
#     inputBed = ""
#     genomeFa = ""
#
#     isoSuppFile = "" if "" else None
#     isoNgsExpFile = "" if "" else None
#     samFile = "" if "" else None
#     junctionFile = "" if "" else None
#
#     baseDir = ""
#     refBed = ""
#     inputBed = inputBed if inputBed else getInputBed(baseDir=baseDir, refBed=refBed)
#     featureData = processInputFile(inputBed, genomeFa, refBed, junctionFile=junctionFile, samFile=samFile)
#     # featureData = pd.read_csv("isoFeatures.txt", sep="\t", index_col=0)
#
#     trainingData, testData = processFeatureData(featureData, crossValid=False)
#     trainingData.replace({False: 0, True: 1}, inplace=True)
#     testData.replace({False: 0, True: 1}, inplace=True)
#     # for scheme in ["bagging", "2step", "None"]:
#     #     for model in ["SVM", "RF", "GB", "ANN", "NB"]:
#     #         sme = SchemeModelEval(featData=featData, model=model, scheme=scheme)
#     #         sme.eval()
#     sme = SchemeModelEval(trainingData=trainingData, model="GB", scheme="bagging")
#     sme.eval()
#     sme.estimatorTest(testData=testData)
#     sme.filterIsoformsByScore(filterScore=0.8, outFile="validIsoforms.lst")
#
#     cvTrainingDataList, cvTestDataList = processFeatureData(featureData, crossValid=True)
#     for i in range(len(cvTrainingDataList)):
#         sme = SchemeModelEval(trainingData=cvTrainingDataList[i], model="GB", scheme="bagging")
#         sme.eval()
#         sme.estimatorTest(testData=cvTestDataList[i])
#         outFile = "validIsoforms.cv_{}.lst".format(i)
#         sme.filterIsoformsByScore(filterScore=0.8, outFile=outFile)

def iso_pu1(dataObj=None, dirSpec=None, refParams=None, hqIsoParams=None, samFile=None, junctionFile=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Start finding high quality novel isoforms for project {} sample {}...".format(projectName, sampleName)
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    prevDir = os.getcwd()
    hqIsoformDir = os.path.join(baseDir, "hqIsoforms")
    resolveDir(hqIsoformDir)

    genomeFa = refParams.ref_genome
    refBed = refParams.ref_bed

    inputBed = None

    if samFile == None and os.path.exists(os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "tmp.bam")):
        samFile = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "tmp.bam")
    if junctionFile == None and os.path.exists(os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")):
        junctionFile = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")

    inputBed = inputBed if inputBed else getInputBed(baseDir=baseDir, refBed=refBed)
    featureData = processInputFile(inputBed, genomeFa, refBed, junctionFile=junctionFile, samFile=samFile)
    # featureData = pd.read_csv("isoFeatures.txt", sep="\t", index_col=0)

    annoIsoData = featureData.loc[featureData.annotation == "annotated",]
    novelIsoData = featureData.loc[featureData.annotation == "novel",]
    posIsoData = annoIsoData.loc[(annoIsoData.flCount >= int(hqIsoParams.pos_fl_coverage)) &
                                 (annoIsoData.minJuncRPKM >= float(hqIsoParams.pos_min_junc_rpkm)),]
    posIsoDataCopy = posIsoData.copy()
    posIsoDataCopy["label"] = 1
    novelIsoDataCopy = novelIsoData.copy()
    novelIsoDataCopy["label"] = 0

    usedFeatures = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength", "orfIntegrity",
                    "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio", "sdJuncCov", "minJuncRPKM",
                    "nIndelsAroundJunc", "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc", "label"]
    dataToClassify = pd.concat([posIsoDataCopy, novelIsoDataCopy])
    dataToClassify = dataToClassify.loc[:, usedFeatures]
    dataToClassify = dataToClassify.sample(frac=1)
    dataToClassify.replace({False: 0, True: 1}, inplace=True)

    bestModel = selectBestModel(featureData, drawAUC=hqIsoParams.draw_auc)
    sme = SchemeModelEval(trainingData=dataToClassify, model=bestModel, scheme="bagging")
    # sme = SchemeModelEval(trainingData=dataToClassify, model="GB", scheme="bagging")
    sme.eval()
    sme.filterIsoformsByScore(filterScore=float(hqIsoParams.filter_score), outFile="validIsoforms.lst")

    cmd = "filter.pl -o validIsoforms.lst {} -2 4 -m i > hq.collapsed.bed12+".format(inputBed)
    subprocess.call(cmd, shell=True)

    os.chdir(prevDir)
    print getCurrentTime() + " Finding high quality novel isoforms for project {} sample {} done!".format(projectName, sampleName)

# def testModels():
#     from sklearn.model_selection import KFold
#     featureData = pd.read_csv("isoFeatures.txt", sep="\t", index_col=0)
#
#     usedFeatures = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength", "orfIntegrity",
#                     "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio", "sdJuncCov", "minJuncRPKM",
#                     "nIndelsAroundJunc", "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc", "label"]
#
#     annoIsoData = featureData.loc[featureData.annotation == "annotated", ]
#     novelIsoData = featureData.loc[featureData.annotation == "novel", ]
#
#
#     ########### construct test dataframe
#     # topN = int(len(novelIsoData) * 0.1)
#     # leastNonCodingIsos = novelIsoData.loc[novelIsoData.codingLabel == "noncoding", ].sort_values(by=["codingP"]).head(topN)
#     leastExpIsos = novelIsoData.loc[(novelIsoData.ratioIsoToGene < 0.05) &
#                                     (novelIsoData.flCount <= 2), ]
#     unreliableJuncIsos = novelIsoData.loc[(novelIsoData.withNovelJunc == True) &
#                                           (novelIsoData.minNovelJuncRPKM < 0.02), ]
#     validNegIndex = list(set(leastExpIsos.index) | set(unreliableJuncIsos.index))
#     selectedNegIndex = np.random.choice(validNegIndex, replace=False, size=800)
#
#     posIsoData = annoIsoData.loc[(annoIsoData.flCount >= 5) & (annoIsoData.ratioIsoToGene >= 0.5) &
#                                  (annoIsoData.minJuncRPKM >= 0.1), ]
#     # validPosIndex = np.random.choice(posIsoData.index, replace=False, size=len(validNegIndex))
#
#     kf = KFold(n_splits=5, random_state=42, shuffle=True)
#     # for model in ["RF", "DT", "GB", "SVM", "NB", "GB", "XGB", "LGB", "ERT", "ANN"]:
#     # for model in ["DT", "GB", "SVM", "NB"]:
#     for model in ["RF"]:
#         fig1 = plt.figure(figsize=[12, 12])
#         ax1 = fig1.add_subplot(111, aspect='equal')
#         ax1.add_patch(
#             patches.Arrow(0.45, 0.5, -0.25, 0.25, width=0.3, color='green', alpha=0.5)
#         )
#         ax1.add_patch(
#             patches.Arrow(0.5, 0.45, 0.25, -0.25, width=0.3, color='red', alpha=0.5)
#         )
#
#         tprs = []
#         aucs = []
#         mean_fpr = np.linspace(0, 1, 100)
#         i = 1
#         for train_i, test_i in kf.split(posIsoData):
#             train_index = posIsoData.index[train_i]
#             test_index = posIsoData.index[test_i]
#             posIsoDataTraining = posIsoData.loc[train_index, ]
#             posIsoDataTest = posIsoData.loc[test_index, ]
#             negIsoDataTest = novelIsoData.loc[selectedNegIndex, ]
#
#             posIsoDataTraining["label"] = 1
#             posIsoDataTest["label"] = 0
#             negIsoDataTest["label"] = 0
#
#             trainingData = pd.concat([posIsoDataTraining, posIsoDataTest, negIsoDataTest])
#             trainingData = trainingData.loc[:, usedFeatures]
#             trainingData = trainingData.sample(frac=1)
#
#             posIsoDataTest["label"] = 1
#             testData = pd.concat([posIsoDataTest, negIsoDataTest])
#             testData = testData.loc[:, usedFeatures]
#             testData = testData.sample(frac=1)
#
#             trainingData.replace({False: 0, True: 1}, inplace=True)
#             testData.replace({False: 0, True: 1}, inplace=True)
#
#             sme = SchemeModelEval(trainingData=trainingData, model=model, scheme="bagging")
#             sme.eval()
#             prediction = sme.finalEstimator.predict_proba(testData.iloc[:, :-1])
#             pu_score = pd.DataFrame({"pu_score": prediction[:, 1]}, index=testData.index)
#             results = pd.DataFrame({
#                 "true_label": testData.label,
#                 "train_label": trainingData.loc[trainingData.label == 0, ].label,
#                 "pu_score": pu_score.pu_score
#             }, columns=["true_label", "train_label", "pu_score"])
#
#             fpr, tpr, t = roc_curve(testData.iloc[:, -1], prediction[:, 1])
#             tprs.append(interp(mean_fpr, fpr, tpr))
#             roc_auc = auc(fpr, tpr)
#             aucs.append(roc_auc)
#             plt.plot(fpr, tpr, lw=2, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
#             i = i + 1
#
#         plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='black')
#         mean_tpr = np.mean(tprs, axis=0)
#         mean_auc = auc(mean_fpr, mean_tpr)
#         plt.plot(mean_fpr, mean_tpr, color='blue', label=r'Mean ROC (AUC = %0.2f )' % (mean_auc), lw=2, alpha=1)
#
#         plt.xlabel('False Positive Rate')
#         plt.ylabel('True Positive Rate')
#         plt.title('ROC')
#         plt.legend(loc="lower right")
#         plt.text(0.32, 0.7, 'More accurate area', fontsize=12)
#         plt.text(0.63, 0.4, 'Less accurate area', fontsize=12)
#         plt.savefig('PU.{}_{}.CV_ROC.pdf'.format(model, "bagging"))
#
#
# def testModels2():
#     from sklearn.model_selection import KFold
#     featureData = pd.read_csv("isoFeatures.txt", sep="\t", index_col=0)
#
#     usedFeatures = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength", "orfIntegrity",
#                     "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio", "sdJuncCov", "minJuncRPKM",
#                     "nIndelsAroundJunc", "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc", "label"]
#
#     annoIsoData = featureData.loc[featureData.annotation == "annotated", ]
#     novelIsoData = featureData.loc[featureData.annotation == "novel", ]
#
#
#     ########### construct test dataframe
#     # topN = int(len(novelIsoData) * 0.1)
#     # leastNonCodingIsos = novelIsoData.loc[novelIsoData.codingLabel == "noncoding", ].sort_values(by=["codingP"]).head(topN)
#     leastExpIsos = novelIsoData.loc[(novelIsoData.ratioIsoToGene < 0.05), ]
#     unreliableJuncIsos = novelIsoData.loc[(novelIsoData.withNovelJunc == True) &
#                                           (novelIsoData.minNovelJuncRPKM < 0.05), ]
#     validNegIndex = list(set(leastExpIsos.index) | set(unreliableJuncIsos.index))
#
#     posIsoData = annoIsoData.loc[(annoIsoData.flCount >= 4) & (annoIsoData.minJuncRPKM >= 0.1), ]
#     posIsoDataInner, posIsoDataOuter = train_test_split(posIsoData, test_size=0.2)
#     validNegIndexInner, validNegOuter = train_test_split(novelIsoData.loc[validNegIndex, ], test_size=0.2)
#
#     outerPosData = posIsoDataInner.copy()
#     posIsoDataOuterCopy = posIsoDataOuter.copy()
#     validNegOuterCopy = validNegOuter.copy()
#     posIsoDataOuterCopy["label"] = 1
#     validNegOuterCopy["label"] = 0
#
#     outerUnlabeledData = pd.concat([posIsoDataOuter, validNegOuter])
#     outerPosData["label"] = 1
#     outerUnlabeledData["label"] = 0
#
#     outerPosData = outerPosData.loc[:, usedFeatures]
#     outerPosData = outerPosData.sample(frac=1)
#     outerUnlabeledData = outerUnlabeledData.loc[:, usedFeatures]
#     outerUnlabeledData = outerUnlabeledData.sample(frac=1)
#
#     outerPosData.replace({False: 0, True: 1}, inplace=True)
#     outerUnlabeledData.replace({False: 0, True: 1}, inplace=True)
#     outerTrainingData = pd.concat([outerPosData, outerUnlabeledData])
#
#     kf = KFold(n_splits=5, random_state=42, shuffle=True)
#     # for model in ["RF", "DT", "GB", "SVM", "NB", "XGB", "LGB", "ERT", "ANN"]:
#     for model in ["RF", "DT", "GB", "SVM", "NB"]:
#     # for model in ["XGB", "LGB", "ERT", "ANN"]:
#     # for model in ["RF"]:
#         for scheme in ["bagging"]:
#             fig1 = plt.figure(figsize=[12, 12])
#             ax1 = fig1.add_subplot(111, aspect='equal')
#             ax1.add_patch(
#                 patches.Arrow(0.45, 0.5, -0.25, 0.25, width=0.3, color='green', alpha=0.5)
#             )
#             ax1.add_patch(
#                 patches.Arrow(0.5, 0.45, 0.25, -0.25, width=0.3, color='red', alpha=0.5)
#             )
#
#             tprs = []
#             aucs = []
#             mean_fpr = np.linspace(0, 1, 100)
#             i = 1
#             for train_i, test_i in kf.split(posIsoDataInner):
#                 train_index = posIsoDataInner.index[train_i]
#                 test_index = posIsoDataInner.index[test_i]
#                 posIsoDataTraining = posIsoDataInner.loc[train_index, ]
#                 posIsoDataTest = posIsoDataInner.loc[test_index, ]
#                 negIsoDataTest = novelIsoData.loc[validNegIndexInner.index, ]
#
#                 posIsoDataTraining["label"] = 1
#                 posIsoDataTest["label"] = 0
#                 negIsoDataTest["label"] = 0
#
#                 trainingData = pd.concat([posIsoDataTraining, posIsoDataTest, negIsoDataTest])
#                 trainingData = trainingData.loc[:, usedFeatures]
#                 trainingData = trainingData.sample(frac=1)
#
#                 posIsoDataTest["label"] = 1
#                 testData = pd.concat([posIsoDataTest, negIsoDataTest])
#                 testData = testData.loc[:, usedFeatures]
#                 testData = testData.sample(frac=1)
#
#                 trainingData.replace({False: 0, True: 1}, inplace=True)
#                 testData.replace({False: 0, True: 1}, inplace=True)
#
#                 sme = SchemeModelEval(trainingData=trainingData, model=model, scheme=scheme)
#                 sme.eval()
#                 prediction = sme.finalEstimator.predict_proba(testData.iloc[:, :-1])
#                 pu_score = pd.DataFrame({"pu_score": prediction[:, 1]}, index=testData.index)
#                 results = pd.DataFrame({
#                     "true_label": testData.label,
#                     "train_label": trainingData.loc[trainingData.label == 0, ].label,
#                     "pu_score": pu_score.pu_score
#                 }, columns=["true_label", "train_label", "pu_score"])
#
#                 fpr, tpr, t = roc_curve(testData.iloc[:, -1], prediction[:, 1])
#                 tprs.append(interp(mean_fpr, fpr, tpr))
#                 roc_auc = auc(fpr, tpr)
#                 aucs.append(roc_auc)
#                 plt.plot(fpr, tpr, lw=2, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
#                 i = i + 1
#
#             plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='black')
#             mean_tpr = np.mean(tprs, axis=0)
#             mean_auc = auc(mean_fpr, mean_tpr)
#             plt.plot(mean_fpr, mean_tpr, color='blue', label=r'Mean ROC (AUC = %0.2f )' % (mean_auc), lw=2, alpha=1)
#
#             plt.xlabel('False Positive Rate')
#             plt.ylabel('True Positive Rate')
#             plt.title('ROC')
#             plt.legend(loc="lower right")
#             plt.text(0.32, 0.7, 'More accurate area', fontsize=12)
#             plt.text(0.63, 0.4, 'Less accurate area', fontsize=12)
#             plt.savefig('PU.{}_{}.inner.CV_ROC.pdf'.format(model, scheme))
#             plt.close()
#
#             ######################
#             fig1 = plt.figure(figsize=[12, 12])
#             ax1 = fig1.add_subplot(111, aspect='equal')
#             ax1.add_patch(
#                 patches.Arrow(0.45, 0.5, -0.25, 0.25, width=0.3, color='green', alpha=0.5)
#             )
#             ax1.add_patch(
#                 patches.Arrow(0.5, 0.45, 0.25, -0.25, width=0.3, color='red', alpha=0.5)
#             )
#
#             sme = SchemeModelEval(trainingData=outerTrainingData, model=model, scheme=scheme)
#             sme.eval()
#             prediction = sme.finalEstimator.predict_proba(outerTrainingData.iloc[:, :-1])
#             pu_score = pd.DataFrame({"pu_score": prediction[:, 1]}, index=outerTrainingData.index)
#             results = pd.DataFrame({
#                 "true_label": pd.concat([outerPosData, posIsoDataOuterCopy, validNegOuterCopy]).label,
#                 "train_label": outerTrainingData.label,
#                 "pu_score": pu_score.pu_score
#             }, columns=["true_label", "train_label", "pu_score"])
#
#             fpr, tpr, t = roc_curve(outerTrainingData.iloc[:, -1], prediction[:, 1])
#             roc_auc = auc(fpr, tpr)
#             plt.plot(fpr, tpr, lw=2, alpha=0.3, label='ROC outer data (AUC = %0.2f )' % (roc_auc))
#             plt.xlabel('False Positive Rate')
#             plt.ylabel('True Positive Rate')
#             plt.title('ROC')
#             plt.legend(loc="lower right")
#             plt.text(0.32, 0.7, 'More accurate area', fontsize=12)
#             plt.text(0.63, 0.4, 'Less accurate area', fontsize=12)
#             plt.savefig('PU.{}_{}.outer.CV_ROC.pdf'.format(model, scheme))
#             plt.close()
#
#
# def testModels3():
#     from sklearn.model_selection import KFold
#     featureData = pd.read_csv("isoFeatures.txt", sep="\t", index_col=0)
#
#     usedFeatures = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength", "orfIntegrity",
#                     "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio", "sdJuncCov", "minJuncRPKM",
#                     "nIndelsAroundJunc", "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc", "label"]
#
#     annoIsoData = featureData.loc[featureData.annotation == "annotated",]
#     novelIsoData = featureData.loc[featureData.annotation == "novel",]
#
#     ########### construct test dataframe
#     # topN = int(len(novelIsoData) * 0.1)
#     # leastNonCodingIsos = novelIsoData.loc[novelIsoData.codingLabel == "noncoding", ].sort_values(by=["codingP"]).head(topN)
#     leastExpIsos = novelIsoData.loc[(novelIsoData.ratioIsoToGene < 0.05),]
#     unreliableJuncIsos = novelIsoData.loc[(novelIsoData.withNovelJunc == True) &
#                                           (novelIsoData.minNovelJuncRPKM < 0.05),]
#     validNegIndex = list(set(leastExpIsos.index) | set(unreliableJuncIsos.index))
#
#     posIsoData = annoIsoData.loc[(annoIsoData.flCount >= 2) & (annoIsoData.minJuncRPKM >= 0.05),]
#     posIsoDataInner, posIsoDataOuter = train_test_split(posIsoData, test_size=0.2, random_state=0)
#     validNegIndexInner, validNegOuter = train_test_split(novelIsoData.loc[validNegIndex,], test_size=0.2,
#                                                          random_state=0)
#
#     outerPosData = posIsoDataInner.copy()
#     outerPosData["label"] = 1
#     posIsoDataOuterCopy = posIsoDataOuter.copy()
#     posIsoDataOuterCopy["label"] = 1
#     validNegOuterCopy = validNegOuter.copy()
#     validNegOuterCopy["label"] = 0
#
#     outerUnlabeledData = pd.concat([posIsoDataOuterCopy, validNegOuterCopy])
#     outerUnlabeledData["label"] = 0
#
#     outerPosData = outerPosData.loc[:, usedFeatures]
#     outerPosData = outerPosData.sample(frac=1)
#     outerUnlabeledData = outerUnlabeledData.loc[:, usedFeatures]
#     outerUnlabeledData = outerUnlabeledData.sample(frac=1)
#
#     outerPosData.replace({False: 0, True: 1}, inplace=True)
#     outerUnlabeledData.replace({False: 0, True: 1}, inplace=True)
#     outerTrainingData = pd.concat([outerPosData, outerUnlabeledData])
#
#     kf = KFold(n_splits=5, random_state=42, shuffle=True)
#     colors = ["#1f497d", "#f79646", "#9bbb59", "#7f7f7f", "#8064a2"]
#     models = ["RF", "GB", "DT", "SVM", "NB"]
#     names = ["RF", "GB", "DT", "SVM", "NB"]
#
#     for scheme in ["bagging"]:
#         fig1 = plt.figure(figsize=[6, 6])
#         ax1 = fig1.add_subplot(111, aspect='equal')
#         for model, color, name in zip(models, colors, names):
#             tprs = []
#             aucs = []
#             mean_fpr = np.linspace(0, 1, 100)
#             i = 1
#             for train_i, test_i in kf.split(posIsoDataInner):
#                 train_index = posIsoDataInner.index[train_i]
#                 test_index = posIsoDataInner.index[test_i]
#                 posIsoDataTraining = posIsoDataInner.loc[train_index,]
#                 posIsoDataTest = posIsoDataInner.loc[test_index,]
#                 negIsoDataTest = novelIsoData.loc[validNegIndexInner.index,]
#
#                 posIsoDataTraining["label"] = 1
#                 posIsoDataTest["label"] = 0
#                 negIsoDataTest["label"] = 0
#
#                 trainingData = pd.concat([posIsoDataTraining, posIsoDataTest, negIsoDataTest])
#                 trainingData = trainingData.loc[:, usedFeatures]
#                 trainingData = trainingData.sample(frac=1)
#
#                 posIsoDataTest["label"] = 1
#                 testData = pd.concat([posIsoDataTest, negIsoDataTest])
#                 testData = testData.loc[:, usedFeatures]
#                 testData = testData.sample(frac=1)
#
#                 trainingData.replace({False: 0, True: 1}, inplace=True)
#                 testData.replace({False: 0, True: 1}, inplace=True)
#
#                 sme = SchemeModelEval(trainingData=trainingData, model=model, scheme=scheme)
#                 sme.eval()
#                 prediction = sme.finalEstimator.predict_proba(testData.iloc[:, :-1])
#                 pu_score = pd.DataFrame({"pu_score": prediction[:, 1]}, index=testData.index)
#                 results = pd.DataFrame({
#                     "true_label": testData.label,
#                     "train_label": trainingData.loc[trainingData.label == 0,].label,
#                     "pu_score": pu_score.pu_score
#                 }, columns=["true_label", "train_label", "pu_score"])
#
#                 fpr, tpr, t = roc_curve(results.true_label, results.pu_score)
#                 tprs.append(interp(mean_fpr, fpr, tpr))
#                 roc_auc = auc(fpr, tpr)
#                 aucs.append(roc_auc)
#                 i = i + 1
#
#             mean_tpr = np.mean(tprs, axis=0)
#             mean_auc = auc(mean_fpr, mean_tpr)
#             plt.plot(mean_fpr, mean_tpr, lw=2, linestyle='--', color=color, label='%s (CV AUC = %0.3f )' % (name, mean_auc))
#
#         ######################
#
#         for model, color, name in zip(models, colors, names):
#             sme = SchemeModelEval(trainingData=outerTrainingData, model=model, scheme=scheme)
#             sme.eval()
#             prediction = sme.finalEstimator.predict_proba(outerUnlabeledData.iloc[:, :-1])
#             pu_score = pd.DataFrame({"pu_score": prediction[:, 1]}, index=outerUnlabeledData.index)
#             results = pd.DataFrame({
#                 "true_label": pd.concat([posIsoDataOuterCopy, validNegOuterCopy]).label,
#                 "train_label": outerUnlabeledData.label,
#                 "pu_score": pu_score.pu_score
#             }, columns=["true_label", "train_label", "pu_score"])
#
#             fpr, tpr, t = roc_curve(results.true_label, results.pu_score)
#             roc_auc = auc(fpr, tpr)
#             plt.plot(fpr, tpr, lw=2, color=color, label='%s (Test AUC = %0.3f )' % (name, roc_auc))
#
#         # plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='grey')
#         plt.xlabel('False Positive Rate')
#         plt.ylabel('True Positive Rate')
#         plt.title('ROC_AUC')
#         plt.legend(loc="lower right")
#
#         plt.savefig('PU.{}.ROC_AUC.pdf'.format(scheme))
#         plt.close()
#
# # testModels()
# # testModels2()
# testModels3()