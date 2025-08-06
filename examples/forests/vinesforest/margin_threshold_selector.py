# import numpy as np
# import pandas as pd
# from sklearn.metrics import confusion_matrix, precision_recall_fscore_support
# from sklearn.model_selection import KFold


# class MarginThresholdSelector:
#   """
#   Select a ΔLL margin threshold τ to control harmful switches via CV precision.

#   X_alt_train, X_base_train: shape (n_models, n_train_obs)
#   X_alt_test,  X_base_test : shape (n_models, n_test_obs)

#   For each fold:
#     - train label (prediction) for model m:  mean(ΔLL_train_fold_m) >= τ
#     - test label (ground truth) for model m: mean(ΔLL_test_fold_m)  >  0
#   Aggregates across folds to estimate precision(τ). Picks τ to meet target_precision
#   (smallest τ achieving it; if none, τ with maximum precision).
#   """

#   def __init__(self, cv=5, target_precision=0.95, n_tau=200, random_state=0):
#     self.cv = cv
#     self.target_precision = target_precision
#     self.n_tau = n_tau
#     self.random_state = random_state
#     self.tau_ = None
#     self.cv_summary_ = None

#   @staticmethod
#   def _fold_means(X):
#     # X shape (n_models, n_obs_fold)
#     # returns per-model scalar mean
#     return X.mean(axis=1)

#   def fit(self, X_alt_train, X_base_train, X_alt_test, X_base_test):
#     rng = np.random.default_rng(self.random_state)

#     # ΔLL matrices
#     D_train = X_alt_train - X_base_train  # (n_models, n_train)
#     D_test = X_alt_test - X_base_test  # (n_models, n_test)

#     # Build CV over *observations* (not models)
#     n_train_obs = D_train.shape[1]
#     kf = KFold(n_splits=self.cv, shuffle=True, random_state=self.random_state)

#     # Candidate taus from quantiles of mean ΔLL on train (robust grid)
#     # Use means per model over entire train as a scale reference.
#     s_train_full = D_train.mean(axis=1)
#     qmin, qmax = np.quantile(s_train_full, [0.01, 0.99])
#     tau_grid = np.linspace(qmin, max(qmax, 0.0), self.n_tau)

#     # Accumulate predictions/labels across folds for each τ
#     precisions = np.zeros_like(tau_grid)
#     recalls = np.zeros_like(tau_grid)
#     supports_pos = np.zeros_like(tau_grid, dtype=int)
#     supports_tot = 0

#     for fold_idx, (tr_idx, te_idx) in enumerate(
#       kf.split(np.arange(n_train_obs))
#     ):
#       # Fold-wise train/test (over observations)
#       Dtr = D_train[:, tr_idx]  # (n_models, n_tr_fold)
#       Dte = D_train[
#         :, te_idx
#       ]  # (n_models, n_te_fold)  -> used as "validation test" for model selection
#       # For ground-truth, we want true generalization; use the held-out *global* test set for labels,
#       # but we still need fold consistency. To avoid leakage, reuse the same global D_test.
#       # Ground-truth label per model is whether it's better on the (global) test set:
#       y_true = self._fold_means(D_test) > 0.0  # (n_models,)

#       s_tr_mean = self._fold_means(Dtr)  # (n_models,)

#       for j, tau in enumerate(tau_grid):
#         y_pred = s_tr_mean >= tau  # predict switch if margin >= tau
#         if y_pred.sum() == 0:
#           # avoid zero-division: precision undefined -> set to 0
#           prec = 0.0
#           rec = 0.0
#           pos_support = int(y_true.sum())
#         else:
#           # Precision/recall against ground truth
#           prec, rec, _, _ = precision_recall_fscore_support(
#             y_true, y_pred, average="binary", zero_division=0
#           )
#           pos_support = int(y_true.sum())

#         precisions[j] += prec
#         recalls[j] += rec
#         supports_pos[j] += pos_support
#       supports_tot += 1

#     # Average over folds
#     precisions /= supports_tot
#     recalls /= supports_tot

#     # Choose τ: smallest τ achieving target precision; otherwise τ with max precision
#     feasible = np.where(precisions >= self.target_precision)[0]
#     if feasible.size > 0:
#       j_star = feasible[0]  # smallest τ meeting precision
#     else:
#       j_star = np.argmax(precisions)

#     self.tau_ = float(tau_grid[j_star])
#     self.cv_summary_ = pd.DataFrame(
#       {"tau": tau_grid, "cv_precision": precisions, "cv_recall": recalls}
#     )
#     return self

#   def decision_from_train(self, X_alt_train, X_base_train):
#     assert self.tau_ is not None, "Call fit() first."
#     D_train = X_alt_train - X_base_train
#     s_train_mean = D_train.mean(axis=1)
#     return s_train_mean >= self.tau_

#   @staticmethod
#   def evaluate_on_test(y_pred_from_train, X_alt_test, X_base_test):
#     D_test = X_alt_test - X_base_test
#     y_true = D_test.mean(axis=1) > 0.0
#     cm = confusion_matrix(y_true, y_pred_from_train)
#     prec, rec, f1, _ = precision_recall_fscore_support(
#       y_true, y_pred_from_train, average="binary", zero_division=0
#     )
#     acc = (y_true == y_pred_from_train).mean()
#     return {
#       "confusion_matrix": cm,
#       "precision": prec,
#       "recall": rec,
#       "f1": f1,
#       "accuracy": acc,
#       "support_pos": int(y_true.sum()),
#       "support_neg": int((~y_true).sum()),
#     }
