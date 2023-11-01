"""
Propensity scores estimation and matching.
"""

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import make_column_transformer
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline


def estimate_ps(df_cov: pd.DataFrame, y: np.ndarray, cat: list[str] | None = None) -> np.ndarray:
    """Estimate propensity scores from covariates and class label.

    Parameters
    ----------
    df_cov: pd.DataFrame of shape=(n_subjects, n_covariates)
        DataFrame holding subjects covariates.
    y: np.ndarray[int], shape=(n_subjects,)
        Binary class label (e.g., diagnosis)
    cat: list of str, default=None
        Categorical covariates in `df_cov` (e.g., site). If None, assumes no categorical
        columns in dataframe.

    Returns
    -------
    ps: np.ndarray of shape (n_subjects,)
        Estimated propensity scores.
    """

    cat = [] if cat is None else cat

    scale_keys = np.setdiff1d(df_cov.columns, cat)
    ct = make_column_transformer((StandardScaler(), scale_keys),
                                 (OneHotEncoder(drop='first'), cat))

    logit = LogisticRegression(penalty=None, random_state=0)
    clf = Pipeline([('ct', ct), ('clf', CalibratedClassifierCV(logit))])
    ps = clf.fit(df_cov, y).predict_proba(df_cov)[:, 1]

    return ps


def get_matches(ps: float, ps_cn: np.ndarray[float], caliper: float | None = .2, n_min: int | None = None,
                n_max: int | None = None):
    """ Get the indices of the closest propensity scores in `pc_cn` to `ps`.

    Parameters
    ----------
    ps: float
        Propensity score of target subject
    ps_cn: np.ndarray, shape = (n_subjects,)
        Propensity scores of control subjects
    caliper: float, default=0.2
        Caliper to use for imperfect matches. If None, no caliper is used.
    n_min: int, default=None
        Minimum number of matches.
    n_max: int default=None
        Maximum number of matches to consider.

    Returns
    -------
    matches: np.ndarray, shape=(n_matches,)
        Indices of the closest matches in `pc_cn`.
    """
    d = np.abs(ps_cn - ps)

    # Indices of closest matches in controls
    idx = np.argsort(d)

    if caliper is not None:
        thresh = caliper * d.std()
        n = np.count_nonzero(d < thresh)
        if n_min is not None:
            n_min = max(n, n_min)
        else:
            n_min = n
    elif n_min is None:
        n_min = d.size

    if n_max is not None:
        n_min = min(n_min, n_max)

    return idx[:n_min]
