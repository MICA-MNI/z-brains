"""
Propensity scores estimation, matching and stratification.
"""


from typing import List, Tuple, Optional

from scipy.optimize import linear_sum_assignment
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import make_column_transformer
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline


def estimate_ps(df_cov: pd.DataFrame, y: np.ndarray,
                cat: Optional[List[str]] = None) -> np.ndarray:
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

    logit = LogisticRegression(penalty='none', random_state=0)
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


def match_ps(y: np.ndarray, ps: np.ndarray, caliper: Optional[float] = None) \
        -> Tuple[pd.DataFrame, np.ndarray]:
    """Match propensity scores.

    Parameters
    ----------
    y: np.ndarray of shape (n_subjects,)
        Binary class label (e.g., diagnosis).
    ps: np.ndarray of shape (n_subjects,)
        Estimated propensity scores.
    caliper: float, default=None
        Caliper to use for imperfect matches. If None, no caliper is used.

    Returns
    -------
    df_match: pd.DataFrame
        DataFrame with matched subjects. Paired subjects share the same value
        in the 'pair' column.
    mask_match: np.ndarray
        Mask of subjects to keep after matching.
    """

    mask = y == 1  # mask positive class (e.g., ASD)
    d = np.abs(ps[mask][:, None] - ps[~mask])

    if caliper is not None:
        thresh = caliper * d.std()
        # set pairs >= caliper to a very large number
        d[d >= thresh] = 10000

    ridx, cidx = linear_sum_assignment(d)
    cost = d[ridx, cidx]

    # Discard matches above threshold
    if caliper is not None:
        m = cost < thresh
        ridx, cidx, cost = ridx[m], cidx[m], cost[m]

    # Convert to original indices
    idx_pos = np.arange(y.size)[mask]
    idx_neg = np.arange(y.size)[~mask]
    ridx = idx_pos[ridx]
    cidx = idx_neg[cidx]

    pair_sel = np.full(y.size, -1, dtype=int)
    pair_sel[ridx] = pair_sel[cidx] = np.arange(ridx.size)
    mask_match = pair_sel >= 0

    df_match = dict(zip(['group', 'ps', 'pair'], [y, ps, pair_sel]))
    df_match = pd.DataFrame(df_match)[mask_match].reset_index(drop=True)
    df_match['group'].replace({0: 'Neg', 1: 'Pos'}, inplace=True)
    df_match.index.name = 'sample_idx'

    return df_match, mask_match
