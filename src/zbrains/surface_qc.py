"""Pre-normalization, control-level outlier detection for surface feature maps.

Detects control subjects whose *raw* surface-sampled feature maps (``.func.gii``)
are likely processing or acquisition outliers, BEFORE intensity normalization
(WhiteStripe, RAVEL), harmonization (ComBat), or normative modelling. Input maps
are never modified.

Method (bilateral-only, per feature)
------------------------------------
1. **Bilateral joint robust scaling.** For each subject a single centre/scale is
   estimated from the *concatenated* left+right map
   (``c = median([xL, xR])``, ``s = 1.4826 * MAD([xL, xR])``) and applied to both
   hemispheres. This removes whole-scan arbitrary scaling while *preserving*
   abnormal relative scaling between hemispheres (independent per-hemisphere
   scaling could conceal a global attenuation affecting only one hemisphere).
2. **Cross-fitted robust PCA** on the bilateral concatenated QC maps. For each
   subject the robust model is fit WITHOUT that subject (K-fold / leave-one-out),
   then the subject is projected in and its **score distance** (robust
   Mahalanobis in PCA score space) and **orthogonal distance** (residual to the
   subspace) are computed. This avoids the self-influence that makes in-sample
   rank scores invalid. (A true ROBPCA, e.g. ``rrcov``, is preferred; when it is
   unavailable this module uses a documented ROBPCA-equivalent: median-centred
   SVD PCA with one reweighting pass and a robust score covariance.)
3. **Raw-intensity QC.** Robust deviations of raw-map summaries (median, MAD,
   dynamic range, fraction-zero, and left-right differences), so gross scaling /
   offset failures that the within-map scaling deliberately removes stay visible.
4. **One combined anomaly statistic.** Each component (orthogonal distance, score
   distance, raw abnormality) is standardized against its robust null, and a
   single statistic ``A = max(Z_OD, Z_SD, Z_raw)`` is formed (max, because a map
   should be flagged if *any* failure mode is extreme).
5. **One jointly-calibrated threshold.** The null distribution of ``A`` is
   obtained once by simulating from the robustly-fitted inlier component
   distribution, which automatically accounts for the correlation between the
   maximized components. A single user parameter ``qc_alpha`` sets the tail:
   ``flagged = A > quantile(A_null, 1 - qc_alpha)``. (This replaces the
   over-conservative empirical-p x Bonferroni x BH scheme, whose minimum
   attainable p was ~ constant/N and could never reach small ``qc_q``.)

If either hemisphere (via the bilateral map) is abnormal, BOTH hemisphere files
are marked for review. Statistical calculations are kept independent from figure
generation and the command-line interface.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

try:  # nibabel needed to read GIFTIs, not to import the pure helpers
    import nibabel as nib
except Exception:  # pragma: no cover
    nib = None

EPSILON = 1e-8
MAD_TO_SIGMA = 1.4826


# ---------------------------------------------------------------------------
# Numeric helpers
# ---------------------------------------------------------------------------
def _median_absolute_deviation(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return float(np.median(np.abs(x - np.median(x))))


def robust_within_map_scale(x: np.ndarray, epsilon: float = EPSILON):
    """Within-map robust standardization ``(x - median) / (1.4826*MAD)`` (utility)."""
    x = np.asarray(x, dtype=float)
    center = float(np.nanmedian(x)) if np.isfinite(x).any() else float("nan")
    mad = _median_absolute_deviation(x)
    scale = MAD_TO_SIGMA * mad if np.isfinite(mad) else float("nan")
    denom = max(scale, epsilon) if np.isfinite(scale) else epsilon
    return (x - center) / denom, center, (scale if np.isfinite(scale) else 0.0)


def robust_bilateral_scale(left: np.ndarray, right: np.ndarray, epsilon: float = EPSILON):
    """Joint bilateral robust standardization (Step 1).

    A single centre/scale from the concatenated left+right values is applied to
    both hemispheres, so global scan scaling is removed but a one-sided global
    attenuation is preserved. Returns ``(zL, zR, center, scale)``.
    """
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    both = np.concatenate([left, right])
    center = float(np.nanmedian(both)) if np.isfinite(both).any() else float("nan")
    mad = _median_absolute_deviation(both)
    scale = MAD_TO_SIGMA * mad if np.isfinite(mad) else float("nan")
    denom = max(scale, epsilon) if np.isfinite(scale) else epsilon
    return (left - center) / denom, (right - center) / denom, center, (scale if np.isfinite(scale) else 0.0)


def robust_z(values: np.ndarray) -> np.ndarray:
    """Signed robust z-score against the (outlier-resistant) cohort null.

    When the cohort MAD is exactly 0 (>50% of controls share a value, common for
    bounded summaries such as ``fraction_zero``) we fall back to the standard
    deviation rather than dividing by ``EPSILON`` -- the latter would amplify a
    tiny deviation by ~1e8 and manufacture false flags. If the values are truly
    constant, all z-scores are 0.
    """
    v = np.asarray(values, dtype=float)
    med = np.nanmedian(v)
    mad = MAD_TO_SIGMA * _median_absolute_deviation(v)
    if np.isfinite(mad) and mad > 0:
        denom = mad
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            std = float(np.nanstd(v))
        if not (np.isfinite(std) and std > 0):
            return np.zeros_like(v)
        denom = std
    return (v - med) / denom


def empirical_upper_tail_pvalue(scores: np.ndarray) -> np.ndarray:
    """Leave-one-out empirical upper-tail p-value (kept for optional reporting)."""
    s = np.asarray(scores, dtype=float)
    out = np.full(s.shape, np.nan, dtype=float)
    finite = np.isfinite(s)
    n = int(finite.sum())
    if n == 0:
        return out
    fv = s[finite]
    for i in np.where(finite)[0]:
        out[i] = min(1.0, (1 + int(np.sum(fv >= s[i])) - 1) / n)
    return out


# ---------------------------------------------------------------------------
# GIFTI loading / validation / raw QC / masking
# ---------------------------------------------------------------------------
def load_gifti_metric(path, multi_array: str = "mean") -> np.ndarray:
    """Load a GIFTI metric as a 1-D float array (mean/first over multiple arrays)."""
    if nib is None:  # pragma: no cover
        raise RuntimeError("nibabel is required to read GIFTI files")
    img = nib.load(str(path))
    if not img.darrays:
        raise ValueError(f"{path} has no GIFTI data arrays")
    arrays = [np.asarray(d.data, dtype=np.float64).reshape(-1) for d in img.darrays]
    if len(arrays) == 1 or multi_array == "first":
        return arrays[0]
    sizes = {a.size for a in arrays}
    if len(sizes) != 1:
        raise ValueError(f"{path} has GIFTI arrays with different vertex counts: {sorted(sizes)}")
    return np.mean(np.vstack(arrays), axis=0)


def _structural_metrics(x: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    n = int(x.size)
    finite = np.isfinite(x)
    vals = x[finite]
    return {
        "n_vertices": n,
        "n_nan": int(np.isnan(x).sum()),
        "n_inf": int(np.isinf(x).sum()),
        "n_zero": int(np.sum(vals == 0.0)) if vals.size else 0,
        "fraction_invalid": float((n - vals.size) / n) if n else 1.0,
        "mad": _median_absolute_deviation(x),
    }


def _structural_flags(metrics: dict, max_invalid_fraction: float = 0.5) -> list:
    flags = []
    if metrics.get("n_vertices", 0) == 0:
        return ["invalid_file"]
    n_finite = metrics["n_vertices"] - metrics["n_nan"] - metrics["n_inf"]
    if n_finite > 0 and metrics["n_zero"] == n_finite:
        flags.append("all_zero")
    if n_finite > 0 and np.isfinite(metrics.get("mad", np.nan)) and metrics["mad"] == 0.0 and "all_zero" not in flags:
        flags.append("constant_map")
    if metrics.get("fraction_invalid", 0.0) > max_invalid_fraction:
        flags.append("missing_values")
    return flags


def calculate_raw_qc_metrics(x: np.ndarray, mask: Optional[np.ndarray] = None) -> dict:
    """Raw-map QC summaries (Step 3) on the ORIGINAL, unnormalized map."""
    x = np.asarray(x, dtype=float)
    if mask is not None:
        x = x[mask]
    finite = np.isfinite(x)
    vals = x[finite]
    n = int(x.size)
    if vals.size == 0:
        nan = float("nan")
        return {"median": nan, "mad": nan, "iqr": nan, "p01": nan, "p05": nan, "p95": nan,
                "p99": nan, "dynamic_range": nan, "fraction_zero": nan,
                "fraction_nonfinite": 1.0, "n_valid_vertices": 0}
    p01, p05, p25, p75, p95, p99 = np.percentile(vals, [1, 5, 25, 75, 95, 99])
    return {
        "median": float(np.median(vals)),
        "mad": _median_absolute_deviation(x),
        "iqr": float(p75 - p25),
        "p01": float(p01), "p05": float(p05), "p95": float(p95), "p99": float(p99),
        "dynamic_range": float(p99 - p01),
        "fraction_zero": float(np.sum(vals == 0.0) / n),
        "fraction_nonfinite": float((n - vals.size) / n),
        "n_valid_vertices": int(vals.size),
    }


def build_common_vertex_mask(maps: Sequence[np.ndarray], drop_low_variance: bool = True,
                             variance_eps: float = 1e-10) -> np.ndarray:
    """Vertices valid for EVERY control (finite everywhere; drop near-zero variance).

    The near-zero-variance drop removes the fsLR medial wall (constant/zero there
    in most feature maps) without a separate atlas file.
    """
    if not maps:
        return np.zeros(0, dtype=bool)
    stacked = np.vstack([np.asarray(m, dtype=float).reshape(-1) for m in maps])
    mask = np.all(np.isfinite(stacked), axis=0)
    if drop_low_variance:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            variances = np.nanvar(np.where(np.isfinite(stacked), stacked, np.nan), axis=0)
        mask = mask & (variances > variance_eps)
    return mask


# ---------------------------------------------------------------------------
# Covariate residualization (optional, cross-fitted)
# ---------------------------------------------------------------------------
def _build_design_matrix(frame: pd.DataFrame, covariates: Sequence[str]):
    n = len(frame)
    cols, names = [np.ones(n)], ["intercept"]
    for cov in covariates:
        if cov not in frame.columns:
            continue
        series = frame[cov]
        if pd.api.types.is_numeric_dtype(series) and series.nunique(dropna=True) > 2:
            values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
            std = np.nanstd(values)
            values = (values - np.nanmean(values)) / (std if std > 0 else 1.0)
            cols.append(np.nan_to_num(values)); names.append(cov)
        else:
            dummies = pd.get_dummies(series.astype("object"), prefix=cov, dummy_na=False, drop_first=True)
            for col in dummies.columns:
                cols.append(dummies[col].to_numpy(dtype=float)); names.append(str(col))
    return np.column_stack(cols), names


def residualize_covariates_cross_fitted(X, frame, covariates, *, n_splits=5,
                                        stratify_col=None, random_seed=42):
    """Cross-fitted covariate residualization (Step: optional). Returns ``(residuals, status)``."""
    usable = [c for c in covariates if c in frame.columns]
    if not usable:
        return X, "no_covariates"
    n = X.shape[0]
    if n < 2 * n_splits:
        return X, f"skipped_cross_fit_small_sample(n={n})"
    try:
        from sklearn.model_selection import KFold, StratifiedKFold
    except Exception:  # pragma: no cover
        return X, "skipped_no_sklearn"
    if stratify_col and stratify_col in frame.columns and frame[stratify_col].nunique(dropna=True) > 1:
        groups = frame[stratify_col].astype("object").fillna("NA").to_numpy()
        if pd.Series(groups).value_counts().min() >= n_splits:
            split_iter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_seed).split(X, groups)
        else:
            split_iter = KFold(n_splits=n_splits, shuffle=True, random_state=random_seed).split(X)
    else:
        split_iter = KFold(n_splits=n_splits, shuffle=True, random_state=random_seed).split(X)
    design, _ = _build_design_matrix(frame, usable)
    residuals = X.copy()
    for tr, te in split_iter:
        beta, _, _, _ = np.linalg.lstsq(design[tr], X[tr], rcond=None)
        residuals[te] = X[te] - design[te] @ beta
    return residuals, "cross_fitted"


# ---------------------------------------------------------------------------
# Robust PCA (ROBPCA-equivalent) and cross-fitted distances
# ---------------------------------------------------------------------------
def _robust_score_cov(scores: np.ndarray, random_seed: int = 42):
    """Robust location/covariance in score space. Returns ``(loc, cov, status)``."""
    scores = np.atleast_2d(scores)
    n, k = scores.shape
    if n > k + 1:
        try:
            from sklearn.covariance import MinCovDet
            mcd = MinCovDet(random_state=random_seed).fit(scores)
            return mcd.location_, mcd.covariance_ + np.eye(k) * 1e-9, "mincovdet"
        except Exception:
            pass
    try:
        from sklearn.covariance import LedoitWolf
        lw = LedoitWolf().fit(scores)
        return np.median(scores, axis=0), lw.covariance_ + np.eye(k) * 1e-9, "ledoitwolf"
    except Exception:
        cov = np.cov(scores, rowvar=False)
        cov = np.atleast_2d(cov) + np.eye(k) * 1e-6
        return np.median(scores, axis=0), cov, "empirical"


@dataclass
class RobustPCA:
    center: np.ndarray
    components: np.ndarray  # k x v
    score_loc: np.ndarray
    score_inv_cov: np.ndarray
    n_components: int
    status: str

    def distances(self, X: np.ndarray):
        X = np.atleast_2d(np.asarray(X, dtype=float))
        centered = X - self.center
        scores = centered @ self.components.T
        d = scores - self.score_loc
        sd = np.sqrt(np.clip(np.einsum("ij,jk,ik->i", d, self.score_inv_cov, d), 0, None))
        recon = scores @ self.components + self.center
        od = np.sqrt(np.mean((X - recon) ** 2, axis=1))
        return od, sd


def fit_robust_pca(X_fit, *, variance_explained=0.90, min_components=2, max_components=50,
                   random_seed=42, reweight=True) -> RobustPCA:
    """Fit an ROBPCA-equivalent model: median-centred SVD PCA + one reweighting
    pass + robust score covariance. (A native ROBPCA is preferred when available.)"""
    def _svd_pca(X, center):
        Xc = X - center
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        return Xc, S, Vt

    center = np.median(X_fit, axis=0)
    Xc, S, Vt = _svd_pca(X_fit, center)
    evr = (S ** 2) / np.sum(S ** 2) if np.sum(S ** 2) > 0 else np.ones_like(S) / max(len(S), 1)
    cum = np.cumsum(evr)
    hard_max = max(1, min(int(max_components), X_fit.shape[0] - 1, X_fit.shape[1]))
    k = int(np.searchsorted(cum, float(variance_explained)) + 1)
    k = min(max(int(min_components), k), hard_max) if hard_max >= min_components else hard_max
    comps = Vt[:k]

    if reweight and X_fit.shape[0] > k + 2:
        scores = (X_fit - center) @ comps.T
        loc_s, cov_s, _ = _robust_score_cov(scores, random_seed)
        d = scores - loc_s
        sd = np.sqrt(np.clip(np.einsum("ij,jk,ik->i", d, np.linalg.pinv(cov_s), d), 0, None))
        recon = scores @ comps + center
        od = np.sqrt(np.mean((X_fit - recon) ** 2, axis=1))
        inliers = (robust_z(sd) < 3.0) & (robust_z(od) < 3.0)
        if inliers.sum() >= max(k + 2, int(0.5 * X_fit.shape[0])) and inliers.sum() < X_fit.shape[0]:
            center = np.median(X_fit[inliers], axis=0)
            _, _, Vt2 = _svd_pca(X_fit[inliers], center)
            k = min(k, Vt2.shape[0])
            comps = Vt2[:k]

    scores = (X_fit - center) @ comps.T
    loc_s, cov_s, cov_status = _robust_score_cov(scores, random_seed)
    return RobustPCA(center=center, components=comps, score_loc=loc_s,
                     score_inv_cov=np.linalg.pinv(cov_s), n_components=k,
                     status=f"robpca_equiv(k={k},cov={cov_status})")


def cross_fitted_bilateral_distances(X, *, variance_explained=0.90, max_components=50,
                                     n_splits=5, minimum_controls=15, random_seed=42):
    """Cross-fitted (K-fold / LOO) orthogonal + score distances on the bilateral
    QC matrix. Returns ``(od, sd, status)`` or ``(None, None, reason)`` if skipped."""
    X = np.asarray(X, dtype=float)
    n = X.shape[0]
    if n < minimum_controls:
        return None, None, f"skipped_pca_few_controls(n={n}<{minimum_controls})"
    from sklearn.model_selection import KFold
    splits = n if n <= 30 else n_splits  # leave-one-out for small cohorts
    folds = list(KFold(n_splits=splits, shuffle=True, random_state=random_seed).split(np.arange(n)))
    od = np.full(n, np.nan)
    sd = np.full(n, np.nan)
    statuses = set()
    for tr, te in folds:
        model = fit_robust_pca(X[tr], variance_explained=variance_explained,
                               max_components=max_components, random_seed=random_seed)
        o, s = model.distances(X[te])
        od[te], sd[te] = o, s
        statuses.add(model.status)
    return od, sd, ";".join(sorted(statuses))


# ---------------------------------------------------------------------------
# Combined anomaly statistic + calibration
# ---------------------------------------------------------------------------
def combined_anomaly_statistic(components: np.ndarray) -> np.ndarray:
    """``A_i = max_j Z_ij`` over finite standardized components (upper tail)."""
    C = np.atleast_2d(np.asarray(components, dtype=float))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmax(C, axis=1)


def calibrate_anomaly_threshold(components: np.ndarray, qc_alpha: float, *,
                                n_boot: int = 20000, random_seed: int = 42,
                                calibration: str = "bootstrap"):
    """Calibrate the max-statistic once. Returns ``(A, threshold, pvalue, status)``.

    ``bootstrap`` simulates the component vector from its robust joint (MVN)
    inlier distribution and takes the max, so the correlation among the maximized
    components is accounted for automatically. ``analytic`` uses a Bonferroni
    normal cutoff over the components (ignores correlation; conservative).
    """
    C = np.atleast_2d(np.asarray(components, dtype=float))
    complete = np.where(np.all(np.isfinite(C), axis=0))[0]
    c = complete.size
    if c == 0:
        return (combined_anomaly_statistic(C), float("inf"),
                np.full(C.shape[0], np.nan), "no_calibration")
    # A must maximize over exactly the components the null is built from, so use
    # the complete columns for both (a per-subject NaN in a component otherwise
    # drops it from the null while A still counts it, deflating the threshold).
    A = combined_anomaly_statistic(C[:, complete])

    if calibration == "analytic":
        from scipy.stats import norm
        thr = float(norm.ppf(1.0 - qc_alpha / c))
        pvalue = np.array([min(1.0, c * float(norm.sf(a))) if np.isfinite(a) else np.nan for a in A])
        return A, thr, pvalue, f"analytic_bonferroni(c={c})"

    Z = C[:, complete]
    loc, cov, cov_status = _robust_score_cov(Z, random_seed)
    # The components are already robust-z standardized (median 0, MAD ~1 for
    # inliers), so the inlier null has unit variance. Simulate from the robust
    # CORRELATION with loc=0 and a unit diagonal: this keeps the max
    # correlation-aware while preventing a contaminated MinCovDet covariance from
    # inflating the diagonal -- and thus the threshold -- at high outlier
    # fractions, which would otherwise silently pass every real anomaly.
    cov = np.atleast_2d(cov)
    d = np.sqrt(np.clip(np.diag(cov), EPSILON, None))
    corr = np.clip(cov / np.outer(d, d), -0.999, 0.999)
    np.fill_diagonal(corr, 1.0)
    rng = np.random.default_rng(random_seed)
    sims = rng.multivariate_normal(np.zeros(c), corr, size=int(n_boot))
    a_null = sims.max(axis=1)
    thr = float(np.quantile(a_null, 1.0 - qc_alpha))
    pvalue = np.array([(1 + int(np.sum(a_null >= a))) / (a_null.size + 1) if np.isfinite(a) else np.nan
                       for a in A])
    return A, thr, pvalue, f"bootstrap(n={int(n_boot)},cov={cov_status})"


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def _normalize_hemisphere(value) -> Optional[str]:
    s = str(value).strip().lower()
    if s in {"l", "lh", "left", "hemi-l"} or s.endswith("l"):
        return "L"
    if s in {"r", "rh", "right", "hemi-r"} or s.endswith("r"):
        return "R"
    return None


@dataclass
class SurfaceQCResult:
    table: pd.DataFrame
    summary: dict
    exclusions: pd.DataFrame
    output_dir: Optional[Path] = None
    figures: list = field(default_factory=list)


def _software_versions() -> dict:
    v = {"numpy": np.__version__, "pandas": pd.__version__}
    try:
        import sklearn
        v["scikit-learn"] = sklearn.__version__
    except Exception:
        pass
    if nib is not None:
        v["nibabel"] = nib.__version__
    return v


def run_surface_qc(
    manifest: pd.DataFrame,
    output_dir=None,
    *,
    qc_alpha: float = 0.01,
    covariates: Sequence[str] = (),
    cross_fit_covariates: bool = False,
    variance_explained: float = 0.90,
    max_components: int = 50,
    minimum_controls: int = 15,
    group_by: Sequence[str] = ("site", "scanner", "protocol"),
    multi_array: str = "mean",
    drop_low_variance: bool = True,
    calibration: str = "bootstrap",
    n_boot: int = 20000,
    random_seed: int = 42,
    make_figures: bool = True,
    write_exclusion_list: bool = False,
) -> SurfaceQCResult:
    """Run bilateral pre-normalization surface QC over a control manifest.

    ``manifest`` needs ``subject, feature, hemisphere, file_path`` (plus optional
    ``site, scanner, protocol, age, sex``). One row per ``subject x feature`` is
    produced with the combined anomaly statistic and a flag at level ``qc_alpha``.
    Input files are never modified.
    """
    required = {"subject", "feature", "hemisphere", "file_path"}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f"manifest missing required columns: {sorted(missing)}")

    man = manifest.copy()
    man["hemisphere"] = man["hemisphere"].map(_normalize_hemisphere)
    optional_cols = [c for c in ("site", "scanner", "protocol", "age", "sex") if c in man.columns]
    warnings_log: list = []
    n_controls = man["subject"].nunique()
    active_group_cols = [c for c in group_by if c in man.columns and man[c].notna().any()]

    rows: list = []
    per_feature_diag: dict = {}
    features = sorted(man["feature"].unique())

    for feature in features:
        fman = man[man["feature"] == feature]
        subjects = sorted(fman["subject"].unique())
        left_files, right_files, covrows = {}, {}, {}
        for subject in subjects:
            srows = fman[fman["subject"] == subject]
            lrow = srows[srows["hemisphere"] == "L"]
            rrow = srows[srows["hemisphere"] == "R"]
            if len(lrow):
                left_files[subject] = str(lrow.iloc[0]["file_path"])
            if len(rrow):
                right_files[subject] = str(rrow.iloc[0]["file_path"])
            covrows[subject] = srows.iloc[0]
        paired = [s for s in subjects if s in left_files and s in right_files]
        for s in subjects:  # surface single-hemisphere inputs instead of dropping silently
            if (s in left_files) != (s in right_files):
                warnings_log.append(
                    f"{s}/{feature}: only one hemisphere present; skipped (QC needs both L and R)")

        # ---- validate + load ----
        loaded = {"L": {}, "R": {}}
        struct_flags = {}
        for subject in paired:
            flags = []
            for hemi, files in (("L", left_files), ("R", right_files)):
                try:
                    arr = load_gifti_metric(files[subject], multi_array=multi_array)
                    loaded[hemi][subject] = arr
                    flags.extend(_structural_flags(_structural_metrics(arr)))
                except Exception as exc:
                    flags.append("invalid_file")
                    warnings_log.append(f"{subject}/{feature}/{hemi}: {exc}")
            struct_flags[subject] = sorted(set(flags))

        # consistent vertex counts per hemisphere
        usable = [s for s in paired if s in loaded["L"] and s in loaded["R"]
                  and "invalid_file" not in struct_flags.get(s, [])]
        for hemi in ("L", "R"):
            sizes = {s: loaded[hemi][s].size for s in usable}
            if sizes:
                modal = pd.Series(list(sizes.values())).mode().iloc[0]
                for s, sz in sizes.items():
                    if sz != modal:
                        struct_flags[s] = sorted(set(struct_flags.get(s, []) + ["unexpected_vertex_count"]))
        usable = [s for s in usable if "unexpected_vertex_count" not in struct_flags.get(s, [])]

        mask_L = build_common_vertex_mask([loaded["L"][s] for s in usable], drop_low_variance) if usable else np.zeros(0, bool)
        mask_R = build_common_vertex_mask([loaded["R"][s] for s in usable], drop_low_variance) if usable else np.zeros(0, bool)

        # ---- raw QC + bilateral joint scaling ----
        raw, bilateral_qc = {}, {}
        constant = set()
        for s in usable:
            rl = calculate_raw_qc_metrics(loaded["L"][s], mask_L)
            rr = calculate_raw_qc_metrics(loaded["R"][s], mask_R)
            raw[s] = {"L": rl, "R": rr,
                      "lr_median_diff": abs(rl["median"] - rr["median"]),
                      "lr_mad_diff": abs(rl["mad"] - rr["mad"])}
            if mask_L.any() and mask_R.any():
                zl, zr, _, scale = robust_bilateral_scale(loaded["L"][s][mask_L], loaded["R"][s][mask_R])
                bilateral_qc[s] = np.concatenate([zl, zr])
                if scale <= EPSILON:
                    constant.add(s)
                    struct_flags[s] = sorted(set(struct_flags.get(s, []) + ["constant_map"]))
            else:
                constant.add(s)
        pca_subjects = [s for s in usable if s not in constant]

        # ---- cross-fitted robust PCA on the bilateral map ----
        od = sd = None
        status = "n/a"
        if len(pca_subjects) >= minimum_controls and mask_L.any() and mask_R.any():
            matrix = np.vstack([bilateral_qc[s] for s in pca_subjects])
            res_status = "none"
            if covariates and cross_fit_covariates:
                fframe = pd.DataFrame([covrows[s] for s in pca_subjects]).reset_index(drop=True)
                matrix, res_status = residualize_covariates_cross_fitted(
                    matrix, fframe, covariates,
                    stratify_col=("site" if "site" in fframe.columns else None),
                    random_seed=random_seed)
            od, sd, pca_status = cross_fitted_bilateral_distances(
                matrix, variance_explained=variance_explained, max_components=max_components,
                minimum_controls=minimum_controls, random_seed=random_seed)
            status = f"residual={res_status};{pca_status}"
        else:
            status = f"skipped_pca_few_controls(n={len(pca_subjects)}<{minimum_controls})"

        # ---- raw abnormality (Step 3 -> one raw component) ----
        raw_abn = {}
        raw_group_status = "feature_only"
        if usable:
            gframe = pd.DataFrame([covrows[s] for s in usable]).reset_index(drop=True)
            grp = np.array(["ALL"] * len(usable))
            if active_group_cols:
                candidate = gframe[active_group_cols].fillna("NA").apply(
                    lambda r: "|".join(map(str, r)), axis=1).to_numpy()
                if pd.Series(candidate).value_counts().min() >= minimum_controls:
                    grp = candidate
                    raw_group_status = "acquisition_group"
            summaries = [
                ("median_L", lambda r: r["L"]["median"]), ("median_R", lambda r: r["R"]["median"]),
                ("mad_L", lambda r: r["L"]["mad"]), ("mad_R", lambda r: r["R"]["mad"]),
                ("dr_L", lambda r: r["L"]["dynamic_range"]), ("dr_R", lambda r: r["R"]["dynamic_range"]),
                ("fz_L", lambda r: r["L"]["fraction_zero"]), ("fz_R", lambda r: r["R"]["fraction_zero"]),
                ("lr_median_diff", lambda r: r["lr_median_diff"]), ("lr_mad_diff", lambda r: r["lr_mad_diff"]),
            ]
            max_abn = np.zeros(len(usable))
            for _name, getter in summaries:
                y = np.array([getter(raw[s]) for s in usable], dtype=float)
                rz = np.full(len(usable), np.nan)
                for g in np.unique(grp):
                    sel = grp == g
                    if sel.sum() >= 3:
                        rz[sel] = np.abs(robust_z(y[sel]))
                max_abn = np.nanmax(np.vstack([max_abn, np.nan_to_num(rz)]), axis=0)
            raw_abn = {s: float(max_abn[i]) for i, s in enumerate(usable)}

        # ---- standardize components + combined anomaly A on pca_subjects ----
        A_by_subject, pval_by_subject = {}, {}
        threshold = float("nan")
        z_od = z_sd = z_raw = {}
        if pca_subjects and od is not None:
            # Gaussianize the (right-skewed, non-negative) distances before robust
            # standardization so the MVN max-calibration is well-behaved: the ROBPCA
            # Wilson-Hilferty transform OD^(2/3), and SD^2 ~ chi^2 -> SD^(2/3).
            z_od_vec = robust_z(np.power(np.clip(od, 0.0, None), 2.0 / 3.0))
            z_sd_vec = robust_z(np.power(np.clip(sd, 0.0, None), 2.0 / 3.0))
            z_raw_vec = robust_z(np.sqrt(np.clip(
                np.array([raw_abn.get(s, np.nan) for s in pca_subjects]), 0.0, None)))
            z_od = {s: float(z_od_vec[i]) for i, s in enumerate(pca_subjects)}
            z_sd = {s: float(z_sd_vec[i]) for i, s in enumerate(pca_subjects)}
            z_raw = {s: float(z_raw_vec[i]) for i, s in enumerate(pca_subjects)}
            comps = np.column_stack([z_od_vec, z_sd_vec, z_raw_vec])
            A, threshold, pval, calib_status = calibrate_anomaly_threshold(
                comps, qc_alpha, n_boot=n_boot, random_seed=random_seed, calibration=calibration)
            status = f"{status};{calib_status}"
            A_by_subject = {s: float(A[i]) for i, s in enumerate(pca_subjects)}
            pval_by_subject = {s: float(pval[i]) for i, s in enumerate(pca_subjects)}

        # ---- emit rows ----
        for s in paired:
            sflags = list(struct_flags.get(s, []))
            structural_fail = bool({"invalid_file", "unexpected_vertex_count", "all_zero",
                                    "constant_map", "missing_values"} & set(sflags))
            rl = raw.get(s, {}).get("L", {})
            rr = raw.get(s, {}).get("R", {})
            A_val = A_by_subject.get(s, np.nan)
            flagged = bool(np.isfinite(A_val) and np.isfinite(threshold) and A_val > threshold) or structural_fail
            rows.append({
                "subject": s, "feature": feature,
                "left_file": left_files.get(s, ""), "right_file": right_files.get(s, ""),
                "site": covrows[s].get("site", "") if "site" in optional_cols else "",
                "scanner": covrows[s].get("scanner", "") if "scanner" in optional_cols else "",
                "protocol": covrows[s].get("protocol", "") if "protocol" in optional_cols else "",
                "n_valid_vertices_left": int(mask_L.sum()), "n_valid_vertices_right": int(mask_R.sum()),
                "raw_median_left": rl.get("median", np.nan), "raw_median_right": rr.get("median", np.nan),
                "raw_mad_left": rl.get("mad", np.nan), "raw_mad_right": rr.get("mad", np.nan),
                "raw_dynamic_range_left": rl.get("dynamic_range", np.nan),
                "raw_dynamic_range_right": rr.get("dynamic_range", np.nan),
                "fraction_zero_left": rl.get("fraction_zero", np.nan),
                "fraction_zero_right": rr.get("fraction_zero", np.nan),
                "left_right_median_difference": raw.get(s, {}).get("lr_median_diff", np.nan),
                "left_right_mad_difference": raw.get(s, {}).get("lr_mad_diff", np.nan),
                "orthogonal_distance_bilateral": float(od[pca_subjects.index(s)]) if (od is not None and s in pca_subjects) else np.nan,
                "score_distance_bilateral": float(sd[pca_subjects.index(s)]) if (sd is not None and s in pca_subjects) else np.nan,
                "z_od": z_od.get(s, np.nan), "z_sd": z_sd.get(s, np.nan), "z_raw": z_raw.get(s, np.nan),
                "raw_abnormality": raw_abn.get(s, np.nan),
                "anomaly_score": A_val,
                "anomaly_threshold": threshold,
                "anomaly_pvalue": pval_by_subject.get(s, np.nan),
                "flagged": flagged,
                "_structural_flags": ";".join(sflags),
                "model_status": status,
            })
        per_feature_diag[feature] = {
            "n_paired": len(paired), "n_usable": len(usable), "n_pca_subjects": len(pca_subjects),
            "n_valid_vertices_left": int(mask_L.sum()), "n_valid_vertices_right": int(mask_R.sum()),
            "anomaly_threshold": threshold, "status": status, "raw_group_status": raw_group_status,
        }

    table = pd.DataFrame(rows)
    if not table.empty:
        table["flag_reason"] = [_flag_reason(r) for _, r in table.iterrows()]
        table = table.drop(columns=["_structural_flags"], errors="ignore")
    else:
        for col in ("flag_reason",):
            table[col] = []
    if qc_alpha < 1.0 / max(n_controls, 1):
        warnings_log.append(
            f"qc_alpha={qc_alpha} is below 1/n_controls={1.0/max(n_controls,1):.4f}; "
            "the bootstrap threshold extrapolates into a sparsely sampled tail.")
    table["warnings"] = ";".join(sorted(set(warnings_log)))[:2000] if warnings_log else ""

    summary = {
        "n_controls": int(n_controls), "n_features": int(len(features)), "n_tests": int(len(table)),
        "qc_alpha": float(qc_alpha), "calibration": calibration,
        "n_flagged": int(table["flagged"].sum()) if "flagged" in table else 0,
        "pca_settings": {"variance_explained": variance_explained, "max_components": max_components,
                         "minimum_controls": minimum_controls},
        "covariates_used": list(covariates), "cross_fit_covariates": bool(cross_fit_covariates),
        "acquisition_group_columns": active_group_cols, "per_feature": per_feature_diag,
        "warnings": sorted(set(warnings_log)), "software_versions": _software_versions(),
        "empirical_pvalue_resolution": 1.0 / max(n_controls, 1),
    }
    exclusions = _build_exclusions(table)

    figures = []
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        generate_qc_tables(table, exclusions, summary, output_dir, write_exclusion_list=write_exclusion_list)
        if make_figures and not table.empty:
            figures = generate_qc_figures(table, output_dir)

    return SurfaceQCResult(table=table, summary=summary, exclusions=exclusions,
                           output_dir=Path(output_dir) if output_dir is not None else None,
                           figures=figures)


def _flag_reason(row: pd.Series) -> str:
    reasons = [t for t in str(row.get("_structural_flags", "") or "").split(";") if t]
    if bool(row.get("flagged", False)):
        comps = {"raw_intensity_outlier": row.get("z_raw", np.nan),
                 "bilateral_spatial_outlier": max(row.get("z_od", np.nan), row.get("z_sd", np.nan))
                 if np.isfinite(row.get("z_od", np.nan)) or np.isfinite(row.get("z_sd", np.nan)) else np.nan}
        finite = {k: v for k, v in comps.items() if np.isfinite(v)}
        if finite:
            best = max(finite.values())
            for k, v in finite.items():
                if v >= best - 1e-12:
                    reasons.append(k)
            reasons.append("combined_statistical_outlier")
    seen = []
    for r in reasons:
        if r and r not in seen:
            seen.append(r)
    return ";".join(seen)


def _build_exclusions(table: pd.DataFrame) -> pd.DataFrame:
    cols = ["subject", "feature", "hemisphere_action", "left_file", "right_file", "anomaly_pvalue", "flag_reason"]
    if table.empty or "flagged" not in table:
        return pd.DataFrame(columns=cols)
    flagged = table[table["flagged"]]
    return pd.DataFrame({
        "subject": flagged["subject"], "feature": flagged["feature"],
        "hemisphere_action": "remove_both",
        "left_file": flagged["left_file"], "right_file": flagged["right_file"],
        "anomaly_pvalue": flagged["anomaly_pvalue"], "flag_reason": flagged["flag_reason"],
    }).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Output tables and figures
# ---------------------------------------------------------------------------
def generate_qc_tables(table, exclusions, summary, output_dir, write_exclusion_list=False):
    output_dir = Path(output_dir)
    tsv = output_dir / "pre_normalization_surface_qc.tsv"
    table.to_csv(tsv, sep="\t", index=False)
    with open(output_dir / "pre_normalization_surface_qc_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, default=str)
    if write_exclusion_list:
        exclusions.to_csv(output_dir / "pre_normalization_surface_qc_exclusions.tsv", sep="\t", index=False)
    return tsv


def generate_qc_figures(table, output_dir) -> list:
    """Per-feature QC report figures. Returns the written paths."""
    import os
    os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib-zbrains"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    figures = []
    for feature, fsub in table.groupby("feature"):
        flagged = fsub["flagged"].fillna(False).to_numpy()
        thr = float(fsub["anomaly_threshold"].dropna().iloc[0]) if fsub["anomaly_threshold"].notna().any() else np.nan
        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        fig.suptitle(f"Pre-normalization surface QC (bilateral): {feature}", fontsize=14, fontweight="bold")

        ax = axes[0, 0]
        A = fsub["anomaly_score"].to_numpy(dtype=float)
        ax.hist(A[np.isfinite(A)], bins=20, color="#4f7d95", edgecolor="#26343d")
        if np.isfinite(thr):
            ax.axvline(thr, color="#c23b22", lw=1.6, label=f"threshold={thr:.2f}")
            ax.legend(frameon=False, fontsize=8)
        ax.set_title("anomaly score A = max(Z)"); ax.set_xlabel("A")

        ax = axes[0, 1]
        od = fsub["orthogonal_distance_bilateral"].to_numpy(dtype=float)
        sd = fsub["score_distance_bilateral"].to_numpy(dtype=float)
        ax.scatter(od[~flagged], sd[~flagged], s=18, color="#6b8fa3", label="ok")
        ax.scatter(od[flagged], sd[flagged], s=44, color="#c23b22", label="flagged")
        for i in np.where(flagged)[0]:
            ax.annotate(str(fsub.iloc[i]["subject"]), (od[i], sd[i]), fontsize=6)
        ax.set_title("bilateral OD vs SD (cross-fitted)"); ax.set_xlabel("orthogonal distance")
        ax.set_ylabel("robust score distance"); ax.legend(frameon=False, fontsize=8)

        ax = axes[0, 2]
        ax.hist(fsub["raw_median_left"].to_numpy(dtype=float), bins=15, alpha=0.6, label="median L", color="#4f7d95")
        ax.hist(fsub["raw_mad_left"].to_numpy(dtype=float), bins=15, alpha=0.6, label="MAD L", color="#c98a3b")
        ax.set_title("raw median / MAD (L)"); ax.legend(frameon=False, fontsize=8)

        ax = axes[1, 0]
        ml = fsub["raw_median_left"].to_numpy(dtype=float); mr = fsub["raw_median_right"].to_numpy(dtype=float)
        ax.scatter(ml[~flagged], mr[~flagged], s=18, color="#6b8fa3")
        ax.scatter(ml[flagged], mr[flagged], s=44, color="#c23b22")
        lims = [np.nanmin([ml, mr]), np.nanmax([ml, mr])]
        if np.all(np.isfinite(lims)):
            ax.plot(lims, lims, "0.6", lw=0.8)
        ax.set_title("raw median L vs R"); ax.set_xlabel("left"); ax.set_ylabel("right")

        ax = axes[1, 1]
        for col, color, lab in (("z_od", "#4f7d95", "Z_OD"), ("z_sd", "#c98a3b", "Z_SD"), ("z_raw", "#7a9e7e", "Z_raw")):
            vals = fsub[col].to_numpy(dtype=float)
            ax.hist(vals[np.isfinite(vals)], bins=15, alpha=0.55, label=lab, color=color)
        ax.set_title("component robust-Z"); ax.legend(frameon=False, fontsize=8)

        ax = axes[1, 2]; ax.set_axis_off()
        ranked = fsub.sort_values("anomaly_score", ascending=False).head(10)
        lines = [f"{r['subject']:>10}  A={r['anomaly_score']:.2f}  {'FLAG' if r['flagged'] else ''}"
                 for _, r in ranked.iterrows()]
        ax.text(0.0, 1.0, "Most abnormal maps\n" + "\n".join(lines), va="top", ha="left",
                family="monospace", fontsize=8)

        fig.tight_layout(rect=[0, 0, 1, 0.96])
        path = output_dir / f"qc_report_feature-{_safe(feature)}.png"
        fig.savefig(path, dpi=150, facecolor="white")
        plt.close(fig)
        figures.append(path)
        ranked.to_csv(output_dir / f"qc_ranked_feature-{_safe(feature)}.csv", index=False)
    return figures


def _safe(value) -> str:
    import re
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", str(value)).strip("_") or "value"


# ---------------------------------------------------------------------------
# Manifest discovery + CLI
# ---------------------------------------------------------------------------
def discover_control_manifest(dataset_root, *, subject_glob="sub-HC*", session="ses-01",
                              label="white", surf="fsLR-32k", smooth=None,
                              features=None, demographics=None) -> pd.DataFrame:
    """Build a QC manifest from ``<root>/<subject>/<session>/maps/cortex/*.func.gii``."""
    import re
    pattern = re.compile(
        r"(?P<subject>sub-[^_]+)_(?P<session>ses-[^_]+)_hemi-(?P<hemi>[LR])_"
        r"surf-(?P<surf>fsLR-[^_]+)_label-(?P<label>[^_]+)_feature-(?P<feature>.+?)_"
        r"smooth-(?P<smooth>[^_]+)\.func\.gii$")
    dataset_root = Path(dataset_root)
    rows = []
    for subject_dir in sorted(dataset_root.glob(subject_glob)):
        cortex = subject_dir / session / "maps" / "cortex"
        if not cortex.exists():
            continue
        for path in sorted(cortex.glob("*.func.gii")):
            m = pattern.match(path.name)
            if not m:
                continue
            g = m.groupdict()
            if g["surf"] != surf or g["label"] != label:
                continue
            if smooth and g["smooth"] != smooth:
                continue
            if features and g["feature"] not in set(features):
                continue
            rows.append({"subject": g["subject"], "feature": g["feature"],
                         "hemisphere": g["hemi"], "file_path": str(path)})
    manifest = pd.DataFrame(rows)
    if demographics is not None and not manifest.empty:
        manifest = manifest.merge(demographics, how="left", on="subject")
    return manifest


def _parse_args(argv=None):
    p = argparse.ArgumentParser(prog="zbrains-qc-surface",
                                description="Pre-normalization bilateral surface-map outlier detection.")
    p.add_argument("--manifest", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--qc-alpha", type=float, default=0.01)
    p.add_argument("--covariates", nargs="*", default=[])
    p.add_argument("--cross-fit-covariates", action="store_true")
    p.add_argument("--variance-explained", type=float, default=0.90)
    p.add_argument("--max-components", type=int, default=50)
    p.add_argument("--minimum-controls", type=int, default=15)
    p.add_argument("--group-by", default="site,scanner,protocol")
    p.add_argument("--calibration", choices=["bootstrap", "analytic"], default="bootstrap")
    p.add_argument("--n-boot", type=int, default=20000)
    p.add_argument("--flag-only", dest="write_exclusion_list", action="store_false", default=False)
    p.add_argument("--write-exclusion-list", dest="write_exclusion_list", action="store_true")
    p.add_argument("--no-figures", dest="make_figures", action="store_false", default=True)
    p.add_argument("--random-seed", type=int, default=42)
    return p.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    sep = "\t" if str(args.manifest).endswith((".tsv", ".txt")) else None
    manifest = pd.read_csv(args.manifest, sep=sep, engine="python")
    group_by = [c.strip() for c in str(args.group_by).split(",") if c.strip()]
    result = run_surface_qc(
        manifest, output_dir=args.output_dir, qc_alpha=args.qc_alpha,
        covariates=args.covariates, cross_fit_covariates=args.cross_fit_covariates,
        variance_explained=args.variance_explained, max_components=args.max_components,
        minimum_controls=args.minimum_controls, group_by=group_by,
        calibration=args.calibration, n_boot=args.n_boot, random_seed=args.random_seed,
        make_figures=args.make_figures, write_exclusion_list=args.write_exclusion_list)
    print(f"Wrote QC outputs to {args.output_dir}")
    print(f"  tests={result.summary['n_tests']}  flagged={result.summary['n_flagged']}  qc_alpha={args.qc_alpha}")
    for w in result.summary["warnings"][:10]:
        print(f"  WARNING: {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
