"""Cross-site feature harmonization for the pooled control cohort.

Harmonizes per-vertex surface features ACROSS sites (e.g. MICs vs NOEL) so a
pooled cross-site control cohort can serve as one normative reference. All
estimators are FIT ON CONTROLS ONLY and then applied, FROZEN, to patients (no
leakage; patient values never enter the fit). Biological covariates (age, sex)
are preserved; only the site (batch) location/scale is removed.

Four transforms (a cheap→rich staircase that decomposes the harmonization effect):
  meancenter   -- remove per-vertex site LOCATION only (anchored to the reference
                  site; ref unchanged). Fewest parameters, most stable at small n.
  sitestd_noeb -- ComBat location+scale WITHOUT empirical Bayes (per-vertex
                  moment matching). Isolates the value of the scale term.
  combat       -- parametric ComBat with empirical-Bayes shrinkage of the site
                  location/scale across the VERTEX dimension (Johnson et al. 2007,
                  Biostatistics 8:118-127; Fortin et al. 2017/2018 NeuroImage).
  mcombat_ref  -- reference-batch (M-ComBat) ComBat: the reference site is left
                  EXACTLY unchanged and the other site is shifted/scaled onto it
                  (Da-ano et al. 2020). Natural for two sites with unequal n.

Design (fit): D = [batch one-hot | covariates]  (batch indicators absorb the
intercept). B = lstsq(D, Y). Standardized residuals Z are computed against the
grand (or reference) mean + covariate effect and pooled variance; per-site
gamma_hat/delta_hat are the mean/var of Z within a site; combat shrinks them via
EB priors pooled across vertices. Adjustment removes gamma/delta and adds the
covariate effect back, so age/sex structure survives and site structure does not.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

EPS = 1e-8
# Below this many controls at a site, its per-vertex SCALE (delta) is too noisy to
# estimate; fall back to location-only for that site (delta forced to 1).
DEFAULT_MIN_SCALE_N = 15


@dataclass(frozen=True)
class HarmonizationParams:
    method: str
    sites: tuple                      # site labels, in the column order used at fit
    ref_site: object                  # reference site label (mcombat_ref/meancenter) or None
    covariate_count: int              # number of covariate columns expected by apply
    grand_mean: np.ndarray            # (V,) standardization mean (ref mean for mcombat_ref)
    beta: np.ndarray                  # (p_cov, V) frozen covariate coefficients
    sigma: np.ndarray                 # (V,) standardization SD = sqrt(pooled/ref variance)
    gamma_star: dict = field(default_factory=dict)   # site -> (V,) additive effect on Z
    delta_star: dict = field(default_factory=dict)   # site -> (V,) multiplicative variance on Z
    location_offset: dict = field(default_factory=dict)  # meancenter: site -> (V,) offset vs ref
    passthrough: np.ndarray = None    # (V,) bool: zero-variance vertices left unchanged


def _one_hot(site, sites):
    idx = {s: i for i, s in enumerate(sites)}
    oh = np.zeros((len(site), len(sites)), dtype=float)
    for r, s in enumerate(site):
        oh[r, idx[s]] = 1.0
    return oh


def _as_cov(cov, n):
    if cov is None:
        return np.zeros((n, 0), dtype=float)
    cov = np.asarray(cov, dtype=float)
    if cov.ndim == 1:
        cov = cov.reshape(n, -1)
    return cov


def fit_site_harmonization(Y, site, cov, method, ref_site="MICs",
                           min_scale_n=DEFAULT_MIN_SCALE_N, eps=EPS):
    """Fit a harmonizer on CONTROLS ONLY. ``Y`` (n_ctrl, V), ``site`` (n_ctrl,),
    ``cov`` (n_ctrl, p) or None. Returns frozen :class:`HarmonizationParams`."""
    Y = np.asarray(Y, dtype=float)
    if Y.ndim != 2:
        raise ValueError("Y must be (n_controls, n_vertices)")
    n, V = Y.shape
    site = np.asarray(site)
    cov = _as_cov(cov, n)
    p = cov.shape[1]
    sites = tuple(dict.fromkeys(site.tolist()))     # preserve first-seen order
    if method in ("mcombat_ref", "meancenter") and ref_site not in sites:
        # anchor to the largest site if the named reference is absent
        counts = {s: int(np.sum(site == s)) for s in sites}
        ref_site = max(counts, key=counts.get)
    ref_idx = sites.index(ref_site) if ref_site in sites else 0

    oh = _one_hot(site, sites)
    design = np.concatenate([oh, cov], axis=1)      # (n, B+p)
    B = np.linalg.lstsq(design, Y, rcond=None)[0]   # (B+p, V)
    batch_int = B[: len(sites)]                     # (B, V) per-site intercepts
    beta = B[len(sites):] if p else np.zeros((0, V))

    # zero-variance vertices (e.g. all-zero medial wall) are passed through
    col_var = np.var(Y, axis=0)
    passthrough = ~(np.isfinite(col_var) & (col_var > eps))

    # ---- meancenter: remove site location only (anchor to ref) --------------
    if method == "meancenter":
        location_offset = {s: (batch_int[i] - batch_int[ref_idx]) for i, s in enumerate(sites)}
        for s in location_offset:
            location_offset[s] = np.where(passthrough, 0.0, location_offset[s])
        return HarmonizationParams(
            method=method, sites=sites, ref_site=ref_site, covariate_count=p,
            grand_mean=np.zeros(V), beta=beta, sigma=np.ones(V),
            location_offset=location_offset, passthrough=passthrough)

    # ---- combat / sitestd_noeb / mcombat_ref: location + scale --------------
    if method == "mcombat_ref":
        grand_mean = batch_int[ref_idx].copy()                      # anchor mean = ref site
    else:
        fractions = np.array([np.mean(site == s) for s in sites])   # sample-size weighted
        grand_mean = fractions @ batch_int                          # (V,)
    stand_mean = grand_mean[None, :] + (cov @ beta if p else 0.0)   # (n, V)

    resid = Y - design @ B
    if method == "mcombat_ref":
        ref_rows = site == ref_site
        var_pooled = np.sum(resid[ref_rows] ** 2, axis=0) / max(int(ref_rows.sum()), 1)
    else:
        var_pooled = np.sum(resid ** 2, axis=0) / n
    sigma = np.sqrt(np.where(passthrough, 1.0, np.maximum(var_pooled, eps)))

    Z = (Y - stand_mean) / sigma[None, :]

    gamma_hat, delta_hat, site_n = {}, {}, {}
    for i, s in enumerate(sites):
        rows = site == s
        nb = int(rows.sum())
        site_n[s] = nb
        gamma_hat[s] = np.mean(Z[rows], axis=0) if nb else np.zeros(V)
        delta_hat[s] = np.var(Z[rows], axis=0, ddof=1) if nb > 1 else np.ones(V)

    gamma_star, delta_star = {}, {}
    for s in sites:
        if method == "mcombat_ref" and s == ref_site:
            gamma_star[s] = np.zeros(V)             # ref passes through unchanged
            delta_star[s] = np.ones(V)
            continue
        gh, dh, nb = gamma_hat[s], delta_hat[s], site_n[s]
        if method == "combat":
            # EB priors pooled ACROSS vertices for this site
            gbar = float(np.mean(gh))
            tau2 = float(np.var(gh, ddof=1)) if gh.size > 1 else 0.0
            m = float(np.mean(dh)); s2 = float(np.var(dh, ddof=1)) if dh.size > 1 else 0.0
            if s2 <= eps or tau2 <= eps:
                gamma_star[s] = gh                  # nothing to shrink toward
                delta_star[s] = dh
            else:
                lam = (2 * s2 + m ** 2) / s2
                theta = (m ** 3 + m * s2) / s2
                g = gh.copy(); d = dh.copy()
                for _ in range(100):
                    g_new = (nb * tau2 * gh + d * gbar) / (nb * tau2 + d)
                    ssq = np.sum((Z[site == s] - g_new[None, :]) ** 2, axis=0)
                    d_new = (theta + 0.5 * ssq) / (nb / 2.0 + lam - 1.0)
                    if np.max(np.abs(g_new - g)) < 1e-4 and np.max(np.abs(d_new - d)) < 1e-4:
                        g, d = g_new, d_new
                        break
                    g, d = g_new, d_new
                gamma_star[s], delta_star[s] = g, d
        else:  # sitestd_noeb: plain moment matching, no shrinkage
            gamma_star[s] = gh
            delta_star[s] = dh
        # small-n / degenerate scale guards -> location only for this site
        if nb < 2 or nb < min_scale_n:
            delta_star[s] = np.ones(V)
        delta_star[s] = np.where(np.isfinite(delta_star[s]) & (delta_star[s] > eps),
                                 delta_star[s], 1.0)
        gamma_star[s] = np.where(np.isfinite(gamma_star[s]), gamma_star[s], 0.0)

    return HarmonizationParams(
        method=method, sites=sites, ref_site=(ref_site if method == "mcombat_ref" else None),
        covariate_count=p, grand_mean=grand_mean, beta=beta, sigma=sigma,
        gamma_star=gamma_star, delta_star=delta_star, passthrough=passthrough)


def apply_site_harmonization(Y, site, cov, params: HarmonizationParams):
    """Apply FROZEN ``params`` to any subjects (controls or patients). Never
    recomputes site statistics. ``Y`` (m, V) or (V,); ``site`` (m,) or scalar;
    ``cov`` (m, p)/(p,)/None. Vertices flagged passthrough or subjects from an
    unknown site are returned unchanged. Non-finite inputs stay non-finite."""
    single = np.asarray(Y).ndim == 1
    Y = np.atleast_2d(np.asarray(Y, dtype=float))
    m, V = Y.shape
    site = np.atleast_1d(np.asarray(site))
    if site.size == 1 and m > 1:
        site = np.repeat(site, m)
    cov = _as_cov(cov, m)
    out = Y.copy()

    for r in range(m):
        s = site[r].item() if hasattr(site[r], "item") else site[r]
        if s not in params.sites:
            continue                                # unknown site -> unchanged
        if params.method == "meancenter":
            off = params.location_offset[s]
            out[r] = Y[r] - off
            continue
        beta_term = (cov[r] @ params.beta) if params.covariate_count else 0.0
        stand_mean = params.grand_mean + beta_term
        z = (Y[r] - stand_mean) / params.sigma
        g = params.gamma_star[s]
        d = np.sqrt(params.delta_star[s])
        adj = (z - g) / d * params.sigma + stand_mean
        out[r] = adj
    # zero-variance vertices are never modified
    if params.passthrough is not None:
        out[:, params.passthrough] = Y[:, params.passthrough]
    return out[0] if single else out


class SiteHarmonizer:
    """Cross-site harmonizer for the pooled control cohort, driven by a two-pass
    export/score protocol so each cohort's OWN analyze_dataset does the loading and
    (site-specific) transforms.

    Pass 1 (export): each cohort's analyze_dataset stashes, per feature-map key
    ``hk``, its control (matrix, subjects, demographics); the driver feeds those in
    via :meth:`add_site_features`. Pass 2 (score): analyze_dataset asks for the
    POOLED harmonized control reference (:meth:`pooled_reference`) and harmonizes
    each patient (:meth:`apply_patient`) using FROZEN params fit on controls only.

    ``method='none'`` is never constructed (the driver skips harmonization). With
    ``site_covariate=True`` no feature transform is applied (identity); the pooled
    reference is returned with an added binary 'site' column for the normative
    design matrix.
    """

    def __init__(self, method, cov_columns, ref_site="MICs", site_covariate=False,
                 min_scale_n=DEFAULT_MIN_SCALE_N):
        self.method = method
        self.cov_columns = list(cov_columns) if cov_columns else []
        self.ref_site = ref_site
        self.site_covariate = bool(site_covariate)
        self.min_scale_n = min_scale_n
        self._exports = {}      # hk -> {site: (matrix, subjects, demo_df)}
        self._params = {}       # hk -> HarmonizationParams (cached)
        self._site_code = {}    # site label -> int code (for the site_covariate column)

    def add_site_features(self, site, export_dict):
        """Store one cohort's exported control features (hk -> (matrix, subjects, demo_df))."""
        if site not in self._site_code:
            self._site_code[site] = len(self._site_code)
        for hk, bundle in export_dict.items():
            self._exports.setdefault(hk, {})[site] = bundle

    def covers(self, hk):
        """True iff >=2 sites contributed controls for this feature-map AND every
        site shares the same feature dimensionality (vertex/depth count). A
        single-site feature, or one whose vertex space differs across sites (e.g.
        hippocampal surfaces from different hippunfold versions), is NOT covered ->
        it keeps its local per-site reference instead of crashing on a shape
        mismatch."""
        exports = self._exports.get(hk, {})
        if len(exports) < 2:
            return False
        shapes = {tuple(np.asarray(m).shape[1:]) for (m, _, _) in exports.values()}
        return len(shapes) == 1

    def harmonizes(self, hk, site):
        """True iff this map is pooled-harmonizable (covers) AND ``site`` itself
        contributed controls to the pool. A site absent from the pool for ``hk``
        (e.g. its controls dropped below the export threshold, or a 3rd site never
        pooled) must NOT be scored against the harmonized reference nor passed
        through unharmonized -- the caller falls back to that site's LOCAL
        reference. Closes the silent-passthrough hole."""
        return self.covers(hk) and site in self._exports.get(hk, {})

    def _cov(self, demo_df, n):
        if not self.cov_columns:
            return np.zeros((n, 0))
        return np.asarray(demo_df[self.cov_columns].to_numpy(dtype=float)).reshape(n, -1)

    def _fit(self, hk):
        if hk in self._params:
            return self._params[hk]
        exports = self._exports[hk]
        mats, sites, covs = [], [], []
        for site, (matrix, subjects, demo_df) in exports.items():
            m = np.asarray(matrix, dtype=float)
            mats.append(m)
            sites.extend([site] * m.shape[0])
            covs.append(self._cov(demo_df, m.shape[0]))
        Y = np.concatenate(mats, axis=0)          # (N, V) or (N, V, d)
        site = np.asarray(sites)
        cov = np.concatenate(covs, axis=0)
        fitter = lambda M: fit_site_harmonization(  # noqa: E731
            M, site, cov, self.method, ref_site=self.ref_site, min_scale_n=self.min_scale_n)
        params = harmonize_multidepth(fitter, Y) if Y.ndim > 2 else fitter(Y)
        self._params[hk] = params
        return params

    def pooled_reference(self, hk):
        """Pooled control reference for scoring pass. Returns
        (harmonized_matrix, subjects, demographics_df) or None (single-site)."""
        if not self.covers(hk):
            return None
        exports = self._exports[hk]
        mats, subjects, demos, sites = [], [], [], []
        for site, (matrix, subs, demo_df) in exports.items():
            m = np.asarray(matrix, dtype=float)
            mats.append(m)
            subjects.extend(list(subs))
            demos.append(demo_df)
            sites.extend([site] * m.shape[0])
        Y = np.concatenate(mats, axis=0)
        site = np.asarray(sites)
        demo = pd.concat(demos, ignore_index=True)
        if self.site_covariate:
            demo = demo.copy()
            demo["site"] = [self._site_code[s] for s in site]
            return Y, subjects, demo               # identity: no feature transform
        params = self._fit(hk)
        cov = self._cov(demo, Y.shape[0])
        applier = lambda M: apply_site_harmonization(M, site, cov, params)  # noqa: E731
        harmonized = harmonize_multidepth(applier, Y) if Y.ndim > 2 else applier(Y)
        return harmonized, subjects, demo

    def apply_patient(self, hk, pat_data, site, pat_demo):
        """Harmonize one patient's feature map with FROZEN params (identity for the
        site_covariate arm). ``pat_demo`` is a row (Series/dict) with cov_columns."""
        if self.site_covariate or not self.covers(hk):
            return pat_data
        params = self._fit(hk)
        cov = (np.array([float(pat_demo[c]) for c in self.cov_columns])
               if self.cov_columns else np.zeros(0))
        pat = np.asarray(pat_data, dtype=float)
        if pat.ndim > 1:   # (V, depth)
            flat = pat.reshape(1, -1)
            out = apply_site_harmonization(flat, site, cov.reshape(1, -1), params)
            return out.reshape(pat.shape)
        return apply_site_harmonization(pat, site, cov, params)


def harmonize_multidepth(fit_or_apply, Y, *args, **kwargs):
    """Helper for 3D (n, V, depth) arrays: reshape to (n, V*depth), run the 2D
    op, reshape back. ``fit_or_apply`` is a callable taking a 2D Y as first arg."""
    Y = np.asarray(Y, dtype=float)
    if Y.ndim <= 2:
        return fit_or_apply(Y, *args, **kwargs)
    n, Vv, d = Y.shape
    flat = Y.reshape(n, Vv * d)
    res = fit_or_apply(flat, *args, **kwargs)
    if isinstance(res, np.ndarray):
        return res.reshape(Y.shape)
    return res
