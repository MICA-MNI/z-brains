"""Greedy staged (coordinate-ascent) optimizer for the z-brains pipeline.

Where ``zbrains_benchmark.run_benchmark`` is a one-factor-at-a-time (OFAT) sweep
that always varies a single axis away from a *fixed* baseline, this driver does
sequential greedy optimization: it starts from the simplest possible pipeline
and, one stage at a time, tries that stage's candidate values, keeps whichever
maximizes the objective, *locks it in*, and carries it forward into the next
stage. Each stage's winner conditions every later stage.

Objective (chosen by the user): the per-subject macro AUROC for lesion detection
from cortical score maps, pooled across the MICs and NOEL cohorts as a single
pipeline (see :mod:`results.greedy_auroc`).

Stage order (the user's requested order):

    0. SIMPLEST_BASELINE -- no intensity normalization (raw nativepro T1w/FLAIR,
       normalization="none"), simplest blur (mean over depths), white-surface
       intensity sampling, and a plain per-vertex z-score with no demographics.
    1. smoothing           -> cortical/hippocampal smoothing FWHM pair
    2. self_normalization  -> wscore_preprocessing  (per-subject spatial rescale)
    3. intensity_sampling  -> intensity_depth_model  (depth-profile -> one value)
    4. blur                -> blur_depth_model        (blur profile -> one value)
    5. normative_model     -> (method, wscore_covariate_model,
                               prediction_variance_percentile):
       zscore -> wscore+demographics -> ... -> Gaussian process kernels
       -> GP + demographics + prediction-uncertainty filtering.

Because normalization is FIXED at "none" for the whole search (the user replaced
WhiteStripe/RAVEL with analysis-stage self-normalization), the heavy image
processing is limited to the visited smoothing x normalization bases; every
candidate after its base exists is a cheap re-analysis of the shared processed maps.

This module is additive: it imports the pure helpers from ``zbrains_benchmark``
and does not modify ``run_benchmark`` or its callers.
"""

from __future__ import annotations

import glob
import hashlib
import functools
import inspect
import warnings
import json
import os
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

import zbrains_benchmark as zb


# ---------------------------------------------------------------------------
# Simplest baseline (Stage 0). Axis values come from zbrains_benchmark.AXES;
# method / prediction_variance_percentile are driver-managed (not axes).
# ---------------------------------------------------------------------------
SIMPLEST_BASELINE = {
    "normalization": "none",                 # raw nativepro T1w/FLAIR; self-normalize downstream
    "wscore_preprocessing": "none",          # no per-subject spatial rescale yet
    "wscore_covariate_model": "linear",      # unused while method="zscore"
    "wscore_distribution": "gaussian",
    "use_curvature_covariates": False,
    "predictive_wscore": False,
    "blur_depth_model": "boundary_gradient",  # perpendicular GM/WM profile slope (Antel 2003/Hong); default + implicit keep for the blur stage
    "intensity_depth_model": "raw",          # resolved from sample_surface by validate()
    "wscore_surface_smoothing_iterations": 0,
    "control_correlation_filter": False,     # labels the correlation-exclusion arm
    # driver-managed (not zbrains_benchmark axes):
    "control_correlation_threshold": zb.DEFAULT_CONTROL_CORRELATION_THRESHOLD,
    # RANK-based pre-dataset-norm exclusion: for each feature, remove that feature
    # from the bottom-correlated fraction of controls. None = no exclusion.
    "control_correlation_quantile": None,    # e.g. 0.05 | 0.10 | 0.20
    "method": "zscore",                      # plain per-vertex z-score, no demographics
    "prediction_variance_percentile": None,  # GP uncertainty filter (None = off)
    # driver-managed feature-engineering, resolved by validate(). sample_surface
    # applies to ALL 5 intensity maps (T1w/FLAIR via intensity_depth_model now;
    # ADC/FA/qT1 sampling backend lands in a later step). Self-norm is T1w/FLAIR only.
    "sample_surface": "white",               # raw | midthickness | white | swm1 (all 5)
    "t1w_flair_self_norm": None,             # swm2 | owncortex (T1w/FLAIR only)
    "quant_sample_surface": None,            # derived: ADC/FA/qT1 plain-sample surface
    # driver-managed normalization composition, resolved into `normalization`:
    "subject_norm": "none",                  # none | whitestripe | wmmean
    "dataset_norm": "none",                  # none | ravel | nyul
    # driver-managed outlier detection. Correlation and ROBPCA both EXCLUDE
    # controls before dataset normalization.
    "outlier_method": "none",                # none | correlation | robpca
    "qc_alpha": None,                        # ROBPCA tail level (0.01 | 0.05 | 0.10)
    # driver-managed processing smoothing. The optimizer starts unsmoothed; the
    # smoothing stages can still opt into positive kernels when they win.
    "cortical_smoothing": 0,
    "hippocampal_smoothing": 0,
    # driver-managed CROSS-SITE feature harmonization (analysis-level; pools both
    # cohorts' controls into one harmonized reference). "none" = per-site
    # references (today's behavior). scoring_site_covariate adds a site dummy to
    # the normative design instead of transforming features (a scoring-stage arm).
    "site_harmonization": "none",            # none | meancenter | sitestd_noeb | combat | mcombat_ref
    "scoring_site_covariate": False,
}

# Keys the driver manages that are NOT zbrains_benchmark axes.
_NON_AXIS_KEYS = (
    "control_correlation_threshold", "control_correlation_quantile",
    "method", "prediction_variance_percentile",
    "sample_surface", "t1w_flair_self_norm", "quant_sample_surface",
    "subject_norm", "dataset_norm",
    "outlier_method", "qc_alpha",
    "cortical_smoothing", "hippocampal_smoothing",
    "site_harmonization", "scoring_site_covariate",
)

# Cross-site harmonization arms (analysis-level feature transforms; fit on pooled
# controls, applied frozen to patients). See src/zbrains/harmonization.py.
HARMONIZATION_ARMS = ("meancenter", "sitestd_noeb", "combat", "mcombat_ref")

# Smoothing search: cortical and hippocampal smoothing are optimized as SEPARATE
# greedy stages (see DEFAULT_STAGES). 0 == None (no smoothing -> the unsmoothed
# sampled metric; processing._smooth_or_copy_metric copies instead of running a
# zero-FWHM kernel). The baseline is 0mm for both compartments and is each
# stage's implicit "keep".
CORTICAL_SMOOTHING_OPTIONS = (0, 2, 5, 10)      # None, 2mm, 5mm, 10mm
HIPPOCAMPAL_SMOOTHING_OPTIONS = (0, 1, 2, 5)    # None, 1mm, 2mm, 5mm

GPR_KERNEL_OPTIONS = {
    "gp_rbf_ard": "gaussian_process",
    "gp_matern32_ard": "gaussian_process_matern32",
    "gp_matern52_ard": "gaussian_process_matern52",
    "gp_rational_quadratic_ard": "gaussian_process_rational_quadratic",
}

GPR_UNCERTAINTY_PERCENTILES = (95.0, 90.0, 80.0)

# Minimum controls for ROBPCA to fit (below this, run_surface_qc uses structural
# flags only). The heavy smoke lowers this for tiny subsets.
SURFACE_QC_MIN_CONTROLS = 15

# Composition of subject-level x dataset-level normalization -> the single
# `normalization` axis value. All 9 (subject_norm, dataset_norm) pairs are now
# buildable by the composable-normalization backend; the canonical desc- label
# for each pair comes straight from the backend so the two never drift.
from zbrains.normalization import compose_normalization_label as _compose_norm


def _resolve_normalization(config):
    """Map (subject_norm, dataset_norm) driver keys onto the ``normalization``
    axis (mutates ``config``) using the backend's canonical composite label."""
    subject_norm = config.get("subject_norm", "none")
    dataset_norm = config.get("dataset_norm", "none")
    config["normalization"] = _compose_norm(subject_norm, dataset_norm)

# T1w/FLAIR feature-engineering resolver tables (driver keys -> intensity_depth_model).
_T1WFLAIR_SURFACE_MODEL = {
    "midthickness": "sample_midthickness", "white": "sample_white", "swm1": "sample_swm1",
}
_T1WFLAIR_SELF_NORM_SUFFIX = {None: "", "none": "", "swm2": "_swm2", "owncortex": "_owncortex"}


def _resolve_intensity_depth_model(config):
    """Resolve the ONE ``sample_surface`` (all 5 maps) + ``t1w_flair_self_norm``
    driver keys into the two analysis params (mutates ``config``):
    ``intensity_depth_model`` (T1w/FLAIR: sampling + optional self-norm) and
    ``quant_sample_surface`` (ADC/FA/qT1: plain sample only). Absent/None/"raw"
    surface leaves raw maps (intensity_depth_model="raw", quant_sample_surface=None)."""
    surface = config.get("sample_surface")
    if not surface or surface in {"none", "raw"}:
        # Raw nativepro maps: reset BOTH derived params so 'raw' is genuinely raw
        # (previously left the baseline's intensity_depth_model in place, so a 'raw'
        # arm silently stayed white-sampled and cache-collided with 'white').
        config["intensity_depth_model"] = "raw"
        config["quant_sample_surface"] = None
        return
    if surface not in _T1WFLAIR_SURFACE_MODEL:
        raise ValueError(f"unknown sample_surface {surface!r}")
    self_norm = config.get("t1w_flair_self_norm")
    if self_norm not in _T1WFLAIR_SELF_NORM_SUFFIX:
        raise ValueError(f"unknown t1w_flair_self_norm {self_norm!r}")
    config["intensity_depth_model"] = (
        _T1WFLAIR_SURFACE_MODEL[surface] + _T1WFLAIR_SELF_NORM_SUFFIX[self_norm]
    )
    config["quant_sample_surface"] = surface  # ADC/FA/qT1: plain sample at this surface


def _validate_smoothing(config):
    # Cortical and hippocampal smoothing are now INDEPENDENT (separate greedy stages),
    # so validate each against its own option set rather than a fixed pair list.
    ctx = config.get("cortical_smoothing")
    hip = config.get("hippocampal_smoothing")
    if ctx not in CORTICAL_SMOOTHING_OPTIONS:
        raise ValueError(
            f"cortical_smoothing={ctx!r} not allowed; must be one of "
            f"{CORTICAL_SMOOTHING_OPTIONS} (0 = None)"
        )
    if hip not in HIPPOCAMPAL_SMOOTHING_OPTIONS:
        raise ValueError(
            f"hippocampal_smoothing={hip!r} not allowed; must be one of "
            f"{HIPPOCAMPAL_SMOOTHING_OPTIONS} (0 = None)"
        )


def _smoothing_label(config):
    ctx = config.get("cortical_smoothing", SIMPLEST_BASELINE["cortical_smoothing"])
    hip = config.get("hippocampal_smoothing", SIMPLEST_BASELINE["hippocampal_smoothing"])
    return f"smoothctx{ctx}hip{hip}"


# Smoothing-kernel CONVENTION tag, ALWAYS present in the processing signature so a
# processed base built under a different convention can never satisfy the
# completeness check. Bumped from implicit-sigma to explicit FWHM when the
# wb_command -metric-smoothing calls gained the `-fwhm` flag (processing.py): the
# on-disk filename token ('smooth-5mm') is unchanged, so WITHOUT this tag a post-fix
# 5mm-FWHM run would silently reuse pre-fix sigma-smoothed (~11.8mm FWHM) maps.
_SMOOTHING_CONVENTION = "smfwhm"


def _processing_signature(config, exclusion_signature=""):
    """Suffix processed/analysis dirs for processing-affecting driver settings.
    The smoothing-convention tag is ALWAYS included so stale (pre-`-fwhm`,
    sigma-smoothed) bases are never reused. The smoothing pair is always included
    too, including the 0/0 baseline, because older runs used the no-pair
    ``smfwhm`` suffix for the previous 5/2 baseline."""
    parts = [_SMOOTHING_CONVENTION, _smoothing_label(config)]
    if exclusion_signature:
        parts.append(exclusion_signature)
    return "_".join(parts)


def _base_signature(config, exclusion_signature=""):
    """Processing signature for the PROCESSED BASE (the per-subject maps).

    Identical to :func:`_processing_signature` EXCEPT the control-exclusion tag is
    dropped when ``dataset_norm == 'none'``: with no dataset-level fit (RAVEL/Nyul),
    every subject's whitestripe/sampled maps are byte-identical regardless of which
    controls are excluded, so all correlation/ROBPCA arms share ONE base instead of
    each re-running
    whitestripe + surface sampling. The exclusion still keys the ANALYSIS dir (whose
    z-reference DOES depend on the excluded set) and re-keys the base once a
    dataset-level fit is active (dataset_norm in {ravel, nyul})."""
    if config.get("dataset_norm", "none") == "none":
        exclusion_signature = ""
    return _processing_signature(config, exclusion_signature)


def _pipeline_settings(cohort, config, normalization=None):
    settings = dict(
        cohort.base_pipeline_settings,
        normalization=normalization or config["normalization"],
    )
    settings["cortical_smoothing"] = config["cortical_smoothing"]
    settings["hippocampal_smoothing"] = config["hippocampal_smoothing"]
    return settings


@dataclass(frozen=True)
class Stage:
    """One greedy stage: try each candidate override, lock the best.

    ``candidates`` is a list of ``(label, overrides)`` where ``overrides`` is a
    partial config applied on top of the current locked config. The current
    config is always re-entered as an implicit "keep" candidate so a stage can
    only ever improve (or tie) the running objective.

    ``condition`` is an optional ``callable(current_config) -> bool``. When it
    returns False the whole stage is SKIPPED (the locked config is carried
    forward unchanged). Used for stages that only make sense given an earlier
    lock -- e.g. the GP uncertainty filter is only swept when a GP kernel won the
    scoring stage, and then it is tuned on THAT winning kernel.
    """

    name: str
    candidates: list = field(default_factory=list)
    condition: object = None


def _default_stages():
    """The greedy stage order (coordinate ascent; each stage locks the arg-max of
    the balanced dataset x disease x feature AUROC, carried forward).

    Ordering follows the processing data flow (upstream image-level decisions
    before downstream surface/analysis ones, so each locked stage conditions the
    next): (1) subject-level norm -> (2) outlier detection [excludes controls] ->
    (3) dataset-level norm [composed on the locked subject norm] -> (4) smoothing
    [the last PROCESSING step -- a surface op applied after volume normalization +
    volume->surface mapping] -> (5) T1w/FLAIR sample surface -> (6) T1w/FLAIR
    self-norm -> (7) blur method (qT1/T1w/FLAIR) -> (8) site harmonization
    [cross-site, analysis-level] -> (9) scoring [normative model] -> (10) GP
    uncertainty filter [CONDITIONAL: only if a GP kernel won stage 9, tuned on
    that kernel]. Stages 1-4 rebuild the reusable base; 5-10 are cheap
    re-analysis. Each stage's implicit "keep" candidate is the current locked
    config, so a stage can only improve or tie.

    Correlation exclusions are feature-specific: a control scan can be omitted
    from one feature's dataset-normalization and normative fits while remaining
    available for every other feature. ROBPCA exclusions remain whole-control.
    """
    return [
        # (1) SUBJECT-LEVEL NORM (volume intensity; most upstream). implicit keep = none.
        Stage("subject_norm", [
            ("whitestripe", {"subject_norm": "whitestripe"}),
            ("wm_mean", {"subject_norm": "wmmean"}),
        ]),
        # (2) OUTLIER DETECTION (excludes flagged observations before dataset-norm).
        Stage("outlier_detection", [
            # Correlation arms: rank controls per feature on the full,
            # subject-normalized, UNSMOOTHED cohort and mask only the affected
            # subject-feature observations BEFORE dataset-norm. The same mask is
            # carried into normative scoring; other features remain available.
            ("corr_q05", {"outlier_method": "correlation", "control_correlation_filter": True,
                          "control_correlation_quantile": 0.05}),
            ("corr_q10", {"outlier_method": "correlation", "control_correlation_filter": True,
                          "control_correlation_quantile": 0.10}),
            ("corr_q20", {"outlier_method": "correlation", "control_correlation_filter": True,
                          "control_correlation_quantile": 0.20}),
            # ROBPCA arms: surface_qc flags controls and EXCLUDES them before
            # dataset-norm (removed from RAVEL/Nyul fit AND the normative reference).
            ("robpca_a0.01", {"outlier_method": "robpca", "qc_alpha": 0.01}),
            ("robpca_a0.05", {"outlier_method": "robpca", "qc_alpha": 0.05}),
            ("robpca_a0.10", {"outlier_method": "robpca", "qc_alpha": 0.10}),
        ]),
        # (3) DATASET-LEVEL NORM, composed on the locked subject norm. implicit keep = none.
        Stage("dataset_norm", [
            ("ravel", {"dataset_norm": "ravel"}),
            ("nyul", {"dataset_norm": "nyul"}),
        ]),
        # (4a) CORTICAL SMOOTHING -- surface smoothing of the cortical maps, applied
        # after the volume intensity normalizations + volume->surface mapping, so it
        # is optimized on the locked normalization. Separate from hippocampal
        # smoothing so each compartment gets its own kernel. 0 = None (unsmoothed).
        # implicit keep = 0mm (SIMPLEST_BASELINE cortical_smoothing).
        Stage("cortical_smoothing", [
            ("cortex_2mm", {"cortical_smoothing": 2}),
            ("cortex_5mm", {"cortical_smoothing": 5}),
            ("cortex_10mm", {"cortical_smoothing": 10}),
        ]),
        # (4b) HIPPOCAMPAL SMOOTHING -- surface smoothing of the hippocampal maps,
        # optimized independently on the locked cortical smoothing. 0 = None.
        # implicit keep = 0mm (SIMPLEST_BASELINE hippocampal_smoothing).
        Stage("hippocampal_smoothing", [
            ("hippocampus_1mm", {"hippocampal_smoothing": 1}),
            ("hippocampus_2mm", {"hippocampal_smoothing": 2}),
            ("hippocampus_5mm", {"hippocampal_smoothing": 5}),
        ]),
        # (5) SAMPLE SURFACE for all 5 intensity maps (T1w/FLAIR/qT1/ADC/FA).
        # implicit keep = white surface (the SIMPLEST_BASELINE default). The 'raw'
        # arm was removed: it did not actually bypass sampling (kept white) and
        # collided in the score cache with the white keep.
        Stage("sample_surface", [
            ("midthickness", {"sample_surface": "midthickness"}),
            ("swm1", {"sample_surface": "swm1"}),
        ]),
        # (6) T1w/FLAIR SELF-NORM on the locked surface. implicit keep = none.
        Stage("t1w_flair_self_norm", [
            ("swm2", {"t1w_flair_self_norm": "swm2"}),
            ("owncortex", {"t1w_flair_self_norm": "owncortex"}),
        ]),
        # (7) BLUR METHOD (*blur companions). Three literature-grounded metrics of
        # gray-white junction integrity, each a distinct construct, all sampling to
        # SWM2 (2mm; together they cover all 4 depths):
        #   boundary_gradient     = OLS intensity-profile slope over all 4 depths (Antel 2003/Hong)  [default/implicit keep]
        #   gray_white_contrast   = pctsurfcon %contrast, GM(mid) vs mean(SWM1,SWM2)
        #   juxtacortical_gradient= mean(SWM1,SWM2)-white sub-boundary extension (transmantle/MELD)
        Stage("blur", [
            ("gray_white_contrast", {"blur_depth_model": "gray_white_contrast"}),
            ("juxtacortical_gradient", {"blur_depth_model": "juxtacortical_gradient"}),
        ]),
        # (8) SITE HARMONIZATION -- cross-site feature harmonization of the POOLED
        # control cohort (analysis-level; fit on pooled controls, applied frozen to
        # patients; both sites scored against one harmonized reference). Checks
        # whether harmonizing helps vs per-site references. implicit keep = none.
        # Cheap->rich staircase decomposes the effect (pooling? additive? scale? EB?):
        Stage("site_harmonization", [
            ("meancenter", {"site_harmonization": "meancenter"}),      # location only
            ("sitestd_noeb", {"site_harmonization": "sitestd_noeb"}),  # location+scale, no EB
            ("combat", {"site_harmonization": "combat"}),              # ComBat EB (Johnson 2007/Fortin)
            ("mcombat_ref", {"site_harmonization": "mcombat_ref"}),    # ref-batch (anchor MICs)
        ]),
        # (9) SCORING. Explicitly score zscore in addition to the implicit keep,
        # so histories show it beside W-score/kNN/GPR. GP kernel families are
        # tested UNFILTERED here; the uncertainty filter is a SEPARATE conditional
        # stage (10) that runs only if a GP kernel wins -- so it is tuned on the
        # winning kernel, not hardcoded to ARD-RBF.
        # NOTE: the site_covariate arm (a site-dummy design-matrix alternative to
        # feature harmonization) is deferred -- it needs a consistent cross-cohort
        # 'site' column across every feature block; the plumbing/keys exist
        # (scoring_site_covariate) but it is not swept until that lands.
        Stage("scoring", [
            ("zscore", {"method": "zscore", "prediction_variance_percentile": None}),
            ("wscore_demographics", {"method": "wscore", "wscore_covariate_model": "linear"}),
            ("knn_half_controls", {"method": "wscore", "wscore_covariate_model": "knn"}),
            *[
                (
                    label,
                    {"method": "wscore", "wscore_covariate_model": model},
                )
                for label, model in GPR_KERNEL_OPTIONS.items()
            ],
        ]),
        # (10) GP UNCERTAINTY FILTER -- CONDITIONAL. Only runs when the scoring
        # stage locked a GP kernel; the candidates set ONLY the variance
        # percentile, so they inherit (are tuned on) the winning kernel. The
        # implicit keep is "no filter", so the stage can only improve or tie.
        # Skipped entirely when zscore/wscore/kNN won (uncertainty is a
        # GP-derived quantity and meaningless for those).
        Stage(
            "gp_uncertainty",
            [
                (f"uncertainty_p{int(pct)}", {"prediction_variance_percentile": pct})
                for pct in GPR_UNCERTAINTY_PERCENTILES
            ],
            condition=lambda c: (
                c.get("method") == "wscore"
                and c.get("wscore_covariate_model") in set(GPR_KERNEL_OPTIONS.values())
            ),
        ),
    ]


DEFAULT_STAGES = _default_stages()


# History rows written with this schema contain enough information to restart a
# stopped greedy run exactly at the last completed candidate/stage. Older CSVs
# are still useful as records, but are intentionally not used for automatic
# resume because they may predate defaults such as white-surface sampling or the
# explicit smoothctx0hip0 baseline signature.
HISTORY_SCHEMA_VERSION = 4


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------
def _axis_config(config):
    """Return just the zbrains_benchmark axis keys (drops driver-managed keys)."""
    return {k: v for k, v in config.items() if k not in _NON_AXIS_KEYS}


def validate(config):
    """Validate/normalize the axis portion of a full driver config.

    First resolves driver-managed feature-engineering keys (T1w/FLAIR sample
    surface + self-norm) into the ``intensity_depth_model`` axis, then validates
    the axis portion against ``zbrains_benchmark``.
    """
    config = dict(config)
    _validate_smoothing(config)
    _resolve_normalization(config)
    _resolve_intensity_depth_model(config)
    outlier_method = config.get("outlier_method", "none")
    if outlier_method not in {"none", "correlation", "robpca"}:
        raise ValueError(f"unknown outlier_method {outlier_method!r}")
    quantile = config.get("control_correlation_quantile")
    if outlier_method == "correlation":
        if not config.get("control_correlation_filter"):
            raise ValueError("correlation exclusion requires control_correlation_filter=True")
        if quantile is None or not 0.0 <= float(quantile) < 1.0:
            raise ValueError(
                "correlation exclusion requires control_correlation_quantile in [0, 1)"
            )
        config["control_correlation_quantile"] = float(quantile)
    axis_only = _axis_config(config)
    normalized, warnings = zb.validate_configuration(axis_only)
    merged = dict(config)
    merged.update(normalized)
    return merged, warnings


def staged_output_directory(config, output_dir_prefix, control_correlation_threshold,
                            exclusion_signature=""):
    """Deterministic analysis output dir, disambiguating the GP variance filter.

    ``output_directory_for`` encodes every axis (incl. normalization -> SELF, the
    normative covariate model, distribution, etc.) but knows nothing about the
    driver-managed ``prediction_variance_percentile``; two GP configs that differ
    only by the filter would otherwise collide, so append a suffix.
    ``exclusion_signature`` keeps distinct pre-dataset-norm exclusion arms apart.
    """
    quantile = config.get("control_correlation_quantile")
    if quantile is None:
        quantile = zb.DEFAULT_CONTROL_CORRELATION_QUANTILE
    base = zb.output_directory_for(
        _axis_config(config), output_dir_prefix, quantile,
        exclusion_signature=exclusion_signature,
    ).rstrip("/")
    pct = config.get("prediction_variance_percentile")
    if pct is not None:
        base = f"{base}_gpvarp{str(pct).replace('.', 'p')}"
    # Cross-site harmonization is analysis-level: it changes the ANALYSIS dir only
    # (never the processed base dir). "none"/False keep the legacy path verbatim.
    harm = config.get("site_harmonization", "none")
    if harm and harm != "none":
        base = f"{base}_harm{harm}"
    if config.get("scoring_site_covariate"):
        base = f"{base}_sitecov"
    return base + "/"


# ---------------------------------------------------------------------------
# Analysis-output provenance manifest: guards against STALE analysis maps from a
# prior run (different code version, feature set, or a config that hashed to the
# same dir) being globbed into the objective. Written after a successful analysis
# write; checked (and stale maps cleared) before the next write to that dir.
# ---------------------------------------------------------------------------
_PROVENANCE_SENTINEL = ".zbrains_provenance.json"
# Map-generating source: any edit here invalidates on-disk analysis maps (catches
# UNCOMMITTED changes -- a content hash, not a git rev).
_CODE_VERSION_FILES = (
    "zbrains_staged.py", "src/zbrains/analysis.py", "src/zbrains/processing.py",
    "src/zbrains/normalization.py", "src/zbrains/dataset.py", "results/greedy_auroc.py",
)
_CODE_VERSION_CACHE = None


def _pipeline_code_version():
    """Short content hash of the map-generating source files, so an analysis dir
    written by a DIFFERENT code version is detected as stale (uncommitted edits too)."""
    global _CODE_VERSION_CACHE
    if _CODE_VERSION_CACHE is None:
        here = os.path.dirname(os.path.abspath(__file__))
        h = hashlib.sha1()
        for rel in _CODE_VERSION_FILES:
            try:
                with open(os.path.join(here, rel), "rb") as fh:
                    h.update(fh.read())
            except OSError:
                h.update(b"?")
        _CODE_VERSION_CACHE = h.hexdigest()[:16]
    return _CODE_VERSION_CACHE


def _analysis_provenance(config, cohort):
    """Provenance record for an analysis dir: code version + config hash + the
    feature list. Any change invalidates (clears) that dir's stale analysis maps."""
    feats = sorted(
        str(f) for f in ((getattr(cohort, "base_pipeline_settings", {}) or {}).get("features") or []))
    return {
        "code_version": _pipeline_code_version(),
        "config": hashlib.sha1(
            json.dumps(config, sort_keys=True, default=str).encode("utf-8")).hexdigest()[:16],
        "features": feats,
    }


def _guard_analysis_dir(output_directory, provenance, *, verbose=False):
    """If ``output_directory``'s provenance manifest is missing or differs from
    ``provenance`` (different code / config / feature-set), CLEAR its stale analysis
    map trees (``{method}_maps``) so a re-run or other-code run cannot mix stale maps
    into the objective. The symlinked base ``maps``/``structural`` are left intact
    (``*_maps`` never matches the ``maps`` symlink). Returns True on a clean match."""
    path = os.path.join(str(output_directory), _PROVENANCE_SENTINEL)
    try:
        with open(path, encoding="utf-8") as fh:
            existing = json.load(fh)
    except (OSError, ValueError):
        existing = None
    if existing == provenance:
        return True
    removed = 0
    for d in glob.glob(os.path.join(str(output_directory), "*", "*", "*_maps")):
        if os.path.isdir(d) and not os.path.islink(d):
            shutil.rmtree(d, ignore_errors=True)
            removed += 1
    if removed and verbose:
        print(f"    provenance mismatch -> cleared {removed} stale analysis map dir(s) "
              f"in {output_directory}")
    return False


def _stamp_analysis_dir(output_directory, provenance):
    """Write the provenance manifest (atomic) after a successful analysis write."""
    d = str(output_directory)
    if not os.path.isdir(d):
        return
    path = os.path.join(d, _PROVENANCE_SENTINEL)
    tmp = os.path.join(d, f".{_PROVENANCE_SENTINEL}.tmp.{os.getpid()}")
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(provenance, fh, sort_keys=True)
        os.replace(tmp, path)
    except OSError:
        try:
            if os.path.exists(tmp):
                os.remove(tmp)
        except OSError:
            pass


def _quant_sample_surface_missing_companions(features):
    lower_features = {str(feature).lower() for feature in (features or [])}
    missing_companions = [
        f"{base}*blur"
        for base in ("adc", "fa", "qt1")
        if base in lower_features and f"{base}*blur" not in lower_features
    ]
    return missing_companions


def _analysis_kwargs(config, control_correlation_threshold, *, features=None, verbose=False):
    """Build the kwargs passed to ``dataset.analyze`` for a config.

    ``prediction_variance_percentile`` is only forwarded when ``analyze`` actually
    accepts it (so the driver still runs before the GP-filter plumbing lands; a
    warning is logged if a filter is requested but unsupported).
    """
    kwargs = {name: config[name] for name in zb.ANALYSIS_AXIS_NAMES}
    kwargs["method"] = config["method"]
    # Staged correlation arms have already restricted the control dataset before
    # dataset normalization. Do not apply the analysis-time per-map filter again.
    if config.get("outlier_method") == "correlation":
        kwargs["control_correlation_filter"] = False
        forward_quantile = False
    else:
        forward_quantile = True
    # Retain support for explicitly requested analysis-time quantiles outside the
    # staged correlation arm, but never forward the removed absolute-threshold API.
    _ccq = config.get("control_correlation_quantile")
    if forward_quantile and _ccq is not None:
        _analyze_params = inspect.signature(_first_cohort_analyze_signature()).parameters
        if "control_correlation_quantile" in _analyze_params:
            kwargs["control_correlation_quantile"] = _ccq
        else:
            print(
                "  WARNING: control_correlation_quantile requested but "
                "dataset.analyze() does not accept it; ignoring it."
            )
    # ADC/FA/qT1 plain-sample surface (driver-managed; derived from sample_surface).
    if config.get("quant_sample_surface") is not None:
        analyze_params = inspect.signature(_first_cohort_analyze_signature()).parameters
        if "quant_sample_surface" in analyze_params:
            missing_companions = _quant_sample_surface_missing_companions(features)
            if missing_companions:
                if verbose:
                    print(
                        "  NOTE: not forwarding quant_sample_surface="
                        f"{config['quant_sample_surface']!r}; feature set is missing "
                        "required companions: " + ", ".join(missing_companions)
                    )
            else:
                kwargs["quant_sample_surface"] = config["quant_sample_surface"]

    pct = config.get("prediction_variance_percentile")
    if pct is not None:
        analyze_params = inspect.signature(_first_cohort_analyze_signature()).parameters
        if "prediction_variance_percentile" in analyze_params:
            kwargs["prediction_variance_percentile"] = pct
        else:
            print(
                "  WARNING: prediction_variance_percentile requested but "
                "dataset.analyze() does not accept it yet; ignoring the GP filter."
            )
    return kwargs


def _first_cohort_analyze_signature():
    """Return zbdataset.analyze for signature introspection (import-light)."""
    from zbrains.dataset import zbdataset

    return zbdataset.analyze


# ---------------------------------------------------------------------------
# Cohort definition and per-candidate evaluation
# ---------------------------------------------------------------------------
@dataclass
class Cohort:
    """One dataset participating in the pooled objective.

    ``score_name`` is the name understood by :mod:`results.greedy_auroc`
    (``"MICs"`` / ``"NOEL"``) used to locate lesion ground truth.
    """

    name: str
    score_name: str
    control_dataset: object
    patient_dataset: object
    output_dir_prefix: str
    base_pipeline_settings: dict


# --- Step 8: pre-dataset-norm control exclusion --------------------------------
# Both outlier methods inspect the full, subject-normalized cohort BEFORE
# dataset-level normalization. ROBPCA removes whole control scans. Correlation
# QC removes only the affected subject-feature observations from RAVEL/Nyul and
# from the normative reference; every unaffected feature remains available.

# Filenames like: sub-X_ses-Y_hemi-L_surf-fsLR-32k_label-white_feature-T1w_smooth-5mm.func.gii
_QC_MAP_RE = re.compile(
    r"_hemi-(?P<hemi>[LR])_surf-fsLR-32k_label-white_feature-(?P<feat>.+?)_smooth-\d+mm\.func\.gii$"
)

# Cache each decision per subject-normalized detector base, method, and parameter
# so QC runs at most once per cohort x arm across the whole greedy search.
_EXCLUSION_CACHE = {}


def _exclusion_signature(exclusions):
    """Stable signature of whole-subject or feature-specific exclusions.

    Whole-subject ROBPCA exclusions are an iterable of ``(pid, ses)`` pairs.
    Correlation exclusions are ``{feature: iterable[(pid, ses)]}``. Empty input
    preserves the legacy path exactly.
    """
    if not exclusions:
        return ""
    if hasattr(exclusions, "items"):
        lines = [
            f"{feature}:{pid}/{ses}"
            for feature, pairs in exclusions.items()
            for pid, ses in pairs
        ]
    else:
        lines = [f"{pid}/{ses}" for pid, ses in exclusions]
    joined = "\n".join(sorted(lines))
    return "excl" + hashlib.sha1(joined.encode("utf-8")).hexdigest()[:10]


def _control_pairs(control_dataset):
    """(pid, ses) string pairs from a control dataset's demographics."""
    data = control_dataset.demographics.data
    return [(str(pid), str(ses)) for pid, ses in
            zip(data["participant_id"], data["session_id"])]


def _control_manifest(base_dir, control_dataset):
    """Build a run_surface_qc manifest from the fsLR-32k white-surface intensity
    maps of the (full) control cohort. One row per subject x feature x hemi;
    ``subject`` encodes 'pid|ses' so multi-session controls never collapse. Only
    surf-fsLR-32k / label-white intensity maps match (composite *blur companions
    are fsnative and are excluded automatically)."""
    rows = []
    for pid, ses in _control_pairs(control_dataset):
        cortex = os.path.join(base_dir, pid, ses, "maps", "cortex")
        seen = set()
        for path in sorted(glob.glob(os.path.join(cortex, "*.func.gii"))):
            m = _QC_MAP_RE.search(os.path.basename(path))
            if not m:
                continue
            feat, hemi = m.group("feat"), m.group("hemi")
            if "*blur" in feat:  # defensive; label-white blur maps don't exist
                continue
            key = (feat, hemi)
            if key in seen:  # one smooth per (feature, hemi)
                continue
            seen.add(key)
            rows.append({"subject": f"{pid}|{ses}", "feature": feat,
                         "hemisphere": hemi, "file_path": path})
    return pd.DataFrame(rows)


def _correlation_quantile_feature_exclusions(
    manifest,
    quantile,
    *,
    minimum_controls=SURFACE_QC_MIN_CONTROLS,
):
    """Return ``{feature: frozenset((pid, ses))}`` for pre-norm masking.

    Each feature is ranked independently. A subject's bilateral left/right maps
    are concatenated, its average correlation to the other controls is computed,
    and the bottom ``quantile`` fraction is masked for that feature only. The
    fixed-count rank selection makes q05/q10/q20 deterministic; ties are broken
    by subject id. At least ``minimum_controls`` rankable observations are kept
    for every feature.
    """
    import math
    import numpy as np
    from zbrains.analysis import _subject_average_correlations
    from zbrains.surface_qc import load_gifti_metric

    q = float(quantile)
    if not 0.0 <= q < 1.0:
        raise ValueError("control_correlation_quantile must be in [0, 1)")
    if manifest is None or manifest.empty or q == 0.0:
        return {}

    exclusions = {}
    for feature in sorted(manifest["feature"].astype(str).unique()):
        feature_rows = manifest[manifest["feature"].astype(str) == feature]
        bilateral = []
        for subject in sorted(feature_rows["subject"].astype(str).unique()):
            subject_rows = feature_rows[feature_rows["subject"].astype(str) == subject]
            paths = {}
            for hemi in ("L", "R"):
                rows = subject_rows[subject_rows["hemisphere"].astype(str) == hemi]
                if len(rows):
                    paths[hemi] = str(rows.iloc[0]["file_path"])
            if set(paths) != {"L", "R"}:
                continue
            try:
                left = np.asarray(load_gifti_metric(paths["L"]), dtype=float).reshape(-1)
                right = np.asarray(load_gifti_metric(paths["R"]), dtype=float).reshape(-1)
            except Exception:
                continue
            bilateral.append((subject, np.concatenate([left, right])))

        if len(bilateral) < 2:
            continue
        # Keep the modal vertex count if malformed inputs disagree in length.
        counts = pd.Series([values.size for _, values in bilateral]).value_counts()
        modal_size = int(counts.index[0])
        bilateral = [(subject, values) for subject, values in bilateral
                     if values.size == modal_size]
        if len(bilateral) < 2:
            continue

        subjects = [subject for subject, _ in bilateral]
        correlations = _subject_average_correlations(
            np.vstack([values for _, values in bilateral])
        )
        floor = min(max(0, int(minimum_controls)), len(subjects))
        n_drop = min(int(math.ceil(q * len(subjects))), len(subjects) - floor)
        if n_drop <= 0:
            continue
        order = sorted(
            range(len(subjects)),
            key=lambda i: (
                float(correlations[i]) if np.isfinite(correlations[i]) else -np.inf,
                subjects[i],
            ),
        )
        removed = []
        for i in order[:n_drop]:
            if "|" not in subjects[i]:
                continue
            pid, ses = subjects[i].split("|", 1)
            removed.append((pid, ses))
        if removed:
            exclusions[feature] = frozenset(removed)
    return exclusions


def _resolve_excluded_pairs(cohort, config, env, *, verbose=True):
    """Resolve pre-dataset-normalization control exclusions for one arm.

    Correlation returns ``{feature: frozenset((pid, ses))}``; ROBPCA returns a
    whole-control ``frozenset((pid, ses))``. Both inspect the SUBJECT-NORM base
    with dataset_norm hardcoded to ``none``, never the composed dataset-norm base.
    """
    method = config.get("outlier_method", "none")
    if method not in {"correlation", "robpca"}:
        return frozenset()
    parameter = (config.get("control_correlation_quantile")
                 if method == "correlation" else config.get("qc_alpha"))
    subject_norm = config.get("subject_norm", "none")
    subj_norm_label = _compose_norm(subject_norm, "none")  # dataset_norm forced none
    # Outlier detection must be a STABLE, smoothing-independent decision, so run it on
    # UNSMOOTHED (0mm) fsLR-32k white-surface maps at the locked subject_norm. Otherwise
    # the surface-QC input is smoothed at the candidate's cortical_smoothing, so every
    # smoothing arm re-runs ROBPCA on differently-smoothed surfaces, yields a different
    # excluded set, and forces a full detector-base + RAVEL-base + reference-CV rebuild
    # (nothing shared across smoothing arms). Forcing 0mm gives ONE detector base
    # (baseline signature) shared by every smoothing/sampling arm.
    detector_config = dict(config)
    detector_config["cortical_smoothing"] = 0
    detector_config["hippocampal_smoothing"] = 0
    processing_signature = _processing_signature(detector_config, "")
    subj_base = zb.processed_base_directory_for(
        subj_norm_label, cohort.output_dir_prefix, processing_signature
    )
    assert "_excl" not in subj_base, "detector base must never be exclusion-keyed"
    cache_key = (subj_base, method, parameter)
    if cache_key in _EXCLUSION_CACHE:
        return _EXCLUSION_CACHE[cache_key]

    # Guarantee the FULL-control subject-norm base exists (never a silent no-op).
    # Built with detector_config so the QC surfaces are UNSMOOTHED (0mm identity copy).
    if not zb.processed_maps_complete(subj_base, [cohort.control_dataset]):
        settings = _pipeline_settings(cohort, detector_config, normalization=subj_norm_label)
        if verbose:
            print(f"[{cohort.name}] Processing controls into {subj_norm_label} base "
                  f"for surface QC (0mm, unsmoothed): {subj_base}")
        cohort.control_dataset.process(
            output_directory=subj_base, env=env, verbose=verbose, **settings
        )

    manifest = _control_manifest(subj_base, cohort.control_dataset)
    control_pairs = set(_control_pairs(cohort.control_dataset))
    if method == "correlation":
        raw_exclusions = _correlation_quantile_feature_exclusions(
            manifest,
            parameter,
            minimum_controls=SURFACE_QC_MIN_CONTROLS,
        )
        excluded = {
            str(feature): frozenset(pair for pair in pairs if pair in control_pairs)
            for feature, pairs in raw_exclusions.items()
        }
        excluded = {feature: pairs for feature, pairs in excluded.items() if pairs}
        if verbose:
            total = sum(len(pairs) for pairs in excluded.values())
            detail = ", ".join(
                f"{feature}={len(pairs)}" for feature, pairs in sorted(excluded.items())
            ) or "none"
            print(
                f"[{cohort.name}] correlation QC (q={float(parameter):g}) masked "
                f"{total} subject-feature observations before dataset normalization "
                f"({detail})."
            )
    else:
        qc_alpha = parameter
        excluded_set = set()
        if not manifest.empty:
            from zbrains.surface_qc import run_surface_qc  # lazy: pulls sklearn/nibabel
            qc_dir = os.path.join(
                cohort.output_dir_prefix, "surface_qc",
                f"{_normlabel(subj_norm_label)}_{processing_signature or _smoothing_label(detector_config)}_a{qc_alpha}",
            )
            result = run_surface_qc(manifest, output_dir=qc_dir, qc_alpha=qc_alpha,
                                    minimum_controls=SURFACE_QC_MIN_CONTROLS,
                                    make_figures=False)
            if result.exclusions is not None and len(result.exclusions):
                for s in result.exclusions["subject"].astype(str).unique():
                    pid, ses = s.split("|", 1)
                    excluded_set.add((pid, ses))
        excluded = frozenset(pair for pair in excluded_set if pair in control_pairs)
        if verbose:
            print(f"[{cohort.name}] surface QC (qc_alpha={qc_alpha}) flagged "
                  f"{len(excluded)}/{len(control_pairs)} controls for exclusion.")
    _EXCLUSION_CACHE[cache_key] = excluded
    return excluded


def _normlabel(normalization):
    return zb._normalization_dir_label(normalization)


def _subset_control_dataset(orig, keep_pairs, name=None):
    """Build a zbdataset restricted to ``keep_pairs`` (a list of (pid,ses)).

    ``name`` overrides the dataset name. This matters for dataset-norm (RAVEL/Nyul):
    ``dataset.process`` FITS the model when the dataset name is a control alias
    (control/controls/hc/...) and APPLIES the frozen model otherwise. So the train-
    fold reference keeps the control name (it FITS the per-fold model), while the
    held-out fold is given a non-control name so it APPLIES that frozen fit (never
    re-fitting on itself) -- the point of per-fold refit."""
    from zbrains.dataset import demographics as _Demographics, zbdataset as _ZBDataset
    keep = sorted(tuple(p) for p in keep_pairs)
    if not keep:
        raise ValueError("empty control subset")
    ctrl_demo = _Demographics(
        orig.demographics.csv_file,
        column_mapping=orig.demographics.column_mapping,
        normative_columns=orig.demographics.normative_columns,
        normative_dtypes=orig.demographics.normative_dtypes,
        reference=orig.demographics.reference,
        subset=[list(p) for p in keep],
    )
    subset_dataset = _ZBDataset(
        name or orig.name, ctrl_demo, orig.micapipe_directory,
        hippunfold_directory=orig.hippunfold_directory,
        freesurfer_directory=orig.freesurfer_directory,
        raw_data_directory=getattr(orig, "raw_data_directory", None),
        cortex=orig.cortex, hippocampus=orig.hippocampus, subcortical=orig.subcortical,
    )
    allowed = set(keep)
    source_exclusions = getattr(orig, "control_feature_exclusions", None) or {}
    subset_dataset.control_feature_exclusions = {
        str(feature): frozenset(
            (str(pid), str(ses)) for pid, ses in pairs
            if (str(pid), str(ses)) in allowed
        )
        for feature, pairs in source_exclusions.items()
    }
    subset_dataset.control_feature_exclusions = {
        feature: pairs
        for feature, pairs in subset_dataset.control_feature_exclusions.items()
        if pairs
    }
    return subset_dataset


def _restricted_control_dataset(cohort, excluded):
    """Apply whole-control or feature-specific pre-normalization exclusions.

    ROBPCA subsets the cohort. Correlation clones the full cohort and attaches a
    feature mask, retaining every scan for its unaffected features. The resulting
    object drives both dataset normalization and normative scoring.
    """
    if not excluded:
        return cohort.control_dataset
    orig = cohort.control_dataset
    if hasattr(excluded, "items"):
        clone = _subset_control_dataset(orig, _control_pairs(orig))
        allowed = set(_control_pairs(orig))
        clone.control_feature_exclusions = {
            str(feature): frozenset(
                (str(pid), str(ses)) for pid, ses in pairs
                if (str(pid), str(ses)) in allowed
            )
            for feature, pairs in excluded.items()
        }
        clone.control_feature_exclusions = {
            feature: pairs for feature, pairs in clone.control_feature_exclusions.items()
            if pairs
        }
        return clone
    keep = [pair for pair in _control_pairs(orig) if pair not in excluded]
    if not keep:
        raise ValueError(f"[{cohort.name}] exclusion removed all controls")
    return _subset_control_dataset(orig, keep)


# Held-out-control scoring for the specificity-aware objective.
SPECIFICITY_CV_FOLDS = 5

# The K control-CV folds are independent analyze() passes; they are scored
# concurrently in THREADS (analyze is CPU+I/O with BLAS/gifti releasing the GIL,
# and is thread-safe -- no chdir/global state). Capped because each concurrent
# fold holds its reference control maps in RAM; the natural ceiling is the fold
# count anyway. Lower for memory, or set env.num_threads=1 to force serial.
CONTROL_FOLD_MAX_WORKERS = 5


def _control_feature_availability(base_dir, control_dataset):
    """{(pid, ses): frozenset(base feature names present)} from the processed base
    cortex maps on disk. Uses the fsLR-32k label-white intensity maps (the same
    maps _control_manifest scans), so it captures which subjects actually have
    each feature (e.g. FLAIR, qT1) -- used to balance the CV folds on feature
    availability. Blur companions follow their base feature and are ignored here."""
    masks = getattr(control_dataset, "control_feature_exclusions", None) or {}
    masked = {}
    for feature, pairs in masks.items():
        key = str(feature).lower().replace("*blur", "")
        if key == "t1map":
            key = "qt1"
        masked.setdefault(key, set()).update(map(tuple, pairs))
    avail = {}
    for pid, ses in _control_pairs(control_dataset):
        feats = set()
        cortex = os.path.join(base_dir, pid, ses, "maps", "cortex")
        for path in glob.glob(os.path.join(cortex, "*.func.gii")):
            m = _QC_MAP_RE.search(os.path.basename(path))
            if m and "*blur" not in m.group("feat"):
                feature = m.group("feat")
                key = feature.lower().replace("*blur", "")
                if key == "t1map":
                    key = "qt1"
                if (pid, ses) not in masked.get(key, set()):
                    feats.add(feature)
        avail[(pid, ses)] = frozenset(feats)
    return avail


def _control_folds(control_ds, k=SPECIFICITY_CV_FOLDS, seed=0, n_age_bins=3,
                   feature_availability=None):
    """Age/sex/feature-availability-BALANCED K-fold partition at SUBJECT level.

    Returns {fold_id: [test (pid,ses) pairs]}. ALL sessions of a held-out subject
    are placed in the SAME fold, so a longitudinal control can never straddle
    train/test and leak into the normative reference (train = the complement). The
    split is over unique subjects, stratified on a composite stratum so each fold
    spans the cohort:
      * SEX  and a coarse AGE bin (the z-score arm doesn't condition on age/sex, so
        unbalanced folds would bias its control false positives), AND
      * a FEATURE-AVAILABILITY signature (``feature_availability`` = {(pid,ses):
        frozenset(features)}) restricted to the features that VARY across subjects
        -- so a sparse feature (e.g. FLAIR present in ~half the cohort) is spread
        evenly and no fold's train reference is starved of it.
    Rare composite strata are collapsed; falls back to plain KFold if the strata
    cannot support K splits."""
    import numpy as np
    from sklearn.model_selection import StratifiedKFold, KFold
    import collections

    data = control_ds.demographics.data
    pairs = _control_pairs(control_ds)                     # aligned to data rows
    # Group row indices by SUBJECT (pid); every session of a subject travels together.
    subj_rows = collections.OrderedDict()
    for i, (pid, _ses) in enumerate(pairs):
        subj_rows.setdefault(pid, []).append(i)
    subjects = list(subj_rows.keys())
    kk = max(2, min(int(k), len(subjects)))

    # One representative stratum per subject (its first row's sex x age-bin).
    parts = []
    if "SEX" in data.columns:
        parts.append(data["SEX"].astype(str).to_numpy())
    if "AGE" in data.columns:
        age = data["AGE"].to_numpy(dtype=float)
        try:
            bins = pd.qcut(age, min(n_age_bins, max(1, len(np.unique(age)))),
                           labels=False, duplicates="drop")
            parts.append(np.asarray(bins).astype(str))
        except Exception:
            pass
    if parts:
        row_strat = ["|".join(str(p[i]) for p in parts) for i in range(len(pairs))]
        subj_strat = [row_strat[subj_rows[s][0]] for s in subjects]
    else:
        subj_strat = ["_" for _ in subjects]

    # Feature-availability dimension: signature over the features that VARY across
    # subjects (universal features carry no balancing information).
    if feature_availability:
        per_subj_feats = {
            s: set().union(*[set(feature_availability.get(pairs[r], frozenset()))
                             for r in subj_rows[s]])
            for s in subjects
        }
        everseen = set().union(*per_subj_feats.values()) if per_subj_feats else set()
        universal = (set.intersection(*[per_subj_feats[s] for s in subjects])
                     if subjects else set())
        variable = everseen - universal
        if variable:
            for idx, s in enumerate(subjects):
                sig = "+".join(sorted(per_subj_feats[s] & variable)) or "none"
                subj_strat[idx] = f"{subj_strat[idx]}|f={sig}"

    subj_strat = np.array(subj_strat)
    if len(set(subj_strat)) > 1:
        counts = collections.Counter(subj_strat)
        subj_strat = np.array([s if counts[s] >= kk else "_rare" for s in subj_strat])
    else:
        subj_strat = None

    def _expand(test_subject_idx):
        out = []
        for si in test_subject_idx:
            out.extend(pairs[r] for r in subj_rows[subjects[si]])
        return out

    folds = {}
    try:
        if subj_strat is not None and len(set(subj_strat)) > 1:
            splitter = StratifiedKFold(n_splits=kk, shuffle=True, random_state=seed)
            it = splitter.split(np.zeros(len(subjects)), subj_strat)
        else:
            it = KFold(n_splits=kk, shuffle=True, random_state=seed).split(subjects)
        for i, (_tr, te) in enumerate(it):
            folds[i] = _expand(te)
    except ValueError:
        it = KFold(n_splits=kk, shuffle=True, random_state=seed).split(subjects)
        for i, (_tr, te) in enumerate(it):
            folds[i] = _expand(te)
    return folds


def _score_heldout_controls(cohort, config, env, control_correlation_threshold,
                            control_dataset, exclusion_signature, *,
                            site_harmonizer=None, verbose=False):
    """Score each control fold against the OTHER controls (leakage-free) so the
    held-out controls get |z| maps that serve as the specificity negative class.
    Reuses the same analysis machinery/base as patient scoring (analysis-level).

    Returns ``{fold_id: [held-out (participant, session) pairs]}`` for the folds
    that were ACTUALLY scored (skipped folds -- too few training controls -- wrote
    no maps and are omitted). The objective uses this to compute one AUROC per fold
    and lock on their mean (reporting the SD), instead of pooling all controls into
    a single value.

    Cost note: this adds K analyze passes per cohort per candidate. There is no
    useful cross-candidate cache -- EVERY greedy config field (normalization,
    sampling, blur, smoothing, method, harmonization, exclusion) changes the shared
    control pipeline, so distinct candidates never share control maps. The expensive
    base processing is already shared across candidates; only the cheap analysis
    layer repeats, and the incumbent is not re-scored across stages (score_cache).
    Lower SPECIFICITY_CV_FOLDS if the GP scoring arm's per-fold refit is too slow."""
    normalization = config["normalization"]
    processing_signature = _processing_signature(config, exclusion_signature)
    base_signature = _base_signature(config, exclusion_signature)   # base: exclusion dropped if dataset_norm=none
    output_directory = staged_output_directory(
        config, cohort.output_dir_prefix, control_correlation_threshold,
        exclusion_signature=processing_signature,
    )
    settings = _pipeline_settings(cohort, config, normalization=normalization)
    akw = _analysis_kwargs(
        config, control_correlation_threshold,
        features=getattr(cohort, "base_pipeline_settings", {}).get("features"),
        verbose=verbose,
    )
    # Per-subject feature availability (from the processed base maps) so folds can
    # be balanced on it -- otherwise a fold could hold out most subjects that have
    # a sparse feature (e.g. FLAIR), starving that feature's reference in a fold.
    base_dir = zb.processed_base_directory_for(
        normalization, cohort.output_dir_prefix, base_signature
    )
    feature_availability = _control_feature_availability(base_dir, control_dataset)

    all_pairs = _control_pairs(control_dataset)

    # Phase 1 (SEQUENTIAL, deterministic): select the folds to score, warning about
    # any that must be skipped. Kept out of the parallel section so warning order and
    # the skip decision are unaffected by concurrent execution.
    fold_jobs = []   # (fold_id, train_pairs, sorted_test_pairs)
    for _fold_id, test_pairs in _control_folds(
        control_dataset, feature_availability=feature_availability
    ).items():
        test_set = set(map(tuple, test_pairs))
        train_pairs = [p for p in all_pairs if p not in test_set]
        if len(train_pairs) < 2 or not test_set:
            # Too few training controls to build a reference for this fold. On a
            # normal-sized cohort this never fires; warn so a tiny cohort silently
            # contributing fewer negatives is visible rather than a silent gap.
            warnings.warn(
                f"{cohort.score_name}: control fold skipped (only {len(train_pairs)} "
                f"training controls) -- {len(test_set)} held-out controls will not "
                f"contribute to the specificity negative pool.",
                RuntimeWarning, stacklevel=2,
            )
            continue
        fold_jobs.append((_fold_id, train_pairs, sorted(test_set)))

    # Phase 2: score each fold. Folds are INDEPENDENT -- each builds its OWN
    # train/heldout subset datasets and analyze() writes |z| maps only for that
    # fold's (disjoint) held-out controls, so they run concurrently in THREADS with
    # output identical to serial (analyze is CPU+I/O -- BLAS/gifti release the GIL --
    # and thread-safe: no chdir/global state).
    def _score_one_fold(job):
        fold_id, train_pairs, test_pairs = job
        train_ds = _subset_control_dataset(control_dataset, train_pairs)
        heldout_ds = _subset_control_dataset(control_dataset, test_pairs)
        # Validate BOTH the held-out set AND the train reference (mirrors the main
        # patient path, which validates patient + control). Without validating
        # train_ds its .features stays [] and analyze() rejects the fold with
        # "Feature sets don't match between patient and reference datasets."
        train_ds.validate(output_directory=output_directory, verbose=False, **settings)
        heldout_ds.validate(output_directory=output_directory, verbose=False, **settings)
        heldout_ds.analyze(
            output_directory=output_directory,
            reference=train_ds,
            site=cohort.score_name,
            site_harmonizer=site_harmonizer,
            scoring_site_covariate=config.get("scoring_site_covariate", False),
            **akw,
        )
        # Only folds that were actually scored contribute negatives; skipped folds
        # wrote no |z| maps, so their held-out controls must NOT be used downstream.
        return fold_id, test_pairs

    n_fold_jobs = min(len(fold_jobs),
                      max(1, int(getattr(env, "num_threads", 1) or 1)),
                      CONTROL_FOLD_MAX_WORKERS)
    if n_fold_jobs > 1:
        from joblib import Parallel, delayed   # heavy import, lazy
        results = Parallel(n_jobs=n_fold_jobs, backend="threading")(
            delayed(_score_one_fold)(job) for job in fold_jobs
        )
    else:
        results = [_score_one_fold(job) for job in fold_jobs]
    return {fold_id: test_pairs for fold_id, test_pairs in results}


def _ensure_fold_base(cohort, config, env, reference_k, heldout_k, exclusion_signature,
                      fold_id, *, verbose=False):
    """Build a PER-FOLD processed base for dataset-norm (RAVEL/Nyul).

    The dataset-norm model is FIT on ``reference_k`` (the train-fold controls, which
    keep a control name so ``dataset.process`` fits) and APPLIED FROZEN to the
    held-out controls + patients (given non-control names, so process applies rather
    than refits). This removes the second-order leak where the held-out fold
    influenced its own normalization. Keyed by ``..._reffold{fold_id}`` so each fold
    gets its own base; sentinel-guarded so it is built once and reused."""
    normalization = config["normalization"]
    base_fold_sig = f"{_base_signature(config, exclusion_signature)}_reffold{fold_id}"
    base_dir = zb.processed_base_directory_for(
        normalization, cohort.output_dir_prefix, base_fold_sig)
    settings = _pipeline_settings(cohort, config, normalization=normalization)
    required = list(_control_pairs(reference_k)) + list(_control_pairs(cohort.patient_dataset))
    if heldout_k is not None:
        required += list(_control_pairs(heldout_k))
    if zb.base_is_marked_complete(base_dir, required):
        return base_dir
    zb.unmark_base_complete(base_dir)
    # 1. FIT dataset-norm on the train-fold reference (control-named -> fits the model).
    if not zb.processed_maps_complete(base_dir, [reference_k]):
        reference_k.process(output_directory=base_dir, env=env, verbose=verbose,
                            skip_existing=False, **settings)
    # 2. APPLY the frozen fit to the held-out controls (non-control name -> applies).
    if heldout_k is not None and not zb.processed_maps_complete(base_dir, [heldout_k]):
        heldout_k.process(output_directory=base_dir, env=env, verbose=verbose,
                          skip_existing=True, **settings)
    # 3. APPLY the frozen fit to patients.
    if not zb.processed_maps_complete(base_dir, [cohort.patient_dataset]):
        cohort.patient_dataset.process(output_directory=base_dir, env=env, verbose=verbose,
                                       skip_existing=True, **settings)
    present = (zb.processed_maps_complete(base_dir, [reference_k])
               and zb.processed_maps_complete(base_dir, [cohort.patient_dataset])
               and (heldout_k is None or zb.processed_maps_complete(base_dir, [heldout_k])))
    if present:
        zb.mark_base_complete(base_dir, required)
    return base_dir


def _export_fold_control_features(cohort, config, control_correlation_threshold,
                                  reference_k, base_k, fold_root, export, *, verbose=False):
    """Load+transform the train-fold reference's control features into ``export`` (for
    the per-fold site harmonizer). Mirrors the two-pass export but fold-aware: the
    base maps come from THIS fold's base ``base_k`` and features are the reference_k
    controls only, so the harmonizer never sees the held-out fold."""
    settings = _pipeline_settings(cohort, config, normalization=config["normalization"])
    akw = _analysis_kwargs(
        config, control_correlation_threshold,
        features=getattr(cohort, "base_pipeline_settings", {}).get("features"), verbose=verbose)
    zb.symlink_processed_outputs(base_k, fold_root,
                                 datasets=[reference_k, cohort.patient_dataset])
    reference_k.validate(output_directory=fold_root, verbose=False, **settings)
    reference_k.analyze(
        output_directory=fold_root, reference=reference_k, site=cohort.score_name,
        export_control_features=export, **akw)


def _analyze_reference_cv(ctx, config, env, control_correlation_threshold, *,
                          needs_controls, verbose=False):
    """Reference cross-validation across cohorts, PER-FOLD refitting BOTH the
    dataset-norm (RAVEL/Nyul, per cohort) AND the cross-site harmonizer (pooled) on
    the train fold only -- so a held-out fold never influences its own normalization,
    harmonization, or z-reference (no second-order leakage).

    ``ctx`` = list of ``(cohort, control_ds, exclusion_sig)``. Returns
    ``(fold_output_dirs, control_subjects_by_name)`` = ``{score_name: {fold: root}}``
    / ``{score_name: {fold: [held-out pairs]}}`` for :func:`reference_cv_score_fn`.
    """
    cohorts = [c for c, _, _ in ctx]
    dataset_norm = config.get("dataset_norm", "none")
    harmonize = config.get("site_harmonization", "none") != "none"

    # --- Per-cohort fold plan (aligned fold ids). dataset_norm != none -> each fold
    # gets its OWN base (fit on that fold's train controls); else all folds share the
    # one exclusion-independent base built by _ensure_base_processed. ---
    plans = {}          # score_name -> {fold_id: {train,test,reference_k,heldout_k,base_k,fold_root}}
    fold_ids = set()
    for cohort, control_ds, sig in ctx:
        processing_signature = _processing_signature(config, sig)
        base_output = staged_output_directory(
            config, cohort.output_dir_prefix, control_correlation_threshold,
            exclusion_signature=processing_signature).rstrip("/")
        shared_base = zb.processed_base_directory_for(
            config["normalization"], cohort.output_dir_prefix, _base_signature(config, sig))
        feat_avail = _control_feature_availability(shared_base, control_ds)
        all_pairs = _control_pairs(control_ds)
        cohort_plan = {}
        for fold_id, test_pairs in _control_folds(
                control_ds, feature_availability=feat_avail).items():
            test_set = set(map(tuple, test_pairs))
            train_pairs = [p for p in all_pairs if p not in test_set]
            if len(train_pairs) < 2 or not test_set:
                warnings.warn(
                    f"{cohort.score_name}: reference fold {fold_id} skipped (only "
                    f"{len(train_pairs)} training controls).", RuntimeWarning, stacklevel=2)
                continue
            reference_k = _subset_control_dataset(control_ds, train_pairs)   # control-named -> FITS
            heldout_k = (_subset_control_dataset(control_ds, sorted(test_set), name="heldout")
                         if needs_controls else None)                       # non-control -> APPLIES frozen
            base_k = (shared_base if dataset_norm == "none"
                      else _ensure_fold_base(cohort, config, env, reference_k, heldout_k,
                                             sig, fold_id, verbose=verbose))
            fold_root = f"{base_output}_cvfold{fold_id}/"
            # Clear stale analysis maps if this fold root was last written by different
            # code / config / feature-set (before any export/scoring writes into it).
            _guard_analysis_dir(fold_root, _analysis_provenance(config, cohort), verbose=verbose)
            cohort_plan[fold_id] = dict(
                train=train_pairs, test=sorted(test_set), reference_k=reference_k,
                heldout_k=heldout_k, base_k=base_k, fold_root=fold_root)
        plans[cohort.score_name] = cohort_plan
        fold_ids |= set(cohort_plan)

    # --- Per-fold cross-site harmonizer, fit on the POOLED train-fold controls only. ---
    harmonizers = {}
    if harmonize:
        for fold_id in sorted(fold_ids):
            h = _build_site_harmonizer(config, cohorts)
            for cohort, control_ds, sig in ctx:
                fp = plans[cohort.score_name].get(fold_id)
                if not fp:
                    continue
                export = {}
                _export_fold_control_features(
                    cohort, config, control_correlation_threshold,
                    fp["reference_k"], fp["base_k"], fp["fold_root"], export, verbose=verbose)
                h.add_site_features(cohort.score_name, export)
            harmonizers[fold_id] = h

    # --- Score each cohort's patients + held-out controls per fold vs reference_k. ---
    fold_output_dirs = {c.score_name: {} for c in cohorts}
    control_subjects_by_name = {c.score_name: {} for c in cohorts}
    for cohort, control_ds, sig in ctx:
        settings = _pipeline_settings(cohort, config, normalization=config["normalization"])
        akw = _analysis_kwargs(
            config, control_correlation_threshold,
            features=getattr(cohort, "base_pipeline_settings", {}).get("features"), verbose=verbose)
        scov = config.get("scoring_site_covariate", False)
        patient_pairs = _control_pairs(cohort.patient_dataset)
        for fold_id, fp in plans[cohort.score_name].items():
            fold_root, reference_k = fp["fold_root"], fp["reference_k"]
            harmonizer_k = harmonizers.get(fold_id)
            # Symlink THIS fold's base maps in; clone the patient dataset (own
            # .features/.valid_subjects) so nothing mutates the shared object.
            zb.symlink_processed_outputs(
                fp["base_k"], fold_root, datasets=[control_ds, cohort.patient_dataset])
            patient_ds = _subset_control_dataset(
                cohort.patient_dataset, patient_pairs, name=cohort.patient_dataset.name)
            reference_k.validate(output_directory=fold_root, verbose=False, **settings)
            patient_ds.validate(output_directory=fold_root, verbose=False, **settings)
            patient_ds.analyze(
                output_directory=fold_root, reference=reference_k, site=cohort.score_name,
                site_harmonizer=harmonizer_k, scoring_site_covariate=scov, **akw)
            if needs_controls:
                heldout_k = fp["heldout_k"]
                heldout_k.validate(output_directory=fold_root, verbose=False, **settings)
                heldout_k.analyze(
                    output_directory=fold_root, reference=reference_k, site=cohort.score_name,
                    site_harmonizer=harmonizer_k, scoring_site_covariate=scov, **akw)
            _stamp_analysis_dir(fold_root, _analysis_provenance(config, cohort))
            fold_output_dirs[cohort.score_name][fold_id] = fold_root
            control_subjects_by_name[cohort.score_name][fold_id] = fp["test"]
    return fold_output_dirs, control_subjects_by_name


def _ensure_base_processed(cohort, config, env, *, control_dataset=None,
                           exclusion_signature="", verbose=True, reprocess_controls=False):
    """Make sure the cohort's exclusion-keyed base for ``normalization`` exists.

    ``normalization`` selects the reusable processed base (zbrains_SELF for
    "none", zbrains_WB for "whitestripe", zbrains_RAVEL for "ravel", plus composed
    labels). ``exclusion_signature`` further keys the base so distinct exclusion
    arms never reuse a stale RAVEL/Nyul fit; ``control_dataset`` (whole-control
    subset for ROBPCA, feature-masked clone for correlation) drives that fit. Runs
    the heavy image
    processing ONCE per (cohort, normalization, exclusion); candidates reuse via
    symlink.
    """
    control_dataset = control_dataset or cohort.control_dataset
    # When dataset_norm=='none' the base is exclusion-INDEPENDENT and SHARED across
    # exclusion arms (see _base_signature), so it must contain ALL controls. The
    # exclusion is then applied logically to the normative reference during
    # validation. Building the shared base from a whole-control subset would leave
    # it missing scans that a later full-control config needs. RAVEL/Nyul
    # (dataset_norm != none) DO fit on the subset/feature-masked control dataset.
    if config.get("dataset_norm", "none") == "none":
        control_dataset = cohort.control_dataset
    normalization = config["normalization"]
    # Base is keyed by _base_signature (exclusion dropped when dataset_norm=='none'),
    # so ROBPCA/exclusion arms reuse the shared base instead of re-processing.
    base_signature = _base_signature(config, exclusion_signature)
    settings = _pipeline_settings(cohort, config, normalization=normalization)
    base_dir = zb.processed_base_directory_for(
        normalization, cohort.output_dir_prefix, base_signature
    )

    # The exact subjects this base must contain (controls actually processed here +
    # patients). Recorded in the sentinel so it certifies THIS subject set -- adding
    # subjects to the participant CSV then re-running correctly re-processes them
    # rather than being masked by a stale 'complete' mark.
    required_pairs = list(_control_pairs(control_dataset)) + list(_control_pairs(cohort.patient_dataset))

    # Fast, RACE-SAFE completion check: the per-base sentinel is written (atomically)
    # only after BOTH controls and patients finish AND covers all required subjects,
    # so if it matches the base is fully done -- trust it in O(1) without re-scanning,
    # and never treat a base a PEER machine is still writing (or one missing newly-added
    # subjects) as complete (shared storage, forward/reverse runs).
    if not reprocess_controls and zb.base_is_marked_complete(base_dir, required_pairs):
        return base_dir
    # We may (re)write this base -> drop any stale sentinel first so a crash mid-build
    # can't leave a base marked complete while it is actually partial.
    zb.unmark_base_complete(base_dir)

    need_controls = reprocess_controls or not zb.processed_maps_complete(
        base_dir, [control_dataset]
    )
    if need_controls:
        if verbose:
            print(f"[{cohort.name}] Processing controls into {normalization} base {base_dir}")
        # reprocess_controls forces a full rebuild; otherwise per-subject reuse
        # skips controls already processed in this base (only new/missing ones run).
        control_dataset.process(
            output_directory=base_dir, env=env, verbose=verbose,
            skip_existing=not reprocess_controls, **settings
        )

    if not zb.processed_maps_complete(base_dir, [cohort.patient_dataset]):
        if verbose:
            print(f"[{cohort.name}] Processing patients into {normalization} base {base_dir}")
        cohort.patient_dataset.process(
            output_directory=base_dir, env=env, verbose=verbose,
            skip_existing=True, **settings
        )

    # Stamp complete ONLY when both datasets are verified fully present -- so a run
    # where some subjects failed to process is NOT falsely trusted (it re-attempts
    # next time). This also ADOPTS a pre-sentinel base that was already fully
    # processed (nothing was rebuilt above; we just mark it -> ZERO reprocessing).
    if (zb.processed_maps_complete(base_dir, [control_dataset])
            and zb.processed_maps_complete(base_dir, [cohort.patient_dataset])):
        zb.mark_base_complete(base_dir, required_pairs)
    return base_dir


def _analyze_candidate(cohort, config, env, control_correlation_threshold, *,
                       control_dataset=None, exclusion_signature="", verbose=False,
                       site_harmonizer=None, export_control_features=None):
    """Symlink the base maps into the candidate's analysis dir and score-map it.

    Returns the analysis ``output_directory`` (where ``{method}_maps`` land). The
    normalization selects which processed base to symlink from; ``control_dataset``
    (the restricted cohort) is the normative reference so the same exclusion drives
    both the dataset-norm fit and the normative fit.

    Cross-site harmonization (two-pass): with ``export_control_features`` (a dict)
    this runs the EXPORT pass -- it loads+transforms this cohort's control features
    and stashes them per feature-map key, skipping patient scoring. With
    ``site_harmonizer`` it runs the SCORING pass -- patients are scored against the
    pooled harmonized control reference.
    """
    control_dataset = control_dataset or cohort.control_dataset
    normalization = config["normalization"]
    # Analysis dir keeps the full exclusion signature (its z-reference depends on the
    # excluded set); the base dir uses _base_signature (exclusion dropped when
    # dataset_norm=='none') so exclusion arms symlink from the shared base.
    processing_signature = _processing_signature(config, exclusion_signature)
    base_signature = _base_signature(config, exclusion_signature)
    output_directory = staged_output_directory(
        config, cohort.output_dir_prefix, control_correlation_threshold,
        exclusion_signature=processing_signature,
    )
    base_dir = zb.processed_base_directory_for(
        normalization, cohort.output_dir_prefix, base_signature
    )
    settings = _pipeline_settings(cohort, config, normalization=normalization)

    # Provenance guard: clear stale {method}_maps if this dir was last written by a
    # different code version / config / feature-set (prevents orphaned maps from a
    # prior run being globbed into the objective). Base maps/structural untouched.
    provenance = _analysis_provenance(config, cohort)
    _guard_analysis_dir(output_directory, provenance, verbose=verbose)

    if Path(output_directory).resolve() != Path(base_dir).resolve():
        zb.symlink_processed_outputs(
            base_dir, output_directory,
            datasets=[control_dataset, cohort.patient_dataset],
        )

    control_dataset.validate(
        output_directory=output_directory, verbose=False, **settings
    )
    cohort.patient_dataset.validate(
        output_directory=output_directory, verbose=False, **settings
    )
    cohort.patient_dataset.analyze(
        output_directory=output_directory,
        reference=control_dataset,
        export_control_features=export_control_features,
        site_harmonizer=site_harmonizer,
        site=cohort.score_name,
        scoring_site_covariate=config.get("scoring_site_covariate", False),
        **_analysis_kwargs(
            config,
            control_correlation_threshold,
            features=getattr(cohort, "base_pipeline_settings", {}).get("features"),
            verbose=verbose,
        ),
    )
    _stamp_analysis_dir(output_directory, provenance)
    return output_directory


# Objective comparison modes (each drives a SEPARATE greedy run):
#   within        -- lesion |z| vs the patient's OWN non-lesion vertices (all disease)
#   within_tle    -- as within, restricted to TLE lesions
#   within_fcd    -- as within, restricted to FCD lesions
#   control       -- lesion |z| vs ALL held-out-control vertices
#   disease_control -- lesion |z| vs the OTHER disease's patients' non-lesional vertices
#   both          -- ALIAS of within (within-subject localization over TLE+FCD together)
OBJECTIVE_MODES = ("within", "within_tle", "within_fcd", "control", "disease_control", "both")

# Modes that need held-out controls scored (the extra K-fold analyze passes).
# ONLY 'control': its negatives ARE controls, so they must be held out of the
# normative fit (else a negative control helped define "normal" -> deflated |z| ->
# inflated specificity). disease_control's negatives are other-disease patients
# (never in the control fit), and within/both compare within-subject, so none of
# them need held-out controls.
_CONTROL_MODES = ("control",)

# Reference cross-validation SCOPE (which objectives re-score patients against each
# control-fold-complement reference and report mean +/- SD, vs a cheap single-shot).
#   "control_disease" (DEFAULT) -- only 'control' is reference-CV'd. It is the ONLY
#                        mode with a leakage risk to cross-validate away: its
#                        negatives are controls, which also train the normative
#                        model, so held-out controls are needed. disease_control
#                        (negatives = other-disease patients, never in the fit) and
#                        within_*/both (within-subject) have no such leakage, so
#                        they stay single-shot -- avoiding the ~5x cost for a
#                        near-zero variance estimate.
#   "full"            -- every objective is reference-CV'd (within_* included).
#   "off"             -- no reference-CV; every objective is single-shot.
CV_SCOPES = ("control_disease", "full", "off")
# The only leakage-prone specificity test -- reference-CV'd under "control_disease".
_REFERENCE_CV_SPECIFICITY_MODES = ("control",)


def _mode_uses_reference_cv(objective_mode, cv_scope):
    """Whether ``objective_mode`` is reference-CV'd under ``cv_scope`` (else single-shot)."""
    if cv_scope not in CV_SCOPES:
        raise ValueError(f"cv_scope must be one of {CV_SCOPES}, got {cv_scope!r}")
    if cv_scope == "full":
        return True
    if cv_scope == "off":
        return False
    return objective_mode in _REFERENCE_CV_SPECIFICITY_MODES   # "control_disease"

# Side channel: reference-CV computes the objective ONCE PER reference fold and the
# greedy locks on their MEAN. The full per-fold vector + mean/SD is stashed here by
# reference_cv_score_fn so evaluate() can record it into the history CSV without
# changing the score_fn's float return contract. Serial greedy -> one slot is safe.
_LAST_CONTROL_FOLD_STATS = None


def _reset_fold_stats():
    global _LAST_CONTROL_FOLD_STATS, _LAST_FEATURE_AUROC
    _LAST_CONTROL_FOLD_STATS = None
    _LAST_FEATURE_AUROC = {}


# Companion side channel: the PER-FEATURE balanced AUROC of the objective's active
# leg ({feature: value}), stashed by default_score_fn so reference_cv_score_fn can
# average it across folds and record one column per feature in the history CSV
# (in addition to the overall mean). Reset with the fold stats each evaluation.
_LAST_FEATURE_AUROC = {}


def _feature_breakdown(detail, features, group_cols):
    """Per-feature balanced AUROC: ``{feature: balanced_macro_auroc(one feature)}``.

    Uses the SAME cell-averaging as the overall objective (``balanced_macro_auroc``)
    but restricted to each feature in turn, so the per-feature columns decompose the
    overall mean consistently. Empty when ``detail`` has no ``feature`` column."""
    from results.greedy_auroc import balanced_macro_auroc
    if detail is None or len(detail) == 0 or "feature" not in getattr(detail, "columns", []):
        return {}
    present = [str(f) for f in detail["feature"].dropna().unique()]
    if features is not None:
        allow = set(features)
        present = [f for f in present if f in allow]
    out = {}
    for feat in sorted(present):
        val = balanced_macro_auroc(detail, features={feat}, group_cols=group_cols)
        if pd.notna(val):
            out[feat] = float(val)
    return out


def _flatten_control_pairs(control_subjects):
    """A cohort's control-subject entry -> a flat list of (pid, ses) pairs.

    Reference-CV passes ONE fold's pairs (already a flat list). Legacy/mock callers
    may pass a ``{fold_id: [pairs]}`` dict (e.g. from the old _score_heldout_controls),
    which is flattened to the pooled set so the single-shot control leg still works."""
    if control_subjects is None:
        return None
    if isinstance(control_subjects, dict):
        seen, flat = set(), []
        for pairs in control_subjects.values():
            for p in (pairs or []):
                t = tuple(p)
                if t not in seen:
                    seen.add(t); flat.append(p)
        return flat
    return control_subjects


def default_score_fn(config, output_dirs_by_score_name, *, control_subjects_by_name=None,
                     objective_mode="both", aggregation="mean_subject_feature",
                     features=None, group_cols=("dataset", "lesion_type", "feature")):
    """Greedy objective: BALANCED macro AUROC over (dataset x disease x feature),
    computable in several comparison MODES (``objective_mode``):

    * ``"within"`` -- WITHIN-SUBJECT localization: each lesion's |z| vs that SAME
      patient's OWN non-lesion vertices (per_subject_feature_auc via
      pooled_macro_auroc). "Does the lesion stand out within this brain?".
    * ``"within_tle"`` / ``"within_fcd"`` -- as ``within`` but restricted to TLE-
      or FCD-type lesions (hippocampal ``hipp-*`` rows count as TLE).
    * ``"control"`` -- CROSS-SUBJECT specificity: lesion |z| vs the |z| of the
      held-out-control vertices at the SAME region (per_patient_vs_control_auc).
      ``control_subjects_by_name[name]`` is this evaluation's control negatives
      (one reference fold's held-out controls under reference-CV, or the pooled set).
    * ``"disease_control"`` -- DISEASE specificity: lesion |z| vs the OTHER
      disease's patients' NON-LESIONAL vertices (per_patient_vs_other_disease_auc);
      TLE vs FCD tissue and vice versa. Uses patients only (no held-out controls).
    * ``"both"`` -- ALIAS of ``within``: within-subject localization over BOTH
      TLE and FCD lesions together (identical to the plain ``within`` leg).

    This is a SINGLE-SHOT objective (one value for the given maps). Reference cross-
    validation -- computing it once per reference fold and averaging (mean +/- SD) --
    is done by :func:`reference_cv_score_fn`, which calls this per fold. In every
    mode each (dataset, lesion_type, feature) cell is averaged, then the cell means
    are averaged, so no cohort/disease/feature dominates the lock.
    """
    if objective_mode not in OBJECTIVE_MODES:
        raise ValueError(f"objective_mode must be one of {OBJECTIVE_MODES}, got {objective_mode!r}")

    global _LAST_FEATURE_AUROC
    _LAST_FEATURE_AUROC = {}   # per-feature breakdown of the active leg (this call)

    from results.greedy_auroc import (  # heavy module
        pooled_macro_auroc, balanced_macro_auroc, per_patient_vs_control_auc,
        per_patient_vs_other_disease_auc,
    )

    def _record(detail):
        """Stash the active leg's per-feature balanced AUROC so it can be recorded
        per feature in the history CSV (in addition to the overall mean)."""
        global _LAST_FEATURE_AUROC
        _LAST_FEATURE_AUROC = _feature_breakdown(detail, features, group_cols)
        return _LAST_FEATURE_AUROC

    def _within_leg(lesion_types=None):
        _scalar, detail = pooled_macro_auroc(
            output_dirs_by_score_name, method=config["method"],
            aggregation=aggregation, features=features, return_details=True,
        )
        if lesion_types is not None and detail is not None and not detail.empty:
            detail = detail[detail["lesion_type"].isin(lesion_types)]
        _record(detail)
        return balanced_macro_auroc(detail, features=features, group_cols=group_cols)

    def _concat_leg(per_cohort_fn, label):
        per_cohort = {name: per_cohort_fn(name, root)
                      for name, root in output_dirs_by_score_name.items()}
        nonempty = {name: p for name, p in per_cohort.items()
                    if p is not None and not p.empty}
        # Warn when ANY expected cohort contributed zero rows -- not only when ALL
        # did. A silently-dropped cohort collapses the "balanced across cohorts"
        # objective to a SUBSET (e.g. MICs-only), misleading the greedy lock.
        missing = [name for name in output_dirs_by_score_name if name not in nonempty]
        if missing:
            warnings.warn(
                f"{label} objective: cohort(s) {missing} produced NO scorable rows -- "
                "the balanced-across-cohorts objective is collapsing to a SUBSET "
                "(GT/layout mismatch, missing held-out-control maps, or a single-disease "
                "cohort for disease_control). Check the data/config.",
                RuntimeWarning, stacklevel=2,
            )
        detail = pd.concat(list(nonempty.values()), ignore_index=True) if nonempty else pd.DataFrame()
        _record(detail)
        return balanced_macro_auroc(detail, features=features, group_cols=group_cols)

    def _control_leg():
        # Single-shot specificity: patient lesion |z| (positives) vs the |z| of THIS
        # evaluation's control negatives at the SAME region, region-matched. The
        # negatives are one reference fold's held-out controls (reference-CV) or the
        # pooled set (legacy). Reference-CV averages this across folds.
        return _concat_leg(
            lambda name, root: per_patient_vs_control_auc(
                root, name, control_root=root,
                control_subjects=_flatten_control_pairs(control_subjects_by_name.get(name)),
                method=config["method"], features=features),
            "control-mode")

    def _disease_control_leg():
        return _concat_leg(
            lambda name, root: per_patient_vs_other_disease_auc(
                root, name, method=config["method"], features=features),
            "disease-control-mode")

    have_controls = bool(control_subjects_by_name)

    if objective_mode == "within":
        return _within_leg()
    if objective_mode == "within_tle":
        return _within_leg({"TLE"})
    if objective_mode == "within_fcd":
        return _within_leg({"FCD"})
    if objective_mode == "disease_control":
        return _disease_control_leg()

    if objective_mode == "control":
        if not have_controls:
            warnings.warn(
                "objective_mode='control' but no held-out controls were scored; "
                "falling back to the within-subject leg for this evaluation.",
                RuntimeWarning, stacklevel=2,
            )
            return _within_leg()
        return _control_leg()

    # objective_mode == "both": ALIAS of within -- within-subject localization over
    # BOTH TLE and FCD lesions together (identical to the plain 'within' leg).
    return _within_leg()


def reference_cv_score_fn(config, fold_output_dirs, *, control_subjects_by_name=None,
                          objective_mode="both", **score_kwargs):
    """Reference cross-validated objective: score EVERY mode once per reference fold
    and lock on the MEAN (mean +/- SD stashed for the history CSV).

    ``fold_output_dirs`` is ``{score_name: {fold_id: analysis_root}}`` where each
    fold's root holds |z| maps computed against reference_k = controls minus fold_k
    (patients always; held-out controls_k too when the mode needs them).
    ``control_subjects_by_name`` is ``{score_name: {fold_id: [held-out pairs]}}``.
    For each fold the single-shot :func:`default_score_fn` is evaluated over that
    fold's roots (its held-out controls as a flat negative list), giving one value
    per fold. Unlike the old design, the PATIENT positives are re-scored per fold
    too, so within/disease/control all vary over the SAME K reference resamplings --
    a uniform mean +/- SD for every objective (the point of reference-CV)."""
    import numpy as np
    global _LAST_CONTROL_FOLD_STATS
    _LAST_CONTROL_FOLD_STATS = None

    fold_ids = sorted({fid for roots in fold_output_dirs.values() for fid in roots})
    per_fold = []
    per_fold_feats = []   # each fold's {feature: balanced AUROC} breakdown
    for k in fold_ids:
        odk = {name: roots[k] for name, roots in fold_output_dirs.items() if k in roots}
        if not odk:
            continue
        csk = None
        if control_subjects_by_name:
            csk = {name: control_subjects_by_name[name][k]
                   for name in odk
                   if name in control_subjects_by_name and k in control_subjects_by_name[name]}
        val = default_score_fn(config, odk, control_subjects_by_name=csk,
                               objective_mode=objective_mode, **score_kwargs)
        if pd.notna(val):
            per_fold.append(float(val))
            per_fold_feats.append(dict(_LAST_FEATURE_AUROC))   # the fold's per-feature split
    if not per_fold:
        return float("nan")
    arr = np.asarray(per_fold, dtype=float)
    mean = float(arr.mean())
    std = float(arr.std(ddof=1)) if arr.size > 1 else 0.0   # sample SD across folds

    # Per-feature: average each feature's balanced AUROC across the folds that scored
    # it (mean + SD), so the history CSV decomposes the overall mean by feature.
    feats = sorted({f for d in per_fold_feats for f in d})
    per_feature, per_feature_std = {}, {}
    for f in feats:
        vals = np.asarray([d[f] for d in per_fold_feats if f in d], dtype=float)
        per_feature[f] = float(vals.mean())
        per_feature_std[f] = float(vals.std(ddof=1)) if vals.size > 1 else 0.0

    _LAST_CONTROL_FOLD_STATS = {
        "aurocs": [round(v, 6) for v in per_fold],
        "mean": mean, "std": std, "n_folds": int(arr.size),
        "per_feature": {f: round(v, 6) for f, v in per_feature.items()},
        "per_feature_std": {f: round(v, 6) for f, v in per_feature_std.items()},
    }
    feat_str = "  ".join(f"{f}={per_feature[f]:.3f}" for f in feats)
    print(f"    [objective {objective_mode}] per-fold AUROC={['%.4f' % v for v in per_fold]} "
          f"-> mean={mean:.4f} SD={std:.4f} (n_folds={arr.size})"
          + (f"  |  per-feature: {feat_str}" if feat_str else ""))
    return mean


def single_shot_score_fn(config, output_dirs_by_score_name, *, control_subjects_by_name=None,
                         objective_mode="both", **score_kwargs):
    """Single-shot objective (no reference-CV): score ONCE over the given maps and
    record it as a 1-"fold" stats entry (mean=value, SD=0, n_folds=1) so the history
    CSV carries the SAME overall + per-feature columns as the reference-CV path.
    Used for the within_* objectives under cv_scope='control_disease'."""
    global _LAST_CONTROL_FOLD_STATS
    _LAST_CONTROL_FOLD_STATS = None
    val = default_score_fn(config, output_dirs_by_score_name,
                           control_subjects_by_name=control_subjects_by_name,
                           objective_mode=objective_mode, **score_kwargs)
    if pd.isna(val):
        return float("nan")
    feats = dict(_LAST_FEATURE_AUROC)   # the leg's per-feature breakdown (this call)
    _LAST_CONTROL_FOLD_STATS = {
        "aurocs": [round(float(val), 6)],
        "mean": float(val), "std": 0.0, "n_folds": 1,
        "per_feature": {f: round(v, 6) for f, v in feats.items()},
        "per_feature_std": {f: 0.0 for f in feats},
    }
    return float(val)


# ---------------------------------------------------------------------------
# Greedy driver
# ---------------------------------------------------------------------------
_BOOL_HISTORY_KEYS = {
    "predictive_wscore",
    "use_curvature_covariates",
    "control_correlation_filter",
}
_INT_HISTORY_KEYS = {
    "wscore_surface_smoothing_iterations",
    "cortical_smoothing",
    "hippocampal_smoothing",
}
_FLOAT_HISTORY_KEYS = {
    "prediction_variance_percentile",
    "control_correlation_threshold",
    "control_correlation_quantile",
    "qc_alpha",
}


def _is_missing_history_value(value):
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _history_bool(value):
    if _is_missing_history_value(value):
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _history_optional_float(value):
    if _is_missing_history_value(value):
        return None
    text = str(value).strip()
    if text.lower() in {"", "none", "nan", "null", "<na>"}:
        return None
    return float(value)


def _history_config(row):
    """Rebuild a validated config from one CSV history row.

    New history rows store ``config_json`` so the resume path is exact. The
    column fallback keeps manually edited/simple histories readable where
    possible, but automatic resume only uses versioned histories.
    """
    raw_json = row.get("config_json")
    if not _is_missing_history_value(raw_json) and str(raw_json).strip():
        config = json.loads(raw_json)
    else:
        config = dict(SIMPLEST_BASELINE)
        for key in (
            *_NON_AXIS_KEYS,
            "normalization",
            "wscore_preprocessing",
            "intensity_depth_model",
            "blur_depth_model",
            "method",
            "wscore_covariate_model",
            "wscore_distribution",
            "predictive_wscore",
            "use_curvature_covariates",
            "wscore_surface_smoothing_iterations",
            "control_correlation_filter",
        ):
            if key in row and not _is_missing_history_value(row[key]):
                config[key] = row[key]

    for key in _BOOL_HISTORY_KEYS:
        if key in config:
            config[key] = _history_bool(config[key])
    for key in _INT_HISTORY_KEYS:
        if key in config and not _is_missing_history_value(config[key]):
            config[key] = int(float(config[key]))
    for key in _FLOAT_HISTORY_KEYS:
        if key in config:
            config[key] = _history_optional_float(config[key])
    return validate(config)[0]


def _history_score(row):
    value = row.get("auroc")
    if _is_missing_history_value(value):
        return float("nan")
    return float(value)


def _history_error(row):
    value = row.get("error")
    if _is_missing_history_value(value):
        return None
    text = str(value).strip()
    return text or None


def _history_path_obj(history_path):
    return None if history_path is None else Path(history_path)


def _write_history(history, history_path):
    path = _history_path_obj(history_path)
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    pd.DataFrame(history).to_csv(tmp, index=False)
    os.replace(tmp, path)


def _load_resume_history(history_path, stages, *, verbose=True):
    path = _history_path_obj(history_path)
    if path is None or not path.exists():
        return None
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return None
    if df.empty:
        return None
    if "history_schema_version" not in df.columns:
        if verbose:
            print(f"Found {path}, but it predates resumable history; starting fresh.")
        return None
    versions = {int(float(v)) for v in df["history_schema_version"].dropna().unique()}
    if versions != {HISTORY_SCHEMA_VERSION}:
        if verbose:
            print(
                f"Found {path}, but its resume schema {sorted(versions)} does not "
                f"match {HISTORY_SCHEMA_VERSION}; starting fresh."
            )
        return None

    history = df.to_dict("records")
    winner_idxs = [i for i, row in enumerate(history) if _history_bool(row.get("is_winner"))]
    if not winner_idxs:
        return None

    last_winner_idx = winner_idxs[-1]
    last_winner = history[last_winner_idx]
    current = _history_config(last_winner)
    baseline_score = _history_score(last_winner)
    stage_names = [stage.name for stage in stages]
    last_stage = str(last_winner.get("stage", ""))
    if last_stage == "baseline":
        start_stage_index = 0
    elif last_stage in stage_names:
        start_stage_index = stage_names.index(last_stage) + 1
    else:
        if verbose:
            print(f"Found {path}, but last winner stage {last_stage!r} is not in this run; starting fresh.")
        return None

    if verbose:
        if start_stage_index >= len(stages):
            print(f"Resuming from {path}: all stages already completed.")
        else:
            print(
                f"Resuming from {path}: last locked stage is {last_stage}; "
                f"continuing at {stages[start_stage_index].name}."
            )
    return {
        "history": history,
        "current": current,
        "baseline_score": baseline_score,
        "start_stage_index": start_stage_index,
    }


def _build_site_harmonizer(config, cohorts):
    """Construct the cross-site harmonizer for a candidate (fit lazily on the
    pooled controls exported in pass 1). Covariates = the control cohort's
    normative columns (age/sex); reference site anchors to MICs when present."""
    from zbrains.harmonization import SiteHarmonizer
    cov = cohorts[0].control_dataset.demographics.normative_columns or []
    names = [c.score_name for c in cohorts]
    ref = "MICs" if "MICs" in names else names[0]
    return SiteHarmonizer(
        config["site_harmonization"], cov, ref_site=ref,
        site_covariate=config.get("scoring_site_covariate", False),
    )


def analyze_config_across_cohorts(config, cohorts, env, control_correlation_threshold, *,
                                  needs_controls, built_bases=None,
                                  reprocess_controls=False, verbose=False,
                                  reference_cv=False):
    """Analyze ONE config across all cohorts.

    Builds the exclusion-keyed base lazily (skipped if already in ``built_bases``)
    and runs the two-pass site harmonization when configured. The return shape
    depends on ``reference_cv``:

    * ``reference_cv=False`` (legacy): ``(output_dirs_by_score_name,
      control_subjects_by_name)`` -- patients scored ONCE vs the full control
      reference, and (when ``needs_controls``) held-out-control folds scored into
      the same dir; ``control_subjects_by_name[name]`` is a ``{fold: pairs}`` map.
      Used by mock/custom-score_fn callers.
    * ``reference_cv=True``: ``(fold_output_dirs, control_subjects_by_name)`` where
      ``fold_output_dirs[name] = {fold: root}`` -- patients (and, when
      ``needs_controls``, held-out controls) re-scored against EACH reference fold,
      for :func:`reference_cv_score_fn`. Every objective then reports mean +/- SD
      over the same K reference resamplings.

    Shared by the greedy ``evaluate()`` and by ``cross_evaluate_configs``.
    """
    if built_bases is None:
        built_bases = set()
    normalization = config["normalization"]
    ctx = []
    for cohort in cohorts:
        excluded = _resolve_excluded_pairs(cohort, config, env, verbose=verbose)
        sig = _exclusion_signature(excluded)
        control_ds = _restricted_control_dataset(cohort, excluded)
        # Cache the built base by its BASE signature (exclusion dropped when
        # dataset_norm=='none'), so exclusion arms that share a base don't rebuild it.
        base_key = (cohort.name, normalization, _base_signature(config, sig))
        if base_key not in built_bases:
            _ensure_base_processed(
                cohort, config, env,
                control_dataset=control_ds, exclusion_signature=sig,
                verbose=verbose, reprocess_controls=reprocess_controls,
            )
            built_bases.add(base_key)
        ctx.append((cohort, control_ds, sig))

    if reference_cv:
        # Reference-CV: re-score patients (and held-out controls when needed) against
        # each reference fold, PER-FOLD refitting the dataset-norm (RAVEL/Nyul) AND
        # the site harmonizer on the train fold only -- no second-order leakage from
        # the held-out fold. Every objective reports mean +/- SD over the folds.
        return _analyze_reference_cv(
            [(c, cds, s) for c, cds, s in ctx if cds is not None],
            config, env, control_correlation_threshold,
            needs_controls=needs_controls, verbose=verbose,
        )

    # Legacy (non-reference-CV) path: cross-site harmonizer fit ONCE on the pooled
    # full control set (two-pass), then patients scored against it.
    harmonizer = None
    if config.get("site_harmonization", "none") != "none":
        harmonizer = _build_site_harmonizer(config, cohorts)
        for cohort, control_ds, sig in ctx:
            export = {}
            _analyze_candidate(
                cohort, config, env, control_correlation_threshold,
                control_dataset=control_ds, exclusion_signature=sig,
                export_control_features=export, verbose=verbose,
            )
            harmonizer.add_site_features(cohort.score_name, export)

    output_dirs = {}
    control_subjects_by_name = {}
    for cohort, control_ds, sig in ctx:
        output_dirs[cohort.score_name] = _analyze_candidate(
            cohort, config, env, control_correlation_threshold,
            control_dataset=control_ds, exclusion_signature=sig,
            site_harmonizer=harmonizer, verbose=verbose,
        )
        if needs_controls and control_ds is not None:
            # The returned fold map ({fold_id: [held-out pairs]}) lets the objective
            # compute a PER-FOLD AUROC (mean +/- SD) rather than one pooled value.
            control_subjects_by_name[cohort.score_name] = _score_heldout_controls(
                cohort, config, env, control_correlation_threshold,
                control_ds, sig, site_harmonizer=harmonizer, verbose=verbose,
            )
    return output_dirs, control_subjects_by_name


def run_staged_optimization(
    *,
    cohorts,
    env,
    stages=DEFAULT_STAGES,
    start_config=None,
    score_fn=None,
    objective_mode="both",
    control_correlation_threshold=zb.DEFAULT_CONTROL_CORRELATION_THRESHOLD,
    reprocess_controls=False,
    verbose=True,
    history_path=None,
    resume=True,
    fail_fast=True,
    cv_scope="control_disease",
):
    """Run the sequential greedy optimization and return (best_config, history).

    Parameters
    ----------
    cohorts : list[Cohort]
        Datasets pooled into the single objective (e.g. MICs and NOEL).
    env : zbenv
        Processing environment.
    stages : list[Stage]
        Ordered greedy stages. Order is significant (unlike OFAT families).
    start_config : dict, optional
        Starting point (defaults to a copy of SIMPLEST_BASELINE).
    score_fn : callable, optional
        ``score_fn(config, output_dirs_by_score_name) -> float``. Defaults to
        :func:`default_score_fn` bound to ``objective_mode``. A custom / mock
        score_fn makes the driver testable without heavy analysis.
    objective_mode : {"within","control","both"}, default "both"
        Which comparison drives the greedy: ``within`` = lesion vs the patient's
        OWN non-lesion vertices; ``control`` = lesion vs held-out-control vertices
        at the SAME region; ``both`` = ALIAS of ``within``. Held-out-control
        scoring (the extra K-fold analyze passes) runs ONLY in ``control`` mode.
        Ignored when a custom ``score_fn`` is supplied.
    history_path : str | Path, optional
        If given, the per-candidate history is checkpointed here as CSV after
        every baseline/candidate/winner row.
    resume : bool, default True
        If ``history_path`` points to a compatible resumable CSV, continue from
        the last locked winner. If the run stopped mid-stage, already recorded
        candidates in that stage are skipped and the remaining candidates are run.
        Candidates recorded as ERRORS are NOT treated as done -- they are retried
        on resume (so a fixed crash re-runs instead of being skipped).
    fail_fast : bool, default True
        If True, a candidate whose evaluation raises (e.g. a RAVEL/ANTs worker
        SIGSEGV, an analysis crash) STOPS the whole run after recording the error,
        rather than silently locking ``keep``/the incumbent. Set False for the old
        best-effort behaviour (skip the failed candidate and continue).
    cv_scope : {"control_disease","full","off"}, default "control_disease"
        Which objectives use reference cross-validation (re-score patients against
        each control-fold-complement reference -> per-fold mean +/- SD). Default
        reference-CVs ONLY the cross-subject specificity tests (control,
        disease_control, both); within_* stay single-shot (cheap, ~1 analyze pass --
        their fold variance is near-zero anyway). "full" reference-CVs every
        objective (~5x cost on within_*); "off" makes all single-shot. Ignored when
        a custom ``score_fn`` is supplied. Every mode still records the overall +
        per-feature history columns (single-shot -> n_folds=1, SD=0).
    """
    if objective_mode not in OBJECTIVE_MODES:
        raise ValueError(f"objective_mode must be one of {OBJECTIVE_MODES}, got {objective_mode!r}")
    if cv_scope not in CV_SCOPES:
        raise ValueError(f"cv_scope must be one of {CV_SCOPES}, got {cv_scope!r}")
    if history_path is not None:
        # Resolve to an ABSOLUTE path now (cwd is correct here); the processing
        # pipeline changes the working directory mid-run, so a relative history
        # path would otherwise scatter the CSV into whatever derivative dir the
        # cwd last landed in instead of the intended location.
        history_path = os.path.abspath(history_path)
    # Reference-CV only drives the DEFAULT objective (a custom/mock score_fn keeps
    # the legacy single-root shape it expects), and only for the modes selected by
    # cv_scope (default: control/disease only; within_* stay single-shot).
    reference_cv = score_fn is None and _mode_uses_reference_cv(objective_mode, cv_scope)
    # Bind the mode into the objective (functools.partial preserves the signature,
    # so evaluate()'s control_subjects_by_name introspection still works). Reference-
    # CV scores once per reference fold and averages; single-shot scores once (both
    # record the same overall + per-feature history columns).
    if score_fn is None:
        objective_fn = reference_cv_score_fn if reference_cv else single_shot_score_fn
        score_fn = functools.partial(objective_fn, objective_mode=objective_mode)
        if verbose:
            print(f"[greedy] cv_scope={cv_scope!r} -> objective_mode={objective_mode!r} "
                  f"uses {'reference-CV (per-fold mean+/-SD)' if reference_cv else 'single-shot'}")
    # Held-out controls are only needed by the control / both legs (within_*
    # and disease_control use patient maps only). Reference-CV re-scores PATIENTS per
    # fold for every mode regardless; this flag only gates scoring held-out controls.
    needs_controls = objective_mode in _CONTROL_MODES
    current = dict(SIMPLEST_BASELINE)
    if start_config is not None:
        current.update(start_config)
    current, warnings = validate(current)
    for w in warnings:
        print(f"  NOTE (baseline): {w}")
    _MODE_DESC = {
        "within": "lesion vs the patient's OWN non-lesion vertices (all disease)",
        "within_tle": "lesion vs own non-lesion vertices (TLE only)",
        "within_fcd": "lesion vs own non-lesion vertices (FCD only)",
        "control": "lesion vs ALL held-out-control vertices",
        "disease_control": "lesion vs the OTHER disease's patients' non-lesional vertices",
        "both": "alias of within (within-subject, TLE+FCD together)",
    }
    print(f"[greedy] objective_mode = {objective_mode!r} ({_MODE_DESC.get(objective_mode, '?')})")

    score_cache = {}
    fold_stats_cache = {}  # config_key -> per-fold control AUROC stats (or None)
    built_bases = set()  # (cohort.name, normalization, processing_sig) already processed

    def evaluate(config):
        """Analyze + score a config across all cohorts (cached by identity key).

        The processed base for this config's normalization is built lazily the
        first time it is needed, so sweeping normalization (none/whitestripe)
        only pays for the bases actually visited.
        """
        key = _config_key(config, control_correlation_threshold)
        if key in score_cache:
            return score_cache[key]
        output_dirs, control_subjects_by_name = analyze_config_across_cohorts(
            config, cohorts, env, control_correlation_threshold,
            needs_controls=needs_controls, built_bases=built_bases,
            reprocess_controls=reprocess_controls, verbose=verbose,
            reference_cv=reference_cv,
        )
        # Pass the held-out-control subjects to the objective only if it accepts
        # them (the real objective does; mock 2-arg score_fns do not).
        score_kwargs = {}
        if control_subjects_by_name and "control_subjects_by_name" in inspect.signature(score_fn).parameters:
            score_kwargs["control_subjects_by_name"] = control_subjects_by_name
        _reset_fold_stats()   # clear any stale stats so a mock score_fn records None
        score = score_fn(config, output_dirs, **score_kwargs)
        # Capture the per-fold stats the objective just stashed (None for a mock
        # score_fn) so they can be recorded to history.
        fold_stats_cache[key] = _LAST_CONTROL_FOLD_STATS
        score_cache[key] = (score, output_dirs)
        return score_cache[key]

    resume_state = _load_resume_history(history_path, stages, verbose=verbose) if resume else None
    if resume_state is not None:
        history = resume_state["history"]
        current = resume_state["current"]
        baseline_score = resume_state["baseline_score"]
        start_stage_index = resume_state["start_stage_index"]
    else:
        # Score and immediately checkpoint the starting baseline.
        baseline_score, baseline_dirs = evaluate(current)
        if not pd.notna(baseline_score):
            raise RuntimeError(
                "Baseline objective is NaN -- no scorable (lesion, feature) cells were "
                "produced. This usually means the discovered lesion/GT layout does not "
                "match the data, or (specificity mode) held-out-control maps were not "
                "written. Refusing to optimize: a NaN baseline would make every "
                "`score > best_score` comparison False and silently lock "
                "SIMPLEST_BASELINE. Fix the data/config first."
            )
        print(f"\n[stage 0] SIMPLEST_BASELINE  ->  pooled macro AUROC = {baseline_score:.4f}")
        history = [_history_row("baseline", "baseline", current, baseline_score,
                                control_correlation_threshold, is_winner=True,
                                fold_stats=fold_stats_cache.get(
                                    _config_key(current, control_correlation_threshold)))]
        _write_history(history, history_path)
        start_stage_index = 0

    for stage in stages[start_stage_index:]:
        # Conditional stages (e.g. gp_uncertainty) are skipped when their
        # precondition on the locked config is not met. Carry the config forward
        # unchanged; record a winner row (once) so resume advances past the stage.
        if stage.condition is not None and not stage.condition(current):
            already = any(row.get("stage") == stage.name and _history_bool(row.get("is_winner"))
                          for row in history)
            print(f"\n[stage {stage.name}] skipped: precondition not met for the locked config")
            if not already:
                history.append(_history_row(
                    stage.name, "keep(skipped:condition-not-met)", current,
                    baseline_score, control_correlation_threshold, is_winner=True))
                _write_history(history, history_path)
            continue

        print(f"\n{'=' * 78}\n[stage] {stage.name}  (current AUROC = {baseline_score:.4f})\n{'=' * 78}")
        best_label, best_config, best_score = "keep", dict(current), baseline_score
        completed_labels = set()

        # A restart can land in the middle of a stage. Treat any row already
        # recorded for this stage as done, and recover the best-so-far from the
        # successful rows so the stage can finish without repeating work.
        for row in history:
            if row.get("stage") != stage.name:
                continue
            label = str(row.get("candidate", ""))
            if label.startswith("WINNER:"):
                continue
            if _history_error(row):
                continue    # an errored candidate is NOT done -> retry it on resume
            completed_labels.add(label)
            score = _history_score(row)
            if not pd.notna(score):
                continue
            candidate_config = _history_config(row)
            if score > best_score:
                best_label, best_config, best_score = label, candidate_config, score

        for label, overrides in stage.candidates:
            if label in completed_labels:
                print(f"  skip {label}: already recorded in {history_path}")
                continue
            candidate = dict(current)
            candidate.update(overrides)
            try:
                candidate, warns = validate(candidate)
            except ValueError as exc:
                print(f"  skip {label}: invalid combination ({exc})")
                continue
            for w in warns:
                print(f"  NOTE ({label}): {w}")
            cand_key = _config_key(candidate, control_correlation_threshold)
            try:
                score, _dirs = evaluate(candidate)
            except Exception as exc:
                print(f"  {label}: FAILED ({type(exc).__name__}: {exc})")
                history.append(_history_row(stage.name, label, candidate, float("nan"),
                                            control_correlation_threshold, error=str(exc)))
                _write_history(history, history_path)
                if fail_fast:
                    # A real evaluation failure (e.g. a RAVEL/ANTs worker SIGSEGV)
                    # must STOP the run -- otherwise the stage silently locks
                    # `keep`/the incumbent as if the failed arm were simply worse,
                    # corrupting the optimization. The error row is checkpointed
                    # above and is RETRIED on resume once the cause is fixed.
                    raise RuntimeError(
                        f"[stage {stage.name}] candidate {label!r} FAILED "
                        f"({type(exc).__name__}: {exc}). fail_fast=True -> stopping the "
                        f"run rather than locking a winner without this arm. Fix the "
                        f"cause then resume (this candidate re-runs); or pass "
                        f"fail_fast=False to skip failed candidates and continue."
                    ) from exc
                continue
            marker = "  <- new winner" if score > best_score else ""
            print(f"  {label}: AUROC = {score:.4f}{marker}")
            history.append(_history_row(stage.name, label, candidate, score,
                                        control_correlation_threshold,
                                        fold_stats=fold_stats_cache.get(cand_key)))
            _write_history(history, history_path)
            if score > best_score:
                best_label, best_config, best_score = label, candidate, score

        current, baseline_score = best_config, best_score
        # mark the locked winner in history
        history.append(_history_row(stage.name, f"WINNER:{best_label}", current,
                                    best_score, control_correlation_threshold, is_winner=True,
                                    fold_stats=fold_stats_cache.get(
                                        _config_key(current, control_correlation_threshold))))
        _write_history(history, history_path)
        print(f"[stage {stage.name}] locked: {best_label}  ->  AUROC = {best_score:.4f}")

    print(f"\n{'#' * 78}\nFINAL pipeline  ->  pooled macro AUROC = {baseline_score:.4f}")
    print(json.dumps(current, indent=2, default=str))

    history_df = pd.DataFrame(history)
    if history_path is not None:
        _write_history(history, history_path)
        print(f"Wrote greedy history: {history_path}")
    return current, history_df


def _hashable_threshold(threshold):
    """Make a threshold hashable: a per-feature mapping -> sorted tuple."""
    if hasattr(threshold, "items"):
        return tuple(sorted((str(k), v) for k, v in threshold.items()))
    return threshold


def _config_key(config, threshold):
    axis = _axis_config(config)
    return (
        tuple(sorted(axis.items())),
        _hashable_threshold(config.get("control_correlation_threshold", threshold)),
        config["method"],
        config.get("prediction_variance_percentile"),
        _hashable_threshold(threshold),
        # Step 8 outlier drivers: deterministically identify the exclusion set (the
        # per-cohort excluded pairs derive from these + normalization, which is in
        # `axis`). Baseline (none/None) hashes identically to a pre-Step-8 config.
        config.get("outlier_method", "none"),
        _hashable_threshold(config.get("control_correlation_quantile")),
        config.get("qc_alpha"),
        config.get("cortical_smoothing", SIMPLEST_BASELINE["cortical_smoothing"]),
        config.get("hippocampal_smoothing", SIMPLEST_BASELINE["hippocampal_smoothing"]),
        # Cross-site harmonization (analysis-level). none/False hash identically to
        # a pre-harmonization config (default-get equivalence).
        config.get("site_harmonization", "none"),
        bool(config.get("scoring_site_covariate", False)),
    )


def _history_row(stage, label, config, score, threshold, *, is_winner=False, error=None,
                 fold_stats=None):
    row = {
        "history_schema_version": HISTORY_SCHEMA_VERSION,
        "stage": stage,
        "candidate": label,
        "auroc": score,
        "is_winner": is_winner,
        "config_json": json.dumps(config, sort_keys=True, default=str),
    }
    # Reference-CV per-fold AUROC (mean +/- SD). `auroc` above already IS the mean
    # (what the greedy locks on); these columns expose the fold spread so a
    # narrow-margin winner with high variance is visible in the history CSV.
    if fold_stats:
        row["control_fold_mean"] = fold_stats.get("mean")
        row["control_fold_std"] = fold_stats.get("std")
        row["control_n_folds"] = fold_stats.get("n_folds")
        row["control_fold_aurocs"] = json.dumps(fold_stats.get("aurocs"))
        # Per-FEATURE decomposition of the overall mean (one column per feature,
        # `auroc_<feature>`, + its across-fold SD `auroc_<feature>_std`) so the CSV
        # shows which features carry the score in addition to the balanced mean.
        for feat, val in (fold_stats.get("per_feature") or {}).items():
            row[f"auroc_{feat}"] = val
        for feat, val in (fold_stats.get("per_feature_std") or {}).items():
            row[f"auroc_{feat}_std"] = val
    row.update({k: config.get(k) for k in (
        "normalization", "wscore_preprocessing", "intensity_depth_model",
        "blur_depth_model", "method", "wscore_covariate_model",
        "wscore_distribution", "predictive_wscore", "use_curvature_covariates",
        "wscore_surface_smoothing_iterations", "control_correlation_filter",
        "control_correlation_quantile",
        "prediction_variance_percentile",
        "sample_surface", "t1w_flair_self_norm", "quant_sample_surface",
        "subject_norm", "dataset_norm", "outlier_method", "qc_alpha",
        "cortical_smoothing", "hippocampal_smoothing",
        "site_harmonization", "scoring_site_covariate",
    )})
    row["control_correlation_threshold"] = config.get("control_correlation_threshold", threshold)
    if error is not None:
        row["error"] = error
    return row


# ---------------------------------------------------------------------------
# Cross-evaluation: score each winning pipeline under EVERY objective metric
# ---------------------------------------------------------------------------
# The five objectives whose winners are cross-evaluated against one another.
CROSS_EVAL_MODES = ("within_tle", "within_fcd", "within", "control", "disease_control")


def cross_evaluate_configs(configs_by_label, cohorts, env, *,
                           control_correlation_threshold=zb.DEFAULT_CONTROL_CORRELATION_THRESHOLD,
                           modes=CROSS_EVAL_MODES, reprocess_controls=False, verbose=False,
                           reference_cv=True):
    """Score every winning pipeline under EVERY objective metric.

    ``configs_by_label`` maps a label (typically the objective a pipeline WON
    under) to its (already validated) config. Each config is analyzed ONCE across
    the cohorts -- under reference-CV (default) patients + held-out controls are
    re-scored against each reference fold -- then scored under each mode in
    ``modes`` (the mean over reference folds when ``reference_cv``). Returns a
    DataFrame indexed by label, columns = ``modes`` + ``mean`` + ``std`` (the
    mean/SD of each pipeline's AUROC ACROSS the metrics -- a robustness score).
    """
    built_bases = set()
    score_one = reference_cv_score_fn if reference_cv else default_score_fn
    rows = {}
    for label, config in configs_by_label.items():
        if verbose:
            print(f"\n[cross-eval] analyzing winner {label!r} under all {len(modes)} metrics")
        output_dirs, control_subjects = analyze_config_across_cohorts(
            config, cohorts, env, control_correlation_threshold,
            needs_controls=True, built_bases=built_bases,
            reprocess_controls=reprocess_controls, verbose=verbose,
            reference_cv=reference_cv,
        )
        scores = {}
        for mode in modes:
            scores[mode] = score_one(
                config, output_dirs,
                control_subjects_by_name=control_subjects, objective_mode=mode,
            )
        rows[label] = scores
        if verbose:
            print("  " + "  ".join(f"{m}={scores[m]:.4f}" for m in modes))

    matrix = pd.DataFrame(rows).T[list(modes)]
    matrix["mean"] = matrix[list(modes)].mean(axis=1)
    matrix["std"] = matrix[list(modes)].std(axis=1)
    return matrix


def plot_cross_evaluation(matrix, output_path, *, modes=CROSS_EVAL_MODES):
    """Render the cross-evaluation: a winners x metrics AUROC heatmap plus a
    mean+/-SD robustness bar chart with the overall winner highlighted. Returns
    the winning label (highest mean AUROC across the metrics)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    modes = list(modes)
    labels = list(matrix.index)
    M = matrix[modes].to_numpy(dtype=float)
    means = matrix["mean"].to_numpy(dtype=float)
    stds = np.nan_to_num(matrix["std"].to_numpy(dtype=float), nan=0.0)

    order = list(np.argsort(-np.nan_to_num(means, nan=-1.0)))
    winner = labels[order[0]]

    fig, (ax0, ax1) = plt.subplots(
        1, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1.35, 1]})

    im = ax0.imshow(M, cmap="viridis", vmin=0.5, vmax=1.0, aspect="auto")
    ax0.set_xticks(range(len(modes))); ax0.set_xticklabels(modes, rotation=45, ha="right")
    ax0.set_yticks(range(len(labels))); ax0.set_yticklabels(labels)
    ax0.set_xlabel("evaluated under metric"); ax0.set_ylabel("pipeline won under")
    ax0.set_title("Cross-evaluation AUROC (winners x metrics)")
    for i in range(len(labels)):
        for j in range(len(modes)):
            v = M[i, j]
            txt = "n/a" if not np.isfinite(v) else f"{v:.3f}"
            on_diag = labels[i] == modes[j]
            ax0.text(j, i, txt, ha="center", va="center", fontsize=8,
                     fontweight="bold" if on_diag else "normal",
                     color="white" if (np.isfinite(v) and v < 0.8) else "black")
    fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04, label="AUROC")

    slabels = [labels[k] for k in order]
    smeans = [means[k] for k in order]
    sstds = [stds[k] for k in order]
    colors = ["#2563eb" if lab == winner else "#94a3b8" for lab in slabels]
    ypos = list(range(len(slabels)))
    ax1.barh(ypos, smeans, xerr=sstds, color=colors, capsize=4)
    ax1.set_yticks(ypos); ax1.set_yticklabels(slabels); ax1.invert_yaxis()
    ax1.set_xlim(0.5, 1.0)
    ax1.set_xlabel("mean AUROC across the metrics (+/- SD)")
    ax1.set_title(f"Robustness -- overall winner: {winner}")
    for y, k in enumerate(order):
        if np.isfinite(means[k]):
            ax1.text(min(means[k] + sstds[k] + 0.01, 0.98), y, f"{means[k]:.3f}",
                     va="center", fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return winner


if __name__ == "__main__":
    print(__doc__)
    print("\nSIMPLEST_BASELINE:")
    print(json.dumps(SIMPLEST_BASELINE, indent=2, default=str))
    print("\nDEFAULT_STAGES (greedy order):")
    for i, stage in enumerate(DEFAULT_STAGES, 1):
        print(f"  {i}. {stage.name}")
        for label, overrides in stage.candidates:
            print(f"       - {label}: {overrides}")
