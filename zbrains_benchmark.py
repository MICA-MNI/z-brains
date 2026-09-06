"""Shared z-brains benchmarking driver.

This module is the single, authoritative description of the benchmark's
engineering-decision space. The example scripts (``examplenew.py``,
``examplenew_noel.py``) only supply dataset-specific configuration and then
call :func:`run_benchmark`; everything about *what* the sweep varies and *how*
those choices may be combined lives here.

Mental model
------------
* Every knob is an **AXIS** (see :data:`AXES`). The values *within* one axis are
  **mutually exclusive** -- a single run picks exactly one of them.
* Different axes are **independent / freely combinable** by default: any value
  of one axis may be paired with any value of another...
* ...*except* for the explicit cross-axis rules in
  :data:`MUTUALLY_EXCLUSIVE` (combinations the analysis code forbids) and
  :data:`REDUNDANT` (combinations that are allowed but silently ignore one
  setting). All of these are driven by whether ``wscore_covariate_model`` is a
  Gaussian-process kernel.

The sweep itself is **one-factor-at-a-time (OFAT)**: every generated variant
differs from the baseline in exactly one axis. Because OFAT never combines two
non-baseline choices, it cannot by itself hit a mutually-exclusive pair -- but
:func:`validate_configuration` still enforces every rule so that hand-edited
baselines, restricted sweeps, or future combined sweeps fail fast here instead
of deep inside ``analysis.py``.

Run ``python zbrains_benchmark.py`` to print the full decision space.
"""

from __future__ import annotations

import hashlib
import os
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from zbrains.hippunfold import hippunfold_cache_tag


# ============================================================================
# 1. ENGINEERING-DECISION AXES
#    Values within an axis are MUTUALLY EXCLUSIVE (pick exactly one per run).
# ============================================================================
@dataclass(frozen=True)
class Axis:
    """One engineering decision.

    Attributes
    ----------
    name : str
        Keyword accepted by ``dataset.process`` / ``dataset.analyze``.
    values : tuple
        Every valid choice. Mutually exclusive: a run selects one.
    baseline : object
        The value used by the OFAT baseline configuration.
    stage : {"processing", "analysis"}
        ``processing`` axes change the reusable processed base directory and
        force a reprocess; ``analysis`` axes only re-score an existing base.
    family : str
        OFAT grouping used by ``enabled_families`` (see :data:`FAMILIES`).
    scope : str
        Human-readable note on which part of the pipeline the axis affects.
    sweep_values : tuple, optional
        Subset of ``values`` actually iterated during the OFAT sweep. Defaults
        to ``values``. Used to keep an axis valid-but-dormant (e.g. ``ravel``).
    """

    name: str
    values: tuple
    baseline: object
    stage: str
    family: str
    scope: str
    sweep_values: tuple = None

    def swept(self):
        return self.values if self.sweep_values is None else self.sweep_values


AXES = (
    # --- Processing stage -----------------------------------------------------
    # The ONLY axis that changes the reusable processed base dir and forces a
    # full reprocess; analysis axes merely symlink-reuse that base. Currently
    # only "whitestripe" is swept ("ravel" is valid but dormant).
    Axis(
        name="normalization",
        values=("none", "whitestripe", "wmmean", "nyul", "ravel",
                "noneRavel", "wmmeanRavel", "whitestripeNyul", "wmmeanNyul"),
        baseline="whitestripe",
        stage="processing",
        family="processing",
        scope=(
            "T1w/FLAIR intensity source; selects zbrains_SELF (raw, self-normalized) "
            "vs zbrains_WB (WhiteStripe) vs zbrains_RAVEL base dir. 'none' applies no "
            "image-level normalization and defers standardization to analysis-stage "
            "self-normalization."
        ),
        sweep_values=("whitestripe",),
    ),
    # --- Analysis: within-subject spatial rescaling ---------------------------
    Axis(
        name="wscore_preprocessing",
        values=("none", "spatial_zscore", "spatial_robust_z"),
        baseline="none",
        stage="analysis",
        family="spatial_preproc",
        scope="Per-subject spatial standardization of every feature map before the normative fit.",
    ),
    # --- Analysis: normative model (the coupled group) ------------------------
    # wscore_covariate_model drives every cross-axis rule below: choosing a
    # gaussian_process_* kernel constrains distribution / curvature / predictive.
    Axis(
        name="wscore_covariate_model",
        values=(
            "linear",
            "quadratic_age",
            "age_sex_interaction",
            "quadratic_age_sex_interaction",
            "knn",
            "gaussian_process",
            "gaussian_process_rbf_isotropic",
            "gaussian_process_matern32",
            "gaussian_process_matern52",
            "gaussian_process_rational_quadratic",
        ),
        baseline="linear",
        stage="analysis",
        family="normative_model",
        scope="OLS demographic design OR Gaussian-process kernel for the normative fit.",
    ),
    Axis(
        name="wscore_distribution",
        values=("gaussian", "gaussian_mad", "gaussian_winsor", "student_t", "empirical", "shash"),
        baseline="gaussian",
        stage="analysis",
        family="normative_model",
        scope="Distribution fit to normative-model residuals.",
    ),
    Axis(
        name="use_curvature_covariates",
        values=(False, True),
        baseline=False,
        stage="analysis",
        family="normative_model",
        scope="Add vertex-wise curvature as an extra cortical covariate.",
    ),
    Axis(
        name="predictive_wscore",
        values=(False, True),
        baseline=False,
        stage="analysis",
        family="normative_model",
        scope="Add patient prediction uncertainty to the W-score denominator.",
    ),
    # --- Analysis: cortical depth reducers (feature-scoped) -------------------
    # Declared and validated, but NOT swept by default (depth_reducers is not in
    # DEFAULT_ENABLED_FAMILIES), matching the original fixed baseline. They act
    # on disjoint feature families (see the REDUNDANT rule for the one overlap).
    Axis(
        name="blur_depth_model",
        # The 3 literature-grounded metrics the greedy sweeps come first; the
        # legacy reducers remain valid-but-dormant for backward-compat/older runs.
        values=("gray_white_contrast", "boundary_gradient", "juxtacortical_gradient",
                "mean_slope_rms", "mean", "gradient_flattening"),
        baseline="boundary_gradient",
        stage="analysis",
        family="depth_reducers",
        scope=("Reduce the 4-depth *blur profile [mid,white,SWM1,SWM2] to one cortical map; "
               "all 3 swept metrics sample to SWM2 (2mm). gray_white_contrast=pctsurfcon "
               "%contrast GM(mid) vs mean(SWM1,SWM2) (Salat 2009); boundary_gradient=OLS "
               "intensity-profile slope over all 4 depths (Antel 2003/Hong 2014,2017); "
               "juxtacortical_gradient=mean(SWM1,SWM2)-white sub-boundary extension (transmantle/MELD)."),
    ),
    Axis(
        name="intensity_depth_model",
        values=(
            "raw",
            # Self-normalized INTENSITY models (feature meaning preserved):
            "white_swm_zscore",
            "white_swm_mode_zscore",
            "white_swm2_zscore",
            "white_swm_ratio",
            "white_gmwm_contrast",
            "white_percentile_scale",
            "white_surface_robust_z",
            # Depth-contrast models (change feature meaning; not in default sweep):
            "white_swm1_direction_cosine",
            "multisurface_median_abs_dominant",
            # Plain intensity SAMPLERS (surface x self-norm), resolved by the greedy
            # driver from t1w_flair_sample_surface + t1w_flair_self_norm:
            "sample_midthickness", "sample_white", "sample_swm1",
            "sample_midthickness_swm2", "sample_white_swm2", "sample_swm1_swm2",
            "sample_midthickness_owncortex", "sample_white_owncortex", "sample_swm1_owncortex",
        ),
        baseline="multisurface_median_abs_dominant",
        stage="analysis",
        family="depth_reducers",
        scope="Local depth normalization of the fsLR-32k white-surface T1w/FLAIR intensity map.",
    ),
    # --- Analysis: post-scoring surface smoothing -----------------------------
    Axis(
        name="wscore_surface_smoothing_iterations",
        values=(0, 16),
        baseline=0,
        stage="analysis",
        family="postproc",
        scope="Post-score one-ring smoothing steps for fsLR-32k cortex maps.",
    ),
    # --- Control-cohort QC ----------------------------------------------------
    Axis(
        name="control_correlation_filter",
        values=(False, True),
        baseline=False,
        stage="analysis",
        family="control_qc",
        scope="Drop the bottom fraction of controls by mean correlation to other controls (rank-based, per feature).",
    ),
)

AXES_BY_NAME = {axis.name: axis for axis in AXES}
BASELINE = {axis.name: axis.baseline for axis in AXES}
ANALYSIS_AXIS_NAMES = tuple(axis.name for axis in AXES if axis.stage == "analysis")


def _build_families():
    families = {}
    for axis in AXES:
        families.setdefault(axis.family, []).append(axis)
    return families


# family name -> [Axis, ...], in declaration order.
FAMILIES = _build_families()

# Families swept by default. Excludes "depth_reducers" so blur_depth_model and
# intensity_depth_model stay fixed at baseline (their original behaviour).
DEFAULT_ENABLED_FAMILIES = (
    "processing",
    "spatial_preproc",
    "normative_model",
    "postproc",
    "control_qc",
)

# Axes with no cross-axis constraint whatsoever: combine with everything.
INDEPENDENT_AXES = (
    "normalization",
    "wscore_preprocessing",
    "wscore_surface_smoothing_iterations",
    "control_correlation_filter",
)

# Control QC drop-fraction for W/Z scoring when control_correlation_filter is on:
# the RANK-based filter drops this bottom fraction of controls by correlation, PER
# feature (adaptive cut that never collapses a low-correlation feature).
DEFAULT_CONTROL_CORRELATION_QUANTILE = 0.10
# Kept only for compatibility with existing staged-driver call signatures and
# history files.  New correlation-QC arms use the quantile above; the absolute
# threshold is no longer forwarded to ``dataset.analyze``.
DEFAULT_CONTROL_CORRELATION_THRESHOLD = 0.60
# Marker file written once an output directory's analysis has completed.
ANALYSIS_COMPLETION_MARKER = ".analysis_complete_regional_qc_v2"


# ============================================================================
# 2. CROSS-AXIS RULES
#    These are the ONLY couplings between axes. Everything else is independent.
# ============================================================================
def _is_gaussian_process(covariate_model):
    """Gaussian-process kernels are the source of every cross-axis rule."""
    return str(covariate_model).startswith("gaussian_process")


@dataclass(frozen=True)
class ExclusionRule:
    """A cross-axis combination that is forbidden (``analysis.py`` raises)."""

    axes: tuple
    reason: str
    evidence: str
    applies: object  # (config) -> bool ; True means the forbidden pair is present


@dataclass(frozen=True)
class RedundancyRule:
    """A combination that is allowed but silently ignores one setting."""

    axes: tuple
    reason: str
    evidence: str
    applies: object  # (config) -> bool
    normalize_to: tuple = None  # (axis_name, value) forced so labels stay honest


# --- MUTUALLY EXCLUSIVE: rejected up front (would raise inside analysis.py) ---
MUTUALLY_EXCLUSIVE = (
    ExclusionRule(
        axes=("wscore_covariate_model", "wscore_distribution"),
        reason=(
            "Gaussian-process models fit their own predictive Gaussian "
            "distribution, so they cannot be paired with a non-Gaussian "
            "wscore_distribution."
        ),
        evidence="analysis.py:713-717 / 1473-1477 raise ValueError",
        applies=lambda c: _is_gaussian_process(c["wscore_covariate_model"])
        and c["wscore_distribution"] != "gaussian",
    ),
    ExclusionRule(
        axes=("wscore_covariate_model", "use_curvature_covariates"),
        reason=(
            "Gaussian-process models support demographic covariates only; "
            "per-vertex curvature covariates are refused (cortical/blur maps "
            "would be skipped)."
        ),
        evidence="gaussian_process.py:636-640 raises ValueError",
        applies=lambda c: _is_gaussian_process(c["wscore_covariate_model"])
        and bool(c["use_curvature_covariates"]),
    ),
)

# --- REDUNDANT: allowed, but one setting has no effect -----------------------
REDUNDANT = (
    RedundancyRule(
        axes=("wscore_covariate_model", "predictive_wscore"),
        reason=(
            "Gaussian-process models are inherently predictive; "
            "predictive_wscore is never forwarded and has no effect."
        ),
        evidence="gaussian_process.py has no use_prediction_uncertainty param",
        applies=lambda c: _is_gaussian_process(c["wscore_covariate_model"])
        and bool(c["predictive_wscore"]),
        # Force it off so directory names / metadata do not advertise "pred".
        normalize_to=("predictive_wscore", False),
    ),
    RedundancyRule(
        axes=("intensity_depth_model", "blur_depth_model"),
        reason=(
            "With intensity_depth_model='multisurface_median_abs_dominant' the "
            "white-surface T1w/FLAIR intensity map ignores blur_depth_model "
            "(intensity reduction takes precedence). The two otherwise act on "
            "disjoint feature families, so this is a map-scoped no-op only."
        ),
        evidence="analysis.py:768 precedence over blur branches at 788/791/817",
        applies=lambda c: c["intensity_depth_model"] == "multisurface_median_abs_dominant"
        and c["blur_depth_model"] != BASELINE["blur_depth_model"],
        normalize_to=None,  # narrow / feature-scoped -> warn only, do not force
    ),
)


def validate_configuration(config, *, strict=True):
    """Validate one configuration against the axis domains and cross-axis rules.

    Parameters
    ----------
    config : dict
        A full configuration (all axis names -> chosen value).
    strict : bool, default True
        If True, raise on any mutually-exclusive combination. If False, such
        combinations are reported in the returned warnings instead of raising.

    Returns
    -------
    (normalized_config, warnings) : (dict, list[str])
        ``normalized_config`` has any redundant settings coerced to their
        no-effect baseline; ``warnings`` describes redundancies that were found.

    Raises
    ------
    ValueError
        On an unknown value for an axis, or (when ``strict``) on a
        mutually-exclusive combination.
    """
    normalized = dict(config)
    warnings = []

    # Within-axis domain: each value must be one of the axis's allowed choices.
    for name, value in config.items():
        axis = AXES_BY_NAME.get(name)
        if axis is not None and value not in axis.values:
            raise ValueError(
                f"{name}={value!r} is not a valid choice; pick one of {list(axis.values)}"
            )

    # Cross-axis mutual exclusions.
    for rule in MUTUALLY_EXCLUSIVE:
        if rule.applies(config):
            message = (
                f"Mutually exclusive: {' + '.join(rule.axes)} -- {rule.reason} "
                f"[{rule.evidence}]"
            )
            if strict:
                raise ValueError(message)
            warnings.append(message)

    # Cross-axis redundancies.
    for rule in REDUNDANT:
        if rule.applies(config):
            warnings.append(f"Redundant: {' + '.join(rule.axes)} -- {rule.reason}")
            if rule.normalize_to is not None:
                axis_name, value = rule.normalize_to
                normalized[axis_name] = value

    return normalized, warnings


# ============================================================================
# 3. ONE-FACTOR-AT-A-TIME CONFIGURATION GENERATOR
# ============================================================================
def processing_configurations(enabled_families=DEFAULT_ENABLED_FAMILIES):
    """Yield the baseline plus one-factor-at-a-time variants.

    Each variant differs from :data:`BASELINE` in exactly one axis. Variants
    are validated (and redundancy-normalized) here; any that violate a
    mutually-exclusive rule are skipped with a printed note.

    Set the ``ZBRAINS_BENCHMARK_OPTIONS`` environment variable to a comma-
    separated list of variant labels (e.g. ``"baseline,wscore_distribution=shash"``)
    to run only those.

    Returns
    -------
    list[(family, label, config, warnings)]
    """
    unknown = set(enabled_families) - set(FAMILIES)
    if unknown:
        raise ValueError(
            f"Unknown toggle families: {sorted(unknown)} (known: {sorted(FAMILIES)})"
        )

    raw = [("baseline", "baseline", dict(BASELINE))]
    for family in enabled_families:
        for axis in FAMILIES[family]:
            for value in axis.swept():
                if value == axis.baseline:
                    continue
                config = dict(BASELINE)
                config[axis.name] = value
                raw.append((family, f"{axis.name}={value}", config))

    configurations = []
    for family, label, config in raw:
        try:
            normalized, warnings = validate_configuration(config)
        except ValueError as exc:
            print(f"[benchmark] skipping invalid variant {label!r}: {exc}")
            continue
        configurations.append((family, label, normalized, warnings))

    requested = os.environ.get("ZBRAINS_BENCHMARK_OPTIONS", "").strip()
    if requested:
        allowed = {value.strip() for value in requested.split(",") if value.strip()}
        configurations = [item for item in configurations if item[1] in allowed]

    return configurations


# ============================================================================
# 4. OUTPUT-DIRECTORY NAMING
# ============================================================================
def quantile_label(quantile):
    """Filesystem-safe directory label for the RANK-based (drop-bottom-fraction)
    control correlation filter. Scalar fraction or per-feature mapping; distinct
    quantile arms get distinct, deterministic labels (unwieldy per-feature maps
    collapse to a stable short hash).
    """
    if isinstance(quantile, Mapping):
        items = sorted(
            (re.sub(r"[^0-9A-Za-z]", "", str(name)),
             "off" if value is None else str(value).replace(".", "p"))
            for name, value in quantile.items()
        )
        label = "corrqpf_" + "_".join(f"{name}{value}" for name, value in items)
        if len(label) > 100:
            digest = hashlib.md5(label.encode("utf-8")).hexdigest()[:12]
            label = f"corrqpf_{digest}"
        return label
    return f"corrq{str(quantile).replace('.', 'p')}"


# Legacy single-mode desc -> processed-base dir token. Composite labels
# (noneRavel, wmmeanRavel, whitestripeNyul, wmmeanNyul) get a unique uppercase
# alphanumeric token so distinct composed arms never share a base directory.
_NORM_DIR_LABEL = {"ravel": "RAVEL", "none": "SELF", "wmmean": "WMMEAN",
                   "nyul": "NYUL", "whitestripe": "WB"}


def _normalization_dir_label(normalization):
    if normalization in _NORM_DIR_LABEL:
        return _NORM_DIR_LABEL[normalization]
    return re.sub(r"[^0-9A-Za-z]+", "", str(normalization)).upper()


def output_directory_for(config, output_dir_prefix,
                         control_correlation_quantile=DEFAULT_CONTROL_CORRELATION_QUANTILE,
                         exclusion_signature=""):
    """Deterministic analysis output directory for a configuration.

    ``exclusion_signature`` (empty when no pre-dataset-norm control exclusion is
    active) disambiguates analysis dirs of distinct exclusion arms.
    ``control_correlation_quantile`` is the RANK-based control-filter drop-fraction
    (scalar or per-feature mapping), encoded when the filter is on so distinct
    quantile arms get distinct dirs. It is passed explicitly (not read from
    ``config``) because the caller strips driver-managed keys before building the
    axis config.
    """
    norm_label = _normalization_dir_label(config["normalization"])
    parts = [f"zbrains_{norm_label}"]
    if config["use_curvature_covariates"]:
        parts.append("curv")
    if config["predictive_wscore"]:
        parts.append("pred")
    if config["wscore_distribution"] != "gaussian":
        parts.append(config["wscore_distribution"])
    if config["wscore_preprocessing"] != "none":
        parts.append(config["wscore_preprocessing"])
    if config["wscore_covariate_model"] != "linear":
        parts.append(config["wscore_covariate_model"])
    if config["wscore_surface_smoothing_iterations"]:
        parts.append(f"wsmooth{config['wscore_surface_smoothing_iterations']}")
    blur_labels = {
        "gray_white_contrast": "blurgwc",
        "boundary_gradient": "blurboundarygrad",
        "juxtacortical_gradient": "blurjuxtagrad",
        "mean_slope_rms": "blurrms",
        "mean": "blurmean",
        "gradient_flattening": "blurgradientflatteningmedabs",
    }
    parts.append(blur_labels[config["blur_depth_model"]])
    # Every non-"raw" intensity model must get a UNIQUE directory suffix, else
    # distinct sampling/self-norm arms collide on one path and reuse stale maps.
    idm = config["intensity_depth_model"]
    if idm != "raw":
        _idm_short = {"white_swm1_direction_cosine": "depthcosine",
                      "multisurface_median_abs_dominant": "multisurfaceintensity"}
        parts.append(_idm_short.get(idm, "idm" + re.sub(r"[^0-9a-z]+", "", idm.lower())))
    if config["control_correlation_filter"]:
        parts.append(quantile_label(control_correlation_quantile))
    else:
        parts.append("allcontrols")
    if exclusion_signature:
        parts.append(exclusion_signature)
    return f"{output_dir_prefix}/{'_'.join(parts)}/"


def processed_base_directory_for(normalization, output_dir_prefix, exclusion_signature=""):
    """Directory of reusable processed maps, independent of analysis options.

    ``exclusion_signature`` (empty for no pre-dataset-norm control exclusion)
    keys the base so distinct exclusion arms never reuse a stale RAVEL/Nyul fit or
    normative model. Empty reproduces the exact legacy path (backward compatible).
    """
    norm_label = _normalization_dir_label(normalization)
    suffix = f"_{exclusion_signature}" if exclusion_signature else ""
    return f"{output_dir_prefix}/zbrains_{norm_label}{suffix}/"


def hippunfold_signature_for_datasets(*datasets):
    """Return one shared HippUnfold source tag, rejecting mixed inputs."""
    signatures = set()
    for dataset in datasets:
        directory = getattr(dataset, "hippunfold_directory", None)
        if (
            not isinstance(directory, (str, os.PathLike))
            or not directory
            or not getattr(dataset, "hippocampus", False)
        ):
            continue
        version = getattr(dataset, "requested_hippunfold_version", None)
        if version is None and getattr(dataset, "_hippunfold_version_detected", False):
            version = getattr(dataset, "hippunfold_version", None)
        signatures.add(hippunfold_cache_tag(directory, version=version))
    if len(signatures) > 1:
        raise ValueError(
            "Control and patient datasets use different HippUnfold sources or versions."
        )
    return next(iter(signatures), "")


# ============================================================================
# 5. PROCESSED-OUTPUT COMPLETENESS CHECKS AND SYMLINK REUSE
# ============================================================================
def subject_output_directories(dataset):
    if not hasattr(dataset, "valid_subjects"):
        dataset.check_directories()
    valid_subjects = dataset.valid_subjects.get("base", [])
    if valid_subjects:
        return [Path(participant) / session for participant, session in valid_subjects]
    subject_dirs = []
    for _, row in dataset.demographics.data.iterrows():
        participant_id = row["participant_id"]
        session_id = row.get("session_id", "ses-01")
        subject_dirs.append(Path(participant_id) / session_id)
    return subject_dirs


def processed_maps_complete(output_directory, datasets):
    output_directory = Path(output_directory)
    expected = 0
    found = 0
    for dataset in datasets:
        for relative_dir in subject_output_directories(dataset):
            expected += 1
            subject_dir = output_directory / relative_dir
            if (subject_dir / "maps").exists() and (subject_dir / "structural").exists():
                found += 1
    return expected > 0 and found == expected


def subjects_missing_non_hippocampal_processed_outputs(
    output_directory,
    dataset,
):
    """Subjects that cannot safely enter hippocampus-only processing.

    A subject is ready when its enabled cortex/subcortex map trees and reusable
    non-HippUnfold structural geometry are present. Hippocampal maps/surfaces are
    deliberately ignored because those are exactly what the V2 migration rebuilds.
    """
    output_directory = Path(output_directory)
    missing = []
    for relative_dir in subject_output_directories(dataset):
        session = output_directory / relative_dir
        maps = session / "maps"
        structural = session / "structural"
        ready = True
        for structure in ("cortex", "subcortical"):
            if not getattr(dataset, structure, False):
                continue
            structure_dir = maps / structure
            if not (
                structure_dir.is_dir()
                and any(path.is_file() for path in structure_dir.rglob("*"))
            ):
                ready = False
        if not (
            structural.is_dir()
            and any(
                path.is_file() and "_label-hipp_" not in path.name
                for path in structural.iterdir()
            )
        ):
            ready = False
        if not ready:
            missing.append(tuple(relative_dir.parts[-2:]))
    return missing


# Per-base completion sentinel: written ATOMICALLY only after a base has finished
# processing BOTH its controls and patients (see zbrains_staged._ensure_base_processed).
# Two machines sharing storage (e.g. forward/reverse objective order) can then trust a
# base as complete in O(1) WITHOUT re-scanning thousands of map files, and -- crucially
# -- never treat a base another machine is still writing as complete (the sentinel is
# the LAST thing written, so its presence means "fully done"). A base that is complete
# on disk but predates sentinels is adopted (stamped) with no reprocessing.
_BASE_COMPLETE_SENTINEL = ".zbrains_maps_complete"


def _sentinel_pairs(pairs):
    return sorted({(str(p[0]), str(p[1])) for p in (pairs or ())})


def _read_complete_subjects(base_dir):
    """The (pid, ses) set the sentinel certifies, or None if there is no sentinel."""
    path = os.path.join(str(base_dir), _BASE_COMPLETE_SENTINEL)
    try:
        with open(path, encoding="utf-8") as fh:
            rows = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    except OSError:
        return None
    return {tuple(r.split("\t", 1)) for r in rows if "\t" in r}


def base_is_marked_complete(base_dir, required_pairs=None):
    """True iff the base sentinel exists AND (when ``required_pairs`` is given) it
    certifies EVERY required (pid, ses). Recording the subject set means a base is
    NOT wrongly trusted after the participant CSV gains subjects -- the new subjects
    aren't in the sentinel, so it re-processes them (and re-stamps)."""
    recorded = _read_complete_subjects(base_dir)
    if recorded is None:
        return False
    if required_pairs is None:
        return True
    return set(_sentinel_pairs(required_pairs)).issubset(recorded)


def mark_base_complete(base_dir, subject_pairs=()):
    """Atomically stamp ``base_dir`` complete, recording the (pid, ses) it covers
    (tmp file + os.replace so a reader never sees a half-written sentinel)."""
    base_dir = str(base_dir)
    if not os.path.isdir(base_dir):
        return
    path = os.path.join(base_dir, _BASE_COMPLETE_SENTINEL)
    tmp = os.path.join(base_dir, f".{_BASE_COMPLETE_SENTINEL}.tmp.{os.getpid()}")
    lines = ["# subjects (pid<TAB>ses) certified fully processed in this base"]
    lines += [f"{pid}\t{ses}" for pid, ses in _sentinel_pairs(subject_pairs)]
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        os.replace(tmp, path)
    except OSError:
        try:
            if os.path.exists(tmp):
                os.remove(tmp)
        except OSError:
            pass


def unmark_base_complete(base_dir):
    """Remove the completion sentinel (call BEFORE (re)writing a base so a crash
    mid-rebuild can't leave a stale 'complete' mark on a now-partial base)."""
    path = os.path.join(str(base_dir), _BASE_COMPLETE_SENTINEL)
    try:
        os.remove(path)
    except OSError:
        pass


def score_maps_complete(output_directory, dataset, method="wscore"):
    output_directory = Path(output_directory)
    expected = 0
    found = 0
    for relative_dir in subject_output_directories(dataset):
        expected += 1
        score_dir = output_directory / relative_dir / f"{method}_maps"
        if score_dir.exists() and any(path.is_file() for path in score_dir.rglob("*")):
            found += 1
    return expected > 0 and found == expected


def analysis_structure_maps_complete(
    output_directory,
    datasets,
    method="wscore",
    structures=("cortex", "hippocampus", "subcortical"),
):
    """Check every enabled requested score-map structure for every subject.

    ``score_maps_complete`` intentionally implements a lightweight historical
    check (any score file per subject). Selective HippUnfold reuse needs a stricter
    check so a partial score tree can never be mistaken for a safely reusable one.
    """
    output_directory = Path(output_directory)
    expected = 0
    found = 0
    for dataset in datasets:
        enabled = [
            structure
            for structure in structures
            if getattr(dataset, structure, False)
        ]
        for relative_dir in subject_output_directories(dataset):
            for structure in enabled:
                expected += 1
                structure_dir = (
                    output_directory
                    / relative_dir
                    / f"{method}_maps"
                    / structure
                )
                if (
                    structure_dir.is_dir()
                    and any(path.is_file() for path in structure_dir.rglob("*"))
                ):
                    found += 1
    return expected > 0 and found == expected


def subjects_missing_analysis_structure(
    output_directory,
    dataset,
    structure,
    method="wscore",
):
    """Return subject/session pairs missing one enabled analysis structure."""
    if not getattr(dataset, structure, False):
        return []
    output_directory = Path(output_directory)
    missing = []
    for relative_dir in subject_output_directories(dataset):
        structure_dir = (
            output_directory
            / relative_dir
            / f"{method}_maps"
            / structure
        )
        if not (
            structure_dir.is_dir()
            and any(path.is_file() for path in structure_dir.rglob("*"))
        ):
            missing.append(tuple(relative_dir.parts[-2:]))
    return missing


def non_hippocampal_analysis_outputs_complete(
    output_directory,
    datasets,
    method="wscore",
):
    """Whether all enabled cortex/subcortex score trees can be safely reused."""
    return analysis_structure_maps_complete(
        output_directory,
        datasets,
        method=method,
        structures=("cortex", "subcortical"),
    )


def output_fully_processed(output_directory, patient_dataset, marker=ANALYSIS_COMPLETION_MARKER):
    output_directory = Path(output_directory)
    marker_path = output_directory / marker
    if not (
        processed_maps_complete(output_directory, [patient_dataset])
        and score_maps_complete(output_directory, patient_dataset)
        and marker_path.exists()
    ):
        return False
    marker_mtime = marker_path.stat().st_mtime
    for relative_dir in subject_output_directories(patient_dataset):
        maps_directory = output_directory / relative_dir / "maps"
        if any(
            path.is_file() and path.stat().st_mtime > marker_mtime
            for path in maps_directory.rglob("*")
        ):
            print(f"Processed maps are newer than {marker_path}; W-scores will be regenerated.")
            return False
    return True


def symlink_path(source, target):
    source = Path(source)
    target = Path(target)
    if not source.exists():
        return False
    if target.is_symlink():
        if target.resolve() == source.resolve():
            return True
        # Stale symlink: a prior run pointed this at a DIFFERENT base (e.g. before the
        # base-signature change repointed exclusion arms at the shared base, or any
        # base-keying evolution). It is a pointer, not data, and this function is
        # always called with the correct current source -- so REPOINT it rather than
        # aborting the run. Do it ATOMICALLY (temp symlink + os.replace) so a peer
        # machine sharing storage never sees the target briefly missing. (A REAL
        # directory at the target is left untouched below.)
        tmp = target.with_name(f".{target.name}.repoint.{os.getpid()}")
        if tmp.is_symlink() or tmp.exists():
            tmp.unlink()
        tmp.symlink_to(source, target_is_directory=source.is_dir())
        os.replace(tmp, target)
        return True
    if target.exists():
        return True
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        target.symlink_to(source, target_is_directory=source.is_dir())
    except FileExistsError:
        # TOCTOU race: a peer process (two-machine shared storage building the SAME
        # config's analysis dir) created the target between the checks above and
        # here. Idempotent -- accept whatever is now there if it already resolves to
        # our source (or a peer materialized a real dir); otherwise repoint it
        # atomically so we never abort the whole run on a benign collision.
        if target.is_symlink() and target.resolve() == source.resolve():
            return True
        if target.exists() and not target.is_symlink():
            return True
        tmp = target.with_name(f".{target.name}.repoint.{os.getpid()}")
        if tmp.is_symlink() or tmp.exists():
            tmp.unlink()
        tmp.symlink_to(source, target_is_directory=source.is_dir())
        os.replace(tmp, target)
    return True


def symlink_processed_outputs(source_directory, target_directory, datasets):
    source_directory = Path(source_directory)
    target_directory = Path(target_directory)
    target_directory.mkdir(parents=True, exist_ok=True)

    if not source_directory.exists():
        raise FileNotFoundError(
            f"Cannot symlink processed outputs because source does not exist: {source_directory}"
        )

    linked_maps = 0
    linked_structural = 0
    for dataset in datasets:
        for relative_dir in subject_output_directories(dataset):
            source_maps = source_directory / relative_dir / "maps"
            target_maps = target_directory / relative_dir / "maps"
            if symlink_path(source_maps, target_maps):
                linked_maps += 1

            source_structural = source_directory / relative_dir / "structural"
            target_structural = target_directory / relative_dir / "structural"
            if symlink_path(source_structural, target_structural):
                linked_structural += 1

    symlink_path(source_directory / "ravel", target_directory / "ravel")
    print(
        f"Reused {linked_maps} map and {linked_structural} structural directories: "
        f"{source_directory} -> {target_directory}"
    )


def seed_non_hippocampal_processed_outputs(
    source_directory,
    target_directory,
    datasets,
):
    """Link legacy non-hippocampal products into a versioned HippUnfold base.

    ``maps/hippocampus`` and every ``label-hipp`` structural surface are excluded.
    The target therefore reuses cortical, subcortical, normalized-volume, SWM, and
    dataset-normalization products without exposing any legacy hippocampal data.
    Returns ``False`` when the legacy source is absent or lacks a required enabled
    non-hippocampal component; callers can then rebuild only the missing subjects.
    """
    from zbrains.dataset import _link_structural_to_cache

    source_directory = Path(source_directory)
    target_directory = Path(target_directory)
    if not source_directory.is_dir():
        return False
    target_directory.mkdir(parents=True, exist_ok=True)

    complete = True
    linked = 0
    for dataset in datasets:
        for relative_dir in subject_output_directories(dataset):
            source_session = source_directory / relative_dir
            source_maps = source_session / "maps"
            target_maps = target_directory / relative_dir / "maps"
            if not source_maps.is_dir():
                complete = False
                continue
            target_maps.mkdir(parents=True, exist_ok=True)
            for source_entry in source_maps.iterdir():
                if (
                    source_entry.name == "hippocampus"
                    or "_label-hipp_" in source_entry.name
                ):
                    continue
                if symlink_path(source_entry, target_maps / source_entry.name):
                    linked += 1

            if (
                getattr(dataset, "cortex", False)
                and not (
                    (target_maps / "cortex").is_dir()
                    and any(
                        path.is_file()
                        for path in (target_maps / "cortex").rglob("*")
                    )
                )
            ):
                complete = False
            if (getattr(dataset, "subcortical", False)
                    and not (
                        (target_maps / "subcortical").is_dir()
                        and any(
                            path.is_file()
                            for path in (target_maps / "subcortical").rglob("*")
                        )
                    )):
                complete = False

            pid, sid = relative_dir.parts[-2:]
            target_structural = _link_structural_to_cache(
                target_directory,
                pid,
                sid,
                hippunfold_directory=getattr(dataset, "hippunfold_directory", None),
                hippunfold_version=(
                    getattr(dataset, "requested_hippunfold_version", None)
                    or getattr(dataset, "hippunfold_version", None)
                ),
            )
            source_structural = source_session / "structural"
            if not source_structural.is_dir() or target_structural is None:
                complete = False
                continue
            target_structural = Path(target_structural)
            for source_entry in source_structural.iterdir():
                if "_label-hipp_" in source_entry.name:
                    continue
                if symlink_path(source_entry, target_structural / source_entry.name):
                    linked += 1
            if not any(target_structural.iterdir()):
                complete = False

    for model_dir in ("ravel", "nyul"):
        source_model = source_directory / model_dir
        if source_model.exists() and symlink_path(
            source_model,
            target_directory / model_dir,
        ):
            linked += 1

    print(
        f"Seeded {linked} non-hippocampal entries: "
        f"{source_directory} -> {target_directory}"
    )
    return complete


def hippocampal_v2_outputs_complete(output_directory, datasets):
    """Whether every subject has all expected V2 maps and bilateral V2 surfaces."""
    output_directory = Path(output_directory)
    output_tokens = {
        "thickness": "thickness",
        "sa": "SA",
        "t1map": "qT1",
        "qt1": "qT1",
        "adc": "ADC",
        "fa": "FA",
        "flair": "FLAIR",
        "t1w": "T1w",
    }
    expected = 0
    found = 0
    for dataset in datasets:
        features = [
            feature
            for feature in (getattr(dataset, "features", None) or ["thickness"])
            if not str(feature).endswith("*blur")
        ]
        source_valid = getattr(dataset, "source_valid_subjects", None) or {}
        for relative_dir in subject_output_directories(dataset):
            expected += 1
            pid, sid = relative_dir.parts[-2:]
            subject = (pid, sid)
            session = output_directory / relative_dir
            maps = session / "maps" / "hippocampus"
            structural = session / "structural"
            expected_features = []
            for feature in features:
                validity = source_valid.get(feature, {})
                valid_hippocampus = validity.get("structures", {}).get(
                    "hippocampus")
                if valid_hippocampus is not None and subject not in valid_hippocampus:
                    continue
                expected_features.append(
                    output_tokens.get(str(feature).lower(), str(feature).lower())
                )
            maps_ok = bool(expected_features) and all(
                any(maps.glob(
                    f"{pid}_{sid}_hemi-{hemi}_den-8k_label-hipp_midthickness_"
                    f"feature-{feature}_smooth-*mm.func.gii"
                ))
                for feature in expected_features
                for hemi in ("L", "R")
            )
            surfaces_ok = all(
                (
                    structural
                    / f"{pid}_{sid}_hemi-{hemi}_space-{space}_den-8k_"
                    "label-hipp_midthickness.surf.gii"
                ).is_file()
                for hemi in ("L", "R")
                for space in ("T1w", "unfold")
            )
            if maps_ok and surfaces_ok:
                found += 1
    return expected > 0 and found == expected


def seed_non_hippocampal_analysis_outputs(
    source_directory,
    target_directory,
    datasets,
    method,
):
    """Link cortex/subcortical score-map trees while excluding hippocampus."""
    source_directory = Path(source_directory)
    target_directory = Path(target_directory)
    if not source_directory.is_dir():
        return False

    complete = True
    linked = 0
    score_name = f"{method}_maps"
    for dataset in datasets:
        for relative_dir in subject_output_directories(dataset):
            source_scores = source_directory / relative_dir / score_name
            target_scores = target_directory / relative_dir / score_name
            if not source_scores.is_dir():
                complete = False
                continue
            target_scores.mkdir(parents=True, exist_ok=True)
            for structure in ("cortex", "subcortical"):
                if structure == "cortex" and not getattr(dataset, "cortex", False):
                    continue
                if (structure == "subcortical"
                        and not getattr(dataset, "subcortical", False)):
                    continue
                source_structure = source_scores / structure
                if (
                    not source_structure.is_dir()
                    or not any(path.is_file() for path in source_structure.rglob("*"))
                ):
                    complete = False
                    continue
                if symlink_path(source_structure, target_scores / structure):
                    linked += 1
    print(
        f"Seeded {linked} non-hippocampal {score_name} directories: "
        f"{source_directory} -> {target_directory}"
    )
    return complete


# ============================================================================
# 6. SUMMARY / DESCRIPTION HELPERS
# ============================================================================
_CONFIG_PRINT_ORDER = (
    "normalization",
    "use_curvature_covariates",
    "predictive_wscore",
    "wscore_distribution",
    "wscore_preprocessing",
    "wscore_covariate_model",
    "wscore_surface_smoothing_iterations",
    "blur_depth_model",
    "intensity_depth_model",
    "control_correlation_filter",
)


def _describe_config(config, quantile):
    parts = [f"{key}={config[key]}" for key in _CONFIG_PRINT_ORDER]
    parts.append(f"control_corr_quantile={quantile}")
    return ", ".join(parts)


def _summary_row(category, option, config, output_directory, quantile, skipped):
    return {
        "category": category,
        "option": option,
        "config": dict(config),
        "output_directory": output_directory,
        "quantile": quantile,
        "skipped": skipped,
    }


def _format_summary(summary):
    config = summary["config"]
    text = (
        f"category={summary['category']}, "
        f"processing_option={summary['option']}, "
        f"normalization={config['normalization']}, "
        f"curvature_covariates={config['use_curvature_covariates']}, "
        f"predictive_wscore={config['predictive_wscore']}, "
        f"wscore_distribution={config['wscore_distribution']} -> "
        f"wscore_preprocessing={config['wscore_preprocessing']} -> "
        f"wscore_covariate_model={config['wscore_covariate_model']} -> "
        f"wscore_surface_smoothing_iterations={config['wscore_surface_smoothing_iterations']} -> "
        f"blur_depth_model={config['blur_depth_model']} -> "
        f"intensity_depth_model={config['intensity_depth_model']} -> "
        f"control_corr_filter={config['control_correlation_filter']} "
        f"(drop-bottom {summary['quantile']}) -> "
        f"{summary['output_directory']}"
    )
    if summary["skipped"]:
        text += " [skipped existing]"
    return text


def format_decision_space():
    """Return a human-readable map of what combines and what is exclusive."""
    lines = [
        "z-brains benchmark decision space",
        "=" * 72,
        "",
        "Each AXIS is one engineering decision. Values WITHIN an axis are",
        "mutually exclusive (a run picks exactly one). Axes are independent",
        "(freely combinable) unless a cross-axis rule below says otherwise.",
        "(*) marks the baseline value. The sweep varies one axis at a time.",
        "",
    ]
    for family, axes in FAMILIES.items():
        swept = "swept" if family in DEFAULT_ENABLED_FAMILIES else "NOT swept by default"
        lines.append(f"[{family}]  ({axes[0].stage} stage, {swept})")
        for axis in axes:
            rendered = ", ".join(
                f"{value}*" if value == axis.baseline else f"{value}" for value in axis.values
            )
            lines.append(f"  - {axis.name}: {rendered}")
            lines.append(f"      {axis.scope}")
        lines.append("")

    lines.append("Fully independent axes (combine with everything):")
    for name in INDEPENDENT_AXES:
        lines.append(f"  - {name}")
    lines.append("")

    lines.append("Mutually exclusive combinations (rejected up front):")
    for rule in MUTUALLY_EXCLUSIVE:
        lines.append(f"  x {' + '.join(rule.axes)}")
        lines.append(f"      {rule.reason}")
        lines.append(f"      enforced by {rule.evidence}")
    lines.append("")

    lines.append("Redundant combinations (allowed; one setting is ignored):")
    for rule in REDUNDANT:
        lines.append(f"  ~ {' + '.join(rule.axes)}")
        lines.append(f"      {rule.reason}")
    return "\n".join(lines)


# ============================================================================
# 7. DRIVER
# ============================================================================
def run_benchmark(
    *,
    control_dataset,
    patient_dataset,
    env,
    output_dir_prefix,
    base_pipeline_settings,
    enabled_families=DEFAULT_ENABLED_FAMILIES,
    control_correlation_quantile=DEFAULT_CONTROL_CORRELATION_QUANTILE,
    reprocess=False,
    reprocess_controls=False,
    dataset_label="",
):
    """Run the one-factor-at-a-time benchmark sweep.

    Parameters
    ----------
    control_dataset, patient_dataset : zbdataset
        Reference (control) and patient datasets.
    env : zbenv
        Processing environment (Workbench settings, thread counts).
    output_dir_prefix : str
        Root directory prefix for derivative outputs.
    base_pipeline_settings : dict
        Shared process/validate settings (features, smoothing) WITHOUT
        ``normalization`` -- it is added per configuration.
    enabled_families : sequence[str]
        Which OFAT families to sweep (see :data:`FAMILIES`).
    control_correlation_quantile : float
        RANK-based drop-fraction used when ``control_correlation_filter`` is on
        (drop the bottom fraction of controls by correlation, per feature).
    reprocess : bool
        Force regeneration of patient derivatives in the base directory.
    reprocess_controls : bool
        Rebuild the control cohort in every normalization-specific base
        directory used by this run before scoring.
    dataset_label : str
        Optional label (e.g. "NOEL") woven into log lines.
    """
    label = f"{dataset_label} " if dataset_label else ""
    run_summaries = []
    reprocessed_control_base_directories = set()
    hippunfold_signature = hippunfold_signature_for_datasets(
        control_dataset,
        patient_dataset,
    )

    for category, option, config, warnings in processing_configurations(enabled_families):
        normalization = config["normalization"]
        output_directory = output_directory_for(
            config,
            output_dir_prefix,
            control_correlation_quantile,
            exclusion_signature=hippunfold_signature,
        )
        pipeline_settings = dict(base_pipeline_settings, normalization=normalization)
        base_output_directory = processed_base_directory_for(
            normalization,
            output_dir_prefix,
            hippunfold_signature,
        )
        is_base_output = (
            Path(output_directory).resolve() == Path(base_output_directory).resolve()
        )

        print("\n" + "=" * 80)
        print(
            f"Running z-brains {label}combination: "
            f"category={category}, processing_option={option}, "
            f"{_describe_config(config, control_correlation_quantile)}"
        )
        print(f"Output directory: {output_directory}")
        for warning in warnings:
            print(f"  NOTE: {warning}")
        print("=" * 80)

        resolved_base_output = str(Path(base_output_directory).resolve())
        if (
            reprocess_controls
            and resolved_base_output not in reprocessed_control_base_directories
        ):
            print(
                f"Reprocessing the original {label}control cohort in the reusable "
                f"base directory: {base_output_directory}"
            )
            control_dataset.process(
                output_directory=base_output_directory,
                env=env,
                verbose=True,
                **pipeline_settings,
            )
            reprocessed_control_base_directories.add(resolved_base_output)

        if not reprocess_controls and output_fully_processed(output_directory, patient_dataset):
            print(
                "Output directory already has complete processed maps and patient "
                "W-score maps; skipping."
            )
            run_summaries.append(
                _summary_row(
                    category, option, config, output_directory,
                    control_correlation_quantile, skipped=True,
                )
            )
            continue

        if is_base_output:
            if reprocess or not processed_maps_complete(output_directory, [patient_dataset]):
                if not reprocess:
                    print(
                        "Patient processed outputs not found; processing patients "
                        "in the base directory."
                    )
                patient_dataset.process(
                    output_directory=output_directory,
                    env=env,
                    verbose=True,
                    **pipeline_settings,
                )
            else:
                print("Base processed outputs found; skipping processing.")
        else:
            if reprocess:
                print(
                    "Derived analysis directory requested with reprocess=True; "
                    "reusing base processed outputs via symlinks."
                )
            if not processed_maps_complete(base_output_directory, [patient_dataset]):
                print(
                    "Patient outputs are missing from the base directory; "
                    f"processing patients first: {base_output_directory}"
                )
                patient_dataset.process(
                    output_directory=base_output_directory,
                    env=env,
                    verbose=True,
                    **pipeline_settings,
                )
            symlink_processed_outputs(
                base_output_directory,
                output_directory,
                datasets=[control_dataset, patient_dataset],
            )

        control_dataset.validate(
            output_directory=output_directory, verbose=False, **pipeline_settings
        )
        patient_dataset.validate(
            output_directory=output_directory, verbose=False, **pipeline_settings
        )

        analysis_kwargs = {name: config[name] for name in ANALYSIS_AXIS_NAMES}
        patient_dataset.analyze(
            output_directory=output_directory,
            reference=control_dataset,
            method="wscore",
            control_correlation_quantile=control_correlation_quantile,
            **analysis_kwargs,
        )
        (Path(output_directory) / ANALYSIS_COMPLETION_MARKER).touch()

        run_summaries.append(
            _summary_row(
                category, option, config, output_directory,
                control_correlation_quantile, skipped=False,
            )
        )

    print(f"\nCompleted z-brains {label}sweep:")
    for summary in run_summaries:
        print("  " + _format_summary(summary))
    print("done")
    return run_summaries


if __name__ == "__main__":
    print(format_decision_space())
