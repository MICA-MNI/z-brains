#!/usr/bin/env python3
"""Per-subject macro AUROC scoring for the greedy staged optimizer.

The greedy driver (``zbrains_benchmark.run_staged_optimization``) needs a single
scalar per candidate configuration: the per-subject macro AUROC for detecting
lesion vertices from cortical W-/Z-score maps, pooled across the MICs and NOEL
cohorts. This module provides exactly that, reusing the low-level primitives of
``results/make_lesion_wscore_overlap.py`` (lesion discovery, GIFTI reading,
tie-corrected Mann-Whitney AUC) rather than its pooled-vertex ROC machinery
(``compute_roc_curves`` pools every subject's vertices together, which is a
different metric).

Directory layout scored (per ``analysis.py``'s ``f"{method}_maps"``):

    <score_root>/<subject>/<session>/<method>_maps/cortex/
        sub-*_ses-*_hemi-{L,R}_surf-fsLR-32k_label-{white,midthickness}_
        feature-*_smooth-*mm_analysis-{regional,asymmetry}.func.gii

``method`` is ``"wscore"`` or ``"zscore"`` -- this module is method-aware so the
simplest-baseline z-score stage and the w-score/GP stages all score identically.

Usage (import):

    from greedy_auroc import pooled_macro_auroc, DEFAULT_ROOTS_FN
    score = pooled_macro_auroc({"MICs": mics_root, "NOEL": noel_root})

Usage (CLI, for verification / standalone scoring of one directory pair):

    python results/greedy_auroc.py \
        --mics-root /path/zbrains_WB_.../ --noel-root /path/zbrains_WB_.../ \
        --method wscore
"""

from __future__ import annotations

import argparse
import importlib.util
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
HEMI_VERTICES = 32492


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# The harness filename is not importable as a package module; load it by path.
_LA = _load_module("lesion_analysis_for_greedy", REPO_ROOT / "results/make_lesion_wscore_overlap.py")
_DS = _load_module("benchmark_datasets_for_greedy", REPO_ROOT / "results/benchmark_datasets.py")

DATASETS_BY_NAME = {d["name"]: d for d in _DS.DATASETS}


def dataset_lesion_dir(dataset_name: str) -> Path:
    try:
        return Path(DATASETS_BY_NAME[dataset_name]["lesion_dir"])
    except KeyError as exc:
        raise KeyError(
            f"Unknown dataset {dataset_name!r}; known: {sorted(DATASETS_BY_NAME)}"
        ) from exc


def dataset_hipp_density(dataset_name: str):
    """Configured hippocampal density for a cohort's feature and cavity maps."""
    try:
        return DATASETS_BY_NAME[dataset_name].get("hipp_density")
    except KeyError as exc:
        raise KeyError(
            f"Unknown dataset {dataset_name!r}; known: {sorted(DATASETS_BY_NAME)}"
        ) from exc


def dataset_default_root(dataset_name: str) -> Path:
    """The cohort's canonical zbrains_WB derivatives root (used for CLI defaults)."""
    return Path(DATASETS_BY_NAME[dataset_name]["derivatives"])


def _requested_maps(maps_by_key, features, label, analysis, include_blur):
    """Filter discover_wscore_maps() output to the detection maps we score.

    Mirrors load_vertex_records: the requested white-surface regional maps plus,
    when ``include_blur``, the midthickness ``*blur`` companions. ``features``
    (base names, e.g. {"T1w","FLAIR"}) optionally restricts which are scored.
    """
    selected = []
    for (feature, map_label, smooth, map_analysis), hemispheres in maps_by_key.items():
        is_surface = map_label == label and map_analysis == analysis
        is_blur = (
            include_blur
            and feature.endswith("*blur")
            and map_label == "midthickness"
            and map_analysis == analysis
        )
        if not (is_surface or is_blur):
            continue
        if features is not None:
            base = feature[:-5] if feature.endswith("*blur") else feature
            if base not in features and feature not in features:
                continue
        if "L" not in hemispheres or "R" not in hemispheres:
            continue
        selected.append((feature, hemispheres))
    return selected


def per_subject_feature_auc(
    score_root,
    dataset_name: str,
    *,
    method: str = "wscore",
    features=None,
    label: str = "white",
    analysis: str = "regional",
    include_blur: bool = True,
    detection_tail: str = "both",
    lesion_threshold: float = 0.5,
    include_hippocampus: bool = True,
    hipp_density=None,
) -> pd.DataFrame:
    """Per-(subject, feature) lesion-detection AUROC for one cohort/root.

    With ``include_hippocampus`` (default True), hippocampal cavity records are
    scored against the hippocampal score maps ({method}_maps/hippocampus) as their
    OWN cells (feature ``hipp-<feature>``, lesion_type ``TLE``) so hippocampal /
    mesial-temporal detection carries non-zero weight in the objective -- TLE
    signal lives in the hippocampus, not only its cortical cavity projection.
    ``hipp_density=None`` uses the density configured for the cohort; a
    cavity/map vertex-count mismatch is skipped rather than crashing.

    Returns a long DataFrame with columns:
        dataset, subject, session, lesion_type, lesion_id, feature, auc, n_inside
    One row per (lesion record, feature map) whose inside/outside classes are both
    non-empty. ``method`` selects the ``<method>_maps`` score directory.
    """
    score_root = Path(score_root)
    lesion_dir = dataset_lesion_dir(dataset_name)
    if hipp_density is None:
        hipp_density = dataset_hipp_density(dataset_name)
    lesion_records = _LA.discover_lesions(lesion_dir, HEMI_VERTICES)

    rows = []
    for rec in lesion_records:
        subject = rec["subject"]
        session = rec.get("session", "ses-01")
        lesion = np.asarray(rec["lesion"], dtype=float)
        inside = lesion > lesion_threshold
        cortex_dir = score_root / subject / session / f"{method}_maps" / "cortex"
        if not cortex_dir.exists():
            continue
        maps_by_key = _LA.discover_wscore_maps(cortex_dir)
        for feature, hemispheres in _requested_maps(
            maps_by_key, features, label, analysis, include_blur
        ):
            data = np.concatenate(
                [
                    _LA.read_gifti_vector(hemispheres["L"], HEMI_VERTICES),
                    _LA.read_gifti_vector(hemispheres["R"], HEMI_VERTICES),
                ]
            ).astype(float)
            if feature in _LA.SIGN_FLIPPED_FEATURES:
                data = -data
            scores = _LA.detection_scores_from_values(data, detection_tail)
            gt_inside = int(inside.sum())          # lesion size in the GT (pre-masking)
            valid = np.isfinite(scores) & np.isfinite(lesion)
            labels = inside[valid]
            n_in, n_out = int(labels.sum()), int((~labels).sum())
            auc = _LA.auc_from_scores(labels, scores[valid])
            masked_out = False
            if not np.isfinite(auc):
                # FIXED-lesion-set fairness: a REAL lesion whose inside-vertices were
                # all masked to NaN (e.g. by the GP prediction-uncertainty filter) is
                # a MISS, not an absent map -- score AUC=0.5 rather than dropping the
                # row, so an arm cannot "win" by shrinking the evaluation set. A
                # genuinely absent map / no-lesion record (gt_inside==0) or one with
                # no valid outside vertices is still skipped.
                if gt_inside > 0 and n_out > 0 and n_in == 0:
                    auc, masked_out = 0.5, True
                else:
                    continue
            rows.append(
                {
                    "dataset": dataset_name,
                    "subject": subject,
                    "session": session,
                    "lesion_type": rec.get("lesion_type", "lesion"),
                    "lesion_id": rec.get("lesion_id", subject),
                    "feature": feature,
                    "auc": float(auc),
                    "n_inside": gt_inside if masked_out else int(labels.sum()),
                    "masked_out": masked_out,
                }
            )

    # --- Hippocampal channel: score hippocampal cavity GT against the hippocampal
    # score maps as its own feature/structure cells (hipp-<feature>, TLE). ---
    if include_hippocampus:
        hipp_cavities = _LA.discover_hipp_cavities(lesion_dir, density=hipp_density)
        for subject, session in sorted({(s, ses) for (s, ses, _h) in hipp_cavities}):
            hipp_dir = score_root / subject / session / f"{method}_maps" / "hippocampus"
            if not hipp_dir.exists():
                continue
            cav_l = hipp_cavities.get((subject, session, "L"))
            cav_r = hipp_cavities.get((subject, session, "R"))
            left = _LA.read_gifti_vector_any(cav_l) if cav_l else None
            right = _LA.read_gifti_vector_any(cav_r) if cav_r else None
            hipp_n = left.size if left is not None else (right.size if right is not None else None)
            if hipp_n is None:
                continue
            if left is None:
                left = np.zeros(hipp_n, dtype=float)
            if right is None:
                right = np.zeros(hipp_n, dtype=float)
            if left.size != hipp_n or right.size != hipp_n:
                continue                                   # density/mesh mismatch
            hlesion = np.concatenate([left, right]).astype(float)
            hinside = hlesion > lesion_threshold
            gt_inside = int(hinside.sum())
            if gt_inside == 0:
                continue
            for (feature, map_analysis), hemis in _LA.discover_hipp_wscore_maps(
                hipp_dir, density=hipp_density
            ).items():
                if map_analysis != analysis or "L" not in hemis or "R" not in hemis:
                    continue
                if features is not None and feature not in features:
                    continue
                if not include_blur and str(feature).endswith("*blur"):
                    continue
                data = np.concatenate([
                    _LA.read_gifti_vector_any(hemis["L"]),
                    _LA.read_gifti_vector_any(hemis["R"]),
                ]).astype(float)
                if data.size != hlesion.size:
                    continue
                if feature in _LA.SIGN_FLIPPED_FEATURES:
                    data = -data
                scores = _LA.detection_scores_from_values(data, detection_tail)
                valid = np.isfinite(scores) & np.isfinite(hlesion)
                labels = hinside[valid]
                n_in, n_out = int(labels.sum()), int((~labels).sum())
                auc = _LA.auc_from_scores(labels, scores[valid])
                masked_out = False
                if not np.isfinite(auc):
                    if gt_inside > 0 and n_out > 0 and n_in == 0:
                        auc, masked_out = 0.5, True
                    else:
                        continue
                rows.append(
                    {
                        "dataset": dataset_name,
                        "subject": subject,
                        "session": session,
                        "lesion_type": "TLE",
                        "lesion_id": f"{subject}-hipp",
                        "feature": f"hipp-{feature}",
                        "auc": float(auc),
                        "n_inside": gt_inside if masked_out else int(labels.sum()),
                        "masked_out": masked_out,
                    }
                )

    return pd.DataFrame(
        rows,
        columns=[
            "dataset", "subject", "session", "lesion_type",
            "lesion_id", "feature", "auc", "n_inside", "masked_out",
        ],
    )


def _read_cortex_detection_vector(cortex_dir, feature, label, analysis, detection_tail):
    """Bilateral detection statistic (length 2*HEMI_VERTICES) for one subject +
    feature, or None if the map is absent. Applies the per-feature sign flip then
    the tail (|z| for 'both') exactly as the positives are computed."""
    cortex_dir = Path(cortex_dir)
    if not cortex_dir.exists():
        return None
    maps_by_key = _LA.discover_wscore_maps(cortex_dir)
    for feat, hemis in _requested_maps(maps_by_key, [feature], label, analysis, True):
        if feat != feature or "L" not in hemis or "R" not in hemis:
            continue
        data = np.concatenate([
            _LA.read_gifti_vector(hemis["L"], HEMI_VERTICES),
            _LA.read_gifti_vector(hemis["R"], HEMI_VERTICES),
        ]).astype(float)
        if feature in _LA.SIGN_FLIPPED_FEATURES:
            data = -data
        return _LA.detection_scores_from_values(data, detection_tail)
    return None


class _RegionMatchedReference:
    """Region-matched negative provider for the specificity objectives.

    Reads each reference subject's detection vectors (cortex + hippocampus) once
    per feature (lazily, cached), optionally masking that subject's OWN
    lesion/cavity vertices to NaN, and returns the pooled finite |z| at a given
    lesion FOOTPRINT (a boolean vertex mask). So the negatives are the SAME
    anatomical vertices as the patient lesion, taken across the reference subjects,
    excluding any reference vertex that is itself lesional there ("if there is no
    lesion there"). For healthy controls no masks are passed, so every footprint
    vertex is a healthy-tissue negative.
    """

    def __init__(self, root, subjects, *, method, label, analysis, include_blur,
                 detection_tail, hipp_density, cortex_lesion=None, hipp_cavity=None):
        self.root = Path(root)
        self.subjects = list(subjects or [])
        self.method = method
        self.label = label
        self.analysis = analysis
        self.include_blur = include_blur
        self.detection_tail = detection_tail
        self.hipp_density = hipp_density
        self.cortex_lesion = cortex_lesion or {}     # (pid,ses) -> bilateral lesion vec
        self.hipp_cavity = hipp_cavity or {}         # (pid,ses) -> {hemi: path}
        self._cortex = {}                            # feature -> (n, 2*HEMI) or None
        self._hipp = {}                              # feature -> (n, hlen) or None

    def _cortex_matrix(self, feature):
        if feature not in self._cortex:
            rows = []
            for (pid, ses) in self.subjects:
                cdir = self.root / pid / ses / f"{self.method}_maps" / "cortex"
                vec = _read_cortex_detection_vector(
                    cdir, feature, self.label, self.analysis, self.detection_tail)
                if vec is None:
                    continue
                les = self.cortex_lesion.get((pid, ses))
                if les is not None and np.asarray(les).size == vec.size:
                    vec = vec.copy()
                    vec[np.asarray(les, dtype=float) > 0.5] = np.nan   # drop own lesion
                rows.append(vec)
            self._cortex[feature] = np.vstack(rows) if rows else None
        return self._cortex[feature]

    def cortex_negatives(self, feature, inside):
        mat = self._cortex_matrix(feature)
        if mat is None or mat.shape[1] != inside.size:
            return np.empty(0)
        vals = mat[:, inside].ravel()
        return vals[np.isfinite(vals)]

    def _hipp_matrix(self, feature):
        if feature not in self._hipp:
            rows = []
            target_len = None
            for (pid, ses) in self.subjects:
                hdir = self.root / pid / ses / f"{self.method}_maps" / "hippocampus"
                if not hdir.exists():
                    continue
                if not self.include_blur and str(feature).endswith("*blur"):
                    continue
                hemis = None
                for (f, manalysis), h in _LA.discover_hipp_wscore_maps(
                    hdir, density=self.hipp_density
                ).items():
                    if f == feature and manalysis == self.analysis and "L" in h and "R" in h:
                        hemis = h
                        break
                if hemis is None:
                    continue
                vec = np.concatenate([
                    _LA.read_gifti_vector_any(hemis["L"]),
                    _LA.read_gifti_vector_any(hemis["R"]),
                ]).astype(float)
                if feature in _LA.SIGN_FLIPPED_FEATURES:
                    vec = -vec
                vec = _LA.detection_scores_from_values(vec, self.detection_tail)
                cav = self.hipp_cavity.get((pid, ses))
                if cav:
                    outside = _hipp_outside_mask(cav, vec.size)
                    if outside is None:
                        continue                      # can't align cavity -> skip subject
                    vec = vec.copy()
                    vec[~outside] = np.nan             # drop own cavity
                if target_len is None:
                    target_len = vec.size
                if vec.size != target_len:
                    continue                          # density mismatch -> skip subject
                rows.append(vec)
            self._hipp[feature] = np.vstack(rows) if rows else None
        return self._hipp[feature]

    def hipp_negatives(self, feature, hinside):
        mat = self._hipp_matrix(feature)
        if mat is None or mat.shape[1] != hinside.size:
            return np.empty(0)
        vals = mat[:, hinside].ravel()
        return vals[np.isfinite(vals)]


def per_patient_vs_control_auc(
    patient_root,
    dataset_name: str,
    *,
    control_root=None,
    control_subjects=None,
    method: str = "wscore",
    features=None,
    label: str = "white",
    analysis: str = "regional",
    include_blur: bool = True,
    include_hippocampus: bool = True,
    detection_tail: str = "both",
    lesion_threshold: float = 0.5,
    hipp_density=None,
) -> pd.DataFrame:
    """Specificity-aware, PARAMETER-FREE, REGION-MATCHED detection objective.

    Per (patient lesion record, feature): a Mann-Whitney AUROC of the lesion's
    detection statistic |z| (POSITIVES) against the |z| of the held-out CONTROLS at
    the SAME anatomical vertices as the lesion footprint (NEGATIVES) -- i.e. the
    homologous region in every control, not the whole brain. Controls have no
    lesion, so every footprint vertex is a healthy-tissue reading. This asks "is the
    lesion more deviant than the same region in healthy brains?", which controls for
    regional variation in |z| calibration. No threshold, no per-subject reduction;
    over-flagging that same region in controls raises the negatives and drops the
    AUROC. A hippocampal cavity scores its ``hipp-<feature>`` cell against the
    control hippocampus at the cavity vertices.
    """
    patient_root = Path(patient_root)
    control_root = Path(control_root) if control_root is not None else patient_root
    lesion_dir = dataset_lesion_dir(dataset_name)
    if hipp_density is None:
        hipp_density = dataset_hipp_density(dataset_name)

    ref = _RegionMatchedReference(
        control_root, control_subjects, method=method, label=label, analysis=analysis,
        include_blur=include_blur, detection_tail=detection_tail, hipp_density=hipp_density,
    )  # controls: no lesion masks -> every footprint vertex is a valid negative

    rows = []

    def _emit(subject, session, lesion_type, lesion_id, feature, pos, neg, gt_inside):
        if neg is None or neg.size == 0 or gt_inside <= 0:
            return                      # no region-matched negatives, or no GT lesion here
        masked_out = pos.size == 0
        if masked_out:
            # FIXED-lesion-set fairness (mirrors per_subject_feature_auc): a REAL
            # lesion whose inside-vertices were ALL masked to NaN (e.g. by the GP
            # uncertainty / control-correlation filter) is a MISS, not an absent
            # map -- score AUC=0.5 so an arm cannot win by shrinking the eval set.
            auc = 0.5
        else:
            allscores = np.concatenate([pos, neg])
            labels = np.concatenate([np.ones(pos.size, bool), np.zeros(neg.size, bool)])
            auc = _LA.auc_from_scores(labels, allscores)
            if not np.isfinite(auc):
                return
        rows.append({
            "dataset": dataset_name, "subject": subject, "session": session,
            "lesion_type": lesion_type, "lesion_id": lesion_id, "feature": feature,
            "auc": float(auc), "n_inside": int(gt_inside if masked_out else pos.size),
            "masked_out": bool(masked_out),
        })

    # --- cortical lesions: lesion-vertex |z| vs control |z| at the SAME vertices ---
    for rec in _LA.discover_lesions(lesion_dir, HEMI_VERTICES):
        subject = rec["subject"]; session = rec.get("session", "ses-01")
        lesion = np.asarray(rec["lesion"], dtype=float)
        inside = lesion > lesion_threshold
        gt_inside = int(inside.sum())          # GT lesion size (pre-masking)
        cortex_dir = patient_root / subject / session / f"{method}_maps" / "cortex"
        if not cortex_dir.exists():
            continue
        maps_by_key = _LA.discover_wscore_maps(cortex_dir)
        for feature, hemis in _requested_maps(maps_by_key, features, label, analysis, include_blur):
            data = np.concatenate([
                _LA.read_gifti_vector(hemis["L"], HEMI_VERTICES),
                _LA.read_gifti_vector(hemis["R"], HEMI_VERTICES),
            ]).astype(float)
            if feature in _LA.SIGN_FLIPPED_FEATURES:
                data = -data
            scores = _LA.detection_scores_from_values(data, detection_tail)
            valid = np.isfinite(scores) & np.isfinite(lesion)
            pos = scores[valid][inside[valid]]
            neg = ref.cortex_negatives(feature, inside)   # controls at the SAME footprint
            _emit(subject, session, rec.get("lesion_type", "lesion"),
                  rec.get("lesion_id", subject), feature, pos, neg, gt_inside)

    # --- hippocampal cavities: cavity-vertex |z| vs control hipp at the SAME verts ---
    if include_hippocampus:
        hipp_cavities = _LA.discover_hipp_cavities(lesion_dir, density=hipp_density)
        for subject, session in sorted({(s, ses) for (s, ses, _h) in hipp_cavities}):
            hipp_dir = patient_root / subject / session / f"{method}_maps" / "hippocampus"
            if not hipp_dir.exists():
                continue
            cav_l = hipp_cavities.get((subject, session, "L"))
            cav_r = hipp_cavities.get((subject, session, "R"))
            left = _LA.read_gifti_vector_any(cav_l) if cav_l else None
            right = _LA.read_gifti_vector_any(cav_r) if cav_r else None
            hipp_n = left.size if left is not None else (right.size if right is not None else None)
            if hipp_n is None:
                continue
            left = np.zeros(hipp_n, dtype=float) if left is None else left
            right = np.zeros(hipp_n, dtype=float) if right is None else right
            if left.size != hipp_n or right.size != hipp_n:
                continue
            hlesion = np.concatenate([left, right]).astype(float)
            hinside = hlesion > lesion_threshold
            gt_inside = int(hinside.sum())
            if gt_inside == 0:
                continue
            for (feature, map_analysis), hemis in _LA.discover_hipp_wscore_maps(
                hipp_dir, density=hipp_density
            ).items():
                if map_analysis != analysis or "L" not in hemis or "R" not in hemis:
                    continue
                if features is not None and feature not in features:
                    continue
                if not include_blur and str(feature).endswith("*blur"):
                    continue
                data = np.concatenate([
                    _LA.read_gifti_vector_any(hemis["L"]),
                    _LA.read_gifti_vector_any(hemis["R"]),
                ]).astype(float)
                if data.size != hlesion.size:
                    continue
                if feature in _LA.SIGN_FLIPPED_FEATURES:
                    data = -data
                scores = _LA.detection_scores_from_values(data, detection_tail)
                valid = np.isfinite(scores) & np.isfinite(hlesion)
                pos = scores[valid][hinside[valid]]
                neg = ref.hipp_negatives(feature, hinside)   # control hipp at the cavity
                _emit(subject, session, "TLE", f"{subject}-hipp", f"hipp-{feature}",
                      pos, neg, gt_inside)

    return pd.DataFrame(
        rows,
        columns=[
            "dataset", "subject", "session", "lesion_type",
            "lesion_id", "feature", "auc", "n_inside", "masked_out",
        ],
    )


def _hipp_outside_mask(cavity_by_hemi, total_size):
    """Bilateral boolean 'non-lesional hippocampus' mask of length ``total_size``.

    Returns all-True when the subject has NO cavity (every hipp vertex is genuinely
    non-lesional). Returns ``None`` when a cavity IS present but cannot be aligned
    to ``total_size`` (unreadable / density mismatch) -- the caller must then SKIP
    that contribution rather than pool the resected (abnormal) hippocampus into the
    'non-lesional' negative class."""
    if not cavity_by_hemi:
        return np.ones(total_size, dtype=bool)
    left = _LA.read_gifti_vector_any(cavity_by_hemi["L"]) if "L" in cavity_by_hemi else None
    right = _LA.read_gifti_vector_any(cavity_by_hemi["R"]) if "R" in cavity_by_hemi else None
    n = left.size if left is not None else (right.size if right is not None else None)
    if n is None:
        return None                              # cavity present but unreadable -> skip
    left = np.zeros(n) if left is None else left
    right = np.zeros(n) if right is None else right
    cav = np.concatenate([left, right])
    if cav.size != total_size:
        return None                              # can't align the cavity -> skip (no contaminate)
    return ~(cav > 0.5)


def per_patient_vs_other_disease_auc(
    patient_root,
    dataset_name: str,
    *,
    method: str = "wscore",
    features=None,
    label: str = "white",
    analysis: str = "regional",
    include_blur: bool = True,
    include_hippocampus: bool = True,
    detection_tail: str = "both",
    lesion_threshold: float = 0.5,
    hipp_density=None,
) -> pd.DataFrame:
    """DISEASE-CONTROL specificity objective (parameter-free, REGION-MATCHED).

    Per (patient lesion record, feature): a Mann-Whitney AUROC of the lesion's
    detection statistic |z| (POSITIVES) against the |z| of the OTHER disease's
    patients at the SAME anatomical vertices as the lesion footprint (NEGATIVES),
    EXCLUDING any other-disease vertex that is itself lesional there. TLE lesions
    are ranked against FCD patients' tissue at the same location and vice versa --
    "is the lesion more deviant than that region in a different patient
    population?", a stronger specificity bar than healthy controls. Both sides come
    from the same analysis dir, so the |z| scales match; no held-out controls are
    needed.
    """
    patient_root = Path(patient_root)
    lesion_dir = dataset_lesion_dir(dataset_name)
    if hipp_density is None:
        hipp_density = dataset_hipp_density(dataset_name)
    records = _LA.discover_lesions(lesion_dir, HEMI_VERTICES)
    hipp_cavities = (
        _LA.discover_hipp_cavities(lesion_dir, density=hipp_density)
        if include_hippocampus else {}
    )

    # Cortical lesion mask + disease per (subject, session).
    lesion_by_subject = {}
    disease_by_subject = {}
    for rec in records:
        key = (rec["subject"], rec.get("session", "ses-01"))
        lesion_by_subject[key] = np.asarray(rec["lesion"], dtype=float)
        disease_by_subject[key] = rec.get("lesion_type", "lesion")
    hipp_cav_by_subject = {}
    for (pid, ses, hemi), path in hipp_cavities.items():
        hipp_cav_by_subject.setdefault((pid, ses), {})[hemi] = path
        disease_by_subject.setdefault((pid, ses), "TLE")

    subjects_by_disease = {}
    for key, disease in disease_by_subject.items():
        subjects_by_disease.setdefault(disease, []).append(key)

    def _other(disease):
        return "FCD" if disease == "TLE" else "TLE"

    # A region-matched reference per POSITIVE disease D, built over the OTHER
    # disease's patients with THEIR own lesion/cavity vertices masked out (so a
    # footprint vertex that is lesional in the other patient is not a negative).
    emitted = {rec.get("lesion_type", "lesion") for rec in records}
    if include_hippocampus and hipp_cavities:
        emitted.add("TLE")
    refs = {}
    for disease in emitted | {"FCD", "TLE"}:
        refs[disease] = _RegionMatchedReference(
            patient_root, subjects_by_disease.get(_other(disease), []),
            method=method, label=label, analysis=analysis, include_blur=include_blur,
            detection_tail=detection_tail, hipp_density=hipp_density,
            cortex_lesion=lesion_by_subject, hipp_cavity=hipp_cav_by_subject,
        )
    # A single-disease cohort has no OTHER-disease patients -> its lesions get no
    # negatives and would be dropped SILENTLY. Surface it loudly.
    for disease in emitted:
        if not subjects_by_disease.get(_other(disease)):
            warnings.warn(
                f"{dataset_name}: disease_control has NO '{_other(disease)}' patients "
                f"to form a negative pool for '{disease}' lesions -- those lesions are "
                "dropped from this objective (single-disease cohort?).",
                RuntimeWarning, stacklevel=2,
            )

    rows = []

    def _emit(subject, session, lesion_type, lesion_id, feature, pos, neg, gt_inside):
        if neg is None or neg.size == 0 or gt_inside <= 0:
            return
        masked_out = pos.size == 0
        if masked_out:
            auc = 0.5
        else:
            allscores = np.concatenate([pos, neg])
            labels = np.concatenate([np.ones(pos.size, bool), np.zeros(neg.size, bool)])
            auc = _LA.auc_from_scores(labels, allscores)
            if not np.isfinite(auc):
                return
        rows.append({
            "dataset": dataset_name, "subject": subject, "session": session,
            "lesion_type": lesion_type, "lesion_id": lesion_id, "feature": feature,
            "auc": float(auc), "n_inside": int(gt_inside if masked_out else pos.size),
            "masked_out": bool(masked_out),
        })

    # --- cortical lesions vs the OTHER disease at the SAME footprint ---
    for rec in records:
        subject = rec["subject"]; session = rec.get("session", "ses-01")
        lesion = np.asarray(rec["lesion"], dtype=float)
        inside = lesion > lesion_threshold
        gt_inside = int(inside.sum())
        lesion_type = rec.get("lesion_type", "lesion")
        cortex_dir = patient_root / subject / session / f"{method}_maps" / "cortex"
        if not cortex_dir.exists():
            continue
        maps_by_key = _LA.discover_wscore_maps(cortex_dir)
        for feature, hemis in _requested_maps(maps_by_key, features, label, analysis, include_blur):
            data = np.concatenate([
                _LA.read_gifti_vector(hemis["L"], HEMI_VERTICES),
                _LA.read_gifti_vector(hemis["R"], HEMI_VERTICES),
            ]).astype(float)
            if feature in _LA.SIGN_FLIPPED_FEATURES:
                data = -data
            scores = _LA.detection_scores_from_values(data, detection_tail)
            valid = np.isfinite(scores) & np.isfinite(lesion)
            pos = scores[valid][inside[valid]]
            neg = refs[lesion_type].cortex_negatives(feature, inside)
            _emit(subject, session, lesion_type, rec.get("lesion_id", subject),
                  feature, pos, neg, gt_inside)

    # --- hippocampal cavities (TLE) vs FCD patients at the SAME cavity verts ---
    if include_hippocampus:
        for subject, session in sorted({(s, ses) for (s, ses, _h) in hipp_cavities}):
            hipp_dir = patient_root / subject / session / f"{method}_maps" / "hippocampus"
            if not hipp_dir.exists():
                continue
            cav_l = hipp_cavities.get((subject, session, "L"))
            cav_r = hipp_cavities.get((subject, session, "R"))
            left = _LA.read_gifti_vector_any(cav_l) if cav_l else None
            right = _LA.read_gifti_vector_any(cav_r) if cav_r else None
            hipp_n = left.size if left is not None else (right.size if right is not None else None)
            if hipp_n is None:
                continue
            left = np.zeros(hipp_n, dtype=float) if left is None else left
            right = np.zeros(hipp_n, dtype=float) if right is None else right
            if left.size != hipp_n or right.size != hipp_n:
                continue
            hlesion = np.concatenate([left, right]).astype(float)
            hinside = hlesion > lesion_threshold
            gt_inside = int(hinside.sum())
            if gt_inside == 0:
                continue
            for (feature, map_analysis), hemis in _LA.discover_hipp_wscore_maps(
                hipp_dir, density=hipp_density
            ).items():
                if map_analysis != analysis or "L" not in hemis or "R" not in hemis:
                    continue
                if features is not None and feature not in features:
                    continue
                if not include_blur and str(feature).endswith("*blur"):
                    continue
                data = np.concatenate([
                    _LA.read_gifti_vector_any(hemis["L"]),
                    _LA.read_gifti_vector_any(hemis["R"]),
                ]).astype(float)
                if data.size != hlesion.size:
                    continue
                if feature in _LA.SIGN_FLIPPED_FEATURES:
                    data = -data
                scores = _LA.detection_scores_from_values(data, detection_tail)
                valid = np.isfinite(scores) & np.isfinite(hlesion)
                pos = scores[valid][hinside[valid]]
                neg = refs["TLE"].hipp_negatives(feature, hinside)
                _emit(subject, session, "TLE", f"{subject}-hipp", f"hipp-{feature}",
                      pos, neg, gt_inside)

    return pd.DataFrame(
        rows,
        columns=[
            "dataset", "subject", "session", "lesion_type",
            "lesion_id", "feature", "auc", "n_inside", "masked_out",
        ],
    )


def _per_subject_pooled_auc(
    score_root,
    dataset_name,
    *,
    method,
    features,
    label,
    analysis,
    include_blur,
    detection_tail,
    lesion_threshold,
) -> pd.DataFrame:
    """One AUC per lesion record, pooling all requested feature vertices together.

    Columns: dataset, subject, lesion_type, lesion_id, auc.
    """
    score_root = Path(score_root)
    lesion_dir = dataset_lesion_dir(dataset_name)
    lesion_records = _LA.discover_lesions(lesion_dir, HEMI_VERTICES)

    rows = []
    for rec in lesion_records:
        subject = rec["subject"]
        session = rec.get("session", "ses-01")
        lesion = np.asarray(rec["lesion"], dtype=float)
        inside = lesion > lesion_threshold
        cortex_dir = score_root / subject / session / f"{method}_maps" / "cortex"
        if not cortex_dir.exists():
            continue
        maps_by_key = _LA.discover_wscore_maps(cortex_dir)
        labels_all, scores_all = [], []
        for feature, hemispheres in _requested_maps(
            maps_by_key, features, label, analysis, include_blur
        ):
            data = np.concatenate(
                [
                    _LA.read_gifti_vector(hemispheres["L"], HEMI_VERTICES),
                    _LA.read_gifti_vector(hemispheres["R"], HEMI_VERTICES),
                ]
            ).astype(float)
            if feature in _LA.SIGN_FLIPPED_FEATURES:
                data = -data
            scores = _LA.detection_scores_from_values(data, detection_tail)
            valid = np.isfinite(scores) & np.isfinite(lesion)
            if not np.any(valid):
                continue
            scores_all.append(scores[valid])
            labels_all.append(inside[valid])
        if not labels_all:
            continue
        auc = _LA.auc_from_scores(np.concatenate(labels_all), np.concatenate(scores_all))
        if not np.isfinite(auc):
            continue
        rows.append(
            {
                "dataset": dataset_name,
                "subject": subject,
                "lesion_type": rec.get("lesion_type", "lesion"),
                "lesion_id": rec.get("lesion_id", subject),
                "auc": float(auc),
            }
        )
    return pd.DataFrame(
        rows, columns=["dataset", "subject", "lesion_type", "lesion_id", "auc"]
    )


AGGREGATIONS = ("mean_subject_feature", "per_subject_pooled", "best_feature_per_subject")


def macro_auroc_for_root(
    score_root,
    dataset_name: str,
    *,
    method: str = "wscore",
    aggregation: str = "mean_subject_feature",
    features=None,
    return_frame: bool = False,
    **kwargs,
):
    """Per-subject macro AUROC for one cohort/root under one of the AGGREGATIONS.

    * ``mean_subject_feature``   -- mean over all (subject, feature) AUCs
      (matches benchmark_processing_method_matrix's groupby-mean convention).
    * ``per_subject_pooled``     -- pool each subject's feature vertices into one
      AUC, then mean over subjects.
    * ``best_feature_per_subject`` -- per subject take the best feature AUC, then
      mean over subjects (optimistic "did any feature detect it").

    Returns the scalar, or ``(scalar, per_subject_feature_frame)`` if
    ``return_frame`` (the returned frame is always the long per-(subject,feature)
    table, for diagnostics, regardless of aggregation).
    """
    if aggregation not in AGGREGATIONS:
        raise ValueError(f"aggregation must be one of {AGGREGATIONS}; got {aggregation!r}")

    frame = per_subject_feature_auc(
        score_root, dataset_name, method=method, features=features, **kwargs
    )

    if aggregation == "per_subject_pooled":
        pooled = _per_subject_pooled_auc(
            score_root, dataset_name, method=method, features=features, **kwargs
        )
        scalar = float(pooled["auc"].mean()) if not pooled.empty else float("nan")
    elif aggregation == "best_feature_per_subject":
        if frame.empty:
            scalar = float("nan")
        else:
            best = frame.groupby(["subject", "lesion_id"])["auc"].max()
            scalar = float(best.mean())
    else:  # mean_subject_feature
        scalar = float(frame["auc"].mean()) if not frame.empty else float("nan")

    return (scalar, frame) if return_frame else scalar


def pooled_macro_auroc(
    roots_by_dataset: dict,
    *,
    method: str = "wscore",
    aggregation: str = "mean_subject_feature",
    features=None,
    return_details: bool = False,
    **kwargs,
):
    """Single pooled per-subject macro AUROC across all provided cohorts.

    ``roots_by_dataset`` maps dataset name (e.g. "MICs", "NOEL") to its score root.
    Pooling concatenates the per-cohort per-subject(-feature) rows and takes a
    single mean, so cohorts contribute in proportion to patient count (the
    intended "single pooled pipeline" objective) rather than being averaged as
    two equally-weighted cohort means.
    """
    if aggregation == "per_subject_pooled":
        parts = [
            _per_subject_pooled_auc(
                root, name, method=method, features=features, **kwargs
            )
            for name, root in roots_by_dataset.items()
        ]
        detail = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        scalar = float(detail["auc"].mean()) if not detail.empty else float("nan")
    else:
        parts = [
            per_subject_feature_auc(
                root, name, method=method, features=features, **kwargs
            )
            for name, root in roots_by_dataset.items()
        ]
        detail = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        if detail.empty:
            scalar = float("nan")
        elif aggregation == "best_feature_per_subject":
            best = detail.groupby(["dataset", "subject", "lesion_id"])["auc"].max()
            scalar = float(best.mean())
        else:  # mean_subject_feature
            scalar = float(detail["auc"].mean())

    return (scalar, detail) if return_details else scalar


def balanced_macro_auroc(detail, features=None, group_cols=("dataset", "lesion_type")):
    """Macro AUROC that weights each (dataset x disease) cell EQUALLY.

    Takes a per-(subject, feature) AUC frame (as returned by
    :func:`per_subject_feature_auc` / :func:`pooled_macro_auroc(..., return_details=True)`),
    averages the AUCs within each ``group_cols`` cell, then averages those cell
    means -- so a larger cohort or a more-common disease does not dominate the
    objective. ``features`` optionally restricts to a set of feature names
    (e.g. ``{"T1w", "FLAIR"}``). Returns NaN if nothing matches.
    """
    if detail is None or len(detail) == 0:
        return float("nan")
    sub = detail
    if features is not None:
        sub = sub[sub["feature"].isin(set(features))]
    if len(sub) == 0:
        return float("nan")
    cols = [c for c in group_cols if c in sub.columns]
    if not cols:
        return float(sub["auc"].mean())
    return float(sub.groupby(cols)["auc"].mean().mean())


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mics-root", type=Path, default=None)
    parser.add_argument("--noel-root", type=Path, default=None)
    parser.add_argument("--method", default="wscore", choices=["wscore", "zscore"])
    parser.add_argument("--aggregation", default="mean_subject_feature", choices=list(AGGREGATIONS))
    parser.add_argument(
        "--features", default=None,
        help="Comma-separated base feature whitelist (e.g. T1w,FLAIR,ADC,FA). Default: all present.",
    )
    parser.add_argument("--no-blur", action="store_true", help="Exclude *blur companion maps.")
    return parser.parse_args()


def main():
    args = _parse_args()
    features = (
        {f.strip() for f in args.features.split(",") if f.strip()}
        if args.features
        else None
    )
    roots = {}
    if args.mics_root is not None:
        roots["MICs"] = args.mics_root
    if args.noel_root is not None:
        roots["NOEL"] = args.noel_root
    if not roots:
        raise SystemExit("Provide at least one of --mics-root / --noel-root")

    scalar, detail = pooled_macro_auroc(
        roots,
        method=args.method,
        aggregation=args.aggregation,
        features=features,
        include_blur=not args.no_blur,
        return_details=True,
    )
    print(f"Pooled per-subject macro AUROC ({args.aggregation}, method={args.method}): {scalar:.4f}")
    if not detail.empty and "feature" in detail.columns:
        print("\nPer-feature mean AUC (pooled cohorts):")
        by_feat = detail.groupby("feature")["auc"].agg(["mean", "count"]).sort_values("mean", ascending=False)
        print(by_feat.to_string())
        print("\nPer-lesion-type mean AUC:")
        by_type = detail.groupby("lesion_type")["auc"].agg(["mean", "count"])
        print(by_type.to_string())
    elif not detail.empty:
        print(f"\nSubjects scored: {detail['subject'].nunique()}")


if __name__ == "__main__":
    main()
