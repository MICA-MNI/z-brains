#!/usr/bin/env python3
"""Compare cortical W-scores inside vs outside lesion masks.

This script expects lesion CSVs with one value per fsLR-32k vertex, ordered as
left hemisphere vertices followed by right hemisphere vertices. It pairs those
masks with each subject's cortical W-score GIFTI maps and writes summary CSVs
plus feature-wise inside/outside plots.
"""

import argparse
import itertools
import os
import re
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Avenir Next'
plt.rcParams['font.weight'] = 'normal'
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['axes.titleweight'] = 'normal'
DEFAULT_LESION_DIR = Path("results/lesiondata/NOEL")
DEFAULT_WSCORE_ROOT = Path("/host/verges/tank/data/ian/BIDS_NOEL/derivatives/zbrains_WB")
DEFAULT_OUTPUT_DIR = Path("results/noel_lesion_wscore_overlap")
DEFAULT_MICS_DERIVATIVES_ROOT = Path("/host/verges/tank/data/ian/BIDS_NOEL/derivatives")
DEFAULT_HEMI_VERTICES = 32492
# Hippunfold surface density whose W-score maps and cavity labels are appended to
# the cortical data. The NOEL optimization now uses HippUnfold v2 den-8k.
DEFAULT_HIPP_DENSITY = "8k"
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SPHERE_L = REPO_ROOT / "src" / "zbrains" / "data" / "fsLR-32k.L.sphere.reg.surf.gii"
DEFAULT_SPHERE_R = REPO_ROOT / "src" / "zbrains" / "data" / "fsLR-32k.R.sphere.reg.surf.gii"

FEATURE_ORDER = [
    "T1w",
    "T1w*blur",
    "FLAIR",
    "FLAIR*blur",
    "ADC",
    "qT1",
    "qT1*blur",
    "thickness",
]
SIGN_FLIPPED_FEATURES = {"T1w*blur"}
FEATURE_ABBREVIATIONS = {
    "ALL_METRICS": "All",
    "T1w": "T1",
    "T1w*blur": "T1b",
    "FLAIR": "FL",
    "FLAIR*blur": "FLb",
    "ADC": "ADC",
    "FA": "FA",
    "qT1": "qT1",
    "qT1*blur": "qT1b",
    "thickness": "CT",
}

WSCORE_PATTERN = re.compile(
    r"(?P<subject>sub-[^_]+)_(?P<session>ses-[^_]+)_hemi-(?P<hemi>[LR])_surf-fsLR-32k_"
    r"label-(?P<label>[^_]+)_feature-(?P<feature>.+?)_smooth-(?P<smooth>[^_]+)_"
    r"analysis-(?P<analysis>[^.]+)\.func\.gii$"
)
TLE_LESION_PATTERN = re.compile(
    r"(?P<subject>sub-[^_]+)_(?P<session>ses-[^_]+)_hemi-(?P<hemi>[LR])_"
    r"surf-fsLR-32k_label-[^_]+_cavity\.func\.gii$"
)
FCD_LESION_PATTERN = re.compile(
    r"(?P<subject>sub-[^_]+)_(?P<session>ses-[^_]+)_hemi-(?P<hemi>[LR])_"
    r"surf-fsLR-32k_label-(?P<label>[^_]+)_lesion\.func\.gii$"
)
# Hippocampal W-score maps written by the pipeline, e.g.
# sub-PX_ses-01_hemi-L_den-8k_label-hipp_midthickness_feature-T1w_smooth-2mm_analysis-regional.func.gii
HIPP_WSCORE_PATTERN = re.compile(
    r"(?P<subject>sub-[^_]+)_(?P<session>ses-[^_]+)_hemi-(?P<hemi>[LR])_"
    r"den-(?P<den>[^_]+)_label-hipp_midthickness_feature-(?P<feature>.+?)_"
    r"smooth-(?P<smooth>[^_]+)_analysis-(?P<analysis>[^.]+)\.func\.gii$"
)
# Hippocampal cavity ground truth produced by make_hippocampal_cavity_labels.py, e.g.
# sub-PX_ses-01_hemi-R_den-8k_label-hipp_cavity.func.gii
HIPP_CAVITY_PATTERN = re.compile(
    r"(?P<subject>sub-[^_]+)_(?P<session>ses-[^_]+)_hemi-(?P<hemi>[LR])_"
    r"den-(?P<den>[^_]+)_label-hipp_cavity\.func\.gii$"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot each cortical feature's W-score inside vs outside lesion masks."
    )
    parser.add_argument("--lesion-dir", type=Path, default=DEFAULT_LESION_DIR)
    parser.add_argument("--wscore-root", type=Path, default=DEFAULT_WSCORE_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--compare-methods",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Run the lesion-overlap analysis for multiple processed z-brains roots. "
            "Use --no-compare-methods to analyze only --wscore-root."
        ),
    )
    parser.add_argument(
        "--method-root",
        action="append",
        default=None,
        metavar="NAME=PATH",
        help=(
            "Processed z-brains output root to include in method comparison. "
            "Can be repeated. Defaults to WB and RAVEL."
        ),
    )
    parser.add_argument("--session", default="ses-01")
    parser.add_argument("--hemi-vertices", type=int, default=DEFAULT_HEMI_VERTICES)
    parser.add_argument("--lesion-threshold", type=float, default=0.5)
    parser.add_argument("--label", default="white", help="Cortical surface label to analyze.")
    parser.add_argument("--analysis", default="regional", help="W-score analysis type to analyze.")
    parser.add_argument(
        "--include-blur",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Also include regional *blur maps, which are written on the midthickness surface.",
    )
    parser.add_argument(
        "--include-hippocampus",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Append hippocampal W-score vertices (and their cavity ground truth) to the "
            "cortical data so detection metrics span cortex + hippocampus."
        ),
    )
    parser.add_argument(
        "--hipp-density",
        default=DEFAULT_HIPP_DENSITY,
        help="Hippunfold surface density of the appended hippocampal maps (e.g. '8k' or '0p5mm').",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR-corrected significance threshold for inside-vs-outside tests.",
    )
    parser.add_argument(
        "--max-combination-size",
        type=int,
        default=3,
        help="Largest feature combination size to test. Set below 2 to skip combinations.",
    )
    parser.add_argument(
        "--run-combinations",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Test feature combinations and optional permutation/spin nulls.",
    )
    parser.add_argument(
        "--min-combination-subjects",
        type=int,
        default=5,
        help="Minimum complete subjects required to test a feature combination.",
    )
    parser.add_argument(
        "--combination-operators",
        default="+",
        help="Comma-separated operators to test between features in combinations.",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=0,
        help="Number of max-statistic permutations for combination-level procedure correction.",
    )
    parser.add_argument(
        "--permutation-seed",
        type=int,
        default=12345,
        help="Random seed for permutation testing.",
    )
    parser.add_argument(
        "--n-spins",
        type=int,
        default=1000,
        help="Number of fsLR spherical spins for combination-level spatial correction.",
    )
    parser.add_argument(
        "--spin-seed",
        type=int,
        default=54321,
        help="Random seed for fsLR spherical spins.",
    )
    parser.add_argument(
        "--detection-threshold",
        type=float,
        default=1.96,
        help="W-score threshold used to count lesion detections and false-positive vertices.",
    )
    parser.add_argument(
        "--detection-tail",
        choices=["positive", "both", "absolute", "negative"],
        default="both",
        help=(
            "Tail used for detection counts. 'positive' treats values above threshold as detections; "
            "'both' counts either sign; 'negative' counts values below -threshold. "
            "'absolute' is accepted as an alias for 'both'."
        ),
    )
    parser.add_argument(
        "--run-roc",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write ROC/AUC CSVs and ROC plots by sweeping W-score detection thresholds.",
    )
    parser.add_argument(
        "--roc-threshold-min",
        type=float,
        default=0.0,
        help="Minimum W-score detection threshold for ROC curve points.",
    )
    parser.add_argument(
        "--roc-threshold-max",
        type=float,
        default=6.0,
        help="Maximum W-score detection threshold for ROC curve points.",
    )
    parser.add_argument(
        "--roc-threshold-steps",
        type=int,
        default=121,
        help="Number of W-score thresholds sampled for ROC curve points.",
    )
    parser.add_argument(
        "--roc-max-curves",
        type=int,
        default=12,
        help="Maximum number of top-AUC ROC curves to draw per plot.",
    )
    parser.add_argument("--sphere-l", type=Path, default=DEFAULT_SPHERE_L)
    parser.add_argument("--sphere-r", type=Path, default=DEFAULT_SPHERE_R)
    args = parser.parse_args()
    if args.detection_tail == "absolute":
        args.detection_tail = "both"
    return args


def safe_method_name(name):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name).strip())
    return safe.strip("_") or "method"


def parse_method_roots(args):
    if not args.compare_methods:
        return [("single", args.wscore_root)]

    if not args.method_root:
        return discover_default_method_roots()

    roots = []
    for spec in args.method_root:
        if "=" not in spec:
            raise ValueError(f"Invalid --method-root value {spec!r}; expected NAME=PATH")
        name, path = spec.split("=", 1)
        name = name.strip()
        if not name:
            raise ValueError(f"Invalid --method-root value {spec!r}; method name is empty")
        roots.append((name, Path(path).expanduser()))
    return roots


def method_name_from_root(path):
    name = Path(path).name
    if name.startswith("zbrains_"):
        name = name[len("zbrains_"):]
    return name


def discover_default_method_roots():
    roots = []
    for pattern in ("zbrains_WB*", "zbrains_RAVEL*"):
        for path in sorted(DEFAULT_MICS_DERIVATIVES_ROOT.glob(pattern)):
            if path.is_dir():
                roots.append((method_name_from_root(path), path))
    if roots:
        return roots
    return [
        ("WB", DEFAULT_MICS_DERIVATIVES_ROOT / "zbrains_WB"),
        ("RAVEL", DEFAULT_MICS_DERIVATIVES_ROOT / "zbrains_RAVEL"),
    ]


def ordered_features(values):
    values = list(dict.fromkeys(values))
    preferred = [feature for feature in FEATURE_ORDER if feature in values]
    remaining = sorted(feature for feature in values if feature not in FEATURE_ORDER)
    return preferred + remaining


def read_lesion_csv(path, hemi_vertices):
    values = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            value = line.strip().strip('"')
            values.append(np.nan if value == "" else float(value))

    lesion = np.asarray(values, dtype=np.float32)
    expected = 2 * hemi_vertices
    if lesion.size != expected:
        raise ValueError(f"{path} has {lesion.size} rows; expected {expected}")
    return lesion


def read_lesion_gifti(path, hemi_vertices):
    data = read_gifti_vector(path, hemi_vertices)
    return np.asarray(data, dtype=np.float32)


def infer_csv_lesion_type(path):
    parts = [part.lower() for part in Path(path).parts]
    if "fcd" in parts:
        return "FCD"
    if "tle" in parts:
        return "TLE"
    return "lesion"


def resolve_lesion_subject(subject, path):
    """Map label-space subject IDs onto processed z-brains subject IDs."""
    parts = {part.lower() for part in Path(path).parts}
    match = re.fullmatch(r"sub-PX(?P<number>\d+)", subject)
    if {"noel", "fcd"}.issubset(parts) and match:
        return f"sub-FCDPX{match.group('number')}"
    return subject


def discover_lesions(lesion_dir, hemi_vertices):
    """Return bilateral lesion vectors from CSV masks and per-hemi GIFTIs."""
    lesion_dir = Path(lesion_dir)
    lesion_records = []

    for path in sorted(lesion_dir.glob("**/sub-*/sub-*_lesion_32k.csv")):
        label_subject = path.parent.name
        subject = resolve_lesion_subject(label_subject, path)
        lesion_type = infer_csv_lesion_type(path)
        lesion_records.append(
            {
                "subject": subject,
                "label_subject": label_subject,
                "session": "ses-01",
                "lesion_type": lesion_type,
                "lesion_id": f"{subject}_ses-01_{lesion_type}",
                "path": path,
                "lesion": read_lesion_csv(path, hemi_vertices),
            }
        )

    fcd_by_subject = {}
    label_priority = {"midthickness": 0, "whiteormidthickness": 1}
    for path in sorted(lesion_dir.glob("**/*_lesion.func.gii")):
        match = FCD_LESION_PATTERN.search(path.name)
        if not match:
            continue
        metadata = match.groupdict()
        label = metadata["label"].lower()
        key = (metadata["subject"], metadata["session"])
        existing = fcd_by_subject.get(key)
        priority = label_priority.get(label, 99)
        if existing is None or priority < existing["priority"]:
            fcd_by_subject[key] = {
                "priority": priority,
                "label": metadata["label"],
                "hemispheres": {},
            }
        elif priority > existing["priority"]:
            continue
        fcd_by_subject[key]["hemispheres"][metadata["hemi"]] = path

    for (label_subject, session), lesion_data in sorted(fcd_by_subject.items()):
        hemispheres = lesion_data["hemispheres"]
        subject = resolve_lesion_subject(label_subject, next(iter(hemispheres.values())))
        lesion_l = np.zeros(hemi_vertices, dtype=np.float32)
        lesion_r = np.zeros(hemi_vertices, dtype=np.float32)
        paths = []
        if "L" in hemispheres:
            lesion_l = read_lesion_gifti(hemispheres["L"], hemi_vertices)
            paths.append(str(hemispheres["L"]))
        if "R" in hemispheres:
            lesion_r = read_lesion_gifti(hemispheres["R"], hemi_vertices)
            paths.append(str(hemispheres["R"]))
        lesion_records.append(
            {
                "subject": subject,
                "label_subject": label_subject,
                "session": session,
                "lesion_type": "FCD",
                "lesion_id": f"{subject}_{session}_FCD",
                "path": ";".join(paths),
                "lesion": np.concatenate([lesion_l, lesion_r]),
            }
        )

    tle_by_subject = {}
    for path in sorted(lesion_dir.glob("**/*_cavity.func.gii")):
        match = TLE_LESION_PATTERN.search(path.name)
        if not match:
            continue
        metadata = match.groupdict()
        subject = resolve_lesion_subject(metadata["subject"], path)
        key = (subject, metadata["subject"], metadata["session"])
        tle_by_subject.setdefault(key, {})[metadata["hemi"]] = path

    for (subject, label_subject, session), hemispheres in sorted(tle_by_subject.items()):
        lesion_l = np.zeros(hemi_vertices, dtype=np.float32)
        lesion_r = np.zeros(hemi_vertices, dtype=np.float32)
        paths = []
        if "L" in hemispheres:
            lesion_l = read_lesion_gifti(hemispheres["L"], hemi_vertices)
            paths.append(str(hemispheres["L"]))
        if "R" in hemispheres:
            lesion_r = read_lesion_gifti(hemispheres["R"], hemi_vertices)
            paths.append(str(hemispheres["R"]))
        lesion_records.append(
            {
                "subject": subject,
                "label_subject": label_subject,
                "session": session,
                "lesion_type": "TLE",
                "lesion_id": f"{subject}_{session}_TLE",
                "path": ";".join(paths),
                "lesion": np.concatenate([lesion_l, lesion_r]),
            }
        )

    return lesion_records


def read_gifti_vector(path, hemi_vertices):
    gifti = nib.load(str(path))
    if not gifti.darrays:
        raise ValueError(f"{path} has no GIFTI data arrays")

    data = np.asarray(gifti.darrays[0].data, dtype=np.float32).reshape(-1)
    if data.size != hemi_vertices:
        raise ValueError(f"{path} has {data.size} vertices; expected {hemi_vertices}")
    return data


def discover_wscore_maps(cortex_dir):
    maps = {}
    for path in cortex_dir.glob("*.func.gii"):
        match = WSCORE_PATTERN.search(path.name)
        if not match:
            continue
        metadata = match.groupdict()
        key = (
            metadata["feature"],
            metadata["label"],
            metadata["smooth"],
            metadata["analysis"],
        )
        maps.setdefault(key, {})[metadata["hemi"]] = path
    return maps


def read_gifti_vector_any(path):
    """Read a 1-D GIFTI metric vector of arbitrary length (hippocampal meshes
    vary by density, so we do not enforce a fixed vertex count here)."""
    gifti = nib.load(str(path))
    if not gifti.darrays:
        raise ValueError(f"{path} has no GIFTI data arrays")
    return np.asarray(gifti.darrays[0].data, dtype=np.float32).reshape(-1)


def discover_hipp_wscore_maps(hipp_dir, density=None):
    """Map (feature, analysis) -> {hemi: path} for hippocampal midthickness maps."""
    maps = {}
    for path in hipp_dir.glob("*.func.gii"):
        match = HIPP_WSCORE_PATTERN.search(path.name)
        if not match:
            continue
        metadata = match.groupdict()
        if density is not None and metadata["den"] != density:
            continue
        key = (metadata["feature"], metadata["analysis"])
        maps.setdefault(key, {})[metadata["hemi"]] = path
    return maps


def discover_hipp_cavities(lesion_dir, density=None):
    """Map (subject, session, hemi) -> hippocampal cavity GIFTI path."""
    cavities = {}
    for path in sorted(Path(lesion_dir).glob("**/*_label-hipp_cavity.func.gii")):
        match = HIPP_CAVITY_PATTERN.search(path.name)
        if not match:
            continue
        metadata = match.groupdict()
        if density is not None and metadata["den"] != density:
            continue
        cavities[(metadata["subject"], metadata["session"], metadata["hemi"])] = path
    return cavities


def load_hippocampal_appendix(
    subject,
    session,
    wscore_root,
    hipp_cavities,
    density,
):
    """Build the per-subject hippocampal data appended to the cortical vectors.

    Returns a dict with:
      ``n``           -- per-hemisphere hippocampal vertex count (None if unavailable)
      ``lesion``      -- bilateral [L, R] cavity ground truth (zeros where no cavity)
      ``data_by_key`` -- {(feature, analysis): bilateral [L, R] W-score vector}

    Hippocampal maps are matched to cortical maps by (feature, analysis). When a
    cortical feature has no hippocampal counterpart (e.g. *blur), callers append
    NaNs of length 2 * ``n`` so those vertices stay invalid for that feature.
    """
    empty = {"n": None, "lesion": None, "data_by_key": {}}
    hipp_dir = Path(wscore_root) / subject / session / "wscore_maps" / "hippocampus"
    wscore_maps = discover_hipp_wscore_maps(hipp_dir, density=density) if hipp_dir.exists() else {}

    data_by_key = {}
    hipp_n = None
    for (feature, analysis), hemispheres in wscore_maps.items():
        if "L" not in hemispheres or "R" not in hemispheres:
            continue
        left = read_gifti_vector_any(hemispheres["L"])
        right = read_gifti_vector_any(hemispheres["R"])
        if left.size != right.size:
            continue
        hipp_n = left.size
        data_by_key[(feature, analysis)] = np.concatenate([left, right])

    cavity_l_path = hipp_cavities.get((subject, session, "L"))
    cavity_r_path = hipp_cavities.get((subject, session, "R"))
    cavity_left = read_gifti_vector_any(cavity_l_path) if cavity_l_path else None
    cavity_right = read_gifti_vector_any(cavity_r_path) if cavity_r_path else None
    for cavity in (cavity_left, cavity_right):
        if cavity is not None:
            hipp_n = cavity.size if hipp_n is None else hipp_n

    if hipp_n is None:
        return empty

    if cavity_left is None:
        cavity_left = np.zeros(hipp_n, dtype=np.float32)
    if cavity_right is None:
        cavity_right = np.zeros(hipp_n, dtype=np.float32)
    if cavity_left.size != hipp_n or cavity_right.size != hipp_n:
        # Cavity density does not match the W-score mesh; skip hippocampus.
        return empty

    return {
        "n": hipp_n,
        "lesion": np.concatenate([cavity_left, cavity_right]),
        "data_by_key": data_by_key,
    }


def hippocampal_map_data(appendix, feature, analysis):
    """Return the bilateral hippocampal vector for a feature, or NaNs when absent."""
    if appendix["n"] is None:
        return None
    data = appendix["data_by_key"].get((feature, analysis))
    if data is not None:
        return data
    return np.full(2 * appendix["n"], np.nan, dtype=np.float32)


def load_sphere_coordinates(path, hemi_vertices):
    img = nib.load(str(path))
    if not img.darrays:
        raise ValueError(f"{path} has no GIFTI data arrays")
    coords = np.asarray(img.darrays[0].data, dtype=np.float64)
    if coords.shape[0] != hemi_vertices or coords.shape[1] != 3:
        raise ValueError(f"{path} has shape {coords.shape}; expected ({hemi_vertices}, 3)")
    norms = np.linalg.norm(coords, axis=1)
    if np.any(norms == 0):
        raise ValueError(f"{path} contains zero-length sphere coordinates")
    return coords / norms[:, None]


def random_rotation_matrix(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def spin_indices(coords, tree, rotation):
    rotated = coords @ rotation.T
    return tree.query(rotated, k=1)[1]


def summarize_values(values):
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "mean": np.nan,
            "median": np.nan,
            "std": np.nan,
            "q25": np.nan,
            "q75": np.nan,
        }

    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
        "q25": float(np.percentile(values, 25)),
        "q75": float(np.percentile(values, 75)),
    }


def detection_mask_from_values(data, threshold, tail):
    data = np.asarray(data, dtype=float)
    if tail == "positive":
        return data > threshold
    if tail in {"both", "absolute"}:
        return np.abs(data) > threshold
    if tail == "negative":
        return data < -threshold
    raise ValueError(f"Unsupported detection tail: {tail}")


def detection_scores_from_values(data, tail):
    data = np.asarray(data, dtype=float)
    if tail == "positive":
        return data
    if tail in {"both", "absolute"}:
        return np.abs(data)
    if tail == "negative":
        return -data
    raise ValueError(f"Unsupported detection tail: {tail}")


def summarize_data_with_mask(
    data,
    lesion,
    lesion_threshold,
    detection_threshold=1.96,
    detection_tail="both",
):
    lesion_valid = np.isfinite(lesion)
    lesion_mask = lesion > lesion_threshold
    valid = lesion_valid & np.isfinite(data)
    inside_mask = valid & lesion_mask
    outside_mask = valid & (~lesion_mask)
    detected_mask = valid & detection_mask_from_values(data, detection_threshold, detection_tail)
    true_positive = detected_mask & inside_mask
    false_positive = detected_mask & outside_mask
    false_negative = inside_mask & (~detected_mask)
    true_negative = outside_mask & (~detected_mask)
    inside_stats = summarize_values(data[inside_mask])
    outside_stats = summarize_values(data[outside_mask])

    tp = int(true_positive.sum())
    fp = int(false_positive.sum())
    fn = int(false_negative.sum())
    tn = int(true_negative.sum())
    n_lesion = int(inside_mask.sum())
    n_outside = int(outside_mask.sum())
    n_detected = int(detected_mask.sum())

    return {
        "n_valid_vertices": int(valid.sum()),
        "n_lesion_vertices": n_lesion,
        "n_outside_vertices": n_outside,
        "lesion_fraction_valid": float(inside_mask.sum() / valid.sum()) if valid.sum() else np.nan,
        "detection_threshold": float(detection_threshold),
        "detection_tail": detection_tail,
        "n_detected_vertices": n_detected,
        "n_true_positive_vertices": tp,
        "n_false_positive_vertices": fp,
        "n_false_negative_vertices": fn,
        "n_true_negative_vertices": tn,
        "sensitivity": float(tp / n_lesion) if n_lesion else np.nan,
        "false_positive_rate": float(fp / n_outside) if n_outside else np.nan,
        "precision": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        "dice": float((2 * tp) / ((2 * tp) + fp + fn)) if ((2 * tp) + fp + fn) else np.nan,
        "false_positive_to_lesion_ratio": float(fp / n_lesion) if n_lesion else np.nan,
        "inside_mean": inside_stats["mean"],
        "outside_mean": outside_stats["mean"],
        "mean_difference_inside_minus_outside": (
            inside_stats["mean"] - outside_stats["mean"]
            if np.isfinite(inside_stats["mean"]) and np.isfinite(outside_stats["mean"])
            else np.nan
        ),
        "inside_median": inside_stats["median"],
        "outside_median": outside_stats["median"],
        "inside_std": inside_stats["std"],
        "outside_std": outside_stats["std"],
        "inside_q25": inside_stats["q25"],
        "inside_q75": inside_stats["q75"],
        "outside_q25": outside_stats["q25"],
        "outside_q75": outside_stats["q75"],
    }


def load_vertex_records(args, wscore_root=None):
    records = []
    skipped = []
    wscore_root = Path(wscore_root) if wscore_root is not None else args.wscore_root
    lesion_records = discover_lesions(args.lesion_dir, args.hemi_vertices)
    include_hipp = getattr(args, "include_hippocampus", False)
    hipp_density = getattr(args, "hipp_density", DEFAULT_HIPP_DENSITY)
    hipp_cavities = discover_hipp_cavities(args.lesion_dir) if include_hipp else {}

    for lesion_record in lesion_records:
        subject = lesion_record["subject"]
        session = lesion_record.get("session", args.session)
        lesion = lesion_record["lesion"]
        cortex_dir = wscore_root / subject / session / "wscore_maps" / "cortex"
        if not cortex_dir.exists():
            skipped.append((subject, "missing cortex dir", str(cortex_dir)))
            continue

        appendix = {"n": None, "lesion": None, "data_by_key": {}}
        if include_hipp:
            appendix = load_hippocampal_appendix(
                subject, session, wscore_root, hipp_cavities, hipp_density
            )
            if appendix["n"] is not None:
                lesion = np.concatenate([lesion, appendix["lesion"]])

        maps = []
        wscore_maps = discover_wscore_maps(cortex_dir)
        for (feature, label, smooth, analysis), hemispheres in sorted(wscore_maps.items()):
            is_requested_surface = label == args.label and analysis == args.analysis
            is_requested_blur = (
                args.include_blur
                and feature.endswith("*blur")
                and label == "midthickness"
                and analysis == args.analysis
            )
            if not (is_requested_surface or is_requested_blur):
                continue
            if "L" not in hemispheres or "R" not in hemispheres:
                skipped.append((subject, f"missing hemisphere for {feature}/{label}/{analysis}", str(cortex_dir)))
                continue
            data = np.concatenate(
                [
                    read_gifti_vector(hemispheres["L"], args.hemi_vertices),
                    read_gifti_vector(hemispheres["R"], args.hemi_vertices),
                ]
            )
            if feature in SIGN_FLIPPED_FEATURES:
                data = -data
            hipp_data = hippocampal_map_data(appendix, feature, analysis)
            if hipp_data is not None:
                data = np.concatenate([data, hipp_data])
            maps.append(
                {
                    "feature": feature,
                    "label": label,
                    "smooth": smooth,
                    "analysis": analysis,
                    "data": data,
                }
            )
        if maps:
            records.append(
                {
                    "subject": subject,
                    "session": session,
                    "lesion_type": lesion_record["lesion_type"],
                    "lesion_id": lesion_record["lesion_id"],
                    "lesion_path": lesion_record["path"],
                    "lesion": lesion,
                    "maps": maps,
                }
            )
    return records, pd.DataFrame(skipped, columns=["subject", "reason", "path"])


def summary_from_vertex_records(
    records,
    lesion_threshold,
    detection_threshold=1.96,
    detection_tail="both",
):
    rows = []
    for record in records:
        for map_record in record["maps"]:
            row = {
                "subject": record["subject"],
                "session": record["session"],
                "lesion_type": record.get("lesion_type", "lesion"),
                "lesion_id": record.get("lesion_id", record["subject"]),
                "lesion_path": record.get("lesion_path", ""),
                "feature": map_record["feature"],
                "label": map_record["label"],
                "smooth": map_record["smooth"],
                "analysis": map_record["analysis"],
            }
            row.update(
                summarize_data_with_mask(
                    map_record["data"],
                    record["lesion"],
                    lesion_threshold,
                    detection_threshold=detection_threshold,
                    detection_tail=detection_tail,
                )
            )
            rows.append(row)
    return pd.DataFrame(rows)


def build_summary(args, wscore_root=None):
    rows = []
    skipped = []
    wscore_root = Path(wscore_root) if wscore_root is not None else args.wscore_root
    lesion_records = discover_lesions(args.lesion_dir, args.hemi_vertices)
    include_hipp = getattr(args, "include_hippocampus", False)
    hipp_density = getattr(args, "hipp_density", DEFAULT_HIPP_DENSITY)
    hipp_cavities = discover_hipp_cavities(args.lesion_dir) if include_hipp else {}

    for lesion_record in lesion_records:
        subject = lesion_record["subject"]
        session = lesion_record.get("session", args.session)
        lesion = lesion_record["lesion"]

        cortex_dir = (
            wscore_root
            / subject
            / session
            / "wscore_maps"
            / "cortex"
        )
        if not cortex_dir.exists():
            skipped.append((subject, "missing cortex dir", str(cortex_dir)))
            continue

        appendix = {"n": None, "lesion": None, "data_by_key": {}}
        if include_hipp:
            appendix = load_hippocampal_appendix(
                subject, session, wscore_root, hipp_cavities, hipp_density
            )
            if appendix["n"] is not None:
                lesion = np.concatenate([lesion, appendix["lesion"]])

        lesion_valid = np.isfinite(lesion)
        lesion_mask = lesion > args.lesion_threshold

        wscore_maps = discover_wscore_maps(cortex_dir)
        for (feature, label, smooth, analysis), hemispheres in sorted(wscore_maps.items()):
            is_requested_surface = label == args.label and analysis == args.analysis
            is_requested_blur = (
                args.include_blur
                and feature.endswith("*blur")
                and label == "midthickness"
                and analysis == args.analysis
            )
            if not (is_requested_surface or is_requested_blur):
                continue

            if "L" not in hemispheres or "R" not in hemispheres:
                skipped.append(
                    (
                        subject,
                        f"missing hemisphere for {feature}/{label}/{analysis}",
                        str(cortex_dir),
                    )
                )
                continue

            data = np.concatenate(
                [
                    read_gifti_vector(hemispheres["L"], args.hemi_vertices),
                    read_gifti_vector(hemispheres["R"], args.hemi_vertices),
                ]
            )
            if feature in SIGN_FLIPPED_FEATURES:
                data = -data
            hipp_data = hippocampal_map_data(appendix, feature, analysis)
            if hipp_data is not None:
                data = np.concatenate([data, hipp_data])
            valid = lesion_valid & np.isfinite(data)
            inside_mask = valid & lesion_mask
            outside_mask = valid & (~lesion_mask)
            inside = data[inside_mask]
            outside = data[outside_mask]
            inside_stats = summarize_values(inside)
            outside_stats = summarize_values(outside)

            rows.append(
                {
                    "subject": subject,
                    "session": session,
                    "lesion_type": lesion_record["lesion_type"],
                    "lesion_id": lesion_record["lesion_id"],
                    "lesion_path": lesion_record["path"],
                    "feature": feature,
                    "label": label,
                    "smooth": smooth,
                    "analysis": analysis,
                    "n_valid_vertices": int(valid.sum()),
                    "n_lesion_vertices": int(inside_mask.sum()),
                    "n_outside_vertices": int(outside_mask.sum()),
                    "lesion_fraction_valid": (
                        float(inside_mask.sum() / valid.sum()) if valid.sum() else np.nan
                    ),
                    "inside_mean": inside_stats["mean"],
                    "outside_mean": outside_stats["mean"],
                    "mean_difference_inside_minus_outside": (
                        inside_stats["mean"] - outside_stats["mean"]
                        if np.isfinite(inside_stats["mean"])
                        and np.isfinite(outside_stats["mean"])
                        else np.nan
                    ),
                    "inside_median": inside_stats["median"],
                    "outside_median": outside_stats["median"],
                    "inside_std": inside_stats["std"],
                    "outside_std": outside_stats["std"],
                    "inside_q25": inside_stats["q25"],
                    "inside_q75": inside_stats["q75"],
                    "outside_q25": outside_stats["q25"],
                    "outside_q75": outside_stats["q75"],
                }
            )

    return pd.DataFrame(rows), pd.DataFrame(skipped, columns=["subject", "reason", "path"])


def benjamini_hochberg(p_values):
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(p_values.shape, np.nan, dtype=float)
    valid = np.isfinite(p_values)
    if not np.any(valid):
        return adjusted

    valid_p = p_values[valid]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    n = ranked.size
    adjusted_ranked = ranked * n / np.arange(1, n + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0, 1)

    valid_indices = np.where(valid)[0]
    adjusted[valid_indices[order]] = adjusted_ranked
    return adjusted


def p_to_stars(p_value):
    if not np.isfinite(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def metric_order_key(metric, group_blur_last=False):
    feature, label, analysis = metric
    try:
        feature_idx = FEATURE_ORDER.index(feature)
    except ValueError:
        feature_idx = len(FEATURE_ORDER)
    blur_group = 1 if group_blur_last and feature.endswith("*blur") else 0
    return blur_group, feature_idx, feature, label, analysis


def ordered_metrics(subset, group_blur_last=False):
    metrics = subset[["feature", "label", "analysis"]].drop_duplicates()
    values = [tuple(row) for row in metrics.to_numpy()]
    return sorted(values, key=lambda metric: metric_order_key(metric, group_blur_last))


def format_metric_label(feature, label, show_surface=False):
    feature_label = FEATURE_ABBREVIATIONS.get(feature, feature)
    if show_surface and feature.endswith("*blur"):
        return feature_label
    if show_surface:
        return feature_label
    return feature_label


def get_test_row(tests, feature, label, analysis):
    if tests is None or tests.empty:
        return None
    match = tests[
        (tests["feature"] == feature)
        & (tests["label"] == label)
        & (tests["analysis"] == analysis)
    ]
    if match.empty:
        return None
    return match.iloc[0]


def compute_inside_outside_tests(summary, alpha=0.05):
    """Run paired subject-level tests for each feature/label/analysis."""
    from scipy import stats

    rows = []
    if summary.empty:
        return pd.DataFrame(rows)

    group_cols = ["label", "analysis", "feature"]
    for (label, analysis, feature), group in summary.groupby(group_cols, sort=False):
        paired = group[
            ["subject", "inside_mean", "outside_mean", "mean_difference_inside_minus_outside"]
        ].dropna()
        paired = paired[np.isfinite(paired["inside_mean"]) & np.isfinite(paired["outside_mean"])]
        n = int(len(paired))

        inside = paired["inside_mean"].to_numpy(dtype=float)
        outside = paired["outside_mean"].to_numpy(dtype=float)
        diff = inside - outside

        mean_inside = float(np.mean(inside)) if n else np.nan
        mean_outside = float(np.mean(outside)) if n else np.nan
        mean_diff = float(np.mean(diff)) if n else np.nan
        median_diff = float(np.median(diff)) if n else np.nan
        sd_diff = float(np.std(diff, ddof=1)) if n > 1 else np.nan
        cohen_dz = mean_diff / sd_diff if np.isfinite(sd_diff) and sd_diff > 0 else np.nan
        pct_inside_greater = float(np.mean(diff > 0)) if n else np.nan

        t_stat = np.nan
        t_p = np.nan
        wilcoxon_stat = np.nan
        wilcoxon_p = np.nan
        if n >= 2:
            t_result = stats.ttest_rel(inside, outside, nan_policy="omit")
            t_stat = float(t_result.statistic)
            t_p = float(t_result.pvalue)
            if not np.allclose(diff, 0):
                try:
                    w_result = stats.wilcoxon(diff, zero_method="wilcox", alternative="two-sided")
                    wilcoxon_stat = float(w_result.statistic)
                    wilcoxon_p = float(w_result.pvalue)
                except ValueError:
                    pass

        rows.append(
            {
                "label": label,
                "analysis": analysis,
                "feature": feature,
                "n_subjects": n,
                "mean_inside": mean_inside,
                "mean_outside": mean_outside,
                "mean_difference_inside_minus_outside": mean_diff,
                "median_difference_inside_minus_outside": median_diff,
                "sd_difference": sd_diff,
                "cohen_dz": cohen_dz,
                "percent_subjects_inside_greater": pct_inside_greater,
                "paired_t_stat": t_stat,
                "paired_t_p": t_p,
                "wilcoxon_stat": wilcoxon_stat,
                "wilcoxon_p": wilcoxon_p,
            }
        )

    tests = pd.DataFrame(rows)
    if tests.empty:
        return tests

    tests["paired_t_p_fdr"] = np.nan
    tests["wilcoxon_p_fdr"] = np.nan
    for (label, analysis), idx in tests.groupby(["label", "analysis"]).groups.items():
        idx = list(idx)
        tests.loc[idx, "paired_t_p_fdr"] = benjamini_hochberg(tests.loc[idx, "paired_t_p"].to_numpy())
        tests.loc[idx, "wilcoxon_p_fdr"] = benjamini_hochberg(tests.loc[idx, "wilcoxon_p"].to_numpy())

    tests["significant_wilcoxon_fdr"] = tests["wilcoxon_p_fdr"] < alpha
    tests["significant_paired_t_fdr"] = tests["paired_t_p_fdr"] < alpha
    tests["direction"] = np.where(
        tests["mean_difference_inside_minus_outside"] > 0,
        "inside_higher",
        np.where(tests["mean_difference_inside_minus_outside"] < 0, "inside_lower", "no_difference"),
    )
    tests = tests.sort_values(
        ["label", "analysis", "wilcoxon_p_fdr", "paired_t_p_fdr", "feature"],
        na_position="last",
    )
    return tests


def compute_detection_performance(summary):
    """Aggregate thresholded lesion-detection and false-positive counts."""
    required = {
        "n_true_positive_vertices",
        "n_false_positive_vertices",
        "n_false_negative_vertices",
        "n_true_negative_vertices",
    }
    if summary.empty or not required.issubset(summary.columns):
        return pd.DataFrame()

    rows = []
    group_cols = ["label", "analysis", "feature"]
    for (label, analysis, feature), group in summary.groupby(group_cols, sort=False):
        tp = int(group["n_true_positive_vertices"].sum())
        fp = int(group["n_false_positive_vertices"].sum())
        fn = int(group["n_false_negative_vertices"].sum())
        tn = int(group["n_true_negative_vertices"].sum())
        lesion_vertices = tp + fn
        outside_vertices = fp + tn
        detected_vertices = tp + fp

        rows.append(
            {
                "label": label,
                "analysis": analysis,
                "feature": feature,
                "n_subjects": int(group["subject"].nunique()),
                "detection_threshold": float(group["detection_threshold"].iloc[0])
                if "detection_threshold" in group
                else np.nan,
                "detection_tail": str(group["detection_tail"].iloc[0])
                if "detection_tail" in group
                else "",
                "total_lesion_vertices": int(lesion_vertices),
                "total_outside_vertices": int(outside_vertices),
                "total_detected_vertices": int(detected_vertices),
                "total_true_positive_vertices": tp,
                "total_false_positive_vertices": fp,
                "total_false_negative_vertices": fn,
                "total_true_negative_vertices": tn,
                "sensitivity": float(tp / lesion_vertices) if lesion_vertices else np.nan,
                "false_positive_rate": float(fp / outside_vertices) if outside_vertices else np.nan,
                "precision": float(tp / detected_vertices) if detected_vertices else np.nan,
                "dice": float((2 * tp) / ((2 * tp) + fp + fn)) if ((2 * tp) + fp + fn) else np.nan,
                "false_positive_to_lesion_ratio": float(fp / lesion_vertices) if lesion_vertices else np.nan,
                "mean_false_positive_vertices_per_subject": float(
                    group.groupby("subject")["n_false_positive_vertices"].sum().mean()
                ),
            }
        )

    performance = pd.DataFrame(rows)
    if performance.empty:
        return performance
    return performance.sort_values(
        ["false_positive_rate", "dice", "feature"],
        ascending=[True, False, True],
        na_position="last",
    )


def auc_from_scores(labels, scores):
    """Mann-Whitney AUC with average ranks for tied scores."""
    from scipy.stats import rankdata

    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    valid = np.isfinite(scores)
    labels = labels[valid]
    scores = scores[valid]
    n_pos = int(labels.sum())
    n_neg = int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan

    ranks = rankdata(scores, method="average")
    pos_rank_sum = float(ranks[labels].sum())
    return float((pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def roc_thresholds_from_args(args):
    if args.roc_threshold_steps <= 1:
        return np.asarray([float(args.roc_threshold_max)], dtype=float)
    return np.linspace(
        float(args.roc_threshold_min),
        float(args.roc_threshold_max),
        int(args.roc_threshold_steps),
        dtype=float,
    )


def compute_roc_curves(records, lesion_threshold, detection_tail="both", thresholds=None, include_pooled=True):
    """Compute vertex-level ROC curves using the lesion mask as ground truth."""
    if thresholds is None:
        thresholds = np.linspace(0.0, 6.0, 121, dtype=float)
    thresholds = np.asarray(thresholds, dtype=float)

    metric_vectors = {}
    for record in records:
        lesion = np.asarray(record["lesion"], dtype=float)
        lesion_valid = np.isfinite(lesion)
        lesion_label = lesion > lesion_threshold
        for map_record in record["maps"]:
            data = np.asarray(map_record["data"], dtype=float)
            valid = lesion_valid & np.isfinite(data)
            if not np.any(valid):
                continue
            key = (map_record["feature"], map_record["label"], map_record["smooth"], map_record["analysis"])
            metric_vectors.setdefault(key, {"labels": [], "scores": []})
            metric_vectors[key]["labels"].append(lesion_label[valid])
            metric_vectors[key]["scores"].append(detection_scores_from_values(data[valid], detection_tail))

    if include_pooled and metric_vectors:
        pooled = {"labels": [], "scores": []}
        for vectors in metric_vectors.values():
            pooled["labels"].extend(vectors["labels"])
            pooled["scores"].extend(vectors["scores"])
        metric_vectors[("ALL_METRICS", "all", "all", "all")] = pooled

    roc_rows = []
    auc_rows = []
    for (feature, label, smooth, analysis), vectors in metric_vectors.items():
        labels = np.concatenate(vectors["labels"]).astype(bool)
        scores = np.concatenate(vectors["scores"]).astype(float)
        valid = np.isfinite(scores)
        labels = labels[valid]
        scores = scores[valid]
        n_pos = int(labels.sum())
        n_neg = int((~labels).sum())
        if n_pos == 0 or n_neg == 0:
            continue

        best_youden = -np.inf
        best_threshold = np.nan
        best_sensitivity = np.nan
        best_fpr = np.nan
        for threshold in thresholds:
            detected = scores > threshold
            tp = int(np.sum(detected & labels))
            fp = int(np.sum(detected & (~labels)))
            fn = n_pos - tp
            tn = n_neg - fp
            sensitivity = float(tp / n_pos)
            fpr = float(fp / n_neg)
            specificity = 1.0 - fpr
            youden = sensitivity - fpr
            precision = float(tp / (tp + fp)) if (tp + fp) else np.nan
            dice = float((2 * tp) / ((2 * tp) + fp + fn)) if ((2 * tp) + fp + fn) else np.nan
            if np.isfinite(youden) and youden > best_youden:
                best_youden = youden
                best_threshold = float(threshold)
                best_sensitivity = sensitivity
                best_fpr = fpr

            roc_rows.append(
                {
                    "feature": feature,
                    "label": label,
                    "smooth": smooth,
                    "analysis": analysis,
                    "detection_tail": detection_tail,
                    "lesion_threshold": float(lesion_threshold),
                    "detection_threshold": float(threshold),
                    "n_lesion_vertices": n_pos,
                    "n_outside_vertices": n_neg,
                    "true_positive_vertices": tp,
                    "false_positive_vertices": fp,
                    "false_negative_vertices": fn,
                    "true_negative_vertices": tn,
                    "sensitivity": sensitivity,
                    "false_positive_rate": fpr,
                    "specificity": specificity,
                    "precision": precision,
                    "dice": dice,
                    "youden_j": youden,
                }
            )

        auc_rows.append(
            {
                "feature": feature,
                "label": label,
                "smooth": smooth,
                "analysis": analysis,
                "detection_tail": detection_tail,
                "lesion_threshold": float(lesion_threshold),
                "n_lesion_vertices": n_pos,
                "n_outside_vertices": n_neg,
                "auc": auc_from_scores(labels, scores),
                "best_youden_j": float(best_youden) if np.isfinite(best_youden) else np.nan,
                "best_detection_threshold": best_threshold,
                "best_sensitivity": best_sensitivity,
                "best_false_positive_rate": best_fpr,
            }
        )

    roc = pd.DataFrame(roc_rows)
    auc = pd.DataFrame(auc_rows)
    if not auc.empty:
        auc = auc.sort_values(["auc", "feature"], ascending=[False, True], na_position="last")
    return roc, auc


def plot_roc_curves(
    roc,
    auc,
    output_path,
    title,
    max_curves=12,
    current_threshold=None,
    pooled_only=False,
):
    if roc.empty or auc.empty:
        return None

    import matplotlib.pyplot as plt

    max_curves = max(1, int(max_curves))
    if pooled_only:
        top_auc = auc[auc["feature"] == "ALL_METRICS"].copy()
    else:
        pooled_auc = auc[auc["feature"] == "ALL_METRICS"].copy()
        regular_auc = auc[auc["feature"] != "ALL_METRICS"].head(max(0, max_curves - len(pooled_auc))).copy()
        top_auc = pd.concat([pooled_auc, regular_auc], ignore_index=True)
    if top_auc.empty:
        if pooled_only:
            return None
        top_auc = auc.head(max_curves).copy()
    fig, ax = plt.subplots(figsize=(7.2, 6.4), facecolor="white")
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(top_auc), 1)))

    for color, (_, row) in zip(colors, top_auc.iterrows()):
        subset = roc[
            (roc["feature"] == row["feature"])
            & (roc["label"] == row["label"])
            & (roc["smooth"] == row["smooth"])
            & (roc["analysis"] == row["analysis"])
        ].copy()
        if "method" in row and "method" in roc.columns:
            subset = subset[subset["method"] == row["method"]]
        if subset.empty:
            continue

        subset = subset.sort_values("false_positive_rate")
        label_prefix = f"{row['method']} " if "method" in row else ""
        is_pooled = row["feature"] == "ALL_METRICS"
        if is_pooled:
            curve_label = f"{label_prefix}All features AUC={row['auc']:.3f}"
        else:
            curve_label = (
                f"{label_prefix}{FEATURE_ABBREVIATIONS.get(row['feature'], row['feature'])}"
                f" [{row['label']}] AUC={row['auc']:.3f}"
            )
        ax.plot(
            subset["false_positive_rate"],
            subset["sensitivity"],
            lw=3.0 if is_pooled else 2.0,
            color=color,
            label=curve_label,
            alpha=1.0 if is_pooled else 0.88,
            linestyle="-" if is_pooled else "--",
        )
        if current_threshold is not None and "detection_threshold" in subset.columns:
            threshold_idx = (subset["detection_threshold"] - current_threshold).abs().idxmin()
            threshold_row = subset.loc[threshold_idx]
            ax.scatter(
                threshold_row["false_positive_rate"],
                threshold_row["sensitivity"],
                s=36,
                color=color,
                edgecolor="white",
                linewidth=0.8,
                zorder=4,
            )

    ax.plot([0, 1], [0, 1], color="#a9a9a9", lw=1.2, linestyle="--", label="Chance")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False-positive rate outside lesion", fontsize=12.5)
    ax.set_ylabel("Sensitivity inside lesion", fontsize=12.5)
    ax.set_title(title, fontsize=15, pad=16)
    ax.grid(color="#e7e9eb", lw=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=8.8, loc="lower right")
    fig.text(
        0.01,
        0.01,
        "ROC points sweep the W-score detection threshold; dots mark the configured detection threshold.",
        ha="left",
        va="bottom",
        fontsize=9,
        color="#6a6a6a",
    )
    fig.tight_layout(rect=[0, 0.035, 1, 1])
    fig.savefig(output_path, dpi=250)
    plt.close(fig)
    return output_path


def plot_roc_feature_subplots(
    roc,
    auc,
    output_path,
    title,
    current_threshold=None,
):
    if roc.empty or auc.empty:
        return None

    import matplotlib.pyplot as plt

    feature_auc = auc[auc["feature"] != "ALL_METRICS"].copy()
    if feature_auc.empty:
        return plot_roc_curves(
            roc,
            auc,
            output_path,
            title,
            current_threshold=current_threshold,
            pooled_only=True,
        )

    features = ordered_features(feature_auc["feature"].dropna().unique())
    n_features = len(features)
    n_cols = min(3, max(1, n_features))
    n_rows = int(np.ceil(n_features / n_cols))
    fig_width = max(6.2, n_cols * 4.8)
    fig_height = max(4.8, n_rows * 4.0)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(fig_width, fig_height),
        facecolor="white",
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    axes = axes.reshape(-1)

    if "method" in auc.columns:
        methods = list(dict.fromkeys(auc["method"].dropna().astype(str)))
    else:
        methods = []
    color_values = plt.cm.tab10(np.linspace(0, 1, max(len(methods), 1)))
    method_colors = {method: color_values[i % len(color_values)] for i, method in enumerate(methods)}
    fallback_colors = plt.cm.tab20(np.linspace(0, 1, 20))

    for feature_idx, feature in enumerate(features):
        ax = axes[feature_idx]
        feature_rows = feature_auc[feature_auc["feature"] == feature].copy()
        if feature_rows.empty:
            ax.set_visible(False)
            continue

        if "method" in feature_rows.columns:
            feature_rows = (
                feature_rows.sort_values(["auc", "method"], ascending=[False, True], na_position="last")
                .drop_duplicates(["method"])
                .sort_values("method")
            )
        else:
            feature_rows = feature_rows.sort_values("auc", ascending=False, na_position="last")

        for row_idx, (_, row) in enumerate(feature_rows.iterrows()):
            subset = roc[
                (roc["feature"] == row["feature"])
                & (roc["label"] == row["label"])
                & (roc["smooth"] == row["smooth"])
                & (roc["analysis"] == row["analysis"])
            ].copy()
            method = ""
            if "method" in row and pd.notna(row["method"]):
                method = str(row["method"])
            if method and "method" in subset.columns:
                subset = subset[subset["method"] == row["method"]]
            if subset.empty:
                continue

            subset = subset.sort_values("false_positive_rate")
            color = method_colors.get(method, fallback_colors[row_idx % len(fallback_colors)])
            curve_label = (
                f"{method} AUC={row['auc']:.3f}"
                if method
                else f"{row['label']} AUC={row['auc']:.3f}"
            )
            ax.plot(
                subset["false_positive_rate"],
                subset["sensitivity"],
                lw=2.3,
                color=color,
                label=curve_label,
                alpha=0.92,
            )
            if current_threshold is not None and "detection_threshold" in subset.columns:
                threshold_idx = (subset["detection_threshold"] - current_threshold).abs().idxmin()
                threshold_row = subset.loc[threshold_idx]
                ax.scatter(
                    threshold_row["false_positive_rate"],
                    threshold_row["sensitivity"],
                    s=30,
                    color=color,
                    edgecolor="white",
                    linewidth=0.75,
                    zorder=4,
                )

        feature_label = FEATURE_ABBREVIATIONS.get(feature, feature)
        ax.set_title(feature_label, fontsize=14, pad=10)
        ax.plot([0, 1], [0, 1], color="#b7b7b7", lw=1.0, linestyle="--", zorder=1)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(color="#e7e9eb", lw=0.8)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, fontsize=8.8, loc="lower right")

    for ax in axes[n_features:]:
        ax.set_visible(False)

    for ax in axes:
        if ax.get_visible():
            ax.tick_params(axis="both", labelsize=10.5)
    for ax in axes[(n_rows - 1) * n_cols:n_rows * n_cols]:
        if ax.get_visible():
            ax.set_xlabel("False-positive rate outside lesion", fontsize=11.5)
    for row_idx in range(n_rows):
        ax = axes[row_idx * n_cols]
        if ax.get_visible():
            ax.set_ylabel("Sensitivity inside lesion", fontsize=11.5)

    fig.suptitle(title, fontsize=17, y=0.985)
    fig.text(
        0.01,
        0.01,
        "Each panel shows one feature across processed datasets; dots mark the configured detection threshold.",
        ha="left",
        va="bottom",
        fontsize=9.5,
        color="#6a6a6a",
    )
    fig.tight_layout(rect=[0, 0.035, 1, 0.955])
    fig.savefig(output_path, dpi=250)
    plt.close(fig)
    return output_path


def aggregate_detection_counts(summary):
    required = {
        "n_true_positive_vertices",
        "n_false_positive_vertices",
        "n_false_negative_vertices",
        "n_true_negative_vertices",
    }
    if summary.empty or not required.issubset(summary.columns):
        return {
            "total_true_positive_vertices": 0,
            "total_false_positive_vertices": 0,
            "total_false_negative_vertices": 0,
            "total_true_negative_vertices": 0,
            "sensitivity": np.nan,
            "false_positive_rate": np.nan,
            "precision": np.nan,
            "dice": np.nan,
        }

    tp = int(summary["n_true_positive_vertices"].sum())
    fp = int(summary["n_false_positive_vertices"].sum())
    fn = int(summary["n_false_negative_vertices"].sum())
    tn = int(summary["n_true_negative_vertices"].sum())
    lesion_vertices = tp + fn
    outside_vertices = fp + tn
    detected_vertices = tp + fp
    return {
        "total_true_positive_vertices": tp,
        "total_false_positive_vertices": fp,
        "total_false_negative_vertices": fn,
        "total_true_negative_vertices": tn,
        "sensitivity": float(tp / lesion_vertices) if lesion_vertices else np.nan,
        "false_positive_rate": float(fp / outside_vertices) if outside_vertices else np.nan,
        "precision": float(tp / detected_vertices) if detected_vertices else np.nan,
        "dice": float((2 * tp) / ((2 * tp) + fp + fn)) if ((2 * tp) + fp + fn) else np.nan,
    }


def signed_feature_label(metric, sign):
    feature, label, analysis = metric
    prefix = "+" if sign > 0 else "-"
    return f"{prefix}{feature}[{label},{analysis}]"


def metric_label(metric):
    feature, label, analysis = metric
    return f"{feature}[{label},{analysis}]"


def parse_combination_operators(operators):
    valid = {"+", "-", "*", "/"}
    parsed = [op.strip() for op in str(operators).split(",") if op.strip()]
    invalid = [op for op in parsed if op not in valid]
    if invalid:
        raise ValueError(f"Unsupported combination operators: {invalid}. Use any of +,-,*,/")
    return parsed or ["+", "-"]


def safe_divide(numerator, denominator, eps=1e-6):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    return np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan, dtype=float),
        where=np.abs(denominator) > eps,
    )


def apply_operator(lhs, rhs, operator):
    if operator == "+":
        return lhs + rhs
    if operator == "-":
        return lhs - rhs
    if operator == "*":
        return lhs * rhs
    if operator == "/":
        return safe_divide(lhs, rhs)
    raise ValueError(f"Unsupported operator: {operator}")


def arithmetic_expression_label(combo, operators, oriented_sign=1):
    expression = metric_label(combo[0])
    for operator, metric in zip(operators, combo[1:]):
        expression = f"({expression} {operator} {metric_label(metric)})"
    if oriented_sign < 0:
        expression = f"-1 * {expression}"
    return expression


def metric_spec(metric):
    return "|".join(metric)


def parse_metric_spec(spec):
    parts = str(spec).split("|")
    if len(parts) != 3:
        raise ValueError(f"Invalid metric spec: {spec}")
    return tuple(parts)


def compute_combination_tests(
    summary,
    max_size=4,
    min_subjects=5,
    alpha=0.05,
    operators=None,
):
    """Test arithmetic feature combinations for stronger lesion detection."""
    from scipy import stats

    if summary.empty or max_size < 2:
        return pd.DataFrame()
    operators = parse_combination_operators(",".join(operators) if isinstance(operators, list) else operators or "+")

    usable = summary[
        np.isfinite(summary["inside_mean"])
        & np.isfinite(summary["outside_mean"])
    ].copy()
    if usable.empty:
        return pd.DataFrame()

    metrics = ordered_metrics(usable, group_blur_last=True)
    metric_tables = {}
    for metric in metrics:
        feature, label, analysis = metric
        metric_tables[metric] = usable[
            (usable["feature"] == feature)
            & (usable["label"] == label)
            & (usable["analysis"] == analysis)
        ].set_index("subject")[["inside_mean", "outside_mean"]]

    rows = []
    max_size = min(max_size, len(metrics))
    for size in range(2, max_size + 1):
        for combo in itertools.combinations(metrics, size):
            common_subjects = set(metric_tables[combo[0]].index)
            for metric in combo[1:]:
                common_subjects &= set(metric_tables[metric].index)
            subjects = sorted(common_subjects)
            if len(subjects) < min_subjects:
                continue

            for operator_pattern in itertools.product(operators, repeat=size - 1):
                first_table = metric_tables[combo[0]].loc[subjects]
                inside = first_table["inside_mean"].to_numpy(dtype=float)
                outside = first_table["outside_mean"].to_numpy(dtype=float)
                for operator, metric in zip(operator_pattern, combo[1:]):
                    table = metric_tables[metric].loc[subjects]
                    inside = apply_operator(
                        inside,
                        table["inside_mean"].to_numpy(dtype=float),
                        operator,
                    )
                    outside = apply_operator(
                        outside,
                        table["outside_mean"].to_numpy(dtype=float),
                        operator,
                    )

                finite = np.isfinite(inside) & np.isfinite(outside)
                inside = inside[finite]
                outside = outside[finite]
                diff = inside - outside
                n = int(len(diff))
                if n < min_subjects:
                    continue
                mean_diff = float(np.mean(diff))

                oriented_sign = 1
                if mean_diff < 0:
                    oriented_sign = -1
                    inside, outside = -inside, -outside
                    diff = -diff
                    mean_diff = -mean_diff

                sd_diff = float(np.std(diff, ddof=1)) if n > 1 else np.nan
                cohen_dz = mean_diff / sd_diff if np.isfinite(sd_diff) and sd_diff > 0 else np.nan
                t_stat = np.nan
                t_p = np.nan
                wilcoxon_stat = np.nan
                wilcoxon_p = np.nan
                if n >= 2:
                    t_result = stats.ttest_rel(inside, outside, nan_policy="omit")
                    t_stat = float(t_result.statistic)
                    t_p = float(t_result.pvalue)
                    if not np.allclose(diff, 0):
                        try:
                            w_result = stats.wilcoxon(diff, zero_method="wilcox", alternative="greater")
                            wilcoxon_stat = float(w_result.statistic)
                            wilcoxon_p = float(w_result.pvalue)
                        except ValueError:
                            pass

                rows.append(
                    {
                        "combination_size": size,
                        "operator_pattern": " ".join(operator_pattern),
                        "metrics": ";".join(metric_spec(metric) for metric in combo),
                        "orientation": oriented_sign,
                        "combination": arithmetic_expression_label(combo, operator_pattern, oriented_sign),
                        "n_subjects": n,
                        "mean_inside": float(np.mean(inside)),
                        "mean_outside": float(np.mean(outside)),
                        "mean_difference_inside_minus_outside": mean_diff,
                        "median_difference_inside_minus_outside": float(np.median(diff)),
                        "sd_difference": sd_diff,
                        "cohen_dz": cohen_dz,
                        "percent_subjects_inside_greater": float(np.mean(diff > 0)),
                        "paired_t_stat": t_stat,
                        "paired_t_p": t_p,
                        "wilcoxon_stat": wilcoxon_stat,
                        "wilcoxon_p": wilcoxon_p,
                    }
                )

    combos = pd.DataFrame(rows)
    if combos.empty:
        return combos

    combos["paired_t_p_fdr"] = benjamini_hochberg(combos["paired_t_p"].to_numpy())
    combos["wilcoxon_p_fdr"] = benjamini_hochberg(combos["wilcoxon_p"].to_numpy())
    combos["significant_wilcoxon_fdr"] = combos["wilcoxon_p_fdr"] < alpha
    combos["significant_paired_t_fdr"] = combos["paired_t_p_fdr"] < alpha
    combos = combos.sort_values(
        [
            "paired_t_p_fdr",
            "paired_t_p",
            "cohen_dz",
            "combination_size",
            "combination",
        ],
        ascending=[True, True, False, True, True],
        na_position="last",
    )
    return combos


def permute_inside_outside_by_subject(summary, rng):
    permuted = summary.copy()
    subjects = sorted(permuted["subject"].dropna().unique())
    flip_by_subject = {
        subject: rng.choice([-1.0, 1.0])
        for subject in subjects
    }

    inside = permuted["inside_mean"].to_numpy(dtype=float)
    outside = permuted["outside_mean"].to_numpy(dtype=float)
    midpoint = 0.5 * (inside + outside)
    half_diff = 0.5 * (inside - outside)
    flips = permuted["subject"].map(flip_by_subject).to_numpy(dtype=float)

    permuted["inside_mean"] = midpoint + flips * half_diff
    permuted["outside_mean"] = midpoint - flips * half_diff
    permuted["mean_difference_inside_minus_outside"] = (
        permuted["inside_mean"] - permuted["outside_mean"]
    )
    return permuted


def add_max_statistic_permutation_pvalues(
    observed_combinations,
    summary,
    max_size=4,
    min_subjects=5,
    alpha=0.05,
    operators=None,
    n_permutations=1000,
    seed=12345,
    verbose=True,
):
    if observed_combinations.empty or n_permutations <= 0:
        observed_combinations = observed_combinations.copy()
        observed_combinations["procedure_permutation_p"] = np.nan
        observed_combinations["significant_procedure_permutation"] = False
        return observed_combinations, pd.DataFrame()

    rng = np.random.default_rng(seed)
    null_rows = []
    null_max_t = np.full(n_permutations, np.nan, dtype=float)
    for i in range(n_permutations):
        permuted_summary = permute_inside_outside_by_subject(summary, rng)
        permuted_combinations = compute_combination_tests(
            permuted_summary,
            max_size=max_size,
            min_subjects=min_subjects,
            alpha=alpha,
            operators=operators,
        )
        if not permuted_combinations.empty:
            t_values = np.abs(permuted_combinations["paired_t_stat"].to_numpy(dtype=float))
            if np.any(np.isfinite(t_values)):
                null_max_t[i] = float(np.nanmax(t_values))
        null_rows.append({"permutation": i + 1, "max_abs_paired_t_stat": null_max_t[i]})
        if verbose and (i + 1) % max(1, n_permutations // 10) == 0:
            print(f"  permutation max-statistic progress: {i + 1}/{n_permutations}")

    observed = observed_combinations.copy()
    observed_t = np.abs(observed["paired_t_stat"].to_numpy(dtype=float))
    valid_null = null_max_t[np.isfinite(null_max_t)]
    p_values = np.full(observed_t.shape, np.nan, dtype=float)
    if valid_null.size:
        for idx, value in enumerate(observed_t):
            if np.isfinite(value):
                p_values[idx] = (1.0 + np.sum(valid_null >= value)) / (valid_null.size + 1.0)

    observed["procedure_permutation_p"] = p_values
    observed["significant_procedure_permutation"] = observed["procedure_permutation_p"] < alpha
    observed = observed.sort_values(
        [
            "procedure_permutation_p",
            "paired_t_p_fdr",
            "paired_t_p",
            "cohen_dz",
            "combination_size",
            "combination",
        ],
        ascending=[True, True, True, False, True, True],
        na_position="last",
    )
    return observed, pd.DataFrame(null_rows)


def spin_vertex_records(records, left_indices, right_indices, hemi_vertices):
    spun = []
    cortex_vertices = 2 * hemi_vertices
    for record in records:
        lesion = record["lesion"]
        # Only the cortical fsLR portion has a matching sphere; any appended
        # hippocampal vertices are carried through unspun.
        spun_lesion = np.concatenate(
            [
                lesion[:hemi_vertices][left_indices],
                lesion[hemi_vertices:cortex_vertices][right_indices],
                lesion[cortex_vertices:],
            ]
        )
        spun.append(
            {
                "subject": record["subject"],
                "session": record["session"],
                "lesion_type": record.get("lesion_type", "lesion"),
                "lesion_id": record.get("lesion_id", record["subject"]),
                "lesion_path": record.get("lesion_path", ""),
                "lesion": spun_lesion,
                "maps": record["maps"],
            }
        )
    return spun


def add_spin_max_statistic_pvalues(
    observed_combinations,
    records,
    args,
    n_spins=1000,
    seed=54321,
    verbose=True,
):
    if observed_combinations.empty or n_spins <= 0:
        observed = observed_combinations.copy()
        observed["spin_procedure_p"] = np.nan
        observed["significant_spin_procedure"] = False
        return observed, pd.DataFrame()

    from scipy.spatial import cKDTree

    coords_l = load_sphere_coordinates(args.sphere_l, args.hemi_vertices)
    coords_r = load_sphere_coordinates(args.sphere_r, args.hemi_vertices)
    tree_l = cKDTree(coords_l)
    tree_r = cKDTree(coords_r)
    rng = np.random.default_rng(seed)

    null_rows = []
    null_max_t = np.full(n_spins, np.nan, dtype=float)
    for i in range(n_spins):
        rotation = random_rotation_matrix(rng)
        left_indices = spin_indices(coords_l, tree_l, rotation)
        right_indices = spin_indices(coords_r, tree_r, rotation)
        spin_records = spin_vertex_records(records, left_indices, right_indices, args.hemi_vertices)
        spin_summary = summary_from_vertex_records(
            spin_records,
            args.lesion_threshold,
            detection_threshold=args.detection_threshold,
            detection_tail=args.detection_tail,
        )
        spin_combinations = compute_combination_tests(
            spin_summary,
            max_size=args.max_combination_size,
            min_subjects=args.min_combination_subjects,
            alpha=args.alpha,
            operators=args.combination_operators,
        )
        if not spin_combinations.empty:
            t_values = np.abs(spin_combinations["paired_t_stat"].to_numpy(dtype=float))
            if np.any(np.isfinite(t_values)):
                null_max_t[i] = float(np.nanmax(t_values))
        null_rows.append({"spin": i + 1, "max_abs_paired_t_stat": null_max_t[i]})
        if verbose and (i + 1) % max(1, n_spins // 10) == 0:
            print(f"  spin max-statistic progress: {i + 1}/{n_spins}")

    observed = observed_combinations.copy()
    observed_t = np.abs(observed["paired_t_stat"].to_numpy(dtype=float))
    valid_null = null_max_t[np.isfinite(null_max_t)]
    p_values = np.full(observed_t.shape, np.nan, dtype=float)
    if valid_null.size:
        for idx, value in enumerate(observed_t):
            if np.isfinite(value):
                p_values[idx] = (1.0 + np.sum(valid_null >= value)) / (valid_null.size + 1.0)

    observed["spin_procedure_p"] = p_values
    observed["significant_spin_procedure"] = observed["spin_procedure_p"] < args.alpha
    observed = observed.sort_values(
        [
            "spin_procedure_p",
            "procedure_permutation_p" if "procedure_permutation_p" in observed.columns else "paired_t_p_fdr",
            "paired_t_p",
            "cohen_dz",
            "combination_size",
            "combination",
        ],
        ascending=[True, True, True, False, True, True],
        na_position="last",
    )
    return observed, pd.DataFrame(null_rows)


def combination_subject_values(summary, combination_row):
    usable = summary[
        np.isfinite(summary["inside_mean"])
        & np.isfinite(summary["outside_mean"])
    ].copy()
    metrics = [parse_metric_spec(spec) for spec in str(combination_row["metrics"]).split(";")]
    operators = str(combination_row["operator_pattern"]).split()
    orientation = int(combination_row.get("orientation", 1))

    metric_tables = {}
    for metric in metrics:
        feature, label, analysis = metric
        metric_tables[metric] = usable[
            (usable["feature"] == feature)
            & (usable["label"] == label)
            & (usable["analysis"] == analysis)
        ].set_index("subject")[["inside_mean", "outside_mean"]]

    common_subjects = set(metric_tables[metrics[0]].index)
    for metric in metrics[1:]:
        common_subjects &= set(metric_tables[metric].index)
    subjects = sorted(common_subjects)
    if not subjects:
        return pd.DataFrame(columns=["subject", "inside_mean", "outside_mean"])

    first_table = metric_tables[metrics[0]].loc[subjects]
    inside = first_table["inside_mean"].to_numpy(dtype=float)
    outside = first_table["outside_mean"].to_numpy(dtype=float)
    for operator, metric in zip(operators, metrics[1:]):
        table = metric_tables[metric].loc[subjects]
        inside = apply_operator(inside, table["inside_mean"].to_numpy(dtype=float), operator)
        outside = apply_operator(outside, table["outside_mean"].to_numpy(dtype=float), operator)

    inside = orientation * inside
    outside = orientation * outside
    valid = np.isfinite(inside) & np.isfinite(outside)
    return pd.DataFrame(
        {
            "subject": np.asarray(subjects)[valid],
            "inside_mean": inside[valid],
            "outside_mean": outside[valid],
        }
    )


def draw_inside_outside_plot(
    subset,
    metrics,
    output_path,
    title,
    tests=None,
    alpha=0.05,
    show_surface=False,
    group_sections=False,
):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    if subset.empty or not metrics:
        return None

    fig_width = max(9.5, len(metrics) * (1.25 if group_sections else 1.15))
    fig, ax = plt.subplots(figsize=(fig_width, 6.1 if group_sections else 5.8), facecolor="white")
    scale = 1.25 if group_sections else 1.0

    outside_color = "#385466"
    inside_color = "#d24d61"
    line_color = "#b8bdc2"
    mean_edge = "#111111"
    grid_color = "#e7e9eb"

    y_min = float(np.nanmin([subset["inside_mean"].min(), subset["outside_mean"].min()]))
    y_max = float(np.nanmax([subset["inside_mean"].max(), subset["outside_mean"].max()]))
    y_span = max(y_max - y_min, 1.0)

    if group_sections:
        blur_indices = [i for i, (feature, _, _) in enumerate(metrics) if feature.endswith("*blur")]
        if blur_indices:
            first_blur = min(blur_indices)
            ax.axvspan(first_blur - 0.5, len(metrics) - 0.5, color="#f8edf0", alpha=0.65, zorder=0)
            ax.axvline(first_blur - 0.5, color="#d8d8d8", lw=1.1, zorder=1)
            label_y = y_max + 0.075 * y_span
            if first_blur > 0:
                ax.text(
                    (first_blur - 1) / 2,
                    label_y,
                    "White surface metrics",
                    ha="center",
                    va="bottom",
                    fontsize=11.5,
                    color="#4b4f56",
                )
            ax.text(
                (first_blur + len(metrics) - 1) / 2,
                label_y,
                "Blur metrics",
                ha="center",
                va="bottom",
                fontsize=11.5,
                color="#4b4f56",
            )

    for i, (feature, label, analysis) in enumerate(metrics):
        feature_df = subset[
            (subset["feature"] == feature)
            & (subset["label"] == label)
            & (subset["analysis"] == analysis)
        ].sort_values("subject")
        x_outside = i - 0.18
        x_inside = i + 0.18

        for _, row in feature_df.iterrows():
            ax.plot(
                [x_outside, x_inside],
                [row["outside_mean"], row["inside_mean"]],
                color=line_color,
                lw=1.1,
                alpha=0.62,
                zorder=2,
            )

        jitter = (
            np.linspace(-0.04, 0.04, len(feature_df))
            if len(feature_df) > 1
            else np.array([0.0])
        )
        ax.scatter(
            np.full(len(feature_df), x_outside) + jitter,
            feature_df["outside_mean"],
            s=34,
            color=outside_color,
            edgecolor="white",
            linewidth=0.55,
            alpha=0.92,
            zorder=4,
        )
        ax.scatter(
            np.full(len(feature_df), x_inside) + jitter,
            feature_df["inside_mean"],
            s=34,
            color=inside_color,
            edgecolor="white",
            linewidth=0.55,
            alpha=0.94,
            zorder=4,
        )

        outside_mean = feature_df["outside_mean"].mean()
        inside_mean = feature_df["inside_mean"].mean()
        ax.plot(
            [x_outside - 0.095, x_outside + 0.095],
            [outside_mean, outside_mean],
            color=mean_edge,
            lw=2.8,
            solid_capstyle="round",
            zorder=5,
        )
        ax.plot(
            [x_inside - 0.095, x_inside + 0.095],
            [inside_mean, inside_mean],
            color=inside_color,
            lw=2.8,
            solid_capstyle="round",
            zorder=5,
        )

        test_row = get_test_row(tests, feature, label, analysis)
        if test_row is not None:
            q_value = float(test_row["paired_t_p_fdr"])
            stars = p_to_stars(q_value)
            if stars and q_value < alpha:
                feature_y_max = float(
                    np.nanmax([feature_df["inside_mean"].max(), feature_df["outside_mean"].max()])
                )
                ax.text(
                    i,
                    feature_y_max + 0.075 * y_span,
                    f"{stars}\nq={q_value:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=12,
                    fontweight="bold",
                    color=inside_color,
                    linespacing=0.9,
                )

    ax.axhline(0, color="#757575", lw=1.0, alpha=0.85, zorder=1)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(
        [format_metric_label(feature, label, show_surface=show_surface) for feature, label, _ in metrics],
        rotation=0 if show_surface else 30,
        ha="center" if show_surface else "right",
        fontsize=12 if group_sections else 10,
    )
    ax.tick_params(axis="y", labelsize=10 * scale)
    ax.set_ylabel("Mean cortical W-score", fontsize=13 if group_sections else 11)
    ax.set_xlabel("")
    ax.set_title(title, fontsize=17 if group_sections else 14, pad=34 if group_sections else 18)
    ax.grid(axis="y", color=grid_color, lw=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#333333")
    ax.spines["bottom"].set_color("#333333")
    ax.set_xlim(-0.6, len(metrics) - 0.4)
    ax.set_ylim(y_min - 0.12 * y_span, y_max + (0.28 if group_sections else 0.20) * y_span)

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=outside_color,
               markeredgecolor="white", markersize=7.5, label="Outside lesion"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=inside_color,
               markeredgecolor="white", markersize=7.5, label="Inside lesion"),
        Line2D([0], [0], color=mean_edge, lw=2.8, label="Group mean"),
    ]
    if group_sections:
        ax.legend(
            handles=legend_handles,
            frameon=False,
            ncol=3,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.10),
            fontsize=11.5,
            columnspacing=1.5,
            handlelength=1.8,
        )
    else:
        ax.legend(handles=legend_handles, frameon=False, ncol=3, loc="upper right", fontsize=10)

    fig.text(
        0.01,
        0.01,
        "Lines connect paired subject means; significance labels show paired t-test FDR q-values.",
        ha="left",
        va="bottom",
        fontsize=10 if group_sections else 8.5,
        color="#6a6a6a",
    )
    fig.tight_layout(rect=[0, 0.035, 1, 0.91 if group_sections else 1])
    fig.savefig(output_path, dpi=250)
    plt.close(fig)
    return output_path


def plot_inside_outside(summary, output_dir, tests=None, alpha=0.05):
    plot_paths = []
    required = {"inside_mean", "outside_mean", "label", "analysis"}
    if summary.empty or not required.issubset(summary.columns):
        return plot_paths

    for label in sorted(summary["label"].dropna().unique()):
        for analysis in sorted(summary["analysis"].dropna().unique()):
            subset = summary[
                (summary["label"] == label)
                & (summary["analysis"] == analysis)
                & np.isfinite(summary["inside_mean"])
                & np.isfinite(summary["outside_mean"])
            ].copy()
            if subset.empty:
                continue

            output_path = output_dir / f"lesion_vs_outside_label-{label}_analysis-{analysis}.png"
            plot_path = draw_inside_outside_plot(
                subset,
                ordered_metrics(subset),
                output_path,
                f"Lesion vs outside-lesion W-scores ({label}, {analysis})",
                tests=tests,
                alpha=alpha,
            )
            if plot_path:
                plot_paths.append(plot_path)

    return plot_paths


def plot_selected_metrics(summary, output_dir, tests=None, alpha=0.05):
    if summary.empty or not {"inside_mean", "outside_mean"}.issubset(summary.columns):
        return None
    subset = summary[
        np.isfinite(summary["inside_mean"])
        & np.isfinite(summary["outside_mean"])
    ].copy()
    if subset.empty:
        return None

    output_path = output_dir / "lesion_vs_outside_selected_metrics.png"
    return draw_inside_outside_plot(
        subset,
        ordered_metrics(subset, group_blur_last=True),
        output_path,
        "Lesion vs outside-lesion W-scores",
        tests=tests,
        alpha=alpha,
        show_surface=True,
        group_sections=True,
    )


def short_combination_label(combination):
    label = str(combination)
    replacements = {
        "[white,regional]": "",
        "[midthickness,regional]": "",
    }
    for feature in sorted(FEATURE_ABBREVIATIONS, key=len, reverse=True):
        replacements[feature] = FEATURE_ABBREVIATIONS[feature]
    for old, new in replacements.items():
        label = label.replace(old, new)
    label = label.replace(" -1 * ", "-")
    label = label.replace("(", "").replace(")", "")
    label = label.replace("-1 * ", "-")
    label = label.replace(" + ", "+").replace(" - ", "-")
    return label


def plot_significant_combinations(summary, combinations, output_dir, alpha=0.05, max_plots=12, spin_null=None):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    if combinations.empty:
        return None

    if "significant_spin_procedure" in combinations.columns:
        significant = combinations[combinations["significant_spin_procedure"]].copy()
    elif "significant_procedure_permutation" in combinations.columns:
        significant = combinations[combinations["significant_procedure_permutation"]].copy()
    else:
        significant = combinations[combinations["significant_paired_t_fdr"]].copy()
    if significant.empty:
        return None

    significant = significant.head(max_plots).reset_index(drop=True)
    combo_values = []
    for _, row in significant.iterrows():
        values = combination_subject_values(summary, row)
        combo_values.append(values)

    nonempty = [not values.empty for values in combo_values]
    if not any(nonempty):
        return None

    significant = significant.loc[nonempty].reset_index(drop=True)
    combo_values = [values for values, keep in zip(combo_values, nonempty) if keep]
    y_values = []
    for values in combo_values:
        y_values.extend(values["inside_mean"].tolist())
        y_values.extend(values["outside_mean"].tolist())
    y_values = np.asarray(y_values, dtype=float)
    y_min = float(np.nanmin(y_values))
    y_max = float(np.nanmax(y_values))
    y_span = max(y_max - y_min, 1.0)

    output_path = output_dir / "lesion_vs_outside_significant_combinations.png"
    fig_width = max(12.5, len(significant) * 2.05)
    has_spin_null = (
        spin_null is not None
        and not spin_null.empty
        and "max_abs_paired_t_stat" in spin_null.columns
        and np.isfinite(spin_null["max_abs_paired_t_stat"].to_numpy(dtype=float)).any()
    )
    if has_spin_null:
        fig, (ax, null_ax) = plt.subplots(
            1,
            2,
            figsize=(fig_width + 3.6, 6.9),
            facecolor="white",
            gridspec_kw={"width_ratios": [5.6, 1.55], "wspace": 0.24},
        )
    else:
        fig, ax = plt.subplots(figsize=(fig_width, 6.8), facecolor="white")
        null_ax = None

    outside_color = "#385466"
    inside_color = "#d24d61"
    line_color = "#b8bdc2"
    mean_edge = "#111111"
    grid_color = "#e7e9eb"
    observed_line_colors = [
        "#d24d61",
        "#2f6f95",
        "#7a5cc2",
        "#c97a2b",
        "#167a6f",
        "#9a4d78",
        "#5f6f30",
        "#3f5873",
    ]

    for i, (row, values) in enumerate(zip(significant.to_dict("records"), combo_values)):
        values = values.sort_values("subject")
        x_outside = i - 0.18
        x_inside = i + 0.18
        for _, subject_row in values.iterrows():
            ax.plot(
                [x_outside, x_inside],
                [subject_row["outside_mean"], subject_row["inside_mean"]],
                color=line_color,
                lw=1.1,
                alpha=0.62,
                zorder=2,
            )
        jitter = (
            np.linspace(-0.04, 0.04, len(values))
            if len(values) > 1
            else np.array([0.0])
        )
        ax.scatter(
            np.full(len(values), x_outside) + jitter,
            values["outside_mean"],
            s=34,
            color=outside_color,
            edgecolor="white",
            linewidth=0.55,
            alpha=0.92,
            zorder=4,
        )
        ax.scatter(
            np.full(len(values), x_inside) + jitter,
            values["inside_mean"],
            s=34,
            color=inside_color,
            edgecolor="white",
            linewidth=0.55,
            alpha=0.94,
            zorder=4,
        )

        outside_mean = values["outside_mean"].mean()
        inside_mean = values["inside_mean"].mean()
        ax.plot(
            [x_outside - 0.095, x_outside + 0.095],
            [outside_mean, outside_mean],
            color=mean_edge,
            lw=2.8,
            solid_capstyle="round",
            zorder=5,
        )
        ax.plot(
            [x_inside - 0.095, x_inside + 0.095],
            [inside_mean, inside_mean],
            color=inside_color,
            lw=2.8,
            solid_capstyle="round",
            zorder=5,
        )

        q_value = row.get("spin_procedure_p", row.get("procedure_permutation_p", row["paired_t_p_fdr"]))
        stars = p_to_stars(q_value)
        ax.text(
            i,
            max(values["inside_mean"].max(), values["outside_mean"].max()) + 0.075 * y_span,
            f"{stars}\np={q_value:.3f}",
            ha="center",
            va="bottom",
            fontsize=13,
            fontweight="bold",
            color=inside_color,
            linespacing=0.9,
        )

    ax.axhline(0, color="#757575", lw=1.0, alpha=0.85, zorder=1)
    ax.set_xticks(range(len(significant)))
    combo_labels = [
        f"{i + 1}. {short_combination_label(row['combination'])}"
        for i, (_, row) in enumerate(significant.iterrows())
    ]
    ax.set_xticklabels(
        combo_labels,
        rotation=0,
        ha="center",
        fontsize=15,
        linespacing=0.95,
    )
    ax.tick_params(axis="y", labelsize=13)
    ax.set_ylabel("Combined metric value", fontsize=15)
    ax.set_title("FDR-significant combined lesion metrics", fontsize=18, pad=22)
    ax.grid(axis="y", color=grid_color, lw=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#333333")
    ax.spines["bottom"].set_color("#333333")
    ax.set_xlim(-0.6, len(significant) - 0.4)
    ax.set_ylim(y_min - 0.12 * y_span, y_max + 0.26 * y_span)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=outside_color,
                   markeredgecolor="white", markersize=7.5, label="Outside lesion"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=inside_color,
                   markeredgecolor="white", markersize=7.5, label="Inside lesion"),
            Line2D([0], [0], color=mean_edge, lw=2.8, label="Group mean"),
        ],
        frameon=False,
        ncol=3,
        loc="upper right",
        fontsize=12.5,
    )

    if has_spin_null:
        null_values = spin_null["max_abs_paired_t_stat"].to_numpy(dtype=float)
        null_values = null_values[np.isfinite(null_values)]
        bins = min(32, max(10, int(np.sqrt(null_values.size) * 1.5)))
        null_ax.hist(
            null_values,
            bins=bins,
            color="#d9e0e7",
            edgecolor="white",
            linewidth=0.8,
            zorder=1,
        )
        observed_handles = []
        for i, (_, row) in enumerate(significant.iterrows()):
            observed_t = abs(float(row["paired_t_stat"]))
            color = observed_line_colors[i % len(observed_line_colors)]
            null_ax.axvline(observed_t, color=color, lw=2.5, alpha=0.95, zorder=3)
            y_top = null_ax.get_ylim()[1]
            null_ax.text(
                observed_t,
                y_top * 0.96,
                str(i + 1),
                ha="center",
                va="top",
                fontsize=12,
                fontweight="bold",
                color=color,
                bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": color, "lw": 1.0},
                zorder=4,
            )
            observed_handles.append(
                Line2D([0], [0], color=color, lw=2.5, label=f"{i + 1}")
            )

        null_ax.set_title("Spin max-|t| null", fontsize=15, pad=18)
        null_ax.set_xlabel("Max |paired t| across all tested combinations", fontsize=12.5)
        null_ax.set_ylabel("Spins", fontsize=12.5)
        null_ax.tick_params(axis="both", labelsize=11.5)
        null_ax.grid(axis="y", color=grid_color, lw=0.8)
        null_ax.set_axisbelow(True)
        null_ax.spines["top"].set_visible(False)
        null_ax.spines["right"].set_visible(False)
        null_ax.spines["left"].set_color("#333333")
        null_ax.spines["bottom"].set_color("#333333")
        null_ax.set_xlim(
            min(float(null_values.min()), float(np.abs(significant["paired_t_stat"]).min())) * 0.92,
            13,
        )
        null_ax.legend(
            handles=observed_handles,
            title="Observed",
            frameon=False,
            ncol=2,
            loc="upper right",
            fontsize=10.5,
            title_fontsize=11,
        )

    fig.text(
        0.01,
        0.01,
        "Only combinations passing the fsLR spin max-statistic test are shown. Right panel: spin max-|t| null; vertical lines are observed combinations.",
        ha="left",
        va="bottom",
        fontsize=11,
        color="#6a6a6a",
    )
    if has_spin_null:
        fig.subplots_adjust(left=0.06, right=0.985, top=0.86, bottom=0.17, wspace=0.24)
    else:
        fig.subplots_adjust(left=0.075, right=0.985, top=0.88, bottom=0.15)
    fig.savefig(output_path, dpi=250)
    plt.close(fig)
    return output_path


def best_subject_average_lesion_effect(summary):
    if summary.empty or "mean_difference_inside_minus_outside" not in summary.columns:
        return None

    subject_effects = (
        summary[np.isfinite(summary["mean_difference_inside_minus_outside"])]
        .groupby("subject", as_index=False)
        .agg(
            average_inside_minus_outside=("mean_difference_inside_minus_outside", "mean"),
            n_metrics=("mean_difference_inside_minus_outside", "count"),
        )
    )
    if subject_effects.empty:
        return None

    return subject_effects.sort_values(
        ["average_inside_minus_outside", "n_metrics", "subject"],
        ascending=[False, False, True],
    ).iloc[0]


def add_method_columns(frame, method, wscore_root):
    frame = frame.copy()
    frame.insert(0, "method", method)
    frame.insert(1, "wscore_root", str(wscore_root))
    return frame


def best_positive_row(frame, stat_column="paired_t_stat"):
    if frame is None or frame.empty or stat_column not in frame.columns:
        return None

    candidates = frame.copy()
    if "mean_difference_inside_minus_outside" in candidates.columns:
        candidates = candidates[candidates["mean_difference_inside_minus_outside"] > 0]
    candidates = candidates[np.isfinite(candidates[stat_column].to_numpy(dtype=float))]
    if candidates.empty:
        return None
    return candidates.sort_values(stat_column, ascending=False).iloc[0]


def count_true(frame, column):
    if frame is None or frame.empty or column not in frame.columns:
        return 0
    return int(frame[column].fillna(False).astype(bool).sum())


def method_detection_assessment(method_results, alpha=0.05):
    rows = []
    for result in method_results:
        method = result["method"]
        wscore_root = result["wscore_root"]
        summary = result["summary"]
        tests = result["tests"]
        combinations = result["combinations"]

        feature_candidates = tests.copy()
        if not feature_candidates.empty and "direction" in feature_candidates.columns:
            feature_candidates = feature_candidates[feature_candidates["direction"] == "inside_higher"]

        combination_candidates = combinations.copy()
        if not combination_candidates.empty and "mean_difference_inside_minus_outside" in combination_candidates.columns:
            combination_candidates = combination_candidates[
                combination_candidates["mean_difference_inside_minus_outside"] > 0
            ]

        best_feature = best_positive_row(feature_candidates)
        best_combination = best_positive_row(combination_candidates)
        best_subject = best_subject_average_lesion_effect(summary)
        detection_counts = aggregate_detection_counts(summary)
        auc = result.get("auc", pd.DataFrame())
        best_auc = None
        if auc is not None and not auc.empty and "auc" in auc.columns:
            best_auc = auc[np.isfinite(auc["auc"].to_numpy(dtype=float))]
            best_auc = best_auc.sort_values("auc", ascending=False).iloc[0] if not best_auc.empty else None

        n_subjects = int(summary["subject"].nunique()) if not summary.empty else 0
        n_features = (
            int(summary[["feature", "label", "analysis"]].drop_duplicates().shape[0])
            if not summary.empty
            else 0
        )
        n_feature_sig = count_true(feature_candidates, "significant_paired_t_fdr")
        n_combo_sig = count_true(combination_candidates, "significant_paired_t_fdr")
        n_perm_sig = count_true(combination_candidates, "significant_procedure_permutation")
        n_spin_sig = count_true(combination_candidates, "significant_spin_procedure")

        rows.append(
            {
                "method": method,
                "wscore_root": str(wscore_root),
                "n_subjects": n_subjects,
                "n_metrics": n_features,
                "summary_rows": int(len(summary)),
                "n_significant_features_fdr": n_feature_sig,
                "n_significant_combinations_fdr": n_combo_sig,
                "n_permutation_significant_combinations": n_perm_sig,
                "n_spin_significant_combinations": n_spin_sig,
                "total_true_positive_vertices": detection_counts["total_true_positive_vertices"],
                "total_false_positive_vertices": detection_counts["total_false_positive_vertices"],
                "total_false_negative_vertices": detection_counts["total_false_negative_vertices"],
                "total_true_negative_vertices": detection_counts["total_true_negative_vertices"],
                "sensitivity": detection_counts["sensitivity"],
                "false_positive_rate": detection_counts["false_positive_rate"],
                "precision": detection_counts["precision"],
                "dice": detection_counts["dice"],
                "best_auc": float(best_auc["auc"]) if best_auc is not None else np.nan,
                "best_auc_feature": best_auc["feature"] if best_auc is not None else "",
                "best_auc_label": best_auc["label"] if best_auc is not None else "",
                "best_auc_threshold": (
                    float(best_auc["best_detection_threshold"]) if best_auc is not None else np.nan
                ),
                "best_feature": best_feature["feature"] if best_feature is not None else "",
                "best_feature_label": best_feature["label"] if best_feature is not None else "",
                "best_feature_t": (
                    float(best_feature["paired_t_stat"]) if best_feature is not None else np.nan
                ),
                "best_feature_dz": (
                    float(best_feature["cohen_dz"]) if best_feature is not None else np.nan
                ),
                "best_feature_q": (
                    float(best_feature["paired_t_p_fdr"]) if best_feature is not None else np.nan
                ),
                "best_combination": (
                    best_combination["combination"] if best_combination is not None else ""
                ),
                "best_combination_t": (
                    float(best_combination["paired_t_stat"]) if best_combination is not None else np.nan
                ),
                "best_combination_dz": (
                    float(best_combination["cohen_dz"]) if best_combination is not None else np.nan
                ),
                "best_combination_q": (
                    float(best_combination["paired_t_p_fdr"]) if best_combination is not None else np.nan
                ),
                "best_combination_permutation_p": (
                    float(best_combination["procedure_permutation_p"])
                    if best_combination is not None
                    and "procedure_permutation_p" in best_combination
                    and np.isfinite(best_combination["procedure_permutation_p"])
                    else np.nan
                ),
                "best_combination_spin_p": (
                    float(best_combination["spin_procedure_p"])
                    if best_combination is not None
                    and "spin_procedure_p" in best_combination
                    and np.isfinite(best_combination["spin_procedure_p"])
                    else np.nan
                ),
                "best_subject": best_subject["subject"] if best_subject is not None else "",
                "best_subject_average_inside_minus_outside": (
                    float(best_subject["average_inside_minus_outside"])
                    if best_subject is not None
                    else np.nan
                ),
                "alpha": alpha,
            }
        )

    assessment = pd.DataFrame(rows)
    if assessment.empty:
        return assessment

    sort_columns = [
        "n_spin_significant_combinations",
        "n_permutation_significant_combinations",
        "n_significant_combinations_fdr",
        "n_significant_features_fdr",
        "best_auc",
        "dice",
        "false_positive_rate",
        "best_combination_t",
        "best_feature_t",
    ]
    assessment = assessment.sort_values(
        sort_columns,
        ascending=[False, False, False, False, False, False, True, False, False],
        na_position="last",
    ).reset_index(drop=True)
    assessment.insert(0, "rank", np.arange(1, len(assessment) + 1))
    return assessment


def subset_method_results_by_lesion_type(method_results, lesion_type, args):
    subset_results = []
    lesion_type = str(lesion_type)
    for result in method_results:
        summary = result["summary"]
        if summary.empty or "lesion_type" not in summary.columns:
            continue
        subset = summary[summary["lesion_type"].astype(str) == lesion_type].copy()
        if subset.empty:
            continue

        tests = compute_inside_outside_tests(subset, alpha=args.alpha)
        detection_performance = compute_detection_performance(subset)
        if args.run_combinations:
            combinations = compute_combination_tests(
                subset,
                max_size=args.max_combination_size,
                min_subjects=args.min_combination_subjects,
                alpha=args.alpha,
                operators=args.combination_operators,
            )
        else:
            combinations = pd.DataFrame()

        subset_results.append(
            {
                **result,
                "summary": subset,
                "tests": tests,
                "significant": tests[tests["significant_paired_t_fdr"]].copy()
                if not tests.empty
                else tests.copy(),
                "detection_performance": detection_performance,
                "auc": pd.DataFrame(),
                "combinations": combinations,
            }
        )
    return subset_results


def write_lesion_type_method_comparisons(method_results, args):
    written = []
    lesion_types = sorted(
        {
            str(lesion_type)
            for result in method_results
            if not result["summary"].empty and "lesion_type" in result["summary"].columns
            for lesion_type in result["summary"]["lesion_type"].dropna().unique()
        }
    )
    for lesion_type in lesion_types:
        subset_results = subset_method_results_by_lesion_type(method_results, lesion_type, args)
        assessment = method_detection_assessment(subset_results, alpha=args.alpha)
        if assessment.empty:
            continue
        safe_type = safe_method_name(lesion_type)
        assessment.insert(1, "lesion_type", lesion_type)
        assessment_path = args.output_dir / f"lesion_wscore_overlap_method_assessment_{safe_type}.csv"
        assessment.to_csv(assessment_path, index=False)
        plot_path = plot_method_comparison(
            assessment,
            args.output_dir,
            filename=f"lesion_wscore_overlap_method_comparison_{safe_type}.png",
            title=f"{lesion_type} lesion-detection strength by processed dataset",
        )
        written.append((lesion_type, assessment_path, plot_path))
    return written


def write_lesion_type_roc_outputs(method_results, args):
    if not args.run_roc:
        return []

    lesion_types = sorted(
        {
            str(record.get("lesion_type", "lesion"))
            for result in method_results
            for record in result.get("vertex_records", [])
        }
    )
    written = []
    thresholds = roc_thresholds_from_args(args)
    for lesion_type in lesion_types:
        roc_frames = []
        auc_frames = []
        for result in method_results:
            records = [
                record
                for record in result.get("vertex_records", [])
                if str(record.get("lesion_type", "lesion")) == lesion_type
            ]
            if not records:
                continue
            roc, auc = compute_roc_curves(
                records,
                args.lesion_threshold,
                detection_tail=args.detection_tail,
                thresholds=thresholds,
            )
            if not roc.empty:
                roc_frames.append(add_method_columns(roc, result["method"], result["wscore_root"]))
            if not auc.empty:
                auc_frames.append(add_method_columns(auc, result["method"], result["wscore_root"]))

        if not roc_frames or not auc_frames:
            continue

        safe_type = safe_method_name(lesion_type)
        combined_roc = pd.concat(roc_frames, ignore_index=True)
        combined_auc = pd.concat(auc_frames, ignore_index=True).sort_values(
            ["auc", "method", "feature"],
            ascending=[False, True, True],
            na_position="last",
        )
        roc_path = args.output_dir / f"lesion_wscore_overlap_roc_{safe_type}_all_methods.csv"
        auc_path = args.output_dir / f"lesion_wscore_overlap_auc_{safe_type}_all_methods.csv"
        combined_roc.to_csv(roc_path, index=False)
        combined_auc.to_csv(auc_path, index=False)
        plot_path = plot_roc_feature_subplots(
            combined_roc,
            combined_auc,
            args.output_dir / f"lesion_wscore_overlap_roc_{safe_type}_all_methods.png",
            f"{lesion_type} ROC curves across processed datasets",
            current_threshold=args.detection_threshold,
        )
        written.append((lesion_type, roc_path, auc_path, plot_path))
    return written


def plot_method_comparison(
    assessment,
    output_dir,
    filename="lesion_wscore_overlap_method_comparison.png",
    title="Lesion-detection strength by processed dataset",
):
    if assessment.empty:
        return None

    methods = assessment["method"].astype(str).tolist()
    x = np.arange(len(methods))
    has_combinations = (
        "best_combination_t" in assessment.columns
        and np.isfinite(assessment["best_combination_t"].to_numpy(dtype=float)).any()
    )
    output_path = output_dir / filename

    if not has_combinations:
        sensitivity = np.nan_to_num(assessment["sensitivity"].to_numpy(dtype=float), nan=0.0)
        dice = np.nan_to_num(assessment["dice"].to_numpy(dtype=float), nan=0.0)
        false_positive_rate = np.nan_to_num(assessment["false_positive_rate"].to_numpy(dtype=float), nan=0.0)

        fig, ax = plt.subplots(figsize=(max(8.0, len(methods) * 1.9), 5.5), facecolor="white")
        width = 0.24
        bars_sens = ax.bar(x - width, sensitivity, width, label="Sensitivity", color="#d24d61")
        bars_dice = ax.bar(x, dice, width, label="Dice", color="#385466")
        bars_fpr = ax.bar(x + width, false_positive_rate, width, label="False-positive rate", color="#9aa6ad")

        for bar, fp_count in zip(bars_fpr, assessment["total_false_positive_vertices"].to_numpy(dtype=int)):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.025,
                f"FP={fp_count}",
                ha="center",
                va="bottom",
                fontsize=10,
                color="#333333",
            )

        ax.set_xticks(x)
        ax.set_xticklabels(methods, fontsize=12)
        ax.tick_params(axis="y", labelsize=11)
        ax.set_ylabel("Thresholded detection metric", fontsize=12.5)
        ax.set_title(title, fontsize=15, pad=16)
        ax.grid(axis="y", color="#e7e9eb", lw=0.8)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, fontsize=11.5)
        ax.set_ylim(0, max(1.0, float(np.nanmax([sensitivity.max(), dice.max(), false_positive_rate.max()])) * 1.22))
        fig.text(
            0.01,
            0.01,
            "Thresholded detections use the configured W-score threshold and tail; FP labels count outside-lesion vertices.",
            ha="left",
            va="bottom",
            fontsize=9.5,
            color="#6a6a6a",
        )
        fig.tight_layout(rect=[0, 0.035, 1, 1])
        fig.savefig(output_path, dpi=250)
        plt.close(fig)
        return output_path

    best_feature_t = assessment["best_feature_t"].to_numpy(dtype=float)
    best_combination_t = assessment["best_combination_t"].to_numpy(dtype=float)
    best_feature_t = np.nan_to_num(best_feature_t, nan=0.0)
    best_combination_t = np.nan_to_num(best_combination_t, nan=0.0)

    fig, ax = plt.subplots(figsize=(max(7.5, len(methods) * 1.8), 5.4), facecolor="white")
    width = 0.36
    feature_color = "#385466"
    combo_color = "#d24d61"
    bars_feature = ax.bar(x - width / 2, best_feature_t, width, label="Best single metric", color=feature_color)
    bars_combo = ax.bar(x + width / 2, best_combination_t, width, label="Best combination", color=combo_color)

    y_max = max(1.0, float(np.nanmax([best_feature_t.max(), best_combination_t.max()])))
    for bars, column in [
        (bars_feature, "n_significant_features_fdr"),
        (bars_combo, "n_significant_combinations_fdr"),
    ]:
        for bar, count in zip(bars, assessment[column].to_numpy(dtype=int)):
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + 0.04 * y_max,
                f"n={count}",
                ha="center",
                va="bottom",
                fontsize=10.5,
                color="#333333",
            )

    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=12)
    ax.tick_params(axis="y", labelsize=11)
    ax.set_ylabel("Best positive paired t-statistic", fontsize=12.5)
    ax.set_title(title, fontsize=15, pad=16)
    ax.grid(axis="y", color="#e7e9eb", lw=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=11.5)
    ax.set_ylim(0, y_max * 1.24)
    fig.text(
        0.01,
        0.01,
        "Higher bars indicate stronger within-lesion vs outside-lesion separation; n labels count FDR-significant tests.",
        ha="left",
        va="bottom",
        fontsize=9.5,
        color="#6a6a6a",
    )
    fig.tight_layout(rect=[0, 0.035, 1, 1])
    fig.savefig(output_path, dpi=250)
    plt.close(fig)
    return output_path


def run_analysis_for_method(args, method, wscore_root, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)

    vertex_records, skipped = load_vertex_records(args, wscore_root=wscore_root)
    summary = summary_from_vertex_records(
        vertex_records,
        args.lesion_threshold,
        detection_threshold=args.detection_threshold,
        detection_tail=args.detection_tail,
    )

    summary_path = output_dir / "lesion_wscore_overlap_summary.csv"
    primary_path = output_dir / "lesion_wscore_overlap_selected_metrics.csv"
    skipped_path = output_dir / "lesion_wscore_overlap_skipped.csv"
    tests_path = output_dir / "lesion_wscore_overlap_feature_tests.csv"
    significant_path = output_dir / "lesion_wscore_overlap_significant_features.csv"
    detection_path = output_dir / "lesion_wscore_overlap_detection_performance.csv"
    roc_path = output_dir / "lesion_wscore_overlap_roc.csv"
    auc_path = output_dir / "lesion_wscore_overlap_auc.csv"
    combinations_path = output_dir / "lesion_wscore_overlap_feature_combinations.csv"
    permutation_null_path = output_dir / "lesion_wscore_overlap_combination_permutation_null.csv"
    spin_null_path = output_dir / "lesion_wscore_overlap_combination_spin_null.csv"

    summary.to_csv(summary_path, index=False)
    skipped.to_csv(skipped_path, index=False)
    tests = compute_inside_outside_tests(summary, alpha=args.alpha)
    tests.to_csv(tests_path, index=False)
    significant = tests[tests["significant_paired_t_fdr"]].copy() if not tests.empty else tests.copy()
    significant.to_csv(significant_path, index=False)
    detection_performance = compute_detection_performance(summary)
    detection_performance.to_csv(detection_path, index=False)
    if args.run_roc:
        roc, auc = compute_roc_curves(
            vertex_records,
            args.lesion_threshold,
            detection_tail=args.detection_tail,
            thresholds=roc_thresholds_from_args(args),
        )
    else:
        roc = pd.DataFrame()
        auc = pd.DataFrame()
    roc.to_csv(roc_path, index=False)
    auc.to_csv(auc_path, index=False)

    if args.run_combinations:
        combinations = compute_combination_tests(
            summary,
            max_size=args.max_combination_size,
            min_subjects=args.min_combination_subjects,
            alpha=args.alpha,
            operators=args.combination_operators,
        )
        combinations, permutation_null = add_max_statistic_permutation_pvalues(
            combinations,
            summary,
            max_size=args.max_combination_size,
            min_subjects=args.min_combination_subjects,
            alpha=args.alpha,
            operators=args.combination_operators,
            n_permutations=args.n_permutations,
            seed=args.permutation_seed,
            verbose=True,
        )
        combinations, spin_null = add_spin_max_statistic_pvalues(
            combinations,
            vertex_records,
            args,
            n_spins=args.n_spins,
            seed=args.spin_seed,
            verbose=True,
        )
    else:
        combinations = pd.DataFrame()
        permutation_null = pd.DataFrame()
        spin_null = pd.DataFrame()
    combinations.to_csv(combinations_path, index=False)
    permutation_null.to_csv(permutation_null_path, index=False)
    spin_null.to_csv(spin_null_path, index=False)

    primary = summary.copy()
    primary.to_csv(primary_path, index=False)

    plot_paths = plot_inside_outside(summary, output_dir, tests=tests, alpha=args.alpha)
    roc_plot = None
    if args.run_roc:
        roc_plot = plot_roc_curves(
            roc,
            auc,
            output_dir / "lesion_wscore_overlap_roc.png",
            f"ROC curves ({method})",
            max_curves=args.roc_max_curves,
            current_threshold=args.detection_threshold,
        )
        if roc_plot:
            plot_paths.insert(0, roc_plot)
    selected_plot = plot_selected_metrics(summary, output_dir, tests=tests, alpha=args.alpha)
    if selected_plot:
        plot_paths.insert(0, selected_plot)
    combination_plot = None
    if args.run_combinations:
        combination_plot = plot_significant_combinations(
            summary,
            combinations,
            output_dir,
            alpha=args.alpha,
            spin_null=spin_null,
        )
    if combination_plot:
        plot_paths.insert(0, combination_plot)

    return {
        "method": method,
        "wscore_root": wscore_root,
        "output_dir": output_dir,
        "summary": summary,
        "skipped": skipped,
        "tests": tests,
        "significant": significant,
        "detection_performance": detection_performance,
        "roc": roc,
        "auc": auc,
        "combinations": combinations,
        "permutation_null": permutation_null,
        "spin_null": spin_null,
        "vertex_records": vertex_records,
        "plot_paths": plot_paths,
        "combination_plot": combination_plot,
        "paths": {
            "summary": summary_path,
            "primary": primary_path,
            "skipped": skipped_path,
            "tests": tests_path,
            "significant": significant_path,
            "detection": detection_path,
            "roc": roc_path,
            "auc": auc_path,
            "combinations": combinations_path,
            "permutation_null": permutation_null_path,
            "spin_null": spin_null_path,
        },
    }


def print_analysis_result(args, result, lesion_count):
    method = result["method"]
    prefix = "" if method == "single" else f"[{method}] "
    summary = result["summary"]
    significant = result["significant"]
    detection_performance = result["detection_performance"]
    auc = result["auc"]
    combinations = result["combinations"]
    combination_plot = result["combination_plot"]
    plot_paths = result["plot_paths"]
    paths = result["paths"]

    print(f"{prefix}lesion masks: {lesion_count}")
    print(f"{prefix}summary rows: {len(summary)}")
    if not summary.empty:
        print(f"{prefix}features: {', '.join(ordered_features(summary['feature'].unique()))}")
    print(f"{prefix}wrote summary: {paths['summary']}")
    print(f"{prefix}wrote primary summary: {paths['primary']}")
    print(f"{prefix}wrote skipped: {paths['skipped']}")
    print(f"{prefix}wrote feature tests: {paths['tests']}")
    print(f"{prefix}wrote significant features: {paths['significant']}")
    print(f"{prefix}wrote detection performance: {paths['detection']}")
    print(f"{prefix}wrote ROC curve points: {paths['roc']}")
    print(f"{prefix}wrote AUC summary: {paths['auc']}")
    print(f"{prefix}wrote feature combinations: {paths['combinations']}")
    print(f"{prefix}wrote combination permutation null: {paths['permutation_null']}")
    print(f"{prefix}wrote combination spin null: {paths['spin_null']}")
    if not significant.empty:
        print(f"{prefix}significant features:")
        for _, row in significant.iterrows():
            print(
                "  "
                f"{row['label']} / {row['analysis']} / {row['feature']}: "
                f"{row['direction']}, "
                f"paired t q={row['paired_t_p_fdr']:.4g}"
            )
    best_subject = best_subject_average_lesion_effect(summary)
    if best_subject is not None:
        print(
            f"{prefix}best average lesion effect: "
            f"{best_subject['subject']} "
            f"(mean inside-outside={best_subject['average_inside_minus_outside']:.4g}, "
            f"metrics={int(best_subject['n_metrics'])})"
        )
    detection_counts = aggregate_detection_counts(summary)
    if np.isfinite(detection_counts["false_positive_rate"]):
        print(
            f"{prefix}false-positive detections: "
            f"{detection_counts['total_false_positive_vertices']} outside-lesion vertices "
            f"(FPR={detection_counts['false_positive_rate']:.4g}, "
            f"sensitivity={detection_counts['sensitivity']:.4g}, "
            f"dice={detection_counts['dice']:.4g})"
        )
    if not detection_performance.empty:
        best_detection = detection_performance.sort_values(
            ["dice", "false_positive_rate", "feature"],
            ascending=[False, True, True],
            na_position="last",
        ).iloc[0]
        print(
            f"{prefix}best thresholded feature: "
            f"{best_detection['label']} / {best_detection['analysis']} / {best_detection['feature']} "
            f"(dice={best_detection['dice']:.4g}, "
            f"FP vertices={int(best_detection['total_false_positive_vertices'])})"
        )
    if not auc.empty:
        best_auc = auc.sort_values("auc", ascending=False, na_position="last").iloc[0]
        print(
            f"{prefix}best ROC/AUC feature: "
            f"{best_auc['label']} / {best_auc['analysis']} / {best_auc['feature']} "
            f"(AUC={best_auc['auc']:.4g}, "
            f"best threshold={best_auc['best_detection_threshold']:.4g})"
        )
    if plot_paths:
        print(f"{prefix}plots:")
        for path in plot_paths:
            print(f"  {path}")
    if combination_plot is None:
        if args.run_combinations:
            print(f"{prefix}no FDR-significant combined metric plot was generated")
        else:
            print(f"{prefix}feature combinations were skipped")
    if not combinations.empty:
        print(f"{prefix}top feature combinations:")
        for _, row in combinations.head(10).iterrows():
            permutation_text = (
                f", permutation p={row['procedure_permutation_p']:.4g}"
                if "procedure_permutation_p" in row and np.isfinite(row["procedure_permutation_p"])
                else ""
            )
            spin_text = (
                f", spin p={row['spin_procedure_p']:.4g}"
                if "spin_procedure_p" in row and np.isfinite(row["spin_procedure_p"])
                else ""
            )
            print(
                "  "
                f"{row['combination']}: "
                f"dz={row['cohen_dz']:.3g}, "
                f"paired t q={row['paired_t_p_fdr']:.4g}"
                f"{permutation_text}"
                f"{spin_text}"
            )


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    mpl_config_dir = args.output_dir / ".matplotlib"
    mpl_config_dir.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))

    import matplotlib

    matplotlib.use("Agg")

    lesion_count = len(discover_lesions(args.lesion_dir, args.hemi_vertices))
    method_roots = parse_method_roots(args)

    if not args.compare_methods:
        result = run_analysis_for_method(args, "single", args.wscore_root, args.output_dir)
        print_analysis_result(args, result, lesion_count)
        return

    method_results = []
    used_output_names = set()
    for method, wscore_root in method_roots:
        safe_name = safe_method_name(method)
        output_name = safe_name
        suffix = 2
        while output_name in used_output_names:
            output_name = f"{safe_name}_{suffix}"
            suffix += 1
        used_output_names.add(output_name)

        print(f"\nRunning lesion-overlap analysis for {method}: {wscore_root}")
        if not Path(wscore_root).exists():
            print(f"  Warning: processed dataset root does not exist: {wscore_root}")
        result = run_analysis_for_method(
            args,
            method,
            Path(wscore_root),
            args.output_dir / output_name,
        )
        method_results.append(result)
        print_analysis_result(args, result, lesion_count)

    combined_specs = [
        ("summary", "lesion_wscore_overlap_summary_all_methods.csv"),
        ("skipped", "lesion_wscore_overlap_skipped_all_methods.csv"),
        ("tests", "lesion_wscore_overlap_feature_tests_all_methods.csv"),
        ("significant", "lesion_wscore_overlap_significant_features_all_methods.csv"),
        ("detection_performance", "lesion_wscore_overlap_detection_performance_all_methods.csv"),
        ("roc", "lesion_wscore_overlap_roc_all_methods.csv"),
        ("auc", "lesion_wscore_overlap_auc_all_methods.csv"),
        ("combinations", "lesion_wscore_overlap_feature_combinations_all_methods.csv"),
        ("permutation_null", "lesion_wscore_overlap_combination_permutation_null_all_methods.csv"),
        ("spin_null", "lesion_wscore_overlap_combination_spin_null_all_methods.csv"),
    ]
    for key, filename in combined_specs:
        frames = [
            add_method_columns(result[key], result["method"], result["wscore_root"])
            for result in method_results
            if not result[key].empty
        ]
        combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        combined_path = args.output_dir / filename
        combined.to_csv(combined_path, index=False)
        print(f"wrote combined {key}: {combined_path}")

    combined_roc_frames = [
        add_method_columns(result["roc"], result["method"], result["wscore_root"])
        for result in method_results
        if not result["roc"].empty
    ]
    combined_auc_frames = [
        add_method_columns(result["auc"], result["method"], result["wscore_root"])
        for result in method_results
        if not result["auc"].empty
    ]
    if args.run_roc and combined_roc_frames and combined_auc_frames:
        combined_roc = pd.concat(combined_roc_frames, ignore_index=True)
        combined_auc = pd.concat(combined_auc_frames, ignore_index=True).sort_values(
            ["auc", "method", "feature"],
            ascending=[False, True, True],
            na_position="last",
        )
        combined_roc_plot = plot_roc_feature_subplots(
            combined_roc,
            combined_auc,
            args.output_dir / "lesion_wscore_overlap_roc_all_methods.png",
            "ROC curves across processed datasets",
            current_threshold=args.detection_threshold,
        )
        if combined_roc_plot:
            print(f"wrote combined ROC plot: {combined_roc_plot}")

    lesion_type_roc_outputs = write_lesion_type_roc_outputs(method_results, args)
    for lesion_type, roc_path, auc_path, plot_path in lesion_type_roc_outputs:
        print(f"wrote {lesion_type} ROC curve points: {roc_path}")
        print(f"wrote {lesion_type} AUC summary: {auc_path}")
        if plot_path:
            print(f"wrote {lesion_type} ROC plot: {plot_path}")

    assessment = method_detection_assessment(method_results, alpha=args.alpha)
    assessment_path = args.output_dir / "lesion_wscore_overlap_method_assessment.csv"
    assessment.to_csv(assessment_path, index=False)
    comparison_plot = plot_method_comparison(assessment, args.output_dir)
    print(f"wrote method assessment: {assessment_path}")
    if comparison_plot:
        print(f"wrote method comparison plot: {comparison_plot}")

    lesion_type_outputs = write_lesion_type_method_comparisons(method_results, args)
    for lesion_type, lesion_assessment_path, lesion_plot_path in lesion_type_outputs:
        print(f"wrote {lesion_type} method assessment: {lesion_assessment_path}")
        if lesion_plot_path:
            print(f"wrote {lesion_type} method comparison plot: {lesion_plot_path}")

    if not assessment.empty:
        best = assessment.iloc[0]
        if args.run_combinations:
            print(
                "best lesion-separation dataset: "
                f"{best['method']} "
                f"(best combo t={best['best_combination_t']:.4g}, "
                f"combo FDR hits={int(best['n_significant_combinations_fdr'])}, "
                f"spin hits={int(best['n_spin_significant_combinations'])}, "
                f"FP vertices={int(best['total_false_positive_vertices'])})"
            )
        else:
            print(
                "best lesion-separation dataset: "
                f"{best['method']} "
                f"(dice={best['dice']:.4g}, "
                f"sensitivity={best['sensitivity']:.4g}, "
                f"FPR={best['false_positive_rate']:.4g}, "
                f"FP vertices={int(best['total_false_positive_vertices'])})"
            )


if __name__ == "__main__":
    main()
