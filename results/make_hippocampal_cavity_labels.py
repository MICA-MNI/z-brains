#!/usr/bin/env python3
"""Project TLE resection-cavity masks onto hippocampal surfaces.

For each subject with a ``*_desc-cavity.nii.gz`` mask (in nativepro / T1w space),
this script samples the hippunfold hippocampal surfaces (inner, midthickness and
outer) at every vertex and marks a vertex as ``1`` whenever *any* of the three
surfaces passes through a labelled cavity voxel. The result is a per-hemisphere
binary ``.func.gii`` on the hippocampal mesh, written with the same naming
convention used by the cortical TLE cavity maps so that
``make_lesion_wscore_overlap*.py`` can pick it up as hippocampal ground truth.

Example (MICs, den-8k):

    python results/make_hippocampal_cavity_labels.py \
        --cavity-dir /host/verges/tank/data/Fatemeh/data/cavity_TLE/01.space-nativepro_tp-preop/BIDS_MICs \
        --hippunfold-dir /host/verges/tank/data/BIDS_MICs/derivatives/hippunfold_v2.0.0 \
        --density 8k \
        --output-dir results/lesiondata/MICs/TLE

Example (NOEL HippUnfold v2, den-8k):

    python results/make_hippocampal_cavity_labels.py \
        --cavity-dir /host/verges/tank/data/Fatemeh/data/cavity_TLE/01.space-nativepro_tp-preop/BIDS-NOEL \
        --hippunfold-dir /host/verges/tank/data/ian/BIDS_NOEL/derivatives/hippunfold_v2.0.0 \
        --density 8k \
        --output-dir results/lesiondata/NOEL/TLE
"""

import argparse
import re
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.gifti import GiftiDataArray, GiftiImage

SURFACES = ("inner", "midthickness", "outer")
CAVITY_PATTERN = re.compile(r"(?P<subject>sub-[^_]+)_(?P<session>ses-[^_]+)_.*cavity\.nii(\.gz)?$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--cavity-dir",
        type=Path,
        required=True,
        help="Directory containing *_desc-cavity.nii.gz masks (searched recursively).",
    )
    parser.add_argument(
        "--hippunfold-dir",
        type=Path,
        required=True,
        help="Root of the hippunfold derivatives (with or without a 'hippunfold/' subfolder).",
    )
    parser.add_argument(
        "--density",
        default="8k",
        help="Hippunfold surface density label, e.g. '8k' for HippUnfold v2.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Where to write per-hemisphere hippocampal cavity .func.gii files.",
    )
    parser.add_argument(
        "--cavity-threshold",
        type=float,
        default=0.0,
        help="A cavity voxel is labelled when its value is strictly above this threshold.",
    )
    parser.add_argument("--hemispheres", default="L,R", help="Comma-separated hemispheres to process.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate outputs even when the target .func.gii already exists.",
    )
    return parser.parse_args()


def discover_cavities(cavity_dir):
    """Return (subject, session, path) tuples for every cavity mask found."""
    records = []
    for path in sorted(Path(cavity_dir).glob("**/*cavity.nii*")):
        match = CAVITY_PATTERN.search(path.name)
        if not match:
            continue
        records.append((match.group("subject"), match.group("session"), path))
    return records


def find_surface(hippunfold_dir, subject, session, hemi, density, surface):
    """Locate a hippunfold T1w-space hippocampal surface, tolerating layout variants."""
    name = (
        f"{subject}_{session}_hemi-{hemi}_space-T1w_den-{density}_label-hipp_{surface}.surf.gii"
    )
    candidates = [
        Path(hippunfold_dir) / subject / session / "surf" / name,
        Path(hippunfold_dir) / "hippunfold" / subject / session / "surf" / name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    # Fall back to a glob in case the session folder is named differently.
    matches = sorted(Path(hippunfold_dir).glob(f"**/{name}"))
    return matches[0] if matches else None


def load_surface_coordinates(path):
    """Return the Nx3 POINTSET vertex coordinates from a .surf.gii.

    Hippunfold versions differ in array order: the coordinate (POINTSET,
    intent 1008) array may come after the triangle (intent 1009) array, so we
    select explicitly by intent rather than taking the first Nx3 array.
    """
    gifti = nib.load(str(path))
    pointset = [da for da in gifti.darrays if int(da.intent) == 1008]
    if pointset:
        return np.asarray(pointset[0].data, dtype=np.float64)
    # Fall back: the float Nx3 array is the coordinates; ints are triangles.
    for darray in gifti.darrays:
        data = np.asarray(darray.data)
        if data.ndim == 2 and data.shape[1] == 3 and np.issubdtype(data.dtype, np.floating):
            return data.astype(np.float64)
    raise ValueError(f"{path} does not contain an Nx3 POINTSET coordinate array")


def vertices_in_cavity(coords, cavity_data, inverse_affine, threshold):
    """Boolean per-vertex flag: nearest cavity voxel value is above threshold."""
    homogeneous = np.column_stack([coords, np.ones(len(coords))])
    voxels = homogeneous @ inverse_affine.T
    ijk = np.rint(voxels[:, :3]).astype(np.int64)

    shape = np.asarray(cavity_data.shape[:3])
    in_bounds = np.all((ijk >= 0) & (ijk < shape), axis=1)

    flags = np.zeros(len(coords), dtype=bool)
    valid = ijk[in_bounds]
    sampled = cavity_data[valid[:, 0], valid[:, 1], valid[:, 2]]
    flags[in_bounds] = sampled > threshold
    return flags


def write_func_gii(path, data):
    image = GiftiImage()
    image.add_gifti_data_array(
        GiftiDataArray(data=np.asarray(data, dtype=np.float32), intent="NIFTI_INTENT_NORMAL")
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    image.to_filename(str(path))


def process_subject(subject, session, cavity_path, args, hemispheres):
    cavity_img = nib.load(str(cavity_path))
    cavity_data = np.asarray(cavity_img.dataobj)
    inverse_affine = np.linalg.inv(cavity_img.affine)

    written = []
    for hemi in hemispheres:
        surface_flags = []
        missing = []
        for surface in SURFACES:
            surf_path = find_surface(
                args.hippunfold_dir, subject, session, hemi, args.density, surface
            )
            if surf_path is None:
                missing.append(surface)
                continue
            coords = load_surface_coordinates(surf_path)
            surface_flags.append(
                vertices_in_cavity(coords, cavity_data, inverse_affine, args.cavity_threshold)
            )

        if not surface_flags:
            print(f"  {subject} {session} hemi-{hemi}: no hippocampal surfaces found; skipping")
            continue
        if missing:
            print(f"  {subject} {session} hemi-{hemi}: missing surfaces {missing}; using {len(surface_flags)}")
        if len({flags.size for flags in surface_flags}) != 1:
            raise ValueError(
                f"{subject} {session} hemi-{hemi}: hippocampal surfaces have mismatched vertex counts"
            )

        # A vertex is positive if ANY surface passes through the cavity.
        cavity_label = np.any(np.vstack(surface_flags), axis=0).astype(np.float32)
        out_name = (
            f"{subject}_{session}_hemi-{hemi}_den-{args.density}_label-hipp_cavity.func.gii"
        )
        out_path = args.output_dir / out_name
        if out_path.exists() and not args.overwrite:
            print(f"  {subject} {session} hemi-{hemi}: exists, skipping ({out_path.name})")
            written.append(out_path)
            continue
        write_func_gii(out_path, cavity_label)
        print(
            f"  {subject} {session} hemi-{hemi}: {int(cavity_label.sum())}/{cavity_label.size} "
            f"vertices in cavity -> {out_path.name}"
        )
        written.append(out_path)
    return written


def main():
    args = parse_args()
    hemispheres = [h.strip() for h in args.hemispheres.split(",") if h.strip()]
    cavities = discover_cavities(args.cavity_dir)
    if not cavities:
        print(f"No cavity masks found under {args.cavity_dir}")
        return

    print(f"Found {len(cavities)} cavity masks under {args.cavity_dir}")
    total = 0
    for subject, session, cavity_path in cavities:
        print(f"{subject} {session}: {cavity_path.name}")
        total += len(process_subject(subject, session, cavity_path, args, hemispheres))
    print(f"\nWrote/verified {total} hippocampal cavity maps in {args.output_dir}")


if __name__ == "__main__":
    main()
