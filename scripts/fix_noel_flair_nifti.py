#!/usr/bin/env python3
"""Repair non-FCDPX NOEL raw FLAIR NIfTI-1 headers in place."""

from __future__ import annotations

import argparse
import gzip
import math
import os
import shutil
import struct
import sys
import tempfile
from pathlib import Path


DEFAULT_ROOT = Path("/host/verges/tank/data/ian/BIDS_NOEL/rawdata")


def repair_flair(src: Path) -> str:
    src_stat = src.stat()
    with gzip.open(src, "rb") as handle:
        raw = bytearray(handle.read())

    if len(raw) < 348:
        raise RuntimeError(f"file is too short for a NIfTI-1 header ({len(raw)} bytes)")

    if struct.unpack_from("<i", raw, 0)[0] == 348:
        endian = "<"
    elif struct.unpack_from(">i", raw, 0)[0] == 348:
        endian = ">"
    else:
        raise RuntimeError("file does not contain a valid NIfTI-1 header")

    magic = bytes(raw[344:348])
    if magic != b"n+1\x00":
        raise RuntimeError(
            f"expected single-file NIfTI magic b'n+1\\x00'; found {magic!r}"
        )

    dim = struct.unpack_from(endian + "8h", raw, 40)
    n_dims = dim[0]
    if not 1 <= n_dims <= 7:
        raise RuntimeError(f"invalid number of dimensions: {n_dims}")

    shape = dim[1 : n_dims + 1]
    bitpix = struct.unpack_from(endian + "h", raw, 72)[0]
    if bitpix <= 0 or bitpix % 8:
        raise RuntimeError(f"unexpected bitpix value: {bitpix}")

    data_bytes = math.prod(shape) * (bitpix // 8)
    expected_at_348 = 348 + data_bytes
    expected_at_352 = 352 + data_bytes

    if len(raw) == expected_at_348:
        raw[348:348] = b"\x00\x00\x00\x00"
        action = "inserted missing extension indicator and corrected vox_offset"
    elif len(raw) == expected_at_352:
        action = "kept extension indicator and corrected vox_offset"
    else:
        raise RuntimeError(
            f"size {len(raw)} does not match the expected simple NIfTI layout "
            f"({expected_at_348} or {expected_at_352} bytes)"
        )

    struct.pack_into(endian + "f", raw, 108, 352.0)

    # Stage beside the source so replacement is atomic on the same filesystem.
    fd, temp_name = tempfile.mkstemp(
        prefix=f".{src.name}.", suffix=".tmp", dir=src.parent
    )
    os.close(fd)
    temp = Path(temp_name)
    try:
        with gzip.open(temp, "wb") as handle:
            handle.write(raw)
        shutil.copystat(src, temp)
        os.chmod(temp, src_stat.st_mode)
        os.replace(temp, src)
    finally:
        temp.unlink(missing_ok=True)

    return action


def eligible_flairs(root: Path) -> tuple[list[Path], int]:
    candidates = sorted(
        path
        for path in root.glob("sub-*/**/anat/*.nii.gz")
        if path.name.lower().endswith("flair.nii.gz")
    )
    excluded = [path for path in candidates if "fcdpx" in str(path).lower()]
    included = [path for path in candidates if path not in excluded]
    return included, len(excluded)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help=f"NOEL rawdata root (default: {DEFAULT_ROOT})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="list eligible files without modifying them",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.root.is_dir():
        print(f"ERROR: rawdata root does not exist: {args.root}", file=sys.stderr)
        return 2

    files, excluded_count = eligible_flairs(args.root)
    print(
        f"Eligible raw FLAIR files: {len(files)}; "
        f"FCDPX files excluded: {excluded_count}",
        flush=True,
    )

    if args.dry_run:
        for src in files:
            print(src)
        return 0

    repaired = 0
    failures: list[tuple[Path, str]] = []
    for index, src in enumerate(files, start=1):
        try:
            action = repair_flair(src)
        except Exception as exc:  # Continue so one malformed input does not stop the batch.
            failures.append((src, str(exc)))
            print(f"[{index}/{len(files)}] ERROR {src}: {exc}", flush=True)
        else:
            repaired += 1
            print(
                f"[{index}/{len(files)}] OK {src.name}: {action}",
                flush=True,
            )

    print(
        f"Summary: repaired={repaired}, failed={len(failures)}, "
        f"fcdpx_excluded={excluded_count}",
        flush=True,
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
