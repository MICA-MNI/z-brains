"""_cortex_bool_mask: build the fsnative cortex/medial-wall mask robustly.

FreeSurfer's ?h.cortex.label indices must be turned into a length-n_vertices
boolean mask aligned to ?h.white. When a subject's label carries an index at or
past the surface's last vertex (off-by-one, or a surface/label recon mismatch),
the mask must be built anyway (dropping the bad index) rather than raising an
IndexError that drops the WHOLE subject from processing.

Run: python tests/test_cortex_mask.py  (or under pytest).
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np

from zbrains.dataset import _cortex_bool_mask


def test_in_range_indices_set_true():
    m = _cortex_bool_mask([0, 2, 4], 5, hemi="L", label="sub/ses")
    assert m.dtype == bool and m.shape == (5,)
    assert m.tolist() == [True, False, True, False, True]


def test_index_one_past_end_is_dropped_not_crash():
    # the actual failure mode: an index == n_vertices (216482 for size 216482)
    m = _cortex_bool_mask([0, 1, 5], 5, hemi="L", label="sub/ses")   # 5 is out of range
    assert m.shape == (5,)                       # stays aligned to the surface
    assert m.tolist() == [True, True, False, False, False]


def test_grossly_out_of_range_and_negative_dropped():
    m = _cortex_bool_mask([-1, 3, 99], 5)
    assert m.shape == (5,) and m.tolist() == [False, False, False, True, False]


def test_empty_label_gives_all_medial_wall():
    m = _cortex_bool_mask([], 4)
    assert m.shape == (4,) and not m.any()


def test_full_cortex_no_medial_wall():
    m = _cortex_bool_mask(np.arange(6), 6)
    assert m.all()


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} CORTEX-MASK TESTS PASSED")
