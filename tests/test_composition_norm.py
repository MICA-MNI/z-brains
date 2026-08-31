"""Composable normalization: subject-level (none/whitestripe/wmmean) composed
with dataset-level (none/ravel/nyul).

Covers the pure/deterministic layer of Step 7: canonical label composition, the
inverse decomposition, legacy single-mode aliases, composite-aware desc/path
resolution, backward-compatible descriptor threading, and unique per-composition
base directories in the benchmark driver. The heavy image-level RAVEL/Nyul run is
validated separately in the end-to-end composed smoke.

Run: python tests/test_composition_norm.py  (or under pytest).
"""
from __future__ import annotations

import inspect
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from zbrains import normalization as nz

_PAIRS = [(s, d) for s in ("none", "whitestripe", "wmmean")
          for d in ("none", "ravel", "nyul")]


def test_compose_decompose_roundtrips_all_nine_pairs():
    labels = []
    for s, d in _PAIRS:
        lbl = nz.compose_normalization_label(s, d)
        assert nz.decompose_normalization_label(lbl) == (s, d), (s, d, lbl)
        labels.append(lbl)
    assert len(set(labels)) == len(_PAIRS)  # every pair has a distinct desc


def test_legacy_single_mode_aliases_preserved():
    # The two historical pairs keep their legacy desc- names so existing
    # desc-ravel / desc-nyul bases and example scripts still work.
    assert nz.compose_normalization_label("whitestripe", "ravel") == "ravel"
    assert nz.compose_normalization_label("none", "nyul") == "nyul"
    assert nz.decompose_normalization_label("ravel") == ("whitestripe", "ravel")
    assert nz.decompose_normalization_label("nyul") == ("none", "nyul")
    # subject-only labels are unchanged
    for s in ("none", "whitestripe", "wmmean"):
        assert nz.compose_normalization_label(s, "none") == s
        assert nz.decompose_normalization_label(s) == (s, "none")


def test_alias_normalization_still_decomposes():
    for alias, expect in [("self", ("none", "none")), ("WS", ("whitestripe", "none")),
                          ("white-stripe", ("whitestripe", "none")),
                          ("tissuemean", ("wmmean", "none"))]:
        assert nz.decompose_normalization_label(alias) == expect


def test_resolve_desc_handles_composite_and_legacy():
    assert nz.resolve_normalization_desc("ravel") == "ravel"
    assert nz.resolve_normalization_desc("noneRavel") == "noneRavel"
    assert nz.resolve_normalization_desc("wmmeanNyul") == "wmmeanNyul"


def test_get_normalized_modality_path_uses_composite_desc():
    p = nz.get_normalized_modality_path("/out/sub-01/ses-01", "/mica", "T1w",
                                        normalization="wmmeanRavel")
    assert p.endswith("sub-01_ses-01_space-nativepro_desc-wmmeanRavel_T1w.nii.gz")
    # legacy label still resolves to its historical desc
    p2 = nz.get_normalized_modality_path("/out/sub-01/ses-01", "/mica", "FLAIR",
                                         normalization="ravel")
    assert p2.endswith("desc-ravel_FLAIR.nii.gz")


def test_invalid_pairs_rejected():
    for bad in [("bogus", "none"), ("none", "bogus")]:
        try:
            nz.compose_normalization_label(*bad)
        except ValueError:
            continue
        raise AssertionError(f"expected ValueError for {bad}")


def test_descriptor_threading_is_backward_compatible():
    # Legacy callers (no input_desc/output_desc) keep the historical defaults, so
    # existing single-mode runs are byte-for-byte unchanged.
    for fn, in_default, out_default in [
        (nz.fit_and_apply_ravel_to_controls, "whitestripe", "ravel"),
        (nz.apply_ravel_model_to_subject, "whitestripe", "ravel"),
        (nz.fit_and_apply_nyul_to_controls, None, "nyul"),
        (nz.apply_nyul_model_to_subject, None, "nyul"),
    ]:
        sig = inspect.signature(fn)
        assert sig.parameters["input_desc"].default == in_default, fn.__name__
        assert sig.parameters["output_desc"].default == out_default, fn.__name__


def test_driver_maps_pairs_to_unique_base_directories():
    import zbrains_staged as zs
    import zbrains_benchmark as zb
    bases = []
    for s, d in _PAIRS:
        cfg = dict(zs.SIMPLEST_BASELINE)
        cfg["subject_norm"] = s
        cfg["dataset_norm"] = d
        merged, warnings = zs.validate(cfg)
        assert not warnings, (s, d, warnings)
        assert merged["normalization"] == nz.compose_normalization_label(s, d)
        bases.append(zb.processed_base_directory_for(merged["normalization"], "PFX"))
    assert len(set(bases)) == len(_PAIRS)  # no base-dir collisions across arms


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} COMPOSITION-NORM TESTS PASSED")
