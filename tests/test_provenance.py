"""Analysis-output provenance manifest: stale {method}_maps from a prior run
(different code version / config / feature set) are cleared before re-use, so
orphaned maps can never be globbed into the objective. Base maps/structural kept.

Run: python tests/test_provenance.py  (or under pytest).
"""
from __future__ import annotations

import os
import sys
import tempfile
import types
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import zbrains_staged as st


def _fake_analysis_dir(tmp_path):
    out = tmp_path / "analysis"
    sub = out / "sub-01" / "ses-01"
    (sub / "maps").mkdir(parents=True)          # base (symlink in real runs) -- must be KEPT
    (sub / "structural").mkdir()                # base -- must be KEPT
    (sub / "zscore_maps" / "cortex").mkdir(parents=True)   # analysis output -- stale
    (sub / "wscore_maps").mkdir()                          # analysis output -- stale
    return str(out), sub


def test_guard_clears_stale_maps_on_mismatch_keeps_base(tmp_path):
    out, sub = _fake_analysis_dir(tmp_path)
    prov1 = {"code_version": "v1", "config": "c1", "features": ["T1w"]}

    # no manifest yet -> treated as stale: clears {method}_maps, keeps base
    assert st._guard_analysis_dir(out, prov1) is False
    assert not (sub / "zscore_maps").exists() and not (sub / "wscore_maps").exists()
    assert (sub / "maps").exists() and (sub / "structural").exists()

    st._stamp_analysis_dir(out, prov1)
    (sub / "zscore_maps").mkdir()               # a fresh analysis write

    # SAME provenance -> match, nothing cleared
    assert st._guard_analysis_dir(out, prov1) is True
    assert (sub / "zscore_maps").exists()

    # DIFFERENT code version -> stale, cleared
    assert st._guard_analysis_dir(out, {"code_version": "v2", "config": "c1", "features": ["T1w"]}) is False
    assert not (sub / "zscore_maps").exists()
    assert (sub / "maps").exists()              # base still kept


def test_analysis_provenance_distinguishes_config_features_and_code_stable():
    c_two = types.SimpleNamespace(base_pipeline_settings={"features": ["T1w", "FLAIR"]})
    c_one = types.SimpleNamespace(base_pipeline_settings={"features": ["T1w"]})
    p_z = st._analysis_provenance({"method": "zscore"}, c_two)
    p_w = st._analysis_provenance({"method": "wscore"}, c_two)
    p_z1 = st._analysis_provenance({"method": "zscore"}, c_one)
    assert p_z["config"] != p_w["config"]           # config change -> new provenance
    assert p_z["features"] != p_z1["features"]       # feature-set change -> new provenance
    assert p_z["features"] == ["FLAIR", "T1w"]       # sorted feature list recorded
    # code version is a stable 16-hex content hash of the map-generating source
    assert p_z["code_version"] == p_w["code_version"] == st._pipeline_code_version()
    assert len(p_z["code_version"]) == 16


def test_stamp_roundtrip_on_missing_dir_is_safe(tmp_path):
    # stamping a nonexistent dir is a no-op; guarding it returns False (no manifest)
    missing = str(tmp_path / "nope")
    st._stamp_analysis_dir(missing, {"code_version": "v", "config": "c", "features": []})
    assert st._guard_analysis_dir(missing, {"code_version": "v", "config": "c", "features": []}) is False


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        fn = globals()[name]
        kwargs = {}
        if "tmp_path" in fn.__code__.co_varnames[:fn.__code__.co_argcount]:
            kwargs["tmp_path"] = Path(tempfile.mkdtemp())
        fn(**kwargs)
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} PROVENANCE TESTS PASSED")
