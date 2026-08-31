"""Pre-dataset-norm whole-control and feature-specific control exclusion.

Covers the deterministic driver logic without heavy processing: exclusion
signature, backward-compatible base/output keying, pooled config-key separation,
the run_surface_qc manifest builder (fsLR-32k white intensity maps only, pid|ses
subject), subject-id round-tripping in _resolve_excluded_pairs (with run_surface_qc
mocked), and the restricted control dataset (subset threading). The heavy ANTs
composed run is validated separately in the Step 9 smoke.

Run: python tests/test_outlier_exclusion.py  (or under pytest).
"""
from __future__ import annotations

import sys
import types
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import pandas as pd

import zbrains_staged as zs
import zbrains_benchmark as zb


# --- signature + keying --------------------------------------------------------
def test_exclusion_signature_empty_stable_orderindependent_distinct():
    assert zs._exclusion_signature(frozenset()) == ""
    a = zs._exclusion_signature({("sub-1", "ses-01"), ("sub-2", "ses-01")})
    b = zs._exclusion_signature({("sub-2", "ses-01"), ("sub-1", "ses-01")})
    c = zs._exclusion_signature({("sub-3", "ses-01")})
    assert a == b                       # order-independent
    assert a != c                       # different sets differ
    assert a.startswith("excl") and len(a) == 14

    fa = zs._exclusion_signature({
        "T1w": {("sub-1", "ses-01")},
        "FLAIR": {("sub-2", "ses-01")},
    })
    fb = zs._exclusion_signature({
        "FLAIR": {("sub-2", "ses-01")},
        "T1w": {("sub-1", "ses-01")},
    })
    fc = zs._exclusion_signature({"T1w": {("sub-2", "ses-01")}})
    assert fa == fb and fa != fc


def test_base_and_output_dirs_backward_compatible_when_empty():
    # empty signature must reproduce the EXACT legacy path strings
    assert zb.processed_base_directory_for("wmmeanRavel", "PFX") == "PFX/zbrains_WMMEANRAVEL/"
    assert zb.processed_base_directory_for("wmmeanRavel", "PFX", "") == "PFX/zbrains_WMMEANRAVEL/"
    sig = zs._exclusion_signature({("sub-1", "ses-01")})
    assert zb.processed_base_directory_for("wmmeanRavel", "PFX", sig) == f"PFX/zbrains_WMMEANRAVEL_{sig}/"
    cfg, _ = zs.validate(dict(zs.SIMPLEST_BASELINE))
    axis = zs._axis_config(cfg)
    assert zb.output_directory_for(axis, "PFX", 0.5) == zb.output_directory_for(axis, "PFX", 0.5, "")
    assert zb.output_directory_for(axis, "PFX", 0.5, sig).rstrip("/").endswith(sig)


def test_config_key_separates_arms_and_preserves_baseline():
    base = dict(zs.SIMPLEST_BASELINE)
    k_base = zs._config_key(base, 0.5)
    k01 = zs._config_key(dict(base, outlier_method="robpca", qc_alpha=0.01), 0.5)
    k05 = zs._config_key(dict(base, outlier_method="robpca", qc_alpha=0.05), 0.5)
    assert len({k_base, k01, k05}) == 3
    # baseline hashes identically to a config that predates the Step-8 keys
    legacy = {k: v for k, v in base.items() if k not in ("outlier_method", "qc_alpha")}
    assert zs._config_key(legacy, 0.5) == k_base


def test_config_key_separates_correlation_quantiles():
    base = dict(zs.SIMPLEST_BASELINE, outlier_method="correlation",
                control_correlation_filter=True)
    keys = {
        zs._config_key(dict(base, control_correlation_quantile=q), 0.5)
        for q in (0.05, 0.10, 0.20)
    }
    assert len(keys) == 3


def test_base_signature_drops_exclusion_when_dataset_norm_none():
    # With dataset_norm='none' (the value while the ROBPCA outlier arms run) the
    # PROCESSED BASE must NOT be exclusion-keyed -- excluding a control changes no
    # subject's maps, so all arms share one base (no re-run of whitestripe/sampling).
    # The ANALYSIS dir must STILL carry the exclusion (its z-reference differs).
    sig = zs._exclusion_signature({("sub-1", "ses-01")})
    sig2 = zs._exclusion_signature({("sub-2", "ses-01")})
    assert sig and sig2 and sig != sig2

    none_cfg = dict(zs.SIMPLEST_BASELINE, dataset_norm="none")
    # base signature: exclusion dropped -> identical for ANY exclusion, == no-exclusion
    assert zs._base_signature(none_cfg, sig) == zs._base_signature(none_cfg, "")
    assert zs._base_signature(none_cfg, sig) == zs._base_signature(none_cfg, sig2)
    assert "excl" not in zs._base_signature(none_cfg, sig)
    # analysis signature: exclusion retained -> distinct per exclusion
    assert zs._processing_signature(none_cfg, sig) != zs._processing_signature(none_cfg, sig2)
    assert sig in zs._processing_signature(none_cfg, sig)
    # two ROBPCA arms therefore resolve to the SAME processed base dir ...
    assert (zb.processed_base_directory_for("none", "PFX", zs._base_signature(none_cfg, sig))
            == zb.processed_base_directory_for("none", "PFX", zs._base_signature(none_cfg, sig2)))

    # ... but with a dataset-level fit (RAVEL/Nyul) the exclusion re-keys the base
    # (the fit depends on which controls are in), so arms get DISTINCT bases.
    rav_cfg = dict(zs.SIMPLEST_BASELINE, dataset_norm="ravel")
    assert zs._base_signature(rav_cfg, sig) == zs._processing_signature(rav_cfg, sig)
    assert "excl" in zs._base_signature(rav_cfg, sig)
    assert zs._base_signature(rav_cfg, sig) != zs._base_signature(rav_cfg, sig2)


def test_ensure_base_uses_full_controls_when_dataset_norm_none(monkeypatch):
    # The dataset_norm='none' base is SHARED across exclusion arms, so it must be
    # built from ALL controls even when a restricted set is passed (exclusion applies
    # at analysis time). RAVEL/Nyul bases DO build from the restricted set.
    import types

    class _DS:
        def __init__(self, tag, pid):
            self.tag = tag
            self.demographics = types.SimpleNamespace(
                data=pd.DataFrame({"participant_id": [pid], "session_id": ["ses-01"]}))
        def process(self, **k): _seen.append(self.tag)
    _seen = []
    full = _DS("full", "c1"); restricted = _DS("restricted", "c1"); patient = _DS("patient", "p1")
    cohort = types.SimpleNamespace(name="MICs", control_dataset=full, patient_dataset=patient,
                                   output_dir_prefix="/pfx", base_pipeline_settings={})
    monkeypatch.setattr(zb, "processed_base_directory_for", lambda *a, **k: "/base")
    monkeypatch.setattr(zb, "processed_maps_complete", lambda *a, **k: False)  # force processing
    monkeypatch.setattr(zb, "mark_base_complete", lambda *a, **k: None)        # /base doesn't exist
    monkeypatch.setattr(zs, "_pipeline_settings", lambda *a, **k: {})

    # dataset_norm='none' -> FULL controls, even though a restricted set is passed
    _seen.clear()
    zs._ensure_base_processed(
        cohort, dict(zs.SIMPLEST_BASELINE, dataset_norm="none", normalization="whitestripe"),
        env=None, control_dataset=restricted, exclusion_signature="excl0badc0ffee", verbose=False)
    assert "full" in _seen and "restricted" not in _seen

    # dataset_norm='ravel' -> the RESTRICTED controls drive the fit
    _seen.clear()
    zs._ensure_base_processed(
        cohort, dict(zs.SIMPLEST_BASELINE, dataset_norm="ravel", normalization="whitestripeRavel"),
        env=None, control_dataset=restricted, exclusion_signature="excl0badc0ffee", verbose=False)
    assert "restricted" in _seen and "full" not in _seen


def test_robpca_arms_present_and_validate():
    stage = {s.name: s for s in zs.DEFAULT_STAGES}["outlier_detection"]
    labels = [lbl for lbl, _ in stage.candidates]
    assert "robpca_a0.01" in labels and "robpca_a0.05" in labels and "robpca_a0.10" in labels
    for arm in ({"outlier_method": "robpca", "qc_alpha": 0.01},
                {"outlier_method": "robpca", "qc_alpha": 0.05},
                {"outlier_method": "robpca", "qc_alpha": 0.10}):
        _, warns = zs.validate(dict(zs.SIMPLEST_BASELINE, **arm))
        assert not warns


# --- manifest builder ----------------------------------------------------------
def _touch(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")


def _fake_control(pairs):
    demo = types.SimpleNamespace(data=pd.DataFrame(
        {"participant_id": [p for p, _ in pairs], "session_id": [s for _, s in pairs]}))
    return types.SimpleNamespace(demographics=demo)


def test_control_manifest_selects_fslr32k_white_intensity_only(tmp_path):
    base = tmp_path / "zbrains_SELF"
    pairs = [("sub-HC1", "ses-01"), ("sub-HC2", "ses-01")]
    for pid, ses in pairs:
        cortex = base / pid / ses / "maps" / "cortex"
        for hemi in ("L", "R"):
            pfx = f"{pid}_{ses}_hemi-{hemi}"
            # intensity maps we WANT (fsLR-32k, label-white)
            for feat in ("T1w", "FLAIR", "qT1", "ADC", "FA"):
                _touch(cortex / f"{pfx}_surf-fsLR-32k_label-white_feature-{feat}_smooth-5mm.func.gii")
            # things we must NOT pick up:
            _touch(cortex / f"{pfx}_surf-fsLR-32k_label-midthickness_feature-T1w_smooth-5mm.func.gii")  # wrong label
            _touch(cortex / f"{pfx}_surf-fsLR-5k_label-white_feature-T1w_smooth-5mm.func.gii")           # wrong surf
            _touch(cortex / f"{pfx}_surf-fsLR-32k_label-thickness.func.gii")                             # no feature
            _touch(cortex / f"{pfx}_feature-T1w*blur_surf-fsnative_desc-raw.func.gii")                   # blur/fsnative
    man = zs._control_manifest(str(base), _fake_control(pairs))
    assert set(man.columns) >= {"subject", "feature", "hemisphere", "file_path"}
    # exactly 2 subjects x 5 features x 2 hemis
    assert len(man) == 2 * 5 * 2
    assert set(man["feature"]) == {"T1w", "FLAIR", "qT1", "ADC", "FA"}
    assert set(man["subject"]) == {"sub-HC1|ses-01", "sub-HC2|ses-01"}   # pid|ses encoding
    # exactly one row per (subject, feature, hemi)
    assert not man.duplicated(["subject", "feature", "hemisphere"]).any()


# --- exclusion resolution (run_surface_qc mocked) ------------------------------
def _cohort_stub(tmp_path, pairs):
    ctrl = _fake_control(pairs)
    return types.SimpleNamespace(
        name="MICs", output_dir_prefix=str(tmp_path), control_dataset=ctrl,
        base_pipeline_settings={"cortex": True}, patient_dataset=None)


def test_resolve_excluded_pairs_nonrobpca_is_empty(tmp_path):
    coh = _cohort_stub(tmp_path, [("sub-HC1", "ses-01")])
    assert zs._resolve_excluded_pairs(coh, dict(zs.SIMPLEST_BASELINE), env=None) == frozenset()
    assert zs._resolve_excluded_pairs(
        coh, dict(zs.SIMPLEST_BASELINE, control_correlation_filter=True), env=None) == frozenset()


def test_resolve_correlation_returns_feature_masks(tmp_path, monkeypatch):
    pairs = [(f"sub-HC{i}", "ses-01") for i in range(1, 5)]
    coh = _cohort_stub(tmp_path, pairs)
    monkeypatch.setattr(zb, "processed_maps_complete", lambda *a, **k: True)
    monkeypatch.setattr(
        zs, "_control_manifest",
        lambda *a, **k: pd.DataFrame([{"feature": "T1w"}]),
    )
    calls = []

    def fake_feature_exclusions(manifest, quantile, minimum_controls=None):
        calls.append((float(quantile), minimum_controls))
        return {
            "T1w": frozenset({("sub-HC2", "ses-01"), ("sub-BOGUS", "ses-99")}),
            "FLAIR": frozenset({("sub-HC3", "ses-01")}),
        }

    monkeypatch.setattr(zs, "_correlation_quantile_feature_exclusions",
                        fake_feature_exclusions)
    cfg = dict(zs.SIMPLEST_BASELINE, outlier_method="correlation",
               control_correlation_filter=True, control_correlation_quantile=0.10)
    zs._EXCLUSION_CACHE.clear()
    excluded = zs._resolve_excluded_pairs(coh, cfg, env=None, verbose=False)
    assert excluded == {
        "T1w": frozenset({("sub-HC2", "ses-01")}),
        "FLAIR": frozenset({("sub-HC3", "ses-01")}),
    }
    assert calls == [(0.10, zs.SURFACE_QC_MIN_CONTROLS)]


def test_correlation_quantiles_drop_exact_counts_per_feature(monkeypatch):
    import numpy as np
    import zbrains.surface_qc as sqc

    rng = np.random.RandomState(42)
    subjects = [f"sub-{i:02d}|ses-01" for i in range(20)]
    rows, arrays = [], {}
    for feature in ("T1w", "FLAIR"):
        base = rng.normal(size=80)
        order = range(20) if feature == "T1w" else reversed(range(20))
        for subject, noise_rank in zip(subjects, order):
            values = base + (0.02 + 0.03 * noise_rank) * rng.normal(size=80)
            for hemi, half in (("L", values[:40]), ("R", values[40:])):
                path = f"/{subject}/{feature}/{hemi}.func.gii"
                arrays[path] = half
                rows.append({"subject": subject, "feature": feature,
                             "hemisphere": hemi, "file_path": path})
    manifest = pd.DataFrame(rows)
    monkeypatch.setattr(sqc, "load_gifti_metric", lambda path: arrays[path])

    for quantile, expected in ((0.05, 1), (0.10, 2), (0.20, 4)):
        excluded = zs._correlation_quantile_feature_exclusions(
            manifest, quantile, minimum_controls=2)
        assert set(excluded) == {"T1w", "FLAIR"}
        assert all(len(pairs) == expected for pairs in excluded.values())
    # The two features rank independently rather than removing whole subjects.
    excluded = zs._correlation_quantile_feature_exclusions(
        manifest, 0.20, minimum_controls=2)
    assert excluded["T1w"] != excluded["FLAIR"]


def test_resolve_excluded_pairs_maps_subjects_and_intersects(tmp_path, monkeypatch):
    pairs = [("sub-HC1", "ses-01"), ("sub-HC2", "ses-01"), ("sub-HC3", "ses-01")]
    # build a base tree so the manifest is non-empty and processed_maps_complete=True.
    # The QC detector base is the subject-norm base at UNSMOOTHED 0mm (so the exclusion
    # is smoothing-independent), keyed by _processing_signature of a 0mm detector config.
    _det = dict(zs.SIMPLEST_BASELINE, outlier_method="robpca", qc_alpha=0.05)
    _det["cortical_smoothing"] = 0
    _det["hippocampal_smoothing"] = 0
    _sig = zs._processing_signature(_det, "")
    base = Path(zb.processed_base_directory_for("none", str(tmp_path), _sig))
    for pid, ses in pairs:
        cortex = base / pid / ses / "maps" / "cortex"
        for hemi in ("L", "R"):
            _touch(cortex / f"{pid}_{ses}_hemi-{hemi}_surf-fsLR-32k_label-white_feature-T1w_smooth-0mm.func.gii")
    monkeypatch.setattr(zb, "processed_maps_complete", lambda *a, **k: True)

    captured = {}
    def fake_run_surface_qc(manifest, output_dir=None, **kw):
        captured["subjects"] = sorted(manifest["subject"].unique())
        captured["qc_alpha"] = kw.get("qc_alpha")
        # flag HC2 (real) and a bogus subject not in the cohort (must be dropped)
        excl = pd.DataFrame({"subject": ["sub-HC2|ses-01", "sub-BOGUS|ses-99"]})
        return types.SimpleNamespace(exclusions=excl, table=None, summary={}, output_dir=output_dir)
    import zbrains.surface_qc as sqc
    monkeypatch.setattr(sqc, "run_surface_qc", fake_run_surface_qc)

    coh = _cohort_stub(tmp_path, pairs)
    cfg = dict(zs.SIMPLEST_BASELINE, outlier_method="robpca", qc_alpha=0.05)
    zs._EXCLUSION_CACHE.clear()
    excluded = zs._resolve_excluded_pairs(coh, cfg, env=None, verbose=False)
    assert excluded == frozenset({("sub-HC2", "ses-01")})   # bogus dropped by intersection
    assert captured["qc_alpha"] == 0.05
    assert captured["subjects"] == ["sub-HC1|ses-01", "sub-HC2|ses-01", "sub-HC3|ses-01"]
    # cached: a second call must NOT re-invoke run_surface_qc
    monkeypatch.setattr(sqc, "run_surface_qc",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError("re-ran QC")))
    assert zs._resolve_excluded_pairs(coh, cfg, env=None, verbose=False) == excluded


# --- restricted control dataset ------------------------------------------------
def test_restricted_control_dataset_identity_when_empty(tmp_path):
    coh = _cohort_stub(tmp_path, [("sub-HC1", "ses-01")])
    assert zs._restricted_control_dataset(coh, frozenset()) is coh.control_dataset


def test_restricted_control_dataset_threads_subset(monkeypatch, tmp_path):
    pairs = [("sub-HC1", "ses-01"), ("sub-HC2", "ses-01"), ("sub-HC3", "ses-01")]
    orig = types.SimpleNamespace(
        name="controls",
        demographics=types.SimpleNamespace(
            data=pd.DataFrame({"participant_id": [p for p, _ in pairs],
                               "session_id": [s for _, s in pairs]}),
            csv_file="ctrl.csv", column_mapping={"participant_id": "participant_id"},
            normative_columns=["AGE"], normative_dtypes=["int"], reference=None),
        micapipe_directory="/mica", hippunfold_directory=None, freesurfer_directory="/fs",
        raw_data_directory=None, cortex=True, hippocampus=False, subcortical=False)
    coh = types.SimpleNamespace(name="MICs", control_dataset=orig)

    seen = {}
    def fake_demo(csv, column_mapping=None, normative_columns=None, normative_dtypes=None,
                  reference=None, subset=None):
        seen["subset"] = subset
        return types.SimpleNamespace(_subset=subset)
    def fake_zbd(name, demo, mica, hippunfold_directory=None, freesurfer_directory=None,
                 raw_data_directory=None, cortex=True, hippocampus=True, subcortical=True):
        seen["name"] = name
        return types.SimpleNamespace(name=name, demographics=demo)
    import zbrains.dataset as ds
    monkeypatch.setattr(ds, "demographics", fake_demo)
    monkeypatch.setattr(ds, "zbdataset", fake_zbd)

    restricted = zs._restricted_control_dataset(coh, frozenset({("sub-HC2", "ses-01")}))
    assert seen["name"] == "controls"                       # control alias preserved
    assert sorted(map(tuple, seen["subset"])) == [("sub-HC1", "ses-01"), ("sub-HC3", "ses-01")]
    assert restricted.name == "controls"


def test_restricted_control_dataset_keeps_subjects_and_attaches_feature_mask(
        monkeypatch, tmp_path):
    pairs = [("sub-HC1", "ses-01"), ("sub-HC2", "ses-01"), ("sub-HC3", "ses-01")]
    orig = types.SimpleNamespace(
        name="controls",
        demographics=types.SimpleNamespace(
            data=pd.DataFrame({"participant_id": [p for p, _ in pairs],
                               "session_id": [s for _, s in pairs]}),
            csv_file="ctrl.csv", column_mapping={}, normative_columns=[],
            normative_dtypes=[], reference=None),
        micapipe_directory="/mica", hippunfold_directory=None,
        freesurfer_directory=None, raw_data_directory=None,
        cortex=True, hippocampus=False, subcortical=False,
        control_feature_exclusions={})
    coh = types.SimpleNamespace(name="MICs", control_dataset=orig)
    seen = {}

    def fake_demo(csv, **kwargs):
        seen["subset"] = kwargs["subset"]
        return types.SimpleNamespace(_subset=kwargs["subset"])

    def fake_zbd(name, demo, *args, **kwargs):
        return types.SimpleNamespace(name=name, demographics=demo,
                                     control_feature_exclusions={})

    import zbrains.dataset as ds
    monkeypatch.setattr(ds, "demographics", fake_demo)
    monkeypatch.setattr(ds, "zbdataset", fake_zbd)
    masks = {
        "T1w": frozenset({("sub-HC2", "ses-01")}),
        "FLAIR": frozenset({("sub-HC3", "ses-01")}),
    }
    restricted = zs._restricted_control_dataset(coh, masks)
    assert sorted(map(tuple, seen["subset"])) == pairs
    assert restricted.control_feature_exclusions == masks


def test_dataset_feature_mask_preserves_other_features_and_base_membership():
    from zbrains.dataset import zbdataset

    ds = zbdataset.__new__(zbdataset)
    ds.features = ["T1w", "T1w*blur", "FLAIR", "FA"]
    ds.control_feature_exclusions = {
        "T1w": {("sub-1", "ses-01")},
        "FLAIR": {("sub-2", "ses-01")},
    }
    pairs = [("sub-1", "ses-01"), ("sub-2", "ses-01"), ("sub-3", "ses-01")]
    valid = {
        "base": list(pairs),
        "structures": {"cortex": list(pairs)},
    }
    for feature in ds.features:
        valid[feature] = {
            "all": list(pairs),
            "structures": {"cortex": list(pairs), "hippocampus": [],
                           "subcortical": []},
        }

    ds._apply_control_feature_exclusions(valid)
    assert valid["base"] == pairs
    assert valid["structures"]["cortex"] == pairs
    assert ("sub-1", "ses-01") not in valid["T1w"]["all"]
    assert ("sub-1", "ses-01") not in valid["T1w*blur"]["all"]
    assert ("sub-1", "ses-01") in valid["FLAIR"]["all"]
    assert ("sub-2", "ses-01") not in valid["FLAIR"]["all"]
    assert valid["FA"]["all"] == pairs


def test_feature_masks_reach_ravel_fit_before_dataset_normalization(
        monkeypatch, tmp_path):
    """The RAVEL fitter receives a different kept-subject list per modality."""
    import os
    import zbrains.dataset as dataset_module
    from zbrains.dataset import zbdataset

    pairs = [("sub-1", "ses-01"), ("sub-2", "ses-01"), ("sub-3", "ses-01")]
    ds = zbdataset.__new__(zbdataset)
    ds.name = "controls"
    ds.features = []
    ds.control_feature_exclusions = {
        "T1w": {pairs[0]},
        "FLAIR": {pairs[1]},
    }
    ds.freesurfer_directory = None
    ds.hippunfold_directory = None
    ds.micapipe_directory = "/unused"
    ds.raw_data_directory = None
    ds.cortex = False
    ds.hippocampus = False
    ds.subcortical = False
    ds.hippunfold_version = 1

    def fake_add_features(*features, verbose=True):
        ds.features = list(features)
        ds.valid_subjects = {
            "base": list(pairs),
            "structures": {"cortex": [], "hippocampus": [], "subcortical": []},
        }
        for feature in features:
            ds.valid_subjects[feature] = {
                "all": list(pairs),
                "structures": {"cortex": [], "hippocampus": [], "subcortical": []},
            }
        ds._apply_control_feature_exclusions(ds.valid_subjects)

    ds.add_features = fake_add_features
    ds._copy_structural_files = lambda *a, **k: None
    ds._create_midline_from_freesurfer = lambda *a, **k: None
    ds.validate = lambda *a, **k: None
    monkeypatch.setattr(dataset_module, "_link_structural_to_cache", lambda *a, **k: None)
    monkeypatch.setattr(dataset_module, "prepare_t1w_flair_identity", lambda *a, **k: None)
    monkeypatch.setattr(dataset_module, "ensure_synthseg_csf", lambda *a, **k: None)
    monkeypatch.setattr(dataset_module, "_norm_cache_dir_for", lambda *a, **k: None)

    fitted = {}

    def fake_fit(*, subjects, modalities, **kwargs):
        fitted[modalities[0]] = list(subjects)

    monkeypatch.setattr(dataset_module, "fit_and_apply_ravel_to_controls", fake_fit)
    previous_cwd = os.getcwd()
    try:
        ds.process(
            output_directory=str(tmp_path / "ravel-base"),
            features=["T1w", "FLAIR"],
            normalization="noneRavel",
            env=types.SimpleNamespace(num_threads=1),
            n_jobs=1,
            verbose=False,
            skip_existing=False,
        )
    finally:
        os.chdir(previous_cwd)

    assert fitted["T1w"] == [pairs[1], pairs[2]]
    assert fitted["FLAIR"] == [pairs[0], pairs[2]]


if __name__ == "__main__":
    class _MP:
        """Restores patched attributes after each test so a monkeypatch can't leak
        into the next (e.g. a patched processed_base_directory_for)."""
        def __init__(self): self._undo = []
        def setattr(self, obj, name, val):
            self._undo.append((obj, name, getattr(obj, name)))
            setattr(obj, name, val)
        def undo(self):
            for obj, name, orig in reversed(self._undo):
                setattr(obj, name, orig)
            self._undo = []
    import tempfile
    names = sorted(n for n in dir() if n.startswith("test_"))
    passed = 0
    for name in names:
        fn = globals()[name]
        varn = fn.__code__.co_varnames[:fn.__code__.co_argcount]
        kwargs = {}
        if "tmp_path" in varn:
            kwargs["tmp_path"] = Path(tempfile.mkdtemp())
        mp = None
        if "monkeypatch" in varn:
            mp = _MP()
            kwargs["monkeypatch"] = mp
        try:
            fn(**kwargs)
        finally:
            if mp is not None:
                mp.undo()
        print(f"[ok] {name}")
        passed += 1
    print(f"\nALL {passed} OUTLIER-EXCLUSION TESTS PASSED")
