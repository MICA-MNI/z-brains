"""Parameter-free specificity-aware objective: per_patient_vs_control_auc.

Positives = the lesion's |z|; negatives = ALL held-out-control vertices' |z| across
cortex + hippocampus (the patient's own non-lesion vertices are NOT used). No
threshold, no per-subject reduction. Verifies: (1) AUC ~1 when the lesion is far
more deviant than controls; (2) AUC collapses toward 0.5 when the pipeline
over-flags controls (the specificity penalty that the old within-subject objective
was blind to); (3) the negative pool includes hippocampal control vertices.

Run: python tests/test_control_negatives_objective.py  (or under pytest).
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np
import nibabel as nib

import results.greedy_auroc as ga

HEMI = ga.HEMI_VERTICES


def _save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = nib.gifti.GiftiDataArray(data=np.asarray(data, dtype=np.float32),
                                   intent="NIFTI_INTENT_NORMAL")
    nib.save(nib.gifti.GiftiImage(darrays=[arr]), str(path))


def _cortex(root, sub, feat, base, lesion_hi=None, method="wscore"):
    L = np.full(HEMI, base, np.float32); R = np.full(HEMI, base, np.float32)
    if lesion_hi is not None:
        L[:30] = lesion_hi                              # first 30 L vertices = lesion
    d = Path(root) / sub / "ses-01" / f"{method}_maps" / "cortex"
    for hemi, vec in (("L", L), ("R", R)):
        _save(d / f"{sub}_ses-01_hemi-{hemi}_surf-fsLR-32k_label-white_feature-{feat}_smooth-5mm_analysis-regional.func.gii", vec)


def _hipp(root, sub, feat, base, n=50, method="wscore"):
    d = Path(root) / sub / "ses-01" / f"{method}_maps" / "hippocampus"
    for hemi in ("L", "R"):
        _save(d / f"{sub}_ses-01_hemi-{hemi}_den-8k_label-hipp_midthickness_feature-{feat}_smooth-2mm_analysis-regional.func.gii",
              np.full(n, base, np.float32))


CTRLS = [("sub-HC1", "ses-01"), ("sub-HC2", "ses-01"), ("sub-HC3", "ses-01")]


def _lesion_record():
    lesion = np.zeros(2 * HEMI, dtype=np.float32); lesion[:30] = 1.0
    return [{"subject": "sub-PX1", "session": "ses-01", "lesion": lesion,
             "lesion_type": "FCD", "lesion_id": "sub-PX1"}]


def test_lesion_beats_controls_and_overflag_is_penalized(monkeypatch):
    with tempfile.TemporaryDirectory() as td:
        _cortex(td, "sub-PX1", "T1w", 0.4, lesion_hi=6.0)        # patient: lesion |z|=6
        for c, _ in CTRLS:
            _cortex(td, c, "T1w", 0.4)                           # healthy controls |z|~0.4
        monkeypatch.setattr(ga, "dataset_lesion_dir", lambda n: Path(td))
        monkeypatch.setattr(ga._LA, "discover_lesions", lambda *a, **k: _lesion_record())

        good = ga.per_patient_vs_control_auc(
            td, "MICs", control_root=td, control_subjects=CTRLS, method="wscore",
            features=["T1w"], include_hippocampus=False)
        assert not good.empty and good.iloc[0]["feature"] == "T1w"
        auc_good = float(good.iloc[0]["auc"])
        assert auc_good > 0.95                                   # lesion >> healthy controls

        # OVER-FLAG the controls (|z|=6 everywhere) -> negatives now as high as the
        # lesion -> AUC collapses. The old within-subject objective could not see this.
        for c, _ in CTRLS:
            _cortex(td, c, "T1w", 6.0)
        bad = ga.per_patient_vs_control_auc(
            td, "MICs", control_root=td, control_subjects=CTRLS, method="wscore",
            features=["T1w"], include_hippocampus=False)
        auc_bad = float(bad.iloc[0]["auc"])
        assert auc_bad < auc_good and auc_bad < 0.65, (auc_good, auc_bad)


def test_region_matched_reference_returns_controls_at_footprint():
    # The negatives are the controls' |z| at the SAME vertices as the lesion
    # footprint (region-matched), NOT the whole brain.
    with tempfile.TemporaryDirectory() as td:
        for c, _ in CTRLS:
            _cortex(td, c, "T1w", 0.4)
            _hipp(td, c, "T1w", 0.4, n=50)
        ref = ga._RegionMatchedReference(
            td, CTRLS, method="wscore", label="white", analysis="regional",
            include_blur=False, detection_tail="both", hipp_density=None)
        inside = np.zeros(2 * HEMI, dtype=bool); inside[:30] = True   # 30-vertex footprint
        neg = ref.cortex_negatives("T1w", inside)
        assert neg.size == len(CTRLS) * 30                # 3 controls x 30 footprint verts
        assert np.allclose(neg, 0.4, atol=1e-5)
        hinside = np.zeros(2 * 50, dtype=bool); hinside[:10] = True
        hneg = ref.hipp_negatives("T1w", hinside)
        assert hneg.size == len(CTRLS) * 10               # controls' hipp at the same verts


def test_hippocampal_cavity_emits_hipp_row_against_bare_key_pool(monkeypatch):
    # Regression: the hipp cavity is LABELLED hipp-<feature> but its negatives live
    # under the BARE <feature> key (cortex+hipp control vertices concatenated). A
    # key mismatch would silently drop every TLE/hipp row.
    with tempfile.TemporaryDirectory() as td:
        n = 50
        for c, _ in CTRLS:                                  # healthy control pool
            _cortex(td, c, "T1w", 0.4)
            _hipp(td, c, "T1w", 0.4, n=n)
        # patient: hippocampal cavity strongly deviant (first 20 L hipp vertices)
        L = np.full(n, 0.4, np.float32); L[:20] = 6.0
        R = np.full(n, 0.4, np.float32)
        d = Path(td) / "sub-PX1" / "ses-01" / "wscore_maps" / "hippocampus"
        for hemi, vec in (("L", L), ("R", R)):
            _save(d / f"sub-PX1_ses-01_hemi-{hemi}_den-8k_label-hipp_midthickness_feature-T1w_smooth-2mm_analysis-regional.func.gii", vec)
        cav_l = np.zeros(n, np.float32); cav_l[:20] = 1.0
        cav_r = np.zeros(n, np.float32)
        cdir = Path(td) / "cavities"; cdir.mkdir()
        _save(cdir / "L.func.gii", cav_l); _save(cdir / "R.func.gii", cav_r)

        monkeypatch.setattr(ga, "dataset_lesion_dir", lambda n: Path(td))
        monkeypatch.setattr(ga._LA, "discover_lesions", lambda *a, **k: [])   # no cortical lesion
        monkeypatch.setattr(ga._LA, "discover_hipp_cavities", lambda *a, **k: {
            ("sub-PX1", "ses-01", "L"): cdir / "L.func.gii",
            ("sub-PX1", "ses-01", "R"): cdir / "R.func.gii",
        })

        out = ga.per_patient_vs_control_auc(
            td, "MICs", control_root=td, control_subjects=CTRLS, method="wscore",
            features=["T1w"], include_hippocampus=True)
        hipp_rows = out[out["feature"] == "hipp-T1w"]
        assert len(hipp_rows) == 1, out           # the row is NOT dropped by a key mismatch
        assert hipp_rows.iloc[0]["lesion_type"] == "TLE"
        assert float(hipp_rows.iloc[0]["auc"]) > 0.9


def test_all_masked_real_lesion_scores_half_not_dropped(monkeypatch):
    # Regression: a REAL lesion whose inside vertices are all NaN (masked by a
    # filter) must be a MISS (AUC=0.5, masked_out=True), not a dropped row -- else
    # an arm could win by shrinking the evaluation set.
    with tempfile.TemporaryDirectory() as td:
        _cortex(td, "sub-PX1", "T1w", 0.4, lesion_hi=np.nan)   # lesion vertices masked
        for c, _ in CTRLS:
            _cortex(td, c, "T1w", 0.4)
        monkeypatch.setattr(ga, "dataset_lesion_dir", lambda n: Path(td))
        monkeypatch.setattr(ga._LA, "discover_lesions", lambda *a, **k: _lesion_record())

        out = ga.per_patient_vs_control_auc(
            td, "MICs", control_root=td, control_subjects=CTRLS, method="wscore",
            features=["T1w"], include_hippocampus=False)
        row = out[out["feature"] == "T1w"]
        assert len(row) == 1, out                  # scored, not dropped
        assert float(row.iloc[0]["auc"]) == 0.5
        assert bool(row.iloc[0]["masked_out"]) is True
        assert int(row.iloc[0]["n_inside"]) == 30  # GT lesion size preserved


def _disease_record(subject, lesion_type, hi):
    lesion = np.zeros(2 * HEMI, dtype=np.float32); lesion[:30] = 1.0
    return {"subject": subject, "session": "ses-01", "lesion": lesion,
            "lesion_type": lesion_type, "lesion_id": subject}


def _cortex_region(td, sub, hi, lo, region):
    """Write a cortex T1w map with |z|=hi over the L-hemi ``region`` slice, lo elsewhere."""
    L = np.full(HEMI, lo, np.float32); L[region] = hi
    R = np.full(HEMI, lo, np.float32)
    d = Path(td) / sub / "ses-01" / "wscore_maps" / "cortex"
    for hemi, vec in (("L", L), ("R", R)):
        _save(d / f"{sub}_ses-01_hemi-{hemi}_surf-fsLR-32k_label-white_feature-T1w_smooth-5mm_analysis-regional.func.gii", vec)


def _region_record(sub, lesion_type, region):
    les = np.zeros(2 * HEMI, dtype=np.float32); les[region] = 1.0
    return {"subject": sub, "session": "ses-01", "lesion": les,
            "lesion_type": lesion_type, "lesion_id": sub}


def test_disease_control_scores_each_disease_vs_the_other(monkeypatch):
    # Region-matched: an FCD lesion at verts 0..30 is ranked against the TLE
    # patient at THOSE SAME verts (healthy there), and a TLE lesion at 100..130
    # against the FCD patient there. Non-overlapping footprints -> valid negatives.
    with tempfile.TemporaryDirectory() as td:
        _cortex_region(td, "sub-FCD1", 6.0, 0.4, slice(0, 30))
        _cortex_region(td, "sub-TLE1", 5.0, 0.4, slice(100, 130))
        records = [_region_record("sub-FCD1", "FCD", slice(0, 30)),
                   _region_record("sub-TLE1", "TLE", slice(100, 130))]
        monkeypatch.setattr(ga, "dataset_lesion_dir", lambda n: Path(td))
        monkeypatch.setattr(ga._LA, "discover_lesions", lambda *a, **k: records)

        out = ga.per_patient_vs_other_disease_auc(
            td, "MICs", method="wscore", features=["T1w"], include_hippocampus=False)
        by_type = {r["lesion_type"]: r for _i, r in out.iterrows()}
        assert set(by_type) == {"FCD", "TLE"}, out
        assert float(by_type["FCD"]["auc"]) > 0.95   # FCD lesion vs TLE patient @ 0..30
        assert float(by_type["TLE"]["auc"]) > 0.95   # TLE lesion vs FCD patient @ 100..130


def test_region_matched_reference_excludes_own_lesion():
    # A reference subject's OWN lesion vertices must be excluded from the negatives
    # (so genuinely abnormal tissue never leaks into the "negative" class).
    with tempfile.TemporaryDirectory() as td:
        _cortex(td, "sub-FCD1", "T1w", 0.4, lesion_hi=6.0)   # lesion at first 30 = 6
        lesion = np.zeros(2 * HEMI, dtype=np.float32); lesion[:30] = 1.0
        ref = ga._RegionMatchedReference(
            td, [("sub-FCD1", "ses-01")], method="wscore", label="white",
            analysis="regional", include_blur=False, detection_tail="both",
            hipp_density=None, cortex_lesion={("sub-FCD1", "ses-01"): lesion})
        own = np.zeros(2 * HEMI, dtype=bool); own[:30] = True
        assert ref.cortex_negatives("T1w", own).size == 0        # own lesion -> no negatives
        elsewhere = np.zeros(2 * HEMI, dtype=bool); elsewhere[100:130] = True
        neg = ref.cortex_negatives("T1w", elsewhere)
        assert neg.size == 30 and float(neg.max()) < 1.0          # healthy tissue (~0.4)


def test_disease_control_single_disease_cohort_warns_and_drops(monkeypatch):
    # A cohort with only ONE disease has no OTHER-disease negative pool -> its
    # lesions must be dropped WITH a warning (not silently vanish).
    import warnings as _w
    with tempfile.TemporaryDirectory() as td:
        _cortex(td, "sub-FCD1", "T1w", 0.4, lesion_hi=6.0)
        records = [_disease_record("sub-FCD1", "FCD", 6.0)]   # FCD only, no TLE
        monkeypatch.setattr(ga, "dataset_lesion_dir", lambda n: Path(td))
        monkeypatch.setattr(ga._LA, "discover_lesions", lambda *a, **k: records)
        with _w.catch_warnings(record=True) as caught:
            _w.simplefilter("always")
            out = ga.per_patient_vs_other_disease_auc(
                td, "MICs", method="wscore", features=["T1w"], include_hippocampus=False)
        assert out.empty                                       # FCD lesions dropped (no TLE pool)
        assert any("no 'tle' patients" in str(x.message).lower() or
                   "no 'TLE' patients" in str(x.message) for x in caught), \
            [str(x.message) for x in caught]


def test_region_match_ignores_control_deviation_outside_footprint(monkeypatch):
    # Controls hot (|z|=6) EVERYWHERE EXCEPT the lesion footprint (0.4 there).
    # Region-matched, the lesion (6 at the footprint) still scores AUC~1 because
    # only the controls' footprint vertices are negatives -- their deviation
    # elsewhere is ignored. (The old whole-brain pool would have collapsed this.)
    with tempfile.TemporaryDirectory() as td:
        _cortex(td, "sub-PX1", "T1w", 0.4, lesion_hi=6.0)        # lesion first-30 = 6
        for c, _ in CTRLS:
            L = np.full(HEMI, 6.0, np.float32); L[:30] = 0.4     # hot except footprint
            R = np.full(HEMI, 6.0, np.float32)
            d = Path(td) / c / "ses-01" / "wscore_maps" / "cortex"
            for hemi, vec in (("L", L), ("R", R)):
                _save(d / f"{c}_ses-01_hemi-{hemi}_surf-fsLR-32k_label-white_feature-T1w_smooth-5mm_analysis-regional.func.gii", vec)
        monkeypatch.setattr(ga, "dataset_lesion_dir", lambda n: Path(td))
        monkeypatch.setattr(ga._LA, "discover_lesions", lambda *a, **k: _lesion_record())
        out = ga.per_patient_vs_control_auc(
            td, "MICs", control_root=td, control_subjects=CTRLS, method="wscore",
            features=["T1w"], include_hippocampus=False)
        assert float(out.iloc[0]["auc"]) > 0.95


def test_hipp_outside_mask_skips_on_mismatch_but_includes_when_no_cavity():
    # no cavity -> all hipp non-lesional (all-True); size mismatch -> None (skip).
    import numpy as np
    assert ga._hipp_outside_mask(None, 100).all()
    assert ga._hipp_outside_mask({}, 100).all()
    with tempfile.TemporaryDirectory() as td:
        cav = Path(td) / "cav_L.func.gii"
        _save(cav, np.ones(10, np.float32))                    # 10 (L) + 10 (R) = 20 != 100
        assert ga._hipp_outside_mask({"L": cav, "R": cav}, 100) is None
        assert ga._hipp_outside_mask({"L": cav, "R": cav}, 20) is not None  # aligns -> mask


if __name__ == "__main__":
    class _MP:
        def __init__(self):
            self._undo = []
        def setattr(self, obj, name, val):
            self._undo.append((obj, name, getattr(obj, name)))
            setattr(obj, name, val)
        def undo(self):
            for obj, name, orig in reversed(self._undo):
                setattr(obj, name, orig)
            self._undo = []
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        fn = globals()[name]
        mp = _MP()
        try:
            fn(mp) if "monkeypatch" in fn.__code__.co_varnames else fn()
        finally:
            mp.undo()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} CONTROL-NEGATIVES-OBJECTIVE TESTS PASSED")
