"""Fix E: the greedy objective must score the HIPPOCAMPAL channel (cavity GT vs
hippocampal score maps) as its own cells, so TLE / mesial-temporal signal carries
non-zero weight. Builds a synthetic hipp cavity + hipp score map on disk and
asserts per_subject_feature_auc emits a hipp-<feature> row.

Run: python tests/test_hippocampal_objective.py  (or under pytest).
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


def _save_gii(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = nib.gifti.GiftiDataArray(data=np.asarray(data, dtype=np.float32),
                                   intent="NIFTI_INTENT_NORMAL")
    nib.save(nib.gifti.GiftiImage(darrays=[arr]), str(path))


def test_hippocampal_channel_scored_as_own_cell(monkeypatch):
    n = 60
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        lesion_dir = td / "lesions"
        score_root = td / "scores"
        sub, ses = "sub-PX999", "ses-01"

        # hippocampal cavity GT (left hemi): first 12 vertices are cavity
        cav = np.zeros(n, dtype=np.float32); cav[:12] = 1.0
        _save_gii(lesion_dir / sub / f"{sub}_{ses}_hemi-L_den-8k_label-hipp_cavity.func.gii", cav)

        # hippocampal score maps: high |z| inside the cavity -> AUC ~ 1
        hipp = score_root / sub / ses / "wscore_maps" / "hippocampus"
        lscore = np.zeros(n, dtype=np.float32); lscore[:12] = 6.0
        _save_gii(hipp / f"{sub}_{ses}_hemi-L_den-8k_label-hipp_midthickness_feature-T1w_smooth-2mm_analysis-regional.func.gii", lscore)
        _save_gii(hipp / f"{sub}_{ses}_hemi-R_den-8k_label-hipp_midthickness_feature-T1w_smooth-2mm_analysis-regional.func.gii", np.zeros(n, dtype=np.float32))

        # isolate the hippocampal path: no cortical lesions, lesion_dir -> synthetic
        monkeypatch.setattr(ga, "dataset_lesion_dir", lambda name: lesion_dir)
        monkeypatch.setattr(ga._LA, "discover_lesions", lambda *a, **k: [])

        frame = ga.per_subject_feature_auc(
            score_root, "MICs", method="wscore", features=["T1w"],
            analysis="regional", include_blur=False, include_hippocampus=True,
        )

    assert not frame.empty, "hippocampal channel produced no rows"
    hipp_rows = frame[frame["feature"] == "hipp-T1w"]
    assert len(hipp_rows) == 1, frame["feature"].tolist()
    row = hipp_rows.iloc[0]
    assert row["lesion_type"] == "TLE"
    assert row["n_inside"] == 12
    assert row["auc"] > 0.9                       # cavity vertices score highest


def test_hippocampus_off_produces_no_hipp_rows(monkeypatch):
    n = 40
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        lesion_dir = td / "lesions"; score_root = td / "scores"
        sub, ses = "sub-PX998", "ses-01"
        cav = np.zeros(n, dtype=np.float32); cav[:8] = 1.0
        _save_gii(lesion_dir / sub / f"{sub}_{ses}_hemi-L_den-8k_label-hipp_cavity.func.gii", cav)
        hipp = score_root / sub / ses / "wscore_maps" / "hippocampus"
        _save_gii(hipp / f"{sub}_{ses}_hemi-L_den-8k_label-hipp_midthickness_feature-T1w_smooth-2mm_analysis-regional.func.gii", cav * 5)
        monkeypatch.setattr(ga, "dataset_lesion_dir", lambda name: lesion_dir)
        monkeypatch.setattr(ga._LA, "discover_lesions", lambda *a, **k: [])
        frame = ga.per_subject_feature_auc(
            score_root, "MICs", method="wscore", features=["T1w"],
            analysis="regional", include_blur=False, include_hippocampus=False,
        )
    assert frame.empty or not (frame["feature"] == "hipp-T1w").any()


if __name__ == "__main__":
    class _MP:
        def setattr(self, obj, name, val): setattr(obj, name, val)
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name](_MP())
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} HIPPOCAMPAL-OBJECTIVE TESTS PASSED")
