"""Per-base completion sentinel: race-safe base reuse across processes/machines.

`_ensure_base_processed` stamps a `.zbrains_maps_complete` marker only after BOTH
controls and patients are verified present, checks it first (O(1) short-circuit), and
never stamps a partially-processed base -- so a peer machine sharing storage can trust
a base as complete without re-scanning and never treats a base still being written as
done. Pre-sentinel bases are adopted with zero reprocessing.

Run: python tests/test_base_sentinel.py  (or under pytest).
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

import pandas as pd

import zbrains_staged as zs
import zbrains_benchmark as zb


# --- sentinel primitives -------------------------------------------------------
def test_sentinel_roundtrip_and_subject_coverage(tmp_path):
    base = str(tmp_path / "zbrains_WB_smfwhm")
    os.makedirs(base)
    assert not zb.base_is_marked_complete(base)
    # stamp certifying two subjects
    zb.mark_base_complete(base, [("s1", "ses-01"), ("s2", "ses-01")])
    assert os.path.exists(os.path.join(base, zb._BASE_COMPLETE_SENTINEL))
    assert zb.base_is_marked_complete(base)                                   # existence
    assert zb.base_is_marked_complete(base, [("s1", "ses-01")])              # subset covered
    assert zb.base_is_marked_complete(base, [("s1", "ses-01"), ("s2", "ses-01")])
    # a subject NOT certified -> not complete for that requirement (CSV grew)
    assert not zb.base_is_marked_complete(base, [("s1", "ses-01"), ("s3", "ses-01")])
    zb.unmark_base_complete(base)
    assert not zb.base_is_marked_complete(base)
    zb.unmark_base_complete(base)                       # idempotent, no crash
    # marking a nonexistent dir is a safe no-op
    zb.mark_base_complete(str(tmp_path / "nope"), [("s1", "ses-01")])
    assert not zb.base_is_marked_complete(str(tmp_path / "nope"))


def test_symlink_path_repoints_stale_symlink(tmp_path):
    # symlink_path must REPOINT a symlink whose target changed (e.g. an analysis
    # dir's `maps` symlink after the base-signature change repointed exclusion arms
    # at the shared base) rather than raising -- but leave a REAL dir untouched.
    (tmp_path / "baseA" / "maps").mkdir(parents=True)
    (tmp_path / "baseB" / "maps").mkdir(parents=True)
    src_a = tmp_path / "baseA" / "maps"
    src_b = tmp_path / "baseB" / "maps"
    tgt = tmp_path / "analysis" / "maps"

    assert zb.symlink_path(src_a, tgt) is True
    assert tgt.is_symlink() and tgt.resolve() == src_a.resolve()
    # target now points to a DIFFERENT base -> repoint (no FileExistsError)
    assert zb.symlink_path(src_b, tgt) is True
    assert tgt.resolve() == src_b.resolve()
    # idempotent: same source again keeps it
    assert zb.symlink_path(src_b, tgt) is True and tgt.resolve() == src_b.resolve()
    # a REAL directory at the target is left as-is (never unlinked)
    real_tgt = tmp_path / "real" / "maps"
    real_tgt.mkdir(parents=True)
    assert zb.symlink_path(src_a, real_tgt) is True
    assert not real_tgt.is_symlink()


def _ds(tag, pairs, on_process):
    d = types.SimpleNamespace(demographics=types.SimpleNamespace(
        data=pd.DataFrame({"participant_id": [p for p, _ in pairs],
                           "session_id": [s for _, s in pairs]})))
    d.process = lambda **k: on_process(tag)
    return d


def _cohort(tmp_path, on_process, controls=(("c1", "ses-01"),), patients=(("p1", "ses-01"),)):
    return types.SimpleNamespace(
        name="MICs", control_dataset=_ds("control", controls, on_process),
        patient_dataset=_ds("patient", patients, on_process),
        output_dir_prefix=str(tmp_path), base_pipeline_settings={})


def _prep(tmp_path, monkeypatch):
    base = str(tmp_path / "zbrains_WB_smfwhm")
    os.makedirs(base, exist_ok=True)
    monkeypatch.setattr(zb, "processed_base_directory_for", lambda *a, **k: base)
    monkeypatch.setattr(zs, "_pipeline_settings", lambda *a, **k: {})
    return base


_NONE_CFG = None  # set per-call below


def test_ensure_base_stamps_after_build_and_adopts_without_reprocess(tmp_path, monkeypatch):
    procs = []
    base = _prep(tmp_path, monkeypatch)
    cohort = _cohort(tmp_path, procs.append)
    cfg = dict(zs.SIMPLEST_BASELINE, dataset_norm="none", normalization="whitestripe")
    # base already complete on disk (files present) but NOT stamped -> ADOPT: no
    # processing, just stamp.
    monkeypatch.setattr(zb, "processed_maps_complete", lambda *a, **k: True)
    zs._ensure_base_processed(cohort, cfg, env=None, verbose=False)
    assert procs == []                                  # nothing rebuilt
    assert zb.base_is_marked_complete(base)             # ... adopted + stamped


def test_ensure_base_shortcircuits_on_sentinel(tmp_path, monkeypatch):
    procs = []
    base = _prep(tmp_path, monkeypatch)
    cohort = _cohort(tmp_path, procs.append)
    cfg = dict(zs.SIMPLEST_BASELINE, dataset_norm="none", normalization="whitestripe")
    zb.mark_base_complete(base, [("c1", "ses-01"), ("p1", "ses-01")])   # a peer finished it
    # a stamped base (covering the current subjects) must be trusted in O(1):
    # processed_maps_complete is NOT consulted
    def boom(*a, **k):
        raise AssertionError("stamped base must not be re-scanned")
    monkeypatch.setattr(zb, "processed_maps_complete", boom)
    zs._ensure_base_processed(cohort, cfg, env=None, verbose=False)
    assert procs == []


def test_ensure_base_reprocesses_when_subject_added(tmp_path, monkeypatch):
    # A base stamped for {c1,p1} must NOT be trusted once the participant CSV gains a
    # subject (c2): the new subject isn't certified, so it re-processes + re-stamps.
    procs = []
    base = _prep(tmp_path, monkeypatch)
    cohort = _cohort(tmp_path, procs.append,
                     controls=(("c1", "ses-01"), ("c2", "ses-01")))   # c2 is new
    cfg = dict(zs.SIMPLEST_BASELINE, dataset_norm="none", normalization="whitestripe")
    zb.mark_base_complete(base, [("c1", "ses-01"), ("p1", "ses-01")])  # stale: no c2
    monkeypatch.setattr(zb, "processed_maps_complete", lambda *a, **k: True)
    zs._ensure_base_processed(cohort, cfg, env=None, verbose=False)
    # not short-circuited -> it ran the control processing path and re-stamped WITH c2
    assert zb.base_is_marked_complete(base, [("c1", "ses-01"), ("c2", "ses-01"), ("p1", "ses-01")])


def test_ensure_base_not_stamped_when_incomplete(tmp_path, monkeypatch):
    procs = []
    base = _prep(tmp_path, monkeypatch)
    cohort = _cohort(tmp_path, procs.append)
    cfg = dict(zs.SIMPLEST_BASELINE, dataset_norm="none", normalization="whitestripe")
    # controls present, patients NEVER complete -> a partially-built base must NOT be
    # stamped (so it is re-attempted next time, never falsely trusted).
    monkeypatch.setattr(zb, "processed_maps_complete",
                        lambda output_directory, datasets: datasets[0] is cohort.control_dataset)
    zs._ensure_base_processed(cohort, cfg, env=None, verbose=False)
    assert "patient" in procs                            # tried to build patients
    assert not zb.base_is_marked_complete(base)          # ... but not stamped (incomplete)


def test_ensure_base_reprocess_unmarks_then_restamps(tmp_path, monkeypatch):
    procs = []
    base = _prep(tmp_path, monkeypatch)
    cohort = _cohort(tmp_path, procs.append)
    cfg = dict(zs.SIMPLEST_BASELINE, dataset_norm="none", normalization="whitestripe")
    zb.mark_base_complete(base, [("c1", "ses-01"), ("p1", "ses-01")])
    monkeypatch.setattr(zb, "processed_maps_complete", lambda *a, **k: True)
    # reprocess_controls=True must NOT short-circuit; it rebuilds and re-stamps.
    zs._ensure_base_processed(cohort, cfg, env=None, verbose=False, reprocess_controls=True)
    assert "control" in procs                            # forced control reprocess
    assert zb.base_is_marked_complete(base)              # re-stamped after rebuild


if __name__ == "__main__":
    class _MP:
        def __init__(self): self._undo = []
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
        varn = fn.__code__.co_varnames[:fn.__code__.co_argcount]
        kwargs = {}
        if "tmp_path" in varn:
            kwargs["tmp_path"] = Path(tempfile.mkdtemp())
        mp = None
        if "monkeypatch" in varn:
            mp = _MP(); kwargs["monkeypatch"] = mp
        try:
            fn(**kwargs)
        finally:
            if mp is not None:
                mp.undo()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} BASE-SENTINEL TESTS PASSED")
