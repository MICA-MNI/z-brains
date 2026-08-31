"""'None' (0mm) smoothing support + separate cortical/hippocampal smoothing stages.

Covers processing._smooth_or_copy_metric (0/None FWHM -> copy the unsmoothed metric
instead of a zero-kernel wb_command that would produce no output) and the greedy
split of smoothing into two independent stages with the requested option sets.

Run: python tests/test_smoothing_none.py  (or under pytest).
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import zbrains_staged as st
from zbrains import processing as zp


# --- processing helper: FWHM=0 copies, FWHM>0 smooths --------------------------
def test_smooth_or_copy_metric_copies_when_fwhm_zero(monkeypatch):
    calls = []
    monkeypatch.setattr(zp.subprocess, "run", lambda *a, **k: calls.append(a))
    with tempfile.TemporaryDirectory() as d:
        src = os.path.join(d, "in.func.gii")
        open(src, "w").write("METRIC")
        # FWHM 0 -> copy the unsmoothed metric, NEVER call wb_command
        out0 = os.path.join(d, "out0.func.gii")
        zp._smooth_or_copy_metric("/wb", "/surf.gii", src, 0, out0)
        assert os.path.exists(out0) and open(out0).read() == "METRIC"
        # None -> also treated as no smoothing
        outN = os.path.join(d, "outN.func.gii")
        zp._smooth_or_copy_metric("/wb", "/surf.gii", src, None, outN)
        assert open(outN).read() == "METRIC"
        assert calls == []                       # no wb_command for the 0/None arms
        # FWHM > 0 -> runs wb_command -metric-smoothing ... -fwhm
        out5 = os.path.join(d, "out5.func.gii")
        zp._smooth_or_copy_metric("/wb", "/surf.gii", src, 5, out5)
        assert len(calls) == 1
        argv = calls[0][0]
        assert "-metric-smoothing" in argv and "5" in argv and "-fwhm" in argv


def test_smooth_or_copy_metric_missing_input_no_crash(monkeypatch):
    monkeypatch.setattr(zp.subprocess, "run", lambda *a, **k: None)
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "out.func.gii")
        zp._smooth_or_copy_metric("/wb", "/surf.gii", os.path.join(d, "nope.gii"), 0, out)
        assert not os.path.exists(out)           # nothing to copy -> no output, no crash


# --- greedy: separate cortical + hippocampal smoothing stages ------------------
def test_cortical_hippocampal_smoothing_separate_stages():
    names = [s.name for s in st.DEFAULT_STAGES]
    assert "smoothing" not in names                         # old combined stage gone
    i, j = names.index("cortical_smoothing"), names.index("hippocampal_smoothing")
    assert j == i + 1                                       # hippocampal right after cortical
    cs = [s for s in st.DEFAULT_STAGES if s.name == "cortical_smoothing"][0]
    hs = [s for s in st.DEFAULT_STAGES if s.name == "hippocampal_smoothing"][0]
    # explicit arms exclude the baseline (implicit "keep": 0mm for both)
    assert [o["cortical_smoothing"] for _, o in cs.candidates] == [2, 5, 10]
    assert [o["hippocampal_smoothing"] for _, o in hs.candidates] == [1, 2, 5]
    # each stage touches ONLY its own key
    assert all(list(o) == ["cortical_smoothing"] for _, o in cs.candidates)
    assert all(list(o) == ["hippocampal_smoothing"] for _, o in hs.candidates)


def test_smoothing_options_validate_including_none():
    assert st.CORTICAL_SMOOTHING_OPTIONS == (0, 2, 5, 10)
    assert st.HIPPOCAMPAL_SMOOTHING_OPTIONS == (0, 1, 2, 5)
    for v in st.CORTICAL_SMOOTHING_OPTIONS:
        _, w = st.validate(dict(st.SIMPLEST_BASELINE, cortical_smoothing=v))
        assert not w
    for v in st.HIPPOCAMPAL_SMOOTHING_OPTIONS:
        _, w = st.validate(dict(st.SIMPLEST_BASELINE, hippocampal_smoothing=v))
        assert not w
    # out-of-range values are rejected
    for bad in (dict(cortical_smoothing=3), dict(hippocampal_smoothing=10)):
        try:
            st.validate(dict(st.SIMPLEST_BASELINE, **bad))
            assert False, "expected ValueError"
        except ValueError:
            pass
    # 0mm ('None') is now the baseline for the whole optimization; positive
    # kernels key distinct processed bases.
    baseline, _ = st.validate(dict(st.SIMPLEST_BASELINE))
    assert baseline["cortical_smoothing"] == 0
    assert baseline["hippocampal_smoothing"] == 0
    assert st._processing_signature(baseline) == f"{st._SMOOTHING_CONVENTION}_smoothctx0hip0"
    smoothed, _ = st.validate(dict(st.SIMPLEST_BASELINE, cortical_smoothing=2))
    assert "smoothctx2hip0" in st._processing_signature(smoothed)


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
        mp = _MP()
        try:
            fn(mp) if "monkeypatch" in fn.__code__.co_varnames else fn()
        finally:
            mp.undo()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} SMOOTHING-NONE TESTS PASSED")
