"""Norm-phase volume cache: reuse WhiteStripe/RAVEL/Nyul nativepro volumes across
smoothing/sampling arms WITHOUT re-running the (RAVEL-dominated) normalization, keyed
so cross-validation folds never share a fit (no leakage).

The volumes (desc-<final>_{T1w,FLAIR}.nii.gz) and the dataset-level fit model are
produced BEFORE any volume-to-surface sampling and on-surface smoothing, so they are
identical across every smoothing/sampling arm of the SAME (normalization, exclusion,
fold). The cache dir is derived from the base path with ONLY the smoothing token
stripped, so:
  * smoothing/sampling arms of one (norm, excl, fold) share ONE cache, but
  * a per-fold `_reffold{k}` base (train-only RAVEL fit) NEVER shares with the
    full-controls run or another fold -> no cross-validation leakage.

Run: python tests/test_norm_cache.py  (or under pytest).
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import zbrains.dataset as d


# --- keying / leakage guard ---------------------------------------------------
def test_cache_key_shared_across_smoothing_distinct_across_fold_and_exclusion():
    f = d._norm_cache_dir_for
    main = f("/p/zbrains_RAVEL_smfwhm_smoothctx0hip0")
    legacy = f("/p/zbrains_RAVEL_smfwhm")                   # old no-token baseline
    ctx10 = f("/p/zbrains_RAVEL_smfwhm_smoothctx10hip5")
    excl = f("/p/zbrains_RAVEL_smfwhm_smoothctx0hip0_exclABCD")
    fold3 = f("/p/zbrains_RAVEL_smfwhm_smoothctx0hip0_reffold3")
    fold4 = f("/p/zbrains_RAVEL_smfwhm_smoothctx2hip1_reffold4")

    assert main == legacy == ctx10                          # every smoothing arm shares
    assert "norm_cache_RAVEL_smfwhm" == os.path.basename(main)
    assert main != excl                                     # exclusion re-keys the fit
    assert main != fold3 and fold3 != fold4                 # each fold isolated (no leak)
    assert os.path.basename(fold3) == "norm_cache_RAVEL_smfwhm_reffold3"
    assert f("/p/some/plain/dir") is None                   # non-zbrains path -> no cache


# --- hydrate / save roundtrip -------------------------------------------------
def _write(path, text="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def _vol(root, pid, ses, desc, mod):
    return os.path.join(root, pid, ses, "maps",
                        d._norm_volume_filename(pid, ses, desc, mod))


def test_save_then_hydrate_symlinks_volumes_and_model(tmp_path):
    base1 = str(tmp_path / "zbrains_RAVEL_smfwhm_smoothctx5hip2")   # arm A (baseline-ish)
    base2 = str(tmp_path / "zbrains_RAVEL_smfwhm_smoothctx0hip2")   # arm B (different smoothing)
    cache = d._norm_cache_dir_for(base1)
    assert cache == d._norm_cache_dir_for(base2)                    # same cache for both arms

    subs = {("HC1", "ses-01"): ["T1w", "FLAIR"], ("HC2", "ses-01"): ["T1w"]}
    # arm A produced real volumes + a RAVEL fit model
    for (pid, ses), mods in subs.items():
        for m in mods:
            _write(_vol(base1, pid, ses, "ravel", m), f"{pid}{m}")
    _write(os.path.join(base1, "ravel", "T1w_ravel_model.npz"), "model")

    d._norm_cache_save(cache, base1, subs, "ravel")
    # cache now holds the volumes + model
    assert os.path.isfile(_vol(cache, "HC1", "ses-01", "ravel", "FLAIR"))
    assert os.path.isfile(os.path.join(cache, "ravel", "T1w_ravel_model.npz"))

    # arm B hydrates: volumes become symlinks into base2, model is copied in
    assert d._norm_cache_hydrate(cache, base2, subs, "ravel") is True
    hydrated = _vol(base2, "HC1", "ses-01", "ravel", "T1w")
    assert os.path.islink(hydrated)
    assert os.path.realpath(hydrated) == os.path.realpath(_vol(cache, "HC1", "ses-01", "ravel", "T1w"))
    assert os.path.isfile(os.path.join(base2, "ravel", "T1w_ravel_model.npz"))   # model restored


def test_hydrate_is_all_or_nothing_on_partial_cache(tmp_path):
    base = str(tmp_path / "zbrains_RAVEL_smfwhm_smoothctx0hip2")
    cache = d._norm_cache_dir_for(base)
    subs = {("HC1", "ses-01"): ["T1w", "FLAIR"]}
    # only ONE of the two required volumes is cached
    _write(_vol(cache, "HC1", "ses-01", "ravel", "T1w"), "t1")
    assert d._norm_cache_hydrate(cache, base, subs, "ravel") is False
    # nothing hydrated into the base (so the norm phase recomputes the full set)
    assert not os.path.exists(_vol(base, "HC1", "ses-01", "ravel", "T1w"))
    assert not os.path.exists(_vol(base, "HC1", "ses-01", "ravel", "FLAIR"))


def test_save_skips_symlinks_and_missing(tmp_path):
    # A base whose volumes are already cache SYMLINKS (i.e. it was hydrated) must not
    # be re-saved (no copying a symlink back over its own target), and missing/failed
    # subjects are skipped rather than erroring.
    base = str(tmp_path / "zbrains_RAVEL_smfwhm_smoothctx0hip2")
    cache = d._norm_cache_dir_for(base)
    subs = {("HC1", "ses-01"): ["T1w"], ("HC2", "ses-01"): ["T1w"]}   # HC2 will be "failed"
    _write(_vol(cache, "HC1", "ses-01", "ravel", "T1w"), "t1")
    # hydrate HC1 (becomes a symlink); HC2 absent -> all-or-nothing returns False,
    # so hydrate nothing; make HC1 a symlink manually to exercise the save skip
    os.makedirs(os.path.dirname(_vol(base, "HC1", "ses-01", "ravel", "T1w")), exist_ok=True)
    d._atomic_symlink(_vol(cache, "HC1", "ses-01", "ravel", "T1w"),
                      _vol(base, "HC1", "ses-01", "ravel", "T1w"))
    d._norm_cache_save(cache, base, subs, "ravel")           # must not raise
    # cache still has exactly HC1's single real file; HC2 never created
    assert os.path.isfile(_vol(cache, "HC1", "ses-01", "ravel", "T1w"))
    assert not os.path.islink(_vol(cache, "HC1", "ses-01", "ravel", "T1w"))
    assert not os.path.exists(_vol(cache, "HC2", "ses-01", "ravel", "T1w"))


def test_atomic_symlink_repoints_and_keeps_real_file(tmp_path):
    a = tmp_path / "a.nii.gz"; a.write_text("A")
    b = tmp_path / "b.nii.gz"; b.write_text("B")
    link = tmp_path / "link.nii.gz"
    assert d._atomic_symlink(str(a), str(link)) is True
    assert os.path.realpath(link) == os.path.realpath(a)
    # a REAL file at the destination is left untouched (never clobbered)
    real = tmp_path / "real.nii.gz"; real.write_text("R")
    assert d._atomic_symlink(str(a), str(real)) is True
    assert not os.path.islink(real) and real.read_text() == "R"


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        fn = globals()[name]
        kwargs = {}
        if "tmp_path" in fn.__code__.co_varnames[:fn.__code__.co_argcount]:
            kwargs["tmp_path"] = Path(tempfile.mkdtemp())
        fn(**kwargs)
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} NORM-CACHE TESTS PASSED")
