"""Base-independent structural cache: the structural/ dir (surfaces + Laplace/SWM)
is pure geometry, so it's generated once and symlinked into every base rather than
re-solved per normalization/smoothing/exclusion arm.

Run: python tests/test_structural_cache.py  (or under pytest).
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from zbrains.dataset import _link_structural_to_cache
from zbrains.hippunfold import hippunfold_cache_tag


def test_structural_symlinked_to_cache_and_reused_across_bases():
    with tempfile.TemporaryDirectory() as prefix:
        base = lambda b: os.path.join(prefix, b)  # noqa: E731

        # fresh base -> structural/ becomes a symlink into <prefix>/struct_cache/...
        a = _link_structural_to_cache(base("zbrains_SELF_smfwhm"), "sub-HC1", "ses-01")
        assert os.path.islink(a)
        cache = os.path.realpath(a)
        assert cache.endswith("struct_cache/sub-HC1/ses-01/structural")

        # a write "through" the symlink lands in the cache
        open(os.path.join(a, "surf.gii"), "w").write("geom")

        # a DIFFERENT base symlinks to the SAME cache -> reuses the geometry
        b = _link_structural_to_cache(base("zbrains_WMMEANRAVEL_smfwhm_exclx"), "sub-HC1", "ses-01")
        assert os.path.realpath(b) == cache
        assert open(os.path.join(b, "surf.gii")).read() == "geom"


def test_existing_real_structural_dir_left_untouched():
    with tempfile.TemporaryDirectory() as prefix:
        real = os.path.join(prefix, "zbrains_OLD", "sub-HC2", "ses-01", "structural")
        os.makedirs(real)
        open(os.path.join(real, "x"), "w").write("real")
        r = _link_structural_to_cache(os.path.join(prefix, "zbrains_OLD"), "sub-HC2", "ses-01")
        assert not os.path.islink(r) and os.path.isdir(r)   # real dir kept as-is


def test_no_zbrains_component_is_noop():
    with tempfile.TemporaryDirectory() as prefix:
        assert _link_structural_to_cache(os.path.join(prefix, "plain"), "sub-HC1", "ses-01") is None


def test_structural_cache_is_separated_by_hippunfold_source_and_version():
    with tempfile.TemporaryDirectory() as prefix:
        v1_source = os.path.join(prefix, "hippunfold_v1.4.1")
        v2_source = os.path.join(prefix, "hippunfold_v2.0.0")
        os.makedirs(v1_source)
        os.makedirs(v2_source)

        v1 = _link_structural_to_cache(
            os.path.join(prefix, "zbrains_V1"),
            "sub-01",
            "ses-01",
            hippunfold_directory=v1_source,
            hippunfold_version=1,
        )
        v2 = _link_structural_to_cache(
            os.path.join(prefix, "zbrains_V2"),
            "sub-01",
            "ses-01",
            hippunfold_directory=v2_source,
            hippunfold_version=2,
        )

        assert os.path.realpath(v1) != os.path.realpath(v2)
        assert hippunfold_cache_tag(v1_source, 1) in os.path.realpath(v1)
        assert hippunfold_cache_tag(v2_source, 2) in os.path.realpath(v2)
        open(os.path.join(v1, "v1-only.surf.gii"), "w").write("v1")
        assert not os.path.exists(os.path.join(v2, "v1-only.surf.gii"))


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} STRUCTURAL-CACHE TESTS PASSED")
