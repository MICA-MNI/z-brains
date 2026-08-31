"""_ThreadLocalStdout: a per-thread stdout proxy that must NEVER crash the compute.

Subject processing routes every print() through this proxy while running subjects
in worker threads. A logging hiccup (a proxy reconstructed without __init__ via
deepcopy/pickle, or read mid-^C-teardown) previously raised a cryptic
AttributeError('_local') out of a print(), which was caught by the per-subject
handler and SILENTLY DROPPED THE WHOLE SUBJECT from the base. The proxy must
degrade gracefully to the real stdout instead.

Run: python tests/test_stdout_proxy.py  (or under pytest).
"""
from __future__ import annotations

import io
import sys
import threading
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from zbrains.dataset import _ThreadLocalStdout


def test_default_routing_and_per_thread_target():
    default = io.StringIO()
    p = _ThreadLocalStdout(default)
    p.write("main-default\n")                    # no target -> default
    assert default.getvalue() == "main-default\n"

    tgt = io.StringIO()
    p.set_target(tgt)
    p.write("routed\n")                          # this thread -> its target
    assert tgt.getvalue() == "routed\n"
    p.clear_target()
    p.write("back-to-default\n")
    assert default.getvalue().endswith("back-to-default\n")


def test_targets_are_isolated_across_threads():
    default = io.StringIO()
    p = _ThreadLocalStdout(default)
    seen = {}

    def worker(name):
        buf = io.StringIO()
        p.set_target(buf)
        p.write(f"{name}\n")
        p.write(f"{name}-again\n")
        seen[name] = buf.getvalue()
        p.clear_target()

    ts = [threading.Thread(target=worker, args=(f"t{i}",)) for i in range(4)]
    for t in ts: t.start()
    for t in ts: t.join()
    # each thread only ever saw its OWN writes -- no cross-thread clobbering
    for name, val in seen.items():
        assert val == f"{name}\n{name}-again\n"


def test_half_built_proxy_does_not_raise_local(monkeypatch=None):
    # Simulate a proxy reconstructed WITHOUT __init__ (deepcopy/pickle) or torn down:
    # __dict__ has neither _default nor _local. write()/flush() must fall back to the
    # real interpreter stdout, NOT raise AttributeError('_local').
    broken = _ThreadLocalStdout.__new__(_ThreadLocalStdout)   # __init__ NOT called
    assert "_local" not in broken.__dict__ and "_default" not in broken.__dict__

    captured = io.StringIO()
    real = sys.__stdout__
    sys.__stdout__ = captured                      # where the graceful fallback should land
    try:
        broken.write("survived\n")                 # must not raise
        broken.flush()                             # must not raise
    finally:
        sys.__stdout__ = real
    assert captured.getvalue() == "survived\n"


def test_delegates_public_attrs_but_guards_private():
    p = _ThreadLocalStdout(sys.__stdout__)
    # public attribute delegation still works (isatty/encoding/etc.)
    assert hasattr(p, "write")
    # underscore misses still raise AttributeError (so deepcopy/pickle probing
    # __deepcopy__/__getstate__ sees "not implemented"), but the message is the name
    raised = False
    try:
        p.__dict__  # fine
        object.__getattribute__(p, "totally_missing_dunder_placeholder")
    except AttributeError:
        pass
    try:
        p._nonexistent_private
    except AttributeError as e:
        raised = True
        assert str(e) == "_nonexistent_private"
    assert raised


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} STDOUT-PROXY TESTS PASSED")
