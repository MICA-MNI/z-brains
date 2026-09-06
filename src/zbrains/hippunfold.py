"""Small helpers for version-safe HippUnfold cache identities."""

import hashlib
import os
import re


_VERSION_RE = re.compile(r"(?:^|[_-])v(?P<major>\d+)(?:[._-]|$)", re.IGNORECASE)
_CACHE_TAG_RE = re.compile(r"_hippunfoldv(?:\d+|unknown)-[0-9a-f]{10}")


def infer_hippunfold_major(hippunfold_directory):
    """Infer the major version from a derivative directory name when possible."""
    if not hippunfold_directory:
        return None
    match = _VERSION_RE.search(os.path.basename(os.path.normpath(os.fspath(hippunfold_directory))))
    return int(match.group("major")) if match else None


def hippunfold_cache_tag(hippunfold_directory, version=None):
    """Return a stable cache tag tied to both major version and source directory."""
    if not hippunfold_directory:
        return ""
    if version is None:
        version = infer_hippunfold_major(hippunfold_directory)
    if version is not None:
        version = int(version)
        if version not in (1, 2):
            raise ValueError(f"Unsupported HippUnfold major version: {version}")
    source = os.path.realpath(os.path.abspath(os.fspath(hippunfold_directory)))
    source_hash = hashlib.sha1(source.encode("utf-8")).hexdigest()[:10]
    major = str(version) if version is not None else "unknown"
    return f"hippunfoldv{major}-{source_hash}"


def without_hippunfold_cache_tag(path):
    """Return the matching legacy output path with its HippUnfold tag removed."""
    return _CACHE_TAG_RE.sub("", os.fspath(path))
