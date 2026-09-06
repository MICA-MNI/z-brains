import os
import sys
import datetime
import copy
from zbrains.hippunfold import hippunfold_cache_tag
from zbrains.processing import apply_blurring, apply_hippocampal_processing, apply_subcortical_processing, apply_cortical_processing, generate_superficial_white_matter
from zbrains.normalization import (
    apply_nyul_model_to_subject,
    apply_ravel_model_to_subject,
    compose_normalization_label,
    decompose_normalization_label,
    ensure_synthseg_csf,
    fit_and_apply_nyul_to_controls,
    fit_and_apply_ravel_to_controls,
    load_ravel_models,
    normalize_normalization_mode,
    prepare_t1w_flair_identity,
    prepare_t1w_flair_whitestripe,
    prepare_t1w_flair_wmmean,
    requested_ravel_modalities,
    resolve_normalization_desc,
)
from zbrains.analysis import analyze_dataset
from zbrains.clinical_reports import generate_clinical_report
import shutil
import re
import threading
import multiprocessing
from joblib import Parallel, delayed
import nibabel as nib
import numpy as np

class LogRedirect:
    """Redirect stdout to a per-subject log file, optionally mirroring to terminal.

    With ``mirror_terminal=False`` it writes ONLY to the log file -- used by the
    parallel subject loop, where mirroring N concurrent subjects onto one terminal
    would interleave into noise (each subject's complete log still lands on disk).
    """
    def __init__(self, log_file, mirror_terminal=True):
        self.log_file = log_file
        # Capture the REAL stdout, unwrapping the thread-local proxy if one is
        # installed, so mirroring can never recurse back through the proxy into
        # this same object.
        self.terminal_stdout = getattr(sys.stdout, "_default", sys.stdout)
        self.mirror_terminal = mirror_terminal

        # Remove old log file if it exists
        if os.path.exists(log_file):
            os.remove(log_file)

        # Open in write mode to create a new file
        self.log = open(self.log_file, "w", encoding="utf-8")

    def __del__(self):
        self.close()

    def close(self):
        if getattr(self, "log", None):
            try:
                self.log.close()
            finally:
                self.log = None

    def write(self, message):
        if self.mirror_terminal:
            self.terminal_stdout.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        if self.mirror_terminal:
            self.terminal_stdout.flush()
        self.log.flush()


class _ThreadLocalStdout:
    """A ``sys.stdout`` stand-in that dispatches writes to a PER-THREAD target.

    Installed once around the parallel subject loop so each worker thread logs to
    its own file via ``set_target`` -- WITHOUT the shared-global ``sys.stdout``
    swap that made concurrent subjects clobber each other's redirect (the reason
    the parallel path was originally disabled). A thread with no target set (e.g.
    the orchestrating main thread) falls back to the real stdout.
    """
    def __init__(self, default):
        self._default = default
        self._local = threading.local()

    def _get_local(self):
        # Read via __dict__ (never __getattr__, which raises for underscore names)
        # and re-create lazily if absent, so a proxy reconstructed WITHOUT __init__
        # (deepcopy/pickle) or read mid-teardown degrades gracefully instead of
        # raising AttributeError('_local') out of a caller's print().
        local = self.__dict__.get("_local")
        if local is None:
            local = threading.local()
            self.__dict__["_local"] = local
        return local

    def _current(self):
        # A logging proxy must NEVER be able to crash the compute: fall back all the
        # way to the real interpreter stdout rather than raising if internal state
        # is missing (half-built instance, teardown, unexpected access order).
        target = getattr(self._get_local(), "target", None)
        return target or self.__dict__.get("_default") or sys.__stdout__

    def set_target(self, stream):
        self._get_local().target = stream

    def clear_target(self):
        self._get_local().target = None

    def write(self, message):
        return self._current().write(message)

    def flush(self):
        cur = self._current()
        flush = getattr(cur, "flush", None)
        if callable(flush):
            flush()

    def __getattr__(self, name):
        # Delegate isatty()/fileno()/encoding/etc. to the current target; guard
        # private names so an early/attribute miss can't recurse via _current().
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self._current(), name)


# Per-subject processing runs each subject in a worker THREAD: the heavy work is
# external wb_command/ANTs subprocesses that release the GIL, so threads give real
# speedup with bit-identical output. A process pool is deliberately avoided --
# loky/fork + this numpy build double-inits ("CPU dispatcher tracer already
# initlized") and SIGSEGVs. Capped because each concurrent subject holds a
# subject's volumes in RAM; raise/lower via env.num_threads if you have the memory.
PROCESSING_MAX_THREADS = 12
NORMALIZATION_MAX_JOBS = 4


class demographics():
    def __init__(self, csv_file, column_mapping=None, normative_columns=None, normative_dtypes=None, reference=None, subset=None):
        # Resolve to an ABSOLUTE path at construction time (cwd is correct here).
        # Later steps of the pipeline change the working directory, and the
        # demographics may be re-read then (e.g. when building a control subset for
        # held-out-control K-fold scoring); a stored relative path would fail.
        self.csv_file = os.path.abspath(csv_file) if isinstance(csv_file, (str, os.PathLike)) else csv_file
        self.data = None
        self.column_mapping = {
            "ID": "participant_id",
            "SES": "session_id",
        } if column_mapping is None else column_mapping
        self.normative_columns = normative_columns
        self.normative_dtypes = normative_dtypes
        self.reference = reference
        self.subset = subset
        self.binary_encodings = {}  # Store binary variable encodings
        self.load_data()

    def __repr__(self):
        return f"Demographics(csv_file={self.csv_file})"
    def __str__(self):
        return f"Demographics data loaded from {self.csv_file}"
    
    def load_data(self):
        import pandas as pd
        import numpy as np
        try:
            self.data = pd.read_csv(self.csv_file)
        except Exception as e:
            raise ValueError(f"Failed to load demographics data: {e}")
        if self.data.empty:
            raise ValueError("Demographics data is empty.")
        # Rename columns based on mapping
        self.data.rename(columns=self.column_mapping, inplace=True)
        print(f"Demographics data loaded from {self.csv_file} with {len(self.data)} entries.")
        if not set(self.column_mapping.values()).issubset(self.data.columns):
            raise ValueError("Demographics data does not contain all required columns.")
        
        if self.subset is not None:
            subset_pairs = []
            for item in self.subset:
                if not isinstance(item, (list, tuple)) or len(item) != 2:
                    raise ValueError(f"Subset entries must be 2-element lists/tuples. Invalid entry: {item}")
                subset_pairs.append(tuple(item))
            subset_set = set(subset_pairs)
            available_pairs = set(zip(self.data['participant_id'], self.data['session_id']))
            missing_pairs = [pair for pair in subset_pairs if pair not in available_pairs]
            if missing_pairs:
                raise ValueError(f"Subset entries not found in demographics data: {missing_pairs}")
            mask = self.data[['participant_id', 'session_id']].apply(tuple, axis=1).isin(subset_set)
            self.data = self.data[mask].reset_index(drop=True)
            if self.data.empty:
                raise ValueError("Subset filtering resulted in empty demographics data.")
            print(f"Filtered demographics to subset with {len(self.data)} entries.")
        
        # Use reference normative columns if available and none are specified
        single_subject = len(self.data) == 1
        using_reference = False
        
        if self.reference is not None:
            if self.normative_columns is None and hasattr(self.reference, 'normative_columns'):
                self.normative_columns = self.reference.normative_columns
                self.normative_dtypes = self.reference.normative_dtypes
                using_reference = True
                print(f"Using normative columns from reference: {self.normative_columns}")
        
        # Validate normative columns and dtypes
        normative_columns = self.normative_columns if self.normative_columns else []
        
        # Check if normative_dtypes is provided when normative_columns is specified
        if normative_columns and self.normative_dtypes is None:
            raise ValueError("normative_dtypes must be provided when normative_columns is specified")
            
        # Check if normative_dtypes has the same length as normative_columns
        if normative_columns and len(normative_columns) != len(self.normative_dtypes):
            raise ValueError(f"normative_dtypes ({len(self.normative_dtypes)}) must have the same length as normative_columns ({len(normative_columns)})")
            
        # Validate each normative column
        for col in normative_columns:
            if col not in self.data.columns:
                raise ValueError(f"Normative column '{col}' not found in demographics data.")
            
            # Check for missing values
            if self.data[col].isnull().any():
                raise ValueError(f"Normative column '{col}' contains missing values.")
        
        # Validate data types if both normative_columns and normative_dtypes are specified
        if normative_columns and self.normative_dtypes:
            for i, (col, dtype) in enumerate(zip(normative_columns, self.normative_dtypes)):
                try:
                    # Validate column values based on specified data type
                    if dtype.lower() == 'int':
                        # Coerce numeric values and round them to integer years.
                        # Many BIDS/site exports store age as decimal years
                        # (e.g. 32.7); downstream normative models expect an
                        # integer column when dtype="int".
                        numeric_values = pd.to_numeric(self.data[col], errors='raise')
                        rounded_values = np.rint(numeric_values).astype(int)
                        changed = ~np.isclose(numeric_values, rounded_values)
                        if np.any(changed):
                            changed_values = [
                                float(v) for v in sorted(pd.Series(numeric_values[changed]).unique())
                            ]
                            print(
                                f"Rounded non-integer values in column '{col}' "
                                f"to nearest integer: {changed_values}"
                            )
                        self.data[col] = rounded_values
                    
                    elif dtype.lower() == 'float':
                        # Try converting to float
                        self.data[col] = self.data[col].astype(float)
                    
                    elif dtype.lower() == 'binary':
                        # Check if we should use reference encoding
                        if self.reference is not None and hasattr(self.reference, 'binary_encodings') and col in self.reference.binary_encodings:
                            # Get encoding from reference
                            ref_value_to_binary = self.reference.binary_encodings[col]
                            
                            # Convert using reference encoding
                            self.data[col] = self.data[col].astype(str).str.lower().map(
                                lambda x: ref_value_to_binary.get(x.lower(), 0)
                            )
                            print(f"Applied reference encoding for binary column '{col}'")
                            
                        else:
                            # Check if the column has valid binary values (1 or 2 unique values)
                            unique_values = self.data[col].astype(str).str.lower().unique()
                            # Sort to ensure deterministic encoding (0 for first alphabetical, 1 for second)
                            unique_values.sort()
                            
                            if len(unique_values) > 2:
                                raise ValueError(f"Column '{col}' has {len(unique_values)} unique values, expected 1 or 2 for binary data type")
                            
                            # Convert to numeric binary values (0/1)
                            # Using sorted values: 1st alphabetical -> 0, 2nd -> 1
                            # This handles len=1 (val -> 0) and len=2 (val1->0, val2->1)
                            value_to_binary = {val: i for i, val in enumerate(unique_values)}
                            
                            # Store the encoding for future reference
                            self.binary_encodings[col] = value_to_binary
            
                            # Apply encoding
                            self.data[col] = self.data[col].astype(str).str.lower().map(value_to_binary)
                    
                    else:
                        raise ValueError(f"Unsupported dtype '{dtype}' for column '{col}'. "
                                        f"Supported types: 'int', 'float', 'binary'")
                        
                except Exception as e:
                    raise ValueError(f"Error validating column '{col}' with dtype '{dtype}': {e}")
            
            print(f"Normative columns: {normative_columns} with types {self.normative_dtypes} successfully validated.")
            if using_reference and single_subject:
                print("Used reference encoding for binary variables to ensure consistency in normative modeling.")
        else:
            print(f"Normative columns: {normative_columns} successfully validated.")
    
        return self


def _link_structural_to_cache(
    output_directory,
    participant_id,
    session_id,
    verbose=False,
    hippunfold_directory=None,
    hippunfold_version=None,
):
    """Point this base's ``structural/`` dir at a base-independent per-subject cache.

    The structural outputs are pure geometry, but the hippocampal files depend on
    the exact HippUnfold source. They are therefore shared across optimization arms
    only within a version/source-specific cache. This prevents v1 0p5mm geometry
    from being exposed to a v2 8k run. No-op when ``output_directory`` has no
    ``zbrains_*`` component (legacy) or a REAL structural dir already exists here.
    Must be called BEFORE any structural files are written.
    """
    parts = os.path.abspath(output_directory).split(os.sep)
    idx = next((i for i, p in enumerate(parts)
                if p.startswith("zbrains_") or p.startswith("zbrains-")), None)
    if idx is None:
        return None
    cache_tag = hippunfold_cache_tag(
        hippunfold_directory,
        version=hippunfold_version,
    )
    parts[idx] = "struct_cache" + (f"_{cache_tag}" if cache_tag else "")
    cache_struct = os.path.join(os.sep.join(parts), participant_id, session_id, "structural")
    base_struct = os.path.join(output_directory, participant_id, session_id, "structural")
    if os.path.isdir(base_struct) and not os.path.islink(base_struct):
        return base_struct                 # real dir already here -> leave it
    os.makedirs(cache_struct, exist_ok=True)
    os.makedirs(os.path.dirname(base_struct), exist_ok=True)
    cache_abs = os.path.abspath(cache_struct)
    if os.path.islink(base_struct):
        if os.path.realpath(base_struct) == cache_abs:
            return base_struct             # already linked correctly
        try:
            os.remove(base_struct)
        except OSError:
            pass
    try:
        os.symlink(cache_abs, base_struct)
        if verbose:
            print(f"  Structural cache: {participant_id}/{session_id} -> {cache_struct}")
    except FileExistsError:
        pass
    return base_struct


# --- Norm-phase volume cache (smoothing-independent) --------------------------
# WhiteStripe/RAVEL/Nyul produce nativepro ``desc-<final>_{T1w,FLAIR}.nii.gz`` volumes
# (and, for RAVEL/Nyul, a dataset-level fit model) BEFORE any volume-to-surface
# sampling; on-surface smoothing is applied later, per-vertex. So those volumes are
# byte-identical across every cortical/hippocampal smoothing (and sampling) arm. The
# processed base directory is keyed by smoothing, so each smoothing arm otherwise
# re-runs the whole (usually RAVEL-dominated) normalization. These helpers reuse the
# volumes from a sibling cache keyed by everything EXCEPT smoothing.
_SMOOTH_TOKEN_RE = re.compile(r"_smoothctx\d+hip\d+")
_HIPPUNFOLD_TOKEN_RE = re.compile(
    r"_hippunfoldv(?:\d+|unknown)-[0-9a-f]{10}"
)
_NORM_MODEL_DIRS = ("ravel", "nyul")


def _norm_cache_dir_for(output_directory):
    """Sibling norm-cache base dir for a processed base, keyed by everything EXCEPT
    on-surface smoothing. Derived from the base path so the per-fold ``_reffold{k}``
    tag and any exclusion tag are PRESERVED -- a cross-validation fold's train-only
    RAVEL/Nyul fit is never shared with the full-controls run or another fold (no
    leakage). Returns None for a non-``zbrains_`` path (legacy/direct)."""
    parts = os.path.abspath(output_directory).split(os.sep)
    idx = next((i for i, p in enumerate(parts)
                if p.startswith("zbrains_") or p.startswith("zbrains-")), None)
    if idx is None:
        return None
    comp = parts[idx]
    sep = "_" if comp.startswith("zbrains_") else "-"
    comp = "norm_cache" + sep + comp.split("zbrains" + sep, 1)[1]
    comp = _SMOOTH_TOKEN_RE.sub("", comp)                 # smoothing arms share ONE cache
    # HippUnfold changes only hippocampal geometry/sampling. Nativepro intensity
    # volumes and their RAVEL/Nyul models are safe to share across HU versions.
    comp = _HIPPUNFOLD_TOKEN_RE.sub("", comp)
    parts[idx] = comp
    return os.sep.join(parts[:idx + 1])


def _atomic_symlink(src, dst):
    """Point ``dst`` at ``src`` atomically (temp symlink + os.replace) so a peer on
    shared storage never sees a half-made link. Leaves a REAL file/dir at ``dst``
    untouched (only creates when absent or already the same link)."""
    if os.path.lexists(dst):
        if os.path.islink(dst) and os.path.realpath(dst) == os.path.realpath(src):
            return True
        return os.path.exists(dst)                        # something real is here -> keep it
    tmp = f"{dst}.link.{os.getpid()}.{threading.get_ident()}"
    try:
        os.symlink(os.path.abspath(src), tmp)
        os.replace(tmp, dst)
    except FileExistsError:
        pass
    except OSError:
        try:
            os.remove(tmp)
        except OSError:
            pass
        return False
    return True


def _copy_atomic(src, dst):
    """Copy ``src`` -> ``dst`` atomically, skipping if ``dst`` already exists (so
    concurrent arms/machines populating one shared cache never clobber each other)."""
    if os.path.exists(dst):
        return True
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = f"{dst}.tmp.{os.getpid()}.{threading.get_ident()}"
    try:
        shutil.copy2(src, tmp)
        os.replace(tmp, dst)
    except OSError:
        try:
            os.remove(tmp)
        except OSError:
            pass
        return False
    return True


def _mirror_model_dirs(src_root, dst_root):
    """Copy any RAVEL/Nyul fit-model files from ``src_root/<model>/`` into
    ``dst_root/<model>/`` (real files only, skip existing). Small .npz/.json; used
    both to save a fit into the cache and to restore it into a base."""
    for name in _NORM_MODEL_DIRS:
        sdir = os.path.join(src_root, name)
        if not os.path.isdir(sdir):
            continue
        try:
            entries = os.listdir(sdir)
        except OSError:
            continue
        for fn in entries:
            sp = os.path.join(sdir, fn)
            if os.path.isfile(sp) and not os.path.islink(sp):
                _copy_atomic(sp, os.path.join(dst_root, name, fn))


def _norm_volume_filename(participant_id, session_id, desc, modality):
    return f"{participant_id}_{session_id}_space-nativepro_desc-{desc}_{modality}.nii.gz"


def _norm_cache_hydrate(cache_dir, base_dir, subjects_modalities, desc, *, verbose=False):
    """Symlink cached desc-<final> volumes into ``base_dir`` maps/ (ALL-OR-NOTHING)
    and restore any cached fit model. Returns True iff EVERY (subject, modality) was
    present in the cache -- a partial cache hydrates nothing (beyond the model) and
    returns False so the norm phase recomputes fully (never leaving some controls
    without the desc-<subject_norm> the dataset fit needs)."""
    if os.path.isdir(cache_dir):
        _mirror_model_dirs(cache_dir, base_dir)          # model first (patient apply may need it)
    pending = []
    for (participant_id, session_id), mods in subjects_modalities.items():
        for modality in mods:
            fn = _norm_volume_filename(participant_id, session_id, desc, modality)
            src = os.path.join(cache_dir, participant_id, session_id, "maps", fn)
            if not os.path.exists(src):
                return False                             # incomplete -> full recompute
            dst = os.path.join(base_dir, participant_id, session_id, "maps", fn)
            pending.append((src, dst))
    if not pending:
        return False                                    # nothing cacheable to reuse
    for src, dst in pending:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if not _atomic_symlink(src, dst):
            return False
    return True


def _norm_cache_save(cache_dir, base_dir, subjects_modalities, desc, *, verbose=False):
    """Persist freshly-computed desc-<final> volumes (+ fit model) into the norm cache
    so sibling smoothing/sampling arms reuse them. Skips symlinks (already cached) and
    missing/failed subjects; atomic + skip-existing so it is safe under concurrency."""
    if not cache_dir:
        return
    os.makedirs(cache_dir, exist_ok=True)
    _mirror_model_dirs(base_dir, cache_dir)
    n = 0
    for (participant_id, session_id), mods in subjects_modalities.items():
        for modality in sorted(mods):
            fn = _norm_volume_filename(participant_id, session_id, desc, modality)
            src = os.path.join(base_dir, participant_id, session_id, "maps", fn)
            if not os.path.isfile(src) or os.path.islink(src):
                continue                                # failed/missing, or already a cache symlink
            if _copy_atomic(src, os.path.join(cache_dir, participant_id, session_id, "maps", fn)):
                n += 1
    if verbose and n:
        print(f"Cached {n} desc-{desc} volume(s) to {cache_dir} for smoothing/sampling reuse.")


def _cortex_bool_mask(cortex_idx, n_vertices, *, hemi="", label=""):
    """Build a length-``n_vertices`` cortex boolean mask from FreeSurfer label indices.

    The medial-wall mask must stay aligned to the fsnative surface (``lh/rh.white``),
    so it is always sized to that surface. FreeSurfer's ``?h.cortex.label`` and
    ``?h.white`` should share a vertex count, but for some subjects the label carries
    an index at/after the surface's last vertex (a 1-past-the-end off-by-one, or a
    surface/label recon mismatch). Rather than crash the whole subject on the
    out-of-bounds fancy-index assignment, drop the offending indices (they name
    vertices that don't exist on this surface anyway) and warn so a gross mismatch
    stays visible instead of silently corrupting the mask.
    """
    idx = np.asarray(cortex_idx, dtype=np.int64).ravel()
    mask = np.zeros(int(n_vertices), dtype=bool)
    if idx.size == 0:
        return mask
    in_range = (idx >= 0) & (idx < n_vertices)
    n_bad = int((~in_range).sum())
    if n_bad:
        print(
            f"  Warning: {label} hemi-{hemi} cortex label has {n_bad} vertex "
            f"index(es) outside the fsnative surface (n={n_vertices}, "
            f"max idx {int(idx.max())}); dropping them from the cortex mask."
        )
    mask[idx[in_range]] = True
    return mask


class zbdataset():
    def __init__(self, name, demographics : demographics, micapipe_directory, hippunfold_directory=None, freesurfer_directory=None, raw_data_directory=None, cortex=True, hippocampus=True, subcortical=True, hippunfold_version=None):
        
        self.name = name
        self.demographics = demographics
        self.micapipe_directory = micapipe_directory
        self.hippunfold_directory = hippunfold_directory
        self.freesurfer_directory = freesurfer_directory
        self.raw_data_directory = raw_data_directory
        self.features = []
        self.cortex = cortex
        self.hippocampus = hippocampus
        self.subcortical = subcortical
        self.valid_dataset = False
        if hippunfold_version is not None and int(hippunfold_version) not in (1, 2):
            raise ValueError("hippunfold_version must be 1, 2, or None")
        self.requested_hippunfold_version = (
            int(hippunfold_version) if hippunfold_version is not None else None
        )
        self.hippunfold_version = self.requested_hippunfold_version or 1
        self._hippunfold_version_detected = False
        # Optional feature-specific control QC mask.  The staged correlation arms
        # populate this as ``{feature: {(participant_id, session_id), ...}}``.
        # A masked scan remains a member of the dataset and stays available for
        # every other feature; it is removed only from the matching feature (and
        # that feature's *blur companion, which shares the same source volume).
        self.control_feature_exclusions = {}
       
    def __repr__(self):
        return f"Dataset(name={self.name})"

    def __str__(self):
        return f"Dataset: {self.name}"

    @staticmethod
    def _feature_exclusion_key(feature):
        """Canonical base-feature key used by feature-specific control masks."""
        key = str(feature).strip().lower().replace("*blur", "")
        return "qt1" if key == "t1map" else key

    def _apply_control_feature_exclusions(self, valid_subjects=None):
        """Prune only masked subject-feature observations from a validity table.

        ``base`` and the dataset-wide structure lists are deliberately untouched:
        a control excluded for FLAIR, for example, must still contribute T1w, FA,
        and every other available feature.  A base feature and its ``*blur``
        companion share one mask because both are extracted from the same volume
        and therefore cannot use different RAVEL/Nyul training cohorts.
        """
        target = self.valid_subjects if valid_subjects is None else valid_subjects
        exclusions = getattr(self, "control_feature_exclusions", None) or {}
        if not exclusions or not target:
            return target

        by_feature = {}
        for feature, pairs in exclusions.items():
            key = self._feature_exclusion_key(feature)
            by_feature.setdefault(key, set()).update(
                (str(pid), str(ses)) for pid, ses in (pairs or ())
            )

        for feature in list(getattr(self, "features", ())) or list(target):
            data = target.get(feature)
            if not isinstance(data, dict) or "all" not in data:
                continue
            excluded = by_feature.get(self._feature_exclusion_key(feature), set())
            if not excluded:
                continue
            data["all"] = [tuple(pair) for pair in data.get("all", ())
                           if tuple(map(str, pair)) not in excluded]
            for structure, pairs in data.get("structures", {}).items():
                data["structures"][structure] = [
                    tuple(pair) for pair in pairs
                    if tuple(map(str, pair)) not in excluded
                ]
        return target
    
    def check_directories(self):
        if self.hippunfold_directory is None:
            print("Warning: hippunfold_directory is not specified, hippocampal data will not be available.")
        if self.freesurfer_directory is None:
            print("Warning: freesurfer_directory is not specified, subcortical data will not be available.")
        if not getattr(self, 'raw_data_directory', None):
            print("Warning: raw_data_directory is not specified, T1w/FLAIR WhiteStripe normalization will use micapipe nativepro volumes.")
        print("Checking agreement between demographics data and micapipe directory...")

        if not os.path.exists(self.micapipe_directory):
            raise ValueError(f"Micapipe directory {self.micapipe_directory} does not exist.")
        if getattr(self, 'raw_data_directory', None) and not os.path.exists(self.raw_data_directory):
            raise ValueError(f"Raw data directory {self.raw_data_directory} does not exist.")
        if self.hippunfold_directory and not os.path.exists(self.hippunfold_directory):
            raise ValueError(f"Hippunfold directory {self.hippunfold_directory} does not exist.")
        if self.freesurfer_directory and not os.path.exists(self.freesurfer_directory):
            raise ValueError(f"Freesurfer directory {self.freesurfer_directory} does not exist.")
        print("All directories exist.")

        # Detect HippUnfold version if applicable
        if self.hippocampus and self.hippunfold_directory:
            # Probe a few directories to find the resolution
            found_version = False
            for _, row in self.demographics.data.iterrows():
                pid = row['participant_id']
                sid = row['session_id']

                if self.requested_hippunfold_version == 2:
                    v2_path = os.path.join(
                        self.hippunfold_directory,
                        pid,
                        sid,
                        f"surf/{pid}_{sid}_hemi-L_space-T1w_den-8k_label-hipp_midthickness.surf.gii",
                    )
                    if os.path.exists(v2_path):
                        found_version = True
                        break
                    continue

                if self.requested_hippunfold_version == 1:
                    v1_path = os.path.join(
                        self.hippunfold_directory,
                        "hippunfold",
                        pid,
                        sid,
                        f"surf/{pid}_{sid}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
                    )
                    if os.path.exists(v1_path):
                        found_version = True
                        break
                    continue
                
                # Check for 0.5mm file (v1 structure: hippunfold_dir/hippunfold/sub-X)
                v1_path = os.path.join(self.hippunfold_directory, "hippunfold", pid, sid, 
                                      f"surf/{pid}_{sid}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii")
                if os.path.exists(v1_path):
                    self.hippunfold_version = 1
                    self._hippunfold_version_detected = True
                    print("Detected HippUnfold resolution: 0.5mm (v1)")
                    found_version = True
                    break
                
                # Check for 8k file (v2 structure: hippunfold_dir/sub-X)
                v2_path = os.path.join(self.hippunfold_directory, pid, sid, 
                                      f"surf/{pid}_{sid}_hemi-L_space-T1w_den-8k_label-hipp_midthickness.surf.gii")
                if os.path.exists(v2_path):
                    self.hippunfold_version = 2
                    self._hippunfold_version_detected = True
                    print("Detected HippUnfold resolution: 8k (v2)")
                    found_version = True
                    break
            
            if self.requested_hippunfold_version is not None:
                self.hippunfold_version = self.requested_hippunfold_version
                self._hippunfold_version_detected = True
                if found_version:
                    resolution = "8k" if self.hippunfold_version == 2 else "0.5mm"
                    print(
                        f"Using configured HippUnfold v{self.hippunfold_version} "
                        f"({resolution}); automatic fallback is disabled."
                    )
                else:
                    print(
                        f"Warning: configured HippUnfold v{self.hippunfold_version}, "
                        "but no matching surface was found while probing subjects. "
                        "Keeping the configured version; v1/v2 fallback is disabled."
                    )
            elif not found_version:
                print("Warning: Could not detect HippUnfold version (checked 8k/v2 and 0.5mm/v1). Defaulting to v1.")

        # Initialize valid_subjects structure
        self.valid_subjects = {
            'base': [],
            'structures': {
                'cortex': [],
                'hippocampus': [],
                'subcortical': []
            }
        }
        
        # Check if demographics data matches micapipe directory
        missing_micapipe = []
        missing_hippunfold = []
        missing_freesurfer = []
        
        for _, row in self.demographics.data.iterrows():
            participant_id = row['participant_id']
            session_id = row['session_id']
            
            # Track if this participant/session has all required directories
            is_valid = True
            
            # Check for directory in micapipe
            participant_session_path = os.path.join(self.micapipe_directory, f"{participant_id}", f"{session_id}")
            if not os.path.exists(participant_session_path):
                missing_micapipe.append((participant_id, session_id))
                is_valid = False
                
            # Check for directory in hippunfold if specified
            if self.hippocampus and self.hippunfold_directory:
                if self.hippunfold_version == 2:
                    hippunfold_path = os.path.join(self.hippunfold_directory, f"{participant_id}", f"{session_id}")
                else:
                    hippunfold_path = os.path.join(self.hippunfold_directory, "hippunfold", f"{participant_id}", f"{session_id}")
                
                if not os.path.exists(hippunfold_path):
                    missing_hippunfold.append((participant_id, session_id))
                    is_valid = False
                    
            # Check for directory in freesurfer if specified
            if self.subcortical and self.freesurfer_directory:
                freesurfer_path = os.path.join(self.freesurfer_directory, f"{participant_id}_{session_id}")
                if not os.path.exists(freesurfer_path):
                    missing_freesurfer.append((participant_id, session_id))
                    is_valid = False
            
            # If all required directories exist, add to valid subjects
            if is_valid:
                self.valid_subjects['base'].append((participant_id, session_id))
        
        # Report missing directories
        if missing_micapipe:
            print(f"Warning: {len(missing_micapipe)} participant/session directories not found in micapipe directory:")
            for p, s in missing_micapipe:
                print(f"  - {p}/{s}")
                
        if missing_hippunfold:
            print(f"Warning: {len(missing_hippunfold)} participant/session directories not found in hippunfold directory:")
            for p, s in missing_hippunfold:
                print(f"  - {p}/{s}")
                
        if missing_freesurfer:
            print(f"Warning: {len(missing_freesurfer)} participant/session directories not found in freesurfer directory:")
            for p, s in missing_freesurfer:
                print(f"  - {p}/{s}")
        
        print(f"Found {len(self.valid_subjects['base'])} valid subjects with complete directory structure.")
        return self
    
    def add_features(self, *features, verbose=True):
        """
        Add features to the dataset.
        Features should be specified as strings, e.g., "FA", "ADC", etc.
        
        This function checks for the existence of required files in the micapipe directory
        for each feature and each valid subject.
        
        Parameters:
        -----------
        *features : str
            Names of features to add
        verbose : bool, default=True
            If True, prints detailed information about missing files and features.
            If False, only prints the summary with percentages.
        """
        
        self.check_directories()
        
        # Standardize feature names consistently - apply mapping first
        feature_mapping = {
            "fa": "FA",
            "adc": "ADC", 
            "thickness": "thickness",
            "sa": "SA",
            "flair": "FLAIR",
            "qt1": "qT1",
            "qt1*blur": "qT1*blur",
            "flair*blur": "FLAIR*blur",
            "fmri": "fMRI",
            "t1w": "T1w",
            "t1w*blur": "T1w*blur",
            "adc*blur": "ADC*blur",
            "fa*blur": "FA*blur",
        }
        
        # Normalize all features to consistent case
        normalized_features = []
        for feature in features:
            # Convert to lowercase for mapping lookup
            feature_lower = str(feature).lower()
            
            # Apply mapping to get standardized case
            if feature_lower in feature_mapping:
                normalized_features.append(feature_mapping[feature_lower])
            else:
                # For unknown features, keep original case but warn user
                normalized_features.append(str(feature))
                if verbose:
                    print(f"Warning: Unknown feature '{feature}' - using as provided")
        
        # Store the normalized features
        self.features = normalized_features
        
        if verbose:
            print(f"Features {self.features} added to the dataset {self.name}.")
        
        # Define file patterns for each feature (using standardized names)
        feature_files = {
            "FA": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_FA.func.gii",
                   "maps/{participant_id}_{session_id}_space-nativepro_model-DTI_map-FA.nii.gz"],
            "ADC": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_ADC.func.gii", 
                    "maps/{participant_id}_{session_id}_space-nativepro_model-DTI_map-ADC.nii.gz"],
            "thickness": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-thickness.func.gii"],
            "FLAIR": ["maps/{participant_id}_{session_id}_space-nativepro_map-flair.nii.gz"],
            "qT1": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_T1map.func.gii", 
                    "maps/{participant_id}_{session_id}_space-nativepro_map-T1map.nii.gz"],
            
            # Fix: Remove surface file requirement for blur features (inputs only)
            "FLAIR*blur": ["maps/{participant_id}_{session_id}_space-nativepro_map-flair.nii.gz"],
            "qT1*blur": ["maps/{participant_id}_{session_id}_space-nativepro_map-T1map.nii.gz"],
            "ADC*blur": ["maps/{participant_id}_{session_id}_space-nativepro_model-DTI_map-ADC.nii.gz"],
            "FA*blur": ["maps/{participant_id}_{session_id}_space-nativepro_model-DTI_map-FA.nii.gz"],

            "fMRI": ["func/desc-se_task-rest_acq-AP_bold/surf/{participant_id}_{session_id}_surf-fsLR-32k_desc-timeseries_clean.shape.gii"],
            "SA": ["surf/{participant_id}_{session_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-{surfacetype}.surf.gii"],
            "T1w": ["anat/{participant_id}_{session_id}_space-nativepro_T1w.nii.gz"],
            "T1w*blur": ["anat/{participant_id}_{session_id}_space-nativepro_T1w.nii.gz"],
        }
        
        # Rest of the method remains the same...
        # Define required files for different structures
        required_cortical_files = [
            "anat/{participant_id}_{session_id}_space-nativepro_T1w.nii.gz",
            # "anat/{participant_id}_{session_id}_space-nativepro_T1w_brain_mask.nii.gz",
            "surf/{participant_id}_{session_id}_hemi-L_surf-fsnative_label-sphere.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_surf-fsnative_label-sphere.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-white.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-white.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_surf-fsnative_label-sphere.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_surf-fsnative_label-sphere.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-pial.surf.gii",
        ]
        
        required_hippocampal_files = [
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_outer.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_outer.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii"
        ]

        required_hippocampal_files_v2 = [
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-8k_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-8k_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-8k_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-8k_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-8k_label-hipp_outer.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-8k_label-hipp_outer.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-unfold_den-8k_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-unfold_den-8k_label-hipp_midthickness.surf.gii"
        ]

        required_freesurfer_files = [
            "{participant_id}_{session_id}/surf/lh.white",
            "{participant_id}_{session_id}/surf/rh.white",
            "{participant_id}_{session_id}/label/lh.cortex.label",
            "{participant_id}_{session_id}/label/rh.cortex.label"
        ]

        required_subcortical_files = [
            "{participant_id}_{session_id}/stats/aseg.stats",
        ]

        # Track missing files and feature availability
        missing_files = {}
        missing_features = {}
        feature_availability = {feature: 0 for feature in self.features}
        subjects_with_complete_data = []
        
        # Initialize feature-specific valid subjects with nested structure
        for feature in self.features:
            self.valid_subjects[feature] = {
                'all': [],  # All structures
                'structures': {
                    'cortex': [],
                    'hippocampus': [],
                    'subcortical': []
                }
            }
        
        # Check for each valid subject from the base list
        for participant_id, session_id in self.valid_subjects['base']:
            subject_missing_files = []
            subject_missing_features = set()
            has_required_cortical = True
            has_required_hippocampal = True
            has_required_subcortical = True
            subject = (participant_id, session_id)
            
            # Check required cortical files
            for file_pattern in required_cortical_files:
                file_path = os.path.join(
                    self.micapipe_directory,
                    participant_id,
                    session_id,
                    file_pattern.format(participant_id=participant_id, session_id=session_id)
                )
                if not os.path.exists(file_path):
                    subject_missing_files.append(file_path)
                    print(f"Missing required cortical file for subject {participant_id}/{session_id}: {file_path}")
                    has_required_cortical = False

            for file_pattern in required_freesurfer_files:
                file_path = os.path.join(
                    self.freesurfer_directory,
                    file_pattern.format(participant_id=participant_id, session_id=session_id)
                )
                if not os.path.exists(file_path):
                    subject_missing_files.append(file_path)
                    print(f"Missing required cortical file for subject {participant_id}/{session_id}: {file_path}")
                    has_required_cortical = False
            
            # Check hippocampal files if enabled
            if self.hippocampus and self.hippunfold_directory:
                
                # Select required files based on detected version
                if self.hippunfold_version == 2:
                    req_patterns = required_hippocampal_files_v2
                    base_dir_for_file = self.hippunfold_directory
                else:
                    req_patterns = required_hippocampal_files
                    base_dir_for_file = os.path.join(self.hippunfold_directory, "hippunfold")

                for file_pattern in req_patterns:
                    file_path = os.path.join(
                        base_dir_for_file,
                        participant_id,
                        session_id,
                        file_pattern.format(participant_id=participant_id, session_id=session_id)
                    )
                    
                    if not os.path.exists(file_path):
                        subject_missing_files.append(file_path)
                        # Only spam print if verbose or just rely on summary
                        # print(f"Missing required hippocampal file for subject {participant_id}/{session_id}: {file_path}")
                        has_required_hippocampal = False

            # Check subcortical files if enabled
            if self.subcortical and self.freesurfer_directory:
                for file_pattern in required_subcortical_files:
                    file_path = os.path.join(
                        self.freesurfer_directory,
                        file_pattern.format(participant_id=participant_id, session_id=session_id)
                    )
                    if not os.path.exists(file_path):
                        subject_missing_files.append(file_path)
                        print(f"Missing required subcortical file for subject {participant_id}/{session_id}: {file_path}")
                        has_required_subcortical = False
            
            # Update structure-specific lists
            if has_required_cortical:
                self.valid_subjects['structures']['cortex'].append(subject)
            
            if has_required_hippocampal:
                self.valid_subjects['structures']['hippocampus'].append(subject)
            
            if has_required_subcortical:
                self.valid_subjects['structures']['subcortical'].append(subject)
            
            # Check feature-specific files and structure availability
            for feature in self.features:
                feature_complete = True
                # Keep track of which structures have this feature
                feature_has_cortical = has_required_cortical
                feature_has_hippocampal = has_required_hippocampal
                feature_has_subcortical = has_required_subcortical
                
                if feature not in feature_files:
                    if verbose:
                        print(f"Warning: Unknown feature '{feature}'. No file checks performed.")
                    continue
                
                for file_pattern in feature_files[feature]:
                    # Handle file patterns with hemisphere placeholders
                    if "{hemi}" in file_pattern:
                        for hemi in ["L", "R"]:
                            # Handle file patterns with surface type placeholders
                            if "{surfacetype}" in file_pattern:
                                for surfacetype in ["midthickness", "white"]:
                                    file_path = os.path.join(
                                        self.micapipe_directory,
                                        participant_id,
                                        session_id,
                                        file_pattern.format(
                                            participant_id=participant_id, 
                                            session_id=session_id,
                                            hemi=hemi,
                                            surfacetype=surfacetype
                                        )
                                    )
                                    if not os.path.exists(file_path):
                                        subject_missing_files.append(file_path)
                                        feature_complete = False
                                        feature_has_cortical = False  # Surface files affect cortical processing
                            else:
                                # Pattern has hemi but not surfacetype
                                file_path = os.path.join(
                                    self.micapipe_directory,
                                    participant_id,
                                    session_id,
                                    file_pattern.format(
                                        participant_id=participant_id, 
                                        session_id=session_id,
                                        hemi=hemi
                                    )
                                )
                                if not os.path.exists(file_path):
                                    subject_missing_files.append(file_path)
                                    feature_complete = False
                                    feature_has_cortical = False
                    else:
                        # Pattern has neither hemi nor surfacetype (volume file)
                        file_path = os.path.join(
                            self.micapipe_directory,
                            participant_id,
                            session_id,
                            file_pattern.format(
                                participant_id=participant_id, 
                                session_id=session_id
                            )
                        )
                        if not os.path.exists(file_path):
                            subject_missing_files.append(file_path)
                            feature_complete = False
                            # Volume files affect all structures
                            feature_has_cortical = False
                            feature_has_hippocampal = False
                            feature_has_subcortical = False
                
                # Update feature availability
                if feature_complete:
                    feature_availability[feature] += 1
                    self.valid_subjects[feature]['all'].append(subject)
                    
                    # Update structure-specific availability for this feature
                    if feature_has_cortical:
                        self.valid_subjects[feature]['structures']['cortex'].append(subject)
                    
                    if feature_has_hippocampal and not feature.lower().endswith('*blur') and feature.lower() != 'fmri':
                        # Blur features and fMRI don't apply to hippocampus
                        self.valid_subjects[feature]['structures']['hippocampus'].append(subject)
                    
                    if feature_has_subcortical and not feature.lower().endswith('*blur') and feature.lower() != 'fmri' and feature.lower() != 'sa':
                        # Blur features, fMRI, and SA don't apply to subcortical
                        self.valid_subjects[feature]['structures']['subcortical'].append(subject)
                else:
                    subject_missing_features.add(feature)
            
            # Track missing features for this subject
            if subject_missing_features:
                missing_features[subject] = subject_missing_features
                
            # Track missing files and complete subjects
            if subject_missing_files:
                missing_files[subject] = subject_missing_files
            else:
                subjects_with_complete_data.append(subject)
        
        # Apply staged correlation QC BEFORE any dataset-level normalization.
        # This changes only the feature-specific lists consulted by
        # subject_normalization_modalities() and the per-feature processing loops;
        # the control remains in ``base`` for all of its unaffected features.
        self._apply_control_feature_exclusions(self.valid_subjects)
        for feature in self.features:
            if feature in self.valid_subjects:
                feature_availability[feature] = len(self.valid_subjects[feature]["all"])

        # Store the results
        self.subjects_with_complete_data = subjects_with_complete_data
        self.missing_files = missing_files
        self.missing_features = missing_features
        self.feature_availability = feature_availability
        self.source_valid_subjects = copy.deepcopy(self.valid_subjects)
        self.source_feature_availability = dict(feature_availability)
        
        # Report findings
        if verbose:
            if subjects_with_complete_data:
                print(f"Found {len(subjects_with_complete_data)} subjects with complete data for all requested features.")
            
            if len(self.valid_subjects['structures']['cortex']) < len(self.valid_subjects['base']):
                print(f"Warning: {len(self.valid_subjects['base']) - len(self.valid_subjects['structures']['cortex'])} subjects are missing required cortical files.")
            
            if self.hippocampus and len(self.valid_subjects['structures']['hippocampus']) < len(self.valid_subjects['base']):
                print(f"Warning: {len(self.valid_subjects['base']) - len(self.valid_subjects['structures']['hippocampus'])} subjects are missing required hippocampal files.")
                
            if self.subcortical and len(self.valid_subjects['structures']['subcortical']) < len(self.valid_subjects['base']):
                print(f"Warning: {len(self.valid_subjects['base']) - len(self.valid_subjects['structures']['subcortical'])} subjects are missing required subcortical files.")
        
        # Always print feature availability summary
        print("\nFeature availability summary:")
        total_subjects = len(self.valid_subjects['base'])
        denom = total_subjects if total_subjects > 0 else 1

        for feature in self.features:
            if feature in feature_availability:
                avail_count = feature_availability[feature]
                print(f"  {feature}: {avail_count}/{total_subjects} subjects ({(avail_count/denom)*100:.1f}%)")
                
                # Print structure-specific availability for this feature
                if verbose:
                    cortex_count = len(self.valid_subjects[feature]['structures']['cortex'])
                    print(f"    - cortex: {cortex_count}/{total_subjects} subjects ({(cortex_count/denom)*100:.1f}%)")
                    
                    if self.hippocampus and not feature.lower().endswith('*blur') and feature.lower() != 'fmri':
                        hippo_count = len(self.valid_subjects[feature]['structures']['hippocampus'])
                        print(f"    - hippocampus: {hippo_count}/{total_subjects} subjects ({(hippo_count/denom)*100:.1f}%)")
                    
                    if self.subcortical and not feature.lower().endswith('*blur') and feature.lower() != 'fmri' and feature.lower() != 'sa':
                        subcort_count = len(self.valid_subjects[feature]['structures']['subcortical'])
                        print(f"    - subcortical: {subcort_count}/{total_subjects} subjects ({(subcort_count/denom)*100:.1f}%)")
        
        # Print structure availability summary
        print("\nStructure availability summary:")
        print(f"  cortex: {len(self.valid_subjects['structures']['cortex'])}/{total_subjects} subjects ({(len(self.valid_subjects['structures']['cortex'])/denom)*100:.1f}%)")
        if self.hippocampus:
            print(f"  hippocampus: {len(self.valid_subjects['structures']['hippocampus'])}/{total_subjects} subjects ({(len(self.valid_subjects['structures']['hippocampus'])/denom)*100:.1f}%)")
        if self.subcortical:
            print(f"  subcortical: {len(self.valid_subjects['structures']['subcortical'])}/{total_subjects} subjects ({(len(self.valid_subjects['structures']['subcortical'])/denom)*100:.1f}%)")
        
        # Only print detailed missing features if verbose
        if verbose and missing_features:
            print(f"\nMissing features by subject:")
            for (pid, sid), missing_feats in missing_features.items():
                print(f"  {pid}/{sid}: missing {', '.join(missing_feats)}")
        
        return self

    def process(
        self,
        output_directory,
        features,
        cortical_smoothing=5,
        hippocampal_smoothing=2,
        env=None,
        verbose=True,
        n_jobs=None,
        normalization="ravel",
        skip_existing=True,
        hippocampus_only=False,
    ):
        """
        Process the dataset with specified features and smoothing parameters using joblib parallelization.
        Logs all output to a log file in each subject's session directory.
        
        Parameters:
        -----------
        output_directory : str
            Directory to store processed data
        features : list or None, default=None
            List of features to process. If None, uses previously added features.
        cortical_smoothing : int, default=5
            Smoothing parameter for cortical features
        hippocampal_smoothing : int, default=2
            Smoothing parameter for hippocampal features
        env : object
            Environment object containing paths to required tools
        verbose : bool, default=True
            If True, prints detailed processing information
        n_jobs : int, optional
            Number of parallel jobs to run. If None, uses all available CPU cores.
            If 1, runs sequentially.
        normalization : {"ravel", "whitestripe"}, default="ravel"
            T1w/FLAIR normalization source used for volume-to-surface and
            subcortical extraction. "ravel" runs WhiteStripe then RAVEL.
            "whitestripe" runs only WhiteStripe and uses those outputs.
        hippocampus_only : bool, default=False
            Regenerate only HippUnfold-dependent structural and feature maps.
            Matching non-hippocampal maps and normalized volumes must already
            exist in ``output_directory``.
            
        Returns:
        --------
        self
            The dataset object
        """
        # Use the same consistent feature mapping
        feature_mapping = {
            "fa": "FA",
            "adc": "ADC",
            "thickness": "thickness",
            "sa": "SA",
            "flair": "FLAIR",
            "qt1": "qT1",
            "flair*blur": "FLAIR*blur",
            "qt1*blur": "qT1*blur",
            "fmri": "fMRI",
            "t1w": "T1w",
            "t1w*blur": "T1w*blur",
            "adc*blur": "ADC*blur",
            "fa*blur": "FA*blur",
        }
        
        # Normalize features consistently
        normalized_features = []
        for feature in features:
            feature_lower = str(feature).lower()
            if feature_lower in feature_mapping:
                normalized_features.append(feature_mapping[feature_lower])
            else:
                normalized_features.append(str(feature))
                if verbose:
                    print(f"Warning: Unknown feature '{feature}' in process method")
        
        features = normalized_features
        # Composable normalization: decompose the (possibly composite) label into
        # a subject-level phase (none/whitestripe/wmmean) and a dataset-level phase
        # (none/ravel/nyul). The subject phase writes desc-<subject_norm>; the
        # dataset phase reads that and writes desc-<composite>; extraction (via
        # get_normalized_modality_path) reads desc-<composite>.
        subject_norm, dataset_norm = decompose_normalization_label(normalization)
        normalization = compose_normalization_label(subject_norm, dataset_norm)
        self.intensity_normalization = normalization
        self.cortical_smoothing = cortical_smoothing
        self.hippocampal_smoothing = hippocampal_smoothing
        self.add_features(*features, verbose=verbose)

        if hippocampus_only and not (self.hippocampus and self.hippunfold_directory):
            raise ValueError(
                "hippocampus_only=True requires an enabled HippUnfold dataset"
            )

        if verbose:
            print(f"Processing dataset {self.name} with cortical smoothing {cortical_smoothing} and hippocampal smoothing {hippocampal_smoothing}.")
            if requested_ravel_modalities(features):
                print(f"T1w/FLAIR extraction will use subject-norm '{subject_norm}' + "
                      f"dataset-norm '{dataset_norm}' (desc-{normalization}).")
        
        if env is None:
            raise ValueError("Environment (zbenv) must be provided for processing.")

        if not os.path.exists(output_directory):
            os.makedirs(output_directory, exist_ok=True)
            print(f"Created output directory: {output_directory}")
        else:
            print(f"Output directory already exists: {output_directory}")

        # Pin the process working directory to a stable, absolute location that
        # persists for the whole run. Long multi-stage runs have been observed to
        # lose their launch CWD (os.getcwd() -> FileNotFoundError, and relative
        # PYTHONPATH breaks child Python startup); anchoring CWD here makes the
        # pipeline immune regardless of the cause. All pipeline paths are absolute.
        try:
            os.chdir(os.path.abspath(output_directory))
        except OSError:
            pass

        # Timestamp for all logs
        timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        
        # NO MAIN LOG FILE CREATION - REMOVED

        # Determine number of jobs
        if n_jobs is None:
            n_jobs = env.num_threads or multiprocessing.cpu_count()

        def _coerce_job_count(requested_jobs, n_items):
            if n_items <= 0:
                return 1
            requested_jobs = int(requested_jobs)
            if requested_jobs == -1:
                requested = multiprocessing.cpu_count()
            else:
                requested = max(1, requested_jobs)
            return min(requested, n_items)
        
        # Use all base valid subjects
        valid_subjects_to_process = self.valid_subjects['base']
        hipp_res = "8k" if self.hippunfold_version == 2 else "0p5mm"

        # Per-subject REUSE: skip subjects whose outputs (maps/ + structural/)
        # already exist in this base, so a re-run or an enlarged cohort only
        # processes the genuinely-missing subjects instead of everything. This is
        # only safe for subject-level normalizations (dataset_norm == "none", i.e.
        # none/whitestripe/wmmean) where each subject is processed independently;
        # RAVEL/Nyul (dataset_norm in {ravel,nyul}) are DATASET-level fits that need
        # every control, so they are never filtered here.
        if skip_existing and hippocampus_only:
            def _hippocampal_outputs_complete(subject):
                pid, sid = subject
                session_dir = os.path.join(output_directory, pid, sid)
                map_dir = os.path.join(session_dir, "maps", "hippocampus")
                structural_dir = os.path.join(session_dir, "structural")
                if not os.path.isdir(map_dir) or not os.path.isdir(structural_dir):
                    return False

                output_tokens = {
                    "thickness": "thickness",
                    "sa": "SA",
                    "t1map": "qT1",
                    "qt1": "qT1",
                    "adc": "ADC",
                    "fa": "FA",
                    "flair": "FLAIR",
                    "t1w": "T1w",
                }
                valid_features = [
                    feature for feature in self.features
                    if not feature.endswith("*blur")
                    and subject
                    in self.valid_subjects[feature]["structures"]["hippocampus"]
                ]
                expected_maps = []
                for feature in valid_features:
                    token = output_tokens.get(feature.lower(), feature.lower())
                    for hemi in ("L", "R"):
                        expected_maps.append(os.path.join(
                            map_dir,
                            f"{pid}_{sid}_hemi-{hemi}_den-{hipp_res}_label-hipp_"
                            f"midthickness_feature-{token}_smooth-"
                            f"{hippocampal_smoothing}mm.func.gii",
                        ))
                expected_surfaces = [
                    os.path.join(
                        structural_dir,
                        f"{pid}_{sid}_hemi-{hemi}_space-{space}_den-{hipp_res}_"
                        "label-hipp_midthickness.surf.gii",
                    )
                    for hemi in ("L", "R")
                    for space in ("T1w", "unfold")
                ]
                return bool(expected_maps) and all(
                    os.path.isfile(path) for path in expected_maps + expected_surfaces
                )

            all_valid = list(valid_subjects_to_process)
            valid_subjects_to_process = [
                subject for subject in all_valid
                if not _hippocampal_outputs_complete(subject)
            ]
            n_skipped = len(all_valid) - len(valid_subjects_to_process)
            if n_skipped and verbose:
                print(
                    f"Hippocampal reuse: skipping {n_skipped}/{len(all_valid)} "
                    f"complete subjects; processing {len(valid_subjects_to_process)}."
                )
        elif skip_existing and dataset_norm == "none":
            def _outputs_complete(subject):
                pid, sid = subject
                sess = os.path.join(output_directory, pid, sid)
                maps_dir, struct_dir = os.path.join(sess, "maps"), os.path.join(sess, "structural")
                return (os.path.isdir(maps_dir) and bool(os.listdir(maps_dir))
                        and os.path.isdir(struct_dir) and bool(os.listdir(struct_dir)))
            all_valid = list(valid_subjects_to_process)
            valid_subjects_to_process = [s for s in all_valid if not _outputs_complete(s)]
            n_skipped = len(all_valid) - len(valid_subjects_to_process)
            if n_skipped and verbose:
                print(f"Per-subject reuse: skipping {n_skipped}/{len(all_valid)} already-processed "
                      f"subjects; processing {len(valid_subjects_to_process)}.")

        # Threads for the per-subject loop: coerce to the item count, then cap for
        # memory. Set once here; the loop below installs a thread-local stdout proxy
        # (see _ThreadLocalStdout) so concurrent subjects don't clobber each other.
        subject_n_jobs = _coerce_job_count(n_jobs, len(valid_subjects_to_process))
        if subject_n_jobs > PROCESSING_MAX_THREADS:
            if verbose:
                print(f"Capping per-subject threads {subject_n_jobs} -> {PROCESSING_MAX_THREADS} "
                      f"(memory); set env.num_threads to override.")
            subject_n_jobs = PROCESSING_MAX_THREADS
        _stdout_proxy = None    # set to a _ThreadLocalStdout only while the parallel loop runs

        if verbose:
            print(f"Using {subject_n_jobs if subject_n_jobs > 1 else 'sequential'} processing "
                  f"for {len(valid_subjects_to_process)} subjects")

        # Identify blur features
        blur_features = [feature for feature in self.features if feature.endswith("*blur")]
        
        normalization_modalities = requested_ravel_modalities(self.features)
        normalization_failed_subjects = set()

        if normalization_modalities and not hippocampus_only:
            dataset_name = str(self.name).lower()
            is_control_dataset = dataset_name in {"control", "controls", "hc", "healthy_controls", "reference", "normative"}
            normalization_total_threads = max(1, int(env.num_threads or 1))

            def _positive_int_setting(name):
                value = getattr(env, name, None)
                if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
                    return max(1, int(value))
                return None

            configured_normalization_jobs = _positive_int_setting("normalization_jobs")
            default_normalization_jobs = min(
                NORMALIZATION_MAX_JOBS,
                normalization_total_threads,
                max(1, len(valid_subjects_to_process)),
            )
            normalization_n_jobs = min(
                configured_normalization_jobs or default_normalization_jobs,
                normalization_total_threads,
                max(1, len(valid_subjects_to_process)),
            )
            configured_threads_per_job = _positive_int_setting(
                "normalization_threads_per_job"
            )
            available_threads_per_job = max(
                1,
                normalization_total_threads // normalization_n_jobs,
            )
            normalization_threads = min(
                configured_threads_per_job or available_threads_per_job,
                available_threads_per_job,
            )
            configured_tmp_dir = getattr(env, "ravel_tmp_dir", None)
            normalization_tmp_dir = (
                os.fspath(configured_tmp_dir)
                if isinstance(configured_tmp_dir, (str, os.PathLike))
                else None
            )
            micapipe_directory = self.micapipe_directory
            raw_data_directory = getattr(self, "raw_data_directory", None)

            if verbose:
                role = "control" if is_control_dataset else "patient"
                print(
                    f"Preparing {role} T1w/FLAIR normalization for "
                    f"{len(valid_subjects_to_process)} subjects: {', '.join(normalization_modalities)}"
                )
                print(
                    f"Running T1w/FLAIR normalization with {normalization_n_jobs} "
                    f"parallel session job(s) ({normalization_threads} "
                    "ANTs/SynthSeg threads per job)."
                )

            def subject_normalization_modalities(subject):
                subject_modalities = []
                for modality in normalization_modalities:
                    modality_features = [modality, f"{modality}*blur"]
                    if any(
                        feature in self.valid_subjects
                        and subject in self.valid_subjects[feature]["all"]
                        for feature in modality_features
                    ):
                        subject_modalities.append(modality)
                return subject_modalities

            # Reuse smoothing-independent normalized volumes across smoothing/sampling
            # arms (see _norm_cache_* helpers). Keyed by the base path with the smoothing
            # token stripped, so the per-fold _reffold{k} tag is preserved -> no CV leak.
            _final_desc = resolve_normalization_desc(normalization)
            _norm_cache_dir = _norm_cache_dir_for(output_directory)
            _subjects_modalities = {
                subject: subject_normalization_modalities(subject)
                for subject in valid_subjects_to_process
            }

            norm_hydrated = False
            if _norm_cache_dir:
                norm_hydrated = _norm_cache_hydrate(
                    _norm_cache_dir, output_directory, _subjects_modalities, _final_desc,
                    verbose=verbose,
                )
                if norm_hydrated and verbose:
                    print(f"Reusing cached desc-{_final_desc} T1w/FLAIR volumes "
                          f"(smoothing-independent) from {_norm_cache_dir}; skipping "
                          f"WhiteStripe/RAVEL/Nyul.")

            def prepare_subject_whitestripe(subject):
                participant_id, session_id = subject
                subject_modalities = subject_normalization_modalities(subject)
                if not subject_modalities:
                    if verbose:
                        print(f"Skipping T1w/FLAIR normalization for {participant_id}/{session_id}; no requested T1w/FLAIR modalities are available.")
                    return participant_id, session_id, True, None, []

                # Cache reuse: the full desc-<norm> volume set was hydrated from a
                # sibling smoothing/sampling arm, so neither WhiteStripe nor the dataset
                # fit needs to run. (All-or-nothing: norm_hydrated is only True when every
                # subject/modality was present, so the dataset fit never lacks an input.)
                if norm_hydrated:
                    return participant_id, session_id, True, None, subject_modalities

                session_output_dir = os.path.join(output_directory, participant_id, session_id)
                maps_dir = os.path.join(session_output_dir, "maps")
                os.makedirs(maps_dir, exist_ok=True)

                # RAVEL needs the SynthSeg CSF control region regardless of the
                # subject-level normalization; the identity/wmmean paths do not
                # write it, so fill it in when dataset_norm == "ravel".
                def _ensure_ravel_csf():
                    if dataset_norm == "ravel":
                        ensure_synthseg_csf(
                            participant_id=participant_id,
                            session_id=session_id,
                            output_dir=maps_dir,
                            micapipe_dir=micapipe_directory,
                            modalities=subject_modalities,
                            threads=normalization_threads,
                            verbose=verbose,
                        )

                if subject_norm == "none":
                    # Self-normalization: link raw nativepro T1w/FLAIR as desc-none.
                    # The dataset-level phase (ravel/nyul) reads this intermediate.
                    try:
                        prepare_t1w_flair_identity(
                            participant_id=participant_id,
                            session_id=session_id,
                            output_dir=maps_dir,
                            micapipe_dir=micapipe_directory,
                            modalities=subject_modalities,
                            verbose=verbose,
                        )
                        _ensure_ravel_csf()
                        return participant_id, session_id, True, None, subject_modalities
                    except Exception as identity_error:
                        import traceback
                        return (
                            participant_id,
                            session_id,
                            False,
                            f"{identity_error}\n{traceback.format_exc()}",
                            subject_modalities,
                        )

                if subject_norm == "wmmean":
                    # WM-mean image normalization: SynthSeg WM mask + whole-WM
                    # tissue mean/SD -> desc-wmmean (writes the SynthSeg label too).
                    try:
                        prepare_t1w_flair_wmmean(
                            participant_id=participant_id,
                            session_id=session_id,
                            output_dir=maps_dir,
                            micapipe_dir=micapipe_directory,
                            modalities=subject_modalities,
                            threads=normalization_threads,
                            verbose=verbose,
                        )
                        _ensure_ravel_csf()
                        return participant_id, session_id, True, None, subject_modalities
                    except Exception as wmmean_error:
                        import traceback
                        return (participant_id, session_id, False,
                                f"{wmmean_error}\n{traceback.format_exc()}", subject_modalities)

                # subject_norm == "whitestripe": WhiteStripe writes desc-whitestripe
                # plus the SynthSeg label and CSF mask (RAVEL-ready).
                session_tmp_dir = os.path.join(
                    session_output_dir,
                    f"tmp_norm_{participant_id}_{session_id}_{os.getpid()}_{threading.get_ident()}")
                os.makedirs(session_tmp_dir, exist_ok=True)

                try:
                    completed_modalities = []
                    modality_errors = []
                    if verbose:
                        print(
                            f"Normalizing {participant_id}/{session_id} with WhiteStripe "
                            f"({', '.join(subject_modalities)})..."
                        )
                    for modality in subject_modalities:
                        try:
                            prepare_t1w_flair_whitestripe(
                                participant_id=participant_id,
                                session_id=session_id,
                                output_dir=maps_dir,
                                micapipe_dir=micapipe_directory,
                                tmp_dir=session_tmp_dir,
                                raw_dir=raw_data_directory,
                                modalities=[modality],
                                threads=normalization_threads,
                                verbose=verbose,
                            )
                            completed_modalities.append(modality)
                        except Exception as modality_error:
                            modality_errors.append(f"{modality}: {modality_error}")
                            if verbose:
                                print(
                                    f"Error during {modality} WhiteStripe normalization for "
                                    f"{participant_id}/{session_id}: {modality_error}"
                                )

                    if completed_modalities:
                        error_text = "\n".join(modality_errors) if modality_errors else None
                        return participant_id, session_id, True, error_text, completed_modalities

                    return (
                        participant_id,
                        session_id,
                        False,
                        "No requested modalities completed WhiteStripe normalization. "
                        + "\n".join(modality_errors),
                        subject_modalities,
                    )
                except Exception as e:
                    import traceback
                    return participant_id, session_id, False, f"{e}\n{traceback.format_exc()}", subject_modalities
                finally:
                    if os.path.exists(session_tmp_dir):
                        shutil.rmtree(session_tmp_dir, ignore_errors=True)

            if normalization_n_jobs == 1:
                whitestripe_results = [
                    prepare_subject_whitestripe(subject)
                    for subject in valid_subjects_to_process
                ]
            else:
                whitestripe_results = Parallel(n_jobs=normalization_n_jobs, backend="threading")(
                    delayed(prepare_subject_whitestripe)(subject)
                    for subject in valid_subjects_to_process
                )

            normalized_modalities_by_subject = {}
            for participant_id, session_id, success, error, subject_modalities in whitestripe_results:
                subject = (participant_id, session_id)
                if not success:
                    normalization_failed_subjects.add(subject)
                    if verbose:
                        print(f"Error during WhiteStripe normalization for {participant_id}/{session_id}: {error}")
                else:
                    normalized_modalities_by_subject[subject] = set(subject_modalities)

            normalization_subjects = [
                subject for subject in valid_subjects_to_process
                if subject not in normalization_failed_subjects
                and normalized_modalities_by_subject.get(subject)
            ]
            if not normalization_subjects:
                raise ValueError("No subjects completed WhiteStripe normalization.")

            if dataset_norm == "none":
                if verbose:
                    source_label = {
                        "none": "raw (self-normalized) nativepro",
                        "wmmean": "WM-mean-normalized",
                    }.get(subject_norm, "WhiteStripe-normalized")
                    print(
                        f"No dataset-level harmonization; downstream extraction will use "
                        f"{source_label} (desc-{normalization}) T1w/FLAIR images."
                    )
            elif dataset_norm == "nyul":
                if is_control_dataset:
                    if norm_hydrated:
                        if verbose:
                            print(f"Reusing cached desc-{normalization} Nyul control volumes; skipping fit.")
                    else:
                        if verbose:
                            print(f"Fitting control Nyul-Udupa standard scale on desc-{subject_norm} "
                                  f"images before extraction (writing desc-{normalization})...")
                        try:
                            fit_and_apply_nyul_to_controls(
                                subjects=normalization_subjects,
                                output_directory=output_directory,
                                micapipe_dir=micapipe_directory,
                                modalities=normalization_modalities,
                                verbose=verbose,
                                input_desc=subject_norm,
                                output_desc=normalization,
                                n_jobs=normalization_threads,
                            )
                        except Exception as e:
                            if verbose:
                                print(f"Error fitting/applying Nyul control model: {e}")
                                import traceback
                                traceback.print_exc()
                            raise
                else:
                    if verbose:
                        print("Applying fitted control Nyul standard scale to patient dataset...")
                    for subject in normalization_subjects:
                        participant_id, session_id = subject
                        subject_modalities = sorted(normalized_modalities_by_subject.get(subject, set()))
                        if not subject_modalities:
                            continue
                        if norm_hydrated:
                            continue                    # cached (hydrated from a sibling arm)
                        try:
                            apply_nyul_model_to_subject(
                                participant_id=participant_id,
                                session_id=session_id,
                                output_directory=output_directory,
                                micapipe_dir=micapipe_directory,
                                modalities=subject_modalities,
                                verbose=verbose,
                                input_desc=subject_norm,
                                output_desc=normalization,
                            )
                        except Exception as e:
                            normalization_failed_subjects.add(subject)
                            if verbose:
                                print(f"Error applying Nyul model for {participant_id}/{session_id}: {e}")
            elif is_control_dataset:
                if verbose:
                    print(f"Fitting control RAVEL model on desc-{subject_norm} images "
                          f"(writing desc-{normalization}) before any volume-to-surface extraction...")
                try:
                    for modality in normalization_modalities:
                        modality_subjects = [
                            subject for subject in normalization_subjects
                            if modality in normalized_modalities_by_subject.get(subject, set())
                        ]
                        if not modality_subjects:
                            if verbose:
                                print(f"Skipping control RAVEL model for {modality}; no subjects have that modality.")
                            continue
                        if norm_hydrated:
                            if verbose:
                                print(f"Reusing cached desc-{normalization} {modality} for controls; skipping RAVEL fit.")
                            continue                    # hydrated from a sibling smoothing arm
                        fit_and_apply_ravel_to_controls(
                            subjects=modality_subjects,
                            output_directory=output_directory,
                            micapipe_dir=micapipe_directory,
                            modalities=[modality],
                            k=1,
                            verbose=verbose,
                            n_jobs=normalization_n_jobs,
                            threads=normalization_threads,
                            input_desc=subject_norm,
                            output_desc=normalization,
                            total_threads=normalization_total_threads,
                            tmp_dir=normalization_tmp_dir,
                        )
                except Exception as e:
                    if verbose:
                        print(f"Error fitting/applying RAVEL control model: {e}")
                        import traceback
                        traceback.print_exc()
                    raise
            else:
                if verbose:
                    print("Applying fitted control RAVEL model to patient dataset...")
                patient_ravel_models = (
                    None
                    if norm_hydrated
                    else load_ravel_models(output_directory, normalization_modalities)
                )

                def apply_subject_ravel(subject):
                    participant_id, session_id = subject
                    subject_modalities = sorted(normalized_modalities_by_subject.get(subject, set()))
                    if not subject_modalities:
                        return participant_id, session_id, True, None
                    if norm_hydrated:
                        return participant_id, session_id, True, None   # cached (hydrated)
                    try:
                        apply_ravel_model_to_subject(
                            participant_id=participant_id,
                            session_id=session_id,
                            output_directory=output_directory,
                            micapipe_dir=micapipe_directory,
                            modalities=subject_modalities,
                            threads=normalization_threads,
                            verbose=verbose,
                            input_desc=subject_norm,
                            output_desc=normalization,
                            models=patient_ravel_models,
                        )
                        return participant_id, session_id, True, None
                    except Exception as e:
                        import traceback
                        return participant_id, session_id, False, f"{e}\n{traceback.format_exc()}"

                if normalization_n_jobs == 1:
                    ravel_results = [apply_subject_ravel(subject) for subject in normalization_subjects]
                else:
                    if verbose:
                        print(f"Applying patient RAVEL normalization with {normalization_n_jobs} parallel jobs...")
                    ravel_results = Parallel(n_jobs=normalization_n_jobs, backend="threading")(
                        delayed(apply_subject_ravel)(subject)
                        for subject in normalization_subjects
                    )

                for participant_id, session_id, success, error in ravel_results:
                    if not success:
                        normalization_failed_subjects.add((participant_id, session_id))
                        if verbose:
                            print(f"Error applying RAVEL model for {participant_id}/{session_id}: {error}")

            # Persist the freshly-computed smoothing-independent volumes (+ fit model)
            # so sibling smoothing/sampling arms of this SAME (norm, exclusion, fold)
            # reuse them. No-op when hydrated from cache (volumes are symlinks) or when
            # the base path is not a zbrains_ dir.
            if _norm_cache_dir and not norm_hydrated:
                _norm_cache_save(_norm_cache_dir, output_directory,
                                 normalized_modalities_by_subject, _final_desc, verbose=verbose)

        # Define a function to process a single subject
        def process_single_subject(subject):
            participant_id, session_id = subject
            if subject in normalization_failed_subjects:
                if verbose:
                    print(f"Skipping {participant_id}/{session_id}; T1w/FLAIR normalization failed.")
                return participant_id, session_id, False
            
            # Create session-specific directories
            session_output_dir = os.path.join(output_directory, participant_id, session_id)
            os.makedirs(session_output_dir, exist_ok=True)
            
            # Create session log file - kept in subject directory
            log_file = os.path.join(session_output_dir, f"{participant_id}_{session_id}_processing_{timestamp}.log")
            
            # Create session-specific tmp directory
            # Per-process/thread-unique scratch dir: two machines sharing storage (or
            # two threads) may build the SAME shared base (e.g. the 0mm ROBPCA detector
            # base) concurrently -- a deterministic tmp path would let one process's
            # rmtree race another's writes ("Directory not empty"). Unique suffix +
            # tolerant cleanup makes it collision-proof.
            session_tmp_dir = os.path.join(
                session_output_dir,
                f"tmp_{participant_id}_{session_id}_{os.getpid()}_{threading.get_ident()}")
            os.makedirs(session_tmp_dir, exist_ok=True)
            
            # Per-subject logging. In the PARALLEL loop a thread-local proxy routes
            # THIS thread's prints to its own log file (no shared sys.stdout swap, so
            # concurrent subjects can't clobber each other); SERIALLY we swap
            # sys.stdout directly so output still mirrors live to the terminal.
            log_redirect = None
            original_stdout = None
            try:
                if _stdout_proxy is not None:
                    # parallel: log to file only (N interleaved terminals = noise)
                    log_redirect = LogRedirect(log_file, mirror_terminal=False)
                    _stdout_proxy.set_target(log_redirect)
                else:
                    original_stdout = sys.stdout
                    log_redirect = LogRedirect(log_file)
                    sys.stdout = log_redirect

                print(f"Processing subject {participant_id}/{session_id}...")
                print(f"Started at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
                print(f"Features: {features}")
                print(f"Cortical smoothing: {cortical_smoothing}mm")
                print(f"Hippocampal smoothing: {hippocampal_smoothing}mm")
                print(f"Log file: {log_file}")
                
                # Add dataset-level information to each subject log
                print(f"Dataset: {self.name}")
                print(f"Total subjects in dataset: {len(valid_subjects_to_process)}")
                print(f"Processing job: {n_jobs if n_jobs > 1 else 'sequential'}")
                print("-" * 50)
                
                try:
                    # Redirect structural/ to a base-independent per-subject cache
                    # (geometry doesn't depend on the base) BEFORE any structural
                    # write, so surfaces + Laplace/SWM are built once and reused.
                    _link_structural_to_cache(
                        output_directory,
                        participant_id,
                        session_id,
                        verbose=verbose,
                        hippunfold_directory=(
                            self.hippunfold_directory if self.hippocampus else None
                        ),
                        hippunfold_version=(
                            self.hippunfold_version
                            if self.hippocampus and self.hippunfold_directory
                            else None
                        ),
                    )
                    # Partial V1 -> V2 migration seeds non-hippocampal geometry
                    # from the legacy base, then copies only the V2 surfaces here.
                    if hippocampus_only:
                        self._copy_hippocampal_structural_files(
                            participant_id,
                            session_id,
                            output_directory,
                            verbose=verbose,
                        )
                    else:
                        self._copy_structural_files(
                            participant_id,
                            session_id,
                            output_directory,
                            verbose=verbose,
                        )
                        self._create_midline_from_freesurfer(
                            participant_id,
                            session_id,
                            output_directory,
                            verbose=verbose,
                        )

                    if self.freesurfer_directory and not hippocampus_only:
                        generate_superficial_white_matter(
                            participant_id=participant_id,
                            session_id=session_id,
                            output_directory=output_directory,
                            workbench_path=env.connectome_workbench_path,
                            micapipe_directory=self.micapipe_directory,
                            freesurfer_directory=self.freesurfer_directory,
                            tmp_dir=session_tmp_dir,
                            verbose=verbose,
                        )
                    
                    # Apply blurring to features that need it
                    if blur_features and not hippocampus_only:
                        # Check which base features are available for this subject
                        available_features = []
                        # Iterate the ORIGINAL blur features to check validity
                        for blur_feature in blur_features:
                            if subject in self.valid_subjects[blur_feature]['all']:
                                # Append the base feature name to available_features for processing
                                available_features.append(blur_feature.replace("*blur", ""))
                        
                        if available_features:  # If ANY features are available, proceed with blurring
                            if verbose:
                                print(f"  Applying additional blur processing for features: {', '.join([f+'*blur' for f in available_features])}")
                            
                            apply_blurring(
                                participant_id=participant_id,
                                session_id=session_id,
                                features=available_features,  # Only process available features
                                output_directory=output_directory,
                                workbench_path=env.connectome_workbench_path,
                                micapipe_directory=self.micapipe_directory,
                                freesurfer_directory=self.freesurfer_directory,
                                tmp_dir=session_tmp_dir,
                                smoothing_fwhm=self.cortical_smoothing,
                                verbose=verbose,
                                normalization=normalization,
                            )

                    # Process cortical features if cortex is enabled
                    if (not hippocampus_only and self.cortex
                            and subject in self.valid_subjects['structures']['cortex']):
                        # Get valid features for cortex for this subject
                        valid_cortical_features = [f for f in self.features if subject in self.valid_subjects[f]['structures']['cortex']]
                        
                        if valid_cortical_features:
                            if verbose:
                                print(f"  Processing cortical data for features: {', '.join(valid_cortical_features)}")
                            
                            apply_cortical_processing(
                                participant_id=participant_id,
                                session_id=session_id,
                                features=valid_cortical_features,
                                output_directory=output_directory,
                                workbench_path=env.connectome_workbench_path,
                                micapipe_directory=self.micapipe_directory,
                                tmp_dir=session_tmp_dir,
                                cortical_smoothing=cortical_smoothing,
                                resolutions=["32k"],#, "5k"],
                                labels=["midthickness", "white"],
                                verbose=verbose,
                                normalization=normalization,
                            )
                
                    # If hippocampus is enabled, process hippocampal data
                    if self.hippocampus and self.hippunfold_directory and subject in self.valid_subjects['structures']['hippocampus']:
                        # Get nonBlur features for hippocampus
                        valid_hipp_features = [f for f in self.features 
                                            if not f.endswith("*blur") 
                                            and subject in self.valid_subjects[f]['structures']['hippocampus']]
                        
                        if valid_hipp_features:
                            apply_hippocampal_processing(
                                participant_id=participant_id,
                                session_id=session_id,
                                features=valid_hipp_features,
                                output_directory=output_directory,
                                workbench_path=env.connectome_workbench_path,
                                micapipe_directory=self.micapipe_directory,
                                hippunfold_directory=self.hippunfold_directory,
                                tmp_dir=session_tmp_dir,
                                smoothing_fwhm=self.hippocampal_smoothing,
                                resolution=hipp_res,
                                hippunfold_version=self.hippunfold_version,
                                verbose=verbose,
                                normalization=normalization,
                            )

                    # If subcortex is enabled, extract subcortical stats
                    if (not hippocampus_only and self.subcortical
                            and self.freesurfer_directory
                            and subject in self.valid_subjects['structures']['subcortical']):
                        # Get nonBlur features for subcortical
                        valid_subcort_features = [f for f in self.features 
                                                if not f.endswith("*blur") 
                                                and subject in self.valid_subjects[f]['structures']['subcortical']]
                        
                        if valid_subcort_features:
                            if verbose:
                                print(f"  Processing subcortical data for features: {', '.join(valid_subcort_features)}")
                            
                            apply_subcortical_processing(
                                participant_id=participant_id,
                                session_id=session_id,
                                features=valid_subcort_features,
                                output_directory=output_directory,
                                micapipe_directory=self.micapipe_directory,
                                freesurfer_directory=self.freesurfer_directory,
                                verbose=verbose,
                                normalization=normalization,
                            )
                    
                    print(f"Completed processing {participant_id}/{session_id} at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
                    return (participant_id, session_id, True)
                
                except Exception as e:
                    if verbose:
                        print(f"Error processing {participant_id}/{session_id}: {e}")
                    import traceback
                    traceback.print_exc()
                    return (participant_id, session_id, False)
                
                finally:
                    # Clean up session-specific tmp directory. Best-effort: a concurrent
                    # peer (two-machine shared storage) writing into a same-named scratch
                    # dir must never crash the run -- a leaked tmp dir is harmless.
                    if os.path.exists(session_tmp_dir):
                        shutil.rmtree(session_tmp_dir, ignore_errors=True)
                        if verbose:
                            print(f"  Cleaned up temporary directory: {session_tmp_dir}")
                    
                    print("-" * 50)
                    print(f"Processing finished at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            finally:
                if _stdout_proxy is not None:
                    _stdout_proxy.clear_target()
                elif original_stdout is not None:
                    sys.stdout = original_stdout
                # Ensure this subject's log is closed
                if log_redirect is not None:
                    log_redirect.close()

            return (participant_id, session_id, True)
        
        # Process subjects. The heavy per-subject work is external wb_command/ANTs
        # subprocesses (GIL released), so THREADS give real speedup with output
        # bit-identical to the serial path -- each subject writes its own maps under
        # its own (subject-keyed) paths and caches, so there is no cross-thread race.
        # A process pool is avoided (loky/fork + numpy double-init -> SIGSEGV).
        if subject_n_jobs > 1 and len(valid_subjects_to_process) > 1:
            print(f"Running parallel processing with {subject_n_jobs} threads for "
                  f"{len(valid_subjects_to_process)} subjects")
            _saved_stdout = sys.stdout
            _stdout_proxy = _ThreadLocalStdout(_saved_stdout)   # visible to process_single_subject
            sys.stdout = _stdout_proxy
            try:
                results = Parallel(n_jobs=subject_n_jobs, backend="threading")(
                    delayed(process_single_subject)(subject)
                    for subject in valid_subjects_to_process
                )
            finally:
                sys.stdout = _saved_stdout
                _stdout_proxy = None
        else:
            print(f"Running sequential processing for {len(valid_subjects_to_process)} subjects")
            results = [process_single_subject(subject) for subject in valid_subjects_to_process]

        # Process results
        failed_subjects = [(pid, sid) for pid, sid, success in results if not success]
        
        # Report results
        if verbose:
            successful_count = len(valid_subjects_to_process) - len(failed_subjects)
            print(f"\nProcessing complete: {successful_count}/{len(valid_subjects_to_process)} subjects successful")
            
            if failed_subjects:
                print(f"Failed subjects ({len(failed_subjects)}):")
                for participant_id, session_id in failed_subjects:
                    print(f"  - {participant_id}/{session_id}")

        # Create a summary file in each subject's directory
        for participant_id, session_id in valid_subjects_to_process:
            session_dir = os.path.join(output_directory, participant_id, session_id)
            summary_file = os.path.join(session_dir, f"processing_summary_{timestamp}.txt")
            
            try:
                with open(summary_file, 'w', encoding='utf-8') as f:
                    f.write(f"===== PROCESSING SUMMARY FOR {participant_id}/{session_id} =====\n")
                    f.write(f"Dataset: {self.name}\n")
                    f.write(f"Completed at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"Total subjects processed: {len(valid_subjects_to_process)}\n")
                    f.write(f"Subjects processed successfully: {len(valid_subjects_to_process) - len(failed_subjects)}/{len(valid_subjects_to_process)}\n")
                    
                    if failed_subjects:
                        f.write(f"Failed subjects ({len(failed_subjects)}):\n")
                        for pid, sid in failed_subjects:
                            f.write(f"  - {pid}/{sid}\n")
                            
                    # Mark this subject's status
                    is_failed = (participant_id, session_id) in failed_subjects
                    f.write(f"\nThis subject status: {'FAILED' if is_failed else 'SUCCESS'}\n")
            except Exception as e:
                if verbose:
                    print(f"Error writing summary for {participant_id}/{session_id}: {e}")

        print(f"Dataset {self.name} processed. Summary files written to each subject directory.")
        
        self.validate(
            features=self.features,
            output_directory=output_directory,
            cortical_smoothing=cortical_smoothing,
            hippocampal_smoothing=hippocampal_smoothing,
            verbose=verbose,
            normalization=normalization,
        )
        return self

    def validate(self, output_directory, features=None, cortical_smoothing=5, hippocampal_smoothing=2,
                 verbose=True, error_threshold=0, warning_threshold=25,
                 only_demographics_subjects=True, normalization=None):
        """
        Validate that all expected processed files exist for each subject.

        The optional normalization argument is accepted so callers can reuse the
        same settings dictionary passed to process(). Validation checks the
        derived files generated from the selected normalization source.

        Returns
        -------
        dict
            {
                "valid_subjects": [(participant_id, session_id), ...],
                "missing_files": { (participant_id, session_id): [paths...] },
                "summary": {...}
            }
        """
        if normalization is not None:
            # Accept composite labels (e.g. noneRavel/whitestripeNyul) as well as
            # the legacy single modes; store the canonical composite desc so
            # downstream path resolution reads the right desc-<composite> volumes.
            subject_norm, dataset_norm = decompose_normalization_label(normalization)
            self.intensity_normalization = compose_normalization_label(subject_norm, dataset_norm)

        feature_case_mapping = {
            "fa": "FA",
            "adc": "ADC",
            "thickness": "thickness",
            "sa": "SA",
            "flair": "FLAIR",
            "qt1": "qT1",
            "t1map": "qT1",
            "flair*blur": "FLAIR*blur",
            "qt1*blur": "qT1*blur",
            "fmri": "fMRI",
            "t1w": "T1w",
            "t1w*blur": "T1w*blur",
            "adc*blur": "ADC*blur",
            "fa*blur": "FA*blur",
        }

        if features is None:
            normalized_features = list(self.features)
        else:
            normalized_features = []
            for feat in features:
                feat_lower = str(feat).lower()
                if feat_lower.endswith("*blur"):
                    base = feat_lower.replace("*blur", "")
                    mapped = feature_case_mapping.get(feat_lower, f"{feature_case_mapping.get(base, feat[:-5])}*blur")
                    normalized_features.append(mapped)
                else:
                    normalized_features.append(feature_case_mapping.get(feat_lower, str(feat)))
        if not normalized_features:
            normalized_features = list(self.features)
        self.features = normalized_features

        if not hasattr(self, "valid_subjects") or "base" not in self.valid_subjects:
            self.check_directories()

        source_valid_subjects = getattr(self, "source_valid_subjects", None)
        if source_valid_subjects is None or any(feature not in source_valid_subjects for feature in self.features):
            self.add_features(*self.features, verbose=False)
            source_valid_subjects = copy.deepcopy(self.valid_subjects)

        def source_has_structure(subject, structure_name):
            if not source_valid_subjects:
                return True
            structures = source_valid_subjects.get("structures", {})
            if structure_name not in structures:
                return True
            return subject in structures[structure_name]

        def source_has_feature_structure(subject, feature_name, structure_name):
            if not source_valid_subjects or feature_name not in source_valid_subjects:
                return True
            feature_data = source_valid_subjects[feature_name]
            structures = feature_data.get("structures", {})
            if structure_name not in structures:
                return True
            return subject in structures[structure_name]

        output_subjects = self.check_output_directories(
            output_directory,
            only_demographics_subjects=only_demographics_subjects,
            verbose=verbose
        )
        if not output_subjects:
            raise ValueError(f"No subject output directories found in {output_directory}. Nothing to validate.")

        if verbose:
            print(f"Validating dataset {self.name} with {len(output_subjects)} output subjects "
                  f"and {len(self.features)} features...")

        feature_output_mapping = {
            "thickness": "thickness",
            "sa": "SA",
            "flair": "FLAIR",
            "adc": "ADC",
            "fa": "FA",
            "qt1": "qT1",
            "t1map": "qT1",
            "fmri": "fMRI",
            "t1w": "T1w",
        }

        feature_meta = []
        for feat in self.features:
            feat_lower = feat.lower()
            is_blur = feat_lower.endswith("*blur")
            base_lower = feat_lower.replace("*blur", "") if is_blur else feat_lower
            cortical_token = feat  # already normalised to processing output naming

            if is_blur:
                if base_lower in {"qt1", "t1map"}:
                    blur_token = "qT1"
                elif base_lower == "flair":
                    blur_token = "FLAIR"
                elif base_lower in {"adc", "fa"}:
                    blur_token = base_lower
                elif base_lower == "t1w":
                    blur_token = "T1w"
                else:
                    blur_token = feat[:-5]
            else:
                blur_token = None

            hippo_token = None if is_blur else feature_output_mapping.get(base_lower, feat)
            subcortical_token = None
            if not is_blur:
                if base_lower in {"thickness", "volume"}:
                     # Map thickness/volume to the "volume" CSV file
                     subcortical_token = "volume"
                elif base_lower not in {"thickness", "volume"}:
                    subcortical_token = feature_output_mapping.get(base_lower, feat.upper())

            feature_meta.append({
                "original": feat,
                "base_lower": base_lower,
                "is_blur": is_blur,
                "cortical_token": cortical_token,
                "blur_token": blur_token,
                "hippo_token": hippo_token,
                "subcortical_token": subcortical_token,
                "requires_cortex": self.cortex,
                "requires_hippocampus": (
                    self.hippocampus and self.hippunfold_directory and not is_blur and hippo_token is not None and base_lower != 'fmri'
                ),
                "requires_subcortical": (
                    self.subcortical and self.freesurfer_directory and not is_blur and subcortical_token is not None and base_lower != 'fmri' and base_lower != 'sa'
                )
            })

        resolutions = ["32k"]#, "5k"]
        surface_labels = ["midthickness", "white"]
        hippocampal_resolution = "8k" if self.hippunfold_version == 2 else "0p5mm"
        blur_suffixes = [
            "_surf-fsnative_desc-raw.func.gii",
            "_surf-fsnative_desc-dist.func.gii",
            "_surf-fsnative_desc-grad.func.gii"
        ]

        processed_valid_subjects = {
            'base': [],
            'structures': {
                'cortex': [],
                'hippocampus': [],
                'subcortical': []
            }
        }
        feature_processed = {
            meta['original']: {
                'all': [],
                'structures': {
                    'cortex': [],
                    'hippocampus': [],
                    'subcortical': []
                }
            } for meta in feature_meta
        }

        missing_files = {}
        missing_features = {}
        feature_availability = {meta['original']: 0 for meta in feature_meta}
        subjects_with_complete_data = []
        all_files_count = 0
        missing_files_count = 0

        for participant_id, session_id in output_subjects:
            subject_missing = []
            subject_missing_features = set()
            bids_id = f"{participant_id}_{session_id}"
            subject_dir = os.path.join(output_directory, participant_id, session_id)
            cortex_dir = os.path.join(subject_dir, "maps", "cortex")
            hippocampus_dir = os.path.join(subject_dir, "maps", "hippocampus")
            subcortical_dir = os.path.join(subject_dir, "maps", "subcortical")
            structural_dir = os.path.join(subject_dir, "structural")

            subject = (participant_id, session_id)
            subject_structures = {
                'cortex': True if self.cortex and source_has_structure(subject, 'cortex') else None,
                'hippocampus': True if self.hippocampus and self.hippunfold_directory and source_has_structure(subject, 'hippocampus') else None,
                'subcortical': True if self.subcortical and self.freesurfer_directory and source_has_structure(subject, 'subcortical') else None
            }
            subject_feature_structures = {}
            for meta in feature_meta:
                feature_name = meta['original']
                subject_feature_structures[feature_name] = {
                    'cortex': (
                        True
                        if meta['requires_cortex']
                        and source_has_feature_structure(subject, feature_name, 'cortex')
                        else None
                    ),
                    'hippocampus': (
                        True
                        if meta['requires_hippocampus']
                        and source_has_feature_structure(subject, feature_name, 'hippocampus')
                        else None
                    ),
                    'subcortical': (
                        True
                        if meta['requires_subcortical']
                        and source_has_feature_structure(subject, feature_name, 'subcortical')
                        else None
                    )
                }

            structural_expected = [
                f"{bids_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
                f"{bids_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
                f"{bids_id}_hemi-L_surf-fsnative_label-sphere.surf.gii",
                f"{bids_id}_hemi-R_surf-fsnative_label-sphere.surf.gii",
                f"{bids_id}_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii",
                f"{bids_id}_hemi-R_space-nativepro_surf-fsnative_label-pial.surf.gii",
                f"{bids_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii",
                f"{bids_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-pial.surf.gii",
                f"{bids_id}_space-nativepro_T1w.nii.gz",
                f"{bids_id}_hemi-L_surf-fsnative_label-medialwall.label.gii",
                f"{bids_id}_hemi-R_surf-fsnative_label-medialwall.label.gii"
            ]
            
            # Only require Laplace file for subjects with at least one
            # source-available blur feature.
            has_subject_blur_features = any(
                meta['is_blur']
                and source_has_feature_structure(subject, meta['original'], 'cortex')
                for meta in feature_meta
            )
            if has_subject_blur_features:
                structural_expected.append(f"{participant_id}_{session_id}-laplace.nii.gz")
            
            if subject_structures['hippocampus'] is True:
                structural_expected.extend([
                    f"{bids_id}_hemi-L_space-unfold_den-{hippocampal_resolution}_label-hipp_midthickness.surf.gii",
                    f"{bids_id}_hemi-R_space-unfold_den-{hippocampal_resolution}_label-hipp_midthickness.surf.gii",
                    f"{bids_id}_hemi-L_space-T1w_den-{hippocampal_resolution}_label-hipp_midthickness.surf.gii",
                    f"{bids_id}_hemi-R_space-T1w_den-{hippocampal_resolution}_label-hipp_midthickness.surf.gii"
                ])

            for filename in structural_expected:
                path = os.path.join(structural_dir, filename)
                all_files_count += 1
                if not os.path.exists(path):
                    subject_missing.append(path)
                    missing_files_count += 1
                    if "hipp" in filename and subject_structures['hippocampus'] is not None:
                        subject_structures['hippocampus'] = False
                    elif "laplace" in filename or filename.endswith(".surf.gii") or filename.endswith(".label.gii"):
                        if subject_structures['cortex'] is not None:
                            subject_structures['cortex'] = False

            if self.cortex:
                for meta in feature_meta:
                    if not meta['requires_cortex']:
                        continue
                    feature_name = meta['original']
                    if subject_feature_structures[feature_name]['cortex'] is not True:
                        continue
                    for hemi in ["L", "R"]:
                        for res in resolutions:
                            for label in (["midthickness"] if meta['is_blur'] or (meta['base_lower'] == 'fmri')  else surface_labels):
                                if meta['base_lower'] == 'fmri':
                                    for fmrifeat in ['rmssd', 'timescales', 'alff', 'falff']:
                                        cortical_path = os.path.join(
                                            cortex_dir,
                                            f"{bids_id}_hemi-{hemi}_surf-fsLR-{res}_label-{label}_feature-"
                                            f"{fmrifeat}_smooth-{cortical_smoothing}mm.func.gii"
                                        )
                                        all_files_count += 1
                                        if not os.path.exists(cortical_path):
                                            subject_missing.append(cortical_path)
                                            missing_files_count += 1
                                            subject_feature_structures[feature_name]['cortex'] = False
                                            subject_missing_features.add(feature_name)
                                else:
                                    cortical_path = os.path.join(
                                        cortex_dir,
                                        f"{bids_id}_hemi-{hemi}_surf-fsLR-{res}_label-{label}_feature-"
                                        f"{meta['cortical_token']}_smooth-{cortical_smoothing}mm.func.gii"
                                    )
                                    all_files_count += 1
                                    if not os.path.exists(cortical_path):
                                        subject_missing.append(cortical_path)                                      
                                        missing_files_count += 1
                                        subject_feature_structures[feature_name]['cortex'] = False
                                        subject_missing_features.add(feature_name)
                    if meta['is_blur'] and meta['blur_token']:
                        for hemi in ["L", "R"]:
                            blur_prefix = os.path.join(
                                cortex_dir,
                                f"{participant_id}_{session_id}_hemi-{hemi}_feature-{meta['blur_token']}*blur"
                            )
                            blur_paths = [f"{blur_prefix}{suffix}" for suffix in blur_suffixes]
                            blur_paths.append(f"{blur_prefix}_surf-fsnative_smooth-{cortical_smoothing}mm.func.gii")
                            for blur_path in blur_paths:
                                all_files_count += 1
                                if not os.path.exists(blur_path):
                                    subject_missing.append(blur_path)
                                    missing_files_count += 1
                                    subject_feature_structures[feature_name]['cortex'] = False
                                    subject_missing_features.add(feature_name)

            if self.hippocampus and self.hippunfold_directory:
                for meta in feature_meta:
                    if not meta['requires_hippocampus']:
                        continue
                    feature_name = meta['original']
                    if subject_feature_structures[feature_name]['hippocampus'] is not True:
                        continue
                    for hemi in ["L", "R"]:
                        hippo_path = os.path.join(
                            hippocampus_dir,
                            f"{participant_id}_{session_id}_hemi-{hemi}_den-{hippocampal_resolution}_label-hipp_"
                            f"midthickness_feature-{meta['hippo_token']}_smooth-{hippocampal_smoothing}mm.func.gii"
                        )
                        all_files_count += 1
                        if not os.path.exists(hippo_path):
                            subject_missing.append(hippo_path)
                            missing_files_count += 1
                            subject_feature_structures[feature_name]['hippocampus'] = False
                            subject_missing_features.add(feature_name)

            if self.subcortical and self.freesurfer_directory:
                # Only check for volume file if thickness or volume is requested as a feature
                check_volume = (
                    subject_structures['subcortical'] is True
                    and any(
                        feat.lower() in ['thickness', 'volume']
                        and source_has_feature_structure(subject, feat, 'subcortical')
                        for feat in self.features
                    )
                )
                
                if check_volume:
                    volume_file = os.path.join(subcortical_dir, f"{bids_id}_feature-volume.csv")
                    all_files_count += 1
                    if not os.path.exists(volume_file):
                        subject_missing.append(volume_file)
                        missing_files_count += 1
                        if subject_structures['subcortical'] is not None:
                            subject_structures['subcortical'] = False

                for meta in feature_meta:
                    if not meta['requires_subcortical']:
                        continue
                    feature_name = meta['original']
                    if subject_feature_structures[feature_name]['subcortical'] is not True:
                        continue
                    subcort_path = os.path.join(
                        subcortical_dir,
                        f"{bids_id}_feature-{meta['subcortical_token']}.csv"
                    )
                    all_files_count += 1
                    if not os.path.exists(subcort_path):
                        subject_missing.append(subcort_path)
                        missing_files_count += 1
                        subject_feature_structures[feature_name]['subcortical'] = False
                        subject_missing_features.add(feature_name)

            if any(flag is True for flag in subject_structures.values() if flag is not None):
                processed_valid_subjects['base'].append((participant_id, session_id))
                for structure_name, flag in subject_structures.items():
                    if flag is True:
                        processed_valid_subjects['structures'][structure_name].append((participant_id, session_id))

            for meta in feature_meta:
                feature_name = meta['original']
                structure_flags = subject_feature_structures[feature_name]
                relevant_flags = [flag for flag in structure_flags.values() if flag is not None]
                if not relevant_flags:
                    continue
                if any(flag is True for flag in relevant_flags):
                    feature_processed[feature_name]['all'].append((participant_id, session_id))
                    for structure_name, flag in structure_flags.items():
                        if flag is True:
                            feature_processed[feature_name]['structures'][structure_name].append((participant_id, session_id))
                else:
                    subject_missing_features.add(feature_name)

            if subject_missing_features:
                missing_features[(participant_id, session_id)] = subject_missing_features

            if subject_missing:
                missing_files[(participant_id, session_id)] = subject_missing
            else:
                subjects_with_complete_data.append((participant_id, session_id))

        for feature_name, data in feature_processed.items():
            data['all'] = list(dict.fromkeys(data['all']))
            for structure_name in data['structures']:
                data['structures'][structure_name] = list(dict.fromkeys(data['structures'][structure_name]))
            feature_availability[feature_name] = len(data['all'])

        for structure_name in processed_valid_subjects['structures']:
            processed_valid_subjects['structures'][structure_name] = list(dict.fromkeys(processed_valid_subjects['structures'][structure_name]))

        processed_valid_subjects.update(feature_processed)

        # Validation may target a shared, already-populated subject-normalized
        # base. Reapply the logical mask so excluded files on disk are not admitted
        # back into the normative reference merely because they physically exist.
        self._apply_control_feature_exclusions(processed_valid_subjects)
        for feature_name in feature_processed:
            feature_availability[feature_name] = len(
                processed_valid_subjects[feature_name]["all"]
            )

        self.valid_subjects = processed_valid_subjects
        self.subjects_with_complete_data = subjects_with_complete_data
        self.missing_files = missing_files
        self.missing_features = missing_features
        self.feature_availability = feature_availability

        total_subjects = len(output_subjects)
        complete_subjects = len(subjects_with_complete_data)
        complete_percentage = (complete_subjects / total_subjects) * 100 if total_subjects else 0.0
        file_presence_percentage = ((all_files_count - missing_files_count) / all_files_count) * 100 if all_files_count else 0.0

        summary = {
            "total_subjects": total_subjects,
            "complete_subjects": complete_subjects,
            "complete_percentage": complete_percentage,
            "total_files_expected": all_files_count,
            "missing_files_count": missing_files_count,
            "file_presence_percentage": file_presence_percentage
        }
        print("\nValidation Summary:")
        print(f"  Total subjects (outputs): {total_subjects}")
        print(f"  Subjects with complete data: {complete_subjects} ({complete_percentage:.1f}%)")
        print(f"  Total expected files: {all_files_count}")
        print(f"  Missing files: {missing_files_count} ({100 - file_presence_percentage:.1f}%)")
        if verbose and missing_files:
            print(f"\nSubjects with missing files: {len(missing_files)}")
            for (pid, sid), files in missing_files.items():
                print(f"  {pid}/{sid}: {len(files)} missing")
                for path in files[:10]:
                    print(f"    - {os.path.basename(path)}")
                if len(files) > 10:
                    print(f"    - ... and {len(files) - 10} more")

        feature_subject_count = sum(len(data['all']) for data in feature_processed.values())
        if feature_subject_count == 0:
            self.valid_dataset = False
            raise ValueError("Validation failed: no subjects have any complete processed feature data.")

        if error_threshold > 0 and complete_percentage < error_threshold:
            self.valid_dataset = False
            raise ValueError(f"Validation failed: {complete_percentage:.1f}% complete "
                             f"(threshold {error_threshold}%).")

        if warning_threshold and complete_percentage < warning_threshold:
            import warnings
            warnings.warn(f"Validation warning: only {complete_percentage:.1f}% of subjects have complete data "
                          f"(warning threshold {warning_threshold}%).")

        self.valid_dataset = True
        self.cortical_smoothing = cortical_smoothing
        self.hippocampal_smoothing = hippocampal_smoothing

        if verbose:
            print("\nFeature availability summary:")
            denom = total_subjects if total_subjects else 1
            for feature_name in self.features:
                avail_count = feature_availability.get(feature_name, 0)
                print(f"  {feature_name}: {avail_count}/{total_subjects} subjects ({(avail_count/denom)*100:.1f}%)")
                if feature_name in feature_processed:
                    for structure_name in ['cortex', 'hippocampus', 'subcortical']:
                        struct_list = feature_processed[feature_name]['structures'][structure_name]
                        if struct_list:
                            print(f"    - {structure_name}: {len(struct_list)}/{total_subjects} subjects ({(len(struct_list)/denom)*100:.1f}%)")

            print("\nStructure availability summary:")
            for structure_name in ['cortex', 'hippocampus', 'subcortical']:
                struct_list = processed_valid_subjects['structures'][structure_name]
                if struct_list:
                    print(f"  {structure_name}: {len(struct_list)}/{total_subjects} subjects ({(len(struct_list)/denom)*100:.1f}%)")

        return {
            "valid_subjects": processed_valid_subjects['base'],
            "missing_files": missing_files,
            "summary": summary
        }
    
    def check_output_directories(self, output_directory, only_demographics_subjects=True, verbose=True):
        """
        Check the output directory structure and identify valid subjects with output data.
        This is used during validation to determine which subjects to validate.
        
        Parameters:
        -----------
        output_directory : str
            Directory where processed data should be stored
        only_demographics_subjects : bool, default=True
            If True, only check directories for subjects in the demographics
        verbose : bool, default=True
            If True, prints detailed information about the directories

        Returns:
        --------
        list
            List of tuples containing (participant_id, session_id) for subjects with output directories
        """
        if not os.path.exists(output_directory):
            if verbose:
                print(f"Output directory {output_directory} does not exist. No subjects to validate.")
            return []
        
        # Find all subject directories in the output directory
        subject_dirs = [d for d in os.listdir(output_directory) 
                    if os.path.isdir(os.path.join(output_directory, d)) and d.startswith('sub-')]
        
        if verbose:
            print(f"Found {len(subject_dirs)} subject directories in output directory.")
        
        # If we're only checking demographics subjects, create a list of valid subject IDs
        valid_demo_subjects = []
        if only_demographics_subjects and self.demographics is not None:
            for _, row in self.demographics.data.iterrows():
                valid_demo_subjects.append((row['participant_id'], row['session_id']))
            
            if verbose:
                print(f"Filtering to {len(valid_demo_subjects)} subjects specified in demographics data.")
        
        valid_output_subjects = []
        missing_cortex = []
        missing_hippocampus = []
        missing_subcortical = []
        skipped_non_demo = []
        
        # Check each subject directory for session directories
        for subject_dir in subject_dirs:
            participant_id = subject_dir
            subject_path = os.path.join(output_directory, subject_dir)
            
            # Find session directories
            session_dirs = [d for d in os.listdir(subject_path) 
                        if os.path.isdir(os.path.join(subject_path, d))]
            
            for session_dir in session_dirs:
                session_id = session_dir
                
                # Skip if only checking demographics subjects and this subject isn't in the demographics
                if only_demographics_subjects and (participant_id, session_id) not in valid_demo_subjects:
                    skipped_non_demo.append((participant_id, session_id))
                    continue
                    
                session_path = os.path.join(subject_path, session_dir)
                
                # Check for required output directories
                has_cortex = os.path.exists(os.path.join(session_path, "maps", "cortex"))
                has_hippocampus = os.path.exists(os.path.join(session_path, "maps", "hippocampus"))
                has_subcortical = os.path.exists(os.path.join(session_path, "maps", "subcortical"))
                
                # Determine if this is a valid subject based on enabled components
                is_valid = True
                
                if self.cortex and not has_cortex:
                    missing_cortex.append((participant_id, session_id))
                    is_valid = False
                    
                if self.hippocampus and not has_hippocampus:
                    missing_hippocampus.append((participant_id, session_id))
                    is_valid = False
                    
                if self.subcortical and not has_subcortical:
                    missing_subcortical.append((participant_id, session_id))
                    is_valid = False
                
                # If all required directories exist, add to valid subjects
                if is_valid:
                    valid_output_subjects.append((participant_id, session_id))
        
        # Report missing directories
        if verbose:
            if skipped_non_demo:
                print(f"Skipped {len(skipped_non_demo)} subject directories not found in demographics data.")
                if len(skipped_non_demo) < 6:
                    for p, s in skipped_non_demo:
                        print(f"  - {p}/{s}")
                else:
                    for p, s in skipped_non_demo[:5]:
                        print(f"  - {p}/{s}")
                    print(f"  - ... and {len(skipped_non_demo) - 5} more")
            
            if missing_cortex and self.cortex:
                print(f"Warning: {len(missing_cortex)} subjects missing cortex output directory:")
                for p, s in missing_cortex[:5]:  # Show only first 5 to avoid excessive output
                    print(f"  - {p}/{s}")
                if len(missing_cortex) > 5:
                    print(f"  - ... and {len(missing_cortex) - 5} more")
                    
            if missing_hippocampus and self.hippocampus:
                print(f"Warning: {len(missing_hippocampus)} subjects missing hippocampus output directory:")
                for p, s in missing_hippocampus[:5]:
                    print(f"  - {p}/{s}")
                if len(missing_hippocampus) > 5:
                    print(f"  - ... and {len(missing_hippocampus) - 5} more")
                    
            if missing_subcortical and self.subcortical:
                print(f"Warning: {len(missing_subcortical)} subjects missing subcortical output directory:")
                for p, s in missing_subcortical[:5]:
                    print(f"  - {p}/{s}")
                if len(missing_subcortical) > 5:
                    print(f"  - ... and {len(missing_subcortical) - 5} more")

        return valid_output_subjects
    
    def analyze(
        self,
        reference,
        method='zscore',
        output_directory=None,
        verbose=True,
        use_curvature_covariates=True,
        gyrification_matching=None,
        predictive_wscore=False,
        wscore_distribution="gaussian",
        wscore_preprocessing="none",
        wscore_covariate_model="linear",
        wscore_surface_smoothing_iterations=0,
        blur_depth_model="mean_slope_rms",
        intensity_depth_model="raw",
        quant_sample_surface=None,
        control_correlation_filter=False,
        control_correlation_quantile=None,
        prediction_variance_percentile=None,
        export_control_features=None,
        site_harmonizer=None,
        site=None,
        scoring_site_covariate=False,
    ):
        """
        Analyze the dataset by comparing it to a reference dataset using specified method.
        
        Parameters:
        -----------
        reference : zbdataset
            Reference dataset to compare against (e.g., control dataset)
        method : str, default='zscore'
            Analysis method to use. Currently supports 'zscore'and 'wscore'
        output_directory : str, optional
            Directory where z-score maps will be saved. If None, uses validation output directory
        verbose : bool, default=True
            If True, prints detailed information about the analysis process
        use_curvature_covariates : bool, default=True
            If True, cortical W-score models include vertex-wise micapipe
            curvature as an additional covariate, and cortical Z-score maps
            residualize each vertex against the same curvature covariates.
        gyrification_matching : bool, optional
            Deprecated alias for use_curvature_covariates.
        predictive_wscore : bool, default=False
            If True, W-score maps include patient prediction uncertainty in
            the denominator.
        wscore_distribution : str, default="gaussian"
            Distribution fitted to normative-model residuals.
        wscore_preprocessing : {"none", "spatial_zscore", "spatial_robust_z"}, default="none"
            Optional within-subject robust spatial normalization before the
            normative regression.
        wscore_covariate_model : str, default="linear"
            Optional quadratic-age and/or age-by-sex demographic design.
        wscore_surface_smoothing_iterations : int, default=0
            Post-score one-ring smoothing steps for fsLR-32k cortex maps.
        blur_depth_model : {"mean_slope_rms", "mean", "gradient_flattening"}, default="mean_slope_rms"
            Multi-depth blur reduction used for cortical W-score maps.
        intensity_depth_model : {"raw", "white_swm1_direction_cosine", "multisurface_median_abs_dominant"}, default="raw"
            Optional local depth normalization for fsLR-32k white-surface
            T1w and FLAIR intensity abnormalities.
        control_correlation_filter : bool, default=False
            If True, each vertex-map feature drops the bottom
            ``control_correlation_quantile`` fraction of controls by their average
            correlation to the other controls for that feature (rank-based).
        control_correlation_quantile : float or dict or None, default=None
            Drop-fraction for the rank-based control filter, applied independently
            per feature/map. A scalar uses the same fraction for every feature; a
            mapping ``feature -> fraction`` sets a per-feature fraction (a feature
            absent from the mapping, or mapped to ``None``, is not filtered).

        Returns:
        --------
        dict
            Dictionary containing analysis results for each feature and region

        """
        
        if gyrification_matching is not None:
            use_curvature_covariates = bool(gyrification_matching)

        # Call the analysis function and store results
        results = analyze_dataset(
            self,
            reference,
            method,
            output_directory,
            verbose,
            use_curvature_covariates=use_curvature_covariates,
            predictive_wscore=predictive_wscore,
            wscore_distribution=wscore_distribution,
            wscore_preprocessing=wscore_preprocessing,
            wscore_covariate_model=wscore_covariate_model,
            wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
            blur_depth_model=blur_depth_model,
            intensity_depth_model=intensity_depth_model,
            quant_sample_surface=quant_sample_surface,
            control_correlation_filter=control_correlation_filter,
            control_correlation_quantile=control_correlation_quantile,
            prediction_variance_percentile=prediction_variance_percentile,
            export_control_features=export_control_features,
            site_harmonizer=site_harmonizer,
            site=site,
            scoring_site_covariate=scoring_site_covariate,
        )
        self.analysis_results = results
        self.reference_demographics = reference.demographics
        return results
    
    def clinical_report(self, output_directory=None, approach='wscore', analyses=['regional','asymmetry'], features=None, 
                    threshold=1.96, threshold_alpha=0.3, color_range=(-3, 3), 
                    cmap='cmo.balance', cmap_asymmetry='cmo.balance_r', 
                    color_bar='bottom', tmp_dir=None, env=None, verbose=True, normalization=None):
        """
        Generate clinical reports for each subject in the dataset.
        
        Parameters:
        -----------
        output_directory : str, optional
            Directory where processed data is stored. If None, uses the directory from analysis
        approach : str, default='wscore'
            Analysis approach used ('zscore' or 'wscore')
        analyses : list, optional
            List of analyses to include in report. If None, uses ['regional']
        features : list, optional
            List of features to include in report. If None, uses dataset features
        threshold : float, default=1.96
            Threshold for significance in maps
        threshold_alpha : float, default=0.3
            Alpha transparency for thresholded regions
        color_range : tuple, default=(-3, 3)
            Color range for visualization
        cmap : str, default='cmo.balance'
            Colormap for regional analysis
        cmap_asymmetry : str, default='cmo.balance_r'
            Colormap for asymmetry analysis
        color_bar : str, default='bottom'
            Position of color bar
        tmp_dir : str, optional
            Base temporary directory. If None, uses session directory
        verbose : bool, default=True
            If True, prints detailed information
        normalization : {"ravel", "whitestripe"}, optional
            Accepted for compatibility with process() settings dictionaries.
            
        Returns:
        --------
        list
            List of generated PDF report file paths
        """

        if env == None:
            raise ValueError("env must be specified to access workbench and other paths")

        if normalization is not None:
            # Accept composite labels (e.g. noneRavel/whitestripeNyul) as well as
            # the legacy single modes; store the canonical composite desc so
            # downstream path resolution reads the right desc-<composite> volumes.
            subject_norm, dataset_norm = decompose_normalization_label(normalization)
            self.intensity_normalization = compose_normalization_label(subject_norm, dataset_norm)

        # Check if analysis has been run
        if not hasattr(self, 'analysis_results'):
            raise ValueError("No analysis results found. Please run dataset.analyze() first.")
            
        # Use features from the dataset if not specified
        if features is None:
            features = list(self.features)

        feature_mapping = {
            "thickness": "thickness",
            "flair": "FLAIR",
            "sa": "SA",
            "adc": "ADC",
            "fa": "FA",
            "qt1": "qT1",
            "qt1*blur": "qT1*blur",
            "flair*blur": "FLAIR*blur",
            "fmri": "fMRI",
            "t1w": "T1w",
            "t1w*blur": "T1w*blur",
            "adc*blur": "ADC*blur",
            "fa*blur": "FA*blur",
        }

        features = [feature_mapping.get(f.lower(), f) for f in features]

        # Default analyses
        if analyses is None:
            analyses = ['regional']  # Default to regional analysis
            
        # Determine hippocampal resolution for report
        res_hip_param = "8k" if self.hippocampus and self.hippunfold_version == 2 else "high"

        # Validate output directory
        if output_directory is None:
            raise ValueError("output_directory must be specified")
            
        # Don't create a separate reports directory - save directly in subject folders
        generated_reports = []
        
        # Extract subjects from analysis results
        valid_subjects = []
        for region_type, region_results in self.analysis_results.items():
            for feature, feature_results in region_results.items():
                if region_type == "subcortical":
                    if f'patient_{approach}s' in feature_results:
                        for result in feature_results[f'patient_{approach}s']:
                            if 'subject' in result:
                                valid_subjects.append(result['subject'])
                else:
                    for map_key, map_results in feature_results.items():
                        if f'patient_{approach}s' in map_results:
                            for result in map_results[f'patient_{approach}s']:
                                if 'subject' in result:
                                    valid_subjects.append(result['subject'])
        
        # Remove duplicates while preserving order
        valid_subjects = list(dict.fromkeys(valid_subjects))
        
        if not valid_subjects:
            if verbose:
                print("No valid subjects found in analysis results. Using all subjects from dataset.")
            valid_subjects = self.valid_subjects['base']
        
        if verbose:
            print(f"Generating clinical reports for {len(valid_subjects)} subjects...")
            
        # Generate report for each subject
        for participant_id, session_id in valid_subjects:
            try:
                # Get demographics for this subject
                subject_demo = self.demographics.data[
                    (self.demographics.data['participant_id'] == participant_id) &
                    (self.demographics.data['session_id'] == session_id)
                ]
                
                if subject_demo.empty:
                    if verbose:
                        print(f"Warning: No demographics found for {participant_id}/{session_id}, skipping...")
                    continue
                    
                # Extract demographics
                age = subject_demo['AGE'].iloc[0] if 'AGE' in subject_demo.columns else None
                sex = subject_demo['SEX'].iloc[0] if 'SEX' in subject_demo.columns else None
                
                # Make sure sex is an int ONLY if it is safely convertible
                if sex is not None:
                    try:
                        sex = int(float(sex))
                    except (ValueError, TypeError):
                        pass # keep sex as string like 'F', 'M'
                        
                # Convert binary sex encoding back to string if needed
                if sex is not None and isinstance(sex, (int, float)):
                    if hasattr(self.reference_demographics, 'binary_encodings'):
                        # Find the original value that maps to this binary code
                        print("ref binary encodings:", self.reference_demographics.binary_encodings)
                        encoding = self.reference_demographics.binary_encodings['SEX']
                        sex = next((k for k, v in encoding.items() if v == sex), sex)
                    else:
                        raise ValueError("Reference demographics must have binary_encodings for 'SEX' to decode numeric values.")
                
                # Subject directory path
                subject_dir = os.path.join(output_directory, participant_id, session_id)
                
                # Create session-specific temporary directory
                session_tmp_dir = os.path.join(subject_dir, "tmp_clinical_report")
                os.makedirs(session_tmp_dir, exist_ok=True)
                
                # Generate report - save directly in subject directory, not in separate reports folder
                report_path = generate_clinical_report(
                    sid=participant_id,
                    ses=session_id,
                    age=age,
                    sex=sex,
                    analyses=analyses,
                    features=features,
                    approach=approach,
                    threshold=threshold,
                    threshold_alpha=threshold_alpha,
                    color_range=color_range,
                    cmap=cmap,
                    cmap_asymmetry=cmap_asymmetry,
                    color_bar=color_bar,
                    tmp_dir=session_tmp_dir,  # Use session-specific tmp directory
                    subject_dir=subject_dir,
                    output_dir=None,  # Don't use separate output directory
                    tag=f"{participant_id}_{session_id}_{approach}_clinical_report",
                    smooth_ctx=self.cortical_smoothing,
                    smooth_hip=self.hippocampal_smoothing,
                    res_hip=res_hip_param,
                    env=env,
                    verbose=verbose,
                    control_data=self.reference_demographics,
                    valid_subjects=self.valid_subjects,
                    hippunfold_version=self.hippunfold_version,
                )
                
                generated_reports.append(report_path)
                
                if verbose:
                    print(f"Generated report for {participant_id}/{session_id}: {report_path}")
                
                # Clean up temporary directory after successful report generation
                import shutil
                if os.path.exists(session_tmp_dir):
                    shutil.rmtree(session_tmp_dir)
                    if verbose:
                        print(f"  Cleaned up temporary directory: {session_tmp_dir}")
                    
            except Exception as e:
                if verbose:
                    print(f"Error generating report for {participant_id}/{session_id}: {e}")
                    import traceback
                    traceback.print_exc()
                # Clean up temporary directory even on error
                import shutil
                session_tmp_dir = os.path.join(output_directory, participant_id, session_id, "tmp_clinical_report")
                if os.path.exists(session_tmp_dir):
                    shutil.rmtree(session_tmp_dir)
                continue
        
        if verbose:
            print(f"Generated {len(generated_reports)} clinical reports in subject directories")
            
        return generated_reports
    
    def _copy_hippocampal_structural_files(
        self,
        participant_id,
        session_id,
        output_directory,
        verbose=False,
    ):
        """Copy only the HippUnfold surfaces needed downstream.

        Used by selective V1 -> V2 migration after the non-hippocampal structural
        files have been linked from the matching legacy base.
        """
        import shutil

        structural_output_dir = os.path.join(
            output_directory,
            participant_id,
            session_id,
            "structural",
        )
        os.makedirs(structural_output_dir, exist_ok=True)
        if self.hippunfold_version == 2:
            hipp_base = self.hippunfold_directory
            den_tag = "8k"
        else:
            hipp_base = os.path.join(self.hippunfold_directory, "hippunfold")
            den_tag = "0p5mm"

        success = True
        for hemi in ("L", "R"):
            for space in ("unfold", "T1w"):
                filename = (
                    f"{participant_id}_{session_id}_hemi-{hemi}_space-{space}_"
                    f"den-{den_tag}_label-hipp_midthickness.surf.gii"
                )
                source_file = os.path.join(
                    hipp_base,
                    participant_id,
                    session_id,
                    "surf",
                    filename,
                )
                target_file = os.path.join(structural_output_dir, filename)
                if not os.path.isfile(source_file):
                    success = False
                    if verbose:
                        print(f"  Warning: structural file not found: {source_file}")
                    continue
                try:
                    shutil.copy2(source_file, target_file)
                    if verbose:
                        print(f"  Copying hippocampal structural file: {filename}")
                except Exception as error:
                    success = False
                    if verbose:
                        print(f"  Error copying {source_file}: {error}")
        return success

    def _copy_structural_files(self, participant_id, session_id, output_directory, verbose=False):
        """
        Copy structural surface files from micapipe directory to output directory.
        
        Parameters:
        -----------
        participant_id : str
            Subject ID
        session_id : str
            Session ID
        output_directory : str
            Base output directory
        verbose : bool
            Whether to print verbose messages
        
        Returns:
        --------
        bool
            Whether the copying was successful
        """
        import shutil
        
        # Define source and target paths
        structural_output_dir = os.path.join(output_directory, participant_id, session_id, "structural")
        os.makedirs(structural_output_dir, exist_ok=True)
        
        # Files to copy
        surface_files = [
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_surf-fsnative_label-sphere.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_surf-fsnative_label-sphere.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-pial.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-pial.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"anat/{participant_id}_{session_id}_space-nativepro_T1w.nii.gz"),
        ]
        
        if self.hippocampus:
            # Determine correct path structure based on version
            if self.hippunfold_version == 2:
                # Flat structure: /hippunfold_dir/sub-X/ses-Y/
                hipp_base = self.hippunfold_directory
                den_tag = "8k"
            else:
                # Nested structure: /hippunfold_dir/hippunfold/sub-X/ses-Y/
                hipp_base = os.path.join(self.hippunfold_directory, "hippunfold")
                den_tag = "0p5mm"

            # Left hemisphere space-unfold
            surface_files.append(os.path.join(hipp_base, participant_id, session_id, 
                f"surf/{participant_id}_{session_id}_hemi-L_space-unfold_den-{den_tag}_label-hipp_midthickness.surf.gii"))

            # Right hemisphere space-unfold
            surface_files.append(os.path.join(hipp_base, participant_id, session_id, 
                f"surf/{participant_id}_{session_id}_hemi-R_space-unfold_den-{den_tag}_label-hipp_midthickness.surf.gii"))

            # Left hemisphere space-T1w
            surface_files.append(os.path.join(hipp_base, participant_id, session_id, 
                f"surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-{den_tag}_label-hipp_midthickness.surf.gii"))

            # Right hemisphere space-T1w
            surface_files.append(os.path.join(hipp_base, participant_id, session_id, 
                f"surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-{den_tag}_label-hipp_midthickness.surf.gii"))

        success = True
        for file_path in surface_files:
            source_file = os.path.join(file_path)
            target_file = os.path.join(structural_output_dir, os.path.basename(file_path))
            
            if os.path.exists(source_file):
                if verbose:
                    print(f"  Copying structural file: {os.path.basename(source_file)}")
                try:
                    shutil.copy2(source_file, target_file)
                except Exception as e:
                    if verbose:
                        print(f"  Error copying {source_file}: {e}")
                    success = False
            else:
                if verbose:
                    print(f"  Warning: Structural file not found: {source_file}")
                success = False
                
        return success
    
    def _create_midline_from_freesurfer(self, participant_id, session_id, output_directory, verbose=False):
        # Pick the surface you consider "fsnative" (coords/triangles don’t actually matter for the mask)
        lh_surf = os.path.join(self.freesurfer_directory,f"{participant_id}_{session_id}/surf/lh.white")
        rh_surf = os.path.join(self.freesurfer_directory,f"{participant_id}_{session_id}/surf/rh.white")
        lh_coords, _ = nib.freesurfer.read_geometry(str(lh_surf))
        rh_coords, _ = nib.freesurfer.read_geometry(str(rh_surf))
        n_lh, n_rh = lh_coords.shape[0], rh_coords.shape[0]

        # Read cortex labels (indices of vertices that are true cortex)
        lh_cortex_idx = nib.freesurfer.read_label(os.path.join(self.freesurfer_directory,f"{participant_id}_{session_id}", "label", "lh.cortex.label"))
        rh_cortex_idx = nib.freesurfer.read_label(os.path.join(self.freesurfer_directory,f"{participant_id}_{session_id}", "label", "rh.cortex.label"))

        # Boolean masks (True = cortex, False = medial wall). Sized to the fsnative
        # white surface; out-of-range label indices are dropped (with a warning)
        # instead of crashing the whole subject -- see _cortex_bool_mask.
        _who = f"{participant_id}/{session_id}"
        lh_is_cortex = _cortex_bool_mask(lh_cortex_idx, n_lh, hemi="L", label=_who)
        rh_is_cortex = _cortex_bool_mask(rh_cortex_idx, n_rh, hemi="R", label=_who)

        # Medial-wall masks (True = medial wall)
        lh_medial_wall = ~lh_is_cortex
        rh_medial_wall = ~rh_is_cortex

        structural_output_dir = os.path.join(output_directory, participant_id, session_id, "structural")
        os.makedirs(structural_output_dir, exist_ok=True)

        lh_medial_wall_file = os.path.join(structural_output_dir, f"{participant_id}_{session_id}_hemi-L_surf-fsnative_label-medialwall.label.gii")
        rh_medial_wall_file = os.path.join(structural_output_dir, f"{participant_id}_{session_id}_hemi-R_surf-fsnative_label-medialwall.label.gii")

        # Save medial wall labels as gifti files
        lh_medial_wall_gifti = nib.gifti.GiftiImage()
        rh_medial_wall_gifti = nib.gifti.GiftiImage()
        lh_medial_wall_gifti.add_gifti_data_array(nib.gifti.GiftiDataArray(data=lh_medial_wall.astype(np.uint8), intent="NIFTI_INTENT_LABEL"))
        rh_medial_wall_gifti.add_gifti_data_array(nib.gifti.GiftiDataArray(data=rh_medial_wall.astype(np.uint8), intent="NIFTI_INTENT_LABEL"))
        nib.save(lh_medial_wall_gifti, lh_medial_wall_file)
        nib.save(rh_medial_wall_gifti, rh_medial_wall_file)
        if verbose:
            print(f"  Created medial wall labels for {participant_id}/{session_id} at {structural_output_dir}")
