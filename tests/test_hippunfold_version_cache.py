"""HippUnfold v1/v2 inputs must never collide in optimization caches."""

import types

import pandas as pd
import pytest

import zbrains_benchmark as zb
import zbrains_staged as st
from results import greedy_auroc as ga
from zbrains.dataset import _norm_cache_dir_for, zbdataset
from zbrains.hippunfold import (
    hippunfold_cache_tag,
    infer_hippunfold_major,
    without_hippunfold_cache_tag,
)


def _cohort(path, version):
    dataset = types.SimpleNamespace(
        hippocampus=True,
        hippunfold_directory=str(path),
        requested_hippunfold_version=version,
        _hippunfold_version_detected=True,
        hippunfold_version=version,
    )
    return types.SimpleNamespace(
        name="NOEL",
        control_dataset=dataset,
        patient_dataset=dataset,
        output_dir_prefix="/derivatives",
        base_pipeline_settings={"features": ["thickness"]},
    )


def test_source_tags_and_processed_bases_differ_between_v1_and_v2(tmp_path):
    v1_path = tmp_path / "hippunfold_v1.4.1"
    v2_path = tmp_path / "hippunfold_v2.0.0"
    v1_path.mkdir()
    v2_path.mkdir()
    assert infer_hippunfold_major(v1_path) == 1
    assert infer_hippunfold_major(v2_path) == 2
    assert hippunfold_cache_tag(v1_path, 1) != hippunfold_cache_tag(v2_path, 2)

    config = dict(st.SIMPLEST_BASELINE)
    v1_signature = st._cohort_base_signature(_cohort(v1_path, 1), config)
    v2_signature = st._cohort_base_signature(_cohort(v2_path, 2), config)
    assert "hippunfoldv1-" in v1_signature
    assert "hippunfoldv2-" in v2_signature
    assert v1_signature != v2_signature
    assert (
        zb.processed_base_directory_for("none", "/derivatives", v1_signature)
        != zb.processed_base_directory_for("none", "/derivatives", v2_signature)
    )


def test_mics_reuse_all_mode_keeps_existing_untagged_paths(tmp_path):
    cohort = _cohort(tmp_path / "hippunfold_v2.0.0", 2)
    cohort.hippunfold_reuse_mode = "all"
    signature = st._cohort_base_signature(cohort, dict(st.SIMPLEST_BASELINE))
    assert "hippunfold" not in signature


def test_nonhippocampal_seed_excludes_all_v1_hippunfold_products(tmp_path):
    pid, sid = "sub-01", "ses-01"
    source = tmp_path / "zbrains_SELF_smfwhm_smoothctx0hip0"
    target = tmp_path / (
        "zbrains_SELF_smfwhm_smoothctx0hip0_hippunfoldv2-aaaaaaaaaa"
    )
    source_maps = source / pid / sid / "maps"
    source_structural = source / pid / sid / "structural"
    for directory in (
        source_maps / "cortex",
        source_maps / "subcortical",
        source_maps / "hippocampus",
        source_structural,
    ):
        directory.mkdir(parents=True, exist_ok=True)
    (source_maps / "cortex" / "cortex.func.gii").touch()
    (source_maps / "subcortical" / "volume.csv").touch()
    (source_maps / "hippocampus" / "v1.func.gii").touch()
    normalized = source_maps / f"{pid}_{sid}_space-nativepro_desc-none_T1w.nii.gz"
    normalized.touch()
    cortical_surface = source_structural / "cortical.surf.gii"
    cortical_surface.touch()
    (source_structural / f"{pid}_{sid}_den-0p5mm_label-hipp_old.surf.gii").touch()

    dataset = types.SimpleNamespace(
        valid_subjects={"base": [(pid, sid)]},
        cortex=True,
        hippocampus=True,
        subcortical=True,
        hippunfold_directory=str(tmp_path / "hippunfold_v2.0.0"),
        requested_hippunfold_version=2,
        hippunfold_version=2,
    )
    assert zb.seed_non_hippocampal_processed_outputs(source, target, [dataset])

    target_maps = target / pid / sid / "maps"
    assert (target_maps / "cortex").resolve() == (source_maps / "cortex").resolve()
    assert (target_maps / "subcortical").resolve() == (source_maps / "subcortical").resolve()
    assert (target_maps / normalized.name).resolve() == normalized.resolve()
    assert not (target_maps / "hippocampus").exists()
    target_structural = target / pid / sid / "structural"
    assert (target_structural / cortical_surface.name).resolve() == cortical_surface.resolve()
    assert not list(target_structural.glob("*_label-hipp_*"))


def test_nonhippocampal_analysis_seed_excludes_hippocampal_scores(tmp_path):
    pid, sid = "sub-01", "ses-01"
    source = tmp_path / "zbrains_SELF_analysis"
    target = tmp_path / "zbrains_SELF_analysis_hippunfoldv2-aaaaaaaaaa"
    source_scores = source / pid / sid / "zscore_maps"
    for structure in ("cortex", "hippocampus", "subcortical"):
        directory = source_scores / structure
        directory.mkdir(parents=True, exist_ok=True)
        (directory / f"{structure}.func.gii").touch()

    dataset = types.SimpleNamespace(
        valid_subjects={"base": [(pid, sid)]},
        cortex=True,
        hippocampus=True,
        subcortical=True,
    )
    assert zb.seed_non_hippocampal_analysis_outputs(
        source,
        target,
        [dataset],
        "zscore",
    )

    target_scores = target / pid / sid / "zscore_maps"
    assert (target_scores / "cortex").resolve() == (source_scores / "cortex").resolve()
    assert (target_scores / "subcortical").resolve() == (
        source_scores / "subcortical"
    ).resolve()
    assert not (target_scores / "hippocampus").exists()
    assert zb.non_hippocampal_analysis_outputs_complete(
        target,
        [dataset],
        method="zscore",
    )
    assert not zb.analysis_structure_maps_complete(
        target,
        [dataset],
        method="zscore",
    )


def test_legacy_base_path_is_recovered_from_v2_tag():
    tagged = "/derivatives/zbrains_SELF_smfwhm_hippunfoldv2-1234567890_exclabc/"
    assert without_hippunfold_cache_tag(tagged) == (
        "/derivatives/zbrains_SELF_smfwhm_exclabc/"
    )


def test_staged_noel_migration_requests_hippocampus_only(monkeypatch, tmp_path):
    calls = []
    processed = set()

    def make_dataset(name, pid):
        dataset = types.SimpleNamespace(
            name=name,
            demographics=types.SimpleNamespace(
                data=pd.DataFrame(
                    {"participant_id": [pid], "session_id": ["ses-01"]}
                )
            ),
            valid_subjects={"base": [(pid, "ses-01")]},
            cortex=True,
            hippocampus=True,
            subcortical=True,
            hippunfold_directory=str(tmp_path / "hippunfold_v2.0.0"),
            requested_hippunfold_version=2,
            _hippunfold_version_detected=True,
            hippunfold_version=2,
        )

        def process(**kwargs):
            calls.append((name, kwargs))
            processed.add(name)

        dataset.process = process
        return dataset

    control = make_dataset("controls", "sub-HC001")
    patient = make_dataset("patients", "sub-PX001")
    cohort = types.SimpleNamespace(
        name="NOEL",
        score_name="NOEL",
        control_dataset=control,
        patient_dataset=patient,
        output_dir_prefix=str(tmp_path),
        base_pipeline_settings={"features": ["thickness"]},
        hippunfold_reuse_mode="nonhippocampal",
    )
    monkeypatch.setattr(st, "_seed_cohort_nonhippocampal_base", lambda *a: True)
    monkeypatch.setattr(zb, "base_is_marked_complete", lambda *a: False)
    monkeypatch.setattr(zb, "processed_maps_complete", lambda *a: True)
    monkeypatch.setattr(
        zb,
        "subjects_missing_non_hippocampal_processed_outputs",
        lambda *a: [],
    )
    monkeypatch.setattr(
        zb,
        "hippocampal_v2_outputs_complete",
        lambda _base, datasets: all(dataset.name in processed for dataset in datasets),
    )
    monkeypatch.setattr(zb, "mark_base_complete", lambda *a: None)
    monkeypatch.setattr(zb, "unmark_base_complete", lambda *a: None)

    st._ensure_base_processed(cohort, dict(st.SIMPLEST_BASELINE), env=object())

    assert [name for name, _ in calls] == ["controls", "patients"]
    assert all(kwargs["hippocampus_only"] is True for _, kwargs in calls)


def test_partial_noel_seed_rebuilds_only_missing_nonhippocampal_subject(
    monkeypatch,
    tmp_path,
):
    calls = []

    def make_dataset(pairs):
        dataset = types.SimpleNamespace(
            name="patients",
            demographics=types.SimpleNamespace(
                data=pd.DataFrame(
                    pairs,
                    columns=["participant_id", "session_id"],
                )
            ),
            valid_subjects={"base": list(pairs)},
            cortex=True,
            hippocampus=True,
            subcortical=True,
        )

        def process(**kwargs):
            calls.append((tuple(pairs), kwargs))

        dataset.process = process
        return dataset

    all_pairs = [("sub-01", "ses-01"), ("sub-02", "ses-01")]
    dataset = make_dataset(all_pairs)
    monkeypatch.setattr(
        zb,
        "subjects_missing_non_hippocampal_processed_outputs",
        lambda *a: [("sub-02", "ses-01")],
    )
    monkeypatch.setattr(zb, "hippocampal_v2_outputs_complete", lambda *a: False)
    monkeypatch.setattr(
        st,
        "_subset_control_dataset",
        lambda _dataset, pairs, name=None: make_dataset(list(pairs)),
    )

    st._process_dataset_for_hippunfold_migration(
        dataset,
        tmp_path,
        env=object(),
        settings={"features": ["thickness"]},
    )

    assert calls[0][0] == (("sub-01", "ses-01"),)
    assert calls[0][1]["hippocampus_only"] is True
    assert calls[1][0] == (("sub-02", "ses-01"),)
    assert "hippocampus_only" not in calls[1][1]


def test_v2_completeness_requires_every_available_hippocampal_feature(tmp_path):
    pid, sid = "sub-01", "ses-01"
    subject = (pid, sid)
    dataset = types.SimpleNamespace(
        valid_subjects={"base": [subject]},
        source_valid_subjects={
            feature: {"structures": {"hippocampus": [subject]}}
            for feature in ("T1w", "FA", "thickness")
        },
        features=["T1w", "FA", "thickness"],
    )
    maps = tmp_path / pid / sid / "maps" / "hippocampus"
    structural = tmp_path / pid / sid / "structural"
    maps.mkdir(parents=True)
    structural.mkdir(parents=True)
    for hemi in ("L", "R"):
        for feature in ("T1w", "thickness"):
            (maps / (
                f"{pid}_{sid}_hemi-{hemi}_den-8k_label-hipp_midthickness_"
                f"feature-{feature}_smooth-0mm.func.gii"
            )).touch()
        for space in ("T1w", "unfold"):
            (structural / (
                f"{pid}_{sid}_hemi-{hemi}_space-{space}_den-8k_"
                "label-hipp_midthickness.surf.gii"
            )).touch()

    assert not zb.hippocampal_v2_outputs_complete(tmp_path, [dataset])
    for hemi in ("L", "R"):
        (maps / (
            f"{pid}_{sid}_hemi-{hemi}_den-8k_label-hipp_midthickness_"
            "feature-FA_smooth-0mm.func.gii"
        )).touch()
    assert zb.hippocampal_v2_outputs_complete(tmp_path, [dataset])


def test_analysis_provenance_changes_with_hippunfold_source(tmp_path):
    v1 = _cohort(tmp_path / "hippunfold_v1.4.1", 1)
    v2 = _cohort(tmp_path / "hippunfold_v2.0.0", 2)
    provenance_v1 = st._analysis_provenance({"method": "zscore"}, v1)
    provenance_v2 = st._analysis_provenance({"method": "zscore"}, v2)
    assert provenance_v1["hippunfold"] != provenance_v2["hippunfold"]


def test_standalone_benchmark_rejects_mixed_hippunfold_sources(tmp_path):
    v1 = _cohort(tmp_path / "hippunfold_v1.4.1", 1).control_dataset
    v2 = _cohort(tmp_path / "hippunfold_v2.0.0", 2).control_dataset

    assert "hippunfoldv2-" in zb.hippunfold_signature_for_datasets(v2, v2)
    with pytest.raises(ValueError, match="different HippUnfold"):
        zb.hippunfold_signature_for_datasets(v1, v2)


def test_normalized_volume_cache_remains_reusable_across_hippunfold_versions(tmp_path):
    shared = "smfwhm_smoothctx0hip0"
    v1 = tmp_path / f"zbrains_NONE_{shared}_hippunfoldv1-aaaaaaaaaa"
    v2 = tmp_path / f"zbrains_NONE_{shared}_hippunfoldv2-bbbbbbbbbb"
    assert _norm_cache_dir_for(v1) == _norm_cache_dir_for(v2)


def test_noel_objective_uses_only_v2_density_cavity_labels(tmp_path):
    stem = "sub-PX001_ses-01_hemi-L_den-{}_label-hipp_cavity.func.gii"
    v1 = tmp_path / stem.format("0p5mm")
    v2 = tmp_path / stem.format("8k")
    v1.touch()
    v2.touch()

    assert ga.dataset_hipp_density("NOEL") == "8k"
    cavities = ga._LA.discover_hipp_cavities(tmp_path, density="8k")
    assert cavities == {("sub-PX001", "ses-01", "L"): v2}


def test_explicit_v2_never_falls_back_to_existing_v1_surface(tmp_path):
    participant, session = "sub-01", "ses-01"
    micapipe = tmp_path / "micapipe"
    (micapipe / participant / session).mkdir(parents=True)
    hippunfold = tmp_path / "hippunfold_v1.4.1"
    v1_surf = hippunfold / "hippunfold" / participant / session / "surf"
    v1_surf.mkdir(parents=True)
    (
        v1_surf
        / f"{participant}_{session}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii"
    ).write_text("v1", encoding="utf-8")
    demographics = types.SimpleNamespace(
        data=pd.DataFrame(
            {"participant_id": [participant], "session_id": [session]}
        )
    )

    dataset = zbdataset(
        "controls",
        demographics,
        str(micapipe),
        hippunfold_directory=str(hippunfold),
        hippunfold_version=2,
        cortex=False,
        hippocampus=True,
        subcortical=False,
    )
    dataset.check_directories()

    assert dataset.hippunfold_version == 2
    assert dataset.valid_subjects["base"] == []


def test_old_optimization_histories_are_not_resume_compatible():
    assert st.HISTORY_SCHEMA_VERSION >= 6
