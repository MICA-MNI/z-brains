import json
import inspect
import os
import sys
import tempfile
import types
import unittest
from unittest.mock import MagicMock, patch

import nibabel as nib
import numpy as np
import pandas as pd

sys.modules.setdefault("skfmm", MagicMock())
clinical_reports = types.ModuleType("zbrains.clinical_reports")
clinical_reports.generate_clinical_report = MagicMock()
sys.modules.setdefault("zbrains.clinical_reports", clinical_reports)

from zbrains.dataset import zbdataset
from zbrains.normalization import (
    get_normalized_modality_path,
    normalize_normalization_mode,
    requested_ravel_modalities,
    whitestripe_normalize_image,
)


class TestWhiteStripeNormalization(unittest.TestCase):
    def test_whitestripe_normalize_image_uses_synthseg_wm_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            image_path = os.path.join(tmp, "image.nii.gz")
            label_path = os.path.join(tmp, "labels.nii.gz")
            output_path = os.path.join(tmp, "whitestripe.nii.gz")
            stats_path = os.path.join(tmp, "stats.json")

            data = np.zeros((7, 7, 7), dtype=np.float32)
            data[1:6, 1:6, 1:6] = np.linspace(90, 110, 125, dtype=np.float32).reshape(5, 5, 5)
            data[0, 0, 0] = 25

            labels = np.zeros((7, 7, 7), dtype=np.int16)
            labels[1:6, 1:6, 1:6] = 2
            labels[0, 0, 0] = 4

            img = nib.Nifti1Image(data, np.eye(4))
            nib.save(img, image_path)
            nib.save(nib.Nifti1Image(labels, np.eye(4)), label_path)

            whitestripe_normalize_image(
                image_path,
                label_path,
                output_path,
                stats_path=stats_path,
                wm_labels=(2,),
            )

            with open(stats_path, "r", encoding="utf-8") as f:
                stats = json.load(f)
            normalized = np.asarray(nib.load(output_path).dataobj)

            nonzero = data != 0
            expected = (data[nonzero] - stats["mean"]) / stats["std"]
            np.testing.assert_allclose(normalized[nonzero], expected, rtol=1e-5, atol=1e-5)
            self.assertEqual(stats["wm_voxels"], 121)
            self.assertGreater(stats["stripe_voxels"], 0)
            self.assertEqual(stats["scale_method"], "trimmed_wm_iqr")
            self.assertGreater(stats["trimmed_wm_iqr_sd"], stats["std"] * 0.99)

    def test_requested_ravel_modalities_only_for_t1w_and_flair(self):
        self.assertEqual(
            requested_ravel_modalities(["FA", "T1w", "FLAIR*blur", "qT1", "t1w*blur"]),
            ["T1w", "FLAIR"],
        )
        self.assertEqual(requested_ravel_modalities(["FA", "T1w"]), ["T1w"])
        self.assertEqual(requested_ravel_modalities(["FA", "FLAIR"]), ["FLAIR"])
        self.assertEqual(requested_ravel_modalities(["FA", "ADC"]), [])

    def test_get_normalized_modality_path_returns_requested_normalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            subject_output = os.path.join(tmp, "sub-01", "ses-01")
            micapipe = os.path.join(tmp, "micapipe", "sub-01", "ses-01")
            maps = os.path.join(subject_output, "maps")
            anat = os.path.join(micapipe, "anat")
            os.makedirs(maps)
            os.makedirs(anat)

            fallback = os.path.join(anat, "sub-01_ses-01_space-nativepro_T1w.nii.gz")
            whitestripe = os.path.join(maps, "sub-01_ses-01_space-nativepro_desc-whitestripe_T1w.nii.gz")
            ravel = os.path.join(maps, "sub-01_ses-01_space-nativepro_desc-ravel_T1w.nii.gz")

            open(fallback, "w", encoding="utf-8").close()
            open(whitestripe, "w", encoding="utf-8").close()
            self.assertEqual(get_normalized_modality_path(subject_output, micapipe, "T1w"), ravel)
            self.assertEqual(
                get_normalized_modality_path(subject_output, micapipe, "T1w", normalization="whitestripe"),
                whitestripe,
            )
            self.assertEqual(normalize_normalization_mode("white-stripe"), "whitestripe")

    def test_downstream_dataset_calls_accept_normalization_setting(self):
        self.assertIn("normalization", inspect.signature(zbdataset.process).parameters)
        self.assertIn("normalization", inspect.signature(zbdataset.validate).parameters)
        self.assertIn("normalization", inspect.signature(zbdataset.clinical_report).parameters)


class TestDatasetNormalizationOrder(unittest.TestCase):
    def test_control_ravel_fits_before_surface_extraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            demo = MagicMock()
            demo.data = pd.DataFrame(
                {
                    "participant_id": ["sub-01", "sub-02"],
                    "session_id": ["ses-01", "ses-01"],
                }
            )
            dataset = zbdataset(
                "controls",
                demographics=demo,
                micapipe_directory=os.path.join(tmp, "micapipe"),
                cortex=True,
                hippocampus=False,
                subcortical=False,
            )
            dataset.features = ["T1w", "FLAIR"]
            subjects = [("sub-01", "ses-01"), ("sub-02", "ses-01")]
            dataset.valid_subjects = {
                "base": subjects,
                "structures": {"cortex": subjects, "hippocampus": [], "subcortical": []},
                "T1w": {"all": subjects, "structures": {"cortex": subjects, "hippocampus": [], "subcortical": []}},
                "FLAIR": {"all": subjects, "structures": {"cortex": subjects, "hippocampus": [], "subcortical": []}},
            }
            env = MagicMock()
            env.num_threads = 1
            env.connectome_workbench_path = "/usr/bin"

            events = []

            def record_whitestripe(participant_id, session_id, **kwargs):
                events.append(("whitestripe", participant_id, session_id, tuple(kwargs["modalities"])))
                return {}

            def record_ravel(**kwargs):
                events.append(("fit_ravel", tuple(kwargs["subjects"]), tuple(kwargs["modalities"])))

            def record_surface(participant_id, session_id, **kwargs):
                events.append(("surface", participant_id, session_id))

            with patch.object(zbdataset, "add_features"), \
                patch.object(zbdataset, "_copy_structural_files"), \
                patch.object(zbdataset, "_create_midline_from_freesurfer"), \
                patch.object(zbdataset, "validate"), \
                patch("zbrains.dataset.prepare_t1w_flair_whitestripe", side_effect=record_whitestripe), \
                patch("zbrains.dataset.fit_and_apply_ravel_to_controls", side_effect=record_ravel), \
                patch("zbrains.dataset.apply_cortical_processing", side_effect=record_surface):
                dataset.process(
                    output_directory=tmp,
                    features=["T1w", "FLAIR"],
                    env=env,
                    verbose=False,
                )

            self.assertEqual(
                events[:4],
                [
                    ("whitestripe", "sub-01", "ses-01", ("T1w",)),
                    ("whitestripe", "sub-01", "ses-01", ("FLAIR",)),
                    ("whitestripe", "sub-02", "ses-01", ("T1w",)),
                    ("whitestripe", "sub-02", "ses-01", ("FLAIR",)),
                ],
            )
            self.assertEqual(events[4][0], "fit_ravel")
            self.assertEqual(events[5][0], "fit_ravel")
            self.assertEqual([event[0] for event in events[6:]], ["surface", "surface"])

    def test_subject_missing_flair_still_processes_t1w(self):
        with tempfile.TemporaryDirectory() as tmp:
            demo = MagicMock()
            demo.data = pd.DataFrame(
                {
                    "participant_id": ["sub-01", "sub-02"],
                    "session_id": ["ses-01", "ses-01"],
                }
            )
            dataset = zbdataset(
                "controls",
                demographics=demo,
                micapipe_directory=os.path.join(tmp, "micapipe"),
                cortex=True,
                hippocampus=False,
                subcortical=False,
            )
            subjects = [("sub-01", "ses-01"), ("sub-02", "ses-01")]
            t1_only_subject = ("sub-02", "ses-01")
            dataset.features = ["T1w", "FLAIR"]
            dataset.valid_subjects = {
                "base": subjects,
                "structures": {"cortex": subjects, "hippocampus": [], "subcortical": []},
                "T1w": {"all": subjects, "structures": {"cortex": subjects, "hippocampus": [], "subcortical": []}},
                "FLAIR": {"all": [subjects[0]], "structures": {"cortex": [subjects[0]], "hippocampus": [], "subcortical": []}},
            }
            env = MagicMock()
            env.num_threads = 1
            env.connectome_workbench_path = "/usr/bin"

            events = []

            def record_whitestripe(participant_id, session_id, **kwargs):
                events.append(("whitestripe", participant_id, session_id, tuple(kwargs["modalities"])))
                return {}

            def record_ravel(**kwargs):
                events.append(("fit_ravel", tuple(kwargs["subjects"]), tuple(kwargs["modalities"])))

            def record_surface(participant_id, session_id, features, **kwargs):
                events.append(("surface", participant_id, session_id, tuple(features)))

            with patch.object(zbdataset, "add_features"), \
                patch.object(zbdataset, "_copy_structural_files"), \
                patch.object(zbdataset, "_create_midline_from_freesurfer"), \
                patch.object(zbdataset, "validate"), \
                patch("zbrains.dataset.prepare_t1w_flair_whitestripe", side_effect=record_whitestripe), \
                patch("zbrains.dataset.fit_and_apply_ravel_to_controls", side_effect=record_ravel), \
                patch("zbrains.dataset.apply_cortical_processing", side_effect=record_surface):
                dataset.process(
                    output_directory=tmp,
                    features=["T1w", "FLAIR"],
                    env=env,
                    verbose=False,
                )

            self.assertIn(("whitestripe", t1_only_subject[0], t1_only_subject[1], ("T1w",)), events)
            self.assertNotIn(("whitestripe", t1_only_subject[0], t1_only_subject[1], ("FLAIR",)), events)
            self.assertIn(("fit_ravel", tuple(subjects), ("T1w",)), events)
            self.assertIn(("fit_ravel", (subjects[0],), ("FLAIR",)), events)
            self.assertIn(("surface", "sub-02", "ses-01", ("T1w",)), events)

    def test_parallel_subject_loop_processes_all_and_isolates_logs(self):
        # env.num_threads>1 with >1 subject takes the THREADED subject loop. It must
        # process every subject (output-identical to serial) and, crucially, each
        # subject must log to its OWN file -- the thread-local stdout proxy replaces
        # the old shared sys.stdout swap that made the parallel path unusable.
        import threading

        with tempfile.TemporaryDirectory() as tmp:
            subjects = [("sub-01", "ses-01"), ("sub-02", "ses-01"), ("sub-03", "ses-01")]
            demo = MagicMock()
            demo.data = pd.DataFrame(
                {"participant_id": [p for p, _ in subjects],
                 "session_id": [s for _, s in subjects]})
            dataset = zbdataset(
                "controls", demographics=demo,
                micapipe_directory=os.path.join(tmp, "micapipe"),
                cortex=True, hippocampus=False, subcortical=False,
            )
            dataset.features = ["T1w", "FLAIR"]
            dataset.valid_subjects = {
                "base": subjects,
                "structures": {"cortex": subjects, "hippocampus": [], "subcortical": []},
                "T1w": {"all": subjects, "structures": {"cortex": subjects, "hippocampus": [], "subcortical": []}},
                "FLAIR": {"all": subjects, "structures": {"cortex": subjects, "hippocampus": [], "subcortical": []}},
            }
            env = MagicMock()
            env.num_threads = 4                       # -> subject_n_jobs = 3 > 1 -> threaded
            env.connectome_workbench_path = "/usr/bin"

            surfaced = []
            threads_seen = set()

            def record_surface(participant_id, session_id, features, **kwargs):
                threads_seen.add(threading.current_thread().name)
                surfaced.append((participant_id, session_id))

            with patch.object(zbdataset, "add_features"), \
                patch.object(zbdataset, "_copy_structural_files"), \
                patch.object(zbdataset, "_create_midline_from_freesurfer"), \
                patch.object(zbdataset, "validate"), \
                patch("zbrains.dataset.prepare_t1w_flair_whitestripe", side_effect=lambda *a, **k: {}), \
                patch("zbrains.dataset.fit_and_apply_ravel_to_controls", side_effect=lambda **k: None), \
                patch("zbrains.dataset.apply_cortical_processing", side_effect=record_surface):
                pre_stdout = sys.stdout
                dataset.process(output_directory=tmp, features=["T1w", "FLAIR"],
                                env=env, verbose=False)
                # the proxy is installed only inside process(); it must be restored
                self.assertIs(sys.stdout, pre_stdout)

            # every subject was processed exactly once ...
            self.assertEqual(sorted(surfaced), sorted(subjects))
            # ... concurrently (more than one worker thread ran) ...
            self.assertGreater(len(threads_seen), 1)
            # ... and each subject's log is its OWN file with only its own
            # "Processing subject" header (no cross-thread bleed).
            for pid, ses in subjects:
                logs = [f for f in os.listdir(os.path.join(tmp, pid, ses))
                        if f.endswith(".log")]
                self.assertEqual(len(logs), 1)
                text = open(os.path.join(tmp, pid, ses, logs[0])).read()
                self.assertIn(f"Processing subject {pid}/{ses}", text)
                for other, oses in subjects:
                    if other != pid:
                        self.assertNotIn(f"Processing subject {other}/{oses}", text)

    def test_process_accepts_preexisting_output_directory_race(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = os.path.join(tmp, "zbrains_SELF_smfwhm_smoothctx0hip0")
            os.makedirs(output_dir)
            demo = MagicMock()
            demo.data = pd.DataFrame({"participant_id": [], "session_id": []})
            dataset = zbdataset(
                "controls",
                demographics=demo,
                micapipe_directory=os.path.join(tmp, "micapipe"),
                cortex=False,
                hippocampus=False,
                subcortical=False,
            )
            dataset.valid_subjects = {
                "base": [],
                "structures": {"cortex": [], "hippocampus": [], "subcortical": []},
            }
            env = MagicMock()
            env.num_threads = 1
            env.connectome_workbench_path = "/usr/bin"

            real_exists = os.path.exists

            def exists_with_creation_race(path):
                if os.path.abspath(path) == os.path.abspath(output_dir):
                    return False
                return real_exists(path)

            with patch.object(zbdataset, "add_features"), \
                patch.object(zbdataset, "validate"), \
                patch("zbrains.dataset.os.chdir"), \
                patch("zbrains.dataset.os.path.exists", side_effect=exists_with_creation_race):
                dataset.process(
                    output_directory=output_dir,
                    features=["thickness"],
                    env=env,
                    verbose=False,
                )

    def test_whitestripe_mode_skips_ravel_and_processes_surfaces(self):
        with tempfile.TemporaryDirectory() as tmp:
            demo = MagicMock()
            demo.data = pd.DataFrame(
                {
                    "participant_id": ["sub-01"],
                    "session_id": ["ses-01"],
                }
            )
            dataset = zbdataset(
                "patients",
                demographics=demo,
                micapipe_directory=os.path.join(tmp, "micapipe"),
                cortex=True,
                hippocampus=False,
                subcortical=False,
            )
            subject = ("sub-01", "ses-01")
            dataset.features = ["T1w"]
            dataset.valid_subjects = {
                "base": [subject],
                "structures": {"cortex": [subject], "hippocampus": [], "subcortical": []},
                "T1w": {"all": [subject], "structures": {"cortex": [subject], "hippocampus": [], "subcortical": []}},
            }
            env = MagicMock()
            env.num_threads = 1
            env.connectome_workbench_path = "/usr/bin"
            events = []

            def record_whitestripe(participant_id, session_id, **kwargs):
                events.append(("whitestripe", tuple(kwargs["modalities"])))
                return {}

            def record_surface(participant_id, session_id, **kwargs):
                events.append(("surface", kwargs["normalization"]))

            with patch.object(zbdataset, "add_features"), \
                patch.object(zbdataset, "_copy_structural_files"), \
                patch.object(zbdataset, "_create_midline_from_freesurfer"), \
                patch.object(zbdataset, "validate"), \
                patch("zbrains.dataset.prepare_t1w_flair_whitestripe", side_effect=record_whitestripe), \
                patch("zbrains.dataset.fit_and_apply_ravel_to_controls") as fit_ravel, \
                patch("zbrains.dataset.apply_ravel_model_to_subject") as apply_ravel, \
                patch("zbrains.dataset.apply_cortical_processing", side_effect=record_surface):
                dataset.process(
                    output_directory=tmp,
                    features=["T1w"],
                    env=env,
                    verbose=False,
                    normalization="whitestripe",
                )

            fit_ravel.assert_not_called()
            apply_ravel.assert_not_called()
            self.assertIn(("whitestripe", ("T1w",)), events)
            self.assertIn(("surface", "whitestripe"), events)


if __name__ == "__main__":
    unittest.main()
