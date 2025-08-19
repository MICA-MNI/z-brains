#!/usr/bin/env python3
import unittest
import sys
import os
import tempfile
from unittest.mock import patch, MagicMock, call
import argparse

from zbrains.cli import parse_args, main
from zbrains.run import run
from zbrains.dataset import demographics, zbdataset

class TestZBrainsArgumentParsing(unittest.TestCase):
    """Test more comprehensive argument parsing scenarios for Z-Brains CLI."""

    def test_quoted_feature_names(self):
        """Test that quoted feature names with spaces/special chars are handled correctly."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC', 'thickness', 'qT1', 'qT1*Blur',
                                'FLAIR', 'FLAIR*Blur',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer']):
            args = parse_args()
            self.assertIn('qT1*Blur', args.features)
            self.assertIn('FLAIR*Blur', args.features)

    def test_output_directory_default_and_override(self):
        """Test default output directory and its override."""
        # Test default value
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer']):
            args = parse_args()
            self.assertEqual(args.output_directory, './zbrains_output')
        
        # Test override
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--output-dir', '/custom/output/path']):
            args = parse_args()
            self.assertEqual(args.output_directory, '/custom/output/path')

    def test_smoothing_parameters(self):
        """Test cortical and hippocampal smoothing parameter handling."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--cortical-smoothing', '15',
                                '--hippocampal-smoothing', '10']):
            args = parse_args()
            self.assertEqual(args.cortical_smoothing, 15)
            self.assertEqual(args.hippocampal_smoothing, 10)

    def test_thread_parameters(self):
        """Test thread count parameters."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--threads', '16',
                                '--wb-threads', '4']):
            args = parse_args()
            self.assertEqual(args.num_threads, 16)
            self.assertEqual(args.num_threads_wb, 4)

    def test_workbench_path(self):
        """Test custom Workbench path setting."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--wb-path', '/custom/path/to/wb_command']):
            args = parse_args()
            self.assertEqual(args.wb_path, '/custom/path/to/wb_command')

    def test_log_file(self):
        """Test log file parameter."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--log-file', '/path/to/logfile.log']):
            args = parse_args()
            self.assertEqual(args.log_file, '/path/to/logfile.log')

    def test_version_action(self):
        """Test that version action exits with version information."""
        with self.assertRaises(SystemExit):
            with patch('sys.argv', ['zbrains', '--version']):
                parse_args()

    def test_force_reprocess_flag(self):
        """Test force reprocess flag."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--force-reprocess']):
            args = parse_args()
            self.assertTrue(args.force_reprocess)

    def test_selective_processing(self):
        """Test flags for selective processing of different brain regions."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--no-cortex',
                                '--no-hippocampus']):
            args = parse_args()
            self.assertFalse(args.cortex)
            self.assertFalse(args.hippocampus)
            self.assertTrue(args.subcortical)  # Default is True

    def test_normative_columns_mismatch(self):
        """Test error handling when normative columns and dtypes don't match."""
        with self.assertRaises(SystemExit):
            with patch('sys.argv', ['zbrains',
                                    '--features', 'FA', 'ADC',
                                    '--control-demographics', 'controls.csv',
                                    '--patient-demographics', 'patients.csv',
                                    '--micapipe-dir', '/path/to/micapipe',
                                    '--hippunfold-dir', '/path/to/hippunfold',
                                    '--freesurfer-dir', '/path/to/freesurfer',
                                    '--normative-columns', 'AGE', 'SEX', 'EDUCATION',
                                    '--normative-dtypes', 'int', 'binary']):  # One fewer dtype than columns
                parse_args()


class TestZBrainsRunFunction(unittest.TestCase):
    """Test the run function in zbrains.run module."""

    @patch('zbrains.run.zbenv')
    @patch('zbrains.run.demographics')
    @patch('zbrains.run.zbdataset')
    @patch('os.makedirs')
    def test_run_initializes_environment(self, mock_makedirs, mock_zbdataset, mock_demographics, mock_zbenv):
        """Test that run correctly initializes environment with provided parameters."""
        # Setup mock return values
        mock_control_demo = MagicMock()
        mock_patient_demo = MagicMock()
        mock_control_dataset = MagicMock()
        mock_patient_dataset = MagicMock()
        mock_env = MagicMock()
        
        # Configure mocks
        mock_demographics.side_effect = [mock_control_demo, mock_patient_demo]
        mock_zbdataset.side_effect = [mock_control_dataset, mock_patient_dataset]
        mock_zbenv.return_value = mock_env
        
        # Configure valid_dataset attributes
        mock_control_dataset.valid_dataset = True
        mock_patient_dataset.valid_dataset = True
        
        # Call run with test parameters
        run(
            features=['FA', 'ADC'], 
            wb_path='/custom/wb_path',
            num_threads=16,
            num_threads_wb=4,
            control_demo_path='control.csv',
            patient_demo_path='patient.csv',
            micapipe_dir='/path/to/micapipe',
            hippunfold_dir='/path/to/hippunfold',
            freesurfer_dir='/path/to/freesurfer',
            output_directory='/output/dir',
            cortical_smoothing=15,
            hippocampal_smoothing=8,
            method='wscore'
        )
        
        # Verify environment was initialized correctly
        mock_zbenv.assert_called_once_with(
            connectome_workbench_path='/custom/wb_path',
            num_threads=16,
            num_threads_wb=4
        )

    @patch('zbrains.run.zbenv')
    @patch('zbrains.run.demographics')
    @patch('zbrains.run.zbdataset')
    @patch('os.makedirs')
    def test_run_creates_output_directory(self, mock_makedirs, mock_zbdataset, mock_demographics, mock_zbenv):
        """Test that run creates the output directory."""
        # Setup mock objects
        mock_demographics.return_value = MagicMock()
        mock_dataset = MagicMock()
        mock_dataset.valid_dataset = True
        mock_zbdataset.return_value = mock_dataset
        
        # Call run with test parameters
        run(
            features=['FA', 'ADC'], 
            control_demo_path='control.csv',
            patient_demo_path='patient.csv',
            micapipe_dir='/path/to/micapipe',
            hippunfold_dir='/path/to/hippunfold',
            freesurfer_dir='/path/to/freesurfer',
            output_directory='/output/dir'
        )
        
        # Verify output directory was created
        mock_makedirs.assert_called_with('/output/dir', exist_ok=True)

    @patch('zbrains.run.zbenv')
    @patch('zbrains.run.demographics')
    @patch('zbrains.run.zbdataset')
    def test_run_validates_features(self, mock_zbdataset, mock_demographics, mock_zbenv):
        """Test that run validates the features parameter."""
        # Test with empty features list
        with self.assertRaises(ValueError):
            run(
                features=[],  # Empty features list should raise error
                control_demo_path='control.csv',
                patient_demo_path='patient.csv',
                micapipe_dir='/path/to/micapipe',
                hippunfold_dir='/path/to/hippunfold',
                freesurfer_dir='/path/to/freesurfer'
            )

    @patch('zbrains.run.zbenv')
    @patch('zbrains.run.demographics')
    @patch('zbrains.run.zbdataset')
    def test_run_validates_required_directories(self, mock_zbdataset, mock_demographics, mock_zbenv):
        """Test that run validates required directories based on processing options."""
        # Test missing micapipe_dir when cortex=True
        with self.assertRaises(ValueError):
            run(
                features=['FA', 'ADC'],
                control_demo_path='control.csv',
                patient_demo_path='patient.csv',
                micapipe_dir=None,  # Missing required directory
                hippunfold_dir='/path/to/hippunfold',
                freesurfer_dir='/path/to/freesurfer',
                cortex=True
            )
        
        # Test missing hippunfold_dir when hippocampus=True
        with self.assertRaises(ValueError):
            run(
                features=['FA', 'ADC'],
                control_demo_path='control.csv',
                patient_demo_path='patient.csv',
                micapipe_dir='/path/to/micapipe',
                hippunfold_dir=None,  # Missing required directory
                freesurfer_dir='/path/to/freesurfer',
                hippocampus=True
            )

    @patch('zbrains.run.zbenv')
    @patch('zbrains.run.demographics')
    @patch('zbrains.run.zbdataset')
    @patch('os.makedirs')
    def test_run_processes_controls_and_patients(self, mock_makedirs, mock_zbdataset, mock_demographics, mock_zbenv):
        """Test that run processes both control and patient datasets."""
        # Setup mock objects
        mock_control_demo = MagicMock()
        mock_patient_demo = MagicMock()
        mock_control_dataset = MagicMock()
        mock_patient_dataset = MagicMock()
        mock_env = MagicMock()
        
        # Configure mocks
        mock_demographics.side_effect = [mock_control_demo, mock_patient_demo]
        mock_zbdataset.side_effect = [mock_control_dataset, mock_patient_dataset]
        mock_zbenv.return_value = mock_env
        
        # Set valid_dataset attribute for both datasets
        mock_control_dataset.valid_dataset = True
        mock_patient_dataset.valid_dataset = True
        
        # Call run with force_reprocess=True to ensure processing is called
        run(
            features=['FA', 'ADC'],
            control_demo_path='control.csv',
            patient_demo_path='patient.csv',
            micapipe_dir='/path/to/micapipe',
            hippunfold_dir='/path/to/hippunfold',
            freesurfer_dir='/path/to/freesurfer',
            force_reprocess=True
        )
        
        # Verify process was called for both datasets
        mock_control_dataset.process.assert_called_once()
        mock_patient_dataset.process.assert_called_once()

    @patch('zbrains.run.zbenv')
    @patch('zbrains.run.demographics')
    @patch('zbrains.run.zbdataset')
    @patch('os.makedirs')
    def test_run_analyzes_and_reports(self, mock_makedirs, mock_zbdataset, mock_demographics, mock_zbenv):
        """Test that run performs analysis and generates clinical reports."""
        # Setup mock objects
        mock_control_demo = MagicMock()
        mock_patient_demo = MagicMock()
        mock_control_dataset = MagicMock()
        mock_patient_dataset = MagicMock()
        mock_env = MagicMock()
        
        # Configure mocks
        mock_demographics.side_effect = [mock_control_demo, mock_patient_demo]
        mock_zbdataset.side_effect = [mock_control_dataset, mock_patient_dataset]
        mock_zbenv.return_value = mock_env
        
        # Set valid_dataset attribute for both datasets
        mock_control_dataset.valid_dataset = True
        mock_patient_dataset.valid_dataset = True
        
        # Call run
        run(
            features=['FA', 'ADC'],
            control_demo_path='control.csv',
            patient_demo_path='patient.csv',
            micapipe_dir='/path/to/micapipe',
            hippunfold_dir='/path/to/hippunfold',
            freesurfer_dir='/path/to/freesurfer',
            method='wscore'
        )
        
        # Verify analyze and clinical_report were called
        mock_patient_dataset.analyze.assert_called_once_with(
            output_directory='./zbrains_output',
            reference=mock_control_dataset,
            method='wscore'
        )
        mock_patient_dataset.clinical_report.assert_called_once()


class TestZBrainsDataset(unittest.TestCase):
    """Test the zbdataset class in zbrains.dataset module."""
    
    @patch('zbrains.dataset.demographics')
    def test_dataset_initialization(self, mock_demographics):
        """Test initialization of zbdataset with various parameters."""
        # Setup mock demographics
        mock_demo = MagicMock()
        mock_demographics.return_value = mock_demo
        
        # Test full initialization
        dataset = zbdataset(
            name="test_dataset",
            demographics=mock_demo,
            micapipe_directory="/path/to/micapipe",
            hippunfold_directory="/path/to/hippunfold",
            freesurfer_directory="/path/to/freesurfer",
            cortex=True,
            hippocampus=True,
            subcortical=True
        )
        
        # Verify attributes
        self.assertEqual(dataset.name, "test_dataset")
        self.assertEqual(dataset.demographics, mock_demo)
        self.assertEqual(dataset.micapipe_directory, "/path/to/micapipe")
        self.assertEqual(dataset.hippunfold_directory, "/path/to/hippunfold")
        self.assertEqual(dataset.freesurfer_directory, "/path/to/freesurfer")
        self.assertTrue(dataset.cortex)
        self.assertTrue(dataset.hippocampus)
        self.assertTrue(dataset.subcortical)
        self.assertFalse(dataset.valid_dataset)
        
        # Test initialization with disabled components
        dataset = zbdataset(
            name="limited_dataset",
            demographics=mock_demo,
            micapipe_directory="/path/to/micapipe",
            hippunfold_directory=None,
            freesurfer_directory=None,
            cortex=True,
            hippocampus=False,
            subcortical=False
        )
        
        # Verify attributes
        self.assertEqual(dataset.name, "limited_dataset")
        self.assertTrue(dataset.cortex)
        self.assertFalse(dataset.hippocampus)
        self.assertFalse(dataset.subcortical)
        self.assertIsNone(dataset.hippunfold_directory)
        self.assertIsNone(dataset.freesurfer_directory)


class TestCliToRunIntegration(unittest.TestCase):
    """Test the integration between CLI parsing and the run function."""

    @patch('zbrains.cli.run')
    @patch('zbrains.cli.parse_args')
    def test_args_passed_to_run(self, mock_parse_args, mock_run):
        """Test that CLI arguments are correctly passed to run function."""
        # Create mock args
        mock_args = MagicMock()
        mock_args.features = ['FA', 'ADC', 'thickness']
        mock_args.control_demo_path = 'data/participants_mics_hc_all.csv'
        mock_args.patient_demo_path = 'data/participants_mics_px_all.csv'
        mock_args.micapipe_dir = '/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0'
        mock_args.hippunfold_dir = '/data/mica3/BIDS_MICs/derivatives/hippunfold_v1.3.0'
        mock_args.freesurfer_dir = '/data/mica3/BIDS_MICs/derivatives/freesurfer'
        mock_args.output_directory = '/host/verges/tank/data/ian/zbrains_outputs'
        mock_args.cortical_smoothing = 10
        mock_args.hippocampal_smoothing = 5
        mock_args.num_threads = 16
        mock_args.num_threads_wb = 4
        mock_args.method = 'wscore'
        mock_args.cortex = True
        mock_args.hippocampus = True
        mock_args.subcortical = True
        mock_args.force_reprocess = False
        mock_args.verbose = True
        mock_args.wb_path = 'wb_command'
        mock_args.log_file = None
        mock_args.normative_columns = ['AGE', 'SEX']
        mock_args.normative_dtypes = ['int', 'binary']
        
        # Set return value for parse_args
        mock_parse_args.return_value = mock_args
        
        # Call main
        main()
        
        # Verify run was called with expected arguments
        mock_run.assert_called_once()
        
        # Check that all arguments were passed correctly
        called_args = mock_run.call_args[1]
        self.assertEqual(called_args['features'], mock_args.features)
        self.assertEqual(called_args['control_demo_path'], mock_args.control_demo_path)
        self.assertEqual(called_args['patient_demo_path'], mock_args.patient_demo_path)
        self.assertEqual(called_args['micapipe_dir'], mock_args.micapipe_dir)
        self.assertEqual(called_args['hippunfold_dir'], mock_args.hippunfold_dir)
        self.assertEqual(called_args['freesurfer_dir'], mock_args.freesurfer_dir)
        self.assertEqual(called_args['output_directory'], mock_args.output_directory)
        self.assertEqual(called_args['cortical_smoothing'], mock_args.cortical_smoothing)
        self.assertEqual(called_args['hippocampal_smoothing'], mock_args.hippocampal_smoothing)
        self.assertEqual(called_args['num_threads'], mock_args.num_threads)
        self.assertEqual(called_args['num_threads_wb'], mock_args.num_threads_wb)
        self.assertEqual(called_args['method'], mock_args.method)


if __name__ == '__main__':
    unittest.main()