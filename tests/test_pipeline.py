#!/usr/bin/env python3
import unittest
import sys
import tempfile
import os
from unittest.mock import patch, MagicMock
from zbrains.cli import parse_args, main

"""
Tests for the Z-Brains CLI module.
"""



class TestZBrainsCliParser(unittest.TestCase):
    """Test the argument parser for the Z-Brains CLI."""

    def test_required_arguments(self):
        """Test that required arguments are enforced."""
        # Test missing --features
        with self.assertRaises(SystemExit):
            with patch('sys.argv', ['zbrains',
                                    '--control-demographics', 'controls.csv',
                                    '--patient-demographics', 'patients.csv',
                                    '--micapipe-dir', '/path/to/micapipe',
                                    '--hippunfold-dir', '/path/to/hippunfold',
                                    '--freesurfer-dir', '/path/to/freesurfer']):
                parse_args()

        # Test missing --control-demographics
        with self.assertRaises(SystemExit):
            with patch('sys.argv', ['zbrains',
                                    '--features', 'FA', 'ADC',
                                    '--patient-demographics', 'patients.csv',
                                    '--micapipe-dir', '/path/to/micapipe',
                                    '--hippunfold-dir', '/path/to/hippunfold',
                                    '--freesurfer-dir', '/path/to/freesurfer']):
                parse_args()

        # Test missing --patient-demographics
        with self.assertRaises(SystemExit):
            with patch('sys.argv', ['zbrains',
                                    '--features', 'FA', 'ADC',
                                    '--control-demographics', 'controls.csv',
                                    '--micapipe-dir', '/path/to/micapipe',
                                    '--hippunfold-dir', '/path/to/hippunfold',
                                    '--freesurfer-dir', '/path/to/freesurfer']):
                parse_args()

    def test_minimal_arguments(self):
        """Test parsing with minimal required arguments."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer']):
            args = parse_args()
            self.assertEqual(args.features, ['FA', 'ADC'])
            self.assertEqual(args.control_demo_path, 'controls.csv')
            self.assertEqual(args.patient_demo_path, 'patients.csv')
            self.assertEqual(args.micapipe_dir, '/path/to/micapipe')
            self.assertEqual(args.hippunfold_dir, '/path/to/hippunfold')
            self.assertEqual(args.freesurfer_dir, '/path/to/freesurfer')
            # Check defaults are set
            self.assertEqual(args.output_directory, './zbrains_output')
            self.assertEqual(args.cortical_smoothing, 10)
            self.assertEqual(args.hippocampal_smoothing, 5)
            self.assertEqual(args.method, 'wscore')
            self.assertEqual(args.num_threads, 4)
            self.assertEqual(args.num_threads_wb, 2)
            self.assertTrue(args.cortex)
            self.assertTrue(args.hippocampus)
            self.assertTrue(args.subcortical)
            self.assertFalse(args.force_reprocess)
            self.assertTrue(args.verbose)
            self.assertEqual(args.normative_columns, ['AGE', 'SEX'])
            self.assertEqual(args.normative_dtypes, ['int', 'binary'])

    def test_full_argument_set(self):
        """Test parsing with the example full command argument set."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC', 'thickness', 'qT1', 
                                'qT1*Blur', 'FLAIR', 'FLAIR*Blur',
                                '--control-demographics', 'data/participants_mics_hc_all.csv',
                                '--patient-demographics', 'data/participants_mics_px_all.csv',
                                '--micapipe-dir', '/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0',
                                '--hippunfold-dir', '/data/mica3/BIDS_MICs/derivatives/hippunfold_v1.3.0',
                                '--freesurfer-dir', '/data/mica3/BIDS_MICs/derivatives/freesurfer',
                                '--output-dir', '/host/verges/tank/data/ian/zbrains_outputs',
                                '--cortical-smoothing', '10',
                                '--hippocampal-smoothing', '5',
                                '--threads', '16',
                                '--wb-threads', '4',
                                '--method', 'wscore']):
            args = parse_args()
            self.assertEqual(args.features, ['FA', 'ADC', 'thickness', 'qT1', 
                                            'qT1*Blur', 'FLAIR', 'FLAIR*Blur'])
            self.assertEqual(args.control_demo_path, 'data/participants_mics_hc_all.csv')
            self.assertEqual(args.patient_demo_path, 'data/participants_mics_px_all.csv')
            self.assertEqual(args.micapipe_dir, '/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0')
            self.assertEqual(args.hippunfold_dir, '/data/mica3/BIDS_MICs/derivatives/hippunfold_v1.3.0')
            self.assertEqual(args.freesurfer_dir, '/data/mica3/BIDS_MICs/derivatives/freesurfer')
            self.assertEqual(args.output_directory, '/host/verges/tank/data/ian/zbrains_outputs')
            self.assertEqual(args.cortical_smoothing, 10)
            self.assertEqual(args.hippocampal_smoothing, 5)
            self.assertEqual(args.num_threads, 16)
            self.assertEqual(args.num_threads_wb, 4)
            self.assertEqual(args.method, 'wscore')

    def test_flag_options(self):
        """Test flag options like --no-cortex, --quiet, etc."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--no-cortex',
                                '--no-hippocampus',
                                '--no-subcortical',
                                '--force-reprocess',
                                '--quiet']):
            args = parse_args()
            self.assertFalse(args.cortex)
            self.assertFalse(args.hippocampus)
            self.assertFalse(args.subcortical)
            self.assertTrue(args.force_reprocess)
            self.assertFalse(args.verbose)

    def test_custom_normative_parameters(self):
        """Test custom normative columns and data types."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--normative-columns', 'AGE', 'SEX', 'EDUCATION',
                                '--normative-dtypes', 'int', 'binary', 'int']):
            args = parse_args()
            self.assertEqual(args.normative_columns, ['AGE', 'SEX', 'EDUCATION'])
            self.assertEqual(args.normative_dtypes, ['int', 'binary', 'int'])

    def test_alternative_statistical_method(self):
        """Test setting alternative statistical method."""
        with patch('sys.argv', ['zbrains',
                                '--features', 'FA', 'ADC',
                                '--control-demographics', 'controls.csv',
                                '--patient-demographics', 'patients.csv',
                                '--micapipe-dir', '/path/to/micapipe',
                                '--hippunfold-dir', '/path/to/hippunfold',
                                '--freesurfer-dir', '/path/to/freesurfer',
                                '--method', 'zscore']):
            args = parse_args()
            self.assertEqual(args.method, 'zscore')

        # Test invalid method
        with self.assertRaises(SystemExit):
            with patch('sys.argv', ['zbrains',
                                    '--features', 'FA', 'ADC',
                                    '--control-demographics', 'controls.csv',
                                    '--patient-demographics', 'patients.csv',
                                    '--micapipe-dir', '/path/to/micapipe',
                                    '--hippunfold-dir', '/path/to/hippunfold',
                                    '--freesurfer-dir', '/path/to/freesurfer',
                                    '--method', 'invalid_method']):
                parse_args()


class TestZBrainsCliMain(unittest.TestCase):
    """Test the main function for the Z-Brains CLI."""

    @patch('zbrains.cli.parse_args')
    @patch('zbrains.cli.run')
    def test_main_successful_execution(self, mock_run, mock_parse_args):
        """Test successful execution of the main function."""
        # Setup mock return values
        mock_args = MagicMock()
        mock_args.features = ['FA', 'ADC']
        mock_args.control_demo_path = 'controls.csv'
        # ... add more arguments as needed
        mock_parse_args.return_value = mock_args
        
        # Call main and check result
        result = main()
        
        # Check if run was called with expected arguments
        mock_run.assert_called_once()
        args_dict = vars(mock_args)
        mock_run.assert_called_with(**args_dict)
        
        # Check that main returns 0 on success
        self.assertEqual(result, 0)

    @patch('zbrains.cli.parse_args')
    @patch('zbrains.cli.run')
    @patch('sys.stderr')
    def test_main_keyboard_interrupt(self, mock_stderr, mock_run, mock_parse_args):
        """Test handling of KeyboardInterrupt in main."""
        # Make run raise KeyboardInterrupt
        mock_run.side_effect = KeyboardInterrupt()
        
        # Call main and check result
        result = main()
        
        # Check error message
        mock_stderr.write.assert_called()
        # Check that main returns 130 on KeyboardInterrupt
        self.assertEqual(result, 130)

    @patch('zbrains.cli.parse_args')
    @patch('zbrains.cli.run')
    @patch('sys.stderr')
    def test_main_general_exception(self, mock_stderr, mock_run, mock_parse_args):
        """Test handling of general exceptions in main."""
        # Make run raise Exception
        mock_run.side_effect = Exception("Test error")
        
        # Call main and check result
        result = main()
        
        # Check error message
        mock_stderr.write.assert_called()
        # Check that main returns 1 on Exception
        self.assertEqual(result, 1)

    @patch('zbrains.cli.parse_args')
    @patch('zbrains.cli.run')
    def test_main_with_log_file(self, mock_run, mock_parse_args):
        """Test execution with log file specified."""
        # Create a temporary directory for log file
        with tempfile.TemporaryDirectory() as temp_dir:
            log_file = os.path.join(temp_dir, "test.log")
            
            # Setup mock return values
            mock_args = MagicMock()
            mock_args.features = ['FA', 'ADC']
            mock_args.log_file = log_file
            mock_parse_args.return_value = mock_args
            
            # Call main
            result = main()
            
            # Check that run was called with log_file parameter
            mock_run.assert_called_once()
            args_dict = vars(mock_args)
            self.assertEqual(args_dict['log_file'], log_file)
            
            # Check that main returns 0 on success
            self.assertEqual(result, 0)


if __name__ == '__main__':
    unittest.main()