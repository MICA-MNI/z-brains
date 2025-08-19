# Z-Brains v2.0.0

**Neuroimaging Processing and Analysis Pipeline**

Z-Brains is a comprehensive neuroimaging pipeline for processing and analyzing structural and diffusion MRI data. It performs patient-control comparisons using normative modeling approaches to identify brain abnormalities.

## Installation

```bash
# Clone repository and install
git clone <repository-url>
cd z-brains-v2.0
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install .
```

Verify installation:
```bash
zbrains --version
python -c "import zbrains; print('Z-Brains installed successfully')"
```

## Quick Start

### View Help
```bash
# Show comprehensive help with examples
zbrains

# Show standard help
zbrains --help
```

### Basic Usage
```bash
zbrains --features FA thickness \
        --control-demographics controls.csv \
        --patient-demographics patients.csv \
        --micapipe-dir /data/micapipe \
        --hippunfold-dir /data/hippunfold \
        --freesurfer-dir /data/freesurfer \
        --wb-path wb_command
```

## Pipeline Workflow

1. **Load Demographics** - Control and patient demographic data
2. **Process Controls** - Create normative models from healthy controls
3. **Process Patients** - Apply same preprocessing to patient data
4. **Statistical Analysis** - Compare patients to controls (zscore/wscore methods)
5. **Clinical Reports** - Generate individual brain maps and reports

## Required Arguments

### Demographics Files
Both control and patient CSV files must include:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `participant_id` | string | Subject identifier | sub-001 |
| `session_id` | string | Session identifier | ses-01 |
| `AGE` | integer | Age in years | 25 |
| `SEX` | binary | Sex (M/F, 0/1, etc.) | M |

**Example CSV format:**
```csv
participant_id,session_id,AGE,SEX
sub-001,ses-01,25,M
sub-002,ses-01,30,F
sub-003,ses-01,28,M
```

Optional columns for custom normative modeling:
- `EDUCATION` - Years of education
- `HANDEDNESS` - Handedness score
- Any other demographic variables

### Derivative Directories (BIDS format)

```
micapipe_dir/
├── sub-001/
│   └── ses-01/
│       ├── anat/
│       ├── surf/
│       └── maps/

hippunfold_dir/
├── hippunfold/
│   └── sub-001/
│       └── ses-01/
│           └── surf/

freesurfer_dir/
├── sub-001/
│   └── ses-01/
│       └── (FreeSurfer outputs)
```

## Available Features

| Feature | Description | Brain Structures |
|---------|-------------|------------------|
| `FA` | Fractional Anisotropy (DTI) | Cortical, Hippocampal, Subcortical |
| `ADC` | Apparent Diffusion Coefficient | Cortical, Hippocampal, Subcortical |
| `thickness` | Cortical thickness | Cortical |
| `qT1` | Quantitative T1 relaxometry | Cortical, Hippocampal, Subcortical |
| `FLAIR` | FLAIR intensity | Cortical, Hippocampal, Subcortical |
| `qT1*blur` | qT1 with spatial blurring | Cortical only |
| `FLAIR*blur` | FLAIR with spatial blurring | Cortical only |

## Command Line Reference

### Required Arguments
```bash
--features FEATURE [FEATURE ...]        # Features to process
--control-demographics PATH             # Control demographics CSV
--patient-demographics PATH             # Patient demographics CSV  
--micapipe-dir PATH                     # Micapipe derivatives directory
--hippunfold-dir PATH                   # Hippunfold derivatives directory
--freesurfer-dir PATH                   # FreeSurfer derivatives directory
--wb-path PATH                          # Path to wb_command
```

### Output Options
```bash
--output-dir PATH                       # Output directory (default: ./zbrains_output)
--log-file PATH                         # Log file path (default: console only)
```

### Processing Parameters
```bash
--cortical-smoothing INT                # Cortical smoothing kernel mm (default: 10)
--hippocampal-smoothing INT             # Hippocampal smoothing kernel mm (default: 5)
--method {wscore,zscore}                # Statistical method (default: wscore)
--threads INT                           # Processing threads (default: 4)
--wb-threads INT                        # Workbench threads (default: 2)
```

### Structure Selection
```bash
--no-cortex                            # Skip cortical processing
--no-hippocampus                       # Skip hippocampal processing
--no-subcortical                       # Skip subcortical processing
```

### Normative Modeling
```bash
--normative-columns COL [COL ...]      # Demographic columns (default: AGE SEX)
--normative-dtypes TYPE [TYPE ...]     # Data types (default: int binary)
                                       # Options: int, float, binary, categorical
```

### Other Options
```bash
--force-reprocess                      # Force reprocessing existing data
--quiet                               # Suppress detailed output
--version                             # Show version
```

## Statistical Methods

### wscore (Weighted Score) - **Recommended**
- Uses normative modeling with demographic covariates
- Accounts for age, sex, and other variables in control distribution
- Provides weighted z-scores based on expected values
- **Recommended for clinical applications**

### zscore (Z-Score)
- Simple standardization using control mean and standard deviation
- Does not account for demographic variables
- Faster computation but less precise normalization

## Example Commands

### Basic Analysis
```bash
zbrains --features FA thickness \
        --control-demographics data/controls.csv \
        --patient-demographics data/patients.csv \
        --micapipe-dir /data/micapipe \
        --hippunfold-dir /data/hippunfold \
        --freesurfer-dir /data/freesurfer \
        --wb-path wb_command
```

### Advanced Multi-Feature Analysis
```bash
zbrains --features FA ADC thickness qT1 FLAIR \
        --control-demographics data/controls.csv \
        --patient-demographics data/patients.csv \
        --micapipe-dir /data/micapipe \
        --hippunfold-dir /data/hippunfold \
        --freesurfer-dir /data/freesurfer \
        --wb-path /usr/local/bin/wb_command \
        --output-dir /results/zbrains \
        --cortical-smoothing 10 \
        --hippocampal-smoothing 5 \
        --method wscore \
        --threads 16 \
        --normative-columns AGE SEX EDUCATION \
        --normative-dtypes int binary int
```

### Cortical-Only Analysis
```bash
zbrains --features thickness FA qT1*blur \
        --control-demographics data/controls.csv \
        --patient-demographics data/patients.csv \
        --micapipe-dir /data/micapipe \
        --freesurfer-dir /data/freesurfer \
        --wb-path wb_command \
        --no-hippocampus \
        --no-subcortical
```

### High-Performance Computing
```bash
zbrains --features FA ADC thickness \
        --control-demographics data/controls.csv \
        --patient-demographics data/patients.csv \
        --micapipe-dir /data/micapipe \
        --hippunfold-dir /data/hippunfold \
        --freesurfer-dir /data/freesurfer \
        --wb-path wb_command \
        --threads 32 \
        --wb-threads 8 \
        --log-file zbrains.log \
        --force-reprocess
```

## Brain Structure Support

### Cortical
- Surface-based analysis on fsLR-32k and fsLR-5k meshes
- Midthickness and white matter surfaces
- Features: thickness, FA, ADC, qT1, FLAIR, blur variants

### Hippocampal
- Hippocampus-specific surface analysis
- Unfolded hippocampal coordinates
- Features: FA, ADC, qT1, FLAIR (no blur variants)

### Subcortical
- Volume-based analysis of deep brain structures
- FreeSurfer segmentation-based
- Features: volumes, FA, ADC, qT1, FLAIR (no blur variants)

## Output Structure

```
zbrains_output/
├── sub-001/
│   └── ses-01/
│       ├── maps/
│       │   ├── cortex/              # Processed cortical maps
│       │   ├── hippocampus/         # Processed hippocampal maps
│       │   └── subcortical/         # Subcortical statistics
│       ├── structural/              # Surface files
│       ├── analysis/                # Statistical comparison results
│       └── *.pdf                    # Clinical reports
├── sub-002/
│   └── ses-01/
│       └── ...
└── logs/                           # Processing logs
```

## Dependencies

### Required Software
- **Connectome Workbench** (`wb_command`) - For surface processing
- **Python 3.8+** - Runtime environment

### Installation Check
```bash
# Verify wb_command is available
wb_command -version

# Check Python environment
python --version
pip list | grep zbrains
```

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| `ModuleNotFoundError: No module named 'zbrains'` | Package not installed | `pip install --force-reinstall .` |
| `wb_command not found` | Workbench not in PATH | Install Workbench or use `--wb-path /full/path/to/wb_command` |
| `FileNotFoundError` for derivatives | Invalid paths | Check directory paths exist and are readable |
| Demographics normalization errors | Missing/mismatched columns | Ensure both CSVs have AGE, SEX with consistent formats |
| Processing very slow | Resource constraints | Reduce `--threads` or smoothing parameters |
| Memory errors | Large datasets | Process subsets or increase available RAM |

### Clean Reinstall
```bash
pip uninstall -y zbrains
pip install .
```

### Debugging
```bash
# Run with detailed logging
zbrains --features FA thickness \
        --control-demographics controls.csv \
        --patient-demographics patients.csv \
        --micapipe-dir /data/micapipe \
        --hippunfold-dir /data/hippunfold \
        --freesurfer-dir /data/freesurfer \
        --wb-path wb_command \
        --log-file debug.log \
        --threads 1

# Monitor progress
tail -f debug.log
```
### Important Notes
- **Control dataset must be processed successfully before patient analysis**
- Missing data is handled gracefully - subjects with incomplete data are skipped
- Processing is parallelized across subjects for efficiency
- All intermediate files are preserved for quality control
- Clinical reports include brain visualizations and statistical summaries

## Batch Processing

### SLURM Example
```bash
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

source .venv/bin/activate

zbrains --features FA ADC thickness \
        --control-demographics controls.csv \
        --patient-demographics patients.csv \
        --micapipe-dir /data/micapipe \
        --hippunfold-dir /data/hippunfold \
        --freesurfer-dir /data/freesurfer \
        --wb-path wb_command \
        --threads 16 \
        --log-file zbrains_${SLURM_JOB_ID}.log
```


**New CLI equivalent:**
```bash
zbrains --features FA ADC thickness \
        --control-demographics controls.csv \
        --patient-demographics patients.csv \
        --micapipe-dir /data/micapipe \
        --hippunfold-dir /data/hippunfold \
        --freesurfer-dir /data/freesurfer \
        --wb-path wb_command \
        --threads 4
```

## Support

For issues, please provide:
- Complete command executed
- Full error traceback
- Git commit hash (`git rev-parse HEAD`)
- Sample demographics CSV headers
- Directory listing of derivatives

**Repository:** https://github.com/MICA-MNI/z-brains  
**License:** BSD-3-Clause

**For MICALab Users**

Here is a working command: 
```bash
zbrains   --features FA ADC thickness qT1 "qT1*Blur" FLAIR "FLAIR*Blur"   --control-demographics /host/verges/tank/data/ian/z-brains-v2.0/data/participants_mics_hc_all.csv   --patient-demographics /host/verges/tank/data/ian/z-brains-v2.0/data/participants_mics_px_all.csv   --micapipe-dir /data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0   --hippunfold-dir /data/mica3/BIDS_MICs/derivatives/hippunfold_v1.3.0   --freesurfer-dir /data/mica3/BIDS_MICs/derivatives/freesurfer   --output-dir /host/verges/tank/data/ian/zbrains_outputs   --cortical-smoothing 10   --hippocampal-smoothing 5   --threads 16   --wb-threads 4   --method wscore --wb-path "/usr/bin/"

```