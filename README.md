![alternate text](./src/data/zbrains_banner.png)
<div align="center">

# Multimodal lesion mapping in focal epilepsy with `zbrains`

![License](https://img.shields.io/badge/license-BSD-brightgreen) [![Version](https://img.shields.io/github/v/tag/MICA-MNI/z-brains)](https://github.com/MICA-MNI/z-brains) [![Documentation Status](https://readthedocs.org/projects/z-brains/badge/?version=latest&color=brightgreen)](https://z-brains.readthedocs.io/en/latest/?badge=latest) [![GitHub issues](https://img.shields.io/github/issues/MICA-MNI/z-brains?color=brightgreen)](https://github.com/MICA-MNI/z-brains/issues) [![GitHub stars](https://img.shields.io/github/stars/MICA-MNI/z-brains.svg?style=flat&label=%E2%9C%A8%EF%B8%8F%20be%20a%20stargazer&color=brightgreen)](https://github.com/MICA-MNI/z-brains/stargazers)

</div>

`zbrains` is developed by [MICA-Lab](https://mica-mni.github.io), affiliated with [McGill University](https://www.mcgill.ca/), the Montreal Neurological Institute and Hospital "[the Neuro](https://www.mcgill.ca/neuro/)," and the McConnell Brain Imaging Center ([BIC](https://www.mcgill.ca/bic/)).

This open access processing and analysis tool aims to identify patient-specific anomalies in brain morphology and microstructure, using features with previously demonstrated potential to accurately localize epileptogenic lesions. `zbrains` uses a set of known software dependencies developed by other groups and aggregated in a published pipeline [micapipe](https://github.com/MICA-MNI/micapipe).

Python Conversion Changelog
- 
### **Functional Changes**
- "**--sub**" argument now accepts the keyword "**all**", which will attempt to run zbrains on all available subjects while processing, and while performing analysis, it will run on all subjects matching the pattern given in the argument (**--patient_prefix**). 
- "**--ses**" argument now also accepts the keyword "**all**", as well as now being optional. If the "**--ses**" argument is not given, it is assumed to be "**all**".
- "**--column_map**" argument has been given a default value of "**participant_id=ID session_id=SES**"
- "**--init**" argument has been removed, in favour of two arguments: "**--wb_path**" and "**--pyinit**"
- "**--wb_path**" specified the path to the workbench bin folder. It is set to "**/data/mica1/01_programs/workbench-1.4.2/bin_linux64**" by default
- "**--pyinit**" specifies the root folder of your desired Python environment. By default it is set to "**/data/mica1/01_programs/miniforge3**". It can be set to any conda-like environment manager that has an environment called "**zbrains**". Setting "**--pyinit=false**" will use the currently activated python environment.
- "**--patient_prefix**" specifies the prefix given to patient cases. Defaults to "**PX**"
- "**--delete_temps**" is a utility function which, if enabled, will not run any analyses, and will instead clean up any loose temp folders in the zbrains dataset derivates folder which may have been leftover from a system crash.

### **How to speed up your analyses with new zbrains**
- zbrain and zbrains_bash are now identical. Batch processing is handled by Python, instead of bash.
- the new "**--n_jobs**" argument will declare how many separate jobs you would like to run in parallel. Speed will increase with higher n_jobs, with larger gains seen when running zbrains on many subjects. Currently, speed increases in the analysis section are limited to **n_jobs**=5. Larger numbers of jobs will have no effect on speed.

### **Backend Changes**
- There is now a testing framework for zbrains. Automated tests will determine if the whole pipeline runs following backend changes. This backend is written with pytest.
- All of the backend functions written in bash have been converted to python. 
- Zbrains is now importable, and importing src.zbrains will enable the programmatic use of zbrains.

## Installation

Make sure to set MICAPIPE and ZBRAINS variables, and add their function to your PATH. For example:

```bash
export MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
export PATH=${PATH}:${MICAPIPE}:${MICAPIPE}/functions
source ${MICAPIPE}/functions/init.sh

export ZBRAINS=/data/mica1/03_projects/jordand/z-brains
export PATH=${PATH}:${ZBRAINS}:${ZBRAINS}/functions
```

## Tutorial

You must be inside the `zbrains` repository to run the following commands.

```bash
cd /path/to/zbrains/repository
ZBRAINS=$(dirname "$(realpath "$0")")
source "${ZBRAINS}/functions/init.sh"
```

`zbrains` requires input and output directories:

- `pth_dataset` points to the BIDS-format dataset path where micapipe, hippunfold, and zbrains directories will be stored
- `micapipe_dir` contains the output of `micapipe` previously run on the BIDS dataset (also BIDS-format)
- `hippunfold_dir` contains the output of `hippunfold` previously run on the BIDS dataset (also BIDS-format)
- `zbrains_dir` points to the directory that will hold `z-brains` outputs

```bash
# path for dataset in BIDS structure
pth_dataset="/path/to/BIDS_dataset"

micapipe_dir="name_of_micapipe_folder"
hippunfold_dir="name_of_hippunfold_folder"
zbrains_dir="name_of_z-brains_folder"
```

### Preparing control data

A `.csv` file with ID and session for healthy controls is required.

```bash
# csv file with ID and session for control participants to be processed
demo_controls='/path/to/control/participants.csv'
```

To process features for healthy controls as a batch, run the following. Note that column names of the healthy controls `.csv` file must be specified under `--column_map`.

```bash
./zbrains_batch --run proc \
   --demo "${demo_controls}" \
   --dataset "${pth_dataset}" \
   --zbrains ${zbrains_dir} \
   --micapipe ${micapipe_dir} \
   --hippunfold ${hippunfold_dir} \
   --column_map participant_id=ID session_id=SES \
   --verbose 2 \
   --scheduler_options "-n 20" #specify threads here
```

### Processing and analyzing patient features

```bash
# specify the list of subject IDs along with corresponding session
subject_ids=(sub-PX001 sub-PX002)
session_ids=(ses-01 ses-01)

# csv file with ID and session for control participants for comparison
demo_controls='/path/to/control/participants.csv'

./zbrains --run analysis \
          --sub "${sid}" \
          --ses "${ses}" \
          --dataset ${pth_dataset} \
          --zbrains ${zbrains_dir} \
          --demo_ref ${demo_controls} \
          --dataset_ref ${pth_dataset} \
          --zbrains_ref ${zbrains_dir} \
          --column_map participant_id=ID session_id=SES \
          --verbose 2
```
