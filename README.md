![alternate text](./src/data/zbrains_banner.png)
<div align="center">

# Multimodal lesion mapping in focal epilepsy with `zbrains`

![License](https://img.shields.io/badge/license-BSD-brightgreen) [![Version](https://img.shields.io/github/v/tag/MICA-MNI/z-brains)](https://github.com/MICA-MNI/z-brains) [![Documentation Status](https://readthedocs.org/projects/z-brains/badge/?version=latest&color=brightgreen)](https://z-brains.readthedocs.io/en/latest/?badge=latest) [![GitHub issues](https://img.shields.io/github/issues/MICA-MNI/z-brains?color=brightgreen)](https://github.com/MICA-MNI/z-brains/issues) [![GitHub stars](https://img.shields.io/github/stars/MICA-MNI/z-brains.svg?style=flat&label=%E2%9C%A8%EF%B8%8F%20be%20a%20stargazer&color=brightgreen)](https://github.com/MICA-MNI/z-brains/stargazers)

</div>

`zbrains` is developed by [MICA-Lab](https://mica-mni.github.io), affiliated with [McGill University](https://www.mcgill.ca/), the Montreal Neurological Institute and Hospital "[the Neuro](https://www.mcgill.ca/neuro/)," and the McConnell Brain Imaging Center ([BIC](https://www.mcgill.ca/bic/)).

This open access processing and analysis tool aims to identify patient-specific anomalies in brain morphology and microstructure, using features with previously demonstrated potential to accurately localize epileptogenic lesions. `zbrains` uses a set of known software dependencies developed by other groups and aggregated in a published pipeline [micapipe](https://github.com/MICA-MNI/micapipe).

Python Conversion Changelog

### **Bugs Fixed**

- qT1 was unused in reports due to a filepath bug, but is now fixed with no changes to paths
- Subcortical assymetry for combined features was broken, but is now fixed.

### **Functional Changes**
- "**--sub**" argument now accepts the keyword "**all**", which will attempt to run zbrains on all available subjects while processing, and while performing analysis, it will run on all subjects matching the pattern given in the argument (**--patient_prefix**). 
- "**--ses**" argument now also accepts the keyword "**all**", as well as now being optional. If the "**--ses**" argument is not given, it is assumed to be "**all**".
- "**--column_map**" argument has been given a default value of "**participant_id=ID session_id=SES**"
- "**--init**" argument has been removed, in favour of two arguments: "**--wb_path**" and "**--pyinit**"
- "**--wb_path**" specified the path to the workbench bin folder. It is set to "**/data/mica1/01_programs/workbench-1.4.2/bin_linux64**" by default
- "**--dicoms**" argument has been added, defaulting to 1. If set to 0, DICOMs will not be generated during analysis. This step can take a long time.
- "**--pyinit**" specifies the root folder of your desired Python environment. By default it is set to "**/data/mica1/01_programs/miniforge3**". It can be set to any conda-like environment manager that has an environment called "**zbrains**". Setting "**--pyinit=false**" will use the currently activated python environment.
- "**--patient_prefix**" specifies the prefix given to patient cases. Defaults to "**PX**"
- "**--delete_temps**" is a utility function which, if enabled, will not run any analyses, and will instead clean up any loose temp folders in the zbrains dataset derivates folder which may have been leftover from a system crash.

### **How to speed up your analyses with new zbrains**
- zbrain and zbrains_bash are now identical. Batch processing is handled by Python, instead of bash.
- the new "**--n_jobs**" argument will declare how many separate jobs you would like to run in parallel. Speed will increase with higher n_jobs, with larger gains seen when running zbrains on many subjects. Currently, speed increases in the analysis section are limited to **n_jobs**=5. Larger numbers of jobs will have no effect on speed.
- an additional argument **n_jobs_wb** dictates the number of jobs used for connectome workbench.

### **Backend Changes**
- There is now a testing framework for zbrains. Automated tests will determine if the whole pipeline runs following backend changes. This backend is written with pytest.
- All of the backend functions written in bash have been converted to python. 
- Zbrains is now importable, and importing src.zbrains will enable the programmatic use of zbrains.

### **Updated Example**

# Windows
Run.bat:
```batch
setlocal enabledelayedexpansion

REM Set the path to the dataset, or the folder containing the 'derivatives' folder
set "pth_dataset=E:\data"

REM Set the directories for micapipe, hippunfold, and zbrains, which will be looked for in the 'derivates' folder
set "micapipe_dir=micapipe"
set "hippunfold_dir=hippunfold"
set "zbrains_dir=zbrains"

REM Set the paths to the demographic control and patient files
REM The demo_controls are only needed for the analysis, and define the control samples compare against. 
REM The demo_patients can be provided if desired, which will run all the patients in the list when the "all" keyword is used, 
REM otherwise the 'all' keyword will run every patient possible, given the micapipe and hippunfold availability, or, for the analysis
REM it will run on all patients who have had zbrains proc run before.
set "demo_controls=E:\HC_participants.csv"
set "demo_patients=E:\PX_participants.csv"

REM Set the subject IDs and session IDs to 'all', using all patients defined in the PX_participants file.
set "sids=all" 
set "sess=all"

REM The code below runs zbrains preserving the old behaviour, with a smooth_ctx of 10, a smooth_hip of 5, and a label_ctx of 'white'
REM The new defaults for this are a smooth_ctx of 5, a smooth_hip of 2, and a label_ctx of 'midthickness'
REM Much of the new volumetric code is dependent on cortical midthickness, so it is recommended.
call conda activate zbrains
call python -m src.zbrains --run "proc analysis" ^
                   --sub "%sids%" ^
                   --micapipe %micapipe_dir% ^
                   --hippunfold %hippunfold_dir% ^
                   --dataset %pth_dataset% ^
                   --zbrains %zbrains_dir% ^
                   --demo_ref %demo_controls% ^
                   --demo %demo_patients% ^
                   --dataset_ref %pth_dataset% ^
                   --zbrains_ref %zbrains_dir% ^
                   --smooth_ctx 10 ^
                   --smooth_hip 5 ^
                   --n_jobs 4 ^
                   --n_jobs_wb 4 ^
		           --label_ctx "white" ^
                   --dicoms 0 ^
                   --wb_path "worbench_path/workbench/bin_windows64" ^
                   --column_map participant_id=ID session_id=SES ^
                   --verbose 2
    )
pause
```


The following example is tailored to use in the MICA Lab, targeting the BIDS_MICs dataset.
by calling the zbrains bash file directly, it will automatically set up environment variable for the workbench path, 
and activate our zbrains environment. We do need to provide a custom conda environment containing updated packages, using pyinit.

# Linux
Run.sh:
```bash
#!/bin/bash

# Enable error handling
set -e

# Set the path to the dataset, or the folder containing the 'derivatives' folder
pth_dataset="/data/mica3/BIDS_MICs"

# Set the directories for micapipe, hippunfold, and zbrains, which will be looked for in the 'derivates' folder
zbrains_dir="zbrains_volumetric"
micapipe_dir="micapipe_v0.2.0"
hippunfold_dir="hippunfold_v1.3.0"

# Set the paths to the demographic control and patient files
# The demo_controls are only needed for the analysis, and define the control samples to compare against.
# The demo_patients can be provided if desired, which will run all the patients in the list when the "all" keyword is used,
# otherwise the 'all' keyword will run every patient possible, given the micapipe and hippunfold availability, or, for the analysis
# it will run on all patients who have had zbrains proc run before.
demo_controls="/host/oncilla/local_raid/oualid/zbrains_csvs/participants_mics_hc.csv"

# Set the subject IDs and session IDs to 'all', using all patients defined in the PX_participants file.
subject_ids="all"
session_ids="all"

# The code below runs zbrains preserving the old behaviour, with a smooth_ctx of 10, a smooth_hip of 5, and a label_ctx of 'white'
# The new defaults for this are a smooth_ctx of 5, a smooth_hip of 2, and a label_ctx of 'midthickness'
# Much of the new volumetric code is dependent on cortical midthickness, so it is recommended.
./zbrains --run "proc analysis"\
        --sub "${subject_ids}" \
        --ses "${session_ids}" \
        --dataset ${pth_dataset} \
        --zbrains ${zbrains_dir} \
        --micapipe ${micapipe_dir} \
        --hippunfold ${hippunfold_dir} \
        --dataset_ref ${pth_dataset} \
        --zbrains_ref ${zbrains_dir} \
        --demo_ref ${demo_controls} \
        --column_map participant_id=ID session_id=SES \
        --smooth_ctx 10 \
        --smooth_hip 5 \
        --n_jobs 4 \
        --n_jobs_wb 4 \
        --label_ctx "white" \
        --dicoms 0 \
        --verbose 2 \
        --pyinit=/data/mica1/03_projects/ian/anaconda3 


# Pause to keep the terminal open (optional, remove if not needed)
read -p "Press any key to continue..."


