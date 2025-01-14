#!/bin/bash

# Enable error handling
set -e

# Set the path to the dataset, or the folder containing the 'derivatives' folder
pth_dataset_ref="/data/mica3/BIDS_MICs"
pth_dataset="/data/mica3/BIDS_PNI"
# Set the directories for micapipe, hippunfold, and zbrains, which will be looked for in the 'derivates' folder
zbrains_dir="zbrains_clinical"
zbrains_dir_ref="zbrains_clinical"
micapipe_dir="micapipe_v0.2.0"
hippunfold_dir="hippunfold_v1.3.0"

# Set the paths to the demographic control and patient files
# The demo_controls are only needed for the analysis, and define the control samples to compare against.
# The demo_patients can be provided if desired, which will run all the patients in the list when the "all" keyword is used,
# otherwise the 'all' keyword will run every patient possible, given the micapipe and hippunfold availability, or, for the analysis
# it will run on all patients who have had zbrains proc run before.
demo_controls="/host/oncilla/local_raid/oualid/zbrains_csvs/participants_mics_hc.csv"

# Set the subject IDs and session IDs to 'all', using all patients defined in the PX_participants file.
subject_ids="sub-PNE002"
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
        --dataset_ref ${pth_dataset_ref} \
        --zbrains_ref ${zbrains_dir_ref} \
        --demo_ref ${demo_controls} \
        --column_map participant_id=ID session_id=SES \
        --smooth_ctx 10 \
        --smooth_hip 5 \
        --n_jobs 4 \
        --n_jobs_wb 4 \
        --label_ctx "white" \
	--wb_path /usr/bin/ \
        --verbose 2 \
        --volumetric 0 \
        --dicoms 0 \
        --pyinit=/data/mica1/03_projects/ian/anaconda3 


# Pause to keep the terminal open (optional, remove if not needed)
read -p "Press any key to continue..."
