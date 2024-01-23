#!/bin/bash

# ANTs, Workbench and python
ZBRAINS=$(dirname "$(realpath "$0")")
source "${ZBRAINS}/functions/init.sh"


pth_dataset="/data_/mica3/BIDS_MICs"
zbrains_dir="zbrains_reborn"
micapipe_dir="micapipe_v0.2.0"
hippunfold_dir="hippunfold_v1.3.0"

demo_controls="/host/oncilla/local_raid/oualid/zbrains_csvs/participants_mics_hc.csv"

subject_ids=(sub-PX122)
session_ids=(ses-01)

#subject_ids=(sub-PX088)
#session_ids=(ses-01)

for i in "${!subject_ids[@]}"
do
  sid=${subject_ids[$i]}
  ses=${session_ids[$i]:-""}
  echo -e "$i\t\tID: $sid, SES: $ses"


  ./zbrains --run proc \
            --sub "${sid}" \
            --ses "${ses}" \
            --dataset ${pth_dataset} \
            --zbrains ${zbrains_dir} \
            --micapipe ${micapipe_dir} \
            --hippunfold ${hippunfold_dir} \
            --struct cortex subcortex \
            --verbose 2
done


for i in "${!subject_ids[@]}"
do
  sid=${subject_ids[$i]}
  ses=${session_ids[$i]:-""}
  echo -e "$i\t\tID: $sid, SES: $ses"


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
done

