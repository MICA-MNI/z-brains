PATH+=:/home/igoodall/anaconda3/envs/zbrains/bin
# path for dataset in BIDS structure
pth_dataset="/mnt/e/data"

micapipe_dir="micapipe"
hippunfold_dir="hippunfold"
zbrains_dir="zbrains_out"

# csv file with ID and session for control participants for comparison
demo_controls="/mnt/e/data/part.csv"

# specify the list of subject IDs along with corresponding session


sids=()
for d in "$pth_dataset"/derivatives/"$micapipe_dir"/*/ ; do
   part2=$(basename "$d")
   # Extract the subject ID without the "sub-" prefix
   sid_part2="${part2#sub-}"
   # Add the "sub-" prefix to the sids array
   sids+=("sub-$sid_part2")
done

echo "${sid}"
echo "${sids}"

sids=(sub-PX070)
for sid in "${sids[@]}" ; do
   echo $sid
   sess=()
   for d in "$pth_dataset"/derivatives/"$micapipe_dir"/"$sid"/*/ ; do
      part2=$(basename "$d")
      # Extract the subject ID without the "sub-" prefix
      ses_part2="${part2#ses-}"
      # Add the "sub-" prefix to the sids array
      sess+=("ses-$ses_part2")
   done
   for ses in "${sess[@]}" ; do
      ./z-brains/zbrains --run proc analysis \
                  --sub "${sid}" \
                  --ses "${ses}" \
                  --micapipe ${micapipe_dir} \
                  --hippunfold ${hippunfold_dir} \
                  --dataset ${pth_dataset} \
                  --zbrains ${zbrains_dir} \
                  --demo_ref ${demo_controls} \
                  --dataset_ref ${pth_dataset} \
                  --zbrains_ref ${zbrains_dir} \
                  --column_map participant_id=ID session_id=SES \
                  --verbose 2
   done
done