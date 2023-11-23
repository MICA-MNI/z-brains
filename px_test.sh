#!/bin/bash


# Path for dataset in BIDS structure
root_path=/data/mica3/BIDS_MICs # /path/to/BIDS/

rawdir=${root_path}/rawdata
micapipedir=${root_path}/derivatives/micapipe_v0.2.0 
hippdir=${root_path}/derivatives/hippunfold_v1.3.0
outdir=${root_path}/derivatives/z-brains_yigutest


# Specify the list of subject IDs along with corresponding session
px_id=(PX001 PX002 PX003)
px_ses=(1 1 1)

# csv file with ID and session for control participants to be processed
PATH_CSV_CONTROLS='/data_/mica1/03_projects/jessica/hackathon2023/lists/control.csv'

i=0
for id in "${px_id[@]}"
do
    ses=${px_ses[$i]}
	
    ./z-brains -sub "$id" -ses "$ses" \
    -rawdir "${rawdir}" \
    -micapipedir "${micapipedir}" \
    -hippdir "${hippdir}" \
    -outdir "${outdir}" \
    -approach "zscore" \
    -demo_cn "${PATH_CSV_CONTROLS}" \
    -mica -verbose 2

    i=$((i+1))
    
done
