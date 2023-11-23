#!/bin/bash


# Path for dataset in BIDS structure
root_path=/data/mica3/BIDS_MICs # /path/to/BIDS/

rawdir=${root_path}/rawdata
micapipedir=${root_path}/derivatives/micapipe_v0.2.0 
hippdir=${root_path}/derivatives/hippunfold_v1.3.0
outdir=${root_path}/derivatives/z-brains_yigutest

# csv file with ID and session for control participants to be processed
PATH_CSV_CONTROLS='/data_/mica1/03_projects/jessica/hackathon2023/lists/control.csv'

while IFS=',' read -r id ses rest
do
    ./z-brains -sub "$id" -ses "$ses" \
    -rawdir "${rawdir}" \
    -micapipedir "${micapipedir}" \
    -hippdir "${hippdir}" \
    -outdir "${outdir}" \
    -run proc \
    -mica -verbose 2

done <<< "$(tail -n +2 "${PATH_CSV_CONTROLS}")"
