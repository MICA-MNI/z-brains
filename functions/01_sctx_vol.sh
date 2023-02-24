#!/bin/bash
#
# Subcortical volume:
#
# Generates subcortical volume outputs:
#
# This workflow makes use of freesurfer outputs and custom python scripts
#
# Atlas and templates are available from:
#
# https://github.com/MICA-MNI/micapipe/tree/master/parcellations
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out Directory
#
umask 003
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
PROC=$8
export OMP_NUM_THREADS=$threads
here=$(pwd)

echo "PROC:"
echo $PROC
#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

#------------------------------------------------------------------------------#
Title "subcortical volume\n\t\tmicapipe-z $Version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Saving temporal dir: $nocleanup"

# Timer
aloita=$(date +%s)
Nsteps=0

# Check inputs: Freesurfer folder
if [ ! -d "$dir_freesurfer" ]; then Error "No freesurfer folder found for idBIDS $id"; exit; fi

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe-z_sctx_vol_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${proc_struct}/surfaces/morphology/sctx_volume/"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir="${dir_freesurfer}/stats/"

#------------------------------------------------------------------------------#
### subcortical volume ###
mkdir -p "${proc_struct}/surfaces/morphology/sctx_volume/"

# Create file header
echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal,ICV" > \
     "${outDir}/${idBIDS}_sctx_volume.csv"
printf "%s,"  "${idBIDS}" >> "${outDir}/${idBIDS}_sctx_volume.csv"

for sub in Left-Accumbens-area Left-Amygdala Left-Caudate Left-Hippocampus Left-Pallidum \
         Left-Putamen Left-Thalamus-Proper Right-Accumbens-area Right-Amygdala \
         Right-Caudate Right-Hippocampus Right-Pallidum Right-Putamen Right-Thalamus-Proper; do
    printf "%g," `grep  ${sub} ${dataDir}/aseg.stats | awk '{print $4}'` >> \
         "${outDir}/${idBIDS}_sctx_volume.csv"
done
    
printf "%g" `cat ${dataDir}/aseg.stats | grep IntraCranialVol | awk -F, '{print $4}'` >> \
     "${outDir}/${idBIDS}_sctx_volume.csv"
echo "" >> "${outDir}/${idBIDS}_sctx_volume.csv"

if [[ -f "${outDir}/${idBIDS}_sctx_volume.csv" ]]; then ((Nsteps++)); fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 1 ]; then status="COMPLETED"; else status="ERROR sctx_vol is missing a processing step"; fi
Title "sctx_vol processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/01
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/sctx_vol_*.txt)"
echo "${id}, ${SES/ses-/}, sctx_vol, $status N=$(printf "%02d" "$Nsteps")/08, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
