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
outDir="${out//micapipe/}/z-brains/scene-nativepro/sub-${id}/${SES}"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir="${dir_freesurfer}/stats/"

#------------------------------------------------------------------------------#
### subcortical volume ###

# Create file header
echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal,ICV" > \
     "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-volume.csv"
printf "%s,"  "${idBIDS}" >> "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-volume.csv"

for sub in Left-Accumbens-area Left-Amygdala Left-Caudate Left-Hippocampus Left-Pallidum \
         Left-Putamen Left-Thalamus-Proper Right-Accumbens-area Right-Amygdala \
         Right-Caudate Right-Hippocampus Right-Pallidum Right-Putamen Right-Thalamus-Proper; do
    printf "%g," `grep  ${sub} ${dataDir}/aseg.stats | awk '{print $4}'` >> \
         "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-volume.csv"
done
    
printf "%g" `cat ${dataDir}/aseg.stats | grep IntraCranialVol | awk -F, '{print $4}'` >> \
     "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-volume.csv"
echo "" >> "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-volume.csv"

if [[ -f "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-volume.csv" ]]; then ((Nsteps++)); fi

#------------------------------------------------------------------------------#
### cortical (and hippocampal) thickness ###

# cortex
if [[ ! -f "$outDir/${idBIDS}_space-conte69_hemi-rh_desc-thickness_10mm.func.gii" ]]; then
    for hemi in lh rh; do

        Do_cmd mris_convert -c \
            ${proc_struct}/surfaces/morphology/${idBIDS}_space-conte69-32k_desc-${hemi}_thickness_10mm.mgh \
            ${proc_struct}/surfaces/conte69/${idBIDS}_space-conte69-32k_desc-${hemi}_midthickness.surf.gii \
            $outDir/${idBIDS}_space-conte69_hemi-${hemi}_desc-thickness_10mm.func.gii

        if [[ -f "$outDir/${idBIDS}_space-conte69_hemi-${hemi}_desc-thickness_10mm.func.gii" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${id} thickness to nativepro"; Nsteps=$((Nsteps + 2))
fi

# hippocampal subfields
dir_hip="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/sub-${id}/"
if [[ ! -f "$outDir/${idBIDS}_space-hipp_hemi-rh_desc-thickness_2mm.func.gii" ]]; then
    for hemi in lh rh; do
        [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])
        Do_cmd wb_command -metric-smoothing \
                          $outDir/${idBIDS}_space-nativepro_desc-hipp_hemi-${hemi}_midthickness.surf.gii \
                          $dir_hip/surf/sub-${id}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_thickness.shape.gii \
                          2 \
                          $outDir/${idBIDS}_space-hipp_hemi-${hemi}_desc-thickness_2mm.func.gii
        
        if [[ -f "$outDir/${idBIDS}_space-hipp_hemi-${hemi}_desc-thickness_2mm.func.gii" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${idBIDS} hippocampal thickness to nativepro"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 05 ]; then status="COMPLETED"; else status="ERROR sctx_vol is missing a processing step"; fi
Title "sctx_vol processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/05
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/proc_morph_*.txt)"
echo "${id}, ${SES/ses-/}, sctx_vol, $status N=$(printf "%02d" "$Nsteps")/05, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"


