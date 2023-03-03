#!/bin/bash
#
# qT1 processing:
#
# Generates vertexwise qT1 intensities (nativepro and conte69) outputs:
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

# Check inputs: Freesurfer space T1
if [ ! -f "$T1freesurfr" ]; then Error "T1 in freesurfer space not found for idBIDS $id : <idBIDSS_DIR>/${id}/mri/T1.mgz"; exit; fi

echo $bids_qT1map
echo $bids_inv2
# Check inputs: qT1 and inv2 images
if [ ! -f "$bids_qT1map" ] || [ ! -f "$bids_inv2" ]; then Error "qT1 or inv2 not found for idBIDS $id : ${idBIDS_bids}/anat/${idBIDS}*_T1map.nii.gz"; exit; fi

#------------------------------------------------------------------------------#
Title "qT1 intensities\n\t\tmicapipe-z $Version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Saving temporal dir: $nocleanup"

# Timer
aloita=$(date +%s)
Nsteps=0

# Freesurfer idBIDSs directory
export SUBJECTS_DIR=${dir_surf}

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe-van_qt1_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${out//micapipe/}/z-brains/scene-nativepro/sub-${id}/${SES}"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir="${dir_freesurfer}/surf"

#------------------------------------------------------------------------------#
### qT1 intensity ###

Info "Subject ${id} qt1 registration to nativepro with antsQuickRigid"
Do_cmd antsRegistrationSyNQuick.sh -d 3 -t r \
                -m ${outDir}/${idBIDS}_space-nativepro_t1w.nii.gz \
                -f ${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz \
                -o ${dir_warp}/${idBIDS}_from-nativepro_t1w_to_qt1_
if [[ -f "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_qt1_0GenericAffine.txt" ]]; then ((Nsteps++)); fi

mkdir "${proc_struct}/qt1"

# Mask t1map with inv2 image
if [[ ! -f "$outDir/${idBIDS}_space-nativepro_qt1.nii.gz" ]]; then
    Do_cmd bet "$bids_inv2" "${tmp}/${idBIDS}" -m -n
    Do_cmd fslmaths "${tmp}/${idBIDS}_mask.nii.gz" -dilD -dilD -dilD -dilD -dilD -dilD -bin "${tmp}/${idBIDS}_mask.nii.gz"
    Do_cmd fslmaths "$bids_qT1map" -mul "${tmp}/${idBIDS}_mask.nii.gz" "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz"

    # Transform qT1 to nativepro
    antsApplyTransforms -d 3 \
                -i ${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz \
                -r ${outDir}/${idBIDS}_space-nativepro_t1w.nii.gz \
                -t [${dir_warp}/${idBIDS}_from-nativepro_t1w_to_qt1_0GenericAffine.txt,1] \
                -o $outDir/${idBIDS}_space-nativepro_qt1.nii.gz
    if [[ -f "$outDir/${idBIDS}_space-nativepro_qt1.nii.gz" ]]; then ((Nsteps++)); fi
else
    Info "idBIDS ${id} qT1 is registered to nativepro"; Nsteps=$((Nsteps + 1))
fi

# Map to surface and apply 10mm smooth
if [[ ! -f "$outDir/${idBIDS}_space-conte69_hemi-rh_midthickness_desc-qt1_10mm.func.gii" ]]; then
    for hemi in lh rh; do
        # Volume to surface    
        Do_cmd wb_command -volume-to-surface-mapping \
                             $outDir/${idBIDS}_space-nativepro_qt1.nii.gz \
                             $outDir/${idBIDS}_space-nativepro_desc-conte69_hemi-${hemi}_midthickness.surf.gii \
                             $tmp/${idBIDS}_space-conte69_hemi-${hemi}_midthickness_desc-qt1.func.gii \
                             -trilinear

        Do_cmd wb_command -metric-smoothing \
                             $outDir/${idBIDS}_space-nativepro_desc-conte69_hemi-${hemi}_midthickness.surf.gii \
                             $tmp/${idBIDS}_space-conte69_hemi-${hemi}_midthickness_desc-qt1.func.gii \
                             10 \
                             $outDir/${idBIDS}_space-conte69_hemi-${hemi}_midthickness_desc-qt1_10mm.func.gii \

        if [[ -f "$outDir/${idBIDS}_space-conte69_hemi-${hemi}_midthickness_desc-qt1_10mm.func.gii" ]]; then ((Nsteps++)); fi
    done
else
    Info "idBIDS ${idBIDS} qT1 is registered to conte69"; Nsteps=$((Nsteps + 2))
fi

# Map qt1 intensities to subcortical structures
if [[ ! -f "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-qt1.csv" ]]; then
    
    echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal" > \
            "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-qt1.csv"
    printf "%s,"  "${idBIDS}" >> "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-qt1.csv"

    for sub in 26 18 11 17 13 12 10 58 54 50 53 52 51 49; do
        if [[ ${sub} == 26 ]]; then sctxname="Left-Accumbens-area"; elif [[ ${sub} == 18 ]]; then sctxname="Left-Amygdala"; \
        elif [[ ${sub} == 11 ]]; then sctxname="Left-Caudate"; elif [[ ${sub} == 17 ]]; then sctxname="Left-Hippocampus"; \
        elif [[ ${sub} == 13 ]]; then sctxname="Left-Pallidum"; elif [[ ${sub} == 12 ]]; then sctxname="Left-Putamen"; \
        elif [[ ${sub} == 10 ]]; then sctxname="Left-Thalamus-Proper"; elif [[ ${sub} == 58 ]]; then sctxname="Right-Accumbens-area"; \
        elif [[ ${sub} == 54 ]]; then sctxname="Right-Amygdala"; elif [[ ${sub} == 50 ]]; then sctxname="Right-Caudate"; \
        elif [[ ${sub} == 53 ]]; then sctxname="Right-Hippocampus"; elif [[ ${sub} == 52 ]]; then sctxname="Right-Pallidum"; \
        elif [[ ${sub} == 51 ]]; then sctxname="Right-Putamen"; elif [[ ${sub} == 49 ]]; then sctxname="Right-Thalamus-Proper"; fi

        # Extract subcortical masks
        Do_cmd mri_binarize --i ${outDir}/${idBIDS}_space-nativepro_aseg.nii.gz \
                            --match "${sub}" \
                            --o "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz"

        # Get qt1 intensities for subcortical mask
        Do_cmd fslmaths $outDir/${idBIDS}_space-nativepro_qt1.nii.gz \
                        -mul "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" \
                        "${tmp}/${idBIDS}_${sctxname}_masked-qt1.nii.gz"

        # Input values in .csv file
        printf "%g," `fslstats "${tmp}/${idBIDS}_${sctxname}_masked-qt1.nii.gz" -M` >> \
            "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-qt1.csv"
        if [[ -f "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-qt1.csv" ]]; then ((Nsteps++)); fi
    done
    echo "" >> "${outDir}/${idBIDS}_space-nativepro_desc-subcortical-qt1.csv"
else
    Info "Subject ${idBIDS} qt1 is registered to conte69"; Nsteps=$((Nsteps + 14))
fi

# Map qt1 intensities to hippocampal subfields
dir_hip="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/sub-${id}/"
if [[ ! -f "$outDir/${idBIDS}_space-hipp_hemi-rh_midthickness_desc-qt1_2mm.func.gii" ]]; then
    for hemi in lh rh; do
        Do_cmd wb_command -volume-to-surface-mapping $outDir/${idBIDS}_space-nativepro_qt1.nii.gz \
                          $outDir/${idBIDS}_space-nativepro_desc-hipp_hemi-${hemi}_midthickness.surf.gii \
                          "${tmp}/${idBIDS}_hemi-${hemi}_qt1_den-0p5mm_label-hipp_midthickness.func.gii" \
			              -trilinear

        Do_cmd wb_command -metric-smoothing \
                           $outDir/${idBIDS}_space-nativepro_desc-hipp_hemi-${hemi}_midthickness.surf.gii \
                          "${tmp}/${idBIDS}_hemi-${hemi}_qt1_den-0p5mm_label-hipp_midthickness.func.gii" \
                          2 \
                          $outDir/${idBIDS}_space-hipp_hemi-${hemi}_midthickness_desc-qt1_2mm.func.gii
        
        if [[ -f "$outDir/${idBIDS}_space-hipp_hemi-${hemi}_midthickness_desc-qt1_2mm.func.gii" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${idBIDS} qT1 is mapped to hippocampus"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 20 ]; then status="COMPLETED"; else status="ERROR qT1 is missing a processing step OR if an extra step was run then registration to nativepro was run de-novo"; fi
Title "proc-qt1 processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/20
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/proc_qt1_*.txt)"
echo "${id}, ${SES/ses-/}, qT1, $status N=$(printf "%02d" "$Nsteps")/20, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
