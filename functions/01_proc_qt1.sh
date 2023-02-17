#!/bin/bash
#
# qT1 processing:
#
# Generates vertexwise qT1 intensities (native, fsa5, and conte69) outputs:
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
outDir="${proc_struct}/surfaces/qt1"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir="${dir_freesurfer}/surf"

#------------------------------------------------------------------------------#
### qT1 intensity ###
mkdir "${proc_struct}/qt1"

# Mask t1map with inv2 image
if [[ ! -f "${proc_struct}/qt1/${idBIDS}_space-fsspace_desc-qt1.nii.gz" ]]; then
    Do_cmd bet "$bids_inv2" "${tmp}/${idBIDS}" -m -n
    Do_cmd fslmaths "${tmp}/${idBIDS}_mask.nii.gz" -dilD -dilD -dilD -dilD -dilD -dilD -bin "${tmp}/${idBIDS}_mask.nii.gz"
    Do_cmd fslmaths "$bids_qT1map" -mul "${tmp}/${idBIDS}_mask.nii.gz" "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz"

    # Register qT1 to fsspace
    Do_cmd bbregister --s "$idBIDS" \
                      --mov "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz" \
                      --reg "${dir_warp}/${idBIDS}_from-qt1_to-fsnative_mode-image_desc-qt1.dat" \
                      --init-coreg --t2 \
                      --o "${proc_struct}/qt1/${idBIDS}_space-fsspace_desc-qt1.nii.gz"
else
    Info "idBIDS ${id} qT1 is registered to fsspace"; Nsteps=$((Nsteps + 1))
fi

# Map to surface (native) and register to fsa5 and apply 10mm smooth
if [[ ! -f "${outDir}/${idBIDS}_space-fsaverage5_desc-rh_qt1_10mm.mgh" ]]; then
    for hemi in lh rh; do
        # Volume to surface    
        Do_cmd mri_vol2surf --mov "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz" \
                            --reg "${dir_warp}/${idBIDS}_from-qt1_to-fsnative_mode-image_desc-qt1.dat" \
                            --surf white --trgsubject "$idBIDS" \
                            --interp trilinear --hemi "$hemi" \
                            --out "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_qt1.mgh"

        Do_cmd mri_surf2surf --hemi "$hemi" \
                             --srcsubject "$idBIDS" \
                             --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_qt1.mgh" \
                             --trgsubject fsaverage5 \
                             --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_qt1.mgh"

        Do_cmd mri_surf2surf --hemi "$hemi" \
                             --fwhm-trg 10 \
                             --srcsubject "$idBIDS" \
                             --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_qt1.mgh" \
                             --trgsubject fsaverage5 \
                             --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_qt1_10mm.mgh"

        if [[ -f "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_qt1_10mm.mgh" ]]; then ((Nsteps++)); fi
    done
else
    Info "idBIDS ${id} qT1 is registered to fsa5"; Nsteps=$((Nsteps + 2))
fi

# Register to conte69 and apply 10mm smooth
if [[ ! -f "${outDir}/${idBIDS}_space-conte69-32k_desc-rh_qt1_10mm.mgh" ]]; then
    for hemi in lh rh; do
        [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

        Do_cmd mri_convert "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_qt1.mgh" \
                           "${tmp}/${idBIDS}_space-fsnative_desc-${hemi}_qt1.func.gii"

        Do_cmd wb_command -metric-resample \
            "${tmp}/${idBIDS}_space-fsnative_desc-${hemi}_qt1.func.gii" \
            "${dir_conte69}/${idBIDS}_${hemi}_sphereReg.surf.gii" \
            "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
            ADAP_BARY_AREA \
            "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_qt1.func.gii" \
            -area-surfs \
            "${dir_freesurfer}/surf/${hemi}.midthickness.surf.gii" \
            "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_midthickness.surf.gii"

        Do_cmd mri_convert "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_qt1.func.gii" \
                           "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_qt1.mgh"

        # Smoothing
        Do_cmd wb_command -metric-smoothing \
            "${util_surface}/fsaverage.${HEMICAP}.midthickness_orig.32k_fs_LR.surf.gii" \
            "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_qt1.func.gii" \
            10 \
            "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_qt1_10mm.func.gii"

        Do_cmd mri_convert "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_qt1_10mm.func.gii" \
                           "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_qt1_10mm.mgh"
        
        if [[ -f "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_qt1_10mm.mgh" ]]; then ((Nsteps++)); fi
    done
else
    Info "idBIDS ${idBIDS} qT1 is registered to conte69"; Nsteps=$((Nsteps + 2))
fi

# Map qt1 intensities to subcortical structures
if [[ ! -f "${outDir}/${idBIDS}_space-qt1_subcortical-intensities.csv" ]]; then
    # Transform to qt1 space
    Do_cmd mri_vol2vol --mov "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz" \
                       --targ "${dir_freesurfer}/mri/aseg.mgz" \
                       --reg "${dir_warp}/${idBIDS}_from-qt1_to-fsnative_mode-image_desc-qt1.dat" \
                       --o "${tmp}/${idBIDS}_space-qt1_atlas-subcortical.nii.gz" \
                       --inv --nearest 
    
    echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal" > \
            "${outDir}/${idBIDS}_space-qt1_subcortical-intensities.csv"
    printf "%s,"  "${idBIDS}" >> "${outDir}/${idBIDS}_space-qt1_subcortical-intensities.csv"

    for sub in 26 18 11 17 13 12 10 58 54 50 53 52 51 49; do
        if [[ ${sub} == 26 ]]; then sctxname="Left-Accumbens-area"; elif [[ ${sub} == 18 ]]; then sctxname="Left-Amygdala"; \
        elif [[ ${sub} == 11 ]]; then sctxname="Left-Caudate"; elif [[ ${sub} == 17 ]]; then sctxname="Left-Hippocampus"; \
        elif [[ ${sub} == 13 ]]; then sctxname="Left-Pallidum"; elif [[ ${sub} == 12 ]]; then sctxname="Left-Putamen"; \
        elif [[ ${sub} == 10 ]]; then sctxname="Left-Thalamus-Proper"; elif [[ ${sub} == 58 ]]; then sctxname="Right-Accumbens-area"; \
        elif [[ ${sub} == 54 ]]; then sctxname="Right-Amygdala"; elif [[ ${sub} == 50 ]]; then sctxname="Right-Caudate"; \
        elif [[ ${sub} == 53 ]]; then sctxname="Right-Hippocampus"; elif [[ ${sub} == 52 ]]; then sctxname="Right-Pallidum"; \
        elif [[ ${sub} == 51 ]]; then sctxname="Right-Putamen"; elif [[ ${sub} == 49 ]]; then sctxname="Right-Thalamus-Proper"; fi

        # Extract subcortical masks
        Do_cmd mri_binarize --i "${tmp}/${idBIDS}_space-qt1_atlas-subcortical.nii.gz" \
                            --match "${sub}" \
                            --o "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz"

        # Get qt1 intensities for subcortical mask
        Do_cmd fslmaths "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz" \
                        -mul "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" \
                        "${tmp}/${idBIDS}_${sctxname}_masked-qt1.nii.gz"

        # Input values in .csv file
        printf "%g," `fslstats "${tmp}/${idBIDS}_${sctxname}_masked-qt1.nii.gz" -M` >> \
            "${outDir}/${idBIDS}_space-qt1_subcortical-intensities.csv"
        if [[ -f "${outDir}/${idBIDS}_space-qt1_subcortical-intensities.csv" ]]; then ((Nsteps++)); fi
    done
    echo "" >> "${outDir}/${idBIDS}_space-qt1_subcortical-intensities.csv"
else
    Info "Subject ${idBIDS} qt1 is registered to conte69"; Nsteps=$((Nsteps + 14))
fi

# Map qt1 intensities to hippocampal subfields
dir_hip="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/sub-${id}/"
if [[ ! -f "${proc_struct}/surfaces/qt1/${idBIDS}_hemi-rh_space-qt1_desc-qt1_den-0p5mm_label-hipp_midthickness_10mm.func.gii" ]]; then
    # Transform transform!
    Do_cmd lta_convert --inreg "${dir_warp}/${idBIDS}_from-qt1_to-fsnative_mode-image_desc-qt1.dat" \
                       --outitk "${tmp}/${idBIDS}_from-qt1_to-fsnative_mode-image_desc-qt1.txt" \
                       --src "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz" \
                       --trg "${dir_freesurfer}/mri/orig.mgz"

    Do_cmd antsApplyTransforms -d 3 \
                               -t ["${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_t1w_0GenericAffine.mat",1] \
                               -t ["${tmp}/${idBIDS}_from-qt1_to-fsnative_mode-image_desc-qt1.txt",1] \
                               -r "${proc_struct}/qt2/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz" \
                               -o Linear["${dir_warp}/${idBIDS}_from-nativepro_t1w_to_qt1_0GenericAffine.mat"] \
                               -v -u float --float                   
    
    Do_cmd ConvertTransformFile 3 "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_qt1_0GenericAffine.mat" \
	                        "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_qt1_0GenericAffine.txt"

    Do_cmd lta_convert --initk "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_qt1_0GenericAffine.txt" \
                       --outfsl "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_qt1_0GenericAffine.mat" \
                       --src "${proc_struct}/${idBIDS}_space-nativepro_t1w.nii.gz" \
                       --trg "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz"

    for hemi in lh rh; do
        [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

        Do_cmd wb_command -surface-apply-affine "${dir_hip}/surf/sub-${id}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                          "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_qt1_0GenericAffine.mat" \
                          "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-qt1_den-0p5mm_label-hipp_midthickness.surf.gii" \
                          -flirt "${proc_struct}/${idBIDS}_space-nativepro_t1w.nii.gz" \
                          "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz"

        Do_cmd wb_command -volume-to-surface-mapping "${proc_struct}/qt1/${idBIDS}_space-qt1_desc-qt1_brain.nii.gz" \
                          "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-qt1_den-0p5mm_label-hipp_midthickness.surf.gii" \
                          "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-qt1_desc-qt1_den-0p5mm_label-hipp_midthickness.func.gii" \
			              -trilinear

        Do_cmd wb_command -metric-smoothing "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-qt1_den-0p5mm_label-hipp_midthickness.surf.gii" \
                          "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-qt1_desc-qt1_den-0p5mm_label-hipp_midthickness.func.gii" \
                          2 "${proc_struct}/surfaces/qt1/${idBIDS}_hemi-${hemi}_space-qt1_desc-qt1_den-0p5mm_label-hipp_midthickness_2mm.func.gii" 
        
        if [[ -f "${proc_struct}/surfaces/qt1/${idBIDS}_hemi-${hemi}_space-qt1_desc-qt1_den-0p5mm_label-hipp_midthickness_2mm.func.gii" ]]; then ((Nsteps++)); fi
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
if [ "$Nsteps" -eq 21 ]; then status="COMPLETED"; else status="ERROR qT1 is missing a processing step"; fi
Title "proc-qt1 processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/21
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/proc_qt1_*.txt)"
echo "${id}, ${SES/ses-/}, qT1, $status N=$(printf "%02d" "$Nsteps")/21, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
