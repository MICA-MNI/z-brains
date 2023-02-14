#!/bin/bash
#
# T2-FLAIR processing:
#
# Generates vertexwise T2-FLAIR intensities (native, fsa5, and conte69) outputs:
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


#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe
    source ${MICAPIPE}/functions/init.sh;
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check inputs: Freesurfer space T1
if [ ! -f "$T1freesurfr" ]; then Error "T1 in freesurfer space not found for Subject $id : <SUBJECTS_DIR>/${id}/mri/T1.mgz"; exit; fi

# Check inputs: T2-FLAIR
if [ ! -f "$bids_flair" ]; then Error "T2-flair not found for Subject $id : ${subject_bids}/anat/${idBIDS}*FLAIR.nii.gz"; exit; fi


#------------------------------------------------------------------------------#
Title "T2-FLAIR intensities\n\t\tmicapipe-z $Version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Saving temporal dir: $nocleanup"

# Timer
aloita=$(date +%s)
Nsteps=0

# Freesurfer SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe-z_flair_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${proc_struct}/surfaces/flair"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir="${dir_freesurfer}/surf"


#------------------------------------------------------------------------------#
### FLAIR intensity correction ###
mkdir "${proc_struct}/flair"

# Bias field correction
flair_N4="${proc_struct}/flair/${idBIDS}_space-flair_desc-flair_N4.nii.gz"
if [[ ! -f "$flair_N4" ]]; then
    Do_cmd N4BiasFieldCorrection -d 3 -i "$bids_flair" -r \
                                -o "$flair_N4" -v
else
    Info "Subject ${id} T2-FLAIR is N4 bias corrected"; Nsteps=$((Nsteps + 1))
fi

# Clamp and rescale intensities
flair_clamp="${proc_struct}/flair/${idBIDS}_space-flair_desc-flair_N4_clamp.nii.gz"
flair_rescale="${proc_struct}/flair/${idBIDS}_space-flair_desc-flair_N4_rescale.nii.gz"
if [[ ! -f "$flair_rescale" ]]; then
    # Clamp intensities
    Do_cmd ImageMath 3 "$flair_clamp" TruncateImageIntensity "$flair_N4" 0.01 0.99 75

    # Rescale intensity [0,100]
    Do_cmd ImageMath 3 "$flair_rescale" RescaleImage "$flair_clamp" 0 100
else
    Info "Subject ${id} T2-FLAIR is intensity corrected"; Nsteps=$((Nsteps + 1))
fi

# Normalize intensities by GM/WM interface, uses 5ttgen
flair_preproc="${proc_struct}/flair/${idBIDS}_space-flair_desc-flair_preproc.nii.gz"
if [[ ! -f "$flair_gmwmi" ]]; then
    # Get gm/wm interface mask
    t1_gmwmi="${tmp}/${idBIDS}_space-nativepro_desc-gmwmi-mask.nii.gz"
    t1_5tt="${tmp}/${idBIDS}_space-nativepro_t1w_5TT.nii.gz"
    if [[ ! -f "$t1_gmwmi" ]]; then
        Info "Calculating Gray matter White matter interface mask"
        Do_cmd 5tt2gmwmi "$t1_5tt" "$t1_gmwmi"; ((Nsteps++))
    else
        Info "Subject ${id} has Gray matter White matter interface mask"; ((Nsteps++))
    fi

    # Register nativepro and flair
    str_flair_affine="${dir_warp}/${idBIDS}_from-flair_to-nativepro_mode-image_desc-affine_"
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$T1nativepro_brain" -m "$flair_rescale" -o "$str_flair_affine" -t a -n "$threads" -p d
    t1_gmwmi_in_flair="${tmp}/flair/${idBIDS}_space-flair_desc-gmwmi-mask.nii.gz"
    Do_cmd antsApplyTransforms -d 3 -i "$t1_gmwmi" -r "$flair_rescale" -t ["$str_flair_affine"0GenericAffine.mat,1] -o "$t1_gmwmi_in_flair" -v -u float

    # binarize mask
    t1_gmwmi_in_flair_thr="${tmp}/flair/${idBIDS}_space-flair_desc-gmwmi-thr.nii.gz"
    fslmaths "$t1_gmwmi_in_flair" -thr 0.5 -bin "$t1_gmwmi_in_flair_thr"

    # compute mean flair intensity in non-zero voxels
    gmwmi_mean=`fslstats "$flair_rescale" -M -k "$t1_gmwmi_in_flair_thr"`

    # Normalize flair
    fslmaths "$flair_rescale" -div $gmwmi_mean "$flair_preproc"
else
    Info "Subject ${id} T2-FLAIR is normalized by GM/WM interface"; Nsteps=$((Nsteps + 1))
fi


#------------------------------------------------------------------------------#
### FLAIR registrations ###
flair_fs="${proc_struct}/flair/${idBIDS}_space-fsspace_desc-flair_preproc.nii.gz"
flair2fs="${dir_warp}/${idBIDS}_from-flair_to-fsnative_mode-image_desc-flair.dat"
if [[ ! -f "$flair_fs" ]]; then

    # Register T2-FLAIR to fsspace
    Do_cmd bbregister --s "$idBIDS" \
                      --mov "$flair_preproc" \
                      --reg "$flair2fs" \
                      --init-coreg --t2 \
                      --o "$flair_fs"

    # Transform transform!
    Do_cmd lta_convert --inreg "$flair2fs" \
                       --outitk "${tmp}/${idBIDS}_from-flair_to-fsnative_mode-image_desc-flair.txt" \
                       --src "$flair_fs" \
                       --trg "${dir_freesurfer}/mri/orig.mgz"

    Do_cmd antsApplyTransforms -d 3 \
                               -t ["${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_t1w_0GenericAffine.mat",1] \
                               -t ["${tmp}/${idBIDS}_from-flair_to-fsnative_mode-image_desc-flair.txt",1] \
                               -r "$flair_fs" \
                               -o Linear["${dir_warp}/${idBIDS}_from-nativepro_t1w_to_flair_0GenericAffine.mat"] \
                               -v -u float --float
    
    Do_cmd ConvertTransformFile 3 "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_flair_0GenericAffine.mat" \
	                        "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_flair_0GenericAffine.txt"

    Do_cmd lta_convert --initk "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_flair_0GenericAffine.txt" \
                       --outfsl "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_flair_0GenericAffine.mat" \
                       --src "$T1nativepro" \
                       --trg "$flair_fs"

else
    Info "Subject ${id} T2-FLAIR is in fsspace"; Nsteps=$((Nsteps + 1))
fi


#------------------------------------------------------------------------------#
### Map intensities to cortex, subcortex, and hippocampus ###

# Map to surface (native) and register to fsa5 and apply 10mm smooth
if [[ ! -f "${outDir}/${idBIDS}_space-fsaverage5_desc-rh_flair_10mm.mgh" ]]; then
    for hemi in lh rh; do
        # Volume to surface    
        Do_cmd mri_vol2surf --mov "$flair_preproc" \
                            --reg "$flair2fs" \
                            --surf white --trgsubject "$idBIDS" \
                            --interp trilinear --hemi "$hemi" \
                            --out "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_flair.mgh"

        Do_cmd mri_surf2surf --hemi "$hemi" \
                             --srcsubject "$idBIDS" \
                             --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_flair.mgh" \
                             --trgsubject fsaverage5 \
                             --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_flair.mgh"

        Do_cmd mri_surf2surf --hemi "$hemi" \
                             --fwhm-trg 10 \
                             --srcsubject "$idBIDS" \
                             --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_flair.mgh" \
                             --trgsubject fsaverage5 \
                             --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_flair_10mm.mgh"

        if [[ -f "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_flair_10mm.mgh" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${id} T2-FLAIR is registered to fsa5"; Nsteps=$((Nsteps + 2))
fi

# Register to conte69 and apply 10mm smooth
if [[ ! -f "${outDir}/${idBIDS}_space-conte69-32k_desc-rh_flair_10mm.mgh" ]]; then
    for hemi in lh rh; do
        [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

        Do_cmd mri_convert "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_flair.mgh" \
                           "${tmp}/${idBIDS}_space-fsnative_desc-${hemi}_flair.func.gii"

        Do_cmd wb_command -metric-resample \
            "${tmp}/${idBIDS}_space-fsnative_desc-${hemi}_flair.func.gii" \
            "${dir_conte69}/${idBIDS}_${hemi}_sphereReg.surf.gii" \
            "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
            ADAP_BARY_AREA \
            "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_flair.func.gii" \
            -area-surfs \
            "${dir_freesurfer}/surf/${hemi}.midthickness.surf.gii" \
            "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_midthickness.surf.gii"

        Do_cmd mri_convert "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_flair.func.gii" \
                           "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_flair.mgh"

        # Smoothing
        Do_cmd wb_command -metric-smoothing \
            "${util_surface}/fsaverage.${HEMICAP}.midthickness_orig.32k_fs_LR.surf.gii" \
            "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_flair.func.gii" \
            10 \
            "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_flair_10mm.func.gii"

        Do_cmd mri_convert "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_flair_10mm.func.gii" \
                           "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_flair_10mm.mgh"
        
        if [[ -f "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_flair_10mm.mgh" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${idBIDS} T2-FLAIR is registered to conte69"; Nsteps=$((Nsteps + 2))
fi

# Map flair intensities to subcortical structures
if [[ ! -f "${outDir}/${idBIDS}_space-flair_subcortical-intensities.csv" ]]; then
    # Transform to flair space
    Do_cmd mri_vol2vol --mov "$flair_preproc" \
                       --targ "${dir_freesurfer}/mri/aseg.mgz" \
                       --reg "${dir_warp}/${idBIDS}_from-flair_to-fsnative_mode-image_desc-flair.dat" \
                       --o "${tmp}/${idBIDS}_space-flair_atlas-subcortical.nii.gz" \
                       --inv --nearest 
    
    echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal" > \
            "${outDir}/${idBIDS}_space-flair_subcortical-intensities.csv"
    printf "%s,"  "${idBIDS}" >> "${outDir}/${idBIDS}_space-flair_subcortical-intensities.csv"

    for sub in 26 18 11 17 13 12 10 58 54 50 53 52 51 49; do
        if [[ ${sub} == 26 ]]; then sctxname="Left-Accumbens-area"; elif [[ ${sub} == 18 ]]; then sctxname="Left-Amygdala"; \
        elif [[ ${sub} == 11 ]]; then sctxname="Left-Caudate"; elif [[ ${sub} == 17 ]]; then sctxname="Left-Hippocampus"; \
        elif [[ ${sub} == 13 ]]; then sctxname="Left-Pallidum"; elif [[ ${sub} == 12 ]]; then sctxname="Left-Putamen"; \
        elif [[ ${sub} == 10 ]]; then sctxname="Left-Thalamus-Proper"; elif [[ ${sub} == 58 ]]; then sctxname="Right-Accumbens-area"; \
        elif [[ ${sub} == 54 ]]; then sctxname="Right-Amygdala"; elif [[ ${sub} == 50 ]]; then sctxname="Right-Caudate"; \
        elif [[ ${sub} == 53 ]]; then sctxname="Right-Hippocampus"; elif [[ ${sub} == 52 ]]; then sctxname="Right-Pallidum"; \
        elif [[ ${sub} == 51 ]]; then sctxname="Right-Putamen"; elif [[ ${sub} == 49 ]]; then sctxname="Right-Thalamus-Proper"; fi

        # Extract subcortical masks
        Do_cmd mri_binarize --i "${tmp}/${idBIDS}_space-flair_atlas-subcortical.nii.gz" \
                            --match "${sub}" \
                            --o "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz"

        # Get flair intensities for subcortical mask
        Do_cmd fslmaths "$flair_preproc" \
                        -mul "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" \
                        "${tmp}/${idBIDS}_${sctxname}_masked-flair.nii.gz"

        # Input values in .csv file
        printf "%g," `fslstats "${tmp}/${idBIDS}_${sctxname}_masked-flair.nii.gz" -M` >> \
            "${outDir}/${idBIDS}_space-flair_subcortical-intensities.csv"
        if [[ -f "${outDir}/${idBIDS}_space-flair_subcortical-intensities.csv" ]]; then ((Nsteps++)); fi
    done
    echo "" >> "${outDir}/${idBIDS}_space-flair_subcortical-intensities.csv"
else
    Info "Subject ${idBIDS} T2-FLAIR is mapped to subcortical areas"; Nsteps=$((Nsteps + 14))
fi

# Map flair intensities to hippocampal subfields
dir_hip="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/sub-${id}/"
if [[ ! -f "${proc_struct}/surfaces/flair/${idBIDS}_hemi-rh_space-flair_desc-flair_N4_den-0p5mm_label-hipp_midthickness_2mm.func.gii" ]]; then

    for hemi in lh rh; do
        [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

        Do_cmd wb_command -surface-apply-affine "${dir_hip}/surf/sub-${id}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                          "${dir_warp}/${idBIDS}_from-nativepro_t1w_to_flair_0GenericAffine.mat" \
                          "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-flair_den-0p5mm_label-hipp_midthickness.surf.gii" \
                          -flirt "$T1nativepro" \
                          "$flair_preproc"

        Do_cmd wb_command -volume-to-surface-mapping "$flair_preproc" \
                          "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-flair_den-0p5mm_label-hipp_midthickness.surf.gii" \
                          "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-flair_desc-flair_N4_den-0p5mm_label-hipp_midthickness.func.gii" \
			              -trilinear

        Do_cmd wb_command -metric-smoothing "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-flair_den-0p5mm_label-hipp_midthickness.surf.gii" \
                          "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-flair_desc-flair_N4_den-0p5mm_label-hipp_midthickness.func.gii" \
                          2 "${proc_struct}/surfaces/flair/${idBIDS}_hemi-${hemi}_space-flair_desc-flair_N4_den-0p5mm_label-hipp_midthickness_2mm.func.gii" 
        
        if [[ -f "${proc_struct}/surfaces/flair/${idBIDS}_hemi-${hemi}_space-flair_desc-flair_N4_den-0p5mm_label-hipp_midthickness_2mm.func.gii" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${idBIDS} T2-FLAIR is mapped to hippocampus"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 21 ]; then status="COMPLETED"; else status="ERROR T2-FLAIR is missing a processing step"; fi
Title "proc-flair processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/21
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/proc_flair_*.txt)"
echo "${id}, ${SES/ses-/}, T2-FLAIR, $status N=$(printf "%02d" "$Nsteps")/21, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
