#!/bin/bash
#
# resampling of conte69 midthickness, pial, and white surfaces to space-nativepro:
#
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

#------------------------------------------------------------------------------#
Title "Surfaces sample to nativepro\n\t\tmicapipe-z $Version, $PROC"
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
tmp="${tmpDir}/${RANDOM}_micapipe-z_surfsamp_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${out//micapipe/}/analysis/scene-nativepro/${idBIDS}"
Do_cmd mkdir -p "$outDir"

# Data locations
#dir_freesurfer="${out//micapipe/}/freesurfer/${idBIDS}"
#dir_conte69="${out}/sub-${id}/${SES}/anat/surfaces/conte69/"

#------------------------------------------------------------------------------#
### resample each surface intensity correction ###

for hemi in lh rh
    do
    # remove offset-to-origin from any gifti surface derived from FS
    Do_cmd python "$ZBRAINS"/functions/removeFSoffset.py ${dir_freesurfer}/surf/${hemi}.midthickness.surf.gii $tmp/no_offset.surf.gii

    # now match target surface to the original fssurf bounding box
    Do_cmd wb_command -surface-match $tmp/no_offset.surf.gii ${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_midthickness.surf.gii $tmp/conte69_space-fsnative.surf.gii
    # note: this works because the conte69 surface still has (nearly) exactly the same shape as the fssurface, just different vertices. The more vertices get removed, the more the shape degrades and the exact bounding box may no longer match. With 32k vertices though, the bounding box stays the same to within <0.01mm

    # apply registration to nativepro
    Do_cmd c3d_affine_tool -itk ${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_t1w_0GenericAffine.mat -inv -o $tmp/aff_fsspacetonativepro.mat
    wb_command -surface-apply-affine $tmp/conte69_space-fsnative.surf.gii $tmp/aff_fsspacetonativepro.mat $outDir/${idBIDS}_space-nativepro_desc-conte69_hemi-${hemi}_midthickness.surf.gii

    ## also apply to pial and white
    Do_cmd mris_convert ${dir_freesurfer}/surf/${hemi}.pial $tmp/${hemi}_pial.surf.gii
    Do_cmd python "$ZBRAINS"/functions/removeFSoffset.py $tmp/${hemi}_pial.surf.gii $tmp/no_offset_pial.surf.gii 
    Do_cmd wb_command -surface-match $tmp/no_offset_pial.surf.gii ${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_pial.surf.gii $tmp/conte69_space-fsnative_pial.surf.gii
    Do_cmd wb_command -surface-apply-affine $tmp/conte69_space-fsnative_pial.surf.gii $tmp/aff_fsspacetonativepro.mat $outDir/${idBIDS}_space-nativepro_desc-conte69_hemi-${hemi}_pial.surf.gii

    Do_cmd mris_convert ${dir_freesurfer}/surf/${hemi}.white $tmp/${hemi}_white.surf.gii
    Do_cmd python "$ZBRAINS"/functions/removeFSoffset.py $tmp/${hemi}_white.surf.gii $tmp/no_offset_white.surf.gii 
    Do_cmd wb_command -surface-match $tmp/no_offset_white.surf.gii ${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_white.surf.gii $tmp/conte69_space-fsnative_white.surf.gii
    Do_cmd wb_command -surface-apply-affine $tmp/conte69_space-fsnative_white.surf.gii $tmp/aff_fsspacetonativepro.mat $outDir/${idBIDS}_space-nativepro_desc-conte69_hemi-${hemi}_white.surf.gii
    done

##### resample relevant images #####

Do_cmd cp ${proc_struct}/${idBIDS}_space-nativepro_t1w.nii.gz $outDir/${idBIDS}_space-nativepro_t1w.nii.gz
Do_cmd mri_convert /${dir_freesurfer}/mri/aseg.mgz $tmp/aseg.nii.gz
Do_cmd antsApplyTransforms -d 3 -n MultiLabel \
    -i $tmp/aseg.nii.gz \
    -r ${outDir}/${idBIDS}_space-nativepro_t1w.nii.gz \
    -t ${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_t1w_0GenericAffine.mat \
    -o ${outDir}/${idBIDS}_space-nativepro_aseg.nii.gz


##### hippocampal surfaces #####

dir_hip="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/sub-${id}"

for hemi in lh rh; do
    [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
    HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])
    Do_cmd  cp ${dir_hip}/surf/sub-${id}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii $outDir/${idBIDS}_space-nativepro_desc-hipp_hemi-${hemi}_midthickness.surf.gii
    done

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 27 ]; then status="COMPLETED"; else status="ERROR surfsamp is missing a processing step"; fi
Title "surfsamp processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/27
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/surfsamp_*.txt)"
echo "${id}, ${SES/ses-/}, surfsamp, $status N=$(printf "%02d" "$Nsteps")/21, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"


