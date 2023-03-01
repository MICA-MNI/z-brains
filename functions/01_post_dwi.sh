#!/bin/bash
#
# DWI post-processing:
#
# Generates vertexwise DWI FA and ADC (native, fsa5, and conte69) outputs:
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
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check inputs: Freesurfer space T1
if [ ! -f "$T1freesurfr" ]; then Error "T1 in freesurfer space not found for Subject $id : <SUBJECTS_DIR>/${id}/mri/T1.mgz"; exit; fi

# Check inputs: DWI FA and ADC maps
dti_FA="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-FA.mif"
dti_ADC="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-ADC.mif"
if [ ! -f "$dti_FA" ]; then Error "FA not found for Subject $id : ${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-FA.mif"; exit; fi
if [ ! -f "$dti_ADC" ]; then Error "ADC not found for Subject $id : ${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-ADC.mif"; exit; fi

#------------------------------------------------------------------------------#
Title "DWI FA and ADC\n\t\tmicapipe-z $Version, $PROC"
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
tmp="${tmpDir}/${RANDOM}_micapipe-z_postdwi_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${out//micapipe/}/analysis/scene-nativepro/${idBIDS}"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir="${dir_freesurfer}/surf"

#------------------------------------------------------------------------------#
# Register DWI space to T1fsnative
Note "Out directory" "$outDir"
# First set all the variables
mat_dwi_affine="${dir_warp}/${idBIDS}_space-dwi_from-dwi_to-nativepro_mode-image_desc-affine_0GenericAffine.mat"
dwi_SyN_str="${dir_warp}/${idBIDS}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_"
dwi_SyN_warp="${dwi_SyN_str}1Warp.nii.gz"
dwi_SyN_Invwarp="${dwi_SyN_str}1InverseWarp.nii.gz"
dwi_SyN_affine="${dwi_SyN_str}0GenericAffine.mat"
fod_wmN="${proc_dwi}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.mif"
fod="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"
T1nativepro_in_dwi=$(ls "${proc_dwi}/${idBIDS}_space-dwi_desc-t1w_nativepro_"*".nii.gz")
mode=$(echo ${T1nativepro_in_dwi/.nii.gz/} | awk -F 'dwi_desc-t1w_nativepro_' '{print $2}')
Do_cmd mrconvert -coord 3 0 "$fod_wmN" "$fod"

# Set the transformation type
Info "Transformation is ${mode}"
if [[ ${mode}  == "SyN" ]]; then
  trans_T12dwi="-t ${dwi_SyN_warp} -t ${dwi_SyN_affine} -t [${mat_dwi_affine},1]" # T1nativepro to DWI
  trans_dwi2T1="-t ${mat_dwi_affine} -t [${dwi_SyN_affine},1] -t ${dwi_SyN_Invwarp}"  # DWI to T1nativepro
elif [[ ${mode}  == "Affine" ]]; then
  trans_T12dwi="-t [${mat_dwi_affine},1]"
  trans_dwi2T1="-t ${mat_dwi_affine}"
fi


function transform_dwi2nativepro() {
  input=$1
  output="${outDir}/${idBIDS}_space-nativepro_model-DTI_map-${2}.nii.gz"
  # Apply registration
  if [[ ! -f "${output}" ]]; then
    Do_cmd antsApplyTransforms -d 3 -i "${input}" \
                        -r ${outDir}/${idBIDS}_space-nativepro_t1w.nii.gz \
                        "$trans_dwi2T1" \
                        -n NearestNeighbor \
                        -o "$output" -v
  if [[ -f "${output}" ]]; then ((Nsteps++)); fi
  else
    Info "Subject ${id} DWI ${2} is registered to fsspace"; ((Nsteps++))
  fi
}
# Create niftis if they don't exist
[[ ! -f ${dti_FA/mif/nii.gz} ]] && mrconvert -force "$dti_FA" ${dti_FA/mif/nii.gz}
[[ ! -f ${dti_ADC/mif/nii.gz} ]] && mrconvert -force "$dti_ADC" ${dti_ADC/mif/nii.gz}

# Apply transformations to a NIFTI file!!!
transform_dwi2nativepro ${dti_FA/mif/nii.gz} "FA"
transform_dwi2nativepro ${dti_ADC/mif/nii.gz} "ADC"

Info "Map to surface and apply 10mm smooth"
if [[ ! -f "$outDir/${idBIDS}_space-conte69_hemi-rh_midthickness_desc-${map}_10mm.func.gii" ]]; then
    for map in FA ADC; do
        for hemi in lh rh; do
            # Volume to surface
            Do_cmd wb_command -volume-to-surface-mapping \
                                 $outDir/${idBIDS}_space-nativepro_model-DTI_map-${map}.nii.gz \
                                 $outDir/${idBIDS}_space-nativepro_desc-conte69_hemi-${hemi}_midthickness.surf.gii \
                                 $outDir/${idBIDS}_space-conte69_hemi-${hemi}_midthickness_desc-${map}.func.gii \
                                 -trilinear

            Do_cmd wb_command -metric-smoothing \
                                 $outDir/${idBIDS}_space-nativepro_desc-conte69_hemi-${hemi}_midthickness.surf.gii \
                                 $outDir/${idBIDS}_space-conte69_hemi-${hemi}_midthickness_desc-${map}.func.gii \
                                 10 \
                                 $outDir/${idBIDS}_space-conte69_hemi-${hemi}_midthickness_desc-${map}_10mm.func.gii \

            if [[ -f "$outDir/${idBIDS}_space-conte69_hemi-${hemi}_midthickness_desc-${map}_10mm.func.gii" ]]; then ((Nsteps++)); fi
        done
    done
else
    Info "Subject ${id} DWI FA and ADC are registered to nativepro"; Nsteps=$((Nsteps + 4))
fi

Info "Map FA/MD to subcortical structures"
if [[ ! -f "${outDir}/${idBIDS}_subcortical-ADC.csv" ]]; then

    for map in FA ADC; do
        echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal" > \
                "${outDir}/${idBIDS}_subcortical-${map}.csv"
        printf "%s,"  "${idBIDS}" >> "${outDir}/${idBIDS}_subcortical-${map}.csv"

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

            # Get FA/MD for subcortical mask
            Do_cmd fslmaths $outDir/${idBIDS}_space-nativepro_model-DTI_map-${map}.nii.gz \
                            -mul "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" \
                            "${tmp}/${idBIDS}_${sctxname}_masked-dwi.nii.gz"

            # Input values in .csv file
            printf "%g," `fslstats "${tmp}/${idBIDS}_${sctxname}_masked-dwi.nii.gz" -M` >> \
                "${outDir}/${idBIDS}_subcortical-${map}.csv"
            if [[ -f "${outDir}/${idBIDS}_subcortical-${map}.csv" ]]; then ((Nsteps++)); fi
        done
        echo "" >> "${outDir}/${idBIDS}_subcortical-${map}.csv"
    done
else
    Info "Subject ${idBIDS} FA/ADC are mapped to subcortical structures"; Nsteps=$((Nsteps + 28))
fi

Info "Map FA/MD to hippocampal subfields"
dir_hip="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/sub-${id}"
#dir_hip="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/${idBIDS/_//}"
affine_xfm="${tmp}/${idBIDS}_from-nativepro_to_dwi_mode-image_desc-affine_0GenericAffine.mat"

if [[ ! -f "$outDir/${idBIDS}_space-hipp_hemi-rh_midthickness_desc-ADC_2mm.func.gii" ]]; then
    # Output FA/MD
    for map in FA ADC; do
    [[ ${map} == "FA" ]] && map2=${dti_FA/.mif/}.nii.gz || map2=${dti_ADC/.mif/}.nii.gz

        for hemi in lh rh; do
            Do_cmd wb_command -volume-to-surface-mapping $outDir/${idBIDS}_space-nativepro_model-DTI_map-${map}.nii.gz \ \
                              $outDir/${idBIDS}_space-nativepro_desc-hipp_hemi-${hemi}_midthickness.surf.gii \
                              "${tmp}/${idBIDS}_hemi-${hemi}_${map}_den-0p5mm_label-hipp_midthickness.func.gii" \
			                  -trilinear

            Do_cmd wb_command -metric-smoothing \
                               $outDir/${idBIDS}_space-nativepro_desc-hipp_hemi-${hemi}_midthickness.surf.gii \
                              "${tmp}/${idBIDS}_hemi-${hemi}_${map}_den-0p5mm_label-hipp_midthickness.func.gii" \
                              2 \
                              $outDir/${idBIDS}_space-hipp_hemi-${hemi}_midthickness_desc-${map}_2mm.func.gii
            
            if [[ -f "$outDir/${idBIDS}_space-hipp_hemi-${hemi}_midthickness_desc-${map}_2mm.func.gii" ]]; then ((Nsteps++)); fi

        done
    done
else
    Info "Subject ${idBIDS} FA/MD are mapped to hippocampal subfields"; Nsteps=$((Nsteps + 4))
fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 42 ]; then status="COMPLETED"; else status="ERROR POST-DWI is missing a processing step"; fi
Title "post-dwi processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/42
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/post_dwi_*.txt)"
echo "${id}, ${SES/ses-/}, POST-DWI, $status N=$(printf "%02d" "$Nsteps")/42, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
