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
outDir="${proc_dwi}/surfaces/"
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
T1fsnative=${proc_struct}/${idBIDS}_space-fsnative_t1w.nii.gz
T1fsnative_brain=${tmp}/${idBIDS}_space-fsspace_T1w_brain.nii.gz
fod="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"
T1nativepro_brain=${proc_struct}/${idBIDS}_space-nativepro_t1w_brain.nii.gz
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

Info "Registering T1nativepro mask to fsnative_space"
Do_cmd antsApplyTransforms -d 3 -e 3 -i ${T1nativepro_brain} -r "${T1fsnative}" \
    -t ["${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_t1w_0GenericAffine.mat",1] \
    -o ${T1fsnative_brain} -v -u int

function transform_dwi2fsnative() {
  input=$1
  output="${tmp}/${idBIDS}_space-fsspace_model-DTI_map-${2}.nii.gz"
  # Apply registration
  if [[ ! -f "${output}" ]]; then
    Do_cmd antsApplyTransforms -d 3 -i "${input}" \
                        -r "$T1fsnative_brain" \
                        -t ["${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_t1w_0GenericAffine.mat",1] \
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
transform_dwi2fsnative ${dti_FA/mif/nii.gz} "FA"
transform_dwi2fsnative ${dti_ADC/mif/nii.gz} "ADC"

Info "Map to surface (native) and register to fsa5 and apply 10mm smooth"
if [[ ! -f "${outDir}/${idBIDS}_space-fsaverage5_desc-rh_model-DTI_map-ADC_10mm.mgh" ]]; then
    for map in FA ADC; do
        for hemi in lh rh; do
            # Volume to surface
            Do_cmd mri_vol2surf --mov "${tmp}/${idBIDS}_space-fsspace_model-DTI_map-${map}.nii.gz" \
                                --regheader "$idBIDS" \
                                --surf white \
                                --interp trilinear --hemi "$hemi" \
                                --out "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_model-DTI_map-${map}.mgh"

            Do_cmd mri_surf2surf --hemi "$hemi" \
                                --srcsubject "$idBIDS" \
                                --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_model-DTI_map-${map}.mgh" \
                                --trgsubject fsaverage5 \
                                --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_model-DTI_map-${map}.mgh"

            Do_cmd mri_surf2surf --hemi "$hemi" \
                                --fwhm-trg 10 \
                                --srcsubject "$idBIDS" \
                                --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_model-DTI_map-${map}.mgh" \
                                --trgsubject fsaverage5 \
                                --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_model-DTI_map-${map}_10mm.mgh"

            if [[ -f "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_model-DTI_map-${map}_10mm.mgh" ]]; then ((Nsteps++)); fi
        done
    done
else
    Info "Subject ${id} DWI FA and ADC are registered to fsa5"; Nsteps=$((Nsteps + 4))
fi

Info "Register to conte69 and apply 10mm smooth"
if [[ ! -f "${outDir}/${idBIDS}_space-conte69-32k_desc-rh_model-DTI_map-ADC_10mm.mgh" ]]; then
    for map in FA ADC; do
        for hemi in lh rh; do
            [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
            HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

            Do_cmd mri_convert "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_model-DTI_map-${map}.mgh" \
                               "${tmp}/${idBIDS}_space-fsnative_desc-${hemi}_model-DTI_map-${map}.func.gii"

            Do_cmd wb_command -metric-resample \
                "${tmp}/${idBIDS}_space-fsnative_desc-${hemi}_model-DTI_map-${map}.func.gii" \
                "${dir_conte69}/${idBIDS}_${hemi}_sphereReg.surf.gii" \
                "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
                ADAP_BARY_AREA \
                "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_model-DTI_map-${map}.func.gii" \
                -area-surfs \
                "${dir_freesurfer}/surf/${hemi}.midthickness.surf.gii" \
                "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_midthickness.surf.gii"

            Do_cmd mri_convert "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_model-DTI_map-${map}.func.gii" \
                               "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_model-DTI_map-${map}.mgh"

            # Smoothing
            Do_cmd wb_command -metric-smoothing \
                "${util_surface}/fsaverage.${HEMICAP}.midthickness_orig.32k_fs_LR.surf.gii" \
                "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_model-DTI_map-${map}.func.gii" \
                10 \
                "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_model-DTI_map-${map}_10mm.func.gii"

            Do_cmd mri_convert "${tmp}/${idBIDS}_space-conte69-32k_desc-${hemi}_model-DTI_map-${map}_10mm.func.gii" \
                               "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_model-DTI_map-${map}_10mm.mgh"

            if [[ -f "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_model-DTI_map-${map}_10mm.mgh" ]]; then ((Nsteps++)); fi
        done
    done
else
    Info "Subject ${idBIDS} DWI FA and ADC maps are registered to conte69"; Nsteps=$((Nsteps + 4))
fi

Info "Map FA/MD to subcortical structures"
if [[ ! -f "${outDir}/${idBIDS}_space-dwi_subcortical-ADC.csv" ]]; then
    # Transform to dwi space
    Do_cmd antsApplyTransforms -d 3 -e 3 \
                               -i "${dir_freesurfer}/mri/aseg.mgz" \
                               -r "${fod}" \
                               -t "${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_t1w_0GenericAffine.mat" \
                               "${trans_T12dwi}" \
                               -o "${tmp}/${idBIDS}_space-dwi_atlas-subcortical.nii.gz" \
                               -v -u float --float --interpolation NearestNeighbor

    for map in FA ADC; do
        echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal" > \
                "${outDir}/${idBIDS}_space-dwi_subcortical-${map}.csv"
        printf "%s,"  "${idBIDS}" >> "${outDir}/${idBIDS}_space-dwi_subcortical-${map}.csv"

        for sub in 26 18 11 17 13 12 10 58 54 50 53 52 51 49; do
            if [[ ${sub} == 26 ]]; then sctxname="Left-Accumbens-area"; elif [[ ${sub} == 18 ]]; then sctxname="Left-Amygdala"; \
            elif [[ ${sub} == 11 ]]; then sctxname="Left-Caudate"; elif [[ ${sub} == 17 ]]; then sctxname="Left-Hippocampus"; \
            elif [[ ${sub} == 13 ]]; then sctxname="Left-Pallidum"; elif [[ ${sub} == 12 ]]; then sctxname="Left-Putamen"; \
            elif [[ ${sub} == 10 ]]; then sctxname="Left-Thalamus-Proper"; elif [[ ${sub} == 58 ]]; then sctxname="Right-Accumbens-area"; \
            elif [[ ${sub} == 54 ]]; then sctxname="Right-Amygdala"; elif [[ ${sub} == 50 ]]; then sctxname="Right-Caudate"; \
            elif [[ ${sub} == 53 ]]; then sctxname="Right-Hippocampus"; elif [[ ${sub} == 52 ]]; then sctxname="Right-Pallidum"; \
            elif [[ ${sub} == 51 ]]; then sctxname="Right-Putamen"; elif [[ ${sub} == 49 ]]; then sctxname="Right-Thalamus-Proper"; fi

            # Extract subcortical masks
            Do_cmd mri_binarize --i "${tmp}/${idBIDS}_space-dwi_atlas-subcortical.nii.gz" \
                                --match "${sub}" \
                                --o "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz"

            # Get FA/MD for subcortical mask
            [[ ${map} == "FA" ]] && map2=${dti_FA/.mif/}.nii.gz || map2=${dti_ADC/.mif/}.nii.gz
            Do_cmd fslmaths ${map2} \
                            -mul "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" \
                            "${tmp}/${idBIDS}_${sctxname}_masked-dwi.nii.gz"

            # Input values in .csv file
            printf "%g," `fslstats "${tmp}/${idBIDS}_${sctxname}_masked-dwi.nii.gz" -M` >> \
                "${outDir}/${idBIDS}_space-dwi_subcortical-${map}.csv"
            if [[ -f "${outDir}/${idBIDS}_space-dwi_subcortical-${map}.csv" ]]; then ((Nsteps++)); fi
        done
        echo "" >> "${outDir}/${idBIDS}_space-dwi_subcortical-${map}.csv"
    done
else
    Info "Subject ${idBIDS} FA/ADC are mapped to subcortical structures"; Nsteps=$((Nsteps + 28))
fi

Info "Map FA/MD to hippocampal subfields"
dir_hip="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/sub-${id}"
#dir_hip="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/${idBIDS/_//}"
affine_xfm="${tmp}/${idBIDS}_from-nativepro_to_dwi_mode-image_desc-affine_0GenericAffine.mat"

if [[ ! -f "${outDir}/${idBIDS}_hemi-rh_space-dwi_desc-dwi-ADC_den-0p5mm_label-hipp_midthickness_10mm.func.gii" ]]; then
    if [[ ${mode} == "Affine" ]]; then
    # Select the transformation file from ANTs
      Do_cmd ConvertTransformFile 3 "${mat_dwi_affine}" "${mat_dwi_affine/.mat/}.txt"
      Do_cmd lta_convert --initk "${mat_dwi_affine/.mat/}.txt" \
                         --outitk "${tmp}/${idBIDS}_from-nativepro_to_dwi_mode-image_desc-affine_0GenericAffine.txt" \
              		       --src "${proc_struct}/${idBIDS}_space-nativepro_t1w_brain.nii.gz" \
              		       --trg "${fod}" \
                         --invert

      Do_cmd lta_convert --initk "${tmp}/${idBIDS}_from-nativepro_to_dwi_mode-image_desc-affine_0GenericAffine.txt" \
                         --outfsl "${affine_xfm}" \
                         --src "${proc_struct}/${idBIDS}_space-nativepro_t1w_brain.nii.gz" \
                         --trg "${fod}"

    elif [[ ${mode} == "SyN" ]]; then
      Info "Convert ${mode} affine transformation"
      SyN_affine="${tmp}/CombinedAffine.mat"
      Do_cmd antsApplyTransforms -v -d 3 -o Linear["${tmp}/CombinedAffine.mat",0] -t "${dwi_SyN_affine}" -t ["${mat_dwi_affine}",1] -r "${fod}"

      Do_cmd ConvertTransformFile 3 "${SyN_affine}" "${SyN_affine/.mat/}.txt"
      Do_cmd lta_convert --initk "${SyN_affine/.mat/}.txt" \
                         --outitk "${tmp}/${idBIDS}_from-nativepro_to_dwi_mode-image_desc-affine_0GenericAffine.txt" \
                         --src "${proc_struct}/${idBIDS}_space-nativepro_t1w_brain.nii.gz" \
                         --trg "${fod}"

      Do_cmd lta_convert --initk "${tmp}/${idBIDS}_from-nativepro_to_dwi_mode-image_desc-affine_0GenericAffine.txt" \
                         --outfsl "${affine_xfm}" \
                         --src "${proc_struct}/${idBIDS}_space-nativepro_t1w_brain.nii.gz" \
                         --trg "${fod}"
    fi

    # Output FA/MD
    for map in FA ADC; do
    [[ ${map} == "FA" ]] && map2=${dti_FA/.mif/}.nii.gz || map2=${dti_ADC/.mif/}.nii.gz

        for hemi in lh rh; do
            [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
            HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])
            #hip_surf="${dir_hip}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" 
            hip_surf="${dir_hip}/surf/sub-${id}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" # Hippunfold surface
            hippsurf_dwi="${outDir}/${idBIDS}_hemi-${HEMICAP}_space-dwi_den-0p5mm_label-hipp_midthickness.surf.gii" # Hippocampal surface in DWI space
            if [[ ${mode} == "SyN" ]]; then
              Info "Apply ${mode} transformations to ${HEMICAP} surface"

              # Affine transformation from T1nativepro to DWI
              Do_cmd wb_command -surface-apply-affine "${hip_surf}" \
                              "${affine_xfm}" \
                              "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-dwi_den-0p5mm_label-hipp_midthickness_tmp.surf.gii" \
                              -flirt "${proc_struct}/${idBIDS}_space-nativepro_t1w.nii.gz" \
                              "${fod}"
              # Warp transformation from spcae-dwi_T1nativepro to DWI
              world_warp=${tmp}/from-t1nativepro_space-dwi_to-dwi_world_warpfield.nii.gz
              Do_cmd wb_command -convert-warpfield -from-itk "${dwi_SyN_Invwarp}" -to-world ${world_warp}
              Do_cmd wb_command -surface-apply-warpfield "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-dwi_den-0p5mm_label-hipp_midthickness_tmp.surf.gii" \
                              ${world_warp} \
                              "${hippsurf_dwi}"
            elif [[ ${mode} == "Affine" ]]; then
              Info "Apply ${mode} transformations to surface"
              #statements
              Do_cmd wb_command -surface-apply-affine "${hip_surf}" \
                              "${affine_xfm}" \
                              "${hippsurf_dwi}" \
                              -flirt "${proc_struct}/${idBIDS}_space-nativepro_t1w.nii.gz" \
                              "${fod}"
            fi
            map_on_surf="${outDir}/${idBIDS}_hemi-${hemi}_space-dwi_desc-dwi-${map}_den-0p5mm_label-hipp_midthickness.func.gii"
            Do_cmd wb_command -volume-to-surface-mapping "${map2}" \
                            "${hippsurf_dwi}" \
                            "${map_on_surf}" \
                            -trilinear

            Do_cmd wb_command -metric-smoothing "${hippsurf_dwi}" \
                            "${map_on_surf}" \
                            2 "${outDir}/${idBIDS}_hemi-${hemi}_space-dwi_desc-dwi-${map}_den-0p5mm_label-hipp_midthickness_2mm.func.gii"

            if [[ -f "${outDir}/${idBIDS}_hemi-${hemi}_space-dwi_desc-dwi-${map}_den-0p5mm_label-hipp_midthickness_2mm.func.gii" ]]; then ((Nsteps++)); fi
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
