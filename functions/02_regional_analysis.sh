#!/bin/bash
#
# Regional analysis: Generates a multivariate feature regional changes report
#
# This workflow makes use of many outputs and custom python scripts
#
# Atlases and templates are avaliable from:
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
thr=$8
demo=$9
PROC=${10}
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
Title "regional analysis\n\t\tmicapipe-z $Version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Saving temporal dir: $nocleanup"

# Timer
aloita=$(date +%s)
Nsteps=0

# Hippunfold directory
hipDir="${out/micapipe/}/hippunfold_v1.0.0/hippunfold/"

# Freesurfer SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Check CORTICAL input features (flair, qt1, thickness, adc); if missing a feature, skip this module
featList_ctx=()
if [ -f "${out//micapipe/}/z-brains/scene-nativepro/${idBIDS}_space-conte69_hemi-rh_midthickness_desc-flair_10mm.func.gii" ]; then
    featList_ctx+=("flair"); fi
if [ -f "${out//micapipe/}/z-brains/scene-nativepro/${idBIDS}_space-conte69_hemi-rh_midthickness_desc-qt1_10mm.func.gii" ]; then
    featList_ctx+=("qt1"); fi
if [ -f "${out//micapipe/}/z-brains/scene-nativepro/${idBIDS}_space-conte69_hemi-rh_midthickness_desc-ADC_10mm.func.gii" ]; then
    featList_ctx+=("adc"); fi
if [ -f "${out//micapipe/}/z-brains/scene-nativepro/${idBIDS}_space-conte69_hemi-rh_desc-thickness_10mm.func.gii" ]; then
    featList_ctx+=("thickness"); fi

# Check SUBCORTICAL input features (flair, qt1, thickness, adc); if missing a feature, skip this module
featList_sctx=()
if [ -f "${out//micapipe/}/${idBIDS}_subcortical-flair.csv" ]; then
    featList_sctx+=("flair"); fi
if [ -f "${out//micapipe/}/${idBIDS}_subcortical-qt1.csv" ]; then
    featList_sctx+=("qt1"); fi
if [ -f "${out//micapipe/}/${idBIDS}_subcortical-ADC.csv" ]; then
    featList_sctx+=("adc"); fi
if [ -f "${out//micapipe/}/${idBIDS}_subcortical-volume.csv" ]; then
    featList_sctx+=("thickness"); fi

# Check HIPPOCAMPAL input features (flair, qt1, thickness, adc); if missing a feature, skip this module
featList_hipp=()
if [ -f "${out//micapipe/}/${idBIDS}_space-hipp_hemi-rh_midthickness_desc-flair_2mm.func.gii" ]; then
    featList_hipp+=("flair"); fi
if [ -f "${out//micapipe/}/${idBIDS}_space-hipp_hemi-rh_midthickness_desc-qt1_2mm.func.gii" ]; then
    featList_hipp+=("qt1"); fi
if [ -f "${out//micapipe/}/${idBIDS}_space-hipp_hemi-rh_midthickness_desc-ADC_2mm.func.gii" ]; then
    featList_hipp+=("adc"); fi
if [ -f "${out//micapipe/}/${idBIDS}_space-hipp_hemi-rh_desc-thickness_2mm.func.gii" ]; then
    featList_hipp+=("thickness"); fi

# Check input feature
Info "Inputs:"
if [[ "$thr" == DEFAULT ]]; then Note "default z-score threshold at |1.96|"; else
Note "thr                   :" "$thr"; fi
if [[ "$demo" == '' ]]; then Note "No demographic file specified"; else
Note "demo                  :" "$demo"; fi

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe-z_asymmetry_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${out/micapipe/}/z-brains/regional/sub-${id}/${SES}/"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"


#------------------------------------------------------------------------------#
### regional analysis ###
Do_cmd python "$ZBRAINS"/functions/02_regional_analysis.py "sub-$id" "$SES" "$out" "$hipDir" \
              "${demo}" "$thr" \
              --featList_ctx "${featList_ctx[@]}" --featList_sctx "${featList_sctx[@]}" \
              --featList_hipp "${featList_hipp[@]}"

# Generate pdf report
Do_cmd python "$ZBRAINS"/functions/02_regional_analysis_report.py "sub-$id" "$SES" "$out" \
              "${demo}" "${thr}" \
              --featList_ctx "${featList_ctx[@]}" --featList_sctx "${featList_sctx[@]}" \
              --featList_hipp "${featList_hipp[@]}"
              
if [[ -f "${outDir}/sub-${id}_regional_analysis.pdf" ]]; then ((Nsteps++)); fi


################################################
# continue here
################################################

# Surface register back to fsspace
for hemi in lh rh; do
  [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
  HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

  Do_cmd mri_convert ${outDir}/sub-${id}_ctx-mz_${hemi}.mgh ${outDir}/sub-${id}_ctx-mz_${hemi}.func.gii

  Do_cmd wb_command -metric-resample \
                    ${outDir}/sub-${id}_ctx-mz_${hemi}.func.gii \
                    "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
                    "${dir_conte69}/${idBIDS}_${hemi}_sphereReg.surf.gii" \
                    ADAP_BARY_AREA \
                    "${tmp}/${idBIDS}_space-fsspace-desc-${hemi}_ctx-mz.func.gii" \
                    -area-surfs \
                    "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_midthickness.surf.gii" \
                    "${dir_freesurfer}/surf/${hemi}.midthickness.surf.gii"

  Do_cmd mri_convert "${tmp}/${idBIDS}_space-fsspace-desc-${hemi}_ctx-mz.func.gii" \
                     "${tmp}/${idBIDS}_space-fsspace-desc-${hemi}_ctx-mz.mgh"
done

# Map surface to volume
Do_cmd mri_surf2vol --o "${outDir}/sub-${id}_ctx-mz.nii.gz" \
                    --subject "$idBIDS" \
                    --so ${dir_freesurfer}/surf/lh.white "${tmp}/${idBIDS}_space-fsspace-desc-lh_ctx-mz.mgh" \
                    --so ${dir_freesurfer}/surf/rh.white "${tmp}/${idBIDS}_space-fsspace-desc-rh_ctx-mz.mgh"

if [[ -f "${outDir}/sub-${id}_ctx-mz.nii.gz" ]]; then ((Nsteps++)); fi

#------------------------------------------------------------------------------#
# QC notification of completion
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completion
if [ "$Nsteps" -eq 2 ]; then status="COMPLETED"; else status="ERROR regional analysis is missing a processing step"; fi
Title "asymmetry processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/02
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/regional_z_*.txt)"
echo "${id}, ${SES/ses-/}, regional_z, $status N=$(printf "%02d" "$Nsteps")/02, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
