#!/bin/bash
#
# hippocampal feature mapping:
#
# Generates hippocampal feature outputs:
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
fsdir=$5
nocleanup=$6
threads=$7
tmpDir=$8
featStr=$9
smoothing=${10}
PROC=${11}
export OMP_NUM_THREADS=$threads
here=$(pwd)

echo "PROC:"
echo $PROC

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs//micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Get the real path of the Inputs
out=${out/micapipe_v0.2.0/}
outz=$(realpath $out)/z-brains
outm=$(realpath $out)/micapipe_v0.2.0
BIDS=$(realpath $BIDS)
id=${id/sub-/}
here=$(pwd)

#------------------------------------------------------------------------------#
Title "hippocampal feature mapping\n\t\tmicapipe-z $this_version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Saving temporary dir: $nocleanup"

# Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe-z_hippo_proc_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Set output directory
if [ "$SES" == "SINGLE" ]; then
  export subject_dirz=$outz/${subject}
  ses=""
else
  export subject_dirz=$outz/${subject}/${SES}
  ses="_${SES}"
fi

# Data location
hippoDir="${out}/hippunfold_v1.3.0/hippunfold/sub-${id}/${SES}/"
subDeriv="${out}/micapipe_v0.2.0/sub-${id}/${SES}/"
mapsDir="${subDeriv}/maps/"
outLogs="${subject_dirz}/logs/"

#------------------------------------------------------------------------------#
# Optional argument handling

# Manage manual inputs: feature processing
featList_hippo=()
if [[ "$featStr" != "all" ]]; then
    IFS=',' read -ra featList_hippo <<< "$featStr"
elif [[ "$featStr" == "all" ]]; then
    if [ -f "${mapsDir}/${idBIDS}_space-nativepro_map-flair.nii.gz" ]; then
        featList_hippo+=("flair"); fi
    if [ -f "${mapsDir}/${idBIDS}_space-nativepro_map-T1map.nii.gz" ]; then
        featList_hippo+=("qt1"); fi
    if [ -f "${mapsDir}/${idBIDS}_space-nativepro_model-DTI_map-ADC.nii.gz" ]; then
        featList_hippo+=("ADC"); fi
    if [ -f "${mapsDir}/${idBIDS}_space-nativepro_model-DTI_map-FA.nii.gz" ]; then
        featList_hippo+=("FA"); fi
    if [ -f "${hippoDir}/surf/${idBIDS}_hemi-R_space-T1w_den-0p5mm_label-hipp_thickness.shape.gii" ]; then
        featList_hippo+=("thickness"); fi
fi

# smoothing
if [[ "$smoothing" == "DEFAULT" ]]; then fwhm=2; else fwhm=${smoothing}; fi

#------------------------------------------------------------------------------#
### Hippocampal thickness ###

Info "Fetching thickness from hippunfold"

hippThickness="${hippoDir}/surf/${idBIDS}_hemi-R_space-T1w_den-0p5mm_label-hipp_thickness.shape.gii"
thickness_out="${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-R_feature-thickness_smooth-${fwhm}mm.func.gii"

if [[ "${featList_hippo[*]}" =~ 'thickness' ]]; then
    N=$((N + 2))

    if [[ -f "${hippThickness}" ]]; then

        # Smooth thickness and output
        for hemi in lh rh; do
            [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
            HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

            Do_cmd wb_command -metric-smoothing "${hippoDir}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                              "${hippThickness}" \
                              ${fwhm} \
                              "${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-${HEMICAP}_feature-thickness_smooth-${fwhm}mm.func.gii"
        done

        if [[ -f "${thickness_out}" ]]; then Nsteps=$((Nsteps + 2)); fi

    else
        Note "Thickness processing requested but did not find hippunfold thickness file: skipping"
    fi
else
    Note "Skipping hippocampal thickness"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
### Hippocampal flair ###

Info "Map T2/FLAIR to hippocampus"

flair_preproc="${mapsDir}/${idBIDS}_space-nativepro_map-flair.nii.gz"
flair_out="${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-R_feature-flair_smooth-${fwhm}mm.func.gii"
if [[ "${featList_hippo[*]}" =~ 'flair' ]]; then
    N=$((N + 2))

    if [[ -f "${flair_preproc}" ]]; then

        # Map flair intensities to hippocampal subfields
        if [[ ! -f "${flair_out}" ]]; then

            for hemi in lh rh; do
                [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
                HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

                Do_cmd wb_command -volume-to-surface-mapping "$flair_preproc" \
                                  "${hippoDir}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                                  "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-T1w_desc-flair_den-0p5mm_feature-hipp_midthickness.func.gii" \
                                  -trilinear

                Do_cmd wb_command -metric-smoothing "${hippoDir}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                                  "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-T1w_desc-flair_den-0p5mm_feature-hipp_midthickness.func.gii" \
                                  ${fwhm} \
                                  "${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-${HEMICAP}_feature-flair_smooth-${fwhm}mm.func.gii"

            done

            if [[ -f "${flair_out}" ]]; then Nsteps=$((Nsteps + 2)); fi

        else
            Note "Subject ${idBIDS} T2-FLAIR is mapped to hippocampus"; Nsteps=$((Nsteps + 2))
        fi

    else
        Note "T2/FLAIR processing requested but did not find pre-processed file: skipping"
    fi
else
    Note "Skipping T2/FLAIR"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
### Hippocampal ADC and FA ###

Info "Map FA/ADC to hippocampus"

# ADC
adc_preproc="${mapsDir}/${idBIDS}_space-nativepro_model-DTI_map-ADC.nii.gz"
adc_out="${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-R_feature-ADC_smooth-${fwhm}mm.func.gii"
if [[ "${featList_hippo[*]}" =~ 'ADC' ]]; then
    N=$((N + 2))

    if [[ -f "${adc_preproc}" ]]; then

        # Map flair intensities to hippocampal subfields
        if [[ ! -f "${adc_out}" ]]; then

            for hemi in lh rh; do
                [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
                HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

                Do_cmd wb_command -volume-to-surface-mapping "$adc_preproc" \
                                  "${hippoDir}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                                  "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-T1w_desc-ADC_den-0p5mm_feature-hipp_midthickness.func.gii" \
                                  -trilinear

                Do_cmd wb_command -metric-smoothing "${hippoDir}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                                  "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-T1w_desc-ADC_den-0p5mm_feature-hipp_midthickness.func.gii" \
                                  ${fwhm} \
                                  "${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-${HEMICAP}_feature-ADC_smooth-${fwhm}mm.func.gii"

            done

            if [[ -f "${adc_out}" ]]; then Nsteps=$((Nsteps + 2)); fi

        else
            Note "Subject ${idBIDS} ADC is mapped to hippocampus"; Nsteps=$((Nsteps + 2))
        fi

    else
        Note "ADC processing requested but did not find pre-processed file: skipping"
    fi
else
    Note "Skipping ADC"; Nsteps=$((Nsteps + 2))
fi

# FA
fa_preproc="${mapsDir}/${idBIDS}_space-nativepro_model-DTI_map-FA.nii.gz"
fa_out="${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-R_feature-FA_smooth-${fwhm}mm.func.gii"
if [[ "${featList_hippo[*]}" =~ 'FA' ]]; then
    N=$((N + 2))

    if [[ -f "${fa_preproc}" ]]; then

        # Map flair intensities to hippocampal subfields
        if [[ ! -f "${fa_out}" ]]; then

            for hemi in lh rh; do
                [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
                HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

                Do_cmd wb_command -volume-to-surface-mapping "$fa_preproc" \
                                  "${hippoDir}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                                  "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-T1w_desc-FA_den-0p5mm_feature-hipp_midthickness.func.gii" \
                                  -trilinear

                Do_cmd wb_command -metric-smoothing "${hippoDir}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                                  "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-T1w_desc-FA_den-0p5mm_feature-hipp_midthickness.func.gii" \
                                  ${fwhm} \
                                  "${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-${HEMICAP}_feature-FA_smooth-${fwhm}mm.func.gii"

            done

            if [[ -f "${fa_out}" ]]; then Nsteps=$((Nsteps + 2)); fi

        else
            Note "Subject ${idBIDS} FA is mapped to hippocampus"; Nsteps=$((Nsteps + 2))
        fi

    else
        Note "FA processing requested but did not find pre-processed file: skipping"
    fi
else
    Note "Skipping FA"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
### Hippocampal qT1 ###

Info "Map qT1 to hippocampus"

qt1_preproc="${mapsDir}/${idBIDS}_space-nativepro_map-T1map.nii.gz"
qt1_out="${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-R_feature-T1map_smooth-${fwhm}mm.func.gii"
if [[ "${featList_hippo[*]}" =~ 'qt1' ]]; then 
    N=$((N + 2))

    if [[ -f "${qt1_preproc}" ]]; then

        # Map flair intensities to hippocampal subfields
        if [[ ! -f "${qt1_out}" ]]; then

            for hemi in lh rh; do
                [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
                HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

                Do_cmd wb_command -volume-to-surface-mapping "$qt1_preproc" \
                                  "${hippoDir}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                                  "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-T1w_desc-T1map_den-0p5mm_feature-hipp_midthickness.func.gii" \
                                  -trilinear

                Do_cmd wb_command -metric-smoothing "${hippoDir}/surf/${idBIDS}_hemi-${HEMICAP}_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii" \
                                  "${tmp}/${idBIDS}_hemi-${HEMICAP}_space-T1w_desc-T1map_den-0p5mm_feature-hipp_midthickness.func.gii" \
                                  ${fwhm} \
                                  "${subject_dirz}/maps/hippocampus/${idBIDS}_hemi-${HEMICAP}_feature-T1map_smooth-${fwhm}mm.func.gii"

            done

            if [[ -f "${qt1_out}" ]]; then Nsteps=$((Nsteps + 2)); fi

        else
            Note "Subject ${idBIDS} qT1 is mapped to hippocampus"; Nsteps=$((Nsteps + 2))
        fi

    else
        Note "qT1 processing requested but did not find pre-processed file: skipping"
    fi

else
    Note "Skipping qT1"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 10 ]; then status="COMPLETED"; else status="ERROR hippo_proc is missing a processing step"; fi
Title "hippo processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/10
\tStatus          : ${status}
\tCheck logs      : $(ls "${outLogs}"/hippo_proc_*.txt)"
#echo "${id}, ${SES/ses-/}, sctx_proc, $status N=$(printf "%02d" "$Nsteps")/08, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${outz}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
