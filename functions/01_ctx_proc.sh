#!/bin/bash
#
# Cortical feature mapping:
#
# Generates cortical volume outputs:
#
# This workflow makes use of freesurfer outputs and custom python scripts
#
umask 003
id=$1
indir=$2
outm=$3
outz=$4
SES=$5
fsdir=$6
nocleanup=$7
threads=$8
tmpDir=$9
featStr=${10}
fwhm=${11}
outRes=${12}
PROC=${13}
export OMP_NUM_THREADS=$threads
here=$(pwd)


#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export ZBRAINS=${$ZBRAINS}
    source "${ZBRAINS}/functions/init.sh" "$threads"
fi

# source utilities
source "$ZBRAINS/functions/utilities.sh"

# Get the real path of the Inputs
outz=$(realpath $outz)
outm=$(realpath $outm)
indir=$(realpath $indir)
id=${id/sub-/}
here=$(pwd)


#------------------------------------------------------------------------------#
Title "Cortical feature mapping\n\t\tmicapipe-z $this_version, $PROC"
zbrains_software
bids_print.variables-ctx
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Saving temporary directory: $nocleanup"

# Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_zbrains_ctx_proc_${idBIDS}"
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
mapsDir="${subject_micapipe}/maps/"


#------------------------------------------------------------------------------#
# Optional argument handling

# Manage manual inputs: 
# Feature processing
featList_ctx=()
if [[ "$featStr" != "all" ]]; then
    IFS=',' read -ra featList_ctx <<< "$featStr"
elif [[ "$featStr" == "all" ]]; then
    if [ -f "${mapsDir}/${idBIDS}_hemi-L_surf-fsLR-32k_label-white_flair.func.gii" ]; then
        featList_ctx+=("flair"); fi
    if [ -f "${mapsDir}/${idBIDS}_hemi-L_surf-fsLR-32k_label-white_T1map.func.gii" ]; then
        featList_ctx+=("qt1"); fi
    if [ -f "${mapsDir}/${idBIDS}_hemi-L_surf-fsLR-32k_label-white_ADC.func.gii" ]; then
        featList_ctx+=("ADC"); fi
    if [ -f "${mapsDir}/${idBIDS}_hemi-L_surf-fsLR-32k_label-white_FA.func.gii" ]; then
        featList_ctx+=("FA"); fi
    if [ -f "${mapsDir}/${idBIDS}_hemi-R_surf-fsLR-32k_label-thickness.func.gii" ]; then
        featList_ctx+=("thickness"); fi
fi


#------------------------------------------------------------------------------#
### Cortical thickness ###

Info "Fetching thickness from cortex"

cortThickness="${subject_micapipe}/maps/${idBIDS}_hemi-R_surf-fsLR-32k_label-thickness.func.gii"
thickness_outHigh="${subject_dirz}/maps/cortex/${idBIDS}_hemi-R_surf-fsLR-32k_feature-thickness_smooth-${fwhm}mm.func.gii"
thickness_outLow="${subject_dirz}/maps/cortex/${idBIDS}_hemi-R_surf-fsLR-5k_feature-thickness_smooth-${fwhm}mm.func.gii" 

if [[ "${featList_ctx[*]}" =~ 'thickness' ]]; then 
    N=$((N + 2))

    if [[ -f "${cortThickness}" ]]; then 
        
        if [ ! -f "${thickness_outHigh}" ] || [ ! -f "${thickness_outLow}" ]; then
        
            # Smooth thickness and output
            for hemi in lh rh; do 
                [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
                HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])
                            
                if [ "$outRes" == "all" ] || [ "$outRes" == "high" ]; then
                
                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-32k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_label-thickness.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_feature-thickness_smooth-${fwhm}mm.func.gii"
                fi
                
                if [ "$outRes" == "all" ] || [ "$outRes" == "low" ]; then
        
                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-5k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_label-thickness.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_feature-thickness_smooth-${fwhm}mm.func.gii"
                fi
                                      
            done
            
            if [ -f "${thickness_outHigh}" ] || [ -f "${thickness_outLow}" ]; then Nsteps=$((Nsteps + 2)); fi
        
        else
            Note "Subject ${idBIDS} thickness is mapped to cortex"; Nsteps=$((Nsteps + 2))
        fi
        
    else
        Note "Thickness processing requested but did not find cortical thickness file: skipping"
    fi
else
    Note "Skipping cortical thickness"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
### Cortical flair ###

Info "Map T2/FLAIR to cortex"

flair_preproc="${mapsDir}/${idBIDS}_hemi-L_surf-fsLR-32k_label-white_flair.func.gii" 
flair_outHigh="${subject_dirz}/maps/cortex/${idBIDS}_hemi-L_surf-fsLR-32k_feature-flair_smooth-${fwhm}mm.func.gii"
flair_outLow="${subject_dirz}/maps/cortex/${idBIDS}_hemi-L_surf-fsLR-5k_feature-flair_smooth-${fwhm}mm.func.gii" 

if [[ "${featList_ctx[*]}" =~ 'flair' ]]; then
    N=$((N + 2))

    if [[ -f "${flair_preproc}" ]]; then

        # Map flair intensities to cortical subfields
        if [ ! -f "${flair_outHigh}" ] || [ ! -f "${flair_outLow}" ]; then

            for hemi in lh rh; do
                [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r 
                HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:]) 
                
                if [ "$outRes" == "all" ] || [ "$outRes" == "high" ]; then
                
                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-32k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_label-white_flair.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_feature-flair_smooth-${fwhm}mm.func.gii"
                fi
                
                if [ "$outRes" == "all" ] || [ "$outRes" == "low" ]; then
                
                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-5k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_label-white_flair.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_feature-flair_smooth-${fwhm}mm.func.gii"
                fi
                
            done

            if [ -f "${flair_outHigh}" ] || [ -f "${flair_outLow}" ]; then Nsteps=$((Nsteps + 2)); fi

        else
            Note "Subject ${idBIDS} T2-FLAIR is mapped to cortex"; Nsteps=$((Nsteps + 2))
        fi

    else
        Note "T2/FLAIR processing requested but did not find pre-processed file: skipping"
    fi
else
    Note "Skipping T2/FLAIR"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
### Cortical ADC and FA ###

Info "Map FA/ADC to cortex"

# ADC
adc_preproc="${mapsDir}/${idBIDS}_hemi-L_surf-fsLR-32k_label-white_ADC.func.gii" 
adc_outHigh="${subject_dirz}/maps/cortex/${idBIDS}_hemi-L_surf-fsLR-32k_feature-ADC_smooth-${fwhm}mm.func.gii" 
adc_outLow="${subject_dirz}/maps/cortex/${idBIDS}_hemi-L_surf-fsLR-5k_feature-ADC_smooth-${fwhm}mm.func.gii" 

if [[ "${featList_ctx[*]}" =~ 'ADC' ]]; then
    N=$((N + 2))

    if [[ -f "${adc_preproc}" ]]; then

        # Map ADC intensities to cortical subfields
        if [ ! -f "${adc_outHigh}" ] || [ ! -f "${adc_outLow}" ]; then

            for hemi in lh rh; do
                [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
                HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])
                
                if [ "$outRes" == "all" ] || [ "$outRes" == "high" ]; then
                
                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-32k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_label-white_ADC.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_feature-ADC_smooth-${fwhm}mm.func.gii"
                fi
                
                if [ "$outRes" == "all" ] || [ "$outRes" == "low" ]; then
                    
                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-5k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_label-white_ADC.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_feature-ADC_smooth-${fwhm}mm.func.gii"
                fi
                
            done

            if [ -f "${adc_outHigh}" ] || [ -f "${adc_outLow}" ]; then Nsteps=$((Nsteps + 2)); fi

        else
            Note "Subject ${idBIDS} ADC is mapped to cortex"; Nsteps=$((Nsteps + 2))
        fi

    else
        Note "ADC processing requested but did not find pre-processed file: skipping"
    fi
else
    Note "Skipping ADC"; Nsteps=$((Nsteps + 2))
fi

# FA
fa_preproc="${mapsDir}/${idBIDS}_hemi-L_surf-fsLR-32k_label-white_FA.func.gii"
fa_outHigh="${subject_dirz}/maps/cortex/${idBIDS}_hemi-L_surf-fsLR-32k_feature-FA_smooth-${fwhm}mm.func.gii"
fa_outLow="${subject_dirz}/maps/cortex/${idBIDS}_hemi-L_surf-fsLR-5k_feature-FA_smooth-${fwhm}mm.func.gii"

if [[ "${featList_ctx[*]}" =~ 'FA' ]]; then
    N=$((N + 2))

    if [[ -f "${fa_preproc}" ]]; then

        # Map FA intensities to cortical subfields
        if [ ! -f "${fa_outHigh}" ] || [ ! -f "${fa_outLow}" ]; then

            for hemi in lh rh; do
                [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
                HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])
                
                if [ "$outRes" == "all" ] || [ "$outRes" == "high" ]; then
                
                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-32k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_label-white_FA.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_feature-FA_smooth-${fwhm}mm.func.gii"
                fi
                
                if [ "$outRes" == "all" ] || [ "$outRes" == "low" ]; then

                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-5k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_label-white_FA.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_feature-FA_smooth-${fwhm}mm.func.gii"
                fi
                
            done

            if [ -f "${fa_outHigh}" ] || [ -f "${fa_outLow}" ]; then Nsteps=$((Nsteps + 2)); fi

        else
            Note "Subject ${idBIDS} FA is mapped to cortex"; Nsteps=$((Nsteps + 2))
        fi

    else
        Note "FA processing requested but did not find pre-processed file: skipping"
    fi
else
    Note "Skipping FA"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
### Cortical qT1 ###

Info "Map qT1 to cortex"

qt1_preproc="${mapsDir}/${idBIDS}_hemi-L_surf-fsLR-32k_label-white_T1map.func.gii"
qt1_outHigh="${subject_dirz}/maps/cortex/${idBIDS}_hemi-L_surf-fsLR-32k_feature-T1map_smooth-${fwhm}mm.func.gii"
qt1_outLow="${subject_dirz}/maps/cortex/${idBIDS}_hemi-L_surf-fsLR-5k_feature-T1map_smooth-${fwhm}mm.func.gii"

if [[ "${featList_ctx[*]}" =~ 'qt1' ]]; then
    N=$((N + 2))

    if [[ -f "${qt1_preproc}" ]]; then

        # Map qT1 intensities to cortical subfields
        if [ ! -f "${qt1_outHigh}" ] || [ ! -f "${qt1_outLow}" ]; then

            for hemi in lh rh; do
                [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
                HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

                if [ "$outRes" == "all" ] || [ "$outRes" == "high" ]; then

                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-32k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_label-white_T1map.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-32k_feature-T1map_smooth-${fwhm}mm.func.gii"
                fi
                
                if [ "$outRes" == "all" ] || [ "$outRes" == "low" ]; then

                    Do_cmd wb_command -metric-smoothing "${subject_micapipe}/surf/${idBIDS}_hemi-${HEMICAP}_space-nativepro_surf-fsLR-5k_label-white.surf.gii" \
                                      "${subject_micapipe}/maps/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_label-white_T1map.func.gii" \
                                      ${fwhm} \
                                      "${subject_dirz}/maps/cortex/${idBIDS}_hemi-${HEMICAP}_surf-fsLR-5k_feature-T1map_smooth-${fwhm}mm.func.gii"
                fi
                
            done

            if [ -f "${qt1_outHigh}" ] || [ -f "${qt1_outLow}" ]; then Nsteps=$((Nsteps + 2)); fi

        else
            Note "Subject ${idBIDS} qT1 is mapped to cortex"; Nsteps=$((Nsteps + 2))
        fi

    else
        Note "qT1 processing requested but did not find pre-processed file: skipping"
    fi

else
    Note "Skipping qT1"; Nsteps=$((Nsteps + 2))
fi

echo $N
echo $Nsteps

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 10 ]; then status="COMPLETED"; else status="ERROR ctx_proc is missing a processing step"; fi
Title "ctx processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/10
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/ctx_proc_*.txt)"
#echo "${id}, ${SES/ses-/}, ctx_proc, $status N=$(printf "%02d" "$Nsteps")/08, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${outz}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
