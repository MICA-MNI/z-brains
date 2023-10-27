#!/bin/bash
#
# Subcortical feature mapping:
#
# Generates subcortical feature outputs:
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
PROC=${10}
export OMP_NUM_THREADS=$threads
here=$(pwd)

echo "PROC:"
echo $PROC

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Get the real path of the Inputs
outz=$(realpath $out)/z-brains
outm=$(realpath $out)/micapipe_v0.2.0
BIDS=$(realpath $BIDS)
id=${id/sub-/}
here=$(pwd)

#------------------------------------------------------------------------------#
Title "subcortical feature mapping\n\t\tmicapipe-z $this_version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Saving temporary dir: $nocleanup"

# Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe-z_sctx_proc_${idBIDS}"
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
volDir="${fsdir}/stats/"
subDeriv="${out}/micapipe_v0.2.0/sub-${id}/${SES}/"
mapsDir="${subDeriv}/maps/"
outLogs="${subject_dirz}/logs/"


#------------------------------------------------------------------------------#
# Optional argument handling

# Manage manual inputs: feature processing
featList_sctx=()
if [[ "$featStr" != "DEFAULT" ]]; then
    IFS=',' read -ra featList_sctx <<< "$featStr"
elif [[ "$featStr" == "DEFAULT" ]]; then
    if [ -f "${mapsDir}/${idBIDS}_space-nativepro_map-flair.nii.gz" ]; then
        featList_sctx+=("flair"); fi
    if [ -f "${mapsDir}/${idBIDS}_space-nativepro_map-T1map.nii.gz" ]; then
        featList_sctx+=("qt1"); fi
    if [ -f "${mapsDir}/${idBIDS}_space-nativepro_model-DTI_map-ADC.nii.gz" ]; then
        featList_sctx+=("ADC"); fi
    if [ -f "${mapsDir}/${idBIDS}_space-nativepro_model-DTI_map-FA.nii.gz" ]; then
        featList_sctx+=("FA"); fi
    if [ -f "${volDir}/aseg.stats" ]; then
        featList_sctx+=("volume"); fi
fi

#------------------------------------------------------------------------------#
### subcortical volume ###

Info "Fetching subcortical volumes from FreeSurfer"

if [[ "${featList_sctx[*]}" =~ 'volume' ]]; then 
    ((N++))

    if [[ -f "${volDir}/aseg.stats" ]]; then
        
        # Create file header
        echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal,ICV" > \
            "${subject_dirz}/maps/subcortex/${idBIDS}_feature-volume.csv" 
        printf "%s,"  "${idBIDS}" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-volume.csv"
                
        for sub in Left-Accumbens-area Left-Amygdala Left-Caudate Left-Hippocampus Left-Pallidum \
        Left-Putamen Left-Thalamus-Proper Right-Accumbens-area Right-Amygdala \
        Right-Caudate Right-Hippocampus Right-Pallidum Right-Putamen Right-Thalamus-Proper; do
            printf "%g," `grep  ${sub} ${volDir}/aseg.stats | awk '{print $4}'` >> \
            "${subject_dirz}/maps/subcortex/${idBIDS}_feature-volume.csv" 
        done
        
        printf "%g" `cat ${volDir}/aseg.stats | grep IntraCranialVol | awk -F, '{print $4}'` >> \
             "${subject_dirz}/maps/subcortex/${idBIDS}_feature-volume.csv" 
        echo "" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-volume.csv" 
    
        if [[ -f "${subject_dirz}/maps/subcortex/${idBIDS}_feature-volume.csv" ]]; then ((Nsteps++)); fi
        
    else
        Note "Volumetric processing requested but did not find subcortical volume file: skipping"
    fi
else
    Note "Skipping subcortical volumetric"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
### subcortical flair ###

Info "Map T2/FLAIR to subcortical structures"

flair_preproc="${mapsDir}/${idBIDS}_space-nativepro_map-flair.nii.gz"
if [[ "${featList_sctx[*]}" =~ 'flair' ]]; then 
    ((N++))
    
    if [[ -f "${flair_preproc}" ]]; then
    
        # Create file header
        echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal" > \
                "${subject_dirz}/maps/subcortex/${idBIDS}_feature-flair.csv" 
        printf "%s,"  "${idBIDS}" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-flair.csv" 
    
        for sub in 26 18 11 17 13 12 10 58 54 50 53 52 51 49; do
            if [[ ${sub} == 26 ]]; then sctxname="Left-Accumbens-area"; elif [[ ${sub} == 18 ]]; then sctxname="Left-Amygdala"; \
            elif [[ ${sub} == 11 ]]; then sctxname="Left-Caudate"; elif [[ ${sub} == 17 ]]; then sctxname="Left-Hippocampus"; \
            elif [[ ${sub} == 13 ]]; then sctxname="Left-Pallidum"; elif [[ ${sub} == 12 ]]; then sctxname="Left-Putamen"; \
            elif [[ ${sub} == 10 ]]; then sctxname="Left-Thalamus-Proper"; elif [[ ${sub} == 58 ]]; then sctxname="Right-Accumbens-area"; \
            elif [[ ${sub} == 54 ]]; then sctxname="Right-Amygdala"; elif [[ ${sub} == 50 ]]; then sctxname="Right-Caudate"; \
            elif [[ ${sub} == 53 ]]; then sctxname="Right-Hippocampus"; elif [[ ${sub} == 52 ]]; then sctxname="Right-Pallidum"; \
            elif [[ ${sub} == 51 ]]; then sctxname="Right-Putamen"; elif [[ ${sub} == 49 ]]; then sctxname="Right-Thalamus-Proper"; fi
    
            # Extract subcortical masks
            Do_cmd ThresholdImage 3 "${subDeriv}/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz" \
                                "${tmp}/${idBIDS}_${sctxname}_mask1.nii.gz" \
                                threshlo "${sub}"
                                
            Do_cmd ThresholdImage 3 "${subDeriv}/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz" \
                                "${tmp}/${idBIDS}_${sctxname}_mask2.nii.gz" \
                                threshhi `awk "BEGIN {print $sub-1}"`
                                
            Do_cmd ImageMath 3 "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" \
                                - "${tmp}/${idBIDS}_${sctxname}_mask1.nii.gz" "${tmp}/${idBIDS}_${sctxname}_mask2.nii.gz"
            
            # Get flair intensities for subcortical mask
            Do_cmd ImageMath 3 "${tmp}/${idBIDS}_${sctxname}_masked-flair.nii.gz" \
                            m "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" "$flair_preproc"
            
            # Get mean intensity in sampled region
            echo `ImageIntensityStatistics 3 "${tmp}/${idBIDS}_${sctxname}_masked-flair.nii.gz"` >> "${tmp}/stats_flair_${sctxname}.txt"
            stats_reg=($(cat "${tmp}/stats_flair_${sctxname}.txt"))
            sum_reg=`echo "${stats_reg[17]}"`
            
            echo `ImageIntensityStatistics 3 "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz"` >> "${tmp}/stats_mask_flair_${sctxname}.txt"
            stats_mask=($(cat "${tmp}/stats_mask_flair_${sctxname}.txt"))
            sum_mask=`echo "${stats_mask[17]}"`
            
            reg_mean=`awk "BEGIN {print $sum_reg/$sum_mask}"`
            
            # Input values in .csv file
            if [[ "${sctxname}" == "Right-Thalamus-Proper" ]]; then
                printf "%g" `echo $reg_mean` >> \
                    "${subject_dirz}/maps/subcortex/${idBIDS}_feature-flair.csv" 
            else
                printf "%g," `echo $reg_mean` >> \
                    "${subject_dirz}/maps/subcortex/${idBIDS}_feature-flair.csv"
            fi 
        done
        
        echo "" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-flair.csv"
        if [[ -f "${subject_dirz}/maps/subcortex/${idBIDS}_feature-flair.csv" ]]; then ((Nsteps++)); fi

    else
        Note "T2/FLAIR processing requested but did not find pre-processed file: skipping"
    fi    
else
    Note "Skipping T2/FLAIR"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
### subcortical ADC and FA ###

Info "Map FA/ADC to subcortical structures"

# ADC
adc_preproc="${mapsDir}/${idBIDS}_space-nativepro_model-DTI_map-ADC.nii.gz"
if [[ "${featList_sctx[*]}" =~ 'ADC' ]]; then 
    ((N++))

    if [[ -f "${adc_preproc}" ]]; then

        echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal" > \
                "${subject_dirz}/maps/subcortex/${idBIDS}_feature-ADC.csv" 
        printf "%s,"  "${idBIDS}" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-ADC.csv" 
    
        for sub in 26 18 11 17 13 12 10 58 54 50 53 52 51 49; do
            if [[ ${sub} == 26 ]]; then sctxname="Left-Accumbens-area"; elif [[ ${sub} == 18 ]]; then sctxname="Left-Amygdala"; \
            elif [[ ${sub} == 11 ]]; then sctxname="Left-Caudate"; elif [[ ${sub} == 17 ]]; then sctxname="Left-Hippocampus"; \
            elif [[ ${sub} == 13 ]]; then sctxname="Left-Pallidum"; elif [[ ${sub} == 12 ]]; then sctxname="Left-Putamen"; \
            elif [[ ${sub} == 10 ]]; then sctxname="Left-Thalamus-Proper"; elif [[ ${sub} == 58 ]]; then sctxname="Right-Accumbens-area"; \
            elif [[ ${sub} == 54 ]]; then sctxname="Right-Amygdala"; elif [[ ${sub} == 50 ]]; then sctxname="Right-Caudate"; \
            elif [[ ${sub} == 53 ]]; then sctxname="Right-Hippocampus"; elif [[ ${sub} == 52 ]]; then sctxname="Right-Pallidum"; \
            elif [[ ${sub} == 51 ]]; then sctxname="Right-Putamen"; elif [[ ${sub} == 49 ]]; then sctxname="Right-Thalamus-Proper"; fi
    
            # Extract subcortical masks
            Do_cmd ThresholdImage 3 "${subDeriv}/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz" \
                                "${tmp}/${idBIDS}_${sctxname}_mask1.nii.gz" \
                                threshlo "${sub}"
                                
            Do_cmd ThresholdImage 3 "${subDeriv}/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz" \
                                "${tmp}/${idBIDS}_${sctxname}_mask2.nii.gz" \
                                threshhi `awk "BEGIN {print $sub-1}"`
                                
            Do_cmd ImageMath 3 "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" \
                                - "${tmp}/${idBIDS}_${sctxname}_mask1.nii.gz" "${tmp}/${idBIDS}_${sctxname}_mask2.nii.gz"
            
            # Get flair intensities for subcortical mask
            Do_cmd ImageMath 3 "${tmp}/${idBIDS}_${sctxname}_masked-ADC.nii.gz" \
                            m "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" "$adc_preproc"
            
            # Get mean intensity in sampled region
            echo `ImageIntensityStatistics 3 "${tmp}/${idBIDS}_${sctxname}_masked-ADC.nii.gz"` >> "${tmp}/stats_ADC_${sctxname}.txt"
            stats_reg=($(cat "${tmp}/stats_ADC_${sctxname}.txt"))
            sum_reg=`echo "${stats_reg[17]}"`
            
            echo `ImageIntensityStatistics 3 "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz"` >> "${tmp}/stats_mask_ADC_${sctxname}.txt"
            stats_mask=($(cat "${tmp}/stats_mask_ADC_${sctxname}.txt"))
            sum_mask=`echo "${stats_mask[17]}"`
            
            reg_mean=`awk "BEGIN {print $sum_reg/$sum_mask}"`
            
            # Input values in .csv file
            if [[ "${sctxname}" == "Right-Thalamus-Proper" ]]; then
                printf "%g" `echo $reg_mean` >> \
                    "${subject_dirz}/maps/subcortex/${idBIDS}_feature-ADC.csv" 
            else
                printf "%g," `echo $reg_mean` >> \
                    "${subject_dirz}/maps/subcortex/${idBIDS}_feature-ADC.csv" 
            fi 
        done

        echo "" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-ADC.csv"
        if [[ -f "${subject_dirz}/maps/subcortex/${idBIDS}_feature-ADC.csv" ]]; then ((Nsteps++)); fi
        
    else
        Note "ADC processing requested but did not find pre-processed file: skipping"
    fi    
else
    Note "Skipping ADC"; ((Nsteps++))
fi

# FA
fa_preproc="${mapsDir}/${idBIDS}_space-nativepro_model-DTI_map-FA.nii.gz"
if [[ "${featList_sctx[*]}" =~ 'FA' ]]; then 
    ((N++))
    
    if [[ -f "${fa_preproc}" ]]; then

        echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal" > \
                "${subject_dirz}/maps/subcortex/${idBIDS}_feature-FA.csv" 
        printf "%s,"  "${idBIDS}" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-FA.csv"  
    
        for sub in 26 18 11 17 13 12 10 58 54 50 53 52 51 49; do
            if [[ ${sub} == 26 ]]; then sctxname="Left-Accumbens-area"; elif [[ ${sub} == 18 ]]; then sctxname="Left-Amygdala"; \
            elif [[ ${sub} == 11 ]]; then sctxname="Left-Caudate"; elif [[ ${sub} == 17 ]]; then sctxname="Left-Hippocampus"; \
            elif [[ ${sub} == 13 ]]; then sctxname="Left-Pallidum"; elif [[ ${sub} == 12 ]]; then sctxname="Left-Putamen"; \
            elif [[ ${sub} == 10 ]]; then sctxname="Left-Thalamus-Proper"; elif [[ ${sub} == 58 ]]; then sctxname="Right-Accumbens-area"; \
            elif [[ ${sub} == 54 ]]; then sctxname="Right-Amygdala"; elif [[ ${sub} == 50 ]]; then sctxname="Right-Caudate"; \
            elif [[ ${sub} == 53 ]]; then sctxname="Right-Hippocampus"; elif [[ ${sub} == 52 ]]; then sctxname="Right-Pallidum"; \
            elif [[ ${sub} == 51 ]]; then sctxname="Right-Putamen"; elif [[ ${sub} == 49 ]]; then sctxname="Right-Thalamus-Proper"; fi
    
            # Extract subcortical masks
            Do_cmd ThresholdImage 3 "${subDeriv}/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz" \
                                "${tmp}/${idBIDS}_${sctxname}_mask1.nii.gz" \
                                threshlo "${sub}"
                                
            Do_cmd ThresholdImage 3 "${subDeriv}/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz" \
                                "${tmp}/${idBIDS}_${sctxname}_mask2.nii.gz" \
                                threshhi `awk "BEGIN {print $sub-1}"`
                                
            Do_cmd ImageMath 3 "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" \
                                - "${tmp}/${idBIDS}_${sctxname}_mask1.nii.gz" "${tmp}/${idBIDS}_${sctxname}_mask2.nii.gz"
            
            # Get flair intensities for subcortical mask
            Do_cmd ImageMath 3 "${tmp}/${idBIDS}_${sctxname}_masked-FA.nii.gz" \
                            m "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" "$fa_preproc"
            
            # Get mean intensity in sampled region
            echo `ImageIntensityStatistics 3 "${tmp}/${idBIDS}_${sctxname}_masked-FA.nii.gz"` >> "${tmp}/stats_FA_${sctxname}.txt"
            stats_reg=($(cat "${tmp}/stats_FA_${sctxname}.txt"))
            sum_reg=`echo "${stats_reg[17]}"`
            
            echo `ImageIntensityStatistics 3 "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz"` >> "${tmp}/stats_mask_FA_${sctxname}.txt"
            stats_mask=($(cat "${tmp}/stats_mask_FA_${sctxname}.txt"))
            sum_mask=`echo "${stats_mask[17]}"`
            
            reg_mean=`awk "BEGIN {print $sum_reg/$sum_mask}"`
            
            # Input values in .csv file
            if [[ "${sctxname}" == "Right-Thalamus-Proper" ]]; then
                printf "%g" `echo $reg_mean` >> \
                    "${subject_dirz}/maps/subcortex/${idBIDS}_feature-FA.csv" 
            else
                printf "%g," `echo $reg_mean` >> \
                    "${subject_dirz}/maps/subcortex/${idBIDS}_feature-FA.csv"
            fi 
        done

        echo "" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-FA.csv"
        if [[ -f "${subject_dirz}/maps/subcortex/${idBIDS}_feature-FA.csv" ]]; then ((Nsteps++)); fi
        
    else
        Note "FA processing requested but did not find pre-processed file: skipping"
    fi
else
    Note "Skipping FA"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
### subcortical qT1 ###

Info "Map qT1 to subcortical structures"

qt1_preproc="${mapsDir}/${idBIDS}_space-nativepro_map-T1map.nii.gz"
if [[ "${featList_sctx[*]}" =~ 'qt1' ]]; then 
    ((N++))
    
    if [[ -f "${qt1_preproc}" ]]; then
    
        # Create file header
        echo "SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal" > \
                "${subject_dirz}/maps/subcortex/${idBIDS}_feature-T1map.csv" 
        printf "%s,"  "${idBIDS}" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-T1map.csv" 
    
        for sub in 26 18 11 17 13 12 10 58 54 50 53 52 51 49; do
            if [[ ${sub} == 26 ]]; then sctxname="Left-Accumbens-area"; elif [[ ${sub} == 18 ]]; then sctxname="Left-Amygdala"; \
            elif [[ ${sub} == 11 ]]; then sctxname="Left-Caudate"; elif [[ ${sub} == 17 ]]; then sctxname="Left-Hippocampus"; \
            elif [[ ${sub} == 13 ]]; then sctxname="Left-Pallidum"; elif [[ ${sub} == 12 ]]; then sctxname="Left-Putamen"; \
            elif [[ ${sub} == 10 ]]; then sctxname="Left-Thalamus-Proper"; elif [[ ${sub} == 58 ]]; then sctxname="Right-Accumbens-area"; \
            elif [[ ${sub} == 54 ]]; then sctxname="Right-Amygdala"; elif [[ ${sub} == 50 ]]; then sctxname="Right-Caudate"; \
            elif [[ ${sub} == 53 ]]; then sctxname="Right-Hippocampus"; elif [[ ${sub} == 52 ]]; then sctxname="Right-Pallidum"; \
            elif [[ ${sub} == 51 ]]; then sctxname="Right-Putamen"; elif [[ ${sub} == 49 ]]; then sctxname="Right-Thalamus-Proper"; fi
    
            # Extract subcortical masks
            Do_cmd ThresholdImage 3 "${subDeriv}/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz" \
                                "${tmp}/${idBIDS}_${sctxname}_mask1.nii.gz" \
                                threshlo "${sub}"
                                
            Do_cmd ThresholdImage 3 "${subDeriv}/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz" \
                                "${tmp}/${idBIDS}_${sctxname}_mask2.nii.gz" \
                                threshhi `awk "BEGIN {print $sub-1}"`
                                
            Do_cmd ImageMath 3 "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" \
                                - "${tmp}/${idBIDS}_${sctxname}_mask1.nii.gz" "${tmp}/${idBIDS}_${sctxname}_mask2.nii.gz"
            
            # Get flair intensities for subcortical mask
            Do_cmd ImageMath 3 "${tmp}/${idBIDS}_${sctxname}_masked-T1map.nii.gz" \
                            m "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz" "$qt1_preproc"
            
            # Get mean intensity in sampled region
            echo `ImageIntensityStatistics 3 "${tmp}/${idBIDS}_${sctxname}_masked-T1map.nii.gz"` >> "${tmp}/stats_T1map_${sctxname}.txt"
            stats_reg=($(cat "${tmp}/stats_T1map_${sctxname}.txt"))
            sum_reg=`echo "${stats_reg[17]}"`
            
            echo `ImageIntensityStatistics 3 "${tmp}/${idBIDS}_${sctxname}_mask.nii.gz"` >> "${tmp}/stats_mask_T1map_${sctxname}.txt"
            stats_mask=($(cat "${tmp}/stats_mask_T1map_${sctxname}.txt"))
            sum_mask=`echo "${stats_mask[17]}"`
            
            reg_mean=`awk "BEGIN {print $sum_reg/$sum_mask}"`
            
            # Input values in .csv file
            if [[ "${sctxname}" == "Right-Thalamus-Proper" ]]; then
                printf "%g" `echo $reg_mean` >> \
                    "${subject_dirz}/maps/subcortex/${idBIDS}_feature-T1map.csv"
            else
                printf "%g," `echo $reg_mean` >> \
                    "${subject_dirz}/maps/subcortex/${idBIDS}_feature-T1map.csv" 
            fi 
        done
        
        echo "" >> "${subject_dirz}/maps/subcortex/${idBIDS}_feature-T1map.csv"
        if [[ -f "${subject_dirz}/maps/subcortex/${idBIDS}_feature-T1map.csv" ]]; then ((Nsteps++)); fi
        
    else
        Note "qT1 processing requested but did not find pre-processed file: skipping"
    fi
else
    Note "Skipping qT1"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 5 ]; then status="COMPLETED"; else status="ERROR sctx_proc is missing a processing step"; fi
Title "sctx processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m.
\tSteps completed : $(printf "%02d" "$Nsteps")/05
\tStatus          : ${status}
\tCheck logs      : $(ls "${outLogs}"/sctx_proc_*.txt)"
#echo "${id}, ${SES/ses-/}, sctx_proc, $status N=$(printf "%02d" "$Nsteps")/05, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${outz}/micapipez_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
