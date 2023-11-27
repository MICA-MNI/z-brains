#!/bin/bash
#
# Feature mapping for cortex, subcortex, and hippocampus


if [[ -z ${ZBRAINS} ]]; then
  echo "ZBRAINS not defined"
  echo "This script should not be run standalone. Please use 'z-brains' directly."
  exit 0;
fi
source "${ZBRAINS}/config.cfg" # Configuration file


# Set umask
umask 003


#------------------------------------------------------------------------------#
while [ "$#" -gt 0 ]; do
  case "$1" in
    --structure)
      if [[ " ${LIST_STRUCTURES[*]} " != *" $2 "* ]]; then SHOW_ERROR "Invalid argument for $1: $2"; exit 1; fi
      structure="$2"
      shift 2
      ;;
    --features)
      comma_separated_features="$2"
      shift 2
      ;;
    --fwhm)
      fwhm="$2"
      shift 2
      ;;
    --tmp)
      tmp_dir="$2"
      shift 2
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    --resolution)
      comma_separated_resolutions="$2"
      shift 2
      ;;
    --logfile)
      logfile="$2"
      shift 2
      ;;
    *)
      SHOW_ERROR "Error: Unknown option '$1'"
      exit 1
      ;;
  esac
done

# Check if mandatory options are provided
declare -A mandatory=([struct]=$structure [features]=$comma_separated_features [tmp]=$tmp_dir [threads]=$threads)
for k in "${!mandatory[@]}"; do
  if [ -z "${mandatory[${k}]}" ]; then
    SHOW_ERROR "Mandatory argument is missing: -$k"
    exit 1
  fi
done

if [[ "$structure" != "subcortex" ]]; then
  [[ -z "$fwhm" ]] && SHOW_ERROR "Mandatory argument is missing: -fwhm" && exit 1;
  [[ -z "$comma_separated_resolutions" ]] && SHOW_ERROR "Mandatory argument is missing: -resolution" && exit 1;
fi


#------------------------------------------------------------------------------#
# Initialization of pertinent scripts and functions
export OMP_NUM_THREADS=$threads


#------------------------------------------------------------------------------#
declare -A map_struct=([cortex]=Cortical [subcortex]=Subcortical [hippocampus]=Hippocampal)

SHOW_TITLE "${map_struct[${structure}]} feature mapping"


#------------------------------------------------------------------------------#
# Initialize timer
SECONDS=0


#------------------------------------------------------------------------------#
declare -A label2struct
label2struct[26]="Left-Accumbens-area"
label2struct[18]="Left-Amygdala"
label2struct[11]="Left-Caudate"
label2struct[17]="Left-Hippocampus"
label2struct[13]="Left-Pallidum"
label2struct[12]="Left-Putamen"
label2struct[10]="Left-Thalamus-Proper"
label2struct[58]="Right-Accumbens-area"
label2struct[54]="Right-Amygdala"
label2struct[50]="Right-Caudate"
label2struct[53]="Right-Hippocampus"
label2struct[52]="Right-Pallidum"
label2struct[51]="Right-Putamen"
label2struct[49]="Right-Thalamus-Proper"


do_cortex_mapping() {
  local feat=$1
  local resolution=$2

  # Input & output locations
  surf_dir=${SUBJECT_MICAPIPE_DIR}/surf
  input_dir=${SUBJECT_MICAPIPE_DIR}/maps
  output_dir=${SUBJECT_OUTPUT_DIR}/${FOLDER_MAPS}/${FOLDER_CTX}

  # Mappings from features names to the way they appear in the input and output filenames
  declare -A map_input=([thickness]=thickness [flair]=white_flair [adc]=white_ADC [fa]=white_FA [qt1]=white_T1map)
  declare -A map_output=([thickness]=thickness [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)

  SHOW_INFO "Map '${feat}' to cortex [resolution=${resolution}]"

  # feat name in filenames
  input_feat="${map_input[${feat}]}"
  output_feat=${map_output[${feat}]}

  n=0
  for h in L R;
  do
    # Set paths
    surf_file="${surf_dir}/${BIDS_ID}_hemi-${h}_space-nativepro_surf-fsLR-${resolution}_label-white.surf.gii"
    input_file="${input_dir}/${BIDS_ID}_hemi-${h}_surf-fsLR-${resolution}_label-${input_feat}.func.gii"
    output_file="${output_dir}/${BIDS_ID}_hemi-${h}_surf-fsLR-${resolution}_feature-${output_feat}_smooth-${fwhm}mm.func.gii"

    # Check if file exists
    if [[ ! -f "${input_file}" ]]; then
      SHOW_WARNING "Subject ${BIDS_ID}: '${feat}' not available. Skipping..."
      break
    fi

    # Perform mapping
    DO_CMD wb_command -metric-smoothing "${surf_file}" "${input_file}" "${fwhm}" "${output_file}";
    [[ -f "$output_file" ]] && n=$((n + 1))
  done

  [[ $n -eq 2 ]] && SHOW_NOTE "Subject ${BIDS_ID}: '${feat}' [resolution=${resolution}] successfully mapped.";
}

do_subcortex_mapping (){
  local feat=$1

  SHOW_INFO "Map '${feat}' to subcortical structures"

  # Input & output locations
  vol_stats_file="${SUBJECT_SURF_DIR}/stats/aseg.stats"
  seg_file=${SUBJECT_MICAPIPE_DIR}/parc/${BIDS_ID}_space-nativepro_T1w_atlas-subcortical.nii.gz
  input_dir=${SUBJECT_MICAPIPE_DIR}/maps
  output_dir=${SUBJECT_OUTPUT_DIR}/${FOLDER_MAPS}/${FOLDER_SCTX}

  # Mappings from features names to the way they appear in the input and output filenames
  declare -A map_input=([flair]='map-flair' [adc]='model-DTI_map-ADC' [fa]='model-DTI_map-FA' [qt1]='map-T1map')
  declare -A map_output=([volume]=volume [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)

  declare -A label2mean_intensity=()
  output_file="${output_dir}/${BIDS_ID}_feature-${map_output[$feat]}.csv"

  # check that input files exist & if not volume, compute mean intensity using and store in label2mean_intensity
  if [[ "$feat" != "volume" ]]; then
    input_file="${input_dir}/${BIDS_ID}_space-nativepro_${map_input[$feat]}.nii.gz"
    if [[ ! -f "${input_file}" ]]; then
      SHOW_WARNING "Subject ${BIDS_ID}: '${feat}' not available. Skipping...";  return
    fi

    mapfile -t a < <(ImageIntensityStatistics 3 "${input_file}" "${seg_file}" | tail -n +2 | awk '{ print $1 "=" $2 }')
    for entry in "${a[@]}"; do
      IFS="=" read -r label mean_intensity <<< "$entry"
      label2mean_intensity+=([$label]=${mean_intensity})
    done
  else
    if [[ ! -f "${vol_stats_file}" ]]; then
      SHOW_WARNING "Subject ${BIDS_ID}: '${feat}' not available. Skipping..."; return
    fi
  fi

  # Write header and subject id to csv
  header="SubjID,Laccumb,Lamyg,Lcaud,Lhippo,Lpal,Lput,Lthal,Raccumb,Ramyg,Rcaud,Rhippo,Rpal,Rput,Rthal"
  [[ "$feat" == "volume" ]] && header+=",ICV"
  echo $header > "${output_file}"
  printf "%s,"  "${BIDS_ID}" >> "${output_file}"

  # get mean for each subcortical structure
  result=()
  for label in "${!label2struct[@]}"
  do
    if [[ "$feat" == "volume" ]]
    then
      struct_name=${label2struct[${label}]}
      value=$(grep "${struct_name}" "${vol_stats_file}" | awk '{ gsub(/^ *| *$/, "", $4); print $4 }')
    else
      value=${label2mean_intensity[${label}]}
    fi
    result+=("${value}")
  done

  # If volume -> add intracranial volume to csv
  if [[ "$feat" == "volume" ]]; then
    value=$(grep IntraCranialVol "${vol_stats_file}" | awk -F, '{ gsub(/^ *| *$/, "", $4); print $4 }')
    result+=("${value}")
  fi

  # write mean result to csv
  printf -v csv_line "%s," "${result[@]}"
  printf "%s\n" "${csv_line%?}" >> "${output_file}" # Remove last comma and add newline

  [[ -f "$output_file" ]] && SHOW_NOTE "Subject ${BIDS_ID}: '${feat}' successfully mapped.";
}

do_hippocampus_mapping() {
  local feat=$1
  local resolution=$2

  SHOW_INFO "Map '${feat}' to hippocampus [resolution=${resolution}]"

  # Input & output locations
  surf_dir=${SUBJECT_HIPPUNFOLD_DIR}/surf
  input_dir=${SUBJECT_MICAPIPE_DIR}/maps
  output_dir=${SUBJECT_OUTPUT_DIR}/${FOLDER_MAPS}/${FOLDER_HIPP}

  # Mappings from features names to the way they appear in the input, intermediate and output filenames
  declare -A map_input=([thickness]=thickness [flair]='map-flair' [adc]='model-DTI_map-ADC' [fa]='model-DTI_map-FA' [qt1]='map-T1map')
  declare -A map_inter=([thickness]=thickness [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)
  declare -A map_output=([thickness]=thickness [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)

  # feat name in filenames
  input_feat="${map_input[${feat}]}"
  inter_feat="${map_inter[${feat}]}"
  output_feat="${map_output[${feat}]}"

  n=0
  for h in L R;
  do

    # Set paths
    surf_file="${surf_dir}/${BIDS_ID}_hemi-${h}_space-T1w_den-${resolution}_label-hipp_midthickness.surf.gii"
    input_file="${input_dir}/${BIDS_ID}_space-nativepro_${input_feat}.nii.gz" # Not used for thickness
    if [[ "$feat" == "thickness" ]]; then
      inter_file="${surf_dir}/${BIDS_ID}_hemi-${h}_space-T1w_den-${resolution}_label-hipp_thickness.shape.gii"
    else
      inter_file="${tmp_dir}/${BIDS_ID}_hemi-${h}_space-T1w_desc-${inter_feat}_den-${resolution}_feature-hipp_midthickness.func.gii"
    fi
    output_file="${output_dir}/${BIDS_ID}_hemi-${h}_den-${resolution}_feature-${output_feat}_smooth-${fwhm}mm.func.gii"

    # Check if file exists
    [[ "$feat" != "thickness" ]] && check_file=${input_file} || check_file=${inter_file}
    if [[ ! -f "${check_file}" ]]; then
      SHOW_WARNING "Subject ${BIDS_ID}: '${feat}' not available. Skipping..."
      break
    fi

    # Perform mapping
    if [[ "$feat" != "thickness" ]]; then
      DO_CMD wb_command -volume-to-surface-mapping "${input_file}" "${surf_file}" "${inter_file}" -trilinear
    fi
    DO_CMD wb_command -metric-smoothing "${surf_file}" "${inter_file}" "${fwhm}" "${output_file}"

    [[ -f "${output_file}" ]] && n=$((n + 1))
  done

  [[ $n -eq 2 ]] && SHOW_NOTE "Subject ${BIDS_ID}: '${feat}' [resolution=${resolution}] successfully mapped.";
}


#------------------------------------------------------------------------------#
# Perform the mapping

# lowercase and split comma-separated into array
IFS=',' read -ra feat_list <<< "${comma_separated_features,,}"

# thickness -> volume
if [[ $structure == "subcortex" ]]; then feat_list=("${feat_list[@]/thickness/volume}"); fi

# do the mapping
for feat in "${feat_list[@]}";
do
  case $structure in
    cortex)
      IFS=',' read -ra resolution_list <<< "${comma_separated_resolutions}"
      for res in "${resolution_list[@]}"; do do_cortex_mapping "${feat}" "${res}"; done
      ;;
    subcortex)
      do_subcortex_mapping "${feat}";
      ;;
    hippocampus)
      IFS=',' read -ra resolution_list <<< "${comma_separated_resolutions}"
      for res in "${resolution_list[@]}"; do do_hippocampus_mapping "${feat}" "${res}"; done
      ;;
  esac
done


#------------------------------------------------------------------------------#
# Wrap up
SHOW_TITLE "${structure^} feature mapping ended in \033[38;5;220m$(printf "%.2f" "$(bc <<< "scale=2; $SECONDS/60")") minutes${COLOR_TITLE}
$([ -n "$logfile" ] && printf "\tCheck logs: %s" "\033[0;32m${logfile}")${COLOR_TITLE}"
