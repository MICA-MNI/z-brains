#!/bin/bash
#
# Feature mapping for cortex, subcortex, and hippocampus


if [[ -z ${ZBRAINS} ]]; then
  echo "ZBRAINS not defined"
  echo "This script must not be run standalone. Please use 'zbrains' directly."
  exit 0;
fi

script_dir=${ZBRAINS}/functions
source "${script_dir}/utilities.sh"


# Set umask
umask 003


#------------------------------------------------------------------------------#
args=("$@")
while (( "${#args[@]}" )); do
  option="${args[0]}"
  case "${option}" in
    --struct)
      PARSE_OPTION_SINGLE_VALUE structure args LIST_STRUCTURES || exit $?
      ;;
    --feat)
      PARSE_OPTION_MULTIPLE_VALUES features args LIST_FEATURES || exit $?
      ;;
    --fwhm)
      PARSE_OPTION_SINGLE_VALUE fwhm args || exit $?
      ;;
    --tmp)
      PARSE_OPTION_SINGLE_VALUE tmp_dir args || exit $?
      ;;
    --resolution)
      PARSE_OPTION_MULTIPLE_VALUES resolutions args LIST_RESOLUTIONS || exit $?
      ;;
    --labels)
      PARSE_OPTION_MULTIPLE_VALUES labels args || exit $?
      ;;
    *)
      SHOW_ERROR "Unknown option '${option}'"
      exit 1
      ;;
  esac
done

# Check if mandatory options are provided
ASSERT_REQUIRED "--struct" "${structure:-}"
ASSERT_REQUIRED "--feat" "${features[@]}"
ASSERT_REQUIRED "--tmp" "${tmp_dir:-}"


if [[ "$structure" != "subcortex" ]]; then
  ASSERT_REQUIRED "--fwhm" "${fwhm:-}"
  ASSERT_REQUIRED "--resolution" "${resolutions[@]:-}"
  ASSERT_REQUIRED "--labels" "${labels[@]:-}"
fi


#------------------------------------------------------------------------------#
declare -A map_struct=([cortex]=Cortical [subcortex]=Subcortical [hippocampus]=Hippocampal)

SHOW_TITLE "${map_struct[${structure}]} feature mapping: ${BIDS_ID}"


#------------------------------------------------------------------------------#
# Initialize timer
SECONDS=0


#------------------------------------------------------------------------------#
map_subcortex() {
  local feat=$1

  SHOW_INFO "${BIDS_ID}: Mapping '${feat}' to subcortical structures"

  # Mappings from features names to the way they appear in the input and output filenames
  declare -A map_input=([flair]='map-flair' [adc]='model-DTI_map-ADC' [fa]='model-DTI_map-FA' [qt1]='map-T1map')
  declare -A map_output=([volume]=volume [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)

  # Input & output locations
  aseg_stats_file="${SUBJECT_SURF_DIR}/stats/aseg.stats"
  seg_file=${SUBJECT_MICAPIPE_DIR}/parc/${BIDS_ID}_space-nativepro_T1w_atlas-subcortical.nii.gz
  input_dir=${SUBJECT_MICAPIPE_DIR}/maps

  #surf_dir=${SUBJECT_HIPPUNFOLD_DIR}/surf
  #input_dir=${SUBJECT_MICAPIPE_DIR}/maps
  if [[ $feat == plugin-* ]]; then 
    #aseg_stats_file="${SUBJECT_SURF_DIR}/stats/aseg.stats"
    #seg_file=${SUBJECT_MICAPIPE_DIR}/parc/${BIDS_ID}_space-nativepro_T1w_atlas-subcortical.nii.gz
    input_dir=${SUBJECT_PLUGIN_DIR}/maps
    feat=${feat:7}; 
  fi

  # feat name in filenames
  if [[ -v map_input[${feat,,}] ]];then
    input_feat="${map_input[${feat,,}]}"
    output_feat=${map_output[${feat,,}]}
  else
    input_feat=$feat
    output_feat=$feat
  fi

  

  #aseg_stats_file="${SUBJECT_SURF_DIR}/stats/aseg.stats"
  #seg_file=${SUBJECT_MICAPIPE_DIR}/parc/${BIDS_ID}_space-nativepro_T1w_atlas-subcortical.nii.gz
  #input_dir=${SUBJECT_MICAPIPE_DIR}/maps
  output_dir=${SUBJECT_OUTPUT_DIR}/${FOLDER_MAPS}/${FOLDER_SCTX}
  output_file="${output_dir}/${BIDS_ID}_feature-${output_feat}.csv"

  # check that input files exist & if not volume
  if [[ "$feat" != "volume" ]]; then
    input_file="${input_dir}/${BIDS_ID}_space-nativepro_${input_feat}.nii.gz"

    for file in ${seg_file} ${input_file}; do
      if [[ ! -f "${file}" ]]; then
        SHOW_WARNING "${BIDS_ID}: cannot map '${feat}' to subcortical structures." "Missing file: ${file}"
        return
      fi
    done

    DO_CMD "python ${script_dir}/subcortical_mapping.py -id ${BIDS_ID} -i ${input_file} -s ${seg_file} -o ${output_file}"

  else
    if [[ ! -f "${aseg_stats_file}" ]]; then
      SHOW_WARNING "${BIDS_ID}: cannot map '${feat}' to subcortical structures." "Missing file ${aseg_stats_file}"; return
    fi

    DO_CMD "python ${script_dir}/subcortical_mapping.py -id ${BIDS_ID} -v ${aseg_stats_file} -o ${output_file}"
  fi

  if [[ -f "$output_file" ]]; then
    SHOW_NOTE "${COLOR_INFO}${BIDS_ID}:" "'${feat}' successfully mapped.";
  else
    SHOW_WARNING "${BIDS_ID}: could not map '${feat}' to subcortical structures."
  fi
}



map_cortex() {
  local feat=$1
  local resol=$2
  local label="${3:-white}"

  root_dir=$SUBJECT_MICAPIPE_DIR  
  if [[ $feat == plugin-* ]]; then root_dir=$SUBJECT_PLUGIN_DIR; feat=${feat:7}; fi
    
  # Input & output locations
  surf_dir=${root_dir}/surf
  input_dir=${root_dir}/maps
  output_dir=${SUBJECT_OUTPUT_DIR}/${FOLDER_MAPS}/${FOLDER_CTX}

  # Mappings from features names to the way they appear in the input and output filenames
  declare -A map_feat=([thickness]=thickness [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)

  SHOW_INFO "${BIDS_ID}: Mapping '${feat}' to cortex [label=${label}, resolution=${resol}]"

  # feat name in filenames
  if [[ -v map_feat[${feat,,}] ]];then
    input_feat="${map_feat[${feat,,}]}"
    output_feat=${map_feat[${feat,,}]}
  else
    input_feat=$feat
    output_feat=$feat
  fi

  n=0
  for h in L R;
  do
    # Set paths
    surf_file="${surf_dir}/${BIDS_ID}_hemi-${h}_space-nativepro_surf-fsLR-${resol}_label-${label}.surf.gii"

    prefix="${BIDS_ID}_hemi-${h}_surf-fsLR-${resol}"
    if [[ "$feat" == "thickness" ]]; then
      input_file="${input_dir}/${prefix}_label-${input_feat}.func.gii"
      output_file="${output_dir}/${prefix}_label-${label}_feature-${output_feat}_smooth-${fwhm}mm.func.gii"
    else
      input_file="${input_dir}/${prefix}_label-${label}_${input_feat}.func.gii"
      output_file="${output_dir}/${prefix}_label-${label}_feature-${output_feat}_smooth-${fwhm}mm.func.gii"
    fi

    # Check if file exists
    for file in ${surf_file} ${input_file}; do
      if [[ ! -f "${file}" ]]; then
        SHOW_WARNING "${BIDS_ID}: cannot map '${feat}' [label=${label}, resolution=${resol}] to cortex." "Missing file: ${file}"
        return
      fi
    done

    # Flair
    # Normalize - req: tissue seg seg from free/fastsurfer and raw Flair, then normalize, the map
    # and continue below

    # Perform mapping
    DO_CMD "${WORKBENCH_PATH}/wb_command -metric-smoothing ${surf_file} ${input_file} ${fwhm} ${output_file}"
    [[ -f "$output_file" ]] && n=$((n + 1))
  done

  if [[ $n -eq 2 ]]; then
    SHOW_NOTE "${COLOR_INFO}${BIDS_ID}:" "'${feat}' [label=${label}, resolution=${resol}] successfully mapped."
  else
    SHOW_WARNING "${BIDS_ID}: could not map '${feat}' [label=${label}, resolution=${resol}] to cortex."
  fi
}


map_hippocampus() {
  local feat=$1
  local resol=$2
  local label="${3:-midthickness}"

  SHOW_INFO "${BIDS_ID}: Mapping '${feat}' to hippocampus [label=${label}, resolution=${resol}]"

  # Input & output locations
  surf_dir=${SUBJECT_HIPPUNFOLD_DIR}/surf
  input_dir=${SUBJECT_MICAPIPE_DIR}/maps
  ISSURF=false
  if [[ $feat == plugin-* ]]; then 
    surf_dir=$SUBJECT_PLUGIN_DIR/surf; 
    input_dir=$SUBJECT_PLUGIN_DIR/maps; 
    feat=${feat:7}; 
  fi
  output_dir=${SUBJECT_OUTPUT_DIR}/${FOLDER_MAPS}/${FOLDER_HIP}

  # Mappings from features names to the way they appear in the input, intermediate and output filenames
  declare -A map_input=([thickness]=thickness [flair]='map-flair' [adc]='model-DTI_map-ADC' [fa]='model-DTI_map-FA'
                        [qt1]='map-T1map')
  declare -A map_inter=([thickness]=thickness [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)
  declare -A map_output=([thickness]=thickness [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)

  # feat name in filenames
  if [[ -v map_input[${feat,,}] ]]; then
    input_feat="${map_input[${feat,,}]}"
    inter_feat="${map_inter[${feat,,}]}"
    output_feat="${map_output[${feat,,}]}"
  else
    input_feat=$feat
    inter_feat=$feat
    output_feat=$feat
  fi

  n=0
  for h in L R;
  do

  # Set paths
  prefix="${BIDS_ID}_hemi-${h}"
  surf_file="${surf_dir}/${prefix}_space-T1w_den-${resol}_label-hipp_${label}.surf.gii"
  input_file="${input_dir}/${BIDS_ID}_space-nativepro_${input_feat}.nii.gz" # Not used for thickness
  output_file="${output_dir}/${prefix}_den-${resol}_label-${label}_feature-${output_feat}_smooth-${fwhm}mm.func.gii"
  # Note the logic here is that if a [shape|func].gii exists, use that. Otherwise map a .nii.gz file
  if [ -f "${surf_dir}/${prefix}_space-T1w_den-${resol}_label-hipp_${feat}.*.gii" ]; then
    inter_file=$(ls "${surf_dir}/${prefix}_space-T1w_den-${resol}_label-hipp_${feat}.*.gii")
    ISSURF=true
  elif [ -f "${input_dir}/${prefix}_space-T1w_den-${resol}_label-hipp_${feat}.*.gii" ]; then
    inter_file=$(ls "${input_dir}/${prefix}_space-T1w_den-${resol}_label-hipp_${feat}.*.gii")
    ISSURF=true
  else
    inter_file="${tmp_dir}/${prefix}_space-T1w_desc-${inter_feat}_den-${resol}_feature-hipp_${label}.func.gii"
    ISSURF=false
  fi

  # Check if file exists
  if $ISSURF; then
    check_file=${inter_file}
  else
    check_file=${input_file}
  fi
  for file in ${surf_file} ${check_file}; do
    if [[ ! -f "${file}" ]]; then
      SHOW_WARNING "${BIDS_ID}: cannot map '${feat}' [label=${label}, resolution=${resol}] to hippocampus." "Missing file: ${file}"
      return
    fi
  done

    # Perform mapping
    if ! $ISSURF; then
      cmd="${WORKBENCH_PATH}/wb_command -volume-to-surface-mapping ${input_file} ${surf_file} ${inter_file} -trilinear"
      DO_CMD "$cmd"
    fi

#    if [[ ! -f "${inter_file}" ]]; then
#      SHOW_WARNING "${BIDS_ID}: cannot map '${feat}' feature." "Missing file: ${inter_file}"
#      return
#    fi
    DO_CMD "${WORKBENCH_PATH}/wb_command -metric-smoothing ${surf_file} ${inter_file} ${fwhm} ${output_file}"

    [[ -f "${output_file}" ]] && n=$((n + 1))
  done

  if [[ $n -eq 2 ]]; then
    SHOW_NOTE "${COLOR_INFO}${BIDS_ID}:" "'${feat}' [label=${label}, resolution=${resol}] successfully mapped."
  else
    SHOW_WARNING "${BIDS_ID}: could not map '${feat}' [label=${label}, resolution=${resol}] to hippocampus."
  fi

}


#------------------------------------------------------------------------------#
# Perform the mapping

# thickness -> volume
[[ $structure == "subcortex" ]] && features=("${features[@]/thickness/volume}")

# do the mapping
for feat in "${features[@]}";
do
  case $structure in
    cortex)
      for res in "${resolutions[@]}"; do
        for lab in "${labels[@]}"; do map_cortex "${feat}" "${res}" "${lab}"; done
      done
      ;;
    subcortex)
      map_subcortex "${feat}";
      ;;
    hippocampus)
      for res in "${resolutions[@]}"; do
        for lab in "${labels[@]}"; do map_hippocampus "${feat}" "${res}" "${lab}"; done
      done
      ;;
  esac
done


#------------------------------------------------------------------------------#
# Wrap up
elapsed=$(printf "%.2f" "$(bc <<< "scale=2; $SECONDS/60")")
SHOW_TITLE "${map_struct[${structure}]} feature mapping for ${BIDS_ID} ended in \033[38;5;220m${elapsed} minutes${NO_COLOR}"
