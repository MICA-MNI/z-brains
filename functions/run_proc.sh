#!/bin/bash
#
# Feature mapping for cortex, subcortex, and hippocampus


if [[ -z ${ZBRAINS} ]]; then
  echo "ZBRAINS not defined"
  echo "This script should not be run standalone. Please use 'zbrains' directly."
  exit 0;
fi


# Set umask
umask 003


#------------------------------------------------------------------------------#
while (( "$#" )); do
  args=("$@")
  case "$1" in
    --struct)
      structure=$(PARSE_OPTION_SINGLE_VALUE args LIST_STRUCTURES)
      shift 2
      ;;
    --feat)
      readarray -t features < <(PARSE_OPTION_MULTIPLE_VALUES args LIST_FEATURES)
      shift $(( ${#features[@]} + 1 ))
      ;;
    --fwhm)
      fwhm=$(PARSE_OPTION_SINGLE_VALUE args)
      shift 2
      ;;
    --tmp)
      tmp_dir=$(PARSE_OPTION_SINGLE_VALUE args)
      shift 2
      ;;
    --resolution)
      readarray -t resolutions < <(PARSE_OPTION_MULTIPLE_VALUES args LIST_RESOLUTIONS)
      shift $(( ${#resolutions[@]} + 1 ))
      ;;
    --logfile)
      logfile=$(PARSE_OPTION_SINGLE_VALUE args)
      shift 2
      ;;
    --labels)
      readarray -t labels < <(PARSE_OPTION_MULTIPLE_VALUES args)
      shift $(( ${#labels[@]} + 1 ))
      ;;
    *)
      SHOW_ERROR "Unknown option '$1'"
      exit 1
      ;;
  esac
done

# Check if mandatory options are provided
ASSERT_REQUIRED "--struct" "$structure"
ASSERT_REQUIRED "--feat" "${features[@]}"
ASSERT_REQUIRED "--tmp" "$tmp_dir"


if [[ "$structure" != "subcortex" ]]; then
  ASSERT_REQUIRED "--fwhm" "$fwhm"
  ASSERT_REQUIRED "--resolution" "${resolutions[@]}"
  ASSERT_REQUIRED "--labels" "${labels[@]}"
fi


#------------------------------------------------------------------------------#
declare -A map_struct=([cortex]=Cortical [subcortex]=Subcortical [hippocampus]=Hippocampal)

SHOW_TITLE "${map_struct[${structure}]} feature mapping"


#------------------------------------------------------------------------------#
# Initialize timer
SECONDS=0


#------------------------------------------------------------------------------#
map_subcortex() {
  local feat=$1

  SHOW_INFO "Map '${feat}' to subcortical structures"

  # Mappings from features names to the way they appear in the input and output filenames
  declare -A map_input=([flair]='map-flair' [adc]='model-DTI_map-ADC' [fa]='model-DTI_map-FA' [qt1]='map-T1map')
  declare -A map_output=([volume]=volume [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)

  # Input & output locations
  aseg_stats_file="${SUBJECT_SURF_DIR}/stats/aseg.stats"
  seg_file=${SUBJECT_MICAPIPE_DIR}/parc/${BIDS_ID}_space-nativepro_T1w_atlas-subcortical.nii.gz
  input_dir=${SUBJECT_MICAPIPE_DIR}/maps
  output_dir=${SUBJECT_OUTPUT_DIR}/${FOLDER_MAPS}/${FOLDER_SCTX}
  output_file="${output_dir}/${BIDS_ID}_feature-${map_output[$feat]}.csv"

  # check that input files exist & if not volume
  if [[ "$feat" != "volume" ]]; then
    input_file="${input_dir}/${BIDS_ID}_space-nativepro_${map_input[$feat]}.nii.gz"
    if [[ ! -f "${input_file}" ]]; then
      SHOW_WARNING "Subject ${BIDS_ID}: '${feat}' not available. Skipping...";  return
    fi

    DO_CMD "python ${ZBRAINS}/functions/subcortical_mapping.py -id ${BIDS_ID} -i ${input_file} -s ${seg_file} \
                                                               -o ${output_file}"

  else
    if [[ ! -f "${aseg_stats_file}" ]]; then
      SHOW_WARNING "Subject ${BIDS_ID}: '${feat}' not available. Skipping..."; return
    fi

    DO_CMD "python ${ZBRAINS}/functions/subcortical_mapping.py -id ${BIDS_ID} -v ${aseg_stats_file} -o ${output_file}"
  fi

  [[ -f "$output_file" ]] && SHOW_NOTE "Subject ${BIDS_ID}: '${feat}' successfully mapped.";
}



map_cortex() {
  local feat=$1
  local resolution=$2
  local label="${3:-white}"

  # Input & output locations
  surf_dir=${SUBJECT_MICAPIPE_DIR}/surf
  input_dir=${SUBJECT_MICAPIPE_DIR}/maps
  output_dir=${SUBJECT_OUTPUT_DIR}/${FOLDER_MAPS}/${FOLDER_CTX}

  # Mappings from features names to the way they appear in the input and output filenames
  declare -A map_feat=([thickness]=thickness [flair]=flair [adc]=ADC [fa]=FA [qt1]=T1map)

  SHOW_INFO "Map '${feat}' to cortex [resolution=${resolution}]"

  # feat name in filenames
  input_feat="${map_feat[${feat}]}"
  output_feat=${map_feat[${feat}]}

  n=0
  for h in L R;
  do
    # Set paths
    surf_file="${surf_dir}/${BIDS_ID}_hemi-${h}_space-nativepro_surf-fsLR-${resolution}_label-${label}.surf.gii"
    if [[ "$feat" == "thickness" ]]; then
      input_file="${input_dir}/${BIDS_ID}_hemi-${h}_surf-fsLR-${resolution}_label-${input_feat}.func.gii"
      output_file="${output_dir}/${BIDS_ID}_hemi-${h}_surf-fsLR-${resolution}_feature-${output_feat}\
      _smooth-${fwhm}mm.func.gii"
    else
      input_file="${input_dir}/${BIDS_ID}_hemi-${h}_surf-fsLR-${resolution}_label-${label}_${input_feat}.func.gii"
      output_file="${output_dir}/${BIDS_ID}_hemi-${h}_surf-fsLR-${resolution}_label-${label}_feature-${output_feat}\
      _smooth-${fwhm}mm.func.gii"
    fi

    # Check if file exists
    if [[ ! -f "${input_file}" ]]; then
      SHOW_WARNING "Subject ${BIDS_ID}: '${feat}' not available. Skipping..."
      break
    fi

    # Flair
    # Normalize - req: tissue seg seg from free/fastsurfer and raw Flair, then normalize, the map
    # and continue below

    # Perform mapping
    DO_CMD "${WORKBENCH_PATH}/wb_command -metric-smoothing ${surf_file} ${input_file} ${fwhm} ${output_file}"
    [[ -f "$output_file" ]] && n=$((n + 1))
  done

  [[ $n -eq 2 ]] && SHOW_NOTE "Subject ${BIDS_ID}: '${feat}' [resolution=${resolution}] successfully mapped.";
}


map_hippocampus() {
  local feat=$1
  local resolution=$2
  local label="${3:-midthickness}"

  SHOW_INFO "Map '${feat}' to hippocampus [resolution=${resolution}]"

  # Input & output locations
  surf_dir=${SUBJECT_HIPPUNFOLD_DIR}/surf
  input_dir=${SUBJECT_MICAPIPE_DIR}/maps
  output_dir=${SUBJECT_OUTPUT_DIR}/${FOLDER_MAPS}/${FOLDER_HIPP}

  # Mappings from features names to the way they appear in the input, intermediate and output filenames
  declare -A map_input=([thickness]=thickness [flair]='map-flair' [adc]='model-DTI_map-ADC' [fa]='model-DTI_map-FA'
                        [qt1]='map-T1map')
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
    surf_file="${surf_dir}/${BIDS_ID}_hemi-${h}_space-T1w_den-${resolution}_label-hipp_${label}.surf.gii"
    input_file="${input_dir}/${BIDS_ID}_space-nativepro_${input_feat}.nii.gz" # Not used for thickness
    if [[ "$feat" == "thickness" ]]; then
      inter_file="${surf_dir}/${BIDS_ID}_hemi-${h}_space-T1w_den-${resolution}_label-hipp_thickness.shape.gii"
    else
      inter_file="${tmp_dir}/${BIDS_ID}_hemi-${h}_space-T1w_desc-${inter_feat}_den-${resolution}\
      _feature-hipp_${label}.func.gii"
    fi
    output_file="${output_dir}/${BIDS_ID}_hemi-${h}_den-${resolution}_label-${label}_feature-${output_feat}\
    _smooth-${fwhm}mm.func.gii"

    # Check if file exists
    [[ "$feat" != "thickness" ]] && check_file=${input_file} || check_file=${inter_file}
    if [[ ! -f "${check_file}" ]]; then
      SHOW_WARNING "Subject ${BIDS_ID}: '${feat}' not available. Skipping..."
      break
    fi

    # Perform mapping
    if [[ "$feat" != "thickness" ]]; then
      DO_CMD "${WORKBENCH_PATH}/wb_command -volume-to-surface-mapping ${input_file} ${surf_file} ${inter_file} \
                                                                      -trilinear"
    fi
    DO_CMD "${WORKBENCH_PATH}/wb_command -metric-smoothing ${surf_file} ${inter_file} ${fwhm} ${output_file}"

    [[ -f "${output_file}" ]] && n=$((n + 1))
  done

  [[ $n -eq 2 ]] && SHOW_NOTE "Subject ${BIDS_ID}: '${feat}' [resolution=${resolution}] successfully mapped.";
}


#------------------------------------------------------------------------------#
# Perform the mapping

# thickness -> volume
if [[ $structure == "subcortex" ]]; then feat_list=("${features[@]/thickness/volume}"); fi

# do the mapping
for feat in "${feat_list[@]}";
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
SHOW_TITLE "${structure^} feature mapping ended in \033[38;5;220m\
$(printf "%.2f" "$(bc <<< "scale=2; $SECONDS/60")") minutes${COLOR_TITLE}
$([ -n "$logfile" ] && printf "\tCheck logs: %s" "\033[0;32m${logfile}")${COLOR_TITLE}"
