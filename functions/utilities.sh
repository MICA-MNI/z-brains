#!/bin/bash

#source ../config.cfg


bids_variables() {
  # This functions assigns variable names acording to:
  #     Subject ID =                        $1
  #     Raw data directory =                $2
  #     Hippunfold directory derivatives =  $3
  #     Micapipe derivatives directory =    $4
  #     Z-brains derivatives directory =    $5
  #     Session ID (optional) =             $6
  id=$1 # id in the for HC001
  raw_dir=$2
  hipp_dir=$3
  micapipe_dir=$4
  out_dir=$5
  session=$6
  umask 002

  # Define utilities/functions directory
  export SCRIPT_DIR=${ZBRAINS}/functions

  export SUBJECT=sub-${id}

  # Handle Single Session
  if [ "$session" == "SINGLE" ]; then
      export SUBJECT_RAW=${raw_dir}/${SUBJECT}           # subject's raw data directory
      export SUBJECT_HIPP=${hipp_dir}/hippundold/${SUBJECT}         # subject's hippunfold directory
      export SUBJECT_MICAPIPE=${micapipe_dir}/${SUBJECT} # subject's micapipe directory
      export SUBJECT_OUT=${out_dir}/${SUBJECT}           # subject's output directory
      session_suffix=""
  else
    predix_session="ses-${session/ses-/}";
    export SUBJECT_RAW=${raw_dir}/${SUBJECT}/${predix_session}           # subject's raw data directory
    export SUBJECT_HIPP=${hipp_dir}/hippunfold/${SUBJECT}/${predix_session}         # subject's hippunfold directory
    export SUBJECT_MICAPIPE=${micapipe_dir}/${SUBJECT}/${predix_session} # subject's micapipe directory
    export SUBJECT_OUT=${out_dir}/${SUBJECT}/${predix_session}           # subject's output directory
    session_suffix="_${predix_session}"
  fi

  export BIDS_ID="${SUBJECT}${session_suffix}"

  # Structural directories derivatives/
  export PROC_STRUCT=$SUBJECT_MICAPIPE/anat  # Structural processing directory
  export DIR_VOLUM=$SUBJECT_MICAPIPE/parc    # Cortical segmentations
  export DIR_CONTE69=$SUBJECT_MICAPIPE/surf  # Conte69 surfaces
  export DIR_MAPS=$SUBJECT_MICAPIPE/maps

  export DIR_WARP=$SUBJECT_MICAPIPE/xfm     # Transformation matrices
  export DIR_LOGS=$SUBJECT_OUT/logs         # directory where log files will be written
  export DIR_QC=$SUBJECT_MICAPIPE/QC        # directory with micapipe qc files

  # post structural Files (the resolution might vary depending on the dataset)
  RES=$(mrinfo "${PROC_STRUCT}"/"${BIDS_ID}"_space-nativepro_T1w.nii.gz -spacing | awk '{printf "%.1f\n", $2}')
  export RES
  export T1NATIVEPRO=${PROC_STRUCT}/${BIDS_ID}_space-nativepro_T1w.nii.gz
  export T1NATIVEPRO_BRAIN=${PROC_STRUCT}/${BIDS_ID}_space-nativepro_T1w_brain.nii.gz
  export T1NATIVEPRO_MASK=${PROC_STRUCT}/${BIDS_ID}_space-nativepro_T1w_brain_mask.nii.gz
  export T1FAST_SEG=${DIR_VOLUM}/${BIDS_ID}_space-nativepro_T1w_atlas-subcortical.nii.gz

  # raw files
  RAW_FLAIR=$(ls "$SUBJECT_RAW"/anat/*FLAIR*.nii* 2>/dev/null)
  export RAW_FLAIR
}

set_surface_directory() {
  local recon=${1}
  DIR_SURF=$(dirname "${out_dir}")/${recon}    # surf
  export DIR_SURF
  export DIR_SUBJSURF=${DIR_SURF}/${BIDS_ID}  # Subject surface dir
  export T1SURF=${DIR_SUBJSURF}/mri/orig.mgz

  # Native midsurface in gifti format
  export LH_MIDSURF="${DIR_CONTE69}/${BIDS_ID}_hemi-L_surf-fsnative_label-midthickness.surf.gii"
  export RH_MIDSURF="${DIR_CONTE69}/${BIDS_ID}_hemi-R_surf-fsnative_label-midthickness.surf.gii"
}

bids_print.variables-ctx() {
  # This functions prints BIDS variables names and files if found
  Info "Cortical processing variables:"
  Note "T1 nativepro         :" "$(find "$T1NATIVEPRO" 2>/dev/null)"
  Note "T1 resolution        :" "$RES"
  Note "fs-LR-32k-lh surface :" "$(find "${DIR_CONTE69}/${BIDS_ID}_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii" 2>/dev/null)"
  Note "fs-LR-32k-rh surface :" "$(find "${DIR_CONTE69}/${BIDS_ID}_hemi-R_space-nativepro_surf-fsLR-32k_label-white.surf.gii" 2>/dev/null)"

}

bids_print() {
  # This functions prints BIDS variables names and files if found

  struct=$1

  case $struct in
    "cortex")
      Info "Cortical processing variables:"
      Note "T1 nativepro         :" "$(find "$T1NATIVEPRO" 2>/dev/null)"
      Note "T1 resolution        :" "$RES"
      Note "fs-LR-32k-lh surface :" "$(find "${DIR_CONTE69}/${BIDS_ID}_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii" 2>/dev/null)"
      Note "fs-LR-32k-rh surface :" "$(find "${DIR_CONTE69}/${BIDS_ID}_hemi-R_space-nativepro_surf-fsLR-32k_label-white.surf.gii" 2>/dev/null)"
      ;;
    "subcortex")
      Info "Subcortical processing variables:"
      Note "T1 nativepro    :" "$(find "$T1NATIVEPRO" 2>/dev/null)"
      Note "T1 sctx seg     :" "$(find "$T1FAST_SEG" 2>/dev/null)"
      Note "Sctx volume     :" "$(find "$DIR_SUBJSURF/stats/aseg.stats" 2>/dev/null)"
      Note "T1 resolution   :" "$RES"
      ;;
    "hippocampus")
      Info "Hippocampal processing variables:"
      Note "T1 nativepro    :" "$(find "$T1NATIVEPRO" 2>/dev/null)"
      Note "T1 resolution   :" "$RES"
      ;;
    *)  # Default case if none of the above patterns match
      echo "Unknown structure in bids_print."
      ;;
  esac
}

bids_print.variables-hipp() {
  # This functions prints BIDS variables names and files if found
  Info "Hippocampal processing variables:"
  Note "T1 nativepro    :" "$(find "$T1NATIVEPRO" 2>/dev/null)"
  Note "T1 resolution   :" "$RES"
}

bids_variables_unset() {
  # This function unsets all environment variables defined by bids_variables
  unset SCRIPT_DIR
  unset SUBJECT
  unset SUBJECT_RAW
  unset SUBJECT_HIPP
  unset SUBJECT_MICAPIPE
  unset SUBJECT_OUT
  unset PROC_STRUCT
  unset DIR_VOLUM
  unset DIR_CONTE69
  unset DIR_MAPS
  unset DIR_WARP
  unset DIR_LOGS
  unset DIR_QC
  unset RES
  unset T1NATIVEPRO
  unset T1NATIVEPRO_BRAIN
  unset T1NATIVEPRO_MASK
  unset T1FAST_SEG
  unset RAW_FLAIR
  unset DIR_SURF
  unset DIR_SUBJSURF
  unset T1SURF
  unset LH_MIDSURF
  unset RH_MIDSURF
}

function zbrains_json() {
  # Name is the name of the raw-BIDS directory
  if [ -f "${BIDS}/dataset_description.json" ]; then
    Name=$(grep Name "${BIDS}"/dataset_description.json | awk -F '"' '{print $4}')
    BIDSVersion=$(grep BIDSVersion "${BIDS}"/dataset_description.json | awk -F '"' '{print $4}')
  else
    Name="BIDS dataset_description NOT found"
    BIDSVersion="BIDS dataset_description NOT found"
  fi

  echo -e "{
    \"Name\": \"${Name}\",
    \"BIDSVersion\": \"${BIDSVersion}\",
    \"DatasetType\": \"derivative\",
    \"GeneratedBy\": [{
        \"Name\": \"z-brains\",
        \"Version\": \"${VERSION}\",
        \"Reference\": \"tbd\",
        \"DOI\": \"tbd\",
        \"URL\": \"https://z-brains.readthedocs.io\",
        \"GitHub\": \"https://github.com/MICA-MNI/z-brains\",
        \"Container\": {
          \"Type\": \"tbd\",
          \"Tag\": \"ello:zbrains:$(echo "${VERSION}" | awk '{print $1}')\"
        },
        \"RunBy\": \"$(whoami)\",
        \"Workstation\": \"$(uname -n)\",
        \"LastRun\": \"$(date)\",
        \"Processing\": \"${PROC}\"
      }]
  }" > "${out_dir}/dataset_description.json"
}

zbrains_software() {
  Info "Z-BRAINS - Software versions"
  Note "ANTS........" "$(antsRegistration --version | awk -F ':' 'NR==1{print $2}')"
  Note "            " "$ANTSPATH"
  Note "WorkBench..." "$(wb_command -version | awk 'NR==3{print $2}')"
  Note "            " "$(which wb_command)"
  Note "python......" "$(python --version | awk 'NR==1{print $2}')"
  Note "            " "$(which python)"
  Note "conda......." "$(conda --version)"
  Note "            " "$(which conda)"
}

function cleanup() {
  # This script will clean the temporary directory
  # and reset the old user path upon interrupt and termination.
  tmp_dir=$1
#  nocleanup=$2

  rm -Rf "${tmp_dir:?}" 2>/dev/null
  bids_variables_unset

#  # Clean temporary directory and temporary fsaverage5
#  if [[ $nocleanup == "FALSE" ]]; then
#      rm -Rf "${tmp_dir:?}" 2>/dev/null
#  else
#      echo -e "z-brains tmp directory was not erased: \n\t\t${tmp_dir}";
#  fi
##  cd "$here" || exit
#  bids_variables_unset
  if [[ -n "$OLD_PATH" ]]; then export PATH=$OLD_PATH; unset OLD_PATH; fi
}


# ----------------------------------------------------------------------------------------------- #
# ------------------------------- Printing and display functions -------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# The following functions are only to print on the terminal colorful messages:
#     Error messages
#     Warning messages
#     Note messages
#     Warn messages
#     Title messages

COLOR_ERROR="\033[38;5;9m"
COLOR_NOTE="\033[38;5;122m"
COLOR_INFO="\033[38;5;75m"
COLOR_WARNING="\033[38;5;184m"
COLOR_TITLE="\033[38;5;141m"
COLOR_NOTE="\033[0;36;10m"
NC="\033[0m" # No color

Error() {
echo -e "${COLOR_ERROR}\n-------------------------------------------------------------\n\n[ ERROR ]..... $1\n
-------------------------------------------------------------${NC}\n"
}
Note() {
if [[ ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then echo -e "\t\t$1\t${COLOR_NOTE}$2 ${NC}"; fi
}
Info() {
if [[ ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then echo  -e "${COLOR_INFO}\n[ INFO ]..... $1 ${NC}"; fi
}
Warning() {
if [[ ${VERBOSE} -gt 0 || ${VERBOSE} -lt 0 ]]; then echo  -e "${COLOR_WARNING}\n[ WARNING ]..... $1 ${NC}"; fi
}
Warn() {
if [[ ${VERBOSE} -gt 0 || ${VERBOSE} -lt 0 ]]; then echo  -e "${COLOR_WARNING}
-------------------------------------------------------------\n
[ WARNING ]..... $1
\n-------------------------------------------------------------${NC}"; fi
}
Title() {
if [[ ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then echo -e "\n${COLOR_TITLE}
-------------------------------------------------------------
\t$1
-------------------------------------------------------------${NC}"; fi
}


# Export
export -f Error Note Info Warning Warn Title


function Do_cmd() {
# do_cmd sends command to stdout before executing it.
str="$(whoami) @ $(uname -n) $(date)"
local l_command=""
local l_sep=" "
local l_index=1

while [ ${l_index} -le $# ]; do
    eval "arg=\${$l_index}"
    if [ "$arg" = "-fake" ]; then
      arg=""
    fi
    if [ "$arg" = "-no_stderr" ]; then
      arg=""
    fi
    if [ "$arg" == "-log" ]; then
      next_arg=$(("${l_index}" + 1))
      eval "logfile=\${${next_arg}}"
      arg=""
      l_index=$((l_index+1))
    fi
    l_command="${l_command}${l_sep}${arg}"
    l_sep=" "
    l_index=$((l_index+1))
   done

#[[ ${QUIET} != true ]] && echo -e "\033[38;5;118m\n${str}:\nCOMMAND -->  \033[38;5;122m${l_command}  \033[0m";
if [[ ${VERBOSE} -gt 2 || ${VERBOSE} -lt 0 ]]; then
  echo -e "\033[38;5;118m\n${str}:\nCOMMAND -->  \033[38;5;122m${l_command}  \033[0m";
fi
if [ -z "$TEST" ]; then $l_command; fi
}

