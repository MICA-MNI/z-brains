#!/bin/bash

export VERBOSE=-1  # Default

function DO_CMD() {
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

if [[ ${VERBOSE} -gt 2 || ${VERBOSE} -lt 0 ]]; then
  echo -e "\033[38;5;118m\n${str}:\nCOMMAND -->  \033[38;5;122m${l_command}  \033[0m";
fi

$l_command
}


#function zbrains_json() {
#  # Name is the name of the raw-BIDS directory
#  if [ -f "${BIDS}/dataset_description.json" ]; then
#    Name=$(grep Name "${BIDS}"/dataset_description.json | awk -F '"' '{print $4}')
#    BIDSVersion=$(grep BIDSVersion "${BIDS}"/dataset_description.json | awk -F '"' '{print $4}')
#  else
#    Name="BIDS dataset_description NOT found"
#    BIDSVersion="BIDS dataset_description NOT found"
#  fi
#
#  echo -e "{
#    \"Name\": \"${Name}\",
#    \"BIDSVersion\": \"${BIDSVersion}\",
#    \"DatasetType\": \"derivative\",
#    \"GeneratedBy\": [{
#        \"Name\": \"zbrains\",
#        \"Version\": \"${VERSION}\",
#        \"Reference\": \"tbd\",
#        \"DOI\": \"tbd\",
#        \"URL\": \"https://z-brains.readthedocs.io\",
#        \"GitHub\": \"https://github.com/MICA-MNI/z-brains\",
#        \"Container\": {
#          \"Type\": \"tbd\",
#          \"Tag\": \"ello:zbrains:$(echo "${VERSION}" | awk '{print $1}')\"
#        },
#        \"RunBy\": \"$(whoami)\",
#        \"Workstation\": \"$(uname -n)\",
#        \"LastRun\": \"$(date)\",
#        \"Processing\": \"${PROC}\"
#      }]
#  }" > "${out_dir}/dataset_description.json"
#}


# ----------------------------------------------------------------------------------------------- #
# ------------------------------- Printing and display functions -------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# The following functions are only to print on the terminal colorful messages:
#     Error messages
#     Warning messages
#     Note messages
#     Info messages
#     Title messages
# VERBOSE is defined by the zbrains main script

COLOR_ERROR="\033[38;5;9m"
COLOR_WARNING="\033[38;5;184m"
COLOR_INFO="\033[38;5;75m"
COLOR_NOTE="\033[0;36;10m"
COLOR_TITLE="\033[38;5;141m"

NO_COLOR="\033[0m" # No color

SHOW_ERROR() {
echo -e "${COLOR_ERROR}\n-------------------------------------------------------------\n\n[ ERROR ]..... $1\n
-------------------------------------------------------------${NO_COLOR}\n"
}

SHOW_WARNING() {
  if [[ -z ${VERBOSE} || ${VERBOSE} -gt 0 || ${VERBOSE} -lt 0 ]]; then echo  -e "${COLOR_WARNING}\n[ WARNING ]..... $1 ${NO_COLOR}"; fi
}

SHOW_NOTE() {
  if [[ -z ${VERBOSE} || ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then echo -e "\t\t$1\t${COLOR_NOTE}$2 ${NO_COLOR}"; fi
}

SHOW_INFO() {
  if [[ -z ${VERBOSE} || ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then echo  -e "${COLOR_INFO}\n[ INFO ]..... $1 ${NO_COLOR}"; fi
}

SHOW_TITLE() {
if [[ -z ${VERBOSE} || ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then echo -e "\n${COLOR_TITLE}
-------------------------------------------------------------
\t$1
-------------------------------------------------------------${NO_COLOR}"
fi
}


function allowed_to_regex() {
local array=("$@")
    array=("${array[@]/#/ }")  # Adds a leading space
    array=("${array[@]/%/ }")  # Adds a trailing space
    IFS='|'; echo "${array[*]}"; unset IFS
}


PARSE_OPTION_SINGLE_VALUE() {
  # PARSE_OPTION_SINGLE_VALUE args allowed_values
  # or
  # PARSE_OPTION_SINGLE_VALUE args
  local -n _args="$1"
  local -n _allowed_values

  [[ $# -gt 1 ]] && _allowed_values="$2" && allowed_regex=$(allowed_to_regex "${_allowed_values[@]}")

  local option="${_args[0]}"
  local value="${_args[1]}"

  # Check if the option has a value and if the value is not another option
  if [[ -z "$value" || "$value" == --* ]]; then
      SHOW_ERROR "${option} option requires a value."
      exit 1
  fi

  # Check if value is in the list of allowed values
  if [[ -v allowed_regex && ! " ${value} " =~ ${allowed_regex} ]]; then
      SHOW_ERROR "Invalid value '$value' for ${option} option. \nAllowed values are: [${_allowed_values[*]}]."
      exit 1
  fi

  echo "$value"
}


PARSE_OPTION_MULTIPLE_VALUES() {
  # PARSE_OPTION_MULTIPLE_VALUES args allowed_values all
  # or
  # PARSE_OPTION_MULTIPLE_VALUES args allowed_values
  # or
  # PARSE_OPTION_MULTIPLE_VALUES args
  local -n _args="$1"
  local -n _allowed_values
  local all  # accept all as a value, additional to those in allowed_values -> denoting all allowed values

  # Check if allowed values and "all" option are provided
  [[ $# -gt 1 ]] && _allowed_values="$2";
  [[ $# -gt 2 ]] && all="$3" && _allowed_values+=("$all");
  [[ -v _allowed_values ]] && allowed_regex=$(allowed_to_regex "${_allowed_values[@]}")

  local option="${_args[0]}"
  _args=("${_args[@]:1}")

  local values=()
  while [[ ${#_args[@]} -gt 0 && "${_args[0]}" != --* ]]; do
    local value="${_args[0]}"

    # Check if the value is in the list of allowed values
    if [[ -v allowed_regex && ! " ${value} " =~ ${allowed_regex} ]]; then
      SHOW_ERROR "Invalid value '$value' for ${option} option. \nAllowed values are: [${_allowed_values[*]}]."
      exit 1
    fi

    values+=("$value")
    _args=("${_args[@]:1}")  # Shift 2 elements from _args
  done

  [[ ${#values[@]} == 0 ]] && SHOW_ERROR "${option} option requires at least one value." && exit 1;

  # If "all" option is set, replace values with all allowed values
  #  if [[ -v all || " ${values[*]} " != *" $all "* ]]; then
  #    values=("${_allowed_values[@]}")
  #  fi

  printf "%s\n" "${values[@]}"
}


ASSERT_REQUIRED() {
  local option="$1"
  local value="$2"
  local error_message="${3:-$option option is required}"

  [[ -z "${value}" ]] && SHOW_ERROR "$error_message" && exit 1;
}

ASSERT_SAME_SIZE() {
  local option1=$1
  local -n list1=$2
  local option2=$3
  local -n list2=$4

  if [[ ${#list1} -ne ${#list2} ]]; then
    SHOW_ERROR "The number of values provided with ${option1} and ${option2} must be the same."
    exit 1
  fi
}

ASSERT_EXISTS() {
  local path="$1"
  local error_message="${2:-$path does not exist}"

  [[ ! -e "$path" ]] && SHOW_ERROR "$error_message" && exit 1;
}


ASSERT_COLUMNS_IN_CSV() {
  local csv="$1"
  local -n _required=$2

  IFS=$([[ $csv == *.tsv ]] && echo -e '\t' || echo ',') read -ra header < "$csv"

  # Check if each required column exists in the CSV file
  for col in "${_required[@]}"; do
    if [[ " ${header[*]} " != *" ${col} "* ]]; then
      SHOW_ERROR "Column '${col}' do not exist in $csv"
      exit 1
    fi
  done
}


SUBMIT_JOB() {
  if [ $# -lt 2 ]; then
    echo "Error: Not enough arguments"
    echo "Usage: submit_job scheduler [scheduler_options] command"
    return 1
  fi

  local scheduler="${1,,}"  # Convert to lower case
  local cmd="${*:2}"
#  shift

  echo scheduler="$scheduler"
  case "$scheduler" in
    "local")
      $cmd
      ;;
    "sge"|"pbs")
      # Check if qsub is installed
      if ! command -v qsub &> /dev/null; then
          echo "Error: qsub could not be found"
          return 1
      fi
      qsub "$cmd"
      ;;
    "slurm")
      # Check if sbatch is installed
      if ! command -v sbatch &> /dev/null; then
          echo "Error: sbatch could not be found"
          return 1
      fi
      sbatch "$cmd"
      ;;
    *)
      echo "Error: Unsupported scheduler '$scheduler'"
      return 1
      ;;
  esac
}


PRINT_VERSION() {
  echo -e "\z-brains April 2023 (Version ${VERSION})\n"
}


# Export
export COLOR_ERROR COLOR_WARNING COLOR_NOTE COLOR_INFO COLOR_TITLE NO_COLOR
export -f DO_CMD SHOW_ERROR SHOW_WARNING SHOW_NOTE SHOW_INFO SHOW_TITLE SUBMIT_JOB PRINT_VERSION \
          PARSE_OPTION_SINGLE_VALUE PARSE_OPTION_MULTIPLE_VALUES ASSERT_REQUIRED ASSERT_SAME_SIZE \
          ASSERT_EXISTS ASSERT_COLUMNS_IN_CSV

