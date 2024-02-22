#!/bin/bash


DO_CMD() {
  # do_cmd sends command to stdout before executing it.
  if [ $# -eq 1 ]; then
    readarray -t array < <(echo "$1" | xargs -n1)
  else
    array=("$@")
  fi

  local str_cmd=""
  for element in "${array[@]}"; do
      [[ $element =~ [[:space:]] ]] && element="\"$element\""
      str_cmd+="$element "
  done
  str_cmd="${str_cmd% }"

#  if [[ -z ${VERBOSE} || ${VERBOSE} -gt 2 || ${VERBOSE} -lt 0 ]]; then
#    str="$(whoami) @ $(uname -n) $(date)"
#    echo -e "\033[38;5;118m${str}:\nCOMMAND -->  \033[38;5;122m${cmd}\n\033[0m";
#  fi

  local header
  header="$(whoami) @ $(uname -n) $(date)"
  log_message 3 "\033[38;5;118m${header}:\nCOMMAND -->  \033[38;5;122m${str_cmd}\n\033[0m"

  "${array[@]}"
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

#SHOW_ERROR() {
##echo -e "${COLOR_ERROR}\n-------------------------------------------------------------\n\n[ ERROR ]..... $1\n
##-------------------------------------------------------------${NO_COLOR}\n"
#  echo ""
#  echo -e "${COLOR_ERROR}[ ERROR ]     $1${NO_COLOR}"
#  for s in "${@:2}"; do echo -e "${COLOR_ERROR}              $s${NO_COLOR}"; done
#  echo ""
#}
#
#SHOW_WARNING() {
#  if [[ -z ${VERBOSE} || ${VERBOSE} -gt 0 || ${VERBOSE} -lt 0 ]]; then
#    echo ""
#    echo -e "${COLOR_WARNING}[ WARNING ]   $1${NO_COLOR}"
#    for s in "${@:2}"; do echo -e "${COLOR_WARNING}              $s${NO_COLOR}"; done
#  fi
#}
#
#SHOW_NOTE() {
#  if [[ -z ${VERBOSE} || ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then
#    echo -e "              $1\t${COLOR_NOTE}$2 ${NO_COLOR}";
#  fi
#}
#
#SHOW_INFO() {
#  if [[ -z ${VERBOSE} || ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then
#    echo ""
#    echo -e "${COLOR_INFO}[ INFO ]      $1${NO_COLOR}"
#    for s in "${@:2}"; do echo -e "${COLOR_INFO}              $s${NO_COLOR}"; done
#  fi
#}
#
#SHOW_TITLE() {
#  if [[ -z ${VERBOSE} || ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then
#    echo -e "\n${COLOR_TITLE}-------------------------------------------------------------${NO_COLOR}"
#    echo -e "${COLOR_TITLE}$1${NO_COLOR}"
#    for s in "${@:2}"; do echo -e "${COLOR_TITLE}$s${NO_COLOR}"; done
#    echo -e "${COLOR_TITLE}-------------------------------------------------------------${NO_COLOR}"
#    echo ""
#  fi
#}

log_message() {
    local level=$1
    local messages="${*:2}"

#    echo -e "${messages[*]}" >> "${LOGFILE}"
    if [[ -n $LOGFILE ]]; then
      echo -e "${messages[*]}" | sed 's/\x1b\[[0-9;]*m//g' >> "${LOGFILE}"
    fi

    if [[ -z $VERBOSE || $VERBOSE -ge $level || $VERBOSE -lt 0 ]]; then
      echo -e "${messages[*]}"
    fi
}

SHOW_ERROR() {
  str="\n${COLOR_ERROR}[ ERROR ]     $1${NO_COLOR}\n"
  for s in "${@:2}"; do str+="${COLOR_ERROR}              $s${NO_COLOR}\n"; done

  log_message 0 "$str"
}

SHOW_WARNING() {
  str="${COLOR_WARNING}[ WARNING ]   $1${NO_COLOR}"
  for s in "${@:2}"; do str+="\n${COLOR_WARNING}              $s${NO_COLOR}"; done

  log_message 1 "$str"
}

SHOW_NOTE() {
  str="              $1\t${COLOR_NOTE}$2 ${NO_COLOR}"

  log_message 2 "$str"
}

SHOW_INFO() {
  str="${COLOR_INFO}[ INFO ]      $1${NO_COLOR}"
  for s in "${@:2}"; do str+="\n${COLOR_INFO}              $s${NO_COLOR}"; done

  log_message 2 "$str"
}

SHOW_TITLE() {
  str="\n${COLOR_TITLE}-------------------------------------------------------------${NO_COLOR}\n"
  str+="${COLOR_TITLE}$1${NO_COLOR}\n"
  for s in "${@:2}"; do str+="${COLOR_TITLE}$s${NO_COLOR}\n"; done
  str+="${COLOR_TITLE}-------------------------------------------------------------${NO_COLOR}\n"

  log_message 2 "$str"
}


function allowed_to_regex() {
  local array=("$@")
#  array=("${array[@]/#/ }")  # Adds a leading space
#  array=("${array[@]/%/ }")  # Adds a trailing space
  IFS='|'; echo "${array[*]}"; unset IFS
}

PARSE_OPTION_SINGLE_VALUE() {
  # PARSE_OPTION_SINGLE_VALUE output_variable args allowed_values
  # or
  # PARSE_OPTION_SINGLE_VALUE output_variable args
  local -n _out=$1
  local -n _args=$2

  local allowed_regex
  if [[ $# -gt 2 ]]; then
    local -n _allowed_values=$3
    allowed_regex=$(allowed_to_regex "${_allowed_values[@]}")
  fi

  local option="${_args[0]}"
  local value="${_args[1]}"

  # Check if the option has a value and if the value is not another option
  if [[ -z "${value}" || "${value}" == --* ]]; then
      SHOW_ERROR "${option} option requires a value."
      return 1
  fi

  # Check if value is in the list of allowed values
  if [[ -v allowed_regex && ! " ${value} " =~ ${allowed_regex} ]]; then
      SHOW_ERROR "Invalid value '${value}' for ${option} option." "Allowed values are: ${allowed_regex//|/, }."
      return 1
  fi

  _out="${value}"
  _args=("${_args[@]:2}")  # shift
  return 0
}


PARSE_OPTION_MULTIPLE_VALUES() {
  # PARSE_OPTION_MULTIPLE_VALUES output_variable args allowed_values all
  # or
  # PARSE_OPTION_MULTIPLE_VALUES output_variable args allowed_values
  # or
  # PARSE_OPTION_MULTIPLE_VALUES output_variable args
  local -n _out=$1
  # shellcheck disable=SC2178
  local -n _args=$2
  [[ $# -gt 2 ]] && local -n _allowed_values=$3
  [[ $# -gt 3 ]] && all=$4;  # accept "all" as an additional value

  local allowed_regex
  if [[ -v "${_allowed_values[*]}" ]]; then
    local list_allowed=("${_allowed_values[@]}")
    [[ -n "$all" ]] && list_allowed+=("$all");
    allowed_regex=$(allowed_to_regex "${list_allowed[@]}")
  fi

  local option="${_args[0]}"
  _args=("${_args[@]:1}")

  local _values=()
  while [[ ${#_args[@]} -gt 0 && "${_args[0]}" != --* ]]; do
    local value="${_args[0]}"

    # Check if the value is in the list of allowed values
    if [[ -v allowed_regex && ! " ${value} " =~ ${allowed_regex} ]]; then
      SHOW_ERROR "Invalid value '$value' for ${option} option." "Allowed values are: ${allowed_regex//|/, }."
      return 1
    fi

    _values+=("$value")
    _args=("${_args[@]:1}")  # shift
  done

  [[ ${#_values[@]} == 0 ]] && SHOW_ERROR "${option} option requires at least one value." && return 1

#  if [[ -v all || " ${_values[*]} " != *" $all "* ]]; then
#    _values=("${_allowed_values[@]}")
#  fi

  _out=("${_values[@]}")
  return 0
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

  if [[ ${#list1[@]} -ne ${#list2[@]} ]]; then
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

