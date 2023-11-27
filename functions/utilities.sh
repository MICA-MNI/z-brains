#!/bin/bash


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
#        \"Name\": \"z-brains\",
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
# VERBOSE is defined by the z-brains main script

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
  if [[ ${VERBOSE} -gt 0 || ${VERBOSE} -lt 0 ]]; then echo  -e "${COLOR_WARNING}\n[ WARNING ]..... $1 ${NO_COLOR}"; fi
}

SHOW_NOTE() {
  if [[ ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then echo -e "\t\t$1\t${COLOR_NOTE}$2 ${NO_COLOR}"; fi
}

SHOW_INFO() {
  if [[ ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then echo  -e "${COLOR_INFO}\n[ INFO ]..... $1 ${NO_COLOR}"; fi
}

SHOW_TITLE() {
if [[ ${VERBOSE} -gt 1 || ${VERBOSE} -lt 0 ]]; then echo -e "\n${COLOR_TITLE}
-------------------------------------------------------------
\t$1
-------------------------------------------------------------${NO_COLOR}"
fi
}


# Export
export COLOR_ERROR COLOR_WARNING COLOR_NOTE COLOR_INFO COLOR_TITLE NO_COLOR
export -f DO_CMD SHOW_ERROR SHOW_WARNING SHOW_NOTE SHOW_INFO SHOW_TITLE

