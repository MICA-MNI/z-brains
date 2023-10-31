#!/bin/bash

export Version="v0.0.2 'Reborn'"

bids_variables() {
  # This functions assigns variable names acording to:
  #     Subject ID =                        $1
  #     Raw data directory =                $2
  #     Hippunfold directory derivatives =  $3
  #     Micapipe derivatives directory =    $4
  #     Z-brains derivatives directory =    $5
  #     Session ID (optional) =             $6
  id=$1
  rawdir=$2
  hippdir=$3
  micapipedir=$4
  outdir=$5
  SES=$6
  umask 002

  # Define utilities/functions directory
  export scriptDir=${ZBRAINS}/functions

  export subject=sub-${id}
  
  # Handle Single Session
  if [ "$SES" == "SINGLE" ]; then
      export subject_raw=${rawdir}/${subject}           # subject's raw data directory
      export subject_hipp=${hippdir}/hippundold/${subject}         # subject's hippunfold directory
      export subject_micapipe=${micapipedir}/${subject} # subject's micapipe directory
      export subject_out=${outdir}/${subject}           # subject's output directory
      ses=""
  else
      export subject_raw=${rawdir}/${subject}/${SES}           # subject's raw data directory
      export subject_hipp=${hippdir}/hippunfold/${subject}/${SES}         # subject's hippunfold directory
      export subject_micapipe=${micapipedir}/${subject}/${SES} # subject's micapipe directory
      export subject_out=${outdir}/${subject}/${SES}           # subject's output directory
      ses="_${SES}"
  fi

  export idBIDS="${subject}${ses}"

  # Structural directories derivatives/
  export proc_struct=$subject_micapipe/anat     # Structural processing directory
  	 export dir_volum=$subject_micapipe/parc    # Cortical segmentations
  	 export dir_conte69=$subject_micapipe/surf  # Conte69 surfaces
     export dir_maps=$subject_micapipe/maps

  export dir_warp=$subject_micapipe/xfm     # Transformation matrices 
  export dir_logs=$subject_out/logs         # directory where log files will be written
  export dir_QC=$subject_micapipe/QC        # directory with micapipe qc files

  # post structural Files (the resolution might vary depending on the dataset)
  export res=$(mrinfo "${proc_struct}"/"${idBIDS}"_space-nativepro_T1w.nii.gz -spacing | awk '{printf "%.1f\n", $2}')
  export T1nativepro=${proc_struct}/${idBIDS}_space-nativepro_T1w.nii.gz
  export T1nativepro_brain=${proc_struct}/${idBIDS}_space-nativepro_T1w_brain.nii.gz
  export T1nativepro_mask=${proc_struct}/${idBIDS}_space-nativepro_T1w_brain_mask.nii.gz
  export T1fast_seg=${subject_micapipe}/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz

  # raw files
  raw_flair=$(ls "$subject_raw"/anat/*FLAIR*.nii* 2>/dev/null)
}

set_surface_directory() {
  local recon=${1}
  export dir_surf=${outdir/\/z-brains/}/${recon}    # surf
  export dir_subjsurf=${dir_surf}/${idBIDS}  # Subject surface dir
  export T1surf=${dir_subjsurf}/mri/orig.mgz

  # Native midsurface in gifti format
  export lh_midsurf="${dir_conte69}/${idBIDS}_hemi-L_surf-fsnative_label-midthickness.surf.gii"
  export rh_midsurf="${dir_conte69}/${idBIDS}_hemi-R_surf-fsnative_label-midthickness.surf.gii"
}

bids_print.variables-ctx() {
  # This functions prints BIDS variables names and files if found
  Info "Cortical processing variables:"
  Note "T1 nativepro         :" "$(find "$T1nativepro" 2>/dev/null)"
  Note "T1 resolution        :" "$res"
  Note "fs-LR-32k-lh surface :" "$(find "${dir_conte69}/${idBIDS}_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii" 2>/dev/null)"
  Note "fs-LR-32k-rh surface :" "$(find "${dir_conte69}/${idBIDS}_hemi-R_space-nativepro_surf-fsLR-32k_label-white.surf.gii" 2>/dev/null)"

}

bids_print.variables-sctx() {
  # This functions prints BIDS variables names and files if found
  Info "Subcortical processing variables:"
  Note "T1 nativepro    :" "$(find "$T1nativepro" 2>/dev/null)"
  Note "T1 sctx seg     :" "$(find "$T1fast_seg" 2>/dev/null)"
  Note "Sctx volume     :" "$(find "$dir_subjsurf/stats/aseg.stats" 2>/dev/null)"
  Note "T1 resolution   :" "$res"
}

bids_print.variables-hipp() {
  # This functions prints BIDS variables names and files if found
  Info "Hippocampal processing variables:"
  Note "T1 nativepro    :" "$(find "$T1nativepro" 2>/dev/null)"
  Note "T1 resolution   :" "$res"
}

bids_variables_unset() {
  # This function unsets all the enviromentalk variables defined by
  # bids_variables
  unset scriptDir
  unset subject
  unset subject_raw
  unset subject_hipp
  unset subject_micapipe
  unset subject_out
  unset proc_struct
  unset dir_volum
  unset dir_conte69
  unset dir_maps
  unset dir_warp
  unset dir_logs
  unset dir_QC
  unset res
  unset T1nativepro
  unset T1nativepro_brain
  unset T1nativepro_mask
  unset T1fast_seg
  unset raw_flair
  unset dir_surf
  unset dir_subjsurf
  unset T1surf
  unset lh_midsurf
  unset rh_midsurf
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
        \"Version\": \"${Version}\",
        \"Reference\": \"tbd\",
        \"DOI\": \"tbd\",
        \"URL\": \"https://z-brains.readthedocs.io\",
        \"GitHub\": \"https://github.com/MICA-MNI/z-brains\",
        \"Container\": {
          \"Type\": \"tbd\",
          \"Tag\": \"ello:zbrains:$(echo ${Version} | awk '{print $1}')\"
        },
        \"RunBy\": \"$(whoami)\",
        \"Workstation\": \"$(uname -n)\",
        \"LastRun\": \"$(date)\",
        \"Processing\": \"${PROC}\"
      }]
  }" > "${outdir}/dataset_description.json"
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
  tmp=$1
  nocleanup=$2
  here=$3
  # Clean temporary directory and temporary fsaverage5
  if [[ $nocleanup == "FALSE" ]]; then
      rm -Rf "$tmp" 2>/dev/null
  else
      echo -e "z-brains tmp directory was not erased: \n\t\t${tmp}";
  fi
  cd "$here"
  bids_variables_unset
  if [[ ! -z "$OLD_PATH" ]]; then  export PATH=$OLD_PATH; unset OLD_PATH; fi
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
Error() {
echo -e "\033[38;5;9m\n-------------------------------------------------------------\n\n[ ERROR ]..... $1\n
-------------------------------------------------------------\033[0m\n"
}
Note(){
# I replaced color \033[38;5;197m to \033[38;5;122m
if [[ ${quiet} != TRUE ]]; then echo -e "\t\t$1\t\033[38;5;122m$2\033[0m"; fi
}
Info() {
Col="38;5;75m" # Color code
if [[ ${quiet} != TRUE ]]; then echo  -e "\033[$Col\n[ INFO ]..... $1 \033[0m"; fi
}
Warning() {
Col="38;5;184m" # Color code
if [[ ${quiet} != TRUE ]]; then echo  -e "\033[$Col\n[ WARNING ]..... $1 \033[0m"; fi
}
Warn() {
Col="38;5;184m" # Color code
if [[ ${quiet} != TRUE ]]; then echo  -e "\033[$Col
-------------------------------------------------------------\n
[ WARNING ]..... $1
\n-------------------------------------------------------------\033[0m"; fi
}
Title() {
if [[ ${quiet} != TRUE ]]; then echo -e "\n\033[38;5;141m
-------------------------------------------------------------
\t$1
-------------------------------------------------------------\033[0m"; fi
}

function Do_cmd() {
# do_cmd sends command to stdout before executing it.
str="$(whoami) @ $(uname -n) $(date)"
local l_command=""
local l_sep=" "
local l_index=1
while [ ${l_index} -le $# ]; do
    eval arg=\${$l_index}
    if [ "$arg" = "-fake" ]; then
      arg=""
    fi
    if [ "$arg" = "-no_stderr" ]; then
      arg=""
    fi
    if [ "$arg" == "-log" ]; then
      nextarg=$(("${l_index}" + 1))
      eval logfile=\${"${nextarg}"}
      arg=""
      l_index=$[${l_index}+1]
    fi
    l_command="${l_command}${l_sep}${arg}"
    l_sep=" "
    l_index=$[${l_index}+1]
   done
if [[ ${quiet} != TRUE ]]; then echo -e "\033[38;5;118m\n${str}:\nCOMMAND -->  \033[38;5;122m${l_command}  \033[0m"; fi
if [ -z "$TEST" ]; then $l_command; fi
}

