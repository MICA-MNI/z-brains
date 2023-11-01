#/bin/bash

# ----------------------------------------------------------------------------------------------- #
# -------------------------------- Software and SGE configuration ------------------------------- #
# ----------------------------------------------------------------------------------------------- #

# Permissions
umask 002

# Save OLD PATH
export OLD_PATH=$PATH

# ----------------------------------------------------------------------------------------------- #
# ------------------------------ Software paths and configuration ------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# User defined PATHS
# AFNI
export AFNIDIR="/data/mica1/01_programs/afni-20.2.06"
# ANTS
export ANTSPATH="/data/mica1/01_programs/ants-2.3.4/bin"
# Workbench
export workbench_path="/data/mica1/01_programs/workbench-1.4.2/bin_linux64"
# FreeSurfer and fastsurfer
export FREESURFER_HOME="/data/mica1/01_programs/freesurfer-7.3.2"
export FASTSURFER_HOME="/data_/mica1/01_programs/fastsurfer"
export fs_licence="/data_/mica1/01_programs/freesurfer-7.3.2/license.txt"
unset TMPDIR

# Remove any other instances of dependencies from the PATH
# ANTS
PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*ants*});IFS=':';echo "${p[*]}";unset IFS)
# Workbench binaries
PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*workbench*});IFS=':';echo "${p[*]}";unset IFS)
# Any other python/conda configuration from the PATH 
PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*conda*});IFS=':';echo "${p[*]}";unset IFS)
# Byebye LD_LIBRARY_PATH variable
LD_LIBRARY_PATH=$(IFS=':';p=($LD_LIBRARY_PATH);unset IFS;p=(${p[@]%%*conda*});IFS=':';echo "${p[*]}";unset IFS)

# FreeSurfer configuration
source "${FREESURFER_HOME}/FreeSurferEnv.sh"
# PYTHON 3.7 configuration
unset PYTHONPATH
unset PYTHONHOME
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
conda3_bin=/data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/bin/
source /data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/etc/profile.d/conda.sh

# Export new PATH with al the necessary binaries
export PATH="${ANTSPATH}:${workbench_path}:${FREESURFER_HOME}/bin/:${FASTSURFER_HOME}:${conda3_bin}:${PATH}"
conda activate /data/mica1/01_programs/micapipe-v0.2.0_conda/micapipe

# Add the number of threads to use here. Note that this is overwritten by
# $NSLOTS if it exists (i.e. when running on SGE).
local_threads="$1"
if [[ -z $local_threads ]]; then export local_threads=10; fi

# Set basic global variables
if [[ ! -z "$NSLOTS" ]]; then
    export threads="$NSLOTS"
else
    export threads="$local_threads"
fi
export OMP_NUM_THREADS="$threads"

# Where processing will run
if [[ -z "$PROC" ]]; then export PROC="LOCAL-MICA"; fi

# Set temporary directory depending on the node
host=$(echo "$HOSTNAME" | awk -F '.' '{print $1}')
case $host in
    fladgate*|yeatman*|oncilla*) export tmpDir="/host/$host/local_raid/temporaryLocalProcessing" ;;
    cassio*|varro*) export tmpDir="/export02/data/temporaryLocalProcessing" ;;
    *) export tmpDir="/data/mica2/temporaryNetworkProcessing" ;;
esac

export SGE_ROOT=/opt/sge
