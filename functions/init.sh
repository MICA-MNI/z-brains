#!/bin/bash

# MICA-specific environment variables required fot zbrains - adjust as needed for your configuration
export ANTSPATH="/data/mica1/01_programs/ants-2.3.4/bin"
# export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1

export WORKBENCH_PATH="/data/mica1/01_programs/workbench-1.4.2/bin_linux64"
# export OMP_NUM_THREADS=1

CONDA_PATH="/data/mica1/01_programs/miniforge3"
source $CONDA_PATH/bin/activate
conda activate zbrains
