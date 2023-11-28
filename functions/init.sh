#!/bin/bash

# MICA-specific environment variables required fot zbrains - adjust as needed for your configuration
export ANTSPATH="/data/mica1/01_programs/ants-2.3.4/bin"
export WORKBENCH_PATH="/data/mica1/01_programs/workbench-1.4.2/bin_linux64"

CONDA_PATH="/data/mica1/01_programs/miniforge3"
source $CONDA_PATH/bin/activate
conda activate zbrains
