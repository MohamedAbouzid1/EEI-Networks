#!/bin/bash
# Script to set up environment before running R script in tmux

# Define the command to run
CMD="$1"

# Source bashrc to set up shell environment first
source ~/.bashrc

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate selenium_env

# Set library path to conda environment
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# Set temporary directory
export TMPDIR=$HOME/tmp
mkdir -p $TMPDIR

# Display environment information for debugging
echo "Running with:"
echo "CONDA_PREFIX=$CONDA_PREFIX"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "TMPDIR=$TMPDIR"

# Execute the command
eval "$CMD"
