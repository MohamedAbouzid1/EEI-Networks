#!/bin/bash
# Wrapper script to set environment variables before running the Python script

# Source bashrc to set up the shell environment first
source ~/.bashrc

# Activate conda environment if available
if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate selenium_env
fi

# Set library path to conda environment
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# Set temporary directory
export TMPDIR=$HOME/tmp
mkdir -p $TMPDIR

# Force enable Wayland support
export MOZ_ENABLE_WAYLAND=1

# Run the Python script with passed arguments
python "$@"
