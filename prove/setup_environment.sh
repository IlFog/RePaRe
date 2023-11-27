#!/bin/bash

# Check if conda is installed and in the PATH
if ! command -v conda &> /dev/null
then
    echo "conda not found. Please install Miniconda or Anaconda."
    exit
fi

# Create the conda environment
conda env create -f snakemake_env.yml

echo "Environment setup complete. Activate with 'conda activate Repare'"