#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --output=DemultiplexOnly-%N-%j.output
#SBATCH --error=DemultiplexOnly-%N-%j.error
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1

# Description: Placeholder script for samples which require demultiplexing
# Author:      AWMGS
# Mode:        BY_SAMPLE
# Use:         sbatch within sample directory

cd "$SLURM_SUBMIT_DIR"

















