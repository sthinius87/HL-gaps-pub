#!/bin/bash

# Script: Run CWL Workflow with Toil (Normal or Restart Mode)
# Description: This script executes a CWL workflow using Toil with Slurm, allowing the user to choose between normal and restart modes.

# Set Slurm arguments for Toil.
export TOIL_SLURM_ARGS="-t 100-00:00:00 --export=ALL,PATH"

# NORMAL
toil-cwl-runner --workDir ./tmp \
--jobStore ./job_store \
--debug \
--batchSystem slurm \
--caching FALSE \
--defaultDisk 1G \
--defaultMemory 1G \
--disableProgress TRUE \
--clean never \
--cleanWorkDir always \
--maxNodes 12 \
--maxJobs 500 \
--defaultCores 4 \
--stats \
scatter_wf.cwl \
scatter_inp.yml \

# RESTART
toil-cwl-runner --workDir ./tmp \
--jobStore ./job_store \
--debug \
--batchSystem slurm \
--caching FALSE \
--defaultDisk 1G \
--defaultMemory 1G \
--disableProgress TRUE \
--clean never \
--cleanWorkDir always \
--maxNodes 12 \
--maxJobs 500 \
--defaultCores 4 \
--stats \
--restart \
scatter_wf.cwl \
scatter_inp.yml \

