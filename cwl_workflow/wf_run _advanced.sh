#!/bin/bash

# Script: Run CWL Workflow with Toil (Normal or Restart Mode)
# Description: This script executes a CWL workflow using Toil with Slurm, allowing the user to choose between normal and restart modes.

# Set Slurm arguments for Toil.
export TOIL_SLURM_ARGS="-t 100-00:00:00 --export=ALL,PATH"

# Define common Toil arguments.
TOIL_COMMON_ARGS=(
    --workDir ./tmp
    --jobStore ./job_store
    --debug
    --batchSystem slurm
    --caching FALSE
    --defaultDisk 1G
    --defaultMemory 1G
    --disableProgress TRUE
    --clean never
    --cleanWorkDir always
    --maxNodes 12
    --maxJobs 500
    --defaultCores 4
    --stats
)

# Define workflow and input files.
WORKFLOW="scatter_wf.cwl"
INPUT_YAML="scatter_inp.yml"

# Prompt the user to choose between normal and restart modes.
read -p "Choose mode (normal/restart): " MODE

# Execute Toil based on the chosen mode.
case "$MODE" in
    normal)
        echo "Running in normal mode..."
        toil-cwl-runner "${TOIL_COMMON_ARGS[@]}" "$WORKFLOW" "$INPUT_YAML"
        ;;
    restart)
        echo "Running in restart mode..."
        toil-cwl-runner "${TOIL_COMMON_ARGS[@]}" --restart "$WORKFLOW" "$INPUT_YAML"
        ;;
    *)
        echo "Invalid mode. Please choose 'normal' or 'restart'."
        exit 1
        ;;
esac

# End of script.