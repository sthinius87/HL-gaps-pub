#!/bin/bash

# Script: Cleanup Temporary Files and Directories
# Description: This script removes temporary files and directories created by previous runs.
#              It ensures a clean working environment for subsequent executions.

# Remove the grep output file, if it exists.
#rm -f grep.out

# Remove the temporary directory and its contents, if it exists.
rm -rf tmp/

# Recreate the temporary directory.
mkdir -p tmp/

# Remove any raw result files, if they exist.
#rm -f results*.raw

# Remove the job store directory and its contents, if it exists.
rm -rf job_store/

# End of script.