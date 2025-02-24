#!/bin/bash

# Base name for the compressed archives
BASE_NAME="db_split_part"
TARGET_DIR="db_split"  # Extract all files to this directory

# Create the target directory if it doesn't exist
if [ ! -d "$TARGET_DIR" ]; then
  mkdir -p "$TARGET_DIR"
fi

# Loop through the archives
for i in {1..5}; do
  TAR_FILE="$BASE_NAME$i.tar.xz"

  # Check if the archive exists
  if [ ! -f "$TAR_FILE" ]; then
    echo "Error: Archive '$TAR_FILE' not found."
    continue  # Skip to the next iteration
  fi

  echo "Extracting '$TAR_FILE' to '$TARGET_DIR'..."

  # Use tar to extract the archive.
  tar -xJvf "$TAR_FILE" -C "$TARGET_DIR" \
    --strip-components=1

  if [ $? -ne 0 ]; then
    echo "Error: Extraction failed for '$TAR_FILE'."
  else
    echo "Extraction complete for '$TAR_FILE'."
  fi
done

echo "All done."
