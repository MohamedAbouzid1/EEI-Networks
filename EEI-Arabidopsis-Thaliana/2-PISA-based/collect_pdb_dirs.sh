#!/bin/bash

# Path to the PISA directory
PISA_DIR="/cosybio/project/EEIcancer/data/PISA"

# Output file name
OUTPUT_FILE="/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/data/pdb_directory_list.txt"

# List all directories in the PISA folder and save names to file
ls -l "$PISA_DIR" | grep "^d" | awk '{print $NF}' > "$OUTPUT_FILE"

# Print summary
echo "Found $(wc -l < "$OUTPUT_FILE") PDB directories"
echo "List saved to $OUTPUT_FILE"