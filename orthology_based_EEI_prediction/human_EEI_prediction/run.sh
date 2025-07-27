#!/bin/bash

# Orthology-based EEI Prediction Pipeline Runner
# This script runs the EEI prediction pipeline with the specified parameters

# Set environment variables - update these paths according to your file locations
EGIO_FILE="/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EGIO/results/ExonGroup_testpro_hsa_mus.txt"
MOUSE_EEI_FILE="/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Mus-Musculus/3-EPPIC-based/results/EPPIC_EEIN_filtered.txt"
HUMAN_EXON_MAP="/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/orthology_based_EEI_prediction/results_CONTACT_BASED/exon_coord_mapping_human.tsv"
MOUSE_EXON_MAP="/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/orthology_based_EEI_prediction/results_CONTACT_BASED/exon_coord_mapping_mouse.tsv"
HUMAN_EEI_FILE="/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Homo-Sapiens-2024/results/3-EPPIC_results/EPPIC_EEIN_filtered.txt" 
OUTPUT_FILE="/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/orthology_based_EEI_prediction/results_EPPIC_Based/predicted_human_eeis_fixed_no_conf_threshold_all_iden_0.tsv"

# Set parameters
IDENTITY_THRESHOLD=0.0

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Run the EEI prediction pipeline
python main.py \
  --egio "$EGIO_FILE" \
  --eei "$MOUSE_EEI_FILE" \
  --human_exon_map "$HUMAN_EXON_MAP" \
  --mouse_exon_map "$MOUSE_EXON_MAP" \
  --human_eei "$HUMAN_EEI_FILE" \
  --identity_threshold "$IDENTITY_THRESHOLD" \
  --output "$OUTPUT_FILE"

# Check if the pipeline completed successfully
if [ $? -eq 0 ]; then
  echo "Pipeline completed successfully. Results saved to $OUTPUT_FILE"
else
  echo "Pipeline failed. Check the error messages above."
  exit 1
fi