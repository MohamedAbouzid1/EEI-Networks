#!/bin/bash

# Example script for running the EEI prediction pipeline with Drosophila data
# This demonstrates how to use the updated species-agnostic pipeline

# Set environment variables for Drosophila-specific files
EGIO_FILE="../EGIO/ExonGroup_testpro_hsa_dme.txt"
ORTHOLOG_EEI_FILE="../data/PISA_networks_filtered/PISA_EEIN_0.5.txt"
HUMAN_EXON_MAP="../../orthology_based_EEI_prediction/results_CONTACT_BASED/exon_coord_mapping_human.tsv"
ORTHOLOG_EXON_MAP="../exon_coord_map/exon_coord_map_dse.tsv"
HUMAN_EEI_FILE="../../EEI-Homo-Sapiens-2024/results/2-PISA_results/PISA_EEIN_0.5_human.txt"
OUTPUT_FILE="../data/predicted/predicted_pisa_human_eeis.tsv"

# Set parameters
IDENTITY_THRESHOLD=0.0

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

echo "Running EEI prediction pipeline with Drosophila data..."
echo "EGIO file: $EGIO_FILE"
echo "Drosophila EEI file: $ORTHOLOG_EEI_FILE"
echo "Human exon map: $HUMAN_EXON_MAP"
echo "Drosophila exon map: $ORTHOLOG_EXON_MAP"
echo "Human EEI file: $HUMAN_EEI_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Identity threshold: $IDENTITY_THRESHOLD"
echo ""

# Run the EEI prediction pipeline
python main.py \
  --egio "$EGIO_FILE" \
  --eei "$ORTHOLOG_EEI_FILE" \
  --human_exon_map "$HUMAN_EXON_MAP" \
  --ortholog_exon_map "$ORTHOLOG_EXON_MAP" \
  --human_eei "$HUMAN_EEI_FILE" \
  --identity_threshold "$IDENTITY_THRESHOLD" \
  --output "$OUTPUT_FILE"

# Check if the pipeline completed successfully
if [ $? -eq 0 ]; then
  echo ""
  echo "Pipeline completed successfully!"
  echo "Results saved to: $OUTPUT_FILE"
  echo ""
  echo "You can now use this same pipeline with other species by updating the file paths:"
  echo "  - Change --egio to point to your species-specific EGIO file"
  echo "  - Change --eei to point to your species-specific EEI file"
  echo "  - Change --ortholog_exon_map to point to your species-specific exon mapping file"
  echo ""
  echo "Example for a different species:"
  echo "  --egio ../EGIO/ExonGroup_testpro_hsa_[species_code].txt"
  echo "  --eei ../data/[species]_EEI_network.txt"
  echo "  --ortholog_exon_map ../exon_coord_map/exon_coord_map_[species].tsv"
else
  echo ""
  echo "Pipeline failed. Check the error messages above."
  exit 1
fi
