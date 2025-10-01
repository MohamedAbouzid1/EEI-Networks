#!/usr/bin/sh

# Subsetting of Drosophila melanogaster GTF file

# select exons only
awk -F '\t' '$3 == "exon"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Drosophila/data/Drosophila_melanogaster.BDGP6.54.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Drosophila/data/processed_exons.tmp"

# select CDS only
awk -F '\t' '$3 == "CDS"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Drosophila/data/Drosophila_melanogaster.BDGP6.54.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Drosophila/data/processed_cds.tmp"
