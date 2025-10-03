#!/usr/bin/sh

# Subsetting of Saccharomyces cerevisiae GTF file

# select exons only
awk -F '\t' '$3 == "exon"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Saccharomyces/Saccharomyces/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Saccharomyces/Saccharomyces/Saccharomyces_cerevisiae/processed_exons.tmp"

# select CDS only
awk -F '\t' '$3 == "CDS"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Saccharomyces/Saccharomyces/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Saccharomyces/Saccharomyces/Saccharomyces_cerevisiae/processed_cds.tmp"
