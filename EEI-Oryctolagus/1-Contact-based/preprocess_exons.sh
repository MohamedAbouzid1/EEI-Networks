#!/usr/bin/sh

# Subsetting of Oryctolagus cuniculus GTF file

# select exons only
awk -F '\t' '$3 == "exon"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Oryctolagus/data/Oryctolagus_cuniculus.OryCun2.0.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Oryctolagus/data/processed_exons.tmp"

# select CDS only
awk -F '\t' '$3 == "CDS"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Oryctolagus/data/Oryctolagus_cuniculus.OryCun2.0.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Oryctolagus/data/processed_cds.tmp"
