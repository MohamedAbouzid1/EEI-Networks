#!/usr/bin/sh

# Subsetting of Bos_taurus GTF file

# select exons only
awk -F '\t' '$3 == "exon"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Bos-Taurus/data/Bos_taurus.ARS-UCD2.0.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Bos-Taurus/data/processed_exons.tmp"

# select CDS only
awk -F '\t' '$3 == "CDS"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Bos-Taurus/data/Bos_taurus.ARS-UCD2.0.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Bos-Taurus/data/processed_cds.tmp"
