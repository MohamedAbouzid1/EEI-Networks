#!/usr/bin/sh

# Subsetting of Rattus norvegicus GTF file

# select exons only
awk -F '\t' '$3 == "exon"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Rattus-Norvegicus/data/Rattus_norvegicus.GRCr8.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Rattus-Norvegicus/data/processed_exons.tmp"

# select CDS only
awk -F '\t' '$3 == "CDS"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Rattus-Norvegicus/data/Rattus_norvegicus.GRCr8.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Rattus-Norvegicus/data/processed_cds.tmp"
