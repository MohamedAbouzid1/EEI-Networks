#!/usr/bin/sh

# Subsetting of Drosophila melanogaster GTF file

# select exons only
awk -F '\t' '$3 == "exon"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Gallus-Gallus/data/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Gallus-Gallus/data/processed_exons.tmp"

# select CDS only
awk -F '\t' '$3 == "CDS"' "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Gallus-Gallus/data/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.114.gtf" > "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Gallus-Gallus/data/processed_cds.tmp"
