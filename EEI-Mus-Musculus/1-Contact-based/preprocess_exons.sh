#!/usr/bin/sh

# Subsetting of Mus Musculus GTF file

# select exons only
awk -F '\t' '$3 == "exon"' ./Mus_musculus/Mus_musculus.GRCm39.113.gtf > ./Mus_musculus/processed_exons.tmp

# select CDS only
awk -F '\t' '$3 == "CDS"' ./Mus_musculus/Mus_musculus.GRCm39.113.gtf > ./Mus_musculus/processed_cds.tmp
