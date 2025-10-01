#!/usr/bin/sh

# Subsetting of Homo Sapiens GTF file

# select exons only
awk -F '\t' '$3 == "exon"' ./data/Homo_sapiens.GRCh38.113.gtf > ./data/processed_exons.tmp

#select CDS only
awk -F '\t' '$3 == "CDS"' ./data/Homo_sapiens.GRCh38.113.gtf > ./data/processed_cds.tmp
