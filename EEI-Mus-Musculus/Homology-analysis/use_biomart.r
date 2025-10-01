###############################################################################
# Purpose: use_biomart.r script to get, from the UniProtID, all the Ensembl ID that are used by EGIO
###############################################################################

rm(list=ls())

# Install biomaRt if not already installed

if (!requireNamespace("biomaRt", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("biomaRt")
}

library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)

species <- args[1]
base_dir <- args[2]
ensembl_dataset <- args[3]


# Use the Ensembl database
ensembl <- useMart("ensembl", dataset=ensembl_dataset) 

# Get Ensembl ID for a given UniProt ID
uniprot_file <- file.path(base_dir, paste0(species,"_uniprotid.txt"))
list <- readLines(uniprot_file)


result <- getBM(attributes=c("uniprotswissprot", "ensembl_gene_id", "ensembl_transcript_id"),
                filters="uniprotswissprot",
                values=list,
                mart=ensembl)

#print(result)

output_file <- file.path(base_dir, paste0(species,"_ensemblid.txt"))


write.table(result, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
