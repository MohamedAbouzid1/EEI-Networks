##############################################################################################
# Purpose: Download CIF data and save fasta 
############################################################################################## 

rm(list=ls())
cat("Starting PDB download script...\n")
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

cat("Loading required libraries...\n")
library('seqinr')
library('data.table')
library('stringr')
cat("Libraries loaded successfully\n")

# NOT turning warnings into errors because conversion of warnign during conversion of 3-letter code to AA symbol
options(warn=0)

################################################################
# all filtered uniprot pdb map
cat("\nReading processed complex file...\n")
uniprot_all <- fread('../data/processed_complex_final.txt', sep='\t', header=TRUE)
keep <- unique(tolower(uniprot_all$ID))
cat(sprintf("Found %d unique PDB structures to process\n", length(keep)))

## Define directories
public_cif_dir <- '/cosybio/project/public_data/PDB_CIF'
base_dir <- '../data'
local_cif_dir <- file.path(base_dir, 'PDB_CIF')
download_list <- file.path(base_dir, 'all_new_CIFs.txt')

# Create local directory if it doesn't exist
if(!dir.exists(base_dir)) {
    cat(sprintf("Creating base directory: %s\n", base_dir))
    dir.create(base_dir, recursive=TRUE)
}

if(!dir.exists(local_cif_dir)) {
    cat(sprintf("Creating local CIF directory: %s\n", local_cif_dir))
    dir.create(local_cif_dir, recursive=TRUE)
}

# Check public directory for existing files
if(dir.exists(public_cif_dir)) {
    cat(sprintf("\nChecking public CIF directory: %s\n", public_cif_dir))
    public_files <- list.files(public_cif_dir, pattern="\\.cif$")
    public_pdbs <- unique(unlist(lapply(strsplit(public_files, '[.]'), '[[', 1)))
    cat(sprintf("Found %d PDB structures in public directory\n", length(public_pdbs)))
    
    # Find which PDBs we need that exist in public directory
    pdbs_to_copy <- intersect(keep, public_pdbs)
    cat(sprintf("Found %d required PDB structures in public directory\n", length(pdbs_to_copy)))
    
    # Copy existing files to local directory
    if(length(pdbs_to_copy) > 0) {
        cat("\nCopying existing PDB files to local directory...\n")
        for(pdb in pdbs_to_copy) {
            src_file <- file.path(public_cif_dir, paste0(pdb, '.cif'))
            dst_file <- file.path(local_cif_dir, paste0(pdb, '.cif'))
            if(!file.exists(dst_file)) {
                cat(sprintf("  Copying %s...\n", pdb))
                file.copy(src_file, dst_file)
            }
        }
        cat("Copying complete\n")
    }
    
    # Find which PDBs we still need to download
    pdbs_to_download <- setdiff(keep, public_pdbs)
    cat(sprintf("\nNeed to download %d new PDB structures\n", length(pdbs_to_download)))
    
    if(length(pdbs_to_download) > 0) {
        cat("PDBs to download:\n")
        cat(paste(head(pdbs_to_download, 5), collapse=", "))
        if(length(pdbs_to_download) > 5) cat(sprintf(" ... and %d more\n", length(pdbs_to_download) - 5))
        
        # Download missing PDBs
        cat("\nPreparing to download missing PDBs...\n")
        allpdbs <- paste(pdbs_to_download, collapse=',')
        cat(sprintf("Saving list of PDBs to download in: %s\n", download_list))
        writeLines(allpdbs, download_list)
        
        cat("\nStarting batch download...\n")
        download_cmd <- paste0('./batch_download.sh -f ', download_list, ' -o ', local_cif_dir, ' -c')
        cat(sprintf("Running command: %s\n", download_cmd))
        system(download_cmd)
        
        cat("\nDecompressing downloaded files...\n")
        decompress_cmd <- paste0('gunzip ', local_cif_dir, '/*.cif.gz')
        cat(sprintf("Running command: %s\n", decompress_cmd))
        system(decompress_cmd)
    }
} else {
    cat(sprintf("\nPublic CIF directory not found: %s\n", public_cif_dir))
    cat("Will download all required PDB structures to local directory\n")
    allpdbs <- paste(keep, collapse=',')
    writeLines(allpdbs, download_list)
    
    cat("\nStarting batch download...\n")
    download_cmd <- paste0('./batch_download.sh -f ', download_list, ' -o ', local_cif_dir, ' -c')
    cat(sprintf("Running command: %s\n", download_cmd))
    system(download_cmd)
    
    cat("\nDecompressing downloaded files...\n")
    decompress_cmd <- paste0('gunzip ', local_cif_dir, '/*.cif.gz')
    cat(sprintf("Running command: %s\n", decompress_cmd))
    system(decompress_cmd)
}

# Verify final results
final_files <- list.files(local_cif_dir, pattern="\\.cif$")
cat(sprintf("\nFinal status: Found %d CIF files in local directory\n", length(final_files)))
missing_pdbs <- setdiff(keep, gsub("\\.cif$", "", final_files))
if(length(missing_pdbs) > 0) {
    cat(sprintf("Warning: %d PDB structures are still missing:\n", length(missing_pdbs)))
    cat(paste(head(missing_pdbs, 5), collapse=", "))
    if(length(missing_pdbs) > 5) cat(sprintf(" ... and %d more\n", length(missing_pdbs) - 5))
} else {
    cat("All required PDB structures are present in local directory\n")
}

cat("\nScript execution completed.\n")
