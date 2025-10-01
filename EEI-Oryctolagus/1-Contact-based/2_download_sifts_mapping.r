##############################################################################################
# Purpose: download all xml files of the cross mapping from SIFTS
##############################################################################################
## WARNING: If some errors in the download occur (some files cannot be downloaded, some are downloaded but not unzipped), please
## use the bash script 'sifts_download.sh' in this same folder. It's able to bypassing errors due to time requests.
## One protein complex resulted missing from the database: 7qe7!

rm(list=ls())
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

library(data.table)
library(gtools)
# turning warnings into errors
options(warn=2)


#################################################
uniprot_all <- fread('../data/processed_complex_final.txt', sep='\t', header=TRUE)
allpdbs <- unique(tolower(uniprot_all$ID))
cat(sprintf("Total number of PDB files needed: %d\n", length(allpdbs)))
######################################################
## check which SIFTS data are already present
# Define directories
public_dir <- '/cosybio/project/public_data/SIFTS'
store_dir <- '../data/SIFTS'  # Local directory where we'll store files

# Create local directory if it doesn't exist
if(!dir.exists(store_dir)) {
    dir.create(store_dir, recursive=TRUE)
    cat(sprintf("Created local directory: %s\n", store_dir))
}

# Get list of files we need
if(dir.exists(public_dir)) {
    # Get list of files from public directory
    public_files <- list.files(public_dir)
    public_pdbs <- unlist(lapply(strsplit(public_files, '[.]'), '[[', 1))
    
    # Find which files we need to copy and which to download
    tocopy <- intersect(allpdbs, public_pdbs)
    todownload <- setdiff(allpdbs, public_pdbs)
    
    cat(sprintf("Found %d existing SIFTS files in public directory\n", length(tocopy)))
    cat(sprintf("Need to download %d additional files\n", length(todownload)))
    
    # Copy existing files from public directory
    if(length(tocopy) > 0) {
        cat("\nCopying existing files from public directory...\n")
        for(pdb in tocopy) {
            # Copy both .xml and .xml.gz files if they exist
            for(ext in c('.xml', '.xml.gz')) {
                src_file <- file.path(public_dir, paste0(pdb, ext))
                if(file.exists(src_file)) {
                    file.copy(src_file, file.path(store_dir, paste0(pdb, ext)), overwrite=TRUE)
                }
            }
        }
        cat(sprintf("Successfully copied %d files\n", length(tocopy)))
    }
} else {
    cat("Warning: Public SIFTS directory not found. Will download all files.\n")
    todownload <- allpdbs
}

if(length(todownload) == 0) {
    cat("All required SIFTS files are available. No downloads needed.\n")
    quit(save="no")
}

## download SIFT data
allsifts <- substr(todownload, 2,3) # extract only the 2nd and 3rd character, that constitute the SIFT

cat(sprintf("\nStarting download of %d SIFTS files...\n", length(todownload)))
download_success <- 0
download_failed <- 0

for(k in 1:length(todownload)){
    cat(sprintf("\nProcessing file %d of %d (%.1f%%): %s\n", 
        k, length(todownload), (k/length(todownload))*100, todownload[k]))
    
    output_name <- paste0(todownload[k],'.xml.gz')
    query <- paste0('https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/',allsifts[k],'/',output_name)
    cmd1 <- paste0('wget -O ',store_dir,'/',output_name,' ',query)
    cmd2 <- paste0("gunzip --force ", store_dir,'/',output_name)
    
    # Try download
    download_status <- system(cmd1)
    if(download_status == 0) {
        cat(sprintf("Successfully downloaded %s\n", output_name))
        # Try decompression
        decompress_status <- system(cmd2)
        if(decompress_status == 0) {
            cat(sprintf("Successfully decompressed %s\n", output_name))
            download_success <- download_success + 1
        } else {
            cat(sprintf("Warning: Failed to decompress %s\n", output_name))
            download_failed <- download_failed + 1
        }
    } else {
        cat(sprintf("Error: Failed to download %s\n", output_name))
        download_failed <- download_failed + 1
    }
}

# Print summary
cat(sprintf("\nFinal Summary:\n"))
cat(sprintf("Total files needed: %d\n", length(allpdbs)))
cat(sprintf("Files copied from public directory: %d\n", length(tocopy)))
cat(sprintf("Files downloaded: %d\n", length(todownload)))
cat(sprintf("Successfully downloaded and decompressed: %d\n", download_success))
cat(sprintf("Failed downloads or decompressions: %d\n", download_failed))

if(download_failed > 0) {
    cat("\nNote: Some files failed to download or decompress. Consider using the sifts_download.sh script for retry.\n")
}
