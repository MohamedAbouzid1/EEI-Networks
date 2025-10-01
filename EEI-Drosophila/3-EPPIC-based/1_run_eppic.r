##############################################################################################
# Purpose: run eppic with the given pdb ids
######################################################################################################

rm(list=ls())

outputDir <- '/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Drosophila/3-EPPIC-based/results'
if(!dir.exists(outputDir)){
  dir.create(outputDir)
}

allpdbs <- tolower(data.table::fread('/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Drosophila/PISA_data/PDBs_specific.txt', header=FALSE)[[1]])

##-- which ids are already present in either of the folders (allspecies or human only)
allids <- list.files(outputDir, pattern='*.scores')

toprocess <- setdiff(allpdbs, unlist(lapply(strsplit(allids, '[.]'), '[[',1)))

##---
## EPPIC for human not giving output for the following eight PDBIDS: ???? c("7uv5", "6l9z", "2j6f", "5jcz", "2j6o", "7aew")

for(k in 1:length(toprocess)){
  
  # Define path to EPPIC JAR file
  eppic_jar <- '/cosybio/project/mabouzid/EEI_networks/software/eppic/eppic-cli/target/uber-eppic-cli-3.4.2-SNAPSHOT.jar'
  
  # Create path to input CIF file
  input_file <- paste0('/cosybio/project/public_data/PDB_CIF/', toprocess[k], '.cif')
  
  # Build complete command
  cmd <- paste0('java -jar "', eppic_jar, '" -i "', input_file, '" -o "', outputDir, '" -s -a 30')
  
  # Print command being executed (helpful for debugging)
  cat("Running:", cmd, "\n")
  
  # Execute command
  system(cmd)
  
  # Print progress
  cat("Processed", k, "of", length(toprocess), ":", toprocess[k], "\n")
}