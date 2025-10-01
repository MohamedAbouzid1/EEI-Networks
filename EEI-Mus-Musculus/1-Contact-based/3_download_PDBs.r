##############################################################################################
# Purpose: Download CIF data and save fasta 
############################################################################################## 

rm(list=ls())
library('seqinr')
library('data.table')
library('stringr')
# NOT turning warnings into errors because conversion of warnign during conversion of 3-letter code to AA symbol
options(warn=0)

################################################################
# all filtered uniprot pdb map
uniprot_all <- fread('Mus_musculus/processed_complex_final.txt', sep='\t', header=TRUE)
keep <- unique(tolower(uniprot_all$ID))

## download CIFs ###
# check if the folder exists
store_cif <- '../public_data/PDB_CIF' 

if(dir.exists(store_cif)){
  pdbids2 <- unique(unlist(lapply(strsplit(list.files(store_cif), '[.]'), '[[', 1)))
  pdbids3 <- setdiff(keep, pdbids2)
  allpdbs <- paste(pdbids3, collapse=',')
}else{
  dir.create(store_cif)
  allpdbs <- paste(keep, collapse=',')
}

if(nchar(allpdbs) != 0){
  writeLines(allpdbs,'../public_data/all_new_CIFs.txt') # save in a txt file all CIFs that are newly download
  system(paste0('Contact-based/batch_download.sh -f ../public_data/all_new_CIFs.txt -o ',store_cif,' -c'))
  system(paste0('gunzip ',store_cif,'/*.cif.gz'))
}
