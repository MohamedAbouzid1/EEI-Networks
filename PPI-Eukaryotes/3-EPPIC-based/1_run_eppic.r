##############################################################################################
# Purpose: run eppic with the given pdb ids
######################################################################################################

rm(list=ls())

outputDir <- 'data/database/EPPIC'
if(!dir.exists(outputDir)){
	dir.create(outputDir)
}

# all species
allpdbs <- data.table::fread('data/database/processed_complex_final.txt', header=TRUE)[[1]]

##-- which ids are already present in either of the folders
all_ids <- list.files(outputDir, pattern='*.scores')

processed_ids <- union(allhuman_ids, all_ids)

toprocess <- setdiff(allpdbs, unlist(lapply(strsplit(processed_ids, '[.]'), '[[',1)))

##---

for(k in 1:length(toprocess)){

	cmd <- paste0('$HOME/eppic/eppic-cli/target/eppic-cli-3.4.2-SNAPSHOT/bin/eppic  -i ',toprocess[k],' -o ',outputDir,' -s -a 30')
	system(cmd)
	
}



