##############################################################################################
# Purpose: download interaction files from PISA main script
##############################################################################################

rm(list=ls())

# human only
#allpdbs <- data.table::fread('../data/uniprot_pdb_Ensembl_finalized.txt', header=TRUE)[[1]]

# store_dir <- '../data/PISA'
# if(!dir.exists(store_dir)){
# 	dir.create(store_dir)
# }

# all species
allpdbs <- data.table::fread('../data/database/all_new_final_file.txt', header=TRUE)[[1]]

store_dir <- '../data/database/PISA'
if(!dir.exists(store_dir)){
	dir.create(store_dir)
}

start <- 1
toprocess <- length(allpdbs)
block <- 100

while(toprocess >= 0){

	end <- start + block - 1

	if(end > length(allpdbs)){
		end <- length(allpdbs)
	}

	#tmp
	cmd <- paste('Rscript runpisa.r',start,end,store_dir)
	#cmd <- paste('screen -dm',tmpcmd)
	system(cmd)

	toprocess <- toprocess-block
	start <- end + 1

}

