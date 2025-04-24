##############################################################################################
# Purpose: download interaction files from PISA main script
##############################################################################################

rm(list=ls())

# all species
allpdbs <- data.table::fread('../data/database/process_complex_final.txt', header=TRUE)[[1]]

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

	tmpcmd <- paste('Rscript runpisa.r',start,end,store_dir)
	cmd <- paste('screen -dm',tmpcmd)
	system(cmd)

	toprocess <- toprocess-block
	start <- end + 1

}

