##############################################################################################
# Purpose: download all xml files of the cross mapping from SIFTS
##############################################################################################
## WARNING: If some errors in the download occur (some files cannot be downloaded, some are downloaded but not unzipped), please
## use the bash script 'sifts_download.sh' in this same folder. It's able to bypassing errors due to time requests.
## One protein complex resulted missing from the database: 7qe7!

rm(list=ls())
library(data.table)
library(gtools)
# turning warnings into errors
options(warn=2)


#################################################
uniprot_all <- fread('Mus_musculus/processed_complex_final.txt', sep='\t', header=TRUE)
allpdbs <- unique(tolower(uniprot_all$ID))
######################################################
## check which SIFTS data are already present
# check if the folder exists
store_dir <- '../public_data/SIFTS'
# to get all names of file we should download
if(dir.exists(store_dir)){
	allfiles <- list.files(store_dir)
	allpresent <- unlist(lapply(strsplit(allfiles, '[.]'), '[[', 1))
	todownload <- setdiff(allpdbs, allpresent)
}else{
  	dir.create(store_dir)
	todownload <- allpdbs
}

## to save in a file the new SIFTS files that should be downloaded
# filepath <- file.path(store_dir, "new_SIFTS.txt")
# fwrite(list(todownload), filepath, sep='\t', quote=FALSE)

## download SIFT data
allsifts <- substr(todownload, 2,3) # extract only the 2nd and 3rd character, that constitute the SIFT

for(k in 1:length(todownload)){

	output_name <- paste0(todownload[k],'.xml.gz')
	query <- paste0('https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/',allsifts[k],'/',output_name)
	cmd1 <- paste0('wget -O ',store_dir,'/',output_name,' ',query)
	cmd2 <- paste0("gunzip --force ", store_dir,'/',output_name)
	system(cmd1)
	system(cmd2)

}
