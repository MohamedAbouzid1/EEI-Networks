##############################################################################################
# Purpose: add CIF number to the files in which Uniprot and PDB are mapped using SIFTS
##############################################################################################

rm(list=ls())
library(data.table)
library(stringr)
# turning warnings into errors
options(warn=2)

#################################################
# all filtered uniprot pdb map
allfiles1 <- list.files('data/database/uniprotPDB_map', full.names=TRUE) # 20 552 protein chains

store_dir <- 'data/database/uniprotPDB_map_final'
if(dir.exists(store_dir)){unlink(store_dir, recursive=TRUE)}
dir.create(store_dir)

## Map CIF numbers ############
#####################################################################
cif_folder <- '../public_data/PDB_CIF'
allfiles <- list.files(cif_folder, full.names=TRUE)

for(k in 1:length(allfiles1)){ # for each of the recently created PDB_map files

	temp <- strsplit(basename(allfiles1[k]),'[_]')[[1]] # split the name of the file
	temp_uni <- temp[1] # save UniProt ID
	temp_pdb <- temp[2] # save PDB ID
	temp_chain <- strsplit(temp[3],'[.]')[[1]][1] # save name of the chain

	wh <- which(allfiles %like% temp_pdb) # save the index of the CIF file referring to that PDB ID
	tfile <- readLines(allfiles[wh]) # read that file

	# Extract start and end positions of the coordinates entries
	wh <- which(tfile == "loop_")+1 # save the index of the lines right after "loop_"
	tfile0 <- trimws(tfile[wh]) # save the content of that lines
	whh1 <- which(tfile0 == "_atom_site.group_PDB") # save the index of the line collecting "_atom_site.group_PDB", with respect all the "first line right after loop_"
	start <- wh[whh1] # save the index of the line collecting "_atom_site.group_PDB", with respect the whole file
	if(whh1 == length(wh)){ # save the end ...?
		end <- length(tfile)
	}else{
		end <- wh[whh1+1]-1-2
	}

	# Extract the coordinates part of the PDB file
	tfile <- tfile[start:end] # save only the line in that interval START:END
	lineID <- word(tfile, 1) # select the first word
	wh <- which(lineID == "ATOM") # save only the indeces of lines having as first word "ATOM"

	# Extract the field entries
	whf <- setdiff(seq(1,length(tfile)), wh) # select the indeces of the lines not starting with "ATOM", that are the lines containing the field name
	fields <- trimws(tfile[whf]) # save the name of the fields

	tfile1 <- trimws(tfile[wh]) # save the lines starting with "ATOM", all in character format
	tfile2 <- read.table(textConnection(tfile1), sep='', colClasses = "character") # fix the visualization, getting a list 

	# author seq num
	auth_seq <- which(fields == "_atom_site.auth_seq_id") # save the index of the column named "_atom_site.auth_seq_id" --> V17

	# cif seq num
	cif_seq <- which(fields == "_atom_site.label_seq_id") # save the index of the column named "_atom_site.label_seq_id" --> V9


	#filter using chain
	chainPosition <- which(fields == "_atom_site.auth_asym_id") # author-based chain ID, save the index of the column named "_atom_site.auth_asym_id" --> V19
	chain <- tfile2[[chainPosition]] # extract all chains, column V19 
	chain[is.na(chain)] <- "NA" # remove from the vector all chain names == "NA"
	wh <- which(chain == temp_chain) # save the index of the lines having as chain the same chain collected from the PDB_map file
	tfile3 <- tfile2[wh, ] # subset tfile2 selecting only the lines with tha same chain collected from the PDB_map file
	tfile4 <- unique(tfile3[, c(cif_seq,auth_seq)]) # subset tfile3: only V9 and V17 columns and selecting only 1 line for each identical couple of cif_seq and auth_seq... no repetitions! 


	# map to the mapping file
	if(file.info(allfiles1[k])$size != 0){ # only if the file is not empty (14 will be excluded)

		allmap <- fread(allfiles1[k]) # read the file
		toadd <- rep('-', length(allmap[[1]])) # create a vector that will be needed later on

		# to catch cases where PDBResNum has a character but that is not always present in the CIF file
		tempnum <- unlist(str_extract_all(as.character(allmap$PDBResNumAuthor), "[0-9]+")) # extract PDBResNumAuthor from the PDB_map file (3rd column)

		# to catch negative numbers
		whneg <- which(substr(as.character(allmap$PDBResNumAuthor), 1, 1) == '-')
		if(length(whneg) != 0){
			tempnum[whneg] <- paste0('-',tempnum[whneg])
		}

		whi <- intersect(tempnum, tfile4[[2]])
		wh1 <- which(tempnum %in% whi)
		wh2 <- which(tfile4[[2]] %in% whi)
		if(length(wh1)!= length(wh2)){next} # skipping cases where PDBResNum has a character but that is not always present in the CIF file, they should be identical!
		toadd[wh1] <- tfile4[[1]][wh2]

		allmap$PDBResNumCIF <- toadd
		fwrite(allmap, paste0(store_dir,'/',basename(allfiles1[k])), sep='\t', quote=FALSE, row.names=FALSE) # create a file identical to the PDB_map but with an additional column referring to CIF number, V9!
	}
	cat('PDB', k, 'of', length(allfiles1), 'done\n')

}
