#######################################################################
# Purpose: mapping between Uniprot data and SIFTS file for each protein chain's residue
#######################################################################

rm(list=ls())
library(data.table)
library(gtools)
#install.packages("XML")
library(XML)
#install.packages("dplyr")
library(dplyr)
# turning warnings into errors
options(warn=0)

#################################################
# all filtered uniprot pdb map
uniprot_pdb1 <- fread('data/database/processed_complex_final.txt', sep='\t', header=TRUE) # upload uniprot_pdb file as a list
pdbids1 <- unique(tolower(uniprot_pdb1$ID)) # take only the IDs: 5732

####################
# PDB ID for which no SIFTS file was available are removed: '7qe7', '7yyq', '8by5'
a <- which(pdbids1 == '7qe7')
b <- which(pdbids1 == '7yyq')
c <- which(pdbids1 == '8by5')

pdbids1 <- pdbids1[-c(a, b, c)]
####################

# chainid from uniprot matches the PDB dbChainId (i.e., the author chain id), which can be 
# different from the corresponding entity.
# The entity corresponds to chain identity in the cif file
########################################################################################
## create uniprot and pdb map based on xmls
store_dir <- 'data/database/uniprotPDB_map'

if(!dir.exists(store_dir)){
	dir.create(store_dir)
}

for(k in 1:length(pdbids1)){  # for each of the 5729 entries
	# get the sub dataframe for this pdb
	tempdata <- uniprot_pdb1[uniprot_pdb1$ID == pdbids1[k], ] # save the line in the list

	temp <- xmlParse(paste0('../public_data/SIFTS/',pdbids1[k],'.xml')) # save the SIFTS file

	residue_set <- getNodeSet(temp, "//rs:residue[@dbSource='PDBe']", "rs") # find matching nodes in an internal XML tree, in particular creates an XMl file version of the SIFTS

	PDBeResNum <- c()
	uniprotId <- c()
	chainId <- c()
	uniprotResId <- c()
	uniprotResNum <- c()

	for(j in 1:length(residue_set)){ # lista dei residues, gli AA... for each residue:

		tempr <- xmlToList(residue_set[[j]])

 		wh <- which(names(tempr) == 'crossRefDb') # get the index of elements which match with 'crossRefDb'
		temprr <- tempr[wh] # save the elements in the given indeces

		# take names of the dbsource
		tempdbs <- unlist(unname(lapply(temprr, function(x) unname(x['dbSource']))))

		# check which one has 'UniProt'
		whu <- which(tempdbs == 'UniProt') # if there exist the UniProt source for that residue

		if(length(whu) != 0){

			# check whether PDB resolution is mapped
			whp <- which(tempdbs == 'PDB') # if there exist the PDB source too
			mappedNum <- unname(temprr[[whp]]['dbResNum']) # save the residue number reported by PDB

			# if some pdb residue number is mapped
			if(mappedNum != 'null'){

				# store author chain id
				chainId <- c(chainId, temprr[[whp]]['dbChainId']) # save the name of the chain (f.e. A) reported by PDB

				# Store author residue number
				PDBeResNum <- c(PDBeResNum, mappedNum) # at the end it would be equal to the list of numbers of all residues

				# UniProt data, for each residue. We got a list of length equal to the sum of all chains' residues
				uniprotId <- c(uniprotId, temprr[[whu]]['dbAccessionId']) # Uniprot IDs of the chain
				uniprotResNum <- c(uniprotResNum, temprr[[whu]]['dbResNum']) # Residue's number
				uniprotResId <- c(uniprotResId, temprr[[whu]]['dbResName']) # Residue's name (which AA)

			}

		}

	}



	chainids1 <- tempdata$CHAIN1 # chain1 names of a specific PDB
	chainids2 <- tempdata$CHAIN2 # chain2 names of a specific PDB

	unids1 <- tempdata$PROTEIN1 # Uniprot names of PROTEIN1 a specific PDB
	unids2 <- tempdata$PROTEIN2 # Uniprot names of PROTEIN2 a specific PDB

	for(j in 1:length(unids1)){ # for each PROTEIN1
									#save info for first protein and chain

		# which chain to keep
		whc <- which(chainId == chainids1[j]) # indeces of residues belonging to that specific chain
		uniprotId1 <- unname(uniprotId[whc]) # removes names or dimnames from R attribute
		uniprotResNum1 <- unname(uniprotResNum[whc])
		uniprotResId1 <- unname(uniprotResId[whc])
		PDBeResNum1 <- unname(PDBeResNum[whc])
		# which uniprot to keep
		whc <- which(uniprotId1 == unids1[j]) # indeces of residues having a specific Uniprot name as PROTEIN1
		uniprotId1 <- unname(uniprotId1[whc])
		uniprotResNum1 <- unname(uniprotResNum1[whc])
		uniprotResId1 <- unname(uniprotResId1[whc])
		PDBeResNum1 <- unname(PDBeResNum1[whc])
		# save file
		Data1 <- data.frame(UNIPROT_SEQ_NUM=uniprotResNum1, UNIPROT=uniprotResId1,PDBResNumAuthor=PDBeResNum1)
		fwrite(Data1, paste0(store_dir,'/',paste0(unids1[j],'_',pdbids1[k],'_',chainids1[j]),'.txt'), sep='\t', quote=FALSE, row.names=FALSE)

									#save info for first protein and chain

		# which chain to keep
		whc <- which(chainId == chainids2[j])
		uniprotId2 <- unname(uniprotId[whc])
		uniprotResNum2 <- unname(uniprotResNum[whc])
		uniprotResId2 <- unname(uniprotResId[whc])
		PDBeResNum2 <- unname(PDBeResNum[whc])
		# which uniprot to keep
		whc <- which(uniprotId2 == unids2[j])
		uniprotId2 <- unname(uniprotId2[whc])
		uniprotResNum2 <- unname(uniprotResNum2[whc])
		uniprotResId2 <- unname(uniprotResId2[whc])
		PDBeResNum2 <- unname(PDBeResNum2[whc])
		# save file
		Data2 <- data.frame(UNIPROT_SEQ_NUM=uniprotResNum2, UNIPROT=uniprotResId2,PDBResNumAuthor=PDBeResNum2)
		fwrite(Data2, paste0(store_dir,'/',paste0(unids2[j],'_',pdbids1[k],'_',chainids2[j]),'.txt'), sep='\t', quote=FALSE, row.names=FALSE)
		
	}

	cat('Protein ', k, ' of ', length(pdbids1), ' done\n')

}
