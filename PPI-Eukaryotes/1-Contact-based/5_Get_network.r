#################################################################
# Purpose: get networks based on distance of chains in PDBs
#################################################################

rm(list=ls())

library(data.table)
library(stringr)
library(plyr)
library(Rcpp)
# install.packages("ggplot2")
library(ggplot2)
library(igraph)
# turning warnings into errors
options(warn=2)

## PDB to unweighted network function
cppFunction("List pdb2net(CharacterVector uniqueids, NumericVector lastPositions, NumericVector xcoord, NumericVector ycoord, NumericVector zcoord, int cutoff){

	int loop1 = uniqueids.size()-1;
	int loop2 = uniqueids.size();
	int esize = loop2*loop2;
	CharacterVector p1(esize);
	CharacterVector p2(esize);
	NumericVector dis(esize);
	int start = 0;
	int counter = 0;

	for(int k=0; k<loop1; k++){

		int i = k+1;

		for(int j=i; j<loop2; j++){

			int startc = lastPositions[j-1]+1;
			int endc = lastPositions[j];
			int startr = start;
			int endr = lastPositions[k];
			double mindist = 100;

			for(int x=startr; x<=endr; x++){

				double xx = xcoord[x];
				double xy = ycoord[x];
				double xz = zcoord[x];

				for(int y=startc; y<=endc; y++){

					double yx = xcoord[y];
					double yy = ycoord[y];
					double yz = zcoord[y];

					double adist = sqrt(pow((yx-xx),2)+pow((yy-xy),2)+pow((yz-xz),2));
					if(adist < mindist){
						mindist = adist;
					}

				}
			}

			if(mindist <= cutoff){
				p1[counter] = uniqueids[k];
				p2[counter] = uniqueids[j];
				dis[counter] = mindist;
				counter = counter+1;
			}

		}
		start = lastPositions[k]+1;

	}

	List L = List::create(p1,p2,dis,counter);
	return L;
  
}")

## all mapped data
afile <- fread('data/database/processed_complex_final.txt', header=TRUE) # 60,516 entries
afile1 <- data.frame(chain1=paste0(afile$ID,'_',afile$CHAIN1),chain2=paste0(afile$ID,'_',afile$CHAIN2)) # dataframe for each chain-chain contact: PDB1_chain1 PDB1_chain2
ig <- graph_from_data_frame(afile1, directed = FALSE) # get a graph with nodes and edges
igg <- simplify(ig, remove.loops=FALSE, remove.multiple=TRUE) # remove repetitions and loops, if any
afile2 <- as.data.frame(as_edgelist(igg)) # get a dataframe from the simplified graph, updated function needed
afile2$pdbid <- tolower(unlist(lapply(strsplit(afile2[[1]], '[_]'), '[[', 1)))
afile2$chain1 <- unlist(lapply(strsplit(afile2[[1]], '[_]'), '[', 2)) # I modified from [[ to [
afile2$chain2 <- unlist(lapply(strsplit(afile2[[2]], '[_]'), '[', 2))

temp_pdb <- unique(afile2[[3]]) # collect pdb ids: 5732

pdbDirectory <- '../public_data/PDB_CIF'
cutoff <- c(6)#c(4, 5, 6, 7, 8)

for(mm in 1:length(cutoff)){

	store_dir <- paste0('data/database/networks_',cutoff[mm])
	if(!dir.exists(store_dir)){
		dir.create(store_dir)
	}

	
	for(k in 1:length(temp_pdb)){

		afile3 <- afile2[afile2$pdbid == temp_pdb[k], ]

		tfile <- readLines(paste0(pdbDirectory,"/",temp_pdb[k],".cif"))
		# Extract start and end positions of the coordinates entries
		wh <- which(tfile == "loop_")+1 # extract the indeces of rows being part of "loop_"
		tfile0 <- trimws(tfile[wh]) # extract all items
		whh1 <- which(tfile0 == "_atom_site.group_PDB") # only the item of our interest
		start <- wh[whh1] # the row from which what we need starts
		if(whh1 == length(wh)){
			end <- length(tfile)
		}else{
			end <- wh[whh1+1]-1-2  # save the end 
		}

		# Extract the coordinates part of the PDB file
		tfile <- tfile[start:end] # subset of CIF file
		lineID <- word(tfile, 1)
		wh <- which(lineID == "ATOM" | lineID == "HETATM") # select only the rows regarding ATOM or HETATM

		# Extract the field entries
		whf <- setdiff(seq(1,length(tfile)), wh)
		fields <- trimws(tfile[whf])

		tfile1 <- trimws(tfile[wh])
		tfile2 <- read.table(textConnection(tfile1), sep='', colClasses = "character")

		# chain id is author defined and not CIF defined
		chainPosition <- which(fields == "_atom_site.auth_asym_id")#_atom_site.label_asym_id
		chain <- tfile2[[chainPosition]]
		chain[is.na(chain)] <- "NA"

		#for each complex from this pdbid
		for(j in 1:length(afile3[[1]])){
			if(k == 5196 && j == 393) { next } else # skip that specific chain because no match between chains 
			{ 
				#filter using chain
				wh <- which(chain == as.character(afile3$chain1[j]) | chain == as.character(afile3$chain2[j])) 
				tfile3 <- tfile2[wh, ] # subset according to the entries regarding the same chain name

				#filter using model number. keep only the first model
				modelPosition <- which(fields == "_atom_site.pdbx_PDB_model_num")
				mdl <- unique(tfile3[[modelPosition]]) # select those tfile3 entries having that specific modelPosition
				wh <- which(tfile3[[modelPosition]] == mdl[1])
				tfile4 <- tfile3[wh, ]

				# extract coordinates for only the heavy atoms (S, O, C, N)
				atomPosition <- which(fields == "_atom_site.type_symbol")
				lineID2 <- tfile4[[atomPosition]]
				wh <- which(lineID2 == "C" | lineID2 == "S" | lineID2 == "O" | lineID2 == "N")
				tfile5 <- tfile4[wh, ]

				# keep only ATOM coordinates
				wh <- which(tfile5[[1]] == "ATOM")
				tfile5 <- tfile5[wh,]
				seqPosition <- which(fields == "_atom_site.label_seq_id") # CIF based, index of this item
				aaPosition <- which(fields == "_atom_site.label_comp_id") # get the AA symbol columnm, index of this item

				### transformation into single chain ###########################
				# get the receptor and ligand information
				temp <- tfile5
				chains <- temp[[chainPosition]]
				uchain <- unique(chains)
				chain1 <- which(chains %in% afile3$chain1[j])
				chain2 <- which(chains %in% afile3$chain2[j])

				# original sequences
				q1 <- temp[[seqPosition]][chain1]
				q2 <- temp[[seqPosition]][chain2]
				temp[chainPosition] <- rep('Z', length(temp[[1]]))
				counter <- 1
				pointer0 <- temp[[seqPosition]][1]
				new_seq <- c(1)

				for(i in 2:length(temp[[1]])){

					pointer1 <- temp[[seqPosition]][i]

					if(pointer1 != pointer0){
						counter <- counter+1
					}

					pointer0 <- pointer1
					new_seq <- c(new_seq, counter)
				}

				temp[seqPosition] <- new_seq

				# extract chain positions of the two proteins
				temp1 <- temp[chain1, ]
				temp2 <- temp[chain2, ]
				p1 <- temp1[[seqPosition]]
				p2 <- temp2[[seqPosition]]
				a1 <- temp1[[aaPosition]]
				a2 <- temp2[[aaPosition]]

				#call to create networks, extract coordinates
				seqq <- temp[[seqPosition]]
				wh <- which(fields == "_atom_site.Cartn_x")
				xcoord <- temp[[wh]]
				wh <- which(fields == "_atom_site.Cartn_y")
				ycoord <- temp[[wh]]
				wh <- which(fields == "_atom_site.Cartn_z")
				zcoord <- temp[[wh]]

				fname <- paste0(temp_pdb[k],'_',afile3$chain1[j],'_',afile3$chain2[j])
				uniqueids <- unique(seqq)
				lastPositions <- length(seqq)-match(unique(seqq),rev(seqq))
				xx <- pdb2net(as.character(uniqueids), as.numeric(lastPositions), as.numeric(xcoord), as.numeric(ycoord), as.numeric(zcoord), as.numeric(cutoff[mm]))


				Data <- data.frame(x=xx[[1]][1:xx[[4]]], y=xx[[2]][1:xx[[4]]])
				fwrite(Data,paste0(store_dir,'/',fname,'.txt'), quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
				
				pp1 <- data.frame(original=q1, new=p1, AA=a1)
				pp2 <- data.frame(original=q2, new=p2, AA=a2)

				fwrite(pp1,paste0(store_dir,'/',fname,'.chain1'), quote=FALSE, sep='\t', row.names=FALSE)
				fwrite(pp2,paste0(store_dir,'/',fname,'.chain2'), quote=FALSE, sep='\t', row.names=FALSE)
			}
		}

		cat('For cutoff:',mm,' :chain ',k, ' of ', length(temp_pdb), ' done\n')

	}

}
