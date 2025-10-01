##############################################################################################
# Purpose: save uniprot entries with other info --> geneid, gene name, pdb chains, resolution of chains, start of chain, end of chain
##############################################################################################

rm(list=ls())
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

library(seqinr)
library(data.table)
library(Rcpp) #for using C++ in R
library(stringr)
library(plyr)
# turning warnings into errors
options(warn=2)

#############################################################################################################
###---- saving uniprot AA sequences
uniprot_dir <- '../data/uniprot_sprot_OC'

allfiles <- list.files(uniprot_dir, pattern = "_Reviewed.txt", full.names=TRUE)
output_file <- '../data/uniprot_sequences_OC.txt'
if(file.exists(output_file)){file.remove(output_file)}

# save all sequences in a vector
# create a unique FASTA file containing all protein AA sequences
for(k in 1:length(allfiles)){

	temp <- readLines(allfiles[k])
	tempseq <- substr(temp, 1,2) #line per line, it consider first and second character: the 2-letters-long category name

	wha <- which(tempseq == 'AC') #Accession Number: unique identifier assigned to each UniProt entry. There could be >1, also older versions.
	# take the canonical uniprot identifier
	temp1 <- unlist(lapply(strsplit(temp[wha[1]], '[;]'), '[[',1)) #save in a vector "AC" and the first element (PRIMARY AC) of the list which collects all Accession Numbers
	temp11 <- lapply(strsplit(temp1,'\\s+'), '[[',2) #save only the PRIMARY AC, as a list
	temp_uniprot <- temp11[[1]] #removing it from the list

	wh1 <- which(tempseq == 'SQ')+1 #select line where sequence starts
	wh2 <- which(tempseq == '//')-1 #select line where sequence ends

	temps1 <- gsub(' ','',paste(temp[wh1:wh2], collapse='')) #save the protein AA sequence
	write.fasta(temps1, temp_uniprot, output_file, open='a', nbchar=60, as.string=FALSE) #sequence, sequence name, the file in which we write, in Append way, putting 60 chars per line, as vectors of single characters (not as strings)
	cat(k,' of ', length(allfiles), ' done\n')

}


###########################################################################################################
# create a mapping of gene name, uniprot ids, geneid, pdb ids, start, and end of the uniprot sequence

# custom Pre-Processing function (C++ Function in R): defined in C++ but integrated in R using Rcpp package.

#############################################################
# Rcpp functions
# This is called: unimappdb1
# It takes as input 2 lists of data (PDB IDs [ex. 3HAO] & UniProt IDs [ex. P46952]) in order to produce a combined "flattened" structure of all combinations that is possible to collect from these lists
# pdbc is a vector of strings - tempuni is a vector of strings - resolution, start, end are vectors of numbers: resolution describes, in particular, the quality of the protein structure... <2Armstrong is often considered quite high - exp is a vector of strings
cppFunction("List unimappdb1(CharacterVector pdbc, CharacterVector tempuni, NumericVector resolution, NumericVector start, NumericVector end, CharacterVector exp){

	int loop1 = pdbc.size();
	int loop2 = tempuni.size();
	int fsize = loop1*loop2; //total number of possible combinations: it defines the size of output vectors
	CharacterVector uniprotid_filt(fsize);
	CharacterVector pdbc_filt(fsize);
	NumericVector resolution_filt(fsize);
	NumericVector start_filt(fsize);
	NumericVector end_filt(fsize);
	CharacterVector exp_filt(fsize);

	int counter = 0;

	for(int k=0; k<loop1; k++){ //iteration over each element of pdbc

		for(int j=0; j<loop2; j++){ //iteration over each element of tempuni
			//for each combination of PDB and UniProt
			uniprotid_filt[counter] = tempuni[j]; //store the following data...
			resolution_filt[counter] = resolution[k];
			pdbc_filt[counter] = pdbc[k];
			start_filt[counter] = start[k];
			end_filt[counter] = end[k];
			exp_filt[counter] = exp[k];
			counter = counter+1;

		}

	}

	List L = List::create(uniprotid_filt,resolution_filt, pdbc_filt, start_filt, end_filt,exp_filt); //at the end a list is returned, that is an R object
	return L;
  
}")

# This is called: unimappdb2
# The same as the unimappdb1 but in this case no information about resolution are stored
cppFunction("List unimappdb2(CharacterVector pdbc, CharacterVector tempuni, NumericVector start, NumericVector end, CharacterVector exp){

	int loop1 = pdbc.size();
	int loop2 = tempuni.size();
	int fsize = loop1*loop2;
	CharacterVector uniprotid_filt(fsize);
	CharacterVector pdbc_filt(fsize);
	NumericVector start_filt(fsize);
	NumericVector end_filt(fsize);
	CharacterVector exp_filt(fsize);

	int counter = 0;

	for(int k=0; k<loop1; k++){

		for(int j=0; j<loop2; j++){

			uniprotid_filt[counter] = tempuni[j];
			pdbc_filt[counter] = pdbc[k];
			start_filt[counter] = start[k];
			end_filt[counter] = end[k];
			exp_filt[counter] = exp[k];
			counter = counter+1;

		}

	}

	List L = List::create(uniprotid_filt, pdbc_filt, start_filt, end_filt,exp_filt);
	return L;
  
}")

#############################################################

#############################################################
# R User-Defined fuctions
# rproc1: using unimappdb1, creates, for each protein, a list of lists for those proteins entries having at least 2 chains and resolution > 3 
rproc1 <- function(temp2, resolution, pdbid, temp_uniprot){

	temp3 <- strsplit(trimws(unlist(lapply(temp2, '[[', 5))),'[=]') # save info about chains with resolved structure
	tempe <- trimws(unlist(lapply(temp2, '[[', 3))) # save experiment info
	## choose temp3 entries with at least two length because that means there is chain information
	temp3_count <- lengths(temp3) # return, for each element, the number of chains it has
	temp3_wh <- which(temp3_count == 2) # only taking the proteins with continuous information

	if(length(temp3_wh) > 0){

		temp3 <- temp3[temp3_wh] # update temp3 considering only prot with at least 2 chains
		tempe <- tempe[temp3_wh] # update tempe considering only prot with at least 2 chains
		resolution <- resolution[temp3_wh] # update resolution considering only prot with at least 2 chains
		pdbid <- pdbid[temp3_wh]  #update pdbid considering only prot with at least 2 chains
		temp4 <- unlist(lapply(temp3, '[[', 1)) # save only chains name "A/B..."
		chain <- unlist(lapply(strsplit(temp4, '[/]'),'[[',1)) # select the first chains for each prot
		temp5 <- unlist(lapply(temp3, '[[', 2)) # save residue numbers resolved "from...to..."
		temp6 <- paste0(unlist(lapply(strsplit(temp5, '[.]'),'[[',1)),'-') # formal modification: . --> -
		start <- unlist(lapply(strsplit(temp6, '[-]'),'[[',1)) # vector of all starting positions
		end <- unlist(lapply(strsplit(temp6, '[-]'),'[[',2)) # vector of all ending positions
		# replacing '' value in start and end vector
		start[start == ''] <- 0 # if empty position, set equal to 0
		end[end == ''] <- 0 # if empty position, set equal to 0
		wh10 <- which(resolution <= 3) # only taking proteins with resolution > 3

		pdbid <- pdbid[wh10] # updating...
		chain <- chain[wh10]
		resolution <- resolution[wh10]
		start <- start[wh10]
		end <- end[wh10]
		tempe <- tempe[wh10]
		pdbc <- paste0(tolower(pdbid), '_', chain)
		## all uniprots
		tempuni <- temp_uniprot
		
		tempxx <- unimappdb1(pdbc, tempuni, as.numeric(resolution), as.numeric(start), as.numeric(end), tempe)

		return(tempxx)

	}else{
		return(NULL)
	}

}

#rproc2: as rproc1 but it doesn't filter based on resolution. It's indeed used for those proteins' entries
## for which resolution wasn't present and that has been further set as 100
rproc2 <- function(temp2, pdbid, temp_uniprot){

	temp3 <- strsplit(trimws(unlist(lapply(temp2, '[[', 5))),'[=]')
	tempe <- trimws(unlist(lapply(temp2, '[[', 3)))

	# choose temp3 entries with at least two length because that means there is chain information
	temp3_count <- lengths(temp3)
	temp3_wh <- which(temp3_count == 2) # only taking the proteins with continuous information

	if(length(temp3_wh) > 0){

		temp3 <- temp3[temp3_wh]
		tempe <- tempe[temp3_wh]
		pdbid <- pdbid[temp3_wh]
		temp4 <- unlist(lapply(temp3, '[[', 1))
		chain <- unlist(lapply(strsplit(temp4, '[/]'),'[[',1))
		temp5 <- unlist(lapply(temp3, '[[', 2))
		temp6 <- paste0(unlist(lapply(strsplit(temp5, '[.]'),'[[',1)),'-') #
		start <- unlist(lapply(strsplit(temp6, '[-]'),'[[',1))
		end <- unlist(lapply(strsplit(temp6, '[-]'),'[[',2))

		# replacing '' value in start and end vector
		start[start == ''] <- 0 
		end[end == ''] <- 0
		
		pdbc <- paste0(tolower(pdbid), '_', chain)
		## all uniprots
		tempuni <- temp_uniprot
		
		tempxx <- unimappdb2(pdbc, tempuni, as.numeric(start), as.numeric(end), tempe)

		return(tempxx)

	}else{
		return(NULL)
	}
	
}

# vectors are initialized
uniprotid1 <- c()
genename1 <- c()
pdbchain1 <- c()
resol1 <- c()
allstart1 <- c()
allend1 <- c()
expr1 <- c()

uniprotid2 <- c()
genename2 <- c()
pdbchain2 <- c()
allstart2 <- c()
allend2 <- c()
expr2 <- c()

pdb_entry <- 0
pdb_str <- 0
pdb_xr <- 0
pdb_nmr <- 0
expr_info <- c()
pdb_files <- c()

for(k in 1:length(allfiles)){ # for each protein-specific file

	temp <- readLines(allfiles[k]) # read each line
	tempcc <- substr(temp, 1,2) # all category names (AC, ID ecc)
	
	# all uniprot ids
	wha <- which(tempcc == 'AC') 
	# take the canonical uniprot identifier
	temp1 <- unlist(lapply(strsplit(temp[wha[1]], '[;]'), '[[',1))
	temp11 <- lapply(strsplit(temp1,'\\s+'), '[[',2)
	temp_uniprot <- temp11[[1]] # saves the primary Accession Number

	# all pdb
	whp <- which(temp %like% ' PDB;') # return the line number of Database Reference lines which contain info about PDBs 

	if(length(whp) == 0){
		next # skip if no pdb entry is present
	}else{
		pdb_entry <- pdb_entry+1 # update the count by 1 if at least one PDB entry exist
		temp_pdb <- temp[whp] # create a list of all lines as strings
		temp2 <- strsplit(temp_pdb,'[;]') # separate strings creating a list of lists
		experiment <- trimws(unlist(lapply(temp2, '[[', 3))) # extract only the 3rd element of each list, the one refering to the tecnique used to collect that info
		expr_info <- union(expr_info, experiment) # save the unique experiment types (no repetitions)

		## only keeping entries with X-ray, EM or NMR
		whx <- which(experiment %in% c('X-ray', 'EM', 'NMR')) # save index of Xray/EM/NMR entry
		if(length(whx) != 0){ pdb_str <- pdb_str+1 } # update the count by 1 if at least one Xray/EM/NMR entry exist
		temp_pdb <- temp_pdb[whx] # DR lines of X-ray/EM/NMR, as a list of string
		temp2 <- strsplit(temp_pdb,'[;]') # split lines in a list of lists
		pdbid <- trimws(unlist(lapply(temp2, '[[', 2))) # save only the corresponding PDB IDs
		pdb_files <- union(pdb_files,unique(pdbid)) # no repetitions
		resolution <- substr(trimws(unlist(lapply(temp2, '[[', 4))),1,4) # select Resolution values (1.50 A --> 1.50 )
		resolution[resolution == '-'] <- 100 # if there's no resolution, save it as 100
		resolution <- as.numeric(resolution)
		whr <- which(resolution <= 3) # save index of lines with resolution < 3

		if(length(whr) != 0){ 

			pdb_xr <- pdb_xr+1 # update the count by 1 if at least one entry with resolution < 3 exist
			temp_pdbx <- temp_pdb[whr] # save the corresponding lines
			temp2 <- strsplit(temp_pdbx,'[;]') #split lines in a list of list

			tempxx <- rproc1(temp2, resolution[whr], pdbid[whr], temp_uniprot)
			if(is.null(tempxx)){next}
			uniprotid1 <- c(uniprotid1, tempxx[[1]]) # save in different variables eache of the sublists
			pdbchain1 <- c(pdbchain1, tempxx[[3]])
			allstart1 <- c(allstart1, tempxx[[4]])
			allend1 <- c(allend1, tempxx[[5]])
			resol1 <- c(resol1, tempxx[[2]])
			expr1 <- c(expr1, tempxx[[6]])

		}else{
			## check if the resolution is 100
			whr <- which(resolution == 100)
			if(length(whr) != 0){

				pdb_nmr <- pdb_nmr+1

				temp_pdbx <- temp_pdb[whr]
				temp2 <- strsplit(temp_pdbx,'[;]')
	
				tempxx <- rproc2(temp2, pdbid[whr], temp_uniprot)
				if(is.null(tempxx)){next}
				uniprotid2 <- c(uniprotid2, tempxx[[1]])
				pdbchain2 <- c(pdbchain2, tempxx[[2]])
				allstart2 <- c(allstart2, tempxx[[3]])
				allend2 <- c(allend2, tempxx[[4]])
				expr2 <- c(expr2, tempxx[[5]])

			}
		}

	}
	cat(k,' of ', length(allfiles), ' done\n')
}

# summarize all information in 2 data frames of 6 different columns: uniprotID, "PDBID"_"chain name", resolution, start, end, experiment 
## the first, id_map1 has entries with resolution < 3
## the second, id_map2 has entries having no resolution (detected using NMR experiments)
id_map1 <- data.frame(uniprotkbac=uniprotid1, pdbchain=pdbchain1, resolution=resol1, start=allstart1, end=allend1, exp=expr1)
id_map2 <- data.frame(uniprotkbac=uniprotid2, pdbchain=pdbchain2, resolution=rep(0,length(uniprotid2)), 
	start=allstart2, end=allend2, exp=expr2)


id_map1$PDBID <- unlist(lapply(strsplit(id_map1$pdbchain, '[_]'), '[[', 1)) # add a new column, "PDBID" only
id_map2$PDBID <- unlist(lapply(strsplit(id_map2$pdbchain, '[_]'), '[[', 1))
id_map1$CHAIN <- unlist(lapply(strsplit(id_map1$pdbchain, '[_]'), '[[', 2)) # add a new column, "chain name" only
id_map2$CHAIN <- unlist(lapply(strsplit(id_map2$pdbchain, '[_]'), '[[', 2))
id_map1$ID <- paste0(id_map1$PDBID) # add a new column, "ID" copying PDBID
id_map2$ID <- paste0(id_map2$PDBID)


##-- all pdb ids ---
# no repetition within each group and in their fusion too
pdbids1 <- unique(id_map1[[7]]) 
pdbids2 <- unique(id_map2[[7]]) 
allpdbids <- union(pdbids1, pdbids2) # their merging


##-- choose the PDBIDs that result in the highest number of complexes ---
allpdbs <- plyr::count(id_map1$ID) # count how many times each PDBID occurs
allpdbsu <- allpdbs[allpdbs$freq > 1, ] # select only those that appear in more than one complex

p1 <- c()
p2 <- c()
pdbid <- c()
chainid1 <- c()
chainid2 <- c()
resol <- c()
st1 <- c()
st2 <- c()
ed1 <- c()
ed2 <- c()

##-- For four of the PDB IDs, the PDB resolution info is not consistant (marginal differences in the values)
idx <- c()
for(k in 1:length(allpdbsu[[1]])){

	temp <- id_map1[id_map1$ID == allpdbsu[[1]][k], ] # extract from the data frame only lines referring to a specific ID
	# if(length(temp[[1]]) > 1){break}
	loop1 <- length(temp[[1]])-1 # measurements in order to know how many loop rounds we should do, for getting all combinations
	loop2 <- length(temp[[1]])
	loop3 <- (loop1*loop2)/2

	for(i in 1:loop1){
		t1 <- temp[i,]
		m <- i+1
		for(j in m:loop2){
			t2 <- temp[j,]
			p1 <- c(p1, t1$uniprotkbac) # save into variables information about the first element in the combination
			chainid1 <- c(chainid1, t1$CHAIN)
			st1 <- c(st1, t1$start)
			ed1 <- c(ed1, t1$end)
			p2 <- c(p2, t2$uniprotkbac) # save into variables information about the second element in the combination
			chainid2 <- c(chainid2, t2$CHAIN)
			st2 <- c(st2, t2$start)
			ed2 <- c(ed2, t2$end)
		}
	}

	pdbid <- c(pdbid, rep(unique(temp$ID), loop3)) # create a vector containing the pdbid (the same), repeated how many times as loop3 suggests
	# if(length(unique(temp$resolution)) > 1){idx <- c(idx,k)}
	resol <- c(resol, rep(unique(temp$resolution)[1], loop3))

	cat(k,' of ', length(allpdbsu[[1]]), ' done\n')

}

# different combinations of proteins, given a specific ID, are formed
cmx_data1 <- data.frame(ID=pdbid, RESOL=resol, PROTEIN1=p1, CHAIN1=chainid1, START1=st1, END1=ed1,
	PROTEIN2=p2, CHAIN2=chainid2, START2=st2, END2=ed2)

##--unique complexes
cmx1 <- igraph::simplify(igraph::graph_from_data_frame(cmx_data1[,c(3,7)], directed=FALSE)) # create a graph with proteins as Nodes and Edges in between couple proteins
cmx11 <- igraph::as_data_frame(cmx1)


##---choose a pdb id for each complex -----
allData <- data.frame(matrix(ncol=length(cmx_data1), nrow=0))
for(k in 1:length(cmx11[[1]])){
	x <- cmx11[[1]][k]
	y <- cmx11[[2]][k]
	wh1 <- which(cmx_data1$PROTEIN1 == x)
	wh2 <- which(cmx_data1$PROTEIN2 == y)
	wha <- intersect(wh1, wh2)
	wh1 <- which(cmx_data1$PROTEIN2 == x)
	wh2 <- which(cmx_data1$PROTEIN1 == y)
	whb <- intersect(wh1, wh2)
	wh <- union(wha, whb)
	temp1 <- cmx_data1[wh,]

	# if multiple
	# compute coverage of each ID-PROTEIN1-PROTEIN2 combinations
	if(nrow(temp1) > 1){
		cov1 <- c()
		cov2 <- c()
		for(j in 1:nrow(temp1)){
			temp2 <- temp1[j,]
			if(x == temp2$PROTEIN1){
				cov1 <- c(cov1, abs(as.numeric(temp2$START1) - as.numeric(temp2$END1)))
				cov2 <- c(cov2, abs(as.numeric(temp2$START2) - as.numeric(temp2$END2)))
			}else{
				cov1 <- c(cov1, abs(as.numeric(temp2$START2) - as.numeric(temp2$END2)))
				cov2 <- c(cov2, abs(as.numeric(temp2$START1) - as.numeric(temp2$END1)))
			}
		}

		# max coverage
		cov1_mx <- which(cov1 == max(cov1))
		cov2_mx <- which(cov2 == max(cov2))
		cov_f <- intersect(cov1_mx, cov2_mx)

		# if no intersection
		if(length(cov_f) == 0){
			if(cov1[cov1_mx[1]] < cov2[cov2_mx[1]]){ # choose the position that has higher coverage for the smaller protein
				# because this says that we choose loss of info in longer protein than the shorter proteins, 
				# potetially leading to less percentage-loss of a given protein
				cov_f <- cov1_mx[1]
			}else{
				cov_f <- cov2_mx[1]
			}
		}
		temp3 <- temp1[cov_f, ]

		# if now more than one row, then select by pdb resolution
		whr <- which(temp3$RESOL == min(temp3$RESOL))
		temp4 <- temp3[whr[1], ]
		allData <- rbind(allData, temp4)
	}else{
		allData <- rbind(allData, temp1)
	}
	cat(k,' of ', length(cmx11[[1]]), ' done\n')
}

netpdb1 <- unique(allData[[1]])
# -- choose the PDBIDs that results in the highest number of complexes ---
allpdbs <- plyr::count(id_map2$ID)
allpdbsu <- allpdbs[allpdbs$freq > 1, ]

p1 <- c()
p2 <- c()
pdbid <- c()
chainid1 <- c()
chainid2 <- c()
resol <- c()
st1 <- c()
st2 <- c()
ed1 <- c()
ed2 <- c()

for(i in 1:loop1){
    t1 <- temp[i,]
    m <- i+1
    for(j in m:loop2){
        t2 <- temp[j,]
        p1 <- c(p1, t1$uniprotkbac)
        chainid1 <- c(chainid1, t1$CHAIN)
        st1 <- c(st1, t1$start)
        ed1 <- c(ed1, t1$end)
        p2 <- c(p2, t2$uniprotkbac)
        chainid2 <- c(chainid2, t2$CHAIN)
        st2 <- c(st2, t2$start)
        ed2 <- c(ed2, t2$end)
        pdbid <- c(pdbid, unique(temp$ID))
        resol <- c(resol, unique(temp$resolution))
    }
}

cmx_data2 <- data.frame(ID=pdbid, RESOL=resol, PROTEIN1=p1, CHAIN1=chainid1, START1=st1, END1=ed1,
    PROTEIN2=p2, CHAIN2=chainid2, START2=st2, END2=ed2)
##--unique complexes
cmx2 <- igraph::simplify(igraph::graph_from_data_frame(cmx_data2[,c(3,7)], directed=FALSE))
cmx11 <- igraph::as_data_frame(cmx2)


##---choose a pdb id for each complex -----
allData2 <- data.frame(matrix(ncol=length(cmx_data2), nrow=0))
for(k in 1:length(cmx11[[1]])){
	x <- cmx11[[1]][k]
	y <- cmx11[[2]][k]
	wh1 <- which(cmx_data2$PROTEIN1 == x)
	wh2 <- which(cmx_data2$PROTEIN2 == y)
	wha <- intersect(wh1, wh2)
	wh1 <- which(cmx_data2$PROTEIN2 == x)
	wh2 <- which(cmx_data2$PROTEIN1 == y)
	whb <- intersect(wh1, wh2)
	wh <- union(wha, whb)
	temp1 <- cmx_data2[wh,]

	# if multiple
	if(nrow(temp1) > 1){
		cov1 <- c()
		cov2 <- c()
		for(j in 1:nrow(temp1)){
			temp2 <- temp1[j,]
			if(x == temp2$PROTEIN1){
				cov1 <- c(cov1, abs(as.numeric(temp2$START1) - as.numeric(temp2$END1)))
				cov2 <- c(cov2, abs(as.numeric(temp2$START2) - as.numeric(temp2$END2)))
			}else{
				cov1 <- c(cov1, abs(as.numeric(temp2$START2) - as.numeric(temp2$END2)))
				cov2 <- c(cov2, abs(as.numeric(temp2$START1) - as.numeric(temp2$END1)))
			}
		}

		# max coverage
		cov1_mx <- which(cov1 == max(cov1))
		cov2_mx <- which(cov2 == max(cov2))
		cov_f <- intersect(cov1_mx, cov2_mx)

		# if no intersection
		if(length(cov_f) == 0){
			if(cov1[cov1_mx[1]] < cov2[cov2_mx[1]]){ # choose the position that has higher coverage for the smaller protein
				# because this says that we choose loss of info in longer protein than the shorter proteins, 
				# potetially leading to less percentage-loss of a given protein
				cov_f <- cov1_mx[1]
			}else{
				cov_f <- cov2_mx[1]
			}
		}
		temp3 <- temp1[cov_f, ]
		# if now more than one row, then select by pdb resolution
		whr <- which(temp3$RESOL == min(temp3$RESOL))
		temp4 <- temp3[whr[1], ]
		allData2 <- rbind(allData2, temp4)
	}else{
		allData2 <- rbind(allData2, temp1)
	}
	cat(k,' of ', length(cmx11[[1]]), ' done\n')
}

netpdb2 <- unique(allData2[[1]])
netpdb <- union(netpdb1, netpdb2)

finaldata <- rbind(allData2, allData)
fwrite(finaldata, '../data/processed_complex_final.txt', sep='\t', row.names=FALSE, quote=FALSE)
