##############################################################################################
# Purpose: map the interfaces to the corrresponding exons
## contact files are in CIF numbering and the interface files are in both CIF and PDB numberings
## Creating EXON-EXON interactions as follows:
## 1. Map interface residues of protein 1 to its exons
## 2. Map interface residues of protein 2 to its exons
## 3. Join each exon in step 1 to each exon in step 2
######################################################################################################
rm(list=ls())
library(Rcpp)
library(data.table)

# Set C++ standard to C++11
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

# Initialize vectors to track skipped PDBs
skipped_no_maps <- character()
skipped_no_scores <- character()
skipped_no_contacts <- character()
skipped_no_interfaces <- character()

cppFunction("List getexons(CharacterVector id1, CharacterVector id2, CharacterVector ex){

    int loop1 = id1.size();
    int loop2 = id2.size();
    CharacterVector exons;

    for(int k=0; k<loop1; k++){
        for(int j=0; j<loop2; j++){

            if((id1[k] == id2[j])){
                exons.push_back(ex[j]);
                break;
            }
        }
    }

    List L = List::create(exons);
    return L;
}")

allmaps <- list.files('/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Oryctolagus/data/uniprot_EnsemblExonPDB_map', full.names=TRUE) # 8970 protein chains

## These pdb ids were not processed by the EPPIC software
not_captured <- c("7uv5", "6l9z", "2j6f", "5jcz", "2j6o", "7aew") # "1n6j", "2as5",  "6xdl" were missing before, not now; "7uv5" is missing now, not before

inputDir <- '/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Oryctolagus/3-EPPIC-based/results/'

allpdbs <- data.table::fread('/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Oryctolagus/data/uniprot_pdb_Ensembl_finalized.txt', header=TRUE)[[1]] # 3592 entries

processed_pdbs <- setdiff(allpdbs, not_captured) 
##--- 2 479  out of 2 485 unique PDB IDs remain

createEEI <- function(tsfile, pdbid, cfile, ifile1, ifile2, allmaps){
	##--- select the map files relevant for this PDB id
	maps <- allmaps[which(allmaps %like% basename(pdbid))] # look for PROT_PDB_chain.tx files
	
	# Debug information
	cat("PDB ID:", pdbid, "\n")
	cat("Number of maps found:", length(maps), "\n")
	cat("Type of maps:", class(maps), "\n")
	
	# Return early if no maps found
	if(length(maps) == 0) {
		cat("No map files found for PDB ID:", pdbid, "- skipping\n")
		return(list(character(0), character(0), numeric(0), numeric(0), numeric(0), numeric(0), 
			numeric(0), numeric(0), character(0), character(0), character(0)))
	}
	
	if(length(maps) > 0) {
		cat("First map file:", maps[1], "\n")
		cat("Type of first map:", class(maps[1]), "\n")
	}
	
	relevant_chains <- trimws(unlist(lapply(strsplit(unlist(lapply(strsplit(basename(maps),'[_]'), '[[', 3)), '[.]'),'[[',1))) # select only the chains
	temp_chains <- strsplit((tsfile$chains),'[+]')[[1]]
	ichains <- intersect(relevant_chains, temp_chains)
	EXON1 <- c()
	EXON2 <- c()
	ASA <- c()
	BSA <- c()
	int1 <- c()
	int2 <- c()
	pro1 <- c()
	pro2 <- c()
	CS_score <- c()
	CR_score <- c()
	PDBID <- c()
	if(length(ichains) == 2 ){ # if both chains of the tsfile is present in my map file
		whx1 <- which(relevant_chains == temp_chains[1]) # assign to each chain an index
		whx2 <- which(relevant_chains == temp_chains[2])
		for(mm in 1:length(whx1)){
			for(nn in 1:length(whx2)){
				map1 <- maps[whx1[mm]] # select the file specific for that chain
				p1 <- strsplit(basename(map1),'[_]')[[1]][1] # save in p1 variable only the protein name
				map2 <- maps[whx2[nn]]
				p2 <- strsplit(basename(map2),'[_]')[[1]][1]
				res1 <- unique(cfile[[1]]) # save residues number for chain 1
				res2 <- unique(cfile[[4]]) # save residues number for chain 2
				temp1 <- data.table::fread(map1) # read file for chain 1
			    whr <- which(temp1$PDBResNumCIF %in% res1) # select only lines of the file referring to the residues involved in the interface
			    temp1 <- temp1[whr, ] # subset the file according to residues involved in the interface
			    temp2 <- data.table::fread(map2)
			    whr <- which(temp2$PDBResNumCIF %in% res2)
			    temp2 <- temp2[whr, ]
			    ifile1 <- ifile1[ifile1$V1 %in% temp1$PDBResNumCIF, ] # filter the interface file of EPPIC too
			    ifile2 <- ifile2[ifile2$V1 %in% temp2$PDBResNumCIF, ]
			    ex1 <- getexons(ifile1[[1]], temp1$PDBResNumCIF, temp1$EXON) # select exons only 
			    ifile1$EXON <- ex1[[1]] # filter the interface file of EPPIC based on exons name
				ex2 <- getexons(ifile2[[1]], temp2$PDBResNumCIF, temp2$EXON)
			    ifile2$EXON <- ex2[[1]]
			    exon1 <- unique(ifile1$EXON) # no repetition in exons name
			    exon2 <- unique(ifile2$EXON)
			    if((length(exon1) != 0) & (length(exon2) != 0)){
				    for(i in 1:length(exon1)){
				    	for(j in 1:length(exon2)){
					    	r1 <- ifile1[ifile1$EXON == exon1[i], ]
					    	r2 <- ifile2[ifile2$EXON == exon2[j], ]
				    		EXON1 <- c(EXON1, exon1[i])
				    		EXON2 <- c(EXON2, exon2[j])
				    		int1 <- c(int1, length(r1[[1]])) 
				    		int2 <- c(int2, length(r2[[1]]))
				    		pro1 <- c(pro1, p1)
				    		pro2 <- c(pro2, p2)
				    		## taking the ASA and BSA values from the interface files
				    		## ignoring contact information --> same as in PISA
				    		ASA <- c(ASA, sum(as.numeric(r1$V4))+sum(as.numeric(r2$V4)))
				    		BSA <- c(BSA, sum(as.numeric(r1$V5))+sum(as.numeric(r2$V5)))
				    		PDBID <- c(PDBID, pdbid)

				    	}
				    }
				    if(tsfile$`call-cr` != 'nopred'){ # the interface should be biologically relevant, "bio"
				    	CR_score <- rep(tsfile$`sc-cr`, length(EXON1))
				    }else{
				    	CR_score <- rep(NA, length(EXON1))
				    }
				    if(tsfile$`call-cs` != 'nopred'){
				    	CS_score <- rep(tsfile$`sc-cs`, length(EXON1))
				    }else{
				    	CS_score <- rep(NA, length(EXON1))
				    }
				}
			}
		}
	} 
    return(list(EXON1, EXON2, ASA, BSA, CS_score, CR_score, int1, int2, pro1, pro2, PDBID))
}


EXON1 <- c()
EXON2 <- c()
ASA <- c()
BSA <- c()
CS_score <- c()
CR_score <- c()
INTFAA1 <- c()
INTFAA2 <- c()
PROTEIN1 <- c()
PROTEIN2 <- c()
PDB_ID <- c()

err_pdbs <- c("8jdk", "8fle", "8fl7", "8p2k", "6zm5", "8idt", "7o80") # PDBs not considered because an error with scan() function arise

for(k in 1:length(processed_pdbs)){
	if(!processed_pdbs[k] %in% err_pdbs) {
		pdbid <- processed_pdbs[k]
		
		# Check if scores file exists
		scores_file <- paste0(inputDir, pdbid, '.scores')
		if(!file.exists(scores_file)) {
			cat("Scores file not found for PDB ID:", pdbid, "- skipping\n")
			skipped_no_scores <- c(skipped_no_scores, pdbid)
			next
		}
		
		# Check if contacts file exists
		contacts_file <- paste0(inputDir, pdbid, '.contacts')
		if(!file.exists(contacts_file)) {
			cat("Contacts file not found for PDB ID:", pdbid, "- skipping\n")
			skipped_no_contacts <- c(skipped_no_contacts, pdbid)
			next
		}
		
		# Check if interfaces file exists
		interfaces_file <- paste0(inputDir, pdbid, '.interfaces')
		if(!file.exists(interfaces_file)) {
			cat("Interfaces file not found for PDB ID:", pdbid, "- skipping\n")
			skipped_no_interfaces <- c(skipped_no_interfaces, pdbid)
			next
		}
		
		score_file <- data.table::fread(scores_file) # read PDB.scores file
		score_file1 <- score_file[score_file$call == 'bio', ] # call specify the origin of a interface score: if == "bio" they are biologically relevant and not artifacts

		if(nrow(score_file1) != 0){

			wh <- which(score_file1$topo %like% '>') # topological type of interface
			if(length(wh > 0)){
				score_file2 <- score_file1[-wh, ] # remove the line if topo == '>'
			}else{
				score_file2 <- score_file1 # keep the file unchanged
			}
			sfile <- score_file2[, c(2,3,6,8,9,10,11,13,14,15,16,18,19,20,21,23,24,25)] # save specific columns

			## --read and process contact file
			contact_file <- readLines(contacts_file)
			ctokens1 <- stringr::word(contact_file, 1)
			ctokens2 <- stringr::word(contact_file, 2)

			## --read and process interface file
			interface_file <- readLines(interfaces_file)
			itokens1 <- stringr::word(interface_file, 1)
			itokens2 <- stringr::word(interface_file, 2)

			for(j in 1:length(sfile[[1]])){ # for each line in the score file
											# collect info from the contact files
				wh1 <- which(ctokens1 == '#') 
				wha1 <- which(ctokens2 == paste0(sfile$id[j],'\t'))
				wha <- intersect(wha1, wh1)+1
				lastp <- which(wh1 == wha1)
				if(lastp == length(wh1)){ ## handling the boundary case
					whb <- length(contact_file)
				}else{
					whb <- wh1[lastp+1]-1
				}
				cfile <- read.table(textConnection(contact_file[wha:whb]), sep='\t', colClasses = "character") # read properly the tokens of contact file, as a dataframe

				wh1 <- which(itokens1 == '#')
				wha1 <- which(itokens2 == paste0(sfile$id[j],'\t'))
				wha <- intersect(wha1, wh1)+4
				lastp <- which(wh1 == wha1)
				if(lastp == length(wh1)){ ## handling the boundary case
					whb <- length(interface_file)
				}else{
					whb <- wh1[lastp+1]-1
				}

				interface_filex <- interface_file[wha:whb]
				## separate the chains
				itokensy <- stringr::word(interface_filex, 1)
				wh <- which(itokensy == '##') # detect the index of ##, that are the characters used to separate strings
				interface_filex_ch1 <- interface_filex[1:(min(wh)-1)]
				interface_filex_ch2 <- interface_filex[(max(wh)+1):length(interface_filex)]

				itokensx1 <- stringr::word(interface_filex_ch1, -1) # select the last column (burial percent, that indicates the amount of protein involved in the interaction)
				itokensx1 <- as.numeric(unlist(lapply(itokensx1, function(x) substr(x, nchar(x)-4, nchar(x)))))
				wh <- which(itokensx1 <= 0.25) # remove the residues not in the interface, with burial percent <= 0.25
				if(length(wh) !=0){
					interface_filex_ch1 <- interface_filex_ch1[-wh]
				}

				itokensx2 <- stringr::word(interface_filex_ch2, -1)# select the last column (burial percent, that indicates the amount of protein involved in the interaction)
				itokensx2 <- as.numeric(unlist(lapply(itokensx2, function(x) substr(x, nchar(x)-4, nchar(x)))))
				wh <- which(itokensx2 <= 0.25) # remove the residues not in the interface, with burial percent <= 0.25
				if(length(wh) !=0){
					interface_filex_ch2 <- interface_filex_ch2[-wh]
				}

				ifile1 <- read.table(textConnection(interface_filex_ch1), sep='\t', colClasses = "character")
				ifile2 <- read.table(textConnection(interface_filex_ch2), sep='\t', colClasses = "character")

				tempcall <- createEEI(sfile[j,], pdbid, cfile, ifile1, ifile2, allmaps) # apply the previously define function

				if(length(tempcall[[1]]) != 0){
					EXON1 <- c(EXON1, tempcall[[1]])
					EXON2 <- c(EXON2, tempcall[[2]])
					ASA <- c(ASA, tempcall[[3]])
					BSA <- c(BSA, tempcall[[4]])
					CS_score <- c(CS_score, tempcall[[5]])
					CR_score <- c(CR_score, tempcall[[6]])
					INTFAA1 <- c(INTFAA1, tempcall[[7]])
					INTFAA2 <- c(INTFAA2, tempcall[[8]])
					PROTEIN1 <- c(PROTEIN1, tempcall[[9]])
					PROTEIN2 <- c(PROTEIN2, tempcall[[10]])
					PDB_ID <- c(PDB_ID, tempcall[[11]])
				} 
			}

		}
		cat('PDB file', k, 'of',length(processed_pdbs), 'done\n')
	} else { 
		cat(processed_pdbs[k], ' skipped (in err_pdbs list)!\n')
	}
}

# Print summary of skipped PDBs
cat("\nSummary of skipped PDBs:\n")
cat("No map files found for:", length(skipped_no_maps), "PDBs\n")
if(length(skipped_no_maps) > 0) cat(paste(skipped_no_maps, collapse=", "), "\n\n")
cat("No scores files found for:", length(skipped_no_scores), "PDBs\n")
if(length(skipped_no_scores) > 0) cat(paste(skipped_no_scores, collapse=", "), "\n\n")
cat("No contacts files found for:", length(skipped_no_contacts), "PDBs\n")
if(length(skipped_no_contacts) > 0) cat(paste(skipped_no_contacts, collapse=", "), "\n\n")
cat("No interfaces files found for:", length(skipped_no_interfaces), "PDBs\n")
if(length(skipped_no_interfaces) > 0) cat(paste(skipped_no_interfaces, collapse=", "), "\n\n")

# Save the results only if we have data
if(length(EXON1) > 0) {
	allData <- data.frame(exon1=EXON1, exon2=EXON2, AA1=INTFAA1, AA2=INTFAA2, protein1=PROTEIN1, protein2=PROTEIN2, 
		BuriedAreaAbs=BSA, SolAccAreaAbs=ASA, CS_SCORE=CS_score, CR_SCORE=CR_score, PDBID=PDB_ID)
	
	data.table::fwrite(allData, paste0('/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Oryctolagus/3-EPPIC-based/results/EPPIC_EEIN.txt'), row.names=FALSE, quote=FALSE, sep='\t')
} else {
	cat("\nNo data was collected. No output file was created.\n")
}



