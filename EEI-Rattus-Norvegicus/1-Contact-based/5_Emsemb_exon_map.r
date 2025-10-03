##############################################################################################
# Purpose: download Ensembl data, use it to get the transcript and sequence of each UniProt ID.
# For each transcript, detect its exons and save all info regarding them.
##############################################################################################

# rm(list=ls())
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

# Set base directory for all paths
base_dir <- "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main"

library(data.table)
### only once
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#BiocManager::install("biomaRt")
#BiocManager::install("pwalign")
###
library(biomaRt)
library(seqinr)
library(stringr)
library(Biostrings)
library(pwalign)



substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

cleanSeq <- function(x){
	last_flag <- substrRight(x, 1)
	if(last_flag == "*"){
		slen <- nchar(x)-1
		x <- substr(x, 1, slen)
	}
	return(x)
}

## Ensemble genome download -------------------------

xx <- "Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz"
system(paste0("wget -O ", file.path(base_dir, "EEI-Rattus-Norvegicus/data/", xx), " https://ftp.ensembl.org/pub/release-114/fasta/rattus_norvegicus/dna/", xx)) # download fasta file
system(paste0("gunzip -f ", file.path(base_dir, "EEI-Rattus-Norvegicus/data/", xx)))

# download GTF file from Ensembl
xx <- "Rattus_norvegicus.GRCr8.114.gtf.gz"
system(paste0("wget -O ", file.path(base_dir, "EEI-Rattus-Norvegicus/data/", xx), " https://ftp.ensembl.org/pub/release-114/gtf/rattus_norvegicus/", xx))
system(paste0("gunzip -f ", file.path(base_dir, "EEI-Rattus-Norvegicus/data/", xx)))

cat('End of downloads!!\n')

# preprocess
system('sh ./preprocess_exons.sh')
exons <- fread(file.path(base_dir, "EEI-Rattus-Norvegicus/data/processed_exons.tmp"),sep='\t') # extract from the GTF only the exons

cds <- fread(file.path(base_dir, "EEI-Rattus-Norvegicus/data/processed_cds.tmp"),sep='\t') # extract from the GTF only the CDS

#### EXONS
xxs <- strsplit(exons$V9, '[;]') # get back the proper GTF structure: a list of lists, each one for each entry

## save gene_id (ENSG...)
wh1 <- lapply(xxs, function(x) which(x %like% 'gene_id')) # get the index of each sublist where the gene_id is contained
gene_id1 <- mapply(function(x, y) x[y], xxs, wh1) # extract the gene_id (ENSG....)
gene_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(gene_id1),'\\s+'), '[[', 2))) # fix some unwanted quotes...

## save transcript_id (ENST...)
wh1 <- lapply(xxs, function(x) which(x %like% 'transcript_id'))
transcript_id1 <- mapply(function(x, y) x[y], xxs, wh1)
transcript_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(transcript_id1),'\\s+'), '[[', 2)))

## save the exon_number
wh1 <- lapply(xxs, function(x) which(x %like% 'exon_number'))
exon_number1 <- mapply(function(x, y) x[y], xxs, wh1)
exon_number <- gsub('\\"', "", unlist(lapply(strsplit(trimws(exon_number1),'\\s+'), '[[', 2)))

## save transcript_biotype ("sRNa", "rRNA", "pseudogene", "protein_coding"...)
wh1 <- lapply(xxs, function(x) which(x %like% 'transcript_biotype'))
transcript_biotype1 <- mapply(function(x, y) x[y], xxs, wh1)
transcript_biotype <- gsub('\\"', "", unlist(lapply(strsplit(trimws(transcript_biotype1),'\\s+'), '[[', 2)))

## save exon_id (ENSE...)
wh1 <- lapply(xxs, function(x) which(x %like% 'exon_id'))
exon_id1 <- mapply(function(x, y) x[y], xxs, wh1)
exon_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(exon_id1),'\\s+'), '[[', 2)))

exons$gene_id <- gene_id # add new columns to exons object, respectively with all these vectors and their proper column name
exons$transcript_id <- transcript_id
exons$transcript_biotype <- transcript_biotype
exons$exon_number <- exon_number
exons$exon_id <- exon_id

exons_f <- exons[,-9] # remove column number 9, the longest, full of string information
fwrite(exons_f, file.path(base_dir, "EEI-Rattus-Norvegicus/data/Ensembl_exons.txt"), sep='\t', quote=FALSE, row.names=FALSE) # and save it in exons_f file
exons_f <- fread(file.path(base_dir, "EEI-Rattus-Norvegicus/data/Ensembl_exons.txt"), header=TRUE)

#repeat everything for CDS

#### CDS 
xxs <- strsplit(cds$V9, '[;]')
# gene_id
wh1 <- lapply(xxs, function(x) which(x %like% 'gene_id'))
gene_id1 <- mapply(function(x, y) x[y], xxs, wh1)
gene_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(gene_id1),'\\s+'), '[[', 2)))

# transcript_id
wh1 <- lapply(xxs, function(x) which(x %like% 'transcript_id'))
transcript_id1 <- mapply(function(x, y) x[y], xxs, wh1)
transcript_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(transcript_id1),'\\s+'), '[[', 2)))

# exon_number
wh1 <- lapply(xxs, function(x) which(x %like% 'exon_number'))
exon_number1 <- mapply(function(x, y) x[y], xxs, wh1)
exon_number <- gsub('\\"', "", unlist(lapply(strsplit(trimws(exon_number1),'\\s+'), '[[', 2)))

# transcript_biotype
wh1 <- lapply(xxs, function(x) which(x %like% 'transcript_biotype'))
transcript_biotype1 <- mapply(function(x, y) x[y], xxs, wh1)
transcript_biotype <- gsub('\\"', "", unlist(lapply(strsplit(trimws(transcript_biotype1),'\\s+'), '[[', 2)))

cds$gene_id <- gene_id
cds$transcript_id <- transcript_id
cds$transcript_biotype <- transcript_biotype
cds$exon_number <- exon_number
cds$nt_len <- (cds$V5-cds$V4)+1

cds_f <- cds[,-9]
fwrite(cds_f, file.path(base_dir, "EEI-Rattus-Norvegicus/data/Ensembl_exon_cds.txt"), sep='\t', quote=FALSE, row.names=FALSE)
cds_f <- fread(file.path(base_dir, "EEI-Rattus-Norvegicus/data/Ensembl_exon_cds.txt"), header=TRUE)

# get the uniprot emsembl mapping ##### USE OF BIOMART #########################
# try several times the connections... they may not work....
ensembl <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", mirror='useast') # useast uswest asia enables one to connect to a specified BioMart database and dataset hosted by Ensembl without having to specify the Ensembl URL.
cat('Ensembl object correctly created!')
# all data info
attrs <- c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id','transcript_biotype', 'uniprotswissprot')
####################################################################################

priassm <- file.path(base_dir, "EEI-Rattus-Norvegicus/data/Rattus_norvegicus.GRCr8.dna.toplevel.fa")

genome <- seqinr::read.fasta(priassm) # get the whole sequence

cat('Genome uploaded correctly!')

# for each uniprot id...
# all uniprot ids
uniprot_pdb <- fread(file.path(base_dir, "EEI-Rattus-Norvegicus/data/processed_complex_final.txt"), sep='\t', header=TRUE) # the file created in 1b

# correct for NA chain recognized as a missing value
uniprot_pdb[is.na(uniprot_pdb)] <- "NA" # the size stays invariate
			      
## remove the complexes where the proteins map to the same chain as per SIFTS mapping
to_keep <- c()
for(i in 1:length(uniprot_pdb[[1]])){ # 22309 pairs
	if(uniprot_pdb$CHAIN1[i] != uniprot_pdb$CHAIN2[i]){
		to_keep <- c(to_keep, i)
	}
}

uniprot_pdb <- uniprot_pdb[to_keep, ] # subset uniprot_pdb keeping only pairs having different chains
upro <- union(uniprot_pdb$PROTEIN1, uniprot_pdb$PROTEIN2) # all proteins
updb <- unique(uniprot_pdb$ID) # all PDB IDs, with no repetitions 
##-- 

uniprot_seqs <- read.fasta(file.path(base_dir, "EEI-Rattus-Norvegicus/data/uniprot_sequences_RN.txt"), seqtype="AA", whole.header=TRUE) # upload protein sequences
uniprotids1 <- uniprot_pdb$PROTEIN1
uniprotids2 <- uniprot_pdb$PROTEIN2
uniprotids <- union(uniprotids1, uniprotids2) 
alld <- getBM(attributes=attrs, filters='uniprotswissprot', values=uniprotids, mart=ensembl) # get from Ensembl the info regarding the UniProt IDs given

bh_transcripts1 <- rep('',length(uniprotids1))

# # countl <- c()
countl <- rep(0, length(uniprotids1))
count_gap <- rep(0, length(uniprotids1))
iden <- rep(0,length(uniprotids1))

tricky_vect <- c("X", "U", "O", "B", "Z", "J", "*")
excluded1 <- c()
excluded2 <- c()
exclud_others <- c()



for(k in 1:length(uniprotids1)){
	name <- uniprotids1[k]
	res <- tricky_vect %in% uniprot_seqs$name
	if(any(res)){
		excluded1 <- c(excluded1, name)
		cat('Protein', name, ' skipped because its UniProt sequence\n')
		next
	}
	if(k == 110 || k == 256){next} # they were generating some errors

	# length of the uniprot seq
	uni_seq <- getLength(uniprot_seqs[uniprotids1[k]]) # get the length of the sequence of the protein
	alldata <- alld[alld$uniprotswissprot == uniprotids1[k], ] # get the Ensembl info regarding that UniProt ID

	# only protein coding
	alldata1 <- alldata[alldata$transcript_biotype == 'protein_coding', ] # select only "protein coding"
	all_transcripts <- alldata1$ensembl_transcript_id # save the name of the transcript of those entries that are protein coding

	if(!identical(all_transcripts, character(0))) { 
        temp_seq <- biomaRt::getSequence(id=all_transcripts, type='ensembl_transcript_id', seqType='peptide', mart=ensembl)
        Sys.sleep(4)

		for(j in 1:length(all_transcripts)) {
			temp_cds <- cds_f[cds_f$transcript_id == all_transcripts[j], ]
			
			if(nrow(temp_cds) > 0) {
				dna_strand <- temp_cds$V7[[1]][1]
				temp_sum <- sum(temp_cds$nt_len)/3
				flags <- 0
				if(temp_sum == uni_seq) {
					flags <- flags+1
					
					## global alignment
					toalign1 <- temp_seq$peptide[which(temp_seq$ensembl_transcript_id == all_transcripts[j])] ## protein seqeunce of first transcript
					toalign1 <- cleanSeq(toalign1) # remove * at the end

					temps <- seqinr::getSequence(genome[temp_cds[[1]][1]])[[1]] 
					allseq <- c()
					for(i in 1:length(temp_cds[[1]])) {
						start <- temp_cds$V4[i]
						end <- temp_cds$V5[i]
						tempseq <- temps[start:end]
						if(dna_strand == "-"){
							tempseq <- rev(seqinr::comp(tempseq))
						}
						allseq <- c(allseq, tempseq) # get the sequence of all transcripts

					}
					t <- seqinr::translate(allseq)
                    if(any(tricky_vect %in% t)) {
                        excluded2 <- c(excluded2, name)
						cat('Protein', name, ' skipped because its biomaRt sequence\n')
                        break
					}
                    else {

						toalign2 <- toupper(paste(seqinr::translate(allseq),collapse=''))

						xxx <- pwalign::pairwiseAlignment(toalign1, toalign2, # specific alignment function
			                     substitutionMatrix = "BLOSUM62", 
			                     gapOpening = -2,
			                     gapExtension = -8, 
			                     scoreOnly = FALSE)

						iso1 <- as.vector(stringr::str_split_fixed(as.character(pwalign::pattern(xxx)), pattern='', n=nchar(as.character(pwalign::pattern(xxx)))))
						iso2 <- as.vector(stringr::str_split_fixed(as.character(pwalign::subject(xxx)), pattern='', n=nchar(as.character(pwalign::subject(xxx)))))

						gap1 <- length(which(iso1 == '-')) # number of gaps in isoform1
						gap2 <- length(which(iso2 == '-')) # number of gaps in isoform 2
						gap <- gap1+gap2
						count_gap[k] <- gap
						iden[k] <- 100-pwalign::pid(xxx, type='PID3') # to calculate pairwise identity of the Pairwise alignment. 100 is maximum
						if(iden[k] == 0 && gap == 0) {
							bh_transcripts1[k] <- all_transcripts[1]
							cat('Protein', name, ':', k,' of ', length(uniprotids1), ' executed successfully!\n')
							break
						} else {
							exclud_others <- c(exclud_others, name)
							cat('Protein', name, ':', k,' of ', length(uniprotids1), ' has no match between # gaps and identity score...\n')
						}
					}
 				} else {
					exclud_others <- c(exclud_others, name)
					cat('Protein', name,':', k, ' of ', length(uniprotids1), ' has "uni_seq" different from "temp_sum"...\n')
				}
			} else {
				exclud_others <- c(exclud_others, name)
				cat('Protein', name,':', k, ' of ', length(uniprotids1), ' has no cds...\n')
			}
		} 
	countl[k] <- flags # put 1 in the vector position of each processed protein
	} else {
		exclud_others <- c(exclud_others, name)
		cat('Protein', name,':', k, ' of ', length(uniprotids1), ' has no transcripts...\n')
	}
cat('IDS1: Protein',k,' of ', length(uniprotids1), ' done\n')
}

countl1 <- count_gap
idenl1 <- iden

bh_transcripts2 <- rep('',length(uniprotids2))
countl <- c()
count_gap <- rep(0, length(uniprotids2))
iden <- rep(0,length(uniprotids2))

# repeat the same for all prot in uniprotids2

for(k in 1:length(uniprotids2)){
	
	name <- uniprotids2[k]
	res <- tricky_vect %in% uniprot_seqs$name
	if(any(res)){
		excluded1 <- c(excluded1, name)
		cat('Protein', name, ' skipped because its UniProt sequence\n')
		next
	}

	# length of the uniprot seq
	uni_seq <- getLength(uniprot_seqs[uniprotids2[k]])
	alldata <- alld[alld$uniprotswissprot == uniprotids2[k], ]

	# only protein coding
	alldata1 <- alldata[alldata$transcript_biotype == 'protein_coding', ]
	all_transcripts <- alldata1$ensembl_transcript_id

 	if(!identical(all_transcripts, character(0))) {
        temp_seq <- biomaRt::getSequence(id=all_transcripts, type='ensembl_transcript_id', seqType='peptide', mart=ensembl)
	Sys.sleep(4)

 		for(j in 1:length(all_transcripts)) {

			temp_cds <- cds_f[cds_f$transcript_id == all_transcripts[j], ]
			
			if(nrow(temp_cds) > 0) {
				dna_strand <- temp_cds$V7[[1]][1]
				temp_sum <- sum(temp_cds$nt_len)/3
				flags <- 0
				if(temp_sum == uni_seq) {
					flags <- flags+1

					## global alignment
					toalign1 <- temp_seq$peptide[which(temp_seq$ensembl_transcript_id == all_transcripts[j])] ## protein seqeunce of first transcript
					toalign1 <- cleanSeq(toalign1)

					temps <- seqinr::getSequence(genome[temp_cds[[1]][1]])[[1]]
					allseq <- c()
					for(i in 1:length(temp_cds[[1]])) {
						start <- temp_cds$V4[i]
						end <- temp_cds$V5[i]
						tempseq <- temps[start:end]
						if(dna_strand == "-"){
							tempseq <- rev(seqinr::comp(tempseq))
						}
						allseq <- c(allseq, tempseq)
					}
					t <- seqinr::translate(allseq)
                    if(any(tricky_vect %in% t)) {
                        excluded2 <- c(excluded2, name)
 						cat('Protein', name, ' skipped because its biomaRt sequence\n')
                        break
					}
                    else {
					
						toalign2 <- toupper(paste(seqinr::translate(allseq),collapse=''))

						xxx <- pwalign::pairwiseAlignment(toalign1, toalign2,
			                     substitutionMatrix = "BLOSUM62", 
			                     gapOpening = -2,
			                     gapExtension = -8, 
			                     scoreOnly = FALSE)

						iso1 <- as.vector(stringr::str_split_fixed(as.character(pwalign::pattern(xxx)), pattern='', n=nchar(as.character(pwalign::pattern(xxx)))))
						iso2 <- as.vector(stringr::str_split_fixed(as.character(pwalign::subject(xxx)), pattern='', n=nchar(as.character(pwalign::subject(xxx)))))

						gap1 <- length(which(iso1 == '-')) # number of gaps in isoform1
						gap2 <- length(which(iso2 == '-')) # number of gaps in isoform 2
						gap <- gap1+gap2
						count_gap[k] <- gap
						iden[k] <- 100-pwalign::pid(xxx, type='PID3')
						if(iden[k] == 0 && gap == 0){
							bh_transcripts2[k] <- all_transcripts[j]
							cat('Protein', name, ':', k,' of ', length(uniprotids1), ' executed successfully!\n')
							break
						} else {
							exclud_others <- c(exclud_others, name)
							cat('Protein', name, ':', k,' of ', length(uniprotids1), ' has no match between # gaps and identity score...\n')
						}
					}
 				} else {
					exclud_others <- c(exclud_others, name)
					cat('Protein', name,':', k, ' of ', length(uniprotids1), ' has "uni_seq" different from "temp_sum"...\n')
				}
			} else {
				exclud_others <- c(exclud_others, name)
				cat('Protein', name,':', k, ' of ', length(uniprotids1), ' has no cds...\n')
			}
		} 
	countl[k] <- flags # put 1 in the vector position of each processed protein
	} else {
		exclud_others <- c(exclud_others, name)
		cat('Protein', name,':', k, ' of ', length(uniprotids1), ' has no transcripts...\n')
	}
cat('IDS2: Protein',k,' of ', length(uniprotids2), ' done\n')
}


uniprot_pdb$bh_ensembl_transcript_id1 <- bh_transcripts1
uniprot_pdb$bh_ensembl_transcript_id2 <- bh_transcripts2

# consider only those uniprot ids that have a BH transcript match
# Here, the BH transcript match is defined as the exact seqeunce length match of uniprot id to a transcript id

uniprot_pdb_final1 <- uniprot_pdb[uniprot_pdb$bh_ensembl_transcript_id1 != '', ] # subset processed_complex_final with only bh trascripts for PROTEIN1 --> FILE
uniprot_pdb_final2 <- uniprot_pdb_final1[uniprot_pdb_final1$bh_ensembl_transcript_id2 != '', ] # subset FILE with only bh transcripts for PRTOEIN2

fwrite(uniprot_pdb_final2, file.path(base_dir, "EEI-Rattus-Norvegicus/data/uniprot_pdb_Ensembl_final.txt"), sep='\t', quote=FALSE, row.names=FALSE) # save in a file

uniprot_pdb_final <- data.table::fread(file.path(base_dir, "EEI-Rattus-Norvegicus/data/uniprot_pdb_Ensembl_final.txt"), header=TRUE)

### do the exon mapping to uniprot sequences
uniprot_uni_ids1 <- uniprot_pdb_final$PROTEIN1 # extract PROTEIN1 names
transcript_uni_ids1 <- uniprot_pdb_final$bh_ensembl_transcript_id1 # extract TRANSCRIPT name (ENST...)
chain_ids1 <- uniprot_pdb_final$CHAIN1 # extract chains name
pdb_ids <- tolower(uniprot_pdb_final$ID) # extract PDB IDs

uniprot_uni_ids2 <- uniprot_pdb_final$PROTEIN2 # extract PROTEIN2 names
transcript_uni_ids2 <- uniprot_pdb_final$bh_ensembl_transcript_id2 # extract TRANSCRIPT name (ENST...)
chain_ids2 <- uniprot_pdb_final$CHAIN2 # extract chains name

#######################################################################################
# map exons to uniprot sequences #####
# store the mappings
store_dir <- file.path(base_dir, "EEI-Rattus-Norvegicus/data/uniprot_Ensembl_Exon_map")

if(dir.exists(store_dir)){
	unlink(store_dir, recursive=TRUE)
}
dir.create(store_dir)

# for one of the interacting protein
for(k in 1:length(transcript_uni_ids1)){
	# get the uniprot sequence
	uniprot_seq <- getSequence(uniprot_seqs[uniprot_uni_ids1[k]])[[1]] # peptide sequence

	# get all exons for this transcript
	allexons <- exons_f[exons_f$transcript_id == transcript_uni_ids1[k], ] # subset exons_f 

	# get exon numbers from cds...denoting the coding exons for this transcript
	cexon_num <- cds_f[cds_f$transcript_id == transcript_uni_ids1[k], ]
	cexon_num <- cexon_num[order(cexon_num$exon_number), ]

	# get exons
	cexon <- allexons[allexons$exon_number %in% cexon_num$exon_number, ]

	# sort exons by genomic start
	cexon <- cexon[order(cexon$exon_number), ]
	if (dim(cexon)[1] == 0) {
		next } 
	else {
	# for each exon
	exon_entry <- rep('',length(uniprot_seq))
	exon_num_entry <- rep('',length(uniprot_seq))
	start_pos <- 1
	end_pos <- 0
	previous <- 0
	for(j in 1:length(cexon[[1]])){

		temp_pos <- (cexon_num$nt_len[j]+previous)/3 # get the length of the aminoacid sequence produce by each specific exon

		if (temp_pos != 0) {

			temp_pos1 <- floor(temp_pos)
			diff <- temp_pos-temp_pos1

			if(diff == 0){ # integer temp_pos
				end_pos <- end_pos+temp_pos1
				exon_entry[start_pos:end_pos] <- cexon$exon_id[j]
				exon_num_entry[start_pos:end_pos] <- cexon$exon_number[j]
				start_pos <- end_pos+1
				previous <- 0
			}else if(diff < 0.5){
				end_pos <- end_pos+temp_pos1
				exon_entry[start_pos:end_pos] <- cexon$exon_id[j]
				exon_num_entry[start_pos:end_pos] <- cexon$exon_number[j]
				start_pos <- end_pos+1
				previous <- 1
			}else{
				end_pos <- end_pos+temp_pos1+1
				exon_entry[start_pos:end_pos] <- cexon$exon_id[j]
				exon_num_entry[start_pos:end_pos] <- cexon$exon_number[j]
				start_pos <- end_pos+1
				previous <- -1
			}
		}
	}
	if (length(uniprot_seq) == length(exon_entry)) {
	Data <- data.frame(UNIPROT_SEQ_NUM=seq(1,length(uniprot_seq)), UNIPROT=uniprot_seq, EXON=exon_entry, EXON_NUM=exon_num_entry)
	
	fwrite(Data, file.path(store_dir, paste0(uniprot_uni_ids1[k],'_',pdb_ids[k],'_',chain_ids1[k],'.txt')), row.names=FALSE, sep='\t', quote=FALSE)
	cat('Protein',k,' of ', length(transcript_uni_ids1), ' done\n')}
	else { cat('There is a difference in length!\n') }
}

}

# for the other interacting protein
for(k in 1:length(transcript_uni_ids2)){

	# get the uniprot sequence
	uniprot_seq <- getSequence(uniprot_seqs[uniprot_uni_ids2[k]])[[1]]

	# get all exons for this transcript
	allexons <- exons_f[exons_f$transcript_id == transcript_uni_ids2[k], ]

	# get exon numbers from cds...denoting the coding exons for this transcript
	cexon_num <- cds_f[cds_f$transcript_id == transcript_uni_ids2[k], ]
	cexon_num <- cexon_num[order(cexon_num$exon_number), ]

	# get exons
	cexon <- allexons[allexons$exon_number %in% cexon_num$exon_number, ]

	# sort exons by genomic start
	cexon <- cexon[order(cexon$exon_number), ]
	if (dim(cexon)[1] == 0) {
		next } 
	else {
	# for each exon
	exon_entry <- rep('',length(uniprot_seq))
	exon_num_entry <- rep('',length(uniprot_seq))
	start_pos <- 1
	end_pos <- 0
	previous <- 0
	for(j in 1:length(cexon[[1]])){

		temp_pos <- (cexon_num$nt_len[j]+previous)/3

		if(temp_pos != 0){

			temp_pos1 <- floor(temp_pos)
			diff <- temp_pos-temp_pos1

			if(diff == 0){ # integer temp_pos
				end_pos <- end_pos+temp_pos1
				exon_entry[start_pos:end_pos] <- cexon$exon_id[j]
				exon_num_entry[start_pos:end_pos] <- cexon$exon_number[j]
				start_pos <- end_pos+1
				previous <- 0
			}else if(diff < 0.5){
				end_pos <- end_pos+temp_pos1
				exon_entry[start_pos:end_pos] <- cexon$exon_id[j]
				exon_num_entry[start_pos:end_pos] <- cexon$exon_number[j]
				start_pos <- end_pos+1
				previous <- 1
			}else{
				end_pos <- end_pos+temp_pos1+1
				exon_entry[start_pos:end_pos] <- cexon$exon_id[j]
				exon_num_entry[start_pos:end_pos] <- cexon$exon_number[j]
				start_pos <- end_pos+1
				previous <- -1
			}

		}
	
	}

	if (length(uniprot_seq) == length(exon_entry)) {
 	Data <- data.frame(UNIPROT_SEQ_NUM=seq(1,length(uniprot_seq)), UNIPROT=uniprot_seq, EXON=exon_entry, EXON_NUM=exon_num_entry)
 	fwrite(Data, file.path(store_dir, paste0(uniprot_uni_ids2[k],'_',pdb_ids[k],'_',chain_ids2[k],'.txt')), row.names=FALSE, sep='\t', quote=FALSE)
 	cat('Protein',k,' of ', length(transcript_uni_ids2), ' done\n') }
	else { cat('There is a difference in length!\n') }
 }
}
