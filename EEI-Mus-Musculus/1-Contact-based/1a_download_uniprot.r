##############################################################################################
# Purpose: generate uniprot files for Mus musculus only
##############################################################################################

rm(list=ls()) # to clear the currentworkspace, removing all objects and variables

# Installation of some missing packages - ONLY ONCE
install.packages("seqinr") #4.2.36
install.packages("data.table") #1.16.2
install.packages("stringr") #1.5.1
install.packages("plyr") #1.8.9
##
library(seqinr)
library(data.table)
library(Rcpp) #1.0.13.1
library(stringr)
library(plyr)

# turning warnings into errors
options(warn=2)

# IF the UniProt file has been already download for Homo Sapiens, leave the following section of lines as comments.
##################################
## create directory to store data
# dir.create('data')

## Proteins from Uniprot: complete UniProtKB/Swiss-Prot data set in flat file format which contains fully annotated curated entries
## It's eight-weekly updated: "UniProtKB/Swiss-Prot Release 2024_05" of 02-Oct-2024
# system("wget -O data/uniprot_sprot.dat.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz")
# system("gunzip --force data/uniprot_sprot.dat.gz")

## save separate files for uniprot entries
# uniprot_dir <- 'data/uniprot'
# dir.create(uniprot_dir)
## create one file per protein ID
# system("awk '/^ID/{if(NR!=1){for(i=0;i<j;i++)print a[i]>\"data/uniprot/file\"k;j=0;k++;}a[j++]=$0;next}{a[j++]=$0;}END{for(i=0;i<j;i++)print a[i]>\"data/uniprot/file\"k}' i=0 k=1  data/uniprot_sprot.dat")

##################################

# create directory to store data
dir.create('Mus_musculus')

# filter to only retain Mus musculus
allfiles <- list.files(uniprot_dir, full.names=TRUE)
output_dir <- '../Mus_musculus/uniprot_sprot_MM'
dir.create(output_dir)
# check for every file if they are Homo Sapiens-related protein. If yes, copy that file in uniprot_sprot_HS directory, properly changing the file name.
counter <- 0
for(k in 1:length(allfiles)){

	tempf <- readLines(allfiles[k])

	for(j in 1:length(tempf)){

		temp <- tempf[j]
		temp2 <- substr(temp, 1,2)

		if((temp2 == 'OS') & (temp %like% 'Mus musculus')){
			tf <- strsplit(tempf[1],'\\s+')
			filename <- paste0(output_dir,'/',tf[[1]][2],'_',strsplit(tf[[1]][3],';')[[1]][1],'.txt')
			file.copy(allfiles[k], filename)
			counter <- counter+1                            
			break
		}

	}

	cat(k,' of ', length(allfiles), ' done\n')

}
