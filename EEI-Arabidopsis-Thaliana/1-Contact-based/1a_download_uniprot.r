##############################################################################################
# Purpose: generate uniprot files 
##############################################################################################

rm(list=ls()) # to clear the currentworkspace, removing all objects and variables

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Set C++ standard to C++11
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# Installation of some missing packages - ONLY ONCE
#install.packages("seqinr") #4.2.36
#install.packages("data.table") #1.16.2
#install.packages("stringr") #1.5.1
#install.packages("plyr") #1.8.9
##
library(seqinr)
library(data.table)
library(Rcpp) #1.0.13.1
library(stringr)
library(plyr)

# turning warnings into errors
options(warn=2)

# IF the UniProt file has been already downloaded, leave the following section of lines as comments.
##################################
## create directory to store data
#if (!dir.exists('/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/data')) dir.create('/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/data')
#data_dir <- '/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/data'

## Proteins from Uniprot: complete UniProtKB/Swiss-Prot data set in flat file format which contains fully annotated curated entries
## It's eight-weekly updated: "UniProtKB/Swiss-Prot Release 2024_05" of 02-Oct-2024
#if (!file.exists(paste0(data_dir, "/uniprot_sprot.dat"))) {
  #system(paste0("wget -O ", data_dir, "/uniprot_sprot.dat.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"))
  #system(paste0("gunzip --force ", data_dir, "/uniprot_sprot.dat.gz"))
#}

## save separate files for uniprot entries
uniprot_dir <- '/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/data/uniprot'
#if (!dir.exists(uniprot_dir)) dir.create(uniprot_dir)
## create one file per protein ID
#system(paste0("awk '/^ID/{if(NR!=1){for(i=0;i<j;i++)print a[i]>\"", data_dir, "/uniprot/file\"k;j=0;k++;}a[j++]=$0;next}{a[j++]=$0;}END{for(i=0;i<j;i++)print a[i]>\"", data_dir, "/uniprot/file\"k}' i=0 k=1  ", data_dir, "/uniprot_sprot.dat"))

##################################

# create directory to store data
if (!dir.exists('../data')) dir.create('../data')

# filter to only retain Arabidopsis thaliana
output_dir <- '../data/uniprot_sprot_AT'
if (!dir.exists(output_dir)) dir.create(output_dir)

# Get list of all files in uniprot directory
allfiles <- list.files(uniprot_dir, full.names=TRUE)

# check for every file if they are Arabidopsis thaliana-related protein. If yes, copy that file in uniprot_sprot_AT directory, properly changing the file name.
counter <- 0
for(k in 1:length(allfiles)){

	# Add error handling for readLines with incomplete final lines
	tryCatch({
		tempf <- readLines(allfiles[k], warn = FALSE)
	}, error = function(e) {
		cat("Error reading file", allfiles[k], ":", e$message, "\n")
		next
	})
	
	# Skip if tempf is NULL or empty
	if(is.null(tempf) || length(tempf) == 0) {
		cat("Skipping file", allfiles[k], "- empty or unreadable\n")
		next
	}

	for(j in 1:length(tempf)){

		temp <- tempf[j]
		temp2 <- substr(temp, 1,2)

		# Check if temp is not empty and has content before processing
		if(length(temp) > 0 && nchar(temp) > 0 && (temp2 == 'OS') && (temp %like% 'Arabidopsis thaliana')){
			tf <- strsplit(tempf[1],'\\s+')
			filename <- paste0(output_dir,'/',tf[[1]][2],'_',strsplit(tf[[1]][3],';')[[1]][1],'.txt')
			if (!file.exists(filename)) file.copy(allfiles[k], filename)
			counter <- counter+1                            
			break
		}

	}

	cat(k,' of ', length(allfiles), ' done\n')

}
