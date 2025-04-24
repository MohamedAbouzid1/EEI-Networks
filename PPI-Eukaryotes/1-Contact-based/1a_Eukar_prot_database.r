#####################################################################
# Purpose: subset UniProtKB/Swiss-prot data set according to the referring species; 
# save in a species-specific folder only proteins belonging to Eukaryotes.
#####################################################################

rm(list=ls()) # get an empty workspace

# Installation of required packages 
# install.packages("seqinr") # 4.2.36
# install.packages("data.table") # 1.16.2
# install.packages("stringr") # 1.5.1
# install.packages("plyr") # 1.8.9
##
# import libraries
library(seqinr)
library(data.table)
library(Rcpp)
library(stringr)
library(plyr)

# Turning warnings into errors
options(warn=2)

# Proteins from Uniprot: complete UniProtKB/Swiss-Prot data set in flat file format which contains fully annotated curated entries
# It's eight-weekly updated: "UniProtKB/Swiss-Prot Release 2024_05" of 02-Oct-2024

system("wget -O data/uniprot_sprot.dat.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz")
system("gunzip --force data/uniprot_sprot.dat.gz")

# Save separate files for uniprot entries
uniprot_dir <- 'data/uniprot' # input directory
db_dir <- 'data/database' # output directory to be created, if not present
if(!dir.exists(db_dir)) {
	dir.create(db_dir)} else
	{cat('Database folder already exists!\n')}

# Create one file per protein ID
system("awk '/^ID/{if(NR!=1){for(i=0;i<j;i++)print a[i]>"data/uniprot/file"k;j=0;k++;}a[j++]=$0;next}{a[j++]=$0;}END{for(i=0;i<j;i++)print a[i]>"data/uniprot/file"k}' i=0 k=1  data/uniprot_sprot.dat")

# Filter according to the belonging species

# List of all files
allfiles <- list.files(uniprot_dir, full.names=TRUE) 

# List of all recorded species, 12,094 different unique species name (Details within parentheses in species names are ignored for grouping)
system("awk '/^OS/ {if (NR>ingroup) print; ingroup=NR+1}' data/uniprot_sprot.dat | sed 's!^OS[[:space:]]*!!g' | sed 's! (.*$!!' | sed 's!\\/!_!'| sed 's!\\.\\s*$!!' | sort | uniq > data/database/species_list.txt")

# For each available species, create a subdirectory in db_dir, if it doesn't already exist, with the name of the specie
species <- readLines(file.path(db_dir, "species_list.txt")) # 12,094 

for(i in 1:length(species)){
    species_path <- file.path(db_dir, species[i])
    if(!dir.exists(species_path)){
        dir.create(species_path)
    }
}


# Filter based on eukaryotes and, if yes, copy that file in the corresponding species folder

counter <- 0

for(k in 1:length(allfiles)){ # for each protein-specific file

	tempf <- readLines(allfiles[k]) # read all its lines

	species_index <- which(substr(tempf, 1, 2) == "OS") # index of the line where the species name is collected
	if (length(species_index) > 0 ){
		temp_species_line <- tempf[species_index[1]] # select whole line
		species_line <- gsub("/", "_", temp_species_line) # substitute any / with _
		direc <- rev(which(sapply(species, function(x) grepl(x, species_line))))[[1]] # select the index of the species in species_list

		classification_index <- which(substr(tempf, 1, 2) == "OC" & grepl("Eukaryota", tempf)) # index of the line where the class "Eukaryota" is collected
		if (length(classification_index) > 0 ){
			classification_line <- tempf[classification_index] # select the line of OC

			tf <- strsplit(tempf[1],'\\s+') # select the first line [ID], divided according the spaces

			filename <- paste0(db_dir, '/', species[direc],'/',tf[[1]][2],'_',strsplit(tf[[1]][3],';')[[1]][1],'.txt') #save the name with the 2nd and the 3rd words of the ID line [CODE_SPECIES_Reviewed]
			if(!file.exists(filename)){
       			file.copy(allfiles[k], filename)
   			}
			counter <- counter+1                            
		}
	}
	cat(k,' of ', length(allfiles), ': done\n')
}
