###################################################################
# Purpose: map homologous exons from Mus musculus EEIs to Homo sapiens homologous exons
# The output file collects pairs of human exons that, based on homology with Mus Musculus, are predicted to interact.
###################################################################
# WARNING: change file paths accordingly
# The script is meant to be run on the high confidence level global networks: those 2, one per species, collecting all the
## the interactions supported by the three approaches. In this case, it's executed on Contact-based networks for both species,
### due to unavailability of PISA- and EPPIC-based networks for Mus Musculus.

rm(list=setdiff(ls()))

library(data.table)

# Upload of the 5 files I need to perform the mapping

# file 1
mus_int <- fread('Mus_musculus/CONTACT_networks/CONTACT_net_6_1.txt', sep='\t') # upload the Contact-based network of Mus musculus

# file 2
mus_gtf <- fread('Mus_musculus/processed_exons.tmp', sep='\t') # upload the exons entries of Mus musculus GTF file
extra_info <- strsplit(mus_gtf$V9, '[;]') # get the column number 9 as a list of lists, one for each entry

# Extract from column number 9, respectively GENE ID and EXON ID
wh1 <- lapply(extra_info, function(x) which(x %like% 'gene_id')) # gene id
gene_id1 <- mapply(function(x, y) x[y], extra_info, wh1)
gene_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(gene_id1),'\\s+'), '[[', 2)))
mus_gtf$gene_id <- gene_id # add it as a new column

wh1 <- lapply(extra_info, function(x) which(x %like% 'exon_id')) # exon id
exon_id1 <- mapply(function(x, y) x[y], extra_info, wh1)
exon_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(exon_id1),'\\s+'), '[[', 2)))
mus_gtf$exon_id <- exon_id # add it as a new column

mus_gtf_f <- mus_gtf[,-9] # remove column number 9, the longest, full of string information. The object is a data.table

# Rename the columns collecting coordinates information
colnames(mus_gtf_f)[1] <- "chr"
colnames(mus_gtf_f)[4] <- "start"
colnames(mus_gtf_f)[5] <- "end"
colnames(mus_gtf_f)[7] <- "strand"

coord_col <- c("chr", "start", "end", "strand")

# file 3
homol_pairs <- fread('EGIO/ExonGroup_testpro_hsa_mus.txt', sep='\t') # upload the output file from EGIO, collecting pairs of homologous exons

# file 4
hsa_gtf_tmp <- fread('data/processed_exons.tmp', sep='\t') # upload the exons entries of Homo Sapiens GTF file
hsa_gtf <- hsa_gtf_tmp[grepl('protein_coding', hsa_gtf_tmp$V9), ]
extra_info <- strsplit(hsa_gtf$V9, '[;]') # get the column number 9 as a list of lists, one for each entry

# file 5
hsa_int <- fread('data/CONTACT_networks/CONTACT_net_6_1.txt', sep='\t') # upload the Contact-based network of Homo Sapiens

# Extract from column number 9, respectively GENE ID and EXON ID
wh1 <- lapply(extra_info, function(x) which(x %like% 'gene_id')) # gene id
gene_id1 <- mapply(function(x, y) x[y], extra_info, wh1)
gene_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(gene_id1),'\\s+'), '[[', 2)))
hsa_gtf$gene_id <- gene_id # add it as a new column

wh1 <- lapply(extra_info, function(x) which(x %like% 'exon_id')) # exon id
exon_id1 <- mapply(function(x, y) x[y], extra_info, wh1)
exon_id <- gsub('\\"', "", unlist(lapply(strsplit(trimws(exon_id1),'\\s+'), '[[', 2)))
hsa_gtf$exon_id <- exon_id # add it as a new column

hsa_gtf_f <- hsa_gtf[,-9] # remove column number 9, the longest, full of string information. The object is a data.table

# Rename the columns collecting coordinates information
colnames(hsa_gtf_f)[1] <- "chr"
colnames(hsa_gtf_f)[4] <- "start"
colnames(hsa_gtf_f)[5] <- "end"
colnames(hsa_gtf_f)[7] <- "strand"

# Add a new column with coordinates in the same EGIO format
hsa_gtf_f$coordinates <- paste0("chr", hsa_gtf_f$chr, ":", hsa_gtf_f$start, ":", hsa_gtf_f$end, ":", hsa_gtf_f$strand, "1")
hsa_gtf_f$coordinates <- gsub("\\+1$", "1", hsa_gtf_f$coordinates)

coord_col <- c("chr", "start", "end", "strand")

# Extract from the network the interacting exon pairs in Mus musculus
exon1 <- mus_int$exon1
exon2 <- mus_int$exon2

new_inter <- list()

for(k in 1:dim(mus_int)[1]){ # for each pair of interacting exons in Mus musculus
    ex1 <- exon1[k]
    ex2 <- exon2[k]

    ex1_ind <- which(grepl(ex1, mus_gtf_f$exon_id)) # there coulbe be one or multiple elements in the vector
    ex2_ind <- which(grepl(ex2, mus_gtf_f$exon_id))

    ex1_gtf <- mus_gtf_f[ex1_ind,] # find GTF entries of the exon of interest. We get a subsset of the data.table
    ex2_gtf <- mus_gtf_f[ex2_ind,]

    if(!all(sapply(ex1_gtf[, ..coord_col], uniqueN) == 1)){
        cat("Multiple coordinates annotated for the same exon1!\n")
    } else if(!all(sapply(ex2_gtf[, ..coord_col], uniqueN) == 1)) {
        cat("Multiple coordinates annotated for the same exon2!\n")
    } else {
        cat("Unique exon coordinates\n")
        ex1_coord_tmp <- paste0("chr", ex1_gtf[1, chr], ":", ex1_gtf[1, start], ":", ex1_gtf[1, end], ":", ex1_gtf[1, strand], "1") # create the string of exon coordinates, in the same format of EGIO file
        ex1_coord <- gsub("\\+1$", "1", ex1_coord_tmp)
        ex2_coord_tmp <- paste0("chr", ex2_gtf[1, chr], ":", ex2_gtf[1, start], ":", ex2_gtf[1, end], ":", ex2_gtf[1, strand], "1")
        ex2_coord <- gsub("\\+1$", "1", ex2_coord_tmp)
    
        coor1_ind <- which(grepl(ex1_coord, homol_pairs$musPos)) # there coulbe be one or multiple elements in the vector
        coor2_ind <- which(grepl(ex2_coord, homol_pairs$musPos))
        if(length(coor1_ind) != 0 && length(coor2_ind) != 0) {
            cat("Homologous found!\n")
            homol_coor1 <- homol_pairs$hsaPos[coor1_ind]
            homol_coor2 <- homol_pairs$hsaPos[coor2_ind]

            homol_ex1_ind <- which(grepl(homol_coor1, hsa_gtf_f$coordinates)) # there coulbe be one or multiple elements in the vector
            homol_ex2_ind <- which(grepl(homol_coor2, hsa_gtf_f$coordinates))

            homol_ex1 <- hsa_gtf_f[homol_ex1_ind,] # subset of Homo Sapiens gtf file
            homol_ex2 <- hsa_gtf_f[homol_ex2_ind,]

            if(!all(sapply(homol_ex1[, exon_id], uniqueN) == 1)) {
                cat("Multiple exon1 names for same coordinates!\n")
                } else if (!all(sapply(homol_ex2[, exon_id], uniqueN) == 1)) {
                    cat("Multiple exon2 names for same coordinates!\n")
                } else {
                homol_pair1 <- c(homol_ex1[1, exon_id], ex1_gtf[1, exon_id]) # create vectors of homologous pairs
                homol_pair2 <- c(homol_ex2[1, exon_id], ex2_gtf[1, exon_id])
                    
                if(!is.na(homol_pair1[1]) && !is.na(homol_pair1[1])) {
                    resA <- hsa_int[(grepl(homol_pair1[1], hsa_int$exon1) & grepl(homol_pair2[1], hsa_int$exon2)), ]
                    resB <- hsa_int[grepl(homol_pair1[1], hsa_int$exon2) & grepl(homol_pair2[1], hsa_int$exon1), ]
                        
                    res <- merge(resA, resB)
                    if(nrow(res) != 0) {
                        cat("Interaction already detected in human...\n")
                    } else {
                    cat("NEW INTERACTION!!!\n")
                    new_inter <- append(new_inter, list(c(homol_ex1[1, exon_id], homol_ex2[1, exon_id])))
                    }
                }   
            }
        } else {
        cat("The exon has NO homologous\n")
        } 
    }
    cat("Exon pair ", k, " out of ", dim(mus_int)[1]," done.\n")
}

# Save new interactions in a dataframe and then in a txt file
results <- do.call(rbind, new_inter)
colnames(results) <- c("exon1", "exon2")

write.table(results, file = "new_interactions.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
