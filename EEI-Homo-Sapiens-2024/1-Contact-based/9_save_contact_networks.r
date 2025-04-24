##############################################################################################
# Purpose: save different contact-based netwroks
##############################################################################################

rm(list=ls())
library(Rcpp)

in_dir <- 'data/CONTACTS'

store_dir <- 'data/CONTACT_networks'
if(!dir.exists(store_dir)){
	dir.create(store_dir)
}

cppFunction("List filtereei(CharacterVector ex1, CharacterVector ex2, CharacterVector exon1, CharacterVector exon2, NumericVector allaa){

    int loop1 = ex1.size();
    int loop2 = exon1.size();
    NumericVector allpositions;

    for(int k=0; k<loop1; k++){

		NumericVector loopx;
		NumericVector tallaa;

        for(int j=0; j<loop2; j++){

            if((ex1[k] == exon1[j]) & (ex2[k] == exon2[j])){
                loopx.push_back(j);
                tallaa.push_back(allaa[j]);
            }

            if((ex1[k] == exon2[j]) & (ex2[k] == exon1[j])){
                loopx.push_back(j);
                tallaa.push_back(allaa[j]);
            }
        }

        int loop3 = tallaa.size();
        int MAX = 0;
        int position = 0;
        for(int j=0; j<loop3; j++){
        	if(tallaa[j] > MAX){
        		MAX = tallaa[j];
        		position = loopx[j];
        	}
            //Rcout << MAX << std::endl;
        }
        //Rcout << loop3 << std::endl;

        allpositions.push_back(position);
    }

    List L = List::create(allpositions);
	return L;
  
}")

cutoff <- c(6)#c(4, 5, 6, 7, 8)
num_aa <- c(1)#c(1,3,5,7,9)

for(k in 1:length(cutoff)){

	temp <- data.table::fread(paste0(in_dir,'/int_exon_pairs',cutoff[k],'.txt'))

	for(j in 1:length(num_aa)){

		wh1 <- which(temp$exon1_coverage >= num_aa[j])
		wh2 <- which(temp$exon2_coverage >= num_aa[j])
		temp1 <- temp[intersect(wh1, wh2), ]
		temp1$allAA <- temp1$exon1_coverage+temp1$exon2_coverage # sum of the coverage of both exons
		gnet <- igraph::as_data_frame(igraph::simplify(igraph::graph_from_data_frame(temp1[,c(3,4)], directed=FALSE))) # get a network
		tempf <- filtereei(gnet[[1]], gnet[[2]], temp1$exon1, temp1$exon2, temp1$allAA) # return the indeces of those pairs that respect the filtering
		ids <- tempf[[1]]+1
		tempff <- temp[ids, ]
		data.table::fwrite(tempff, paste0(store_dir,'/CONTACT_net_',cutoff[k],'_',num_aa[j],'.txt'), sep='\t', row.names=FALSE, quote=FALSE)

	}
}
