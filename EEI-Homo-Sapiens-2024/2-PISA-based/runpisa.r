.##############################################################################################
# Purpose: download interaction files from PISA
##############################################################################################

args=commandArgs(TRUE)

allpdbs <- data.table::fread('../data/allPDB_IDs.txt', header=FALSE)[[1]] # file created in script 1-Contact-based/6_all_map.r

start <- args[1]
end <- args[2]
store_dir <- args[3]

for (k in start:end) {
    if ( allpdbs[k] == "1jbj" ) { # it was generating an error
        next
    }

    if (!dir.exists(paste0(store_dir, '/', allpdbs[k]))) {
        tmpcmd <- paste('python PisaAuto_id.py', allpdbs[k], store_dir)
        system(tmpcmd)
    } else {
        cat(allpdbs[k], 'already present\n')
    }
}
