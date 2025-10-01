##############################################################################################
# Purpose: download interaction files from PISA
##############################################################################################

args=commandArgs(TRUE)

if (length(args) != 3) {
    stop("Error: Required 3 arguments: start_index end_index store_dir")
}

cat("Starting runpisa.r with arguments:", paste(args, collapse=" "), "\n")

# Check if input file exists
input_file <- '../../data/pdbs_redownload.txt'
if (!file.exists(input_file)) {
    stop("Error: Input file '", input_file, "' does not exist!")
}

allpdbs <- try(data.table::fread(input_file, header=FALSE)[[1]])
if (inherits(allpdbs, "try-error")) {
    stop("Error reading the PDB IDs file")
}

start <- as.numeric(args[1])
end <- as.numeric(args[2])
store_dir <- args[3]

if (is.na(start) || is.na(end)) {
    stop("Error: start and end must be numbers")
}

if (start > length(allpdbs) || end > length(allpdbs)) {
    stop("Error: start or end index exceeds number of available PDB IDs")
}

cat(sprintf("Processing PDBs from index %d to %d\n", start, end))

# Create a wrapper script to ensure environment variables are properly set
wrapper_script <- "run_with_env.sh"
cat("Creating environment wrapper script...\n")

# Write the wrapper script
writeLines(
    c("#!/bin/bash",
      "# Wrapper script to set environment variables before running the Python script",
      "",
      "# Source bashrc to set up the shell environment first",
      "source ~/.bashrc",
      "",
      "# Activate conda environment if available",
      "if [ -f \"$(conda info --base)/etc/profile.d/conda.sh\" ]; then",
      "    source \"$(conda info --base)/etc/profile.d/conda.sh\"",
      "    conda activate selenium_env",
      "fi",
      "",
      "# Set library path to conda environment",
      "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH",
      "",
      "# Set temporary directory",
      "export TMPDIR=$HOME/tmp",
      "mkdir -p $TMPDIR",
      "",
      "# Force enable Wayland support",
      "export MOZ_ENABLE_WAYLAND=1",
      "",
      "# Run the Python script with passed arguments",
      "python \"$@\""
    ),
    wrapper_script
)

# Make the wrapper script executable
system(paste("chmod +x", wrapper_script))

for (k in start:end) {
    cat(sprintf("Processing PDB %d: %s\n", k, allpdbs[k]))
    
    if (allpdbs[k] == "1jbj") {
        cat("Skipping 1jbj due to known error\n")
        next
    }
    
    pdb_dir <- paste0(store_dir, '/', allpdbs[k])
    if (!dir.exists(pdb_dir)) {
        # Use the wrapper script to run the Python script
        tmpcmd <- paste('./run_with_env.sh PisaAuto_id_marina.py', allpdbs[k], store_dir)
        cat("Running command:", tmpcmd, "\n")
        result <- try(system(tmpcmd))
        if (inherits(result, "try-error") || result != 0) {
            cat("Error processing PDB:", allpdbs[k], "\n")
        }
    } else {
        cat(allpdbs[k], 'already present\n')
    }
}