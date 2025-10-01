##############################################################################################
# Purpose: download interaction files from PISA main script
##############################################################################################

rm(list=ls())

# Check if input file exists
input_file <- '../PISA_data/PDBs_specific.txt'
if (!file.exists(input_file)) {
	stop("Error: Input file '", input_file, "' does not exist!")
}

# Upload of the file
cat("Reading PDB IDs from file...\n")
allpdbs <- try(data.table::fread(input_file, header=FALSE)[[1]])
if (inherits(allpdbs, "try-error")) {
	stop("Error reading the input file. Check if file is properly formatted.")
}

if (length(allpdbs) == 0) {
	stop("No PDB IDs found in the input file!")
}
cat("Total number of PDBs to process:", length(allpdbs), "\n\n")

store_dir <- '../PISA_data/PISA_results_Drosophila'
if(!dir.exists(store_dir)){
	cat("Creating storage directory:", store_dir, "\n")
	dir.create(store_dir, recursive = TRUE)
}

# Create a script that will be used to set up the proper environment in tmux
cat("Creating tmux setup script...\n")
setup_script <- file.path(store_dir, "tmux_setup.sh")
writeLines(
    c("#!/bin/bash",
      "# Script to set up environment before running R script in tmux",
      "",
      "# Define the command to run",
      "CMD=\"$1\"",
      "",
      "# Source bashrc to set up shell environment first",
      "source ~/.bashrc",
      "",
      "# Activate conda environment",
      "source \"$(conda info --base)/etc/profile.d/conda.sh\"",
      "conda activate selenium_env",
      "",
      "# Set library path to conda environment",
      "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH",
      "",
      "# Set temporary directory",
      "export TMPDIR=$HOME/tmp",
      "mkdir -p $TMPDIR",
      "",
      "# Display environment information for debugging",
      "echo \"Running with:\"",
      "echo \"CONDA_PREFIX=$CONDA_PREFIX\"",
      "echo \"LD_LIBRARY_PATH=$LD_LIBRARY_PATH\"",
      "echo \"TMPDIR=$TMPDIR\"",
      "",
      "# Execute the command",
      "eval \"$CMD\""
    ),
    setup_script
)

# Make the script executable
system(paste("chmod +x", setup_script))

start <- 1
toprocess <- length(allpdbs)
block <- 100

batch_number <- 1
cat("Starting PISA processing...\n")

while(toprocess >= 0){

	end <- start + block - 1

	if(end > length(allpdbs)){
		end <- length(allpdbs)
	}

	cat(sprintf("Creating batch %d: Processing PDBs %d to %d\n", batch_number, start, end))
	
	# Create a log directory for each batch
	log_dir <- file.path(store_dir, "logs")
	dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
	
	# Create tmux session name
	session_name <- sprintf("pisa_batch_drosophila_%d", batch_number)
	
	# Create the R command to be run inside the tmux session
	r_cmd <- sprintf('Rscript runpisa.r %d %d %s > %s/batch_%d.log 2>&1', 
					 start, end, store_dir, log_dir, batch_number)
	
	# Use the setup script to create a properly configured tmux session
	cmd <- sprintf('tmux new-session -d -s %s "%s \'%s\'"', 
               session_name, setup_script, r_cmd)
	
	# Execute and check if tmux session was created
	system_result <- system(cmd)
	if (system_result != 0) {
		cat(sprintf("Warning: Failed to create tmux session for batch %d\n", batch_number))
	}
	
	cat(sprintf("Tmux session '%s' created\n", session_name))
	Sys.sleep(2)  # Add small delay between creating sessions

	toprocess <- toprocess-block
	start <- end + 1
	batch_number <- batch_number + 1

}

cat("\nAll tmux sessions created. To view running sessions use: tmux ls\n")
cat("To attach to a specific session use: tmux attach -t pisa_batch_N (where N is the batch number)\n")
cat("To detach from a session: press Ctrl+B then D\n")
cat("Logs can be found in:", log_dir, "\n")

# Check if any tmux sessions were actually created
system("tmux ls")

# Only log essential information
cat(sprintf("Created %d batch(es)\n", batch_number - 1))