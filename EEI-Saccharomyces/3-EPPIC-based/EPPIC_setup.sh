#!/usr/bin/bash

# Create a directory for software and databases
INSTALL_DIR="/cosybio_project/mabouzid/EEI_networks/software"
DATABASE_DIR="/cosybio_project/mabouzid/EEI_networks/databases"

# Create necessary directories
mkdir -p $INSTALL_DIR
mkdir -p $DATABASE_DIR

## download the git source
cd $INSTALL_DIR
git clone https://github.com/eppic-team/eppic

## setup the databases
mkdir -p $DATABASE_DIR/Uniprot100
cd $DATABASE_DIR/Uniprot100
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
gunzip -c uniref100.fasta.gz | makeblastdb -dbtype prot -logfile makeblastdb.log -parse_seqids -out uniref100.fasta -title uniref100.fasta

## mmseq2 software via conda
# Option 1: If conda is installed
# conda install -c conda-forge -c bioconda mmseqs2

# Option 2: Direct download (recommended for more control)
cd $INSTALL_DIR
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
tar xvfz mmseqs-linux-avx2.tar.gz
MMSEQS_BIN="$INSTALL_DIR/mmseqs/bin/mmseqs"

## clustalomega installation
cd $INSTALL_DIR
# First install argtable2 (dependency)
wget http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz
tar -xf argtable2-13.tar.gz
cd argtable2-13
./configure --prefix=$INSTALL_DIR/argtable2-13
make
make install

# Now install clustal omega
cd $INSTALL_DIR
wget http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz
tar -xf clustal-omega-1.2.4.tar.gz
cd clustal-omega-1.2.4
./configure CFLAGS="-I$INSTALL_DIR/argtable2-13/include" LDFLAGS="-L$INSTALL_DIR/argtable2-13/lib"
make
make install prefix=$INSTALL_DIR/clustal-omega-1.2.4
CLUSTALO_BIN="$INSTALL_DIR/clustal-omega-1.2.4/bin/clustalo"

## EPPIC software
cd $INSTALL_DIR
wget http://eppic-web.org/downloads/eppic.zip
unzip eppic.zip

# Get the path to blastp
BLASTP_BIN=$(which blastp)
if [ -z "$BLASTP_BIN" ]; then
    echo "BLASTP not found in PATH. Please install BLAST+ tools."
    exit 1
fi

# Create the configuration file
cat > ~/.eppic.conf << EOF
BLAST_DB_DIR=$DATABASE_DIR/Uniprot100
BLAST_DB=uniref100.fasta
MMSEQS_BIN=$MMSEQS_BIN
BLASTP_BIN=$BLASTP_BIN
CLUSTALO_BIN=$CLUSTALO_BIN
EOF

echo "EPPIC setup complete. Configuration file created at ~/.eppic.conf"