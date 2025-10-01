#!/bin/bash

# Define the path of working folders, create them if not existing
sql_dir="./SQLtables"
if [ ! -d "$sql_dir" ]; then
    mkdir -p "$sql_dir"
fi


# Define TaxID and short-name of the 2 species of interest: always Homo Sapiens vs another organism
specie1="hsa"
full_name1="homo_sapiens"
taxIdA=9606 # Homo sapiens
sp1_dir="./SQLtables/$specie1"
if [ ! -d "$sp1_dir" ]; then
    mkdir -p "$sp1_dir"
fi

specie2="mus"
full_name2="mus_musculus"
taxIdB=10090 # Mus musculus
sp2_dir="./SQLtables/$specie2"
if [ ! -d "$sp2_dir" ]; then
    mkdir -p "$sp2_dir"
fi

#  Download the SQLtable for homologous genes between the two species

url="https://inparanoidb.sbc.su.se/download/sqltable/${taxIdA}&${taxIdB}&prot"

if [ ! -f "$sp2_dir/SQLtable.$taxIdA.fa-$taxIdB.fa" ]; then
    wget -O "$sp2_dir/SQLtable.$taxIdA.fa-$taxIdB.fa" "$url"
fi


# Extract UniProtID from the SQLtable downloaded -- the order of rows is fundamental, not change it!

sqltable="$sp2_dir/SQLtable.$taxIdA.fa-$taxIdB.fa"

cut -f5 $sqltable | sed -n '1~2p' > $sp2_dir/"$specie2"_uniprotid.txt # extraction of odd rows that belong to specie2
cut -f5 $sqltable | sed -n 'n;p' > $sp2_dir/"$specie1"_uniprotid.txt # extraction of even rows that belong to specie1


# Run use_biomart.r script to get, from the UniProtID, all the Ensembl ID that are used by EGIO. Run it twice, once per organism
sp1_ensembl_dataset="hsapiens_gene_ensembl"
Rscript use_biomart.r "$specie1" "$sp2_dir" "$sp1_ensembl_dataset"

## look for the name of the ensembl dataset and modify it 
sp2_ensembl_dataset="mmusculus_gene_ensembl"  
Rscript use_biomart.r "$specie2" "$sp2_dir" "$sp2_ensembl_dataset"


# Find the Ensembl ID for each of the genes part of the homologous pairs

for s in $specie1 $specie2
do
#echo $s
    uniprot_ids="$sp2_dir/"$s"_uniprotid.txt"
    ensembl_ids="$sp2_dir/"$s"_ensemblid.txt"
    output_file="$sp2_dir/"$s"_homol.txt"

    while IFS= read -r l; do
        #echo "$l"
        ensembl_matching=$(grep -m1 "$l" "$ensembl_ids" | cut -f2)
        if [ -n "$ensembl_matching" ]; then
            echo -e "$l\t$ensembl_matching" >> "$output_file"
        else
            echo -e "$l\tNo match" >> "$output_file"
        fi
    done < "$uniprot_ids"
done

# Merge in one file, each column collecting homologous genes of one of the species
paste $sp2_dir/"$specie1"_homol.txt $sp2_dir/"$specie2"_homol.txt > $sp2_dir/homologous_"$specie1"_"$specie2".tmp

# Remove lines having at least one "No match" between UniProt and Ensembl ID
grep -vE "No match" $sp2_dir/homologous_"$specie1"_"$specie2".tmp | cut -f2,4 > $sp2_dir/homologous_"$specie1"_"$specie2".txt

# Remove tmp file
rm $sp2_dir/homologous_"$specie1"_"$specie2".tmp 

# Add as first row the species' names
sed -i "1i $specie1\t$specie2" $sp2_dir/homologous_"$specie1"_"$specie2".txt

homol_file="$sp2_dir/homologous_"$specie1"_"$specie2".txt"


# Download required files from Ensembl, if not already present, and unzip

for S in $full_name1 $full_name2
do
    if [ "$S" = "$full_name1" ]; then
        sp_dir="$sp1_dir"
    else [ "$S" = "$full_name2" ]
        sp_dir="$sp2_dir"
    fi
    
    # cDNA fasta
    if ! ls "$sp_dir"/*cdna* 1> /dev/null 2>&1; then
        rsync -av --progress rsync://ftp.ensembl.org/ensembl/pub/release-113/fasta/$S/cdna/*all.fa.gz "$sp_dir/"
        gunzip -f "$sp_dir"/*cdna*.gz
    else
        echo "cDNA fasta already present"
    fi

    # CDS fasta
    if ! ls "$sp_dir"/*cds* 1> /dev/null 2>&1; then
        rsync -av --progress rsync://ftp.ensembl.org/ensembl/pub/release-113/fasta/$S/cds/*all.fa.gz "$sp_dir/"
        gunzip -f "$sp_dir"/*cds*.gz
    else
        echo "CDS fasta already present"
    fi

    # reference transcriptome gtf 
    if ! ls "$sp_dir"/*gtf* 1> /dev/null 2>&1; then
        rsync -av --progress rsync://ftp.ensembl.org/ensembl/pub/release-113/gtf/$S/*113.gtf.gz "$sp_dir/"
        gunzip -f "$sp_dir"/*gtf*.gz    
    else
        echo "gtf already present"
    fi
done

# Assign files easy names
cdna1_file=$(ls "$sp1_dir"/*cdna*.fa)
cds1_file=$(ls "$sp1_dir"/*cds*.fa)
gtf1_file=$(ls "$sp1_dir"/*.gtf)

cdna2_file=$(ls "$sp2_dir"/*cdna*.fa)
cds2_file=$(ls "$sp2_dir"/*cds*.fa)
gtf2_file=$(ls "$sp2_dir"/*.gtf)


# Run EGIO
cd ./EGIO
./_RUN_egio.sh -s "$specie1" -S "$specie2" -r ../"$gtf1_file" -R ../"$gtf2_file" -e ../"$cdna1_file" -E ../"$cdna2_file" -o ../"$cds1_file" -O ../"$cds2_file" -h ../"$homol_file" -p 6 -i 0.8 -c 0.8 -m 2 -n -2 -g -1
