###################################################################
# Purpose: get the name of the species for which Protein-Protein Interactions have been detected
## in PPI-Eukaryotes analysis
##################################################################

# select all UniprotIDs in the Contact-based PPI network
grep -vE "protein" ./data/database/CONTACT_networks/CONTACT_net_6_1.txt | cut -f1 | sed 's!_.*$!!g' > ./data/database/CONTACT_networks/p1.tmp
grep -vE "protein" ./data/database/CONTACT_networks/CONTACT_net_6_1.txt | cut -f2 | sed 's!_.*$!!g' > ./data/database/CONTACT_networks/p2.tmp

cat p1.tmp p2.tmp | sort -u > ./data/database/CONTACT_uniq_prot.txt # 7 823 proteins

rm ./data/database/CONTACT_networks/*.tmp

list="./data/database/CONTACT_uniq_prot.txt"

# Get, for each protein, the species to which it belongs
find ./data/database/ -mindepth 1 -not -empty -type d -print0| while IFS= read -r -d '' f; do
    echo "Processing folder: $f"
    for file in "$f"/*Reviewed.txt; 
        do
        sec_line=$(sed -n '2p' "$file")
        #echo "whole line is " $sec_line
        names=$(echo "$sec_line" | cut -d " " -f4- | tr ';' '\n') 
        #echo $names
        while IFS= read -r name; do
            if grep -Fxq "$name" "$list"; then
                echo -e "$name\t${f##*/}" >> ./data/database/CONTACT_prot_species_pairs.txt
            fi
        done <<< "$names"
    done
done
