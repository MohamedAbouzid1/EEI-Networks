#!/bin/bash
################################################################
# Purpose: Generate global EEI networks based on the three EEI definition approaches detecting which interactions are
# supported by more than one approach. The resulting file collects all the EEI detected by all 3 methods.
################################################################

# Define the path of output folder, create them if not existing
out_dir="./data/global_EEINs"
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi

# Upload the files for each of the three EEINs created using the three approaches
cont_netw="./data/CONTACT_networks/CONTACT_net_6_1.txt"
pisa_netw="./data/PISA_networks_filtered/PISA_EEIN_0.5.txt"
eppic_netw="./data/EPPIC_EEIN_filtered.txt"

# Extract the list of interacting exons for each of the three files

# Contact-based interactions
cont_list=$(mktemp) || { echo "Error: mktemp didn't work properly"; exit 1; }
grep -vE "^prot" "$cont_netw" | cut -f3,4 > "$cont_list"
awk -F'\t' 'NF == 2 && $1 != "" && $2 != "" {print toupper($0)}' "$cont_list" > "$cont_list.tmp" && mv "$cont_list.tmp" "$cont_list"
# Pisa-based interactions
pisa_list=$(mktemp) || { echo "Error: mktemp didn't work properly"; exit 1; }
grep -vE "^exon" "$pisa_netw" | cut -f1,2 > "$pisa_list"
awk -F'\t' 'NF == 2 && $1 != "" && $2 != "" {print toupper($0)}' "$pisa_list" > "$pisa_list.tmp" && mv "$pisa_list.tmp" "$pisa_list"
# EPPIC-based interactions
eppic_list=$(mktemp) || { echo "Errore: mktemp non ha funzionato"; exit 1; }
grep -vE "^exon" "$eppic_netw" | cut -f1,2 > "$eppic_list"
awk '{print toupper($0)}' "$eppic_list" > "$eppic_list.tmp" && mv "$eppic_list.tmp" "$eppic_list"


# INTERSECTION between networks from Contact & PISA

for ex1 in $(cut -f1 $cont_list | sort | uniq)
do
ex1_2=$(grep $ex1 $cont_list | cut -f2)
        ex2_temp=$(grep "$ex1" "$pisa_list" | sed -e 's!'$ex1'!!' | sed 's/^[[:space:]]*//')
        if [[ -n "$ex2_temp" ]]; then
            while IFS=$'\n' read -r i;
            do
            printf -v ex2 '%s' $i
        
            if [[ "$ex1_2" =~ "$ex2" ]]; then
                echo -e "$ex1\t$ex2" >> "$out_dir/match_CP.txt"
            fi
            done <<< $ex2_temp
        fi
done



# INTERSECTION between networks from Contact & EPPIC

for en1 in $(cut -f1 $cont_list | sort | uniq)
do
en1_2=$(grep $en1 $cont_list | cut -f2)
        en2_temp=$(grep "$en1" "$eppic_list" | sed -e 's!'$en1'!!' | sed 's/^[[:space:]]*//')
        if [[ -n "$en2_temp" ]]; then
            while IFS=$'\n' read -r i;
            do
            printf -v en2 '%s' $i
        
            if [[ "$en1_2" =~ "$en2" ]]; then
                echo -e "$en1\t$en2" >> "$out_dir/match_CE.txt"
            fi
            done <<< $en2_temp
        fi
done



# INTERSECTION between networks from PISA & EPPIC

for exon1 in $(cut -f1 $pisa_list | sort | uniq)
do
exon1_2=$(grep $exon1 $pisa_list | cut -f2)
        exon2_temp=$(grep "$exon1" "$eppic_list" | sed -e 's!'$exon1'!!' | sed 's/^[[:space:]]*//')
        if [[ -n "$exon2_temp" ]]; then
            while IFS=$'\n' read -r i;
            do
            printf -v exon2 '%s' $i
        
            if [[ "$exon1_2" =~ "$exon2" ]]; then
                echo -e "$exon1\t$exon2" >> "$out_dir/match_PE.txt"
            fi
            done <<< $exon2_temp
        fi
done


rm -f "$cont_list" "$pisa_list" "$eppic_list"


# Get the interactions supported by all 3 approaches

file1="$out_dir/match_CP.txt"
file2="$out_dir/match_CE.txt"

l1=$(mktemp) || { echo "Error: mktemp didn't work properly"; exit 1; }
cat $file1 > $l1
awk '{print toupper($0)}' "$l1" > "$l1.tmp" && mv "$l1.tmp" "$l1"

l2=$(mktemp) || { echo "Error: mktemp didn't work properly"; exit 1; }
cat $file2 > $l2
awk '{print toupper($0)}' "$l2" > "$l2.tmp" && mv "$l2.tmp" "$l2"


while IFS= read -r line; do
num_fields=$(echo "$line" | awk '{print NF}')
    if [[ "$num_fields" -eq 2 ]]; then # check that each row has exactly 2 exons
        read -r E1 E2 <<< "$line"
        m1=$(grep "$E1" "$l2" | sed -e 's!'$E1'!!'| sed 's/^[[:space:]]*//')
        while IFS=$'\n' read -r i; do
        printf -v m2 '%s' $i
        
        if [[ "$E2" == "$m2" ]]; then
            echo -e "$E1\t$E2" >> high_global_network.txt # a high confidence network with edges supported by all three methods
        fi
    done <<< $m1
    fi

done < "$f1"

rm -f "$l2"
