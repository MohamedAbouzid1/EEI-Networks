#!/usr/bin/sh

########################################################################
# Purpose: create, for each species, a txt file listing PDBs in which at least one protein of that species is contained
########################################################################

file="./data/database/process_complex_file.txt"

list=$(mktemp)
grep -vE "^ID" "$file" | cut -f1 | sort | uniq > "$list"
awk '{print toupper($0)}' "$list" > "$list.tmp" && mv "$list.tmp" "$list"
#echo $list
#cat $list
find ./data/database/ -mindepth 1 -not -empty -type d -print0| while IFS= read -r -d '' f; do
    echo "Processing folder: $f"
    matches=$(find "$f" -type f -print0 | xargs -0 grep -E "PDB;")
    #echo $x
    while IFS= read -r i; do
        #echo "Match found: $i"
        x=${i##*PDB;}
        y=${x%%;*}
        pdb_code=$(echo "$y" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        #echo "Checking $pdb_code"
        if [[ -n "$pdb_code" ]] && grep -q "^$pdb_code" "$list"; then
            #echo $pdb_code
            echo $pdb_code >> "$f/PDBs_specific.tmp"
        fi
    done <<< "$matches"
    if [[ -f "$f/PDBs_specific.tmp" ]]; then
        sort -u "$f/PDBs_specific.tmp" -o "$f/PDBs_specific.tmp"
        species=${f##*/}
        #echo $s
        awk -v S="$species" '{print $0, S}' "$f/PDBs_specific.tmp" > "$f/PDBs_specific.txt"
    fi
done
