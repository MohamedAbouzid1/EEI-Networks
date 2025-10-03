#!/usr/bin/env python3
"""
Debug script to test coordinate-to-gene mapping for significant EEIs
"""

import json
import os

def test_coordinate_mapping():
    """Test the coordinate mapping for the 4 significant EEIs"""
    
    # The 4 significant EEI coordinates from your results
    test_coordinates = [
        'chr15:101618188:101618313:-1',
        'chr11:100094934:100095154:-1', 
        'chr15:101617311:101617531:-1',
        'chr11:100095250:100095375:-1'
    ]
    
    # Load the mapping file
    mapping_file = 'coord_gene_mapping_files/combined_mapping.json'
    if not os.path.exists(mapping_file):
        print(f"Error: Mapping file not found: {mapping_file}")
        return
    
    with open(mapping_file, 'r') as f:
        mappings = json.load(f)
    
    coord_mapping = mappings['coordinate_to_gene']
    available_genes = mappings['available_genes']
    
    print("=== COORDINATE MAPPING DEBUG ===")
    print(f"Total coordinates in mapping: {len(coord_mapping)}")
    print(f"Available genes: {len(available_genes)}")
    print()
    
    # Check what chromosomes are available
    available_chroms = set()
    for coord_key, gene_info in coord_mapping.items():
        available_chroms.add(gene_info['chrom'])
    
    print(f"Available chromosomes: {sorted(available_chroms)}")
    print()
    
    # Test each coordinate
    for i, coord in enumerate(test_coordinates, 1):
        print(f"EEI {i}: {coord}")
        
        # Parse the coordinate
        parts = coord.split(':')
        chrom = parts[0]  # chr15, chr11
        start = int(parts[1])
        end = int(parts[2])
        strand = parts[3]
        
        print(f"  Parsed: chrom={chrom}, start={start}, end={end}, strand={strand}")
        
        # Convert to mapping file format
        chrom_num = chrom.replace('chr', '')
        mapping_key = f"{chrom_num}:{start}-{end}"
        
        print(f"  Mapping key: {mapping_key}")
        
        # Check if exact match exists
        if mapping_key in coord_mapping:
            gene_info = coord_mapping[mapping_key]
            print(f"  ✓ EXACT MATCH: {gene_info['gene_name']}")
        else:
            print(f"  ✗ No exact match")
            
            # Check for overlapping coordinates
            overlapping_genes = []
            for coord_key, gene_info in coord_mapping.items():
                if gene_info['chrom'] == chrom_num:
                    coord_start = gene_info['start']
                    coord_end = gene_info['end']
                    
                    # Check for overlap
                    if not (end < coord_start or start > coord_end):
                        overlapping_genes.append(gene_info['gene_name'])
            
            if overlapping_genes:
                print(f"  ✓ OVERLAPPING GENES: {overlapping_genes}")
            else:
                print(f"  ✗ No overlapping genes found")
        
        print()

if __name__ == "__main__":
    test_coordinate_mapping()
