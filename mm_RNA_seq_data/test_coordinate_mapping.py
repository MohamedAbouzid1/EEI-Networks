#!/usr/bin/env python3
"""
Test script to debug coordinate mapping issues.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'utils'))
from utils.coord_gene_mapping import coordinate_to_gene_mapping_with_files
import pandas as pd

def test_coordinate_mapping():
    """Test coordinate mapping with sample coordinates."""
    
    # Test coordinates from the EEIs file
    test_coordinates = [
        "chr4:148568767:148568879:1",
        "chr15:6803290:6803455:1", 
        "chr11:105987817:105987938:1",
        "chr11:6544148:6544286:1"
    ]
    
    mapping_dir = "coord_gene_mapping_files"
    
    print("=== TESTING COORDINATE MAPPING ===")
    
    for coord in test_coordinates:
        print(f"\nTesting coordinate: {coord}")
        
        try:
            mapped_gene = coordinate_to_gene_mapping_with_files(coord, mapping_dir)
            print(f"  Mapped to: {mapped_gene}")
            
            if mapped_gene:
                # Check if the gene exists in expression data
                expression_file = "geo_data/mouse_expression_transposed.tsv"
                expression_df = pd.read_csv(expression_file, sep='\t', index_col=0)
                
                if mapped_gene in expression_df.index:
                    print(f"  ✓ Gene '{mapped_gene}' found in expression data")
                else:
                    print(f"  ✗ Gene '{mapped_gene}' NOT found in expression data")
            else:
                print("  ✗ No mapping found")
                
        except Exception as e:
            print(f"  Error: {e}")
    
    print("\n=== CHECKING MAPPING FILES ===")
    
    # Check if mapping files exist
    combined_file = os.path.join(mapping_dir, "combined_mapping.json")
    if os.path.exists(combined_file):
        print(f"✓ Combined mapping file exists: {combined_file}")
        
        # Load and check the mapping
        import json
        with open(combined_file, 'r') as f:
            mappings = json.load(f)
        
        coord_mapping = mappings['coordinate_to_gene']
        available_genes = mappings['available_genes']
        
        print(f"  Coordinate mappings: {len(coord_mapping)}")
        print(f"  Available genes: {len(available_genes)}")
        
        # Check for chromosome 4 coordinates
        chr4_coords = [k for k in coord_mapping.keys() if k.startswith('4:')]
        print(f"  Chromosome 4 coordinates: {len(chr4_coords)}")
        
        if chr4_coords:
            print(f"  Sample chr4 coordinates: {chr4_coords[:3]}")
            
    else:
        print(f"✗ Combined mapping file not found: {combined_file}")

if __name__ == "__main__":
    test_coordinate_mapping() 