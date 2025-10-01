#!/usr/bin/env python3
"""
Script to create high confidence EEI network by finding EEIs that exist 
in all three networks: Contact, PISA, and EPPIC.
"""

import pandas as pd
import sys
from pathlib import Path
from data_loader import load_contact_network, load_eppic_network, load_pisa_network
from create_nets import find_high_confidence_eeis, create_high_confidence_network


def main():
    """Main function to create high-confidence EEI network."""
    
    # File paths - update these to match your file locations
    contact_file = '../EEI-Mus-Musculus/results/Mus_musculus/CONTACT_networks/CONTACT_net_6_1.txt'
    eppic_file = '../EEI-Mus-Musculus/3-EPPIC-based/results/EPPIC_EEIN_filtered.txt'
    pisa_file = '../EEI-Mus-Musculus/PISA-based/PISA_networks_filtered/PISA_EEIN_0.5.txt'
    output_file = '../EEI-Mus-Musculus/high_confidence_network/high_confidence_eei_network_mm.txt'
    
    # Check if files exist
    files = [contact_file, eppic_file, pisa_file]
    for f in files:
        if not Path(f).exists():
            print(f"Error: File '{f}' not found!")
            print("Please update the file paths in the script.")
            sys.exit(1)
    
    # Load networks
    contact_pairs, contact_pmap, contact_df = load_contact_network(contact_file)
    eppic_pairs, eppic_pmap, eppic_df = load_eppic_network(eppic_file)
    pisa_pairs, pisa_pmap, pisa_df = load_pisa_network(pisa_file)
    
    # Find high-confidence EEIs
    high_conf_pairs = find_high_confidence_eeis(contact_pairs, eppic_pairs, pisa_pairs)
    
    if len(high_conf_pairs) == 0:
        print("\nNo EEIs found in all three networks!")
        print("This might be due to:")
        print("  1. Different exon naming conventions between networks")
        print("  2. Networks covering different protein interactions")
        print("  3. Networks being from different species or databases")
        
        # Show some examples from each network for debugging
        print("\nExample EEIs from each network:")
        print("Contact:", list(contact_pairs)[:3] if contact_pairs else "Empty")
        print("EPPIC:", list(eppic_pairs)[:3] if eppic_pairs else "Empty")
        print("PISA:", list(pisa_pairs)[:3] if pisa_pairs else "Empty")
        
        return
    
    # Create comprehensive high-confidence network
    print("\nCreating comprehensive high-confidence network...")
    high_conf_df = create_high_confidence_network(
        high_conf_pairs, contact_df, eppic_df, pisa_df,
        contact_pmap, eppic_pmap, pisa_pmap
    )
    
    # Sort by proteins and exons for better organization
    high_conf_df = high_conf_df.sort_values(['protein1', 'protein2', 'exon1', 'exon2'])
    
    # Save to file
    high_conf_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nHigh-confidence network saved to: {output_file}")
    
    # Print summary statistics
    print(f"\nSummary:")
    print(f"  Total high-confidence EEIs: {len(high_conf_df)}")
    print(f"  Unique protein pairs: {high_conf_df[['protein1', 'protein2']].drop_duplicates().shape[0]}")
    
    # Show first few rows
    print(f"\nFirst 5 high-confidence EEIs:")
    print(high_conf_df.head().to_string())

if __name__ == "__main__":
    main()