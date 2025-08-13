import pandas as pd
import numpy as np
from collections import defaultdict
import re

def map_human_crpes_to_mouse_orthologs(egio_file, human_crpe_file, 
                                       human_exon_mapping_file,
                                       mouse_exon_mapping_file,
                                       identity_threshold=0.8, 
                                       output_file=None):
    """
    Map human Cancer-Relevant Perturbed EEIs (CRPEs) to their mouse orthologs 
    using EGIO output and exon coordinate mapping files.
    
    Parameters:
    -----------
    egio_file : str
        Path to EGIO output file (ExonGroup_testpro_hsa_mus.txt)
    human_crpe_file : str  
        Path to human CRPE file (e.g., 3_NETHIGH_CRPES.txt)
    human_exon_mapping_file : str
        Path to human exon ID -> coordinate mapping file
    mouse_exon_mapping_file : str
        Path to mouse exon ID -> coordinate mapping file
    identity_threshold : float
        Minimum sequence identity to consider orthology (default: 0.8)
    output_file : str, optional
        Path to save the mapped mouse EEIs
        
    Returns:
    --------
    tuple: (mapped_mouse_eeis, mapping_stats)
        mapped_mouse_eeis : pd.DataFrame
            DataFrame with columns ['mouse_exon1', 'mouse_exon2', 'human_exon1', 'human_exon2', 'identity1', 'identity2']
        mapping_stats : dict
            Dictionary with mapping statistics
    """
    
    print("Loading EGIO ortholog mapping data...")
    # Load EGIO data
    egio_df = pd.read_csv(egio_file, sep='\t')
    print(f"Loaded {len(egio_df)} ortholog mappings from EGIO")
    
    # Filter by identity threshold
    egio_filtered = egio_df[egio_df['Iden'] >= identity_threshold].copy()
    print(f"Filtered to {len(egio_filtered)} mappings with identity >= {identity_threshold}")
    
    # Load exon coordinate mapping files
    print("Loading exon coordinate mapping files...")
    human_exon_map = load_exon_mapping_file(human_exon_mapping_file)
    mouse_exon_map = load_exon_mapping_file(mouse_exon_mapping_file)
    print(f"Loaded {len(human_exon_map)} human exon mappings")
    print(f"Loaded {len(mouse_exon_map)} mouse exon mappings")
    
    # Create human to mouse coordinate mapping
    print("Creating coordinate mappings...")
    human_to_mouse_map = {}
    mouse_to_human_map = {}
    
    for _, row in egio_filtered.iterrows():
        human_coord = row['hsaPos']
        mouse_coord = row['musPos'] 
        identity = row['Iden']
        
        # Store mapping with identity score
        human_to_mouse_map[human_coord] = {
            'mouse_coord': mouse_coord,
            'identity': identity
        }
        mouse_to_human_map[mouse_coord] = {
            'human_coord': human_coord, 
            'identity': identity
        }
    
    print(f"Created {len(human_to_mouse_map)} human->mouse coordinate mappings")
    
    # Load human CRPEs
    print(f"Loading human CRPEs from {human_crpe_file}...")
    human_crpes = []
    
    with open(human_crpe_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and ',' in line:
                exon1, exon2 = line.split(',')
                human_crpes.append((exon1.strip(), exon2.strip()))
    
    print(f"Loaded {len(human_crpes)} human CRPE pairs")
    
    # Map CRPEs to mouse orthologs using exon coordinate mappings
    print("Mapping human CRPEs to mouse orthologs...")
    mapped_eeis = []
    mapping_stats = {
        'total_human_crpes': len(human_crpes),
        'both_exons_mapped': 0,
        'one_exon_mapped': 0,
        'no_exons_mapped': 0,
        'successful_mappings': 0
    }
    
    for human_exon1, human_exon2 in human_crpes:
        # Get coordinates for human exons using the mapping file
        human_coord1 = human_exon_map.get(human_exon1)
        human_coord2 = human_exon_map.get(human_exon2)
        
        exons_found = 0
        if human_coord1:
            exons_found += 1
        if human_coord2:
            exons_found += 1
            
        # Try to map to mouse coordinates
        mouse_coord1 = None
        mouse_coord2 = None
        identity1 = None
        identity2 = None
        
        if human_coord1 and human_coord1 in human_to_mouse_map:
            mouse_mapping1 = human_to_mouse_map[human_coord1]
            mouse_coord1 = mouse_mapping1['mouse_coord']
            identity1 = mouse_mapping1['identity']
            
        if human_coord2 and human_coord2 in human_to_mouse_map:
            mouse_mapping2 = human_to_mouse_map[human_coord2]
            mouse_coord2 = mouse_mapping2['mouse_coord']
            identity2 = mouse_mapping2['identity']
        
        # Count mapping success
        if mouse_coord1 and mouse_coord2:
            mapping_stats['both_exons_mapped'] += 1
            mapping_stats['successful_mappings'] += 1
            
            mapped_eeis.append({
                'mouse_exon1': mouse_coord1,
                'mouse_exon2': mouse_coord2,
                'human_exon1': human_exon1,
                'human_exon2': human_exon2,
                'human_coord1': human_coord1,
                'human_coord2': human_coord2,
                'identity1': identity1,
                'identity2': identity2,
                'avg_identity': (identity1 + identity2) / 2
            })
        elif mouse_coord1 or mouse_coord2:
            mapping_stats['one_exon_mapped'] += 1
        else:
            mapping_stats['no_exons_mapped'] += 1
    
    # Convert to DataFrame
    mapped_df = pd.DataFrame(mapped_eeis)
    
    # Calculate additional statistics
    if len(mapped_df) > 0:
        mapping_stats['avg_identity'] = mapped_df['avg_identity'].mean()
        mapping_stats['min_identity'] = mapped_df[['identity1', 'identity2']].min().min()
        mapping_stats['max_identity'] = mapped_df[['identity1', 'identity2']].max().max()
    else:
        mapping_stats['avg_identity'] = 0
        mapping_stats['min_identity'] = 0
        mapping_stats['max_identity'] = 0
    
    mapping_stats['mapping_rate'] = mapping_stats['successful_mappings'] / mapping_stats['total_human_crpes']
    
    # Print statistics
    print("\n=== MAPPING STATISTICS ===")
    print(f"Total human CRPEs: {mapping_stats['total_human_crpes']}")
    print(f"Successfully mapped: {mapping_stats['successful_mappings']} ({mapping_stats['mapping_rate']:.2%})")
    print(f"Both exons mapped: {mapping_stats['both_exons_mapped']}")
    print(f"One exon mapped: {mapping_stats['one_exon_mapped']}")
    print(f"No exons mapped: {mapping_stats['no_exons_mapped']}")
    if mapping_stats['successful_mappings'] > 0:
        print(f"Average identity: {mapping_stats['avg_identity']:.3f}")
        print(f"Identity range: {mapping_stats['min_identity']:.3f} - {mapping_stats['max_identity']:.3f}")
    
    # Save results if output file specified
    if output_file and len(mapped_df) > 0:
        mapped_df.to_csv(output_file, sep='\t', index=False)
        print(f"\nSaved {len(mapped_df)} mapped mouse EEIs to {output_file}")
        
        # Also save a simplified version with just mouse coordinates
        simple_output = output_file.replace('.txt', '_simple.txt').replace('.tsv', '_simple.tsv')
        simple_df = mapped_df[['mouse_exon1', 'mouse_exon2']].copy()
        simple_df.to_csv(simple_output, sep='\t', index=False, header=False)
        print(f"Saved simplified mouse EEI list to {simple_output}")
    
    return mapped_df, mapping_stats

def load_exon_mapping_file(mapping_file):
    """
    Load exon mapping file (exon_id -> coordinate mapping).
    
    Expected format:
    exon_id    coord
    ENSE00001724518    chr1:1000:2000:+1
    """
    try:
        df = pd.read_csv(mapping_file, sep='\t')
        mapping = {}
        for _, row in df.iterrows():
            mapping[row['exon_id']] = row['coord']
        return mapping
    except Exception as e:
        print(f"Warning: Could not load exon mapping file {mapping_file}: {e}")
        return {}

if __name__ == "__main__":
    egio_file = "EGIO/ExonGroup_testpro_hsa_mus.txt"
    human_crpe_file = "CRPE_human/1_NETLOW_CRPES.txt"
    human_exon_mapping_file = "exon_mapping_files/exon_coord_mapping_human.tsv"
    mouse_exon_mapping_file = "exon_mapping_files/exon_coord_mapping_mouse.tsv"
    output_file = "outputs/mapped_mouse_crpes_low.tsv"
    
    # Run the mapping
    mapped_eeis, stats = map_human_crpes_to_mouse_orthologs(
        egio_file=egio_file,
        human_crpe_file=human_crpe_file,
        human_exon_mapping_file=human_exon_mapping_file,
        mouse_exon_mapping_file=mouse_exon_mapping_file,
        identity_threshold=0.0,
        output_file=output_file
    )
    
    print(f"\nMapping completed. Found {len(mapped_eeis)} orthologous mouse EEI pairs.")