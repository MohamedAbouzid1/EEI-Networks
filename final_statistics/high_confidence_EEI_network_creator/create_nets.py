import pandas as pd
import sys
from pathlib import Path
from data_loader import load_contact_network, load_eppic_network, load_pisa_network


def find_high_confidence_eeis(contact_pairs, eppic_pairs, pisa_pairs):
    """Find EEI pairs that exist in all three networks."""
    print("\nFinding high-confidence EEIs (present in all three networks)...")
    
    # Find intersection of all three sets
    high_confidence = contact_pairs & eppic_pairs & pisa_pairs
    
    print(f"  Found {len(high_confidence)} high-confidence EEI pairs")
    
    # Also calculate pairwise intersections for statistics
    contact_eppic = contact_pairs & eppic_pairs
    contact_pisa = contact_pairs & pisa_pairs
    eppic_pisa = eppic_pairs & pisa_pairs
    
    print(f"\nPairwise overlaps:")
    print(f"  Contact ∩ EPPIC: {len(contact_eppic)} EEIs")
    print(f"  Contact ∩ PISA: {len(contact_pisa)} EEIs")
    print(f"  EPPIC ∩ PISA: {len(eppic_pisa)} EEIs")
    
    return high_confidence

def create_high_confidence_network(high_conf_pairs, contact_df, eppic_df, pisa_df, 
                                  contact_pmap, eppic_pmap, pisa_pmap):
    """Create a comprehensive high-confidence network with data from all sources."""
    
    results = []
    
    for exon_pair in high_conf_pairs:
        exon1, exon2 = exon_pair
        
        # Get protein information (should be consistent across networks)
        proteins = contact_pmap.get(exon_pair, eppic_pmap.get(exon_pair, pisa_pmap.get(exon_pair)))
        
        # Find matching rows in each network
        # Handle both directions of the interaction
        contact_rows = contact_df[
            ((contact_df['exon1'] == exon1) & (contact_df['exon2'] == exon2)) |
            ((contact_df['exon1'] == exon2) & (contact_df['exon2'] == exon1))
        ]
        
        eppic_rows = eppic_df[
            ((eppic_df['exon1'] == exon1) & (eppic_df['exon2'] == exon2)) |
            ((eppic_df['exon1'] == exon2) & (eppic_df['exon2'] == exon1))
        ]
        
        pisa_rows = pisa_df[
            ((pisa_df['exon1'] == exon1) & (pisa_df['exon2'] == exon2)) |
            ((pisa_df['exon1'] == exon2) & (pisa_df['exon2'] == exon1))
        ]
        
        # For each combination, create a comprehensive record
        # We'll take the first match from each network for simplicity
        if not contact_rows.empty and not eppic_rows.empty and not pisa_rows.empty:
            contact_row = contact_rows.iloc[0]
            eppic_row = eppic_rows.iloc[0]
            pisa_row = pisa_rows.iloc[0]
            
            result = {
                'exon1': exon1,
                'exon2': exon2,
                'protein1': proteins[0] if proteins else '',
                'protein2': proteins[1] if proteins else '',
                
                # Contact-based metrics
                'contact_jaccard_percent': contact_row.get('jaccard_percent', None),
                'contact_exon1_coverage_percent': contact_row.get('exon1_coverage_percent', None),
                'contact_exon2_coverage_percent': contact_row.get('exon2_coverage_percent', None),
                
                # EPPIC metrics
                'eppic_buried_area': eppic_row.get('BuriedAreaAbs', None),
                'eppic_cs_score': eppic_row.get('CS_SCORE', None),
                'eppic_cr_score': eppic_row.get('CR_SCORE', None),
                'eppic_pdbid': eppic_row.get('PDBID', None),
                
                # PISA metrics
                'pisa_free_energy': pisa_row.get('FreeEnergy', None),
                'pisa_buried_area': pisa_row.get('BuriedAreaAbs', None),
                'pisa_hydrogen_bonds': pisa_row.get('Hydrogen', None),
                'pisa_salt_bridges': pisa_row.get('Saltbridge', None),
                'pisa_pdbid': pisa_row.get('PDBID', None)
            }
            
            results.append(result)
    
    return pd.DataFrame(results)