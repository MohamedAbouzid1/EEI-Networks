#!/usr/bin/env python3
"""
Script to generate Venn diagrams showing overlap between three EEI detection methods:
1. CONTACT
2. PISA  
3. EPPIC

This script analyzes the exon-exon interactions detected by each method and creates
Venn diagrams to visualize their overlaps.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import os
from pathlib import Path

def load_contact_data(file_path):
    """Load and process CONTACT results."""
    print("Loading CONTACT data...")
    df = pd.read_csv(file_path, sep='\t')
    
    # Create unique interaction identifiers using exon1 and exon2
    df['interaction'] = df['exon1'] + '_' + df['exon2']
    
    # Remove any rows with NaN values in interaction column
    df = df.dropna(subset=['interaction'])
    
    # Get unique interactions (filter out any remaining NaN values)
    interactions = set(df['interaction'].dropna().unique())
    print(f"CONTACT: Found {len(interactions)} unique exon-exon interactions")
    
    return interactions, df

def load_pisa_data(file_path):
    """Load and process PISA results."""
    print("Loading PISA data...")
    df = pd.read_csv(file_path, sep='\t')
    
    # Create unique interaction identifiers using exon1 and exon2
    df['interaction'] = df['exon1'] + '_' + df['exon2']
    
    # Remove any rows with NaN values in interaction column
    df = df.dropna(subset=['interaction'])
    
    # Get unique interactions (filter out any remaining NaN values)
    interactions = set(df['interaction'].dropna().unique())
    print(f"PISA: Found {len(interactions)} unique exon-exon interactions")
    
    return interactions, df

def load_eppic_data(file_path):
    """Load and process EPPIC results."""
    print("Loading EPPIC data...")
    df = pd.read_csv(file_path, sep='\t')
    
    # Create unique interaction identifiers using exon1 and exon2
    df['interaction'] = df['exon1'] + '_' + df['exon2']
    
    # Remove any rows with NaN values in interaction column
    df = df.dropna(subset=['interaction'])
    
    # Get unique interactions (filter out any remaining NaN values)
    interactions = set(df['interaction'].dropna().unique())
    print(f"EPPIC: Found {len(interactions)} unique exon-exon interactions")
    
    return interactions, df

def create_venn_diagram(contact_set, pisa_set, eppic_set, output_dir):
    """Create Venn diagram showing overlaps between the three methods."""
    
    # Calculate overlaps
    contact_only = contact_set - pisa_set - eppic_set
    pisa_only = pisa_set - contact_set - eppic_set
    eppic_only = eppic_set - contact_set - pisa_set
    
    contact_pisa = (contact_set & pisa_set) - eppic_set
    contact_eppic = (contact_set & eppic_set) - pisa_set
    pisa_eppic = (pisa_set & eppic_set) - contact_set
    
    all_three = contact_set & pisa_set & eppic_set
    
    print("\nOverlap Statistics:")
    print(f"CONTACT only: {len(contact_only)}")
    print(f"PISA only: {len(pisa_only)}")
    print(f"EPPIC only: {len(eppic_only)}")
    print(f"CONTACT + PISA: {len(contact_pisa)}")
    print(f"CONTACT + EPPIC: {len(contact_eppic)}")
    print(f"PISA + EPPIC: {len(pisa_eppic)}")
    print(f"All three methods: {len(all_three)}")
    
    # Create Venn diagram
    plt.figure(figsize=(12, 10))
    
    # Create the Venn diagram
    v = venn3([contact_set, pisa_set, eppic_set], 
              set_labels=('CONTACT', 'PISA', 'EPPIC'),
              alpha=0.7)
    
    # Customize colors
    if v.get_patch_by_id('100'):
        v.get_patch_by_id('100').set_color('#ff9999')  # Red for CONTACT only
    if v.get_patch_by_id('010'):
        v.get_patch_by_id('010').set_color('#99ff99')  # Green for PISA only
    if v.get_patch_by_id('001'):
        v.get_patch_by_id('001').set_color('#9999ff')  # Blue for EPPIC only
    if v.get_patch_by_id('110'):
        v.get_patch_by_id('110').set_color('#ffff99')  # Yellow for CONTACT+PISA
    if v.get_patch_by_id('101'):
        v.get_patch_by_id('101').set_color('#ff99ff')  # Magenta for CONTACT+EPPIC
    if v.get_patch_by_id('011'):
        v.get_patch_by_id('011').set_color('#99ffff')  # Cyan for PISA+EPPIC
    if v.get_patch_by_id('111'):
        v.get_patch_by_id('111').set_color('#ffcc99')  # Orange for all three
    
    # Add circles for better visibility
    c = venn3_circles([contact_set, pisa_set, eppic_set], linestyle='solid', linewidth=2)
    
    # Set title and labels
    plt.title('Overlap of Exon-Exon Interactions Detected by Three Methods\n(CONTACT, PISA, EPPIC)', 
              fontsize=16, fontweight='bold', pad=20)
    
    # Add statistics text
    stats_text = f"""Statistics:
Total unique interactions: {len(contact_set | pisa_set | eppic_set)}
CONTACT: {len(contact_set)} | PISA: {len(pisa_set)} | EPPIC: {len(eppic_set)}
Overlap (all three): {len(all_three)} ({len(all_three)/len(contact_set | pisa_set | eppic_set)*100:.1f}%)"""
    
    plt.figtext(0.02, 0.02, stats_text, fontsize=10, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(output_dir, 'venn_diagram_eei_overlap.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nVenn diagram saved to: {output_path}")
    
    plt.show()
    
    return {
        'contact_only': contact_only,
        'pisa_only': pisa_only,
        'eppic_only': eppic_only,
        'contact_pisa': contact_pisa,
        'contact_eppic': contact_eppic,
        'pisa_eppic': pisa_eppic,
        'all_three': all_three
    }

def save_detailed_results(overlaps, output_dir):
    """Save detailed results to text files."""
    
    results_dir = os.path.join(output_dir, 'overlap_analysis')
    os.makedirs(results_dir, exist_ok=True)
    
    # Save each overlap category
    for category, interactions in overlaps.items():
        if interactions:
            # Filter out any NaN values and convert to strings
            valid_interactions = [str(interaction) for interaction in interactions if pd.notna(interaction)]
            
            if valid_interactions:
                filename = os.path.join(results_dir, f'{category}_interactions.txt')
                with open(filename, 'w') as f:
                    f.write(f"# {category.replace('_', ' ').title()} Interactions\n")
                    f.write(f"# Count: {len(valid_interactions)}\n\n")
                    for interaction in sorted(valid_interactions):
                        f.write(f"{interaction}\n")
                print(f"Saved {len(valid_interactions)} {category} interactions to {filename}")

def generate_summary_statistics(contact_set, pisa_set, eppic_set, overlaps, output_dir):
    """Generate and save summary statistics."""
    
    total_unique = len(contact_set | pisa_set | eppic_set)
    
    stats = {
        'Method': ['CONTACT', 'PISA', 'EPPIC', 'All Methods'],
        'Total_Interactions': [len(contact_set), len(pisa_set), len(eppic_set), total_unique],
        'Unique_Only': [len(overlaps['contact_only']), len(overlaps['pisa_only']), len(overlaps['eppic_only']), 0],
        'Shared_with_Others': [
            len(contact_set) - len(overlaps['contact_only']),
            len(pisa_set) - len(overlaps['pisa_only']),
            len(eppic_set) - len(overlaps['eppic_only']),
            0
        ]
    }
    
    stats_df = pd.DataFrame(stats)
    
    # Add percentages
    stats_df['Unique_Percentage'] = (stats_df['Unique_Only'] / stats_df['Total_Interactions'] * 100).round(1)
    stats_df['Shared_Percentage'] = (stats_df['Shared_with_Others'] / stats_df['Total_Interactions'] * 100).round(1)
    
    # Save statistics
    stats_file = os.path.join(output_dir, 'overlap_statistics.csv')
    stats_df.to_csv(stats_file, index=False)
    print(f"Statistics saved to: {stats_file}")
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(stats_df.to_string(index=False))
    
    # Calculate additional metrics
    print(f"\nAdditional Metrics:")
    print(f"Jaccard similarity (CONTACT vs PISA): {len(contact_set & pisa_set) / len(contact_set | pisa_set):.3f}")
    print(f"Jaccard similarity (CONTACT vs EPPIC): {len(contact_set & eppic_set) / len(contact_set | eppic_set):.3f}")
    print(f"Jaccard similarity (PISA vs EPPIC): {len(pisa_set & eppic_set) / len(pisa_set | eppic_set):.3f}")
    print(f"Consensus interactions (all three): {len(overlaps['all_three'])} ({len(overlaps['all_three'])/total_unique*100:.1f}%)")

def main():
    """Main function to run the analysis."""
    
    # Define paths
    base_dir = "/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Homo-Sapiens-2024/results"
    contact_file = os.path.join(base_dir, "1-CONTACT_results", "human_network_final.txt")
    pisa_file = os.path.join(base_dir, "2-PISA_results", "PISA_EEIN_0.5_human.txt")
    eppic_file = os.path.join(base_dir, "3-EPPIC_network", "EPPIC_EEIN_filtered.txt")
    
    # Create output directory
    output_dir = "/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/final_statistics/venn_analysis_results"
    os.makedirs(output_dir, exist_ok=True)
    
    print("EEI Detection Methods Overlap Analysis")
    print("="*50)
    
    # Load data
    contact_interactions, contact_df = load_contact_data(contact_file)
    pisa_interactions, pisa_df = load_pisa_data(pisa_file)
    eppic_interactions, eppic_df = load_eppic_data(eppic_file)
    
    # Create Venn diagram
    overlaps = create_venn_diagram(contact_interactions, pisa_interactions, eppic_interactions, output_dir)
    
    # Save detailed results
    save_detailed_results(overlaps, output_dir)
    
    # Generate summary statistics
    generate_summary_statistics(contact_interactions, pisa_interactions, eppic_interactions, overlaps, output_dir)
    
    print(f"\nAnalysis complete! Results saved to: {output_dir}")

if __name__ == "__main__":
    main()
