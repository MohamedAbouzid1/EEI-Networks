#!/usr/bin/env python3
"""
Script to create phylogenetic tree visualization with EEI counts mapped to branch tips.

This script creates a phylogenetic tree showing the evolutionary relationships between
species and maps the EEI network sizes (Contact, PISA, EPPIC) to the terminal branches
for easy comparison across species and detection methods.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os
from pathlib import Path

# Try to import ete3, if not available, use alternative approach
try:
    from ete3 import Tree, TreeStyle, TextFace, add_face_to_node, NodeStyle, faces, AttrFace
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False
    print("Warning: ete3 not available. Using matplotlib-only approach.")

def create_species_data():
    """Create the species EEI data from Table 3.1."""
    
    species_data = {
        'Mus_musculus': {
            'common_name': 'House Mouse',
            'contact': 2322,
            'pisa': 1387,
            'eppic': 1952
        },
        'Saccharomyces_cerevisiae': {
            'common_name': 'Baker\'s Yeast',
            'contact': 1465,
            'pisa': None,  # No data
            'eppic': None  # No data
        },
        'Rattus_norvegicus': {
            'common_name': 'Brown Rat',
            'contact': 160,
            'pisa': 127,
            'eppic': 106
        },
        'Bos_taurus': {
            'common_name': 'Cattle',
            'contact': 425,
            'pisa': None,  # No data
            'eppic': 223
        },
        'Drosophila_melanogaster': {
            'common_name': 'Fruit Fly',
            'contact': 354,
            'pisa': 125,
            'eppic': 94
        },
        'Gallus_gallus': {
            'common_name': 'Chicken',
            'contact': 51,
            'pisa': 44,
            'eppic': 107
        },
        'Oryctolagus_cuniculus': {
            'common_name': 'European Rabbit',
            'contact': 175,
            'pisa': None,  # No data
            'eppic': 53
        }
    }
    
    return species_data

def create_newick_tree():
    """Create a simplified phylogenetic tree in Newick format."""
    
    # Simplified tree topology based on evolutionary relationships
    # This is a simplified representation for visualization purposes
    newick_tree = "((((Mus_musculus,Rattus_norvegicus),Oryctolagus_cuniculus),(Bos_taurus,Gallus_gallus)),Drosophila_melanogaster,Saccharomyces_cerevisiae);"
    
    return newick_tree

def create_tree_with_ete3(species_data, output_dir):
    """Create phylogenetic tree using ete3 library."""
    
    print("Creating phylogenetic tree with ete3...")
    
    # Create tree from Newick format
    newick_tree = create_newick_tree()
    tree = Tree(newick_tree)
    
    # Set up tree style
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.rotation = 0
    ts.scale = 20
    ts.mode = "r"  # rectangular mode
    ts.arc_start = -180
    ts.arc_span = 180
    
    # Define colors for different methods
    colors = {
        'contact': '#FF6B6B',    # Red
        'pisa': '#4ECDC4',       # Teal
        'eppic': '#45B7D1'       # Blue
    }
    
    # Add EEI count faces to leaf nodes
    for leaf in tree.iter_leaves():
        species = leaf.name
        if species in species_data:
            data = species_data[species]
            
            # Create a container for the faces
            container = faces.SeqGroupFace(
                width=300,
                height=60,
                text_size=8,
                col_width=100
            )
            
            # Add method faces
            methods = ['contact', 'pisa', 'eppic']
            for method in methods:
                count = data.get(method)
                if count is not None:
                    face = TextFace(f"{count}", fsize=8, fgcolor='black', bgcolor=colors[method])
                    face.margin_left = 2
                    face.margin_right = 2
                    container.add_face(face, column=methods.index(method))
                else:
                    face = TextFace("-", fsize=8, fgcolor='gray')
                    face.margin_left = 2
                    face.margin_right = 2
                    container.add_face(face, column=methods.index(method))
            
            # Add the container to the leaf
            add_face_to_node(container, leaf, column=1, position="branch-right")
            
            # Add common name
            common_name = data['common_name']
            name_face = TextFace(f"{common_name}", fsize=10, fgcolor='black')
            add_face_to_node(name_face, leaf, column=2, position="branch-right")
    
    # Render the tree
    output_path = os.path.join(output_dir, 'phylogenetic_tree_eei_ete3.png')
    tree.render(output_path, tree_style=ts, dpi=300, w=800, h=600)
    print(f"ETE3 tree saved to: {output_path}")
    
    return tree, ts

def create_tree_with_matplotlib(species_data, output_dir):
    """Create phylogenetic tree using matplotlib (improved version)."""
    
    print("Creating improved phylogenetic tree with matplotlib...")
    
    # Create figure with better proportions
    fig, ax = plt.subplots(1, 1, figsize=(18, 12))
    
    # Define improved tree structure with better spacing
    # More realistic phylogenetic relationships and branch lengths
    tree_structure = {
        'root': {'x': 0, 'y': 5},
        'eukaryote_split': {'x': 1.5, 'y': 5},
        'animal_branch': {'x': 2.5, 'y': 4.5},
        'vertebrate_branch': {'x': 3.5, 'y': 4.2},
        'mammal_branch': {'x': 4.5, 'y': 3.8},
        'eutherian_branch': {'x': 5.5, 'y': 3.5},
        'rodent_branch': {'x': 6.5, 'y': 4.0},
        'Mus_musculus': {'x': 7.5, 'y': 4.5},
        'Rattus_norvegicus': {'x': 7.5, 'y': 3.5},
        'lagomorpha_branch': {'x': 6.5, 'y': 2.8},
        'Oryctolagus_cuniculus': {'x': 7.5, 'y': 2.8},
        'artiodactyl_branch': {'x': 6.5, 'y': 2.0},
        'Bos_taurus': {'x': 7.5, 'y': 2.0},
        'bird_branch': {'x': 5.5, 'y': 1.5},
        'Gallus_gallus': {'x': 7.5, 'y': 1.5},
        'invertebrate_branch': {'x': 3.5, 'y': 1.0},
        'Drosophila_melanogaster': {'x': 5.0, 'y': 1.0},
        'fungi_branch': {'x': 2.5, 'y': 0.5},
        'Saccharomyces_cerevisiae': {'x': 4.0, 'y': 0.5}
    }
    
    # Define improved colors with better contrast
    colors = {
        'contact': '#E74C3C',    # Vibrant Red
        'pisa': '#3498DB',       # Vibrant Blue  
        'eppic': '#2ECC71'       # Vibrant Green
    }
    
    # Draw improved tree branches with varying thickness
    branch_width = 3
    
    # Main branches
    ax.plot([0, 1.5], [5, 5], 'k-', linewidth=branch_width, alpha=0.8)  # Root to eukaryote split
    ax.plot([1.5, 2.5], [5, 4.5], 'k-', linewidth=branch_width-1, alpha=0.8)  # To animal branch
    ax.plot([2.5, 3.5], [4.5, 4.2], 'k-', linewidth=branch_width-1, alpha=0.8)  # To vertebrate
    ax.plot([3.5, 4.5], [4.2, 3.8], 'k-', linewidth=branch_width-1, alpha=0.8)  # To mammal
    ax.plot([4.5, 5.5], [3.8, 3.5], 'k-', linewidth=branch_width-1, alpha=0.8)  # To eutherian
    ax.plot([1.5, 2.5], [5, 0.5], 'k-', linewidth=branch_width-1, alpha=0.8)  # To fungi
    ax.plot([2.5, 3.5], [0.5, 1.0], 'k-', linewidth=branch_width-1, alpha=0.8)  # Fungi to invertebrate
    ax.plot([3.5, 5.0], [1.0, 1.0], 'k-', linewidth=branch_width-1, alpha=0.8)  # To fly
    
    # Rodent branches
    ax.plot([5.5, 6.5], [3.5, 4.0], 'k-', linewidth=branch_width-2, alpha=0.8)  # To rodent branch
    ax.plot([6.5, 7.5], [4.0, 4.5], 'k-', linewidth=branch_width-2, alpha=0.8)  # To mouse
    ax.plot([6.5, 7.5], [4.0, 3.5], 'k-', linewidth=branch_width-2, alpha=0.8)  # To rat
    
    # Other mammal branches
    ax.plot([5.5, 6.5], [3.5, 2.8], 'k-', linewidth=branch_width-2, alpha=0.8)  # To lagomorpha
    ax.plot([6.5, 7.5], [2.8, 2.8], 'k-', linewidth=branch_width-2, alpha=0.8)  # To rabbit
    ax.plot([5.5, 6.5], [3.5, 2.0], 'k-', linewidth=branch_width-2, alpha=0.8)  # To artiodactyl
    ax.plot([6.5, 7.5], [2.0, 2.0], 'k-', linewidth=branch_width-2, alpha=0.8)  # To cattle
    
    # Bird branch
    ax.plot([4.5, 5.5], [3.8, 1.5], 'k-', linewidth=branch_width-2, alpha=0.8)  # To bird branch
    ax.plot([5.5, 7.5], [1.5, 1.5], 'k-', linewidth=branch_width-2, alpha=0.8)  # To chicken
    
    # Add terminal branches and species labels
    for species in ['Mus_musculus', 'Rattus_norvegicus', 'Oryctolagus_cuniculus', 
                   'Bos_taurus', 'Gallus_gallus', 'Drosophila_melanogaster', 'Saccharomyces_cerevisiae']:
        if species in tree_structure:
            pos = tree_structure[species]
            
            # Draw terminal branch
            ax.plot([pos['x'], pos['x'] + 0.3], [pos['y'], pos['y']], 'k-', linewidth=2, alpha=0.8)
            
            # Add species name with better formatting
            common_name = species_data[species]['common_name']
            ax.text(pos['x'] + 0.4, pos['y'], common_name, 
                   fontsize=11, va='center', ha='left', fontweight='bold')
            
            # Add scientific name in italics
            sci_name = species.replace('_', ' ')
            ax.text(pos['x'] + 0.4, pos['y'] - 0.15, f"({sci_name})", 
                   fontsize=9, va='center', ha='left', style='italic', color='gray')
            
            # Add EEI counts as improved colored bars
            if species in species_data:
                data = species_data[species]
                methods = ['contact', 'pisa', 'eppic']
                
                # Calculate bar positions with better spacing
                bar_height = 0.08
                bar_start = pos['x'] + 2.0
                bar_spacing = 0.25
                
                # Find maximum count for normalization
                max_count = 0
                for method in methods:
                    method_max = max([d.get(method, 0) for d in species_data.values() if d.get(method) is not None])
                    max_count = max(max_count, method_max)
                
                for i, method in enumerate(methods):
                    count = data.get(method)
                    if count is not None:
                        # Normalize count for visualization
                        normalized_count = (count / max_count) * 1.5 if max_count > 0 else 0
                        
                        # Draw bar with better styling
                        bar_x = bar_start + i * bar_spacing
                        
                        # Create rounded rectangle for bars
                        from matplotlib.patches import FancyBboxPatch
                        bar_patch = FancyBboxPatch(
                            (bar_x, pos['y'] - bar_height/2), 
                            normalized_count, bar_height,
                            boxstyle="round,pad=0.01",
                            facecolor=colors[method], 
                            edgecolor='white',
                            linewidth=1,
                            alpha=0.8
                        )
                        ax.add_patch(bar_patch)
                        
                        # Add count label with better positioning
                        if normalized_count > 0.1:  # Only add label if bar is visible
                            ax.text(bar_x + normalized_count/2, pos['y'], str(count), 
                                   fontsize=9, ha='center', va='center', 
                                   fontweight='bold', color='white')
                    else:
                        # Draw dash for missing data
                        bar_x = bar_start + i * bar_spacing
                        ax.text(bar_x + 0.075, pos['y'], '-', 
                               fontsize=14, ha='center', va='center', color='gray', fontweight='bold')
    
    # Add improved method labels with better positioning
    method_labels = ['CONTACT', 'PISA', 'EPPIC']
    label_y = 6.0
    
    for i, (method, label) in enumerate(zip(['contact', 'pisa', 'eppic'], method_labels)):
        bar_x = 9.5 + i * 0.25
        
        # Add colored box for method label
        from matplotlib.patches import FancyBboxPatch
        label_box = FancyBboxPatch(
            (bar_x - 0.08, label_y - 0.05), 0.16, 0.1,
            boxstyle="round,pad=0.01",
            facecolor=colors[method], 
            edgecolor='white',
            linewidth=1,
            alpha=0.8
        )
        ax.add_patch(label_box)
        
        ax.text(bar_x, label_y, label, fontsize=10, ha='center', va='center', 
               fontweight='bold', color='white')
    
    # Add evolutionary time scale (simplified)
    ax.text(0.75, 5.3, 'Common Ancestor', fontsize=10, ha='center', va='center', 
           style='italic', color='darkblue', fontweight='bold')
    
    # Add major evolutionary groups
    group_labels = [
        (1.5, 4.8, 'Eukaryotes', 'darkgreen'),
        (3.0, 3.9, 'Animals', 'darkblue'),
        (4.0, 3.3, 'Vertebrates', 'purple'),
        (5.0, 2.9, 'Mammals', 'brown'),
        (2.5, 1.8, 'Invertebrates', 'orange'),
        (2.0, 0.8, 'Fungi', 'darkred')
    ]
    
    for x, y, label, color in group_labels:
        ax.text(x, y, label, fontsize=9, ha='center', va='center', 
               fontweight='bold', color=color, 
               bbox=dict(boxstyle="round,pad=0.2", facecolor='white', edgecolor=color, alpha=0.7))
    
    # Customize plot with better limits
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(-0.5, 6.5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Add improved title
    plt.title('Phylogenetic Tree with Exon-Exon Interaction (EEI) Network Sizes\nAcross Species and Detection Methods', 
              fontsize=18, fontweight='bold', pad=30)
    
    # Add improved legend
    legend_elements = [
        mpatches.Patch(color=colors['contact'], label='CONTACT', alpha=0.8),
        mpatches.Patch(color=colors['pisa'], label='PISA', alpha=0.8),
        mpatches.Patch(color=colors['eppic'], label='EPPIC', alpha=0.8)
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.02, 0.98),
             frameon=True, fancybox=True, shadow=True, fontsize=12)
    
    # Add informative note
    note_text = ('Note: Bar lengths are proportional to EEI counts. "-" indicates missing data.\n'
                'Species are arranged according to evolutionary relationships.')
    ax.text(0.5, -0.3, note_text, fontsize=10, ha='left', style='italic', 
           color='gray', bbox=dict(boxstyle="round,pad=0.5", facecolor='lightgray', alpha=0.5))
    
    # Add grid lines for better readability (optional)
    ax.grid(True, alpha=0.1, linestyle='--')
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(output_dir, 'phylogenetic_tree_eei_matplotlib.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Matplotlib tree saved to: {output_path}")
    
    plt.show()
    
    return fig

def create_comparative_plot(species_data, output_dir):
    """Create an improved comparative bar plot of EEI counts across species."""
    
    print("Creating improved comparative bar plot...")
    
    # Prepare data for plotting
    species_names = []
    contact_counts = []
    pisa_counts = []
    eppic_counts = []
    
    for species, data in species_data.items():
        species_names.append(data['common_name'])
        contact_counts.append(data['contact'] if data['contact'] is not None else 0)
        pisa_counts.append(data['pisa'] if data['pisa'] is not None else 0)
        eppic_counts.append(data['eppic'] if data['eppic'] is not None else 0)
    
    # Create figure with better proportions
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    
    # Define improved colors
    colors = ['#E74C3C', '#3498DB', '#2ECC71']  # Match phylogenetic tree colors
    methods = ['CONTACT', 'PISA', 'EPPIC']
    
    # Set up bar positions
    x = np.arange(len(species_names))
    width = 0.28
    
    # Create bars with better styling
    bars1 = ax.bar(x - width, contact_counts, width, label='CONTACT', 
                  color=colors[0], alpha=0.8, edgecolor='white', linewidth=1.5)
    bars2 = ax.bar(x, pisa_counts, width, label='PISA', 
                  color=colors[1], alpha=0.8, edgecolor='white', linewidth=1.5)
    bars3 = ax.bar(x + width, eppic_counts, width, label='EPPIC', 
                  color=colors[2], alpha=0.8, edgecolor='white', linewidth=1.5)
    
    # Add value labels on bars with better positioning
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:  # Only add labels for non-zero values
                ax.text(bar.get_x() + bar.get_width()/2., height + 25,
                       f'{int(height):,}', ha='center', va='bottom', 
                       fontweight='bold', fontsize=9)
    
    # Customize plot with better styling
    ax.set_xlabel('Species', fontsize=14, fontweight='bold')
    ax.set_ylabel('Number of EEI Interactions', fontsize=14, fontweight='bold')
    ax.set_title('Comparative Analysis: EEI Network Sizes Across Species and Detection Methods', 
                fontsize=16, fontweight='bold', pad=25)
    ax.set_xticks(x)
    ax.set_xticklabels(species_names, rotation=45, ha='right', fontsize=11)
    
    # Improve legend
    ax.legend(loc='upper right', fontsize=12, frameon=True, 
             fancybox=True, shadow=True)
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    ax.set_axisbelow(True)
    
    # Add statistical summary
    total_contact = sum(contact_counts)
    total_pisa = sum(pisa_counts)
    total_eppic = sum(eppic_counts)
    
    summary_text = f'Total EEIs: CONTACT={total_contact:,} | PISA={total_pisa:,} | EPPIC={total_eppic:,}'
    ax.text(0.5, 0.95, summary_text, transform=ax.transAxes, 
           fontsize=11, ha='center', va='top', 
           bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.7))
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(output_dir, 'eei_counts_comparative_improved.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Improved comparative plot saved to: {output_path}")
    
    plt.show()
    
    return fig

def save_data_summary(species_data, output_dir):
    """Save the EEI data as a CSV file for reference."""
    
    # Prepare data for DataFrame
    data_rows = []
    for species, data in species_data.items():
        data_rows.append({
            'Species': species,
            'Common_Name': data['common_name'],
            'CONTACT': data['contact'],
            'PISA': data['pisa'],
            'EPPIC': data['eppic']
        })
    
    df = pd.DataFrame(data_rows)
    
    # Save to CSV
    output_path = os.path.join(output_dir, 'eei_counts_by_species.csv')
    df.to_csv(output_path, index=False)
    print(f"Data summary saved to: {output_path}")
    
    return df

def main():
    """Main function to create phylogenetic tree visualization."""
    
    # Create output directory
    output_dir = "/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/final_statistics/phylogenetic_analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    print("EEI Phylogenetic Tree Analysis")
    print("="*50)
    
    # Create species data
    species_data = create_species_data()
    
    # Save data summary
    df = save_data_summary(species_data, output_dir)
    
    # Create phylogenetic tree visualization
    if ETE3_AVAILABLE:
        try:
            create_tree_with_ete3(species_data, output_dir)
        except Exception as e:
            print(f"ETE3 tree creation failed: {e}")
            print("Falling back to matplotlib approach...")
            create_tree_with_matplotlib(species_data, output_dir)
    else:
        print("ETE3 not available, using matplotlib approach...")
        create_tree_with_matplotlib(species_data, output_dir)
    
    # Create comparative plot
    create_comparative_plot(species_data, output_dir)
    
    print(f"\nAnalysis complete! All visualizations saved to: {output_dir}")
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print("Species EEI Network Sizes:")
    for species, data in species_data.items():
        print(f"\n{data['common_name']} ({species}):")
        print(f"  CONTACT: {data['contact'] if data['contact'] is not None else 'N/A'}")
        print(f"  PISA: {data['pisa'] if data['pisa'] is not None else 'N/A'}")
        print(f"  EPPIC: {data['eppic'] if data['eppic'] is not None else 'N/A'}")

if __name__ == "__main__":
    main()
