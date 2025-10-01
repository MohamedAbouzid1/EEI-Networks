#!/usr/bin/env python3
"""
Script to create network visualization figures for representative EEI networks.

This script identifies representative protein complexes/networks from each detection method
(CONTACT, PISA, EPPIC) and creates network visualizations showing exon-exon interactions.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from collections import defaultdict, Counter
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style for better plots
plt.style.use('default')
sns.set_palette("husl")

def load_contact_data(file_path):
    """Load and process CONTACT results."""
    print("Loading CONTACT data...")
    df = pd.read_csv(file_path, sep='\t')
    df = df.dropna(subset=['exon1', 'exon2'])
    print(f"CONTACT: Loaded {len(df)} interactions")
    return df

def load_pisa_data(file_path):
    """Load and process PISA results."""
    print("Loading PISA data...")
    df = pd.read_csv(file_path, sep='\t')
    df = df.dropna(subset=['exon1', 'exon2'])
    print(f"PISA: Loaded {len(df)} interactions")
    return df

def load_eppic_data(file_path):
    """Load and process EPPIC results."""
    print("Loading EPPIC data...")
    df = pd.read_csv(file_path, sep='\t')
    df = df.dropna(subset=['exon1', 'exon2'])
    print(f"EPPIC: Loaded {len(df)} interactions")
    return df

def identify_protein_complexes(df, method_name):
    """Identify protein complexes based on PDB structures."""
    print(f"\nIdentifying protein complexes for {method_name}...")
    
    # Handle different column names for different methods
    if method_name == 'CONTACT':
        protein2_col = 'protein1.1'
        # Extract PDB ID from protein names for CONTACT
        df['PDBID'] = df['protein1'].str.extract(r'_(\w+)_')
    else:
        protein2_col = 'protein2'
        if 'PDBID' not in df.columns:
            # Extract PDB ID from protein names if not present
            df['PDBID'] = df['protein1'].str.extract(r'_(\w+)_')
    
    # Group by PDB ID to identify complexes
    pdb_groups = df.groupby('PDBID')
    
    complex_stats = []
    for pdb_id, group in pdb_groups:
        if len(group) < 3:  # Skip small complexes
            continue
            
        # Get unique proteins and exons
        proteins = set(group['protein1'].unique()) | set(group[protein2_col].unique())
        exons = set(group['exon1'].unique()) | set(group['exon2'].unique())
        
        complex_stats.append({
            'PDBID': pdb_id,
            'num_interactions': len(group),
            'num_proteins': len(proteins),
            'num_exons': len(exons),
            'proteins': proteins,
            'exons': exons,
            'data': group,
            'protein2_col': protein2_col  # Store the column name for later use
        })
    
    # Sort by number of interactions and select top complexes
    complex_stats.sort(key=lambda x: x['num_interactions'], reverse=True)
    
    print(f"{method_name}: Found {len(complex_stats)} protein complexes")
    for i, complex_info in enumerate(complex_stats[:5]):
        print(f"  {i+1}. PDB {complex_info['PDBID']}: {complex_info['num_interactions']} interactions, "
              f"{complex_info['num_proteins']} proteins, {complex_info['num_exons']} exons")
    
    return complex_stats

def create_network_graph(complex_data, method_name, pdb_id, protein2_col='protein2'):
    """Create a NetworkX graph from complex data."""
    G = nx.Graph()
    
    # Add edges (exon-exon interactions)
    for _, row in complex_data.iterrows():
        exon1 = row['exon1']
        exon2 = row['exon2']
        
        # Create edge with interaction data
        edge_data = {
            'weight': 1,
            'method': method_name,
            'pdb_id': pdb_id,
            'protein1': row['protein1'],
            'protein2': row[protein2_col]
        }
        
        # Add method-specific attributes
        if method_name == 'CONTACT':
            edge_data.update({
                'jaccard_percent': row.get('jaccard_percent', 0),
                'exon1_coverage_percent': row.get('exon1_coverage_percent', 0),
                'exon2_coverage_percent': row.get('exon2_coverage_percent', 0)
            })
        elif method_name == 'PISA':
            edge_data.update({
                'free_energy': row.get('FreeEnergy', 0),
                'buried_area': row.get('BuriedArea', 0),
                'hydrogen_bonds': row.get('Hydrogen', 0)
            })
        elif method_name == 'EPPIC':
            edge_data.update({
                'cs_score': row.get('CS_SCORE', 0),
                'cr_score': row.get('CR_SCORE', 0),
                'buried_area_abs': row.get('BuriedAreaAbs', 0)
            })
        
        G.add_edge(exon1, exon2, **edge_data)
    
    return G

def visualize_network(G, complex_info, method_name, output_dir, layout_type='spring'):
    """Create network visualization."""
    
    # Set up the plot
    fig, ax = plt.subplots(1, 1, figsize=(15, 12))
    
    # Choose layout
    if layout_type == 'spring':
        pos = nx.spring_layout(G, k=3, iterations=50)
    elif layout_type == 'circular':
        pos = nx.circular_layout(G)
    elif layout_type == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.spring_layout(G)
    
    # Calculate node sizes based on degree
    degrees = dict(G.degree())
    node_sizes = [max(50, degrees[node] * 100) for node in G.nodes()]
    
    # Calculate edge widths based on interaction strength
    edge_weights = []
    for edge in G.edges():
        edge_data = G[edge[0]][edge[1]]
        if method_name == 'CONTACT':
            weight = edge_data.get('jaccard_percent', 1)
        elif method_name == 'PISA':
            weight = abs(edge_data.get('free_energy', 1))
        elif method_name == 'EPPIC':
            weight = abs(edge_data.get('cs_score', 1))
        else:
            weight = 1
        edge_weights.append(max(0.5, weight / 10))
    
    # Draw the network
    nx.draw_networkx_nodes(G, pos, 
                          node_size=node_sizes,
                          node_color='lightblue',
                          alpha=0.8,
                          ax=ax)
    
    nx.draw_networkx_edges(G, pos,
                          width=edge_weights,
                          alpha=0.6,
                          edge_color='gray',
                          ax=ax)
    
    # Add labels for high-degree nodes only
    high_degree_nodes = [node for node, degree in degrees.items() if degree > 3]
    labels = {node: node.split('_')[-1] if len(node.split('_')) > 1 else node 
              for node in high_degree_nodes}
    nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)
    
    # Set title and labels
    title = f"{method_name} Network - PDB {complex_info['PDBID']}\n"
    title += f"Interactions: {complex_info['num_interactions']} | "
    title += f"Proteins: {complex_info['num_proteins']} | "
    title += f"Exons: {complex_info['num_exons']}"
    
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.axis('off')
    
    # Add legend
    legend_elements = [
        plt.scatter([], [], s=100, c='lightblue', alpha=0.8, label='Exon'),
        plt.Line2D([0], [0], color='gray', alpha=0.6, linewidth=2, label='EEI')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    
    # Save the plot
    filename = f"{method_name}_network_{complex_info['PDBID']}_{layout_type}.png"
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved network visualization: {filename}")
    
    plt.show()
    
    return fig

def create_comparative_network_plot(complexes_data, output_dir):
    """Create a comparative plot showing top complexes from all methods."""
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    
    methods = ['CONTACT', 'PISA', 'EPPIC']
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
    
    for i, (method, color) in enumerate(zip(methods, colors)):
        ax = axes[i]
        complexes = complexes_data[method]
        
        if not complexes:
            ax.text(0.5, 0.5, f'No data for {method}', 
                   ha='center', va='center', transform=ax.transAxes)
            continue
        
        # Create bar plot of top 5 complexes
        top_5 = complexes[:5]
        pdb_ids = [comp['PDBID'] for comp in top_5]
        interactions = [comp['num_interactions'] for comp in top_5]
        
        bars = ax.bar(range(len(pdb_ids)), interactions, color=color, alpha=0.7)
        
        # Add value labels on bars
        for bar, value in zip(bars, interactions):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                   str(value), ha='center', va='bottom', fontweight='bold')
        
        ax.set_title(f'{method} - Top Protein Complexes', fontweight='bold')
        ax.set_xlabel('PDB ID')
        ax.set_ylabel('Number of Interactions')
        ax.set_xticks(range(len(pdb_ids)))
        ax.set_xticklabels(pdb_ids, rotation=45)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save comparative plot
    output_path = os.path.join(output_dir, 'comparative_complexes.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved comparative plot: {output_path}")
    
    plt.show()

def analyze_network_properties(complexes_data, output_dir):
    """Analyze and visualize network properties."""
    
    print("\nAnalyzing network properties...")
    
    # Collect statistics
    all_stats = []
    for method, complexes in complexes_data.items():
        for complex_info in complexes:
            G = create_network_graph(complex_info['data'], method, complex_info['PDBID'], 
                                   complex_info.get('protein2_col', 'protein2'))
            
            # Calculate network metrics
            stats = {
                'method': method,
                'pdb_id': complex_info['PDBID'],
                'num_nodes': G.number_of_nodes(),
                'num_edges': G.number_of_edges(),
                'density': nx.density(G),
                'avg_clustering': nx.average_clustering(G),
                'avg_shortest_path': nx.average_shortest_path_length(G) if nx.is_connected(G) else np.nan
            }
            all_stats.append(stats)
    
    stats_df = pd.DataFrame(all_stats)
    
    # Create comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Network density comparison
    sns.boxplot(data=stats_df, x='method', y='density', ax=axes[0,0])
    axes[0,0].set_title('Network Density by Method')
    axes[0,0].set_ylabel('Density')
    
    # Average clustering comparison
    sns.boxplot(data=stats_df, x='method', y='avg_clustering', ax=axes[0,1])
    axes[0,1].set_title('Average Clustering by Method')
    axes[0,1].set_ylabel('Average Clustering')
    
    # Number of edges vs nodes
    sns.scatterplot(data=stats_df, x='num_nodes', y='num_edges', hue='method', ax=axes[1,0])
    axes[1,0].set_title('Network Size: Nodes vs Edges')
    axes[1,0].set_xlabel('Number of Nodes')
    axes[1,0].set_ylabel('Number of Edges')
    
    # Network size distribution
    sns.histplot(data=stats_df, x='num_edges', hue='method', alpha=0.7, ax=axes[1,1])
    axes[1,1].set_title('Distribution of Network Sizes')
    axes[1,1].set_xlabel('Number of Edges')
    
    plt.tight_layout()
    
    # Save network properties plot
    output_path = os.path.join(output_dir, 'network_properties.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved network properties analysis: {output_path}")
    
    plt.show()
    
    # Save statistics
    stats_path = os.path.join(output_dir, 'network_statistics.csv')
    stats_df.to_csv(stats_path, index=False)
    print(f"Saved network statistics: {stats_path}")
    
    return stats_df

def main():
    """Main function to run the network visualization analysis."""
    
    # Define paths
    base_dir = "/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/EEI-Homo-Sapiens-2024/results"
    contact_file = os.path.join(base_dir, "1-CONTACT_results", "human_network_final.txt")
    pisa_file = os.path.join(base_dir, "2-PISA_results", "PISA_EEIN_0.5_human.txt")
    eppic_file = os.path.join(base_dir, "3-EPPIC_network", "EPPIC_EEIN_filtered.txt")
    
    # Create output directory
    output_dir = "/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/final_statistics/network_visualizations"
    os.makedirs(output_dir, exist_ok=True)
    
    print("EEI Network Visualization Analysis")
    print("="*50)
    
    # Load data
    contact_df = load_contact_data(contact_file)
    pisa_df = load_pisa_data(pisa_file)
    eppic_df = load_eppic_data(eppic_file)
    
    # Identify representative complexes
    contact_complexes = identify_protein_complexes(contact_df, 'CONTACT')
    pisa_complexes = identify_protein_complexes(pisa_df, 'PISA')
    eppic_complexes = identify_protein_complexes(eppic_df, 'EPPIC')
    
    complexes_data = {
        'CONTACT': contact_complexes,
        'PISA': pisa_complexes,
        'EPPIC': eppic_complexes
    }
    
    # Create comparative plot
    create_comparative_network_plot(complexes_data, output_dir)
    
    # Analyze network properties
    stats_df = analyze_network_properties(complexes_data, output_dir)
    
    # Visualize top 3 complexes from each method
    print("\nCreating detailed network visualizations...")
    
    for method, complexes in complexes_data.items():
        print(f"\nVisualizing {method} networks...")
        
        # Select top 3 complexes for visualization
        top_complexes = complexes[:3]
        
        for i, complex_info in enumerate(top_complexes):
            print(f"  Creating visualization for {complex_info['PDBID']}...")
            
            # Create network graph
            G = create_network_graph(complex_info['data'], method, complex_info['PDBID'],
                                   complex_info.get('protein2_col', 'protein2'))
            
            # Create visualizations with different layouts
            for layout in ['spring', 'circular', 'kamada_kawai']:
                visualize_network(G, complex_info, method, output_dir, layout)
    
    print(f"\nAnalysis complete! All visualizations saved to: {output_dir}")
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for method, complexes in complexes_data.items():
        print(f"{method}: {len(complexes)} protein complexes identified")
        if complexes:
            top_complex = complexes[0]
            print(f"  Top complex: PDB {top_complex['PDBID']} with {top_complex['num_interactions']} interactions")

if __name__ == "__main__":
    main()
