#!/usr/bin/env python3
"""
Script to create Sankey diagram visualization for orthology relationships 
between human and mouse exons from EGIO output data.

This script parses the EGIO output file and creates interactive Sankey diagrams
showing the different types of orthology relationships between species.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.offline import plot
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style for better plots
plt.style.use('default')
sns.set_palette("husl")

def load_egio_data(file_path):
    """Load and process EGIO output data."""
    print("Loading EGIO data...")
    
    # Read the tab-separated file
    df = pd.read_csv(file_path, sep='\t')
    
    # Remove header row if it exists
    if df.iloc[0, 0] == 'Group':
        df = df.iloc[1:].reset_index(drop=True)
    
    # Rename columns for clarity
    df.columns = ['Group', 'Human_Gene', 'Human_Position', 'Mouse_Gene', 'Mouse_Position', 'Identity', 'Type']
    
    # Convert data types
    df['Group'] = df['Group'].astype(int)
    df['Identity'] = pd.to_numeric(df['Identity'], errors='coerce')
    
    print(f"Loaded {len(df)} orthology relationships")
    print(f"Relationship types found: {df['Type'].value_counts().to_dict()}")
    
    return df

def analyze_orthology_types(df):
    """Analyze the different types of orthology relationships."""
    
    print("\nAnalyzing orthology relationship types...")
    
    # Count different relationship types
    type_counts = df['Type'].value_counts()
    
    relationship_info = {
        '1-1': {
            'description': 'One-to-One orthology',
            'count': type_counts.get('1-1', 0),
            'color': '#2ECC71'  # Green
        },
        '1-0': {
            'description': 'Human-specific (lost in mouse)',
            'count': type_counts.get('1-0', 0),
            'color': '#E74C3C'  # Red
        },
        '0-1': {
            'description': 'Mouse-specific (lost in human)',
            'count': type_counts.get('0-1', 0),
            'color': '#3498DB'  # Blue
        },
        '1-N': {
            'description': 'One-to-Many (human to multiple mouse)',
            'count': type_counts.get('1-N', 0),
            'color': '#F39C12'  # Orange
        },
        'N-1': {
            'description': 'Many-to-One (multiple human to mouse)',
            'count': type_counts.get('N-1', 0),
            'color': '#9B59B6'  # Purple
        }
    }
    
    # Print summary
    total_relationships = len(df)
    print(f"\nOrthology Relationship Summary:")
    print(f"Total relationships: {total_relationships:,}")
    
    for rel_type, info in relationship_info.items():
        if info['count'] > 0:
            percentage = (info['count'] / total_relationships) * 100
            print(f"{rel_type}: {info['count']:,} ({percentage:.1f}%) - {info['description']}")
    
    return relationship_info

def create_sankey_diagram(df, relationship_info, output_dir):
    """Create interactive Sankey diagram using Plotly."""
    
    print("\nCreating Sankey diagram...")
    
    # Prepare data for Sankey diagram
    # We'll create nodes for different relationship types and flows between human and mouse
    
    # Define nodes
    nodes = [
        {'label': 'Human Exons', 'color': '#E74C3C'},
        {'label': 'Mouse Exons', 'color': '#3498DB'},
        {'label': '1:1 Orthologs', 'color': '#2ECC71'},
        {'label': 'Human-Specific', 'color': '#E67E22'},
        {'label': 'Mouse-Specific', 'color': '#2980B9'},
        {'label': '1:N Relations', 'color': '#F39C12'},
        {'label': 'N:1 Relations', 'color': '#9B59B6'}
    ]
    
    # Define links (flows)
    links = []
    type_counts = df['Type'].value_counts()
    
    # 1:1 relationships
    if '1-1' in type_counts:
        count_1_1 = type_counts['1-1']
        links.extend([
            {'source': 0, 'target': 2, 'value': count_1_1, 'color': 'rgba(46, 204, 113, 0.3)'},
            {'source': 2, 'target': 1, 'value': count_1_1, 'color': 'rgba(46, 204, 113, 0.3)'}
        ])
    
    # Human-specific (1-0)
    if '1-0' in type_counts:
        count_1_0 = type_counts['1-0']
        links.extend([
            {'source': 0, 'target': 3, 'value': count_1_0, 'color': 'rgba(231, 76, 60, 0.3)'}
        ])
    
    # Mouse-specific (0-1)
    if '0-1' in type_counts:
        count_0_1 = type_counts['0-1']
        links.extend([
            {'source': 4, 'target': 1, 'value': count_0_1, 'color': 'rgba(52, 152, 219, 0.3)'}
        ])
    
    # 1:N relationships
    if '1-N' in type_counts:
        count_1_N = type_counts['1-N']
        links.extend([
            {'source': 0, 'target': 5, 'value': count_1_N, 'color': 'rgba(243, 156, 18, 0.3)'},
            {'source': 5, 'target': 1, 'value': count_1_N, 'color': 'rgba(243, 156, 18, 0.3)'}
        ])
    
    # N:1 relationships
    if 'N-1' in type_counts:
        count_N_1 = type_counts['N-1']
        links.extend([
            {'source': 0, 'target': 6, 'value': count_N_1, 'color': 'rgba(155, 89, 182, 0.3)'},
            {'source': 6, 'target': 1, 'value': count_N_1, 'color': 'rgba(155, 89, 182, 0.3)'}
        ])
    
    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=[node['label'] for node in nodes],
            color=[node['color'] for node in nodes]
        ),
        link=dict(
            source=[link['source'] for link in links],
            target=[link['target'] for link in links],
            value=[link['value'] for link in links],
            color=[link['color'] for link in links]
        )
    )])
    
    # Update layout
    fig.update_layout(
        title_text="Orthology Relationships Between Human and Mouse Exons<br><sub>Based on EGIO Analysis</sub>",
        title_x=0.5,
        title_font_size=20,
        font_size=12,
        width=1000,
        height=600,
        margin=dict(l=50, r=50, t=80, b=50)
    )
    
    # Save as HTML
    html_path = os.path.join(output_dir, 'sankey_orthology_interactive.html')
    plot(fig, filename=html_path, auto_open=False)
    print(f"Interactive Sankey diagram saved to: {html_path}")
    
    return fig

def create_matplotlib_sankey(df, relationship_info, output_dir):
    """Create static Sankey-like diagram using matplotlib."""
    
    print("\nCreating static Sankey-like diagram...")
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    
    # Define colors
    colors = {
        '1-1': '#2ECC71',  # Green
        '1-0': '#E74C3C',  # Red
        '0-1': '#3498DB',  # Blue
        '1-N': '#F39C12',  # Orange
        'N-1': '#9B59B6'   # Purple
    }
    
    # Calculate counts
    type_counts = df['Type'].value_counts()
    total_count = len(df)
    
    # Create stacked bar chart representation
    y_positions = np.arange(len(relationship_info))
    bar_width = 0.8
    
    # Draw bars for each relationship type
    bottom = 0
    for i, (rel_type, info) in enumerate(relationship_info.items()):
        if info['count'] > 0:
            percentage = (info['count'] / total_count) * 100
            height = percentage
            
            # Draw bar
            ax.barh(i, height, left=0, height=bar_width, 
                   color=info['color'], alpha=0.7, edgecolor='white', linewidth=1)
            
            # Add count label
            ax.text(height/2, i, f"{info['count']:,}\n({percentage:.1f}%)", 
                   ha='center', va='center', fontweight='bold', fontsize=10)
            
            # Add relationship type label
            ax.text(-2, i, f"{rel_type}: {info['description']}", 
                   ha='right', va='center', fontsize=11, fontweight='bold')
    
    # Customize plot
    ax.set_xlim(-25, 100)
    ax.set_ylim(-0.5, len(relationship_info) - 0.5)
    ax.set_xlabel('Percentage of Total Relationships', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of Orthology Relationship Types\nBetween Human and Mouse Exons', 
                fontsize=16, fontweight='bold', pad=20)
    
    # Remove y-axis ticks and labels
    ax.set_yticks([])
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='x')
    ax.set_axisbelow(True)
    
    # Add summary statistics
    summary_text = f"""Summary Statistics:
Total Orthology Relationships: {total_count:,}
One-to-One Orthologs: {type_counts.get('1-1', 0):,} ({type_counts.get('1-1', 0)/total_count*100:.1f}%)
Human-Specific Exons: {type_counts.get('1-0', 0):,} ({type_counts.get('1-0', 0)/total_count*100:.1f}%)
Mouse-Specific Exons: {type_counts.get('0-1', 0):,} ({type_counts.get('0-1', 0)/total_count*100:.1f}%)"""
    
    ax.text(0.98, 0.98, summary_text, transform=ax.transAxes, 
           fontsize=10, va='top', ha='right',
           bbox=dict(boxstyle="round,pad=0.5", facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(output_dir, 'sankey_orthology_static.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Static Sankey diagram saved to: {output_path}")
    
    plt.show()
    
    return fig

def create_identity_distribution_plot(df, output_dir):
    """Create plot showing identity score distribution for different relationship types."""
    
    print("\nCreating identity score distribution plot...")
    
    # Filter out non-numeric identity values and 0-1, 1-0 relationships (no identity)
    df_filtered = df[df['Identity'].notna() & (df['Identity'] > 0)].copy()
    
    if len(df_filtered) == 0:
        print("No identity scores available for plotting")
        return None
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: Box plot of identity scores by relationship type
    relationship_types = df_filtered['Type'].unique()
    box_data = [df_filtered[df_filtered['Type'] == rt]['Identity'].values for rt in relationship_types]
    box_labels = [rt for rt in relationship_types]
    
    bp = ax1.boxplot(box_data, labels=box_labels, patch_artist=True)
    
    # Color the boxes
    colors = ['#2ECC71', '#F39C12', '#9B59B6']  # Green, Orange, Purple
    for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax1.set_title('Identity Score Distribution by Relationship Type', fontweight='bold')
    ax1.set_xlabel('Relationship Type')
    ax1.set_ylabel('Identity Score')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Histogram of all identity scores
    ax2.hist(df_filtered['Identity'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    ax2.axvline(df_filtered['Identity'].mean(), color='red', linestyle='--', 
               label=f'Mean: {df_filtered["Identity"].mean():.3f}')
    ax2.axvline(df_filtered['Identity'].median(), color='green', linestyle='--', 
               label=f'Median: {df_filtered["Identity"].median():.3f}')
    
    ax2.set_title('Distribution of Identity Scores', fontweight='bold')
    ax2.set_xlabel('Identity Score')
    ax2.set_ylabel('Frequency')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(output_dir, 'identity_distribution.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Identity distribution plot saved to: {output_path}")
    
    plt.show()
    
    return fig

def create_gene_level_summary(df, output_dir):
    """Create summary statistics at the gene level."""
    
    print("\nCreating gene-level summary...")
    
    # Count unique genes
    human_genes = df['Human_Gene'].dropna().nunique()
    mouse_genes = df['Mouse_Gene'].dropna().nunique()
    
    # Count genes with different relationship types
    gene_summary = []
    
    for rel_type in ['1-1', '1-0', '0-1', '1-N', 'N-1']:
        subset = df[df['Type'] == rel_type]
        if len(subset) > 0:
            human_gene_count = subset['Human_Gene'].dropna().nunique()
            mouse_gene_count = subset['Mouse_Gene'].dropna().nunique()
            
            gene_summary.append({
                'Relationship_Type': rel_type,
                'Human_Genes': human_gene_count,
                'Mouse_Genes': mouse_gene_count,
                'Total_Relationships': len(subset)
            })
    
    gene_summary_df = pd.DataFrame(gene_summary)
    
    # Save summary
    summary_path = os.path.join(output_dir, 'gene_level_summary.csv')
    gene_summary_df.to_csv(summary_path, index=False)
    print(f"Gene-level summary saved to: {summary_path}")
    
    # Print summary
    print(f"\nGene-Level Summary:")
    print(f"Total unique human genes: {human_genes:,}")
    print(f"Total unique mouse genes: {mouse_genes:,}")
    print(f"\nGene counts by relationship type:")
    print(gene_summary_df.to_string(index=False))
    
    return gene_summary_df

def main():
    """Main function to create orthology relationship visualizations."""
    
    # Define paths
    egio_file = "/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/mm_RNA_seq_data/EGIO/ExonGroup_testpro_hsa_mus.txt"
    
    # Create output directory
    output_dir = "/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/final_statistics/EGIO_analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    print("EGIO Orthology Relationship Analysis")
    print("="*50)
    
    # Load EGIO data
    df = load_egio_data(egio_file)
    
    # Analyze relationship types
    relationship_info = analyze_orthology_types(df)
    
    # Create visualizations
    create_sankey_diagram(df, relationship_info, output_dir)
    create_matplotlib_sankey(df, relationship_info, output_dir)
    create_identity_distribution_plot(df, output_dir)
    
    # Create gene-level summary
    gene_summary_df = create_gene_level_summary(df, output_dir)
    
    print(f"\nAnalysis complete! All visualizations saved to: {output_dir}")
    
    # Print final summary
    print("\n" + "="*60)
    print("FINAL SUMMARY")
    print("="*60)
    total_relationships = len(df)
    print(f"Total orthology relationships analyzed: {total_relationships:,}")
    
    type_counts = df['Type'].value_counts()
    for rel_type in ['1-1', '1-0', '0-1', '1-N', 'N-1']:
        if rel_type in type_counts:
            count = type_counts[rel_type]
            percentage = (count / total_relationships) * 100
            description = relationship_info[rel_type]['description']
            print(f"{rel_type}: {count:,} ({percentage:.1f}%) - {description}")

if __name__ == "__main__":
    main()
