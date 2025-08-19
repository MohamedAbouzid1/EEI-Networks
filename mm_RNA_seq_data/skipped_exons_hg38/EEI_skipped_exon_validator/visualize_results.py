#!/usr/bin/env python3
"""
Generate visualizations for EEI validation results
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

def create_validation_summary_plot(stats_df, output_dir):
    """
    Create a comprehensive summary figure
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('EEI Network Validation Summary', fontsize=16, fontweight='bold')
    
    # 1. Validation rates by method
    ax1 = axes[0, 0]
    methods = stats_df['method']
    validation_rates = stats_df['percentage_eeis_validated']
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
    bars = ax1.bar(methods, validation_rates, color=colors, alpha=0.8)
    ax1.set_ylabel('Validation Rate (%)', fontweight='bold')
    ax1.set_title('A. EEI Validation Rates by Method')
    ax1.set_ylim(0, 15)
    
    # Add value labels on bars
    for bar, val in zip(bars, validation_rates):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # 2. Exon matching rates
    ax2 = axes[0, 1]
    exon_rates = stats_df['percentage_exons_matched']
    bars2 = ax2.bar(methods, exon_rates, color=colors, alpha=0.8)
    ax2.set_ylabel('Exons Matched (%)', fontweight='bold')
    ax2.set_title('B. Exon Matching Rates by Method')
    ax2.set_ylim(0, 30)
    
    for bar, val in zip(bars2, exon_rates):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # 3. Confidence score comparison
    ax3 = axes[0, 2]
    x = np.arange(len(methods))
    width = 0.35
    
    conf_all = stats_df['mean_confidence_all']
    conf_matched = stats_df['mean_confidence_matched']
    
    bars1 = ax3.bar(x - width/2, conf_all, width, label='All EEIs', alpha=0.6)
    bars2 = ax3.bar(x + width/2, conf_matched, width, label='Validated EEIs', alpha=0.8)
    
    ax3.set_ylabel('Mean Confidence Score', fontweight='bold')
    ax3.set_title('C. Confidence Score Comparison')
    ax3.set_xticks(x)
    ax3.set_xticklabels(methods)
    ax3.legend()
    ax3.set_ylim(0, 1)
    
    # 4. Number of validated events
    ax4 = axes[1, 0]
    unique_events = stats_df['unique_skipped_events']
    bars4 = ax4.bar(methods, unique_events, color=colors, alpha=0.8)
    ax4.set_ylabel('Unique Skipped Events', fontweight='bold')
    ax4.set_title('D. Unique Splicing Events Captured')
    
    for bar, val in zip(bars4, unique_events):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(val)}', ha='center', va='bottom', fontweight='bold')
    
    # 5. Total EEIs vs Validated
    ax5 = axes[1, 1]
    total_eeis = stats_df['total_eeis']
    validated_eeis = stats_df['eeis_with_skipped_exons']
    
    x = np.arange(len(methods))
    width = 0.35
    
    bars1 = ax5.bar(x - width/2, total_eeis, width, label='Total EEIs', alpha=0.6)
    bars2 = ax5.bar(x + width/2, validated_eeis, width, label='Validated EEIs', alpha=0.8)
    
    ax5.set_ylabel('Number of EEIs', fontweight='bold')
    ax5.set_title('E. Total vs Validated EEIs')
    ax5.set_xticks(x)
    ax5.set_xticklabels(methods)
    ax5.legend()
    
    # 6. Confidence improvement
    ax6 = axes[1, 2]
    conf_diff = stats_df['confidence_diff']
    bars6 = ax6.bar(methods, conf_diff * 100, color=colors, alpha=0.8)
    ax6.set_ylabel('Confidence Improvement (%)', fontweight='bold')
    ax6.set_title('F. Confidence Score Improvement\nin Validated EEIs')
    ax6.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    for bar, val in zip(bars6, conf_diff * 100):
        height = bar.get_height()
        ax6.text(bar.get_x() + bar.get_width()/2., height,
                f'+{val:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    
    # Save figure
    output_file = output_dir / 'validation_summary.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'validation_summary.pdf', bbox_inches='tight')
    print(f"Saved summary figure to {output_file}")
    
    return fig

def create_psi_distribution_plot(validation_results, output_dir):
    """
    Create PSI distribution plot for validated exons
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('PSI Distribution of Validated Exons', fontsize=16, fontweight='bold')
    
    methods = ['PISA', 'EPPIC', 'Contact']
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
    
    for idx, (method, color) in enumerate(zip(methods, colors)):
        ax = axes[idx]
        
        # Load validation results
        file_path = output_dir.parent / f"{method.lower()}_validation_results.tsv"
        if file_path.exists():
            df = pd.read_csv(file_path, sep='\t')
            
            # Filter for entries with PSI data
            psi_data = df[df['mean_psi'].notna()]
            
            if len(psi_data) > 0:
                # Create violin plot
                parts = ax.violinplot([psi_data['mean_psi']], positions=[0.5], 
                                      widths=0.7, showmeans=True, showmedians=True)
                
                for pc in parts['bodies']:
                    pc.set_facecolor(color)
                    pc.set_alpha(0.7)
                
                # Add scatter points
                y = psi_data['mean_psi']
                x = np.random.normal(0.5, 0.04, size=len(y))
                ax.scatter(x, y, color=color, alpha=0.5, s=50)
                
                # Statistics
                mean_psi = psi_data['mean_psi'].mean()
                median_psi = psi_data['mean_psi'].median()
                
                ax.text(0.5, -0.15, f'Mean: {mean_psi:.3f}\nMedian: {median_psi:.3f}',
                       ha='center', transform=ax.transData, fontsize=10)
            
            ax.set_xlim(0, 1)
            ax.set_ylim(-0.05, 1.05)
            ax.set_xticks([0.5])
            ax.set_xticklabels([method])
            ax.set_ylabel('Mean PSI Value' if idx == 0 else '')
            ax.set_title(f'{method} Method\n(n={len(psi_data)} with PSI data)')
            ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
            ax.axhline(y=0.9, color='red', linestyle='--', alpha=0.3)
            
            # Add text annotations
            ax.text(0.95, 0.92, 'Constitutive\n(PSI > 0.9)', 
                   transform=ax.transAxes, fontsize=8, ha='right', color='red')
    
    plt.tight_layout()
    
    # Save figure
    output_file = output_dir / 'psi_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'psi_distribution.pdf', bbox_inches='tight')
    print(f"Saved PSI distribution to {output_file}")
    
    return fig

def create_overlap_venn_diagram(output_dir):
    """
    Create a Venn diagram showing overlap between methods
    """
    from matplotlib_venn import venn3, venn2
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    
    # Load validation results for each method
    methods = ['pisa', 'eppic', 'contact']
    validated_eeis = {}
    
    for method in methods:
        file_path = output_dir.parent / f"{method}_validation_results.tsv"
        if file_path.exists():
            df = pd.read_csv(file_path, sep='\t')
            # Create unique EEI identifiers
            validated_eeis[method] = set(df[['eei_exon1', 'eei_exon2']].apply(
                lambda x: f"{x['eei_exon1']}_{x['eei_exon2']}", axis=1))
    
    if len(validated_eeis) == 3:
        # Create Venn diagram
        venn = venn3([validated_eeis['pisa'], validated_eeis['eppic'], validated_eeis['contact']],
                     ('PISA', 'EPPIC', 'Contact'),
                     set_colors=('#FF6B6B', '#4ECDC4', '#45B7D1'),
                     alpha=0.5)
        
        # Customize the diagram
        for text in venn.set_labels:
            text.set_fontsize(14)
            text.set_fontweight('bold')
        
        for text in venn.subset_labels:
            if text:
                text.set_fontsize(12)
    
    ax.set_title('Overlap of Validated EEIs Between Methods', 
                fontsize=16, fontweight='bold', pad=20)
    
    # Save figure
    output_file = output_dir / 'method_overlap_venn.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved Venn diagram to {output_file}")
    
    return fig

def generate_summary_table(stats_df, output_dir):
    """
    Generate a formatted summary table for presentation
    """
    # Create presentation-ready table
    summary = {
        'Method': stats_df['method'].tolist(),
        'Total EEIs': stats_df['total_eeis'].tolist(),
        'Validated EEIs': [f"{int(x)} ({y:.1f}%)" for x, y in 
                          zip(stats_df['eeis_with_skipped_exons'], 
                              stats_df['percentage_eeis_validated'])],
        'Unique Exons Matched': [f"{int(x)} ({y:.1f}%)" for x, y in 
                                zip(stats_df['unique_exons_matched'],
                                    stats_df['percentage_exons_matched'])],
        'Mean Confidence (All)': [f"{x:.3f}" for x in stats_df['mean_confidence_all']],
        'Mean Confidence (Validated)': [f"{x:.3f}" for x in stats_df['mean_confidence_matched']],
        'Improvement': [f"+{x:.1f}%" for x in stats_df['confidence_diff'] * 100]
    }
    
    summary_df = pd.DataFrame(summary)
    
    # Save as formatted table
    output_file = output_dir / 'summary_table.csv'
    summary_df.to_csv(output_file, index=False)
    print(f"Saved summary table to {output_file}")
    
    # Also save as LaTeX for paper
    latex_file = output_dir / 'summary_table.tex'
    summary_df.to_latex(latex_file, index=False, caption='Validation of predicted EEI networks using alternative splicing data')
    print(f"Saved LaTeX table to {latex_file}")
    
    return summary_df

def main():
    """
    Generate all visualizations
    """
    # Set paths
    base_dir = Path("/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/mm_RNA_seq_data/skipped_exons_hg38")
    output_dir = base_dir / "validation_results" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load summary statistics
    stats_file = base_dir / "validation_results" / "validation_statistics_summary.tsv"
    stats_df = pd.read_csv(stats_file, sep='\t')
    
    print("Generating visualizations...")
    print("=" * 60)
    
    # Generate all plots
    create_validation_summary_plot(stats_df, output_dir)
    create_psi_distribution_plot(None, output_dir)
    
    try:
        create_overlap_venn_diagram(output_dir)
    except ImportError:
        print("Note: Install matplotlib-venn for Venn diagram: pip install matplotlib-venn")
    
    summary_table = generate_summary_table(stats_df, output_dir)
    
    print("\n" + "=" * 60)
    print("SUMMARY TABLE")
    print("=" * 60)
    print(summary_table.to_string(index=False))
    
    print("\n" + "=" * 60)
    print("All visualizations saved to:", output_dir)
    print("=" * 60)

if __name__ == "__main__":
    main()