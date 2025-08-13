#!/usr/bin/env python3
"""
Script to plot and visualize the results of the EEI survival analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('default')
sns.set_palette("husl")

def load_results(results_dir="results"):
    """Load all results from the analysis."""
    
    # Load main results
    all_results = pd.read_csv(f"{results_dir}/all_eei_survival_results.tsv", sep='\t')
    significant_results = pd.read_csv(f"{results_dir}/significant_eei_survival_results.tsv", sep='\t')
    
    # Load summary statistics
    with open(f"{results_dir}/survival_analysis_summary.json", 'r') as f:
        summary_stats = json.load(f)
    
    return all_results, significant_results, summary_stats

def plot_p_value_distribution(all_results, significant_results, output_dir="results/plots"):
    """Plot p-value distribution."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    plt.figure(figsize=(12, 8))
    
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: P-value distribution
    ax1.hist(all_results['p_value'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.axvline(x=0.05, color='red', linestyle='--', label='p=0.05 threshold')
    ax1.set_xlabel('P-value')
    ax1.set_ylabel('Frequency')
    ax1.set_title('P-value Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Log-transformed p-values
    log_pvals = -np.log10(all_results['p_value'])
    ax2.hist(log_pvals, bins=20, alpha=0.7, color='lightgreen', edgecolor='black')
    ax2.axvline(x=-np.log10(0.05), color='red', linestyle='--', label='p=0.05 threshold')
    ax2.set_xlabel('-log10(P-value)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Log-transformed P-value Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/p_value_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ P-value distribution plot saved to {output_dir}/p_value_distribution.png")

def plot_significance_summary(significant_results, summary_stats, output_dir="results/plots"):
    """Plot summary of significant results."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    plt.figure(figsize=(15, 10))
    
    # Create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Prognostic types
    if len(significant_results) > 0:
        prognostic_counts = significant_results['prognostic_type'].value_counts()
        colors = ['green' if x == 'Favorable' else 'red' for x in prognostic_counts.index]
        ax1.bar(prognostic_counts.index, prognostic_counts.values, color=colors, alpha=0.7)
        ax1.set_title('Prognostic Types of Significant EEIs')
        ax1.set_ylabel('Count')
        for i, v in enumerate(prognostic_counts.values):
            ax1.text(i, v + 0.1, str(v), ha='center', va='bottom')
    
    # Plot 2: P-values of significant EEIs
    if len(significant_results) > 0:
        significant_results_sorted = significant_results.sort_values('p_value')
        colors = ['green' if x == 'Favorable' else 'red' for x in significant_results_sorted['prognostic_type']]
        bars = ax2.bar(range(len(significant_results_sorted)), -np.log10(significant_results_sorted['p_value']), 
                       color=colors, alpha=0.7)
        ax2.set_xlabel('EEI Index')
        ax2.set_ylabel('-log10(P-value)')
        ax2.set_title('Significance of Significant EEIs')
        ax2.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05 threshold')
        ax2.legend()
        
        # Add EEI labels
        for i, (idx, row) in enumerate(significant_results_sorted.iterrows()):
            ax2.text(i, -np.log10(row['p_value']) + 0.1, 
                    f"{row['mouse_exon1'].split(':')[1][:6]}...", 
                    rotation=45, ha='right', va='bottom', fontsize=8)
    
    # Plot 3: Sample sizes
    if len(significant_results) > 0:
        ax3.scatter(significant_results['n_present'], significant_results['n_absent'], 
                   c=significant_results['p_value'], cmap='viridis', s=100, alpha=0.7)
        ax3.set_xlabel('Samples with EEI Present')
        ax3.set_ylabel('Samples with EEI Absent')
        ax3.set_title('Sample Sizes vs P-value')
        plt.colorbar(ax3.collections[0], ax=ax3, label='P-value')
    
    # Plot 4: Identity vs significance
    if len(significant_results) > 0:
        ax4.scatter(significant_results['avg_identity'], -np.log10(significant_results['p_value']), 
                   c=[1 if x == 'Favorable' else 0 for x in significant_results['prognostic_type']], 
                   cmap='RdYlGn', s=100, alpha=0.7)
        ax4.set_xlabel('Average Identity')
        ax4.set_ylabel('-log10(P-value)')
        ax4.set_title('Identity vs Significance')
        ax4.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/significance_summary.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Significance summary plot saved to {output_dir}/significance_summary.png")

def plot_survival_curves(significant_results, expression_file, survival_file, output_dir="results/plots"):
    """Plot Kaplan-Meier survival curves for significant EEIs."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    expression_df = pd.read_csv(expression_file, sep='\t', index_col=0)
    survival_df = pd.read_csv(survival_file, index_col=0)
    
    # Align samples
    common_samples = set(expression_df.columns) & set(survival_df.index)
    expression_df = expression_df[list(common_samples)]
    survival_df = survival_df.loc[list(common_samples)]
    
    # Plot each significant EEI
    for idx, eei_row in significant_results.iterrows():
        try:
            # Get coordinates
            exon1_coord = eei_row['mouse_exon1']
            exon2_coord = eei_row['mouse_exon2']
            
            # Map to genes
            from utils.coord_gene_mapping import coordinate_to_gene_mapping_with_files
            gene1 = coordinate_to_gene_mapping_with_files(exon1_coord, "coord_gene_mapping_files")
            gene2 = coordinate_to_gene_mapping_with_files(exon2_coord, "coord_gene_mapping_files")
            
            if gene1 and gene2 and gene1 in expression_df.index and gene2 in expression_df.index:
                # Get expression data
                gene1_expr = expression_df.loc[gene1]
                gene2_expr = expression_df.loc[gene2]
                
                # Define EEI presence (both genes expressed > threshold)
                threshold = 0.1  # CPM threshold
                eei_present = (gene1_expr > threshold) & (gene2_expr > threshold)
                
                # Create groups
                present_samples = eei_present[eei_present == True].index.tolist()
                absent_samples = eei_present[eei_present == False].index.tolist()
                
                if len(present_samples) >= 2 and len(absent_samples) >= 2:
                    # Get survival data for each group
                    present_survival = survival_df.loc[present_samples]
                    absent_survival = survival_df.loc[absent_samples]
                    
                    # Create plot
                    plt.figure(figsize=(10, 6))
                    
                    # Fit Kaplan-Meier curves
                    kmf_present = KaplanMeierFitter()
                    kmf_absent = KaplanMeierFitter()
                    
                    # Use 'survival_days' instead of 'survival_time'
                    kmf_present.fit(present_survival['survival_days'], present_survival['event_status'], 
                                  label=f'EEI Present (n={len(present_samples)})')
                    kmf_absent.fit(absent_survival['survival_days'], absent_survival['event_status'], 
                                 label=f'EEI Absent (n={len(absent_samples)})')
                    
                    # Plot curves
                    kmf_present.plot(ci_show=True, color='green')
                    kmf_absent.plot(ci_show=True, color='red')
                    
                    # Perform log-rank test
                    logrank_result = logrank_test(
                        present_survival['survival_days'], absent_survival['survival_days'],
                        present_survival['event_status'], absent_survival['event_status']
                    )
                    
                    plt.title(f'EEI: {gene1}-{gene2}\nP-value: {eei_row["p_value"]:.3e} ({eei_row["prognostic_type"]})')
                    plt.xlabel('Time (days)')
                    plt.ylabel('Survival Probability')
                    plt.grid(True, alpha=0.3)
                    
                    # Add median survival times
                    median_present = kmf_present.median_survival_time_
                    median_absent = kmf_absent.median_survival_time_
                    plt.text(0.05, 0.95, f'Median Present: {median_present:.1f} days\nMedian Absent: {median_absent:.1f} days', 
                            transform=plt.gca().transAxes, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                    
                    # Save plot
                    safe_name = f"{gene1}_{gene2}".replace(':', '_').replace('-', '_')
                    plt.savefig(f"{output_dir}/survival_curve_{safe_name}.png", dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    print(f"✓ Survival curve saved for {gene1}-{gene2}")
                    
        except Exception as e:
            print(f"Error plotting EEI {exon1_coord}-{exon2_coord}: {e}")

def plot_analysis_overview(all_results, significant_results, summary_stats, output_dir="results/plots"):
    """Create an overview plot of the entire analysis."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Analysis summary
    total_eeis = summary_stats['total_eeis_tested']
    significant_eeis = summary_stats['significant_eeis']
    tested_eeis = len(all_results)
    
    labels = ['Tested EEIs', 'Significant EEIs', 'Non-significant']
    sizes = [tested_eeis, significant_eeis, tested_eeis - significant_eeis]
    colors = ['lightblue', 'red', 'lightgray']
    
    ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax1.set_title('EEI Analysis Overview')
    
    # Plot 2: Significance rate
    if len(all_results) > 0:
        p_values = all_results['p_value']
        significant_pct = (p_values <= 0.05).sum() / len(p_values) * 100
        ax2.bar(['Significant (p≤0.05)', 'Non-significant'], 
                [significant_pct, 100-significant_pct], 
                color=['red', 'lightgray'], alpha=0.7)
        ax2.set_ylabel('Percentage')
        ax2.set_title('Significance Rate')
        ax2.text(0, significant_pct + 1, f'{significant_pct:.1f}%', ha='center', va='bottom')
    
    # Plot 3: Identity distribution
    if len(all_results) > 0:
        ax3.hist(all_results['avg_identity'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax3.set_xlabel('Average Identity')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Identity Distribution')
        ax3.axvline(all_results['avg_identity'].mean(), color='red', linestyle='--', 
                    label=f'Mean: {all_results["avg_identity"].mean():.3f}')
        ax3.legend()
    
    # Plot 4: Sample size distribution
    if len(all_results) > 0:
        total_samples = all_results['n_present'] + all_results['n_absent']
        ax4.hist(total_samples, bins=15, alpha=0.7, color='lightgreen', edgecolor='black')
        ax4.set_xlabel('Total Samples')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Sample Size Distribution')
        ax4.axvline(total_samples.mean(), color='red', linestyle='--', 
                    label=f'Mean: {total_samples.mean():.1f}')
        ax4.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/analysis_overview.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Analysis overview plot saved to {output_dir}/analysis_overview.png")

def create_results_report(significant_results, summary_stats, output_dir="results"):
    """Create a text report of the results."""
    
    report_file = f"{output_dir}/analysis_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("EEI SURVIVAL ANALYSIS REPORT\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("SUMMARY STATISTICS:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Total EEIs tested: {summary_stats['total_eeis_tested']}\n")
        f.write(f"Significant EEIs (p≤0.05): {summary_stats['significant_eeis']}\n")
        f.write(f"Significance rate: {summary_stats['significance_rate']:.1%}\n")
        f.write(f"Median p-value: {summary_stats['median_p_value']:.4f}\n")
        f.write(f"Minimum p-value: {summary_stats['min_p_value']:.2e}\n\n")
        
        if summary_stats['significant_eeis'] > 0:
            f.write("PROGNOSTIC TYPES:\n")
            f.write("-" * 30 + "\n")
            f.write(f"Favorable EEIs: {summary_stats['favorable_eeis']}\n")
            f.write(f"Unfavorable EEIs: {summary_stats['unfavorable_eeis']}\n")
            f.write(f"Average identity of significant EEIs: {summary_stats.get('avg_identity_significant', 0):.3f}\n\n")
            
            f.write("TOP SIGNIFICANT EEIs:\n")
            f.write("-" * 30 + "\n")
            for idx, row in significant_results.head().iterrows():
                f.write(f"• {row['mouse_exon1']} - {row['mouse_exon2']}\n")
                f.write(f"  P-value: {row['p_value']:.2e}\n")
                f.write(f"  Prognostic type: {row['prognostic_type']}\n")
                f.write(f"  Average identity: {row['avg_identity']:.3f}\n")
                f.write(f"  Sample sizes: {row['n_present']} present, {row['n_absent']} absent\n\n")
        
        f.write("=" * 60 + "\n")
        f.write("Report generated by plot_survival_results.py\n")
        f.write("=" * 60 + "\n")
    
    print(f"✓ Analysis report saved to {report_file}")

def main():
    """Main function to generate all plots."""
    
    print("=== EEI SURVIVAL ANALYSIS PLOTTING ===")
    
    # Load results
    try:
        all_results, significant_results, summary_stats = load_results()
        print(f"✓ Loaded {len(all_results)} total results")
        print(f"✓ Loaded {len(significant_results)} significant results")
    except Exception as e:
        print(f"Error loading results: {e}")
        return
    
    # Create plots directory
    plots_dir = "results/plots"
    os.makedirs(plots_dir, exist_ok=True)
    
    # Generate plots
    print("\nGenerating plots...")
    
    # 1. P-value distribution
    plot_p_value_distribution(all_results, significant_results, plots_dir)
    
    # 2. Significance summary
    if len(significant_results) > 0:
        plot_significance_summary(significant_results, summary_stats, plots_dir)
    
    # 3. Survival curves for significant EEIs
    if len(significant_results) > 0:
        plot_survival_curves(significant_results, "geo_data/mouse_expression_transposed.tsv", 
                           "outputs/survival_df.csv", plots_dir)
    
    # 4. Analysis overview
    plot_analysis_overview(all_results, significant_results, summary_stats, plots_dir)
    
    # 5. Create report
    create_results_report(significant_results, summary_stats)
    
    print(f"\n=== PLOTTING COMPLETE ===")
    print(f"All plots saved to: {plots_dir}")
    print(f"Report saved to: results/analysis_report.txt")

if __name__ == "__main__":
    main() 