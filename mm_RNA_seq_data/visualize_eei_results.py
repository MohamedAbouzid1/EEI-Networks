#!/usr/bin/env python3
"""
EEI Treatment Association Results Visualization
Creates publication-quality plots for supervisor presentation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def load_results(results_dir):
    """Load results from the specified directory."""
    results_file = Path(results_dir) / "eei_treatment_binary_association_results.tsv"
    significant_file = Path(results_dir) / "significant_eei_treatment_binary_associations.tsv"
    
    # Load all results
    results_df = pd.read_csv(results_file, sep='\t')
    print(f"Loaded {len(results_df)} total EEIs from {results_dir}")
    
    # Load significant results if available
    if significant_file.exists():
        significant_df = pd.read_csv(significant_file, sep='\t')
        print(f"Found {len(significant_df)} significant EEIs")
    else:
        significant_df = pd.DataFrame()
    
    return results_df, significant_df

def create_volcano_plot(results_df, output_dir, title_suffix=""):
    """Create volcano plot showing significance vs effect size."""
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'EEI Treatment Association Volcano Plots{title_suffix}', fontsize=16, fontweight='bold')
    
    # Plot 1: Exon1 t-test
    ax1 = axes[0, 0]
    plot_volcano_subplot(results_df, 'exon1_ttest_p', 'exon1_fold_change', 
                        'Exon 1 Expression Changes', ax1)
    
    # Plot 2: Exon2 t-test
    ax2 = axes[0, 1]
    plot_volcano_subplot(results_df, 'exon2_ttest_p', 'exon2_fold_change', 
                        'Exon 2 Expression Changes', ax2)
    
    # Plot 3: Combined t-test
    ax3 = axes[1, 0]
    plot_volcano_subplot(results_df, 'combined_ttest_p', 'combined_mean_diff', 
                        'Combined Expression Changes', ax3)
    
    # Plot 4: Product t-test
    ax4 = axes[1, 1]
    plot_volcano_subplot(results_df, 'product_ttest_p', 'combined_mean_diff', 
                        'Product Expression Changes', ax4)
    
    plt.tight_layout()
    
    # Save plot
    output_path = Path(output_dir) / f"volcano_plots{title_suffix.replace(' ', '_')}.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved volcano plots to {output_path}")

def plot_volcano_subplot(results_df, p_col, effect_col, title, ax):
    """Create a single volcano plot subplot."""
    
    # Filter out NaN values
    mask = ~(results_df[p_col].isna() | results_df[effect_col].isna())
    data = results_df[mask].copy()
    
    # Calculate -log10(p-value)
    data['neg_log_p'] = -np.log10(data[p_col])
    
    # Define significance threshold
    sig_threshold = -np.log10(0.05)
    
    # Color points by significance
    colors = ['red' if p < 0.05 else 'blue' for p in data[p_col]]
    sizes = [50 if p < 0.05 else 20 for p in data[p_col]]
    
    # Create scatter plot
    ax.scatter(data[effect_col], data['neg_log_p'], 
               c=colors, s=sizes, alpha=0.6, edgecolors='black', linewidth=0.5)
    
    # Add significance line
    ax.axhline(y=sig_threshold, color='red', linestyle='--', alpha=0.7, label='p=0.05')
    
    # Add zero effect line
    if 'fold_change' in effect_col:
        ax.axvline(x=1.0, color='gray', linestyle='--', alpha=0.5, label='No change')
    else:
        ax.axvline(x=0.0, color='gray', linestyle='--', alpha=0.5, label='No change')
    
    ax.set_xlabel('Effect Size')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add counts
    sig_count = (data[p_col] < 0.05).sum()
    total_count = len(data)
    ax.text(0.02, 0.98, f'Significant: {sig_count}/{total_count}', 
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

def create_correlation_analysis_plot(results_df, output_dir, title_suffix=""):
    """Create correlation analysis plots."""
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'EEI Correlation Analysis{title_suffix}', fontsize=16, fontweight='bold')
    
    # Plot 1: Correlation in Control vs Treated
    ax1 = axes[0, 0]
    plot_correlation_comparison(results_df, 'corr_control', 'corr_treated', 
                               'Control', 'Treated', 'Exon Correlation Comparison', ax1)
    
    # Plot 2: Correlation difference distribution
    ax2 = axes[0, 1]
    plot_correlation_difference(results_df, 'corr_difference', ax2)
    
    # Plot 3: Correlation vs significance
    ax3 = axes[1, 0]
    plot_correlation_significance(results_df, ax3)
    
    # Plot 4: Correlation change by test type
    ax4 = axes[1, 1]
    plot_correlation_by_test(results_df, ax4)
    
    plt.tight_layout()
    
    # Save plot
    output_path = Path(output_dir) / f"correlation_analysis{title_suffix.replace(' ', '_')}.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved correlation analysis to {output_path}")

def plot_correlation_comparison(results_df, col1, col2, label1, label2, title, ax):
    """Plot correlation comparison between two groups."""
    
    # Filter out NaN values
    mask = ~(results_df[col1].isna() | results_df[col2].isna())
    data = results_df[mask].copy()
    
    # Create scatter plot
    ax.scatter(data[col1], data[col2], alpha=0.6, s=30, edgecolors='black', linewidth=0.3)
    
    # Add diagonal line
    min_val = min(data[col1].min(), data[col2].min())
    max_val = max(data[col1].max(), data[col2].max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.7, label='y=x')
    
    ax.set_xlabel(f'{label1} Correlation')
    ax.set_ylabel(f'{label2} Correlation')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add correlation coefficient
    corr_coef = data[col1].corr(data[col2])
    ax.text(0.02, 0.98, f'r = {corr_coef:.3f}', 
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

def plot_correlation_difference(results_df, diff_col, ax):
    """Plot distribution of correlation differences."""
    
    # Filter out NaN values
    data = results_df[~results_df[diff_col].isna()][diff_col]
    
    # Create histogram
    ax.hist(data, bins=30, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.axvline(x=data.mean(), color='red', linestyle='--', alpha=0.8, label=f'Mean: {data.mean():.3f}')
    ax.axvline(x=data.median(), color='orange', linestyle='--', alpha=0.8, label=f'Median: {data.median():.3f}')
    
    ax.set_xlabel('Correlation Difference (|Control - Treated|)')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of Correlation Changes')
    ax.legend()
    ax.grid(True, alpha=0.3)

def plot_correlation_significance(results_df, ax):
    """Plot correlation difference vs significance."""
    
    # Use minimum p-value for significance
    mask = ~(results_df['corr_difference'].isna() | results_df['min_p_value'].isna())
    data = results_df[mask].copy()
    
    # Calculate -log10(p-value)
    data['neg_log_p'] = -np.log10(data['min_p_value'])
    
    # Color by significance
    colors = ['red' if p < 0.05 else 'blue' for p in data['min_p_value']]
    sizes = [50 if p < 0.05 else 20 for p in data['min_p_value']]
    
    ax.scatter(data['corr_difference'], data['neg_log_p'], 
               c=colors, s=sizes, alpha=0.6, edgecolors='black', linewidth=0.5)
    
    ax.set_xlabel('Correlation Difference')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Correlation Changes vs Statistical Significance')
    ax.grid(True, alpha=0.3)

def plot_correlation_by_test(results_df, ax):
    """Plot correlation differences by test type."""
    
    # Filter out NaN values
    mask = ~(results_df['corr_difference'].isna() | results_df['best_test'].isna())
    data = results_df[mask].copy()
    
    # Create box plot
    test_order = ['exon1_ttest', 'exon2_ttest', 'combined_ttest', 'product_ttest']
    available_tests = [test for test in test_order if test in data['best_test'].values]
    
    if available_tests:
        plot_data = [data[data['best_test'] == test]['corr_difference'] for test in available_tests]
        ax.boxplot(plot_data, labels=available_tests)
        ax.set_ylabel('Correlation Difference')
        ax.set_title('Correlation Changes by Most Significant Test')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3)

def create_summary_statistics_plot(results_df, output_dir, title_suffix=""):
    """Create summary statistics plots."""
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'EEI Analysis Summary Statistics{title_suffix}', fontsize=16, fontweight='bold')
    
    # Plot 1: P-value distribution
    ax1 = axes[0, 0]
    plot_pvalue_distribution(results_df, ax1)
    
    # Plot 2: Test type success rates
    ax2 = axes[0, 1]
    plot_test_success_rates(results_df, ax2)
    
    # Plot 3: Fold change distribution
    ax3 = axes[1, 0]
    plot_fold_change_distribution(results_df, ax3)
    
    # Plot 4: Effect size summary
    ax4 = axes[1, 1]
    plot_effect_size_summary(results_df, ax4)
    
    plt.tight_layout()
    
    # Save plot
    output_path = Path(output_dir) / f"summary_statistics{title_suffix.replace(' ', '_')}.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved summary statistics to {output_path}")

def plot_pvalue_distribution(results_df, ax):
    """Plot distribution of p-values."""
    
    # Get minimum p-values
    data = results_df['min_p_value'].dropna()
    
    # Create histogram
    ax.hist(data, bins=30, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.axvline(x=0.05, color='red', linestyle='--', alpha=0.8, label='p=0.05')
    
    ax.set_xlabel('P-value')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of P-values')
    ax.set_xscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add counts
    sig_count = (data < 0.05).sum()
    total_count = len(data)
    ax.text(0.02, 0.98, f'Significant: {sig_count}/{total_count}', 
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

def plot_test_success_rates(results_df, ax):
    """Plot success rates of different test types."""
    
    if 'best_test' in results_df.columns:
        test_counts = results_df['best_test'].value_counts()
        
        # Create bar plot
        bars = ax.bar(range(len(test_counts)), test_counts.values, 
                      color=sns.color_palette("husl", len(test_counts)))
        
        ax.set_xlabel('Test Type')
        ax.set_ylabel('Number of EEIs')
        ax.set_title('Most Successful Test Types')
        ax.set_xticks(range(len(test_counts)))
        ax.set_xticklabels(test_counts.index, rotation=45, ha='right')
        
        # Add value labels on bars
        for i, bar in enumerate(bars):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                   f'{int(height)}', ha='center', va='bottom')
        
        ax.grid(True, alpha=0.3)

def plot_fold_change_distribution(results_df, ax):
    """Plot distribution of fold changes."""
    
    # Combine fold changes from both exons
    fold_changes = []
    for col in ['exon1_fold_change', 'exon2_fold_change']:
        if col in results_df.columns:
            data = results_df[col].dropna()
            fold_changes.extend(data)
    
    if fold_changes:
        # Create histogram
        ax.hist(fold_changes, bins=30, alpha=0.7, edgecolor='black', linewidth=0.5)
        ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.8, label='No change')
        ax.axvline(x=np.mean(fold_changes), color='orange', linestyle='--', alpha=0.8, 
                   label=f'Mean: {np.mean(fold_changes):.2f}')
        
        ax.set_xlabel('Fold Change (Treated/Control)')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Expression Fold Changes')
        ax.legend()
        ax.grid(True, alpha=0.3)

def plot_effect_size_summary(results_df, ax):
    """Plot summary of effect sizes."""
    
    # Calculate effect size statistics
    effect_stats = []
    effect_labels = []
    
    for col in ['exon1_mean_diff', 'exon2_mean_diff', 'combined_mean_diff']:
        if col in results_df.columns:
            data = results_df[col].dropna()
            if len(data) > 0:
                effect_stats.append([data.mean(), data.std()])
                effect_labels.append(col.replace('_', ' ').title())
    
    if effect_stats:
        # Create bar plot with error bars
        means = [stat[0] for stat in effect_stats]
        stds = [stat[1] for stat in effect_stats]
        
        bars = ax.bar(range(len(means)), means, yerr=stds, 
                      capsize=5, color=sns.color_palette("husl", len(means)))
        
        ax.set_xlabel('Effect Size Type')
        ax.set_ylabel('Mean Difference (Treated - Control)')
        ax.set_title('Effect Size Summary')
        ax.set_xticks(range(len(means)))
        ax.set_xticklabels(effect_labels, rotation=45, ha='right')
        
        # Add horizontal line at zero
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        
        ax.grid(True, alpha=0.3)

def create_top_hits_table(results_df, output_dir, title_suffix="", top_n=20):
    """Create a summary table of top hits."""
    
    # Sort by significance
    sorted_df = results_df.sort_values('min_p_value').head(top_n)
    
    # Select key columns for display
    display_cols = ['mouse_exon1', 'mouse_exon2', 'min_p_value', 'best_test', 
                    'exon1_fold_change', 'exon2_fold_change', 'corr_difference']
    
    # Filter to available columns
    available_cols = [col for col in display_cols if col in sorted_df.columns]
    display_df = sorted_df[available_cols].copy()
    
    # Round numeric columns
    for col in display_df.columns:
        if display_df[col].dtype in ['float64', 'float32']:
            if 'p_value' in col:
                display_df[col] = display_df[col].apply(lambda x: f"{x:.2e}" if pd.notna(x) else "NA")
            elif 'fold_change' in col:
                display_df[col] = display_df[col].apply(lambda x: f"{x:.2f}" if pd.notna(x) else "NA")
            elif 'corr' in col:
                display_df[col] = display_df[col].apply(lambda x: f"{x:.3f}" if pd.notna(x) else "NA")
    
    # Save table
    output_path = Path(output_dir) / f"top_hits_table{title_suffix.replace(' ', '_')}.tsv"
    display_df.to_csv(output_path, sep='\t', index=False)
    print(f"Saved top hits table to {output_path}")
    
    # Print summary
    print(f"\nTop {top_n} Most Significant EEIs:")
    print("=" * 80)
    print(display_df.to_string(index=False))
    
    return display_df

def main():
    """Main function to generate all visualizations."""
    
    # Define output directory for plots
    output_dir = "eei_quantitative_visualization_plots"
    Path(output_dir).mkdir(exist_ok=True)
    
    # List of result directories to analyze
    result_dirs = [
        "results_treatment_quantitative_high",
        "results_treatment_quantitative_low", 
        "results_treatment_quantitative_medium"
    ]
    
    for result_dir in result_dirs:
        if Path(result_dir).exists():
            print(f"\n{'='*60}")
            print(f"Analyzing results from: {result_dir}")
            print(f"{'='*60}")
            
            # Load results
            results_df, significant_df = load_results(result_dir)
            
            # Create title suffix
            title_suffix = f" - {result_dir.replace('results_treatment_quantitative', '').replace('_', ' ').title()}"
            if not title_suffix.strip(' -'):
                title_suffix = " - All EEIs"
            
            # Generate plots
            try:
                create_volcano_plot(results_df, output_dir, title_suffix)
                create_correlation_analysis_plot(results_df, output_dir, title_suffix)
                create_summary_statistics_plot(results_df, output_dir, title_suffix)
                create_top_hits_table(results_df, output_dir, title_suffix)
                
                print(f"✓ Successfully generated all plots for {result_dir}")
                
            except Exception as e:
                print(f"✗ Error generating plots for {result_dir}: {e}")
                continue
    
    print(f"\n{'='*60}")
    print("VISUALIZATION COMPLETE!")
    print(f"All plots saved to: {output_dir}")
    print(f"{'='*60}")
    
    # Create a summary report
    create_summary_report(result_dirs, output_dir)

def create_summary_report(result_dirs, output_dir):
    """Create a summary report of all analyses."""
    
    report_lines = []
    report_lines.append("EEI TREATMENT ASSOCIATION ANALYSIS - SUMMARY REPORT")
    report_lines.append("=" * 60)
    report_lines.append("")
    
    for result_dir in result_dirs:
        if Path(result_dir).exists():
            try:
                results_df, _ = load_results(result_dir)
                
                # Calculate summary statistics
                total_eeis = len(results_df)
                significant_eeis = (results_df['min_p_value'] < 0.05).sum() if 'min_p_value' in results_df.columns else 0
                significance_rate = (significant_eeis / total_eeis * 100) if total_eeis > 0 else 0
                
                # Get test type distribution
                if 'best_test' in results_df.columns:
                    test_dist = results_df['best_test'].value_counts()
                    most_successful = test_dist.index[0] if len(test_dist) > 0 else "None"
                else:
                    most_successful = "Unknown"
                
                report_lines.append(f"Dataset: {result_dir}")
                report_lines.append(f"  Total EEIs: {total_eeis}")
                report_lines.append(f"  Significant EEIs: {significant_eeis}")
                report_lines.append(f"  Significance Rate: {significance_rate:.1f}%")
                report_lines.append(f"  Most Successful Test: {most_successful}")
                report_lines.append("")
                
            except Exception as e:
                report_lines.append(f"Dataset: {result_dir} - ERROR: {e}")
                report_lines.append("")
    
    # Save report
    report_path = Path(output_dir) / "analysis_summary_report.txt"
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Summary report saved to: {report_path}")

if __name__ == "__main__":
    main()
