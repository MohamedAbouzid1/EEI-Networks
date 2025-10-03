import matplotlib.pyplot as plt
import numpy as np
import os


def generate_validation_plots(results_df, output_dir):
    """
    Generate plots validating the approach using paper's findings.
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Response scores distribution
    ax = axes[0, 0]
    if 'response_score' in results_df.columns:
        results_df['response_score'].hist(bins=30, ax=ax)
        ax.set_xlabel('Response Score (Good - Poor)')
        ax.set_ylabel('Number of EEIs')
        ax.set_title('Distribution of EEI Response Scores')
        ax.axvline(x=0, color='red', linestyle='--', alpha=0.5)
    
    # Plot 2: P-value distribution (prefer adjusted)
    ax = axes[0, 1]
    if 'p_adj_fdr_bh' in results_df.columns and results_df['p_adj_fdr_bh'].notna().any():
        (-np.log10(results_df['p_adj_fdr_bh'].dropna())).hist(bins=30, ax=ax)
        ax.set_xlabel('-log10(FDR-adjusted p)')
        ax.set_ylabel('Number of EEIs')
        ax.set_title('Significance of EEI-Response Association (FDR)')
        ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='FDR=0.05')
        ax.legend()
    elif 'p_value' in results_df.columns:
        (-np.log10(results_df['p_value'].dropna())).hist(bins=30, ax=ax)
        ax.set_xlabel('-log10(p-value)')
        ax.set_ylabel('Number of EEIs')
        ax.set_title('Significance of EEI-Response Association')
        ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
        ax.legend()
    
    # Plot 3: Model comparison (4T1 vs E0771)
    ax = axes[1, 0]
    if all(col in results_df.columns for col in ['good_responder_mean', 'moderate_responder_mean']):
        ax.scatter(results_df['moderate_responder_mean'], 
                  results_df['good_responder_mean'],
                  alpha=0.5)
        ax.set_xlabel('4T1+CPA Expression (Moderate Response)')
        ax.set_ylabel('E0771+CPA Expression (Good Response)')
        ax.set_title('EEI Expression: E0771 vs 4T1 Response')
        
        # Add diagonal line
        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),
            np.max([ax.get_xlim(), ax.get_ylim()]),
        ]
        ax.plot(lims, lims, 'k-', alpha=0.3, zorder=0)
    
    # Plot 4: Top significant EEIs (prefer adjusted)
    ax = axes[1, 1]
    if 'p_adj_fdr_bh' in results_df.columns and results_df['p_adj_fdr_bh'].notna().any():
        top_eeis = results_df.nsmallest(20, 'p_adj_fdr_bh')
        if len(top_eeis) > 0:
            y_pos = np.arange(len(top_eeis))
            ax.barh(y_pos, -np.log10(top_eeis['p_adj_fdr_bh']))
            ax.set_yticks(y_pos)
            ax.set_yticklabels([f"EEI_{i+1}" for i in range(len(top_eeis))], fontsize=8)
            ax.set_xlabel('-log10(FDR-adjusted p)')
            ax.set_title('Top 20 Significant EEIs (FDR)')
            ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
    elif 'p_value' in results_df.columns and len(results_df) > 0:
        top_eeis = results_df.nsmallest(20, 'p_value')
        if len(top_eeis) > 0:
            y_pos = np.arange(len(top_eeis))
            ax.barh(y_pos, -np.log10(top_eeis['p_value']))
            ax.set_yticks(y_pos)
            ax.set_yticklabels([f"EEI_{i+1}" for i in range(len(top_eeis))], fontsize=8)
            ax.set_xlabel('-log10(p-value)')
            ax.set_title('Top 20 Significant EEIs')
            ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/eei_response_validation_plots.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Validation plots saved to {output_dir}/eei_response_validation_plots.png")