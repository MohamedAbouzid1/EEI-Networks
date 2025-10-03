import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def analyze_coverage_confidence_distributions(predicted_eeis_file):
    """
    Analyze the distributions of coverage and confidence scores in predicted EEIs.
    
    Args:
        predicted_eeis_file: Path to the TSV file containing predicted EEIs
    """
    # Read the predicted EEIs
    df = pd.read_csv(predicted_eeis_file, sep='\t')
    
    # Basic statistics for coverage percentages
    coverage_stats = {
        'exon1_coverage_percent': df['exon1_coverage_percent'].describe(),
        'exon2_coverage_percent': df['exon2_coverage_percent'].describe(),
        'jaccard_percent': df['jaccard_percent'].describe()
    }
    
    # Basic statistics for confidence scores
    confidence_stats = df['confidence'].describe()
    
    # Create a figure with multiple subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Distribution of confidence scores
    sns.histplot(data=df, x='confidence', bins=30, ax=ax1)
    ax1.set_title('Distribution of Confidence Scores')
    ax1.set_xlabel('Confidence Score')
    ax1.set_ylabel('Count')
    
    # Plot 2: Distribution of coverage percentages
    sns.histplot(data=df, x='jaccard_percent', bins=30, ax=ax2)
    ax2.set_title('Distribution of Jaccard Coverage Percentages')
    ax2.set_xlabel('Jaccard Coverage (%)')
    ax2.set_ylabel('Count')
    
    # Plot 3: Scatter plot of confidence vs coverage
    sns.scatterplot(data=df, x='confidence', y='jaccard_percent', ax=ax3)
    ax3.set_title('Confidence vs Coverage')
    ax3.set_xlabel('Confidence Score')
    ax3.set_ylabel('Jaccard Coverage (%)')
    
    # Plot 4: Box plot of coverage percentages
    coverage_data = pd.melt(df[['exon1_coverage_percent', 'exon2_coverage_percent', 'jaccard_percent']])
    sns.boxplot(data=coverage_data, x='variable', y='value', ax=ax4)
    ax4.set_title('Coverage Percentages Distribution')
    ax4.set_xlabel('Coverage Type')
    ax4.set_ylabel('Coverage (%)')
    ax4.set_xticklabels(['Exon 1', 'Exon 2', 'Jaccard'])
    
    plt.tight_layout()
    plt.savefig('coverage_confidence_analysis.png')
    plt.close()
    
    # Calculate correlations
    correlations = {
        'confidence_vs_exon1_coverage': df['confidence'].corr(df['exon1_coverage_percent']),
        'confidence_vs_exon2_coverage': df['confidence'].corr(df['exon2_coverage_percent']),
        'confidence_vs_jaccard': df['confidence'].corr(df['jaccard_percent']),
        'exon1_vs_exon2_coverage': df['exon1_coverage_percent'].corr(df['exon2_coverage_percent'])
    }
    
    # Calculate percentiles for coverage
    coverage_percentiles = {
        'exon1_coverage': np.percentile(df['exon1_coverage_percent'], [25, 50, 75, 90, 95]),
        'exon2_coverage': np.percentile(df['exon2_coverage_percent'], [25, 50, 75, 90, 95]),
        'jaccard_coverage': np.percentile(df['jaccard_percent'], [25, 50, 75, 90, 95])
    }
    
    # Calculate percentiles for confidence
    confidence_percentiles = np.percentile(df['confidence'], [25, 50, 75, 90, 95])
    
    # Print results
    print("\nCoverage Statistics:")
    for metric, stats in coverage_stats.items():
        print(f"\n{metric}:")
        print(stats)
    
    print("\nConfidence Statistics:")
    print(confidence_stats)
    
    print("\nCorrelations:")
    for metric, value in correlations.items():
        print(f"{metric}: {value:.3f}")
    
    print("\nCoverage Percentiles:")
    for metric, percentiles in coverage_percentiles.items():
        print(f"\n{metric}:")
        print(f"25th percentile: {percentiles[0]:.2f}%")
        print(f"Median: {percentiles[1]:.2f}%")
        print(f"75th percentile: {percentiles[2]:.2f}%")
        print(f"90th percentile: {percentiles[3]:.2f}%")
        print(f"95th percentile: {percentiles[4]:.2f}%")
    
    print("\nConfidence Percentiles:")
    print(f"25th percentile: {confidence_percentiles[0]:.3f}")
    print(f"Median: {confidence_percentiles[1]:.3f}")
    print(f"75th percentile: {confidence_percentiles[2]:.3f}")
    print(f"90th percentile: {confidence_percentiles[3]:.3f}")
    print(f"95th percentile: {confidence_percentiles[4]:.3f}")
    
    # Additional analysis: High confidence vs low confidence predictions
    high_conf_threshold = np.percentile(df['confidence'], 75)  # 75th percentile
    high_coverage_threshold = np.percentile(df['jaccard_percent'], 75)  # 75th percentile
    
    high_conf_high_cov = df[(df['confidence'] >= high_conf_threshold) & 
                           (df['jaccard_percent'] >= high_coverage_threshold)]
    
    print(f"\nHigh confidence (≥{high_conf_threshold:.3f}) and high coverage (≥{high_coverage_threshold:.2f}%) predictions:")
    print(f"Count: {len(high_conf_high_cov)}")
    print(f"Percentage of total: {(len(high_conf_high_cov) / len(df) * 100):.2f}%")
    
    return {
        'coverage_stats': coverage_stats,
        'confidence_stats': confidence_stats,
        'correlations': correlations,
        'coverage_percentiles': coverage_percentiles,
        'confidence_percentiles': confidence_percentiles,
        'high_confidence_high_coverage': high_conf_high_cov
    }

if __name__ == "__main__":
    predicted_eeis_file = "results_CONTACT_BASED/predicted_human_eeis_fixed_no_conf_threshold_all_relationships_iden_0.tsv"
    results = analyze_coverage_confidence_distributions(predicted_eeis_file) 