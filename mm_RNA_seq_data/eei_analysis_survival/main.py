import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from data_loader import load_expression_data, load_metadata, load_exon_mappings
from generate_plots import generate_validation_plots
from analysis import create_response_based_survival_proxy, analyze_eei_ifn_correlation, identify_time_dependent_eeis, create_pseudo_survival_analysis




def main_analysis_pipeline(mapped_eeis_file, expression_file, metadata_file, 
                          mapping_files_dir, output_dir):
    """
    Main pipeline using treatment response as survival proxy.
    """
    
    print("=== EEI VALIDATION USING TREATMENT RESPONSE ===")
    
    # Load data
    mapped_eeis = pd.read_csv(mapped_eeis_file, sep='\t')
    expression_df = load_expression_data(expression_file)
    metadata_df = load_metadata(metadata_file)
    coord_to_exon = load_exon_mappings(mapping_files_dir)
    
    # Create response-based survival proxy
    metadata_df = create_response_based_survival_proxy(
        metadata_df, expression_df, mapped_eeis, coord_to_exon
    )
    
    print(f"\nResponse group distribution:")
    print(metadata_df['response_group'].value_counts())
    
    # Analyze EEI-response associations
    print("\nAnalyzing EEI associations with treatment response...")
    response_results = analyze_eei_ifn_correlation(
        mapped_eeis, expression_df, metadata_df, coord_to_exon
    )
    
    # Identify time-dependent EEIs
    print("\nIdentifying time-dependent EEIs...")
    time_dependent = identify_time_dependent_eeis(
        mapped_eeis, expression_df, metadata_df, coord_to_exon
    )
    
    # Pseudo-survival analysis
    print("\nPerforming pseudo-survival analysis...")
    survival_results = create_pseudo_survival_analysis(response_results, metadata_df)
    
    # Generate validation plots
    print("\nGenerating validation plots...")
    generate_validation_plots(response_results, output_dir)
    
    # Save results
    if len(response_results) > 0:
        response_results.to_csv(f"{output_dir}/eei_response_associations.tsv", sep='\t', index=False)
        print(f"Saved response associations to {output_dir}/eei_response_associations.tsv")
    
    if len(time_dependent) > 0:
        time_dependent.to_csv(f"{output_dir}/time_dependent_eeis.tsv", sep='\t', index=False)
        print(f"Saved time-dependent EEIs to {output_dir}/time_dependent_eeis.tsv")
    
    # Summary statistics
    print("\n=== ANALYSIS SUMMARY ===")
    print(f"Total EEIs analyzed: {len(mapped_eeis)}")
    print(f"EEIs with expression data: {len(response_results)}")
    # Prefer adjusted p-values if available
    if 'p_adj_fdr_bh' in response_results.columns:
        sig_count = (response_results['p_adj_fdr_bh'] < 0.05).sum()
        print(f"Significant EEIs (FDR<0.05): {sig_count}")
    elif 'p_value' in response_results.columns:
        sig_count = (response_results['p_value'] < 0.05).sum()
        print(f"Significant EEIs (p<0.05): {sig_count}")
    print(f"Time-dependent EEIs: {len(time_dependent)}")
    
    if survival_results:
        print(f"Top response-associated EEIs: {survival_results['significant_count']}")
    
    return response_results, time_dependent, survival_results


if __name__ == "__main__":
    # Run the analysis
    results = main_analysis_pipeline(
        mapped_eeis_file="../outputs/mapped_mouse_crpes_brca.tsv",
        expression_file="../crc2022_exon_expression/crc2022_exon_cpm.tsv",
        metadata_file="../crc2022_exon_expression/samplesheet.csv",
        mapping_files_dir="../exon_mapping_files",
        output_dir="../results_survival_proxy_mouse_crpes_brca_corrected_pvalues"
    )