import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), 'utils'))
from utils.coord_gene_mapping import coordinate_to_gene_mapping_with_files

def check_mapping_coverage(mapped_eeis, expression_df, mapping_files_dir):
    """
    Check how many mapped EEIs have both exons/genes present in expression data.
    """
    print(f"\n=== MAPPING COVERAGE ANALYSIS ===")
    
    total_eeis = len(mapped_eeis)
    both_exons_mapped = 0
    one_exon_mapped = 0
    neither_exon_mapped = 0
    
    # Sample a subset for analysis (first 50 EEIs to avoid slowing down)
    sample_size = min(50, total_eeis)
    sample_eeis = mapped_eeis.head(sample_size)
    
    for idx, eei_row in sample_eeis.iterrows():
        mouse_exon1 = eei_row['mouse_exon1']
        mouse_exon2 = eei_row['mouse_exon2']
        
        # Check if exons can be mapped to genes in expression data
        exon1_mapped = find_expression_for_coordinate(mouse_exon1, expression_df, mapping_files_dir) is not None
        exon2_mapped = find_expression_for_coordinate(mouse_exon2, expression_df, mapping_files_dir) is not None
        
        if exon1_mapped and exon2_mapped:
            both_exons_mapped += 1
        elif exon1_mapped or exon2_mapped:
            one_exon_mapped += 1
        else:
            neither_exon_mapped += 1
    
    print(f"Sample analysis (first {sample_size} EEIs):")
    print(f"  Both exons mapped: {both_exons_mapped} ({both_exons_mapped/sample_size*100:.1f}%)")
    print(f"  One exon mapped: {one_exon_mapped} ({one_exon_mapped/sample_size*100:.1f}%)")
    print(f"  Neither exon mapped: {neither_exon_mapped} ({neither_exon_mapped/sample_size*100:.1f}%)")
    
    # Estimate for full dataset
    estimated_both = int(both_exons_mapped * total_eeis / sample_size)
    print(f"\nEstimated for all {total_eeis} EEIs:")
    print(f"  Both exons mapped: ~{estimated_both} ({estimated_both/total_eeis*100:.1f}%)")
    
    return both_exons_mapped, one_exon_mapped, neither_exon_mapped

def test_orthologous_eeis_survival(mapped_eeis_file, expression_file, survival_file, 
                                   cpm_threshold=0.5, p_value_threshold=0.05,
                                   min_group_size=10, output_dir="outputs/", mapping_files_dir=None):
    """
    Test mapped mouse EEIs for survival associations using the methodology 
    from the original cancer paper.
    
    Parameters:
    -----------
    mapped_eeis_file : str
        Path to mapped mouse EEIs file (output from ortholog mapping)
    expression_file : str
        Path to mouse expression data (samples x exons/coordinates)
    survival_file : str
        Path to mouse survival data with columns: [sample_id, survival_time, event_status]
    cpm_threshold : float
        CPM threshold to define exon expression (default: 0.5)
    p_value_threshold : float
        P-value threshold for significance (default: 0.05)
    min_group_size : int
        Minimum number of samples per group for testing (default: 10)
    output_dir : str
        Directory to save results
    mapping_files_dir : str, optional
        Directory containing mapping files for coordinate-to-gene mapping
        
    Returns:
    --------
    tuple: (significant_eeis, all_results, summary_stats)
    """
    
    print("=== MOUSE EEI SURVIVAL ANALYSIS ===")
    print(f"Loading data...")
    
    # Load mapped EEIs
    mapped_eeis = pd.read_csv(mapped_eeis_file, sep='\t')
    print(f"Loaded {len(mapped_eeis)} mapped mouse EEIs")
    
    # Load expression data
    expression_df = load_expression_data(expression_file)
    print(f"Loaded expression data: {expression_df.shape}")
    
    # Load survival data
    survival_df = load_survival_data(survival_file)
    print(f"Loaded survival data: {len(survival_df)} samples")
    
    # Align samples between expression and survival data
    expression_df, survival_df = align_samples(expression_df, survival_df)
    print(f"Aligned data: {len(survival_df)} samples with both expression and survival data")
    
    # Check mapping coverage
    check_mapping_coverage(mapped_eeis, expression_df, mapping_files_dir)
    
    # Test each EEI for survival associations
    print(f"\nTesting {len(mapped_eeis)} EEIs for survival associations...")
    
    results = []
    significant_count = 0
    
    # Debug counters
    missing_expression_count = 0
    group_size_too_small_count = 0
    successful_mapping_count = 0
    
    for idx, eei_row in mapped_eeis.iterrows():
        mouse_exon1 = eei_row['mouse_exon1']
        mouse_exon2 = eei_row['mouse_exon2']
        
        # Test this EEI
        result = test_single_eei_survival(
            mouse_exon1, mouse_exon2, expression_df, survival_df,
            cpm_threshold, min_group_size, eei_row, mapping_files_dir=mapping_files_dir
        )
        
        if result is not None:
            results.append(result)
            successful_mapping_count += 1
            if result['p_value'] <= p_value_threshold:
                significant_count += 1
        else:
            # Track why this EEI failed
            # Check if it's due to missing expression data
            exon1_expr = find_expression_for_coordinate(mouse_exon1, expression_df, mapping_files_dir)
            exon2_expr = find_expression_for_coordinate(mouse_exon2, expression_df, mapping_files_dir)
            
            if exon1_expr is None or exon2_expr is None:
                missing_expression_count += 1
            else:
                # Check if it's due to group size
                eei_present = (exon1_expr > cpm_threshold) & (exon2_expr > cpm_threshold)
                present_samples = eei_present[eei_present == True].index.tolist()
                absent_samples = eei_present[eei_present == False].index.tolist()
                
                if len(present_samples) < min_group_size or len(absent_samples) < min_group_size:
                    group_size_too_small_count += 1
        
        # Progress update
        if (idx + 1) % 50 == 0:
            print(f"Processed {idx + 1}/{len(mapped_eeis)} EEIs...")
    
    # Print debug summary
    print(f"\n=== DEBUG SUMMARY ===")
    print(f"Total mapped EEIs: {len(mapped_eeis)}")
    print(f"EEIs with missing expression data: {missing_expression_count}")
    print(f"EEIs with insufficient group sizes: {group_size_too_small_count}")
    print(f"Successfully mapped and testable EEIs: {successful_mapping_count}")
    print(f"Success rate: {successful_mapping_count}/{len(mapped_eeis)} ({successful_mapping_count/len(mapped_eeis)*100:.1f}%)")
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    if len(results_df) == 0:
        print("No testable EEIs found!")
        return None, None, None
    
    # Filter significant results
    significant_eeis = results_df[results_df['p_value'] <= p_value_threshold].copy()
    significant_eeis = significant_eeis.sort_values('p_value')
    
    # Calculate summary statistics
    summary_stats = calculate_summary_statistics(results_df, significant_eeis, p_value_threshold)
    
    # Print results
    print_results_summary(summary_stats, results_df, significant_eeis)
    
    # Save results
    save_results(results_df, significant_eeis, summary_stats, output_dir)
    
    # Generate plots
    generate_survival_plots(significant_eeis[:5], expression_df, survival_df, cpm_threshold, output_dir)
    
    return significant_eeis, results_df, summary_stats

def load_expression_data(expression_file):
    """Load and process expression data."""
    print(f"Loading expression data from {expression_file}...")
    
    # Try to determine file format
    if expression_file.endswith('.csv'):
        df = pd.read_csv(expression_file, index_col=0)
    else:
        df = pd.read_csv(expression_file, sep='\t', index_col=0)
    
    print(f"Original expression data shape: {df.shape}")
    print(f"First few column names: {list(df.columns[:5])}")
    print(f"First few row names: {list(df.index[:5])}")
    
    # Check if data needs to be transposed
    # If the first column contains gene names and column headers also look like gene names,
    # then the data is transposed and needs to be fixed
    if len(df.columns) > 0 and len(df.index) > 0:
        first_col_name = str(df.columns[0])
        first_row_name = str(df.index[0])
        
        # If both column headers and row names look like gene names (contain common gene patterns)
        gene_patterns = ['LOC', 'ENS', 'Gm', 'Rik', 'Xkr', 'Rp', 'Sox', 'Mrp', 'Lypla', 'Tcea']
        col_looks_like_gene = any(pattern in first_col_name for pattern in gene_patterns)
        row_looks_like_gene = any(pattern in first_row_name for pattern in gene_patterns)
        
        if col_looks_like_gene and row_looks_like_gene:
            print("Detected transposed expression data. Transposing...")
            df = df.T  # Transpose the data
            print(f"After transposition: {df.shape}")
    
    # Ensure sample IDs are strings for matching
    df.columns = df.columns.astype(str)
    
    print(f"Final expression data shape: {df.shape}")
    print(f"Sample IDs (first 5): {list(df.columns[:5])}")
    
    return df

def load_survival_data(survival_file):
    """Load survival data with standard column names."""
    print(f"Loading survival data from {survival_file}...")
    
    if survival_file.endswith('.csv'):
        df = pd.read_csv(survival_file)
    else:
        df = pd.read_csv(survival_file, sep='\t')
    
    print(f"Survival data columns: {list(df.columns)}")
    print(f"Survival data shape: {df.shape}")
    
    # Handle different column name formats
    if 'mouse_id' in df.columns:
        df = df.rename(columns={'mouse_id': 'sample_id'})
    elif 'sample_id' not in df.columns and len(df.columns) >= 1:
        # Assume first column is sample ID
        df = df.rename(columns={df.columns[0]: 'sample_id'})
    
    if 'survival_days' in df.columns:
        df = df.rename(columns={'survival_days': 'survival_time'})
    elif 'survival_time' not in df.columns and len(df.columns) >= 3:
        # Assume third column is survival time
        df = df.rename(columns={df.columns[2]: 'survival_time'})
    
    if 'event_status' not in df.columns and len(df.columns) >= 4:
        # Assume fourth column is event status
        df = df.rename(columns={df.columns[3]: 'event_status'})
    
    # Ensure we have the required columns
    required_cols = ['sample_id', 'survival_time', 'event_status']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"Warning: Missing required columns: {missing_cols}")
        print(f"Available columns: {list(df.columns)}")
        # Try to use available columns
        if len(df.columns) >= 3:
            df = df.iloc[:, :3]  # Take first 3 columns
            df.columns = required_cols
        else:
            raise ValueError(f"Survival data must have at least 3 columns. Found: {list(df.columns)}")
    
    # Ensure proper data types
    df['survival_time'] = pd.to_numeric(df['survival_time'], errors='coerce')
    df['event_status'] = pd.to_numeric(df['event_status'], errors='coerce')
    
    # Remove rows with missing survival data
    df = df.dropna(subset=['survival_time', 'event_status'])
    
    print(f"Processed survival data: {df.shape}")
    print(f"Final columns: {list(df.columns)}")
    
    return df

def align_samples(expression_df, survival_df):
    """Align samples between expression and survival data."""
    expr_samples = set(expression_df.columns)
    surv_samples = set(survival_df['sample_id'])
    
    common_samples = expr_samples.intersection(surv_samples)
    print(f"Found {len(common_samples)} common samples")
    
    # Filter to common samples
    expression_aligned = expression_df[list(common_samples)]
    survival_aligned = survival_df[survival_df['sample_id'].isin(common_samples)].copy()
    
    # Sort both by sample ID for consistency
    expression_aligned = expression_aligned.reindex(sorted(expression_aligned.columns), axis=1)
    survival_aligned = survival_aligned.sort_values('sample_id').reset_index(drop=True)
    
    return expression_aligned, survival_aligned

def test_single_eei_survival(exon1_coord, exon2_coord, expression_df, survival_df, 
                             cpm_threshold, min_group_size, eei_metadata, mapping_files_dir=None):
    """
    Test a single EEI for survival associations using the methodology from the cancer paper.
    
    An EEI is present when both exons are expressed (>threshold).
    An EEI is absent when one or both exons are not expressed (<=threshold).
    """
    
    # Try to find expression data for these coordinates using mapping files
    exon1_expr = find_expression_for_coordinate(exon1_coord, expression_df, mapping_files_dir)
    exon2_expr = find_expression_for_coordinate(exon2_coord, expression_df, mapping_files_dir)
    
    if exon1_expr is None or exon2_expr is None:
        return None
    
    # Define EEI presence: both exons expressed above threshold
    eei_present = (exon1_expr > cpm_threshold) & (exon2_expr > cpm_threshold)
    
    # Create groups for survival analysis
    present_samples = eei_present[eei_present == True].index.tolist()
    absent_samples = eei_present[eei_present == False].index.tolist()
    
    # Check minimum group sizes
    if len(present_samples) < min_group_size or len(absent_samples) < min_group_size:
        return None
    
    # Get survival data for each group
    present_survival = survival_df[survival_df['sample_id'].isin(present_samples)]
    absent_survival = survival_df[survival_df['sample_id'].isin(absent_samples)]
    
    if len(present_survival) == 0 or len(absent_survival) == 0:
        return None
    
    # Perform log-rank test
    try:
        logrank_result = logrank_test(
            present_survival['survival_time'], absent_survival['survival_time'],
            present_survival['event_status'], absent_survival['event_status']
        )
        p_value = logrank_result.p_value
        
        # Calculate median survival times
        kmf = KaplanMeierFitter()
        
        kmf.fit(present_survival['survival_time'], present_survival['event_status'])
        median_present = kmf.median_survival_time_
        
        kmf.fit(absent_survival['survival_time'], absent_survival['event_status'])
        median_absent = kmf.median_survival_time_
        
        # Determine if EEI presence is favorable or unfavorable
        if median_present > median_absent:
            prognostic_type = "Favorable"  # EEI presence associated with longer survival
        else:
            prognostic_type = "Unfavorable"  # EEI presence associated with shorter survival
        
        # Compile results
        result = {
            'mouse_exon1': exon1_coord,
            'mouse_exon2': exon2_coord,
            'human_exon1': eei_metadata.get('human_exon1', ''),
            'human_exon2': eei_metadata.get('human_exon2', ''),
            'p_value': p_value,
            'median_survival_present': median_present,
            'median_survival_absent': median_absent,
            'n_present': len(present_samples),
            'n_absent': len(absent_samples),
            'prognostic_type': prognostic_type,
            'avg_identity': eei_metadata.get('avg_identity', 0),
            'logrank_statistic': logrank_result.test_statistic
        }
        
        return result
        
    except Exception as e:
        print(f"Error testing EEI {exon1_coord}-{exon2_coord}: {e}")
        return None

def find_expression_for_coordinate(coordinate, expression_df, mapping_files_dir=None):
    """
    Find expression data for a genomic coordinate using mapping files if provided.
    """
    if mapping_files_dir is not None:
        mapped_gene = coordinate_to_gene_mapping_with_files(coordinate, mapping_files_dir)
        if mapped_gene:
            # For transposed files, genes are in the index (rows), not columns
            if mapped_gene in expression_df.index:
                return expression_df.loc[mapped_gene]
            elif mapped_gene in expression_df.columns:
                # Fallback for non-transposed files
                return expression_df[mapped_gene]
        return None
    # Fallback to old logic if mapping_files_dir is not provided
    if coordinate in expression_df.columns:
        return expression_df[coordinate]
    return None

def calculate_summary_statistics(results_df, significant_eeis, p_threshold):
    """Calculate summary statistics for the analysis."""
    
    total_tested = len(results_df)
    total_significant = len(significant_eeis)
    
    stats = {
        'total_eeis_tested': total_tested,
        'significant_eeis': total_significant,
        'significance_rate': total_significant / total_tested if total_tested > 0 else 0,
        'median_p_value': results_df['p_value'].median(),
        'min_p_value': results_df['p_value'].min(),
        'favorable_eeis': len(significant_eeis[significant_eeis['prognostic_type'] == 'Favorable']),
        'unfavorable_eeis': len(significant_eeis[significant_eeis['prognostic_type'] == 'Unfavorable']),
    }
    
    if total_significant > 0:
        stats['avg_identity_significant'] = significant_eeis['avg_identity'].mean()
        stats['median_logrank_stat'] = significant_eeis['logrank_statistic'].median()
    
    return stats

def print_results_summary(summary_stats, results_df, significant_eeis):
    """Print a summary of the analysis results."""
    
    print("\n" + "="*50)
    print("SURVIVAL ANALYSIS RESULTS SUMMARY")
    print("="*50)
    
    print(f"Total EEIs tested: {summary_stats['total_eeis_tested']}")
    print(f"Significant EEIs (p≤0.05): {summary_stats['significant_eeis']}")
    print(f"Significance rate: {summary_stats['significance_rate']:.1%}")
    print(f"Median p-value: {summary_stats['median_p_value']:.4f}")
    print(f"Minimum p-value: {summary_stats['min_p_value']:.2e}")
    
    if summary_stats['significant_eeis'] > 0:
        print(f"\nPrognostic Types:")
        print(f"  Favorable EEIs: {summary_stats['favorable_eeis']}")
        print(f"  Unfavorable EEIs: {summary_stats['unfavorable_eeis']}")
        print(f"  Average identity of significant EEIs: {summary_stats.get('avg_identity_significant', 0):.3f}")
        
        print(f"\nTop 5 most significant EEIs:")
        for idx, row in significant_eeis.head().iterrows():
            print(f"  {row['mouse_exon1']} - {row['mouse_exon2']}: p={row['p_value']:.2e} ({row['prognostic_type']})")
    
    print("\n" + "="*50)

def save_results(results_df, significant_eeis, summary_stats, output_dir):
    """Save analysis results to files."""
    
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # Save all results
    results_df.to_csv(f"{output_dir}/all_eei_survival_results.tsv", sep='\t', index=False)
    print(f"Saved all results to {output_dir}/all_eei_survival_results.tsv")
    
    # Save significant results
    if len(significant_eeis) > 0:
        significant_eeis.to_csv(f"{output_dir}/significant_eei_survival_results.tsv", sep='\t', index=False)
        print(f"Saved significant results to {output_dir}/significant_eei_survival_results.tsv")
    
    # Save summary statistics
    import json
    with open(f"{output_dir}/survival_analysis_summary.json", 'w') as f:
        json.dump(summary_stats, f, indent=2)
    print(f"Saved summary statistics to {output_dir}/survival_analysis_summary.json")

def generate_survival_plots(top_eeis, expression_df, survival_df, cpm_threshold, output_dir):
    """Generate Kaplan-Meier survival plots for top significant EEIs."""
    
    import os
    os.makedirs(f"{output_dir}/survival_plots", exist_ok=True)
    
    for idx, eei_row in top_eeis.iterrows():
        try:
            plot_single_eei_survival(eei_row, expression_df, survival_df, cpm_threshold, output_dir)
        except Exception as e:
            print(f"Error plotting EEI {eei_row['mouse_exon1']}-{eei_row['mouse_exon2']}: {e}")

def plot_single_eei_survival(eei_row, expression_df, survival_df, cmp_threshold, output_dir):
    """Plot Kaplan-Meier survival curve for a single EEI."""
    
    exon1_coord = eei_row['mouse_exon1']
    exon2_coord = eei_row['mouse_exon2']
    
    # Get expression data (same logic as in test_single_eei_survival)
    exon1_expr = find_expression_for_coordinate(exon1_coord, expression_df)
    exon2_expr = find_expression_for_coordinate(exon2_coord, expression_df)
    
    if exon1_expr is None or exon2_expr is None:
        return
    
    # Define EEI presence
    eei_present = (exon1_expr > cmp_threshold) & (exon2_expr > cmp_threshold)
    
    # Get survival data for each group
    present_samples = eei_present[eei_present == True].index.tolist()
    absent_samples = eei_present[eei_present == False].index.tolist()
    
    present_survival = survival_df[survival_df['sample_id'].isin(present_samples)]
    absent_survival = survival_df[survival_df['sample_id'].isin(absent_samples)]
    
    # Create Kaplan-Meier plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    kmf = KaplanMeierFitter()
    
    # Plot present group
    kmf.fit(present_survival['survival_time'], present_survival['event_status'], label='EEI Present')
    kmf.plot_survival_function(ax=ax, color='red')
    
    # Plot absent group
    kmf.fit(absent_survival['survival_time'], absent_survival['event_status'], label='EEI Absent')
    kmf.plot_survival_function(ax=ax, color='blue')
    
    # Formatting
    ax.set_title(f"EEI: {exon1_coord} - {exon2_coord}\np-value: {eei_row['p_value']:.2e}")
    ax.set_xlabel('Survival Time (days)')
    ax.set_ylabel('Survival Probability')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Save plot
    filename = f"eei_survival_{exon1_coord.replace(':', '_')}_{exon2_coord.replace(':', '_')}.png"
    filepath = f"{output_dir}/survival_plots/{filename}"
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()

# Main testing function adapted for your specific data
def test_orthologous_eeis_only(mapped_eeis_file, mouse_expression_df, survival_df, 
                               cpm_threshold=0.5, output_dir="outputs/"):
    """
    Test only mouse EEIs that are orthologs of human CRPEs.
    
    This is the main function that implements your requested functionality.
    """
    print("=== Testing Orthologous Mouse EEIs for Survival Associations ===")
    
    # Load mapped EEIs
    mapped_eeis = pd.read_csv(mapped_eeis_file, sep='\t')
    print(f"Loaded {len(mapped_eeis)} mapped mouse EEIs")
    
    # Prepare survival data with proper column names
    survival_data = survival_df.copy()
    if 'mouse_id' in survival_data.columns:
        survival_data = survival_data.rename(columns={'mouse_id': 'sample_id'})
    
    # Align sample IDs between expression and survival data
    expr_samples = set(mouse_expression_df.index)  # Mouse IDs are in index
    surv_samples = set(survival_data['sample_id'].astype(str))
    
    # Convert mouse IDs to string for matching
    mouse_expression_df.index = mouse_expression_df.index.astype(str)
    
    common_samples = expr_samples.intersection(surv_samples)
    print(f"Found {len(common_samples)} samples with both expression and survival data")
    
    if len(common_samples) < 20:
        print("Warning: Very few samples for analysis. Results may not be reliable.")
    
    # Filter to common samples
    expression_aligned = mouse_expression_df.loc[list(common_samples)]
    survival_aligned = survival_data[survival_data['sample_id'].astype(str).isin(common_samples)].copy()
    
    print(f"Proceeding with {len(common_samples)} samples")
    
    # Test each mapped EEI
    results = []
    testable_count = 0
    
    for idx, eei_row in mapped_eeis.iterrows():
        # Extract coordinates
        mouse_exon1 = eei_row['mouse_exon1']
        mouse_exon2 = eei_row['mouse_exon2']
        
        # For this implementation, we'll need to map coordinates to genes
        # This is a placeholder - you'll need to implement coordinate->gene mapping
        gene1 = coordinate_to_gene_mapping_with_files(mouse_exon1, "utils/coord_gene_mapping/")
        gene2 = coordinate_to_gene_mapping_with_files(mouse_exon2, "utils/coord_gene_mapping/")
        
        if gene1 is None or gene2 is None:
            continue
            
        testable_count += 1
        
        # Get expression data for these genes
        expr1 = expression_aligned[gene1]
        expr2 = expression_aligned[gene2]
        
        # Define EEI presence: both exons expressed above threshold
        eei_present = (expr1 > cpm_threshold) & (expr2 > cpm_threshold)
        
        # Create groups
        present_samples = eei_present[eei_present == True].index.tolist()
        absent_samples = eei_present[eei_present == False].index.tolist()
        
        # Check minimum group sizes (at least 5 per group for small dataset)
        min_size = max(5, len(common_samples) // 10)
        if len(present_samples) < min_size or len(absent_samples) < min_size:
            continue
        
        # Get survival data
        present_survival = survival_aligned[survival_aligned['sample_id'].astype(str).isin(present_samples)]
        absent_survival = survival_aligned[survival_aligned['sample_id'].astype(str).isin(absent_samples)]
        
        # Perform log-rank test
        try:
            logrank_result = logrank_test(
                present_survival['survival_days'], absent_survival['survival_days'],
                present_survival['event_status'], absent_survival['event_status']
            )
            
            p_value = logrank_result.p_value
            
            # Calculate median survival
            kmf = KaplanMeierFitter()
            kmf.fit(present_survival['survival_days'], present_survival['event_status'])
            median_present = kmf.median_survival_time_
            
            kmf.fit(absent_survival['survival_days'], absent_survival['event_status'])
            median_absent = kmf.median_survival_time_
            
            # Determine prognostic type
            prognostic_type = "Favorable" if median_present > median_absent else "Unfavorable"
            
            # Store results
            result = {
                'mouse_exon1': mouse_exon1,
                'mouse_exon2': mouse_exon2,
                'mapped_gene1': gene1,
                'mapped_gene2': gene2,
                'human_exon1': eei_row.get('human_exon1', ''),
                'human_exon2': eei_row.get('human_exon2', ''),
                'p_value': p_value,
                'median_survival_present': median_present,
                'median_survival_absent': median_absent,
                'n_present': len(present_samples),
                'n_absent': len(absent_samples),
                'prognostic_type': prognostic_type,
                'avg_identity': eei_row.get('avg_identity', 0),
                'logrank_statistic': logrank_result.test_statistic
            }
            
            results.append(result)
            
        except Exception as e:
            print(f"Error testing EEI {mouse_exon1}-{mouse_exon2}: {e}")
            continue
        
        # Progress update
        if (idx + 1) % 50 == 0:
            print(f"Processed {idx + 1}/{len(mapped_eeis)} EEIs...")
    
    # Process results
    if len(results) == 0:
        print("No testable EEIs found! Check coordinate->gene mapping.")
        return None, None
    
    results_df = pd.DataFrame(results)
    significant_eeis = results_df[results_df['p_value'] <= 0.05].sort_values('p_value')
    
    # Print summary
    print(f"\nSUMMARY:")
    print(f"Total mapped EEIs: {len(mapped_eeis)}")
    print(f"Testable EEIs: {len(results_df)}")
    print(f"Significant EEIs (p≤0.05): {len(significant_eeis)}")
    print(f"Significance rate: {len(significant_eeis)/len(results_df):.1%}")
    
    if len(significant_eeis) > 0:
        print(f"\nTop 5 significant EEIs:")
        for idx, row in significant_eeis.head().iterrows():
            print(f"  {row['mapped_gene1']}-{row['mapped_gene2']}: p={row['p_value']:.3f} ({row['prognostic_type']})")
    
    # Save results
    os.makedirs(output_dir, exist_ok=True)
    results_df.to_csv(f"{output_dir}/orthologous_eei_survival_results.tsv", sep='\t', index=False)
    if len(significant_eeis) > 0:
        significant_eeis.to_csv(f"{output_dir}/significant_orthologous_eeis.tsv", sep='\t', index=False)
    
    return significant_eeis, results_df

def main():
    """
    Main function to run the EEI survival analysis pipeline.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Run EEI survival analysis with mapping files')
    parser.add_argument('--mapped-eeis', required=True, 
                       help='Path to mapped mouse EEIs file (output from ortholog mapping)')
    parser.add_argument('--expression', required=True,
                       help='Path to mouse expression data file')
    parser.add_argument('--survival', required=True,
                       help='Path to mouse survival data file')
    parser.add_argument('--mapping-files-dir', required=True,
                       help='Directory containing mapping files (created by coord_gene_mapping.py)')
    parser.add_argument('--output-dir', default='outputs',
                       help='Output directory for results (default: outputs)')
    parser.add_argument('--cpm-threshold', type=float, default=0.5,
                       help='CPM threshold for exon expression (default: 0.5)')
    parser.add_argument('--p-threshold', type=float, default=0.05,
                       help='P-value threshold for significance (default: 0.05)')
    parser.add_argument('--min-group-size', type=int, default=10,
                       help='Minimum samples per group for testing (default: 10)')
    
    args = parser.parse_args()
    
    # Check if mapping files exist
    combined_mapping_file = os.path.join(args.mapping_files_dir, "combined_mapping.json")
    if not os.path.exists(combined_mapping_file):
        print(f"Error: Mapping files not found in {args.mapping_files_dir}")
        print("Please run the mapping creation script first:")
        print("python utils/coord_gene_mapping.py --gtf <gtf_file> --expression <expression_file> --output-dir <mapping_dir>")
        return
    
    print("=== EEI SURVIVAL ANALYSIS PIPELINE ===")
    print(f"Mapped EEIs file: {args.mapped_eeis}")
    print(f"Expression file: {args.expression}")
    print(f"Survival file: {args.survival}")
    print(f"Mapping files directory: {args.mapping_files_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"CPM threshold: {args.cpm_threshold}")
    print(f"P-value threshold: {args.p_threshold}")
    print(f"Minimum group size: {args.min_group_size}")
    print("=" * 50)
    
    # Run the analysis
    try:
        significant_eeis, all_results, summary_stats = test_orthologous_eeis_survival(
            mapped_eeis_file=args.mapped_eeis,
            expression_file=args.expression,
            survival_file=args.survival,
            cpm_threshold=args.cpm_threshold,
            p_value_threshold=args.p_threshold,
            min_group_size=args.min_group_size,
            output_dir=args.output_dir,
            mapping_files_dir=args.mapping_files_dir
        )
        
        if significant_eeis is not None:
            print(f"\nAnalysis completed successfully!")
            print(f"Results saved to: {args.output_dir}")
        else:
            print("\nAnalysis completed but no testable EEIs were found.")
            
    except Exception as e:
        print(f"Error running analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()