import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')
import os
import sys

def load_exon_mappings(mapping_files_dir):
    """
    Load exon coordinate mappings once and cache them.
    """
    mapping_file = os.path.join(mapping_files_dir, 'exon_coord_mapping_mouse.tsv')
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    # Create dictionaries for fast lookup
    coord_to_exon = {}
    
    for _, row in mapping_df.iterrows():
        coord = row['coord']
        exon_id = row['exon_id']

        # Skip missing values
        if pd.isna(coord) or pd.isna(exon_id):
            continue

        # Normalize to string for consistent lookups
        coord = str(coord)
        exon_id = str(exon_id)
        
        # Also create version without strand for flexible matching
        parts = coord.split(':')
        if len(parts) == 4:
            coord_no_strand = f"{parts[0]}:{parts[1]}:{parts[2]}"
            coord_to_exon[coord_no_strand] = exon_id
            
            # Also store with dash format
            coord_dash = f"{parts[0]}:{parts[1]}-{parts[2]}"
            coord_to_exon[coord_dash] = exon_id
        
        coord_to_exon[coord] = exon_id
    
    return coord_to_exon

def find_expression_with_mapping(coordinate, expression_df, coord_to_exon):
    """
    Find expression data using pre-loaded mappings.
    """
    # Guard against non-string or missing coordinates
    if coordinate is None or (isinstance(coordinate, float) and np.isnan(coordinate)):
        return None
    coordinate = str(coordinate).strip()
    if coordinate == '' or coordinate.lower() == 'nan':
        return None

    # Direct ENSMUSE ID lookup
    if coordinate.startswith('ENSMUSE'):
        if coordinate in expression_df.index:
            return expression_df.loc[coordinate]
        return None
    
    # Use mapping dictionary
    if coord_to_exon and coordinate in coord_to_exon:
        exon_id = coord_to_exon[coordinate]
        if exon_id in expression_df.index:
            return expression_df.loc[exon_id]
    
    # Try alternative coordinate formats
    if '-' in coordinate and coord_to_exon:
        # Convert chr1:12345-67890 to chr1:12345:67890
        parts = coordinate.replace('-', ':').split(':')
        if len(parts) == 3:
            # Try both strands
            for strand in ['1', '-1']:
                alt_coord = f"{parts[0]}:{parts[1]}:{parts[2]}:{strand}"
                if alt_coord in coord_to_exon:
                    exon_id = coord_to_exon[alt_coord]
                    if exon_id in expression_df.index:
                        return expression_df.loc[exon_id]
    
    return None

def load_expression_data(expression_file):
    """Load and process expression data."""
    print(f"Loading expression data from {expression_file}...")
    
    if expression_file.endswith('.csv'):
        df = pd.read_csv(expression_file, index_col=0)
    else:
        df = pd.read_csv(expression_file, sep='\t', index_col=0)
    
    print(f"Expression data shape: {df.shape}")
    print(f"First few exon IDs: {list(df.index[:5])}")
    print(f"Sample IDs (first 5): {list(df.columns[:5])}")
    
    # Ensure sample IDs are strings for matching
    df.columns = df.columns.astype(str)
    
    return df

def load_metadata(metadata_file):
    """Load sample metadata with treatment information."""
    print(f"Loading metadata from {metadata_file}...")
    
    if metadata_file.endswith('.csv'):
        df = pd.read_csv(metadata_file)
    else:
        df = pd.read_csv(metadata_file, sep='\t')
    
    print(f"Metadata shape: {df.shape}")
    print(f"Metadata columns: {list(df.columns)}")
    
    # Use 'Run' column as sample_id to match expression data
    if 'Run' in df.columns:
        df['sample_id'] = df['Run'].astype(str)
    else:
        print("Warning: 'Run' column not found in metadata")
    
    # Create simplified binary treatment column
    if 'treatment' in df.columns:
        df['treatment_binary'] = df['treatment'].apply(
            lambda x: 'Control' if x == 'Vehicle' else 'Treated'
        )
        # Optional: basic counts for sanity check
        try:
            ctrl_n = (df['treatment_binary'] == 'Control').sum()
            trt_n = (df['treatment_binary'] == 'Treated').sum()
            print(f"treatment_binary counts -> Control: {ctrl_n}, Treated: {trt_n}")
        except Exception:
            pass
    
    return df

def test_eei_treatment_association(mapped_eeis_file, expression_file, metadata_file,
                                  cpm_threshold=0.5, p_value_threshold=0.05,
                                  min_group_size=3, output_dir="outputs/", 
                                  mapping_files_dir=None, test_variable='treatment'):
    """
    Test associations between EEI presence and treatment/phenotype variables.
    
    Parameters:
    -----------
    mapped_eeis_file : str
        Path to mapped mouse EEIs file
    expression_file : str
        Path to exon expression data (CPM)
    metadata_file : str
        Path to sample metadata with treatment info
    cpm_threshold : float
        CPM threshold to define exon expression
    p_value_threshold : float
        P-value threshold for significance
    min_group_size : int
        Minimum samples per group for testing
    output_dir : str
        Directory to save results
    mapping_files_dir : str
        Directory containing exon mapping files
    test_variable : str
        Which metadata variable to test ('treatment', 'model', 'timepoint')
    """
    
    print("=== EEI TREATMENT ASSOCIATION ANALYSIS ===")
    print(f"Testing associations with: {test_variable}")
    
    # Load mappings
    coord_to_exon = None
    if mapping_files_dir:
        print(f"Loading exon coordinate mappings from {mapping_files_dir}...")
        coord_to_exon = load_exon_mappings(mapping_files_dir)
        print(f"Loaded {len(coord_to_exon)} coordinate mappings")
    
    # Load data
    mapped_eeis = pd.read_csv(mapped_eeis_file, sep='\t')
    print(f"Loaded {len(mapped_eeis)} mapped mouse EEIs")
    
    expression_df = load_expression_data(expression_file)
    metadata_df = load_metadata(metadata_file)
    
    # Align samples
    expr_samples = set(expression_df.columns)
    meta_samples = set(metadata_df['sample_id'])
    common_samples = list(expr_samples.intersection(meta_samples))
    
    print(f"Found {len(common_samples)} samples with both expression and metadata")
    
    if len(common_samples) < 6:
        print("Error: Too few samples for meaningful analysis")
        return None, None, None
    
    # Filter to common samples
    expression_aligned = expression_df[common_samples]
    metadata_aligned = metadata_df[metadata_df['sample_id'].isin(common_samples)].copy()
    metadata_aligned = metadata_aligned.set_index('sample_id')

    # Quick sanity check: print expression ranges for first few EEIs
    try:
        print("\nExpression range sanity check for first 5 EEIs:")
        for eei in mapped_eeis.head(5).itertuples(index=False):
            mouse_exon1 = getattr(eei, 'mouse_exon1')
            mouse_exon2 = getattr(eei, 'mouse_exon2')
            exon1_expr = find_expression_with_mapping(mouse_exon1, expression_aligned, coord_to_exon)
            exon2_expr = find_expression_with_mapping(mouse_exon2, expression_aligned, coord_to_exon)
            print(f"EEI: {mouse_exon1} - {mouse_exon2}")
            if exon1_expr is not None:
                print(f"  Exon1 range: {exon1_expr.min():.3f} - {exon1_expr.max():.3f}")
            else:
                print("  Exon1: not found in expression data")
            if exon2_expr is not None:
                print(f"  Exon2 range: {exon2_expr.min():.3f} - {exon2_expr.max():.3f}")
            else:
                print("  Exon2: not found in expression data")
    except Exception as e:
        print(f"Warning during expression range sanity check: {e}")
    
    # Diagnostic: Analyze EEI presence patterns across treatment groups
    print(f"\nDiagnostic: EEI presence patterns with CPM threshold {cpm_threshold}")
    print("Analyzing first 20 EEIs for presence patterns...")
    
    presence_patterns = []
    for eei in mapped_eeis.head(20).itertuples(index=False):
        mouse_exon1 = getattr(eei, 'mouse_exon1')
        mouse_exon2 = getattr(eei, 'mouse_exon2')
        
        exon1_expr = find_expression_with_mapping(mouse_exon1, expression_aligned, coord_to_exon)
        exon2_expr = find_expression_with_mapping(mouse_exon2, expression_aligned, coord_to_exon)
        
        if exon1_expr is None or exon2_expr is None:
            continue
            
        # Define EEI presence
        eei_present = (exon1_expr > cpm_threshold) & (exon2_expr > cpm_threshold)
        
        # Calculate presence rates per group
        presence_by_group = eei_present.groupby(metadata_aligned[test_variable]).agg(['sum', 'count'])
        presence_by_group['rate'] = presence_by_group['sum'] / presence_by_group['count']
        
        # Store pattern info
        pattern_info = {
            'eei': f"{mouse_exon1} - {mouse_exon2}",
            'total_present': eei_present.sum(),
            'total_samples': len(eei_present),
            'presence_rate': eei_present.mean(),
            'control_rate': presence_by_group.loc['Control', 'rate'] if 'Control' in presence_by_group.index else 0,
            'treated_rate': presence_by_group.loc['Treated', 'rate'] if 'Treated' in presence_by_group.index else 0
        }
        presence_patterns.append(pattern_info)
    
    # Print summary of presence patterns
    if presence_patterns:
        print(f"\nEEI Presence Pattern Summary (first 20):")
        print(f"{'EEI':<50} {'Total':<8} {'Rate':<8} {'Control':<8} {'Treated':<8}")
        print("-" * 90)
        for pattern in presence_patterns:
            print(f"{pattern['eei'][:49]:<50} {pattern['total_present']:<8} {pattern['presence_rate']:.3f} {pattern['control_rate']:.3f} {pattern['treated_rate']:.3f}")
        
        # Calculate overall statistics
        avg_presence_rate = np.mean([p['presence_rate'] for p in presence_patterns])
        avg_control_rate = np.mean([p['control_rate'] for p in presence_patterns])
        avg_treated_rate = np.mean([p['treated_rate'] for p in presence_patterns])
        
        print(f"\nAverage presence rates:")
        print(f"  Overall: {avg_presence_rate:.3f}")
        print(f"  Control: {avg_control_rate:.3f}")
        print(f"  Treated: {avg_treated_rate:.3f}")
        print(f"  Difference: {abs(avg_control_rate - avg_treated_rate):.3f}")
        
        # Check if most EEIs are present in all/most samples
        high_presence = sum(1 for p in presence_patterns if p['presence_rate'] > 0.8)
        low_presence = sum(1 for p in presence_patterns if p['presence_rate'] < 0.2)
        print(f"\nEEI presence distribution:")
        print(f"  High presence (>80% samples): {high_presence}")
        print(f"  Low presence (<20% samples): {low_presence}")
        print(f"  Medium presence (20-80% samples): {len(presence_patterns) - high_presence - low_presence}")
    
    # Test each EEI
    results = []
    
    for idx, eei_row in mapped_eeis.iterrows():
        mouse_exon1 = eei_row['mouse_exon1']
        mouse_exon2 = eei_row['mouse_exon2']
        
        # Find expression for both exons
        exon1_expr = find_expression_with_mapping(mouse_exon1, expression_aligned, coord_to_exon)
        exon2_expr = find_expression_with_mapping(mouse_exon2, expression_aligned, coord_to_exon)
        
        if exon1_expr is None or exon2_expr is None:
            continue
        
        # Define EEI presence: both exons above threshold
        eei_present = (exon1_expr > cpm_threshold) & (exon2_expr > cpm_threshold)
        
        # Test association with treatment variable
        result = test_single_eei_quantitative(
            eei_row, expression_aligned, metadata_aligned, test_variable,
            min_group_size, coord_to_exon
        )
        
        if result is not None:
            results.append(result)
        
        # Progress update
        if (idx + 1) % 100 == 0:
            print(f"Processed {idx + 1}/{len(mapped_eeis)} EEIs...")
    
    # Process results
    if len(results) == 0:
        print("No testable EEIs found!")
        return None, None, None
    
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    from statsmodels.stats.multitest import multipletests
    results_df['p_adj'] = multipletests(results_df['min_p_value'], method='fdr_bh')[1]
    
    # Filter significant results
    significant_eeis = results_df[results_df['min_p_value'] <= p_value_threshold].copy()
    significant_eeis = significant_eeis.sort_values('min_p_value')
    
    # Calculate summary statistics
    summary_stats = {
        'total_eeis_tested': len(results_df),
        'significant_eeis': len(significant_eeis),
        'significance_rate': len(significant_eeis) / len(results_df) if len(results_df) > 0 else 0,
        'median_p_value': results_df['min_p_value'].median(),
        'min_p_value': results_df['min_p_value'].min(),
        'test_variable': test_variable
    }
    
    # Print summary
    print_results_summary(summary_stats, results_df, significant_eeis)
    
    # Save results
    save_results(results_df, significant_eeis, summary_stats, output_dir, test_variable)
    
    # Generate plots for top hits
    if len(significant_eeis) > 0:
        generate_association_plots(
            significant_eeis.head(5), expression_aligned, metadata_aligned,
            cpm_threshold, test_variable, output_dir, coord_to_exon
        )
    
    return significant_eeis, results_df, summary_stats

def test_single_eei_association(eei_present, metadata_df, test_variable,
                               min_group_size, eei_metadata):
    """
    Test association between EEI expression and a metadata variable using quantitative analysis.
    """
    # Get the test variable values for each sample
    test_values = metadata_df[test_variable]
    
    # Create contingency table
    contingency_df = pd.crosstab(eei_present, test_values)
    
    # Check minimum group sizes
    if contingency_df.min().min() < min_group_size:
        return None
    
    # Perform appropriate statistical test
    if len(contingency_df.columns) == 2:
        # Binary variable - use Fisher's exact test
        from scipy.stats import fisher_exact
        oddsratio, p_value = fisher_exact(contingency_df)
        test_type = 'fisher_exact'
        effect_size = oddsratio
    else:
        # Multiple categories - use chi-square test
        from scipy.stats import chi2_contingency
        chi2, p_value, dof, expected = chi2_contingency(contingency_df)
        test_type = 'chi_square'
        # Calculate Cramér's V as effect size
        n = contingency_df.sum().sum()
        min_dim = min(contingency_df.shape[0] - 1, contingency_df.shape[1] - 1)
        effect_size = np.sqrt(chi2 / (n * min_dim)) if min_dim > 0 else 0
    
    # Calculate EEI presence rate per group
    presence_rates = {}
    for col in contingency_df.columns:
        total = contingency_df[col].sum()
        present = contingency_df.loc[True, col] if True in contingency_df.index else 0
        presence_rates[col] = present / total if total > 0 else 0
    
    result = {
        'mouse_exon1': eei_metadata['mouse_exon1'],
        'mouse_exon2': eei_metadata['mouse_exon2'],
        'human_exon1': eei_metadata.get('human_exon1', ''),
        'human_exon2': eei_metadata.get('human_exon2', ''),
        'p_value': p_value,
        'test_type': test_type,
        'effect_size': effect_size,
        'n_eei_present': eei_present.sum(),
        'n_eei_absent': (~eei_present).sum(),
        'presence_rates': str(presence_rates),
        'avg_identity': eei_metadata.get('avg_identity', 0)
    }
    
    return result

def test_single_eei_quantitative(eei_row, expression_df, metadata_df, test_variable,
                                min_group_size, coord_to_exon):
    """
    Test quantitative association between EEI expression and treatment using multiple approaches.
    """
    mouse_exon1 = eei_row['mouse_exon1']
    mouse_exon2 = eei_row['mouse_exon2']
    
    # Get expression data for both exons
    exon1_expr = find_expression_with_mapping(mouse_exon1, expression_df, coord_to_exon)
    exon2_expr = find_expression_with_mapping(mouse_exon2, expression_df, coord_to_exon)
    
    if exon1_expr is None or exon2_expr is None:
        return None
    
    # Get treatment groups
    control_samples = metadata_df[metadata_df[test_variable] == 'Control'].index
    treated_samples = metadata_df[metadata_df[test_variable] == 'Treated'].index
    
    # Check minimum group sizes
    if len(control_samples) < min_group_size or len(treated_samples) < min_group_size:
        return None
    
    # Extract expression values for each group
    exon1_control = exon1_expr[control_samples]
    exon1_treated = exon1_expr[treated_samples]
    exon2_control = exon2_expr[control_samples]
    exon2_treated = exon2_expr[treated_samples]
    
    results = {}
    
    # Test 1: Individual exon expression differences
    try:
        from scipy.stats import ttest_ind, mannwhitneyu
        
        # T-test for exon1
        t_stat1, p_value1 = ttest_ind(exon1_control, exon1_treated)
        results['exon1_ttest_p'] = p_value1
        results['exon1_ttest_stat'] = t_stat1
        
        # T-test for exon2
        t_stat2, p_value2 = ttest_ind(exon2_control, exon2_treated)
        results['exon2_ttest_p'] = p_value2
        results['exon2_ttest_stat'] = t_stat2
        
        # Mann-Whitney U test (non-parametric alternative)
        u_stat1, mw_p1 = mannwhitneyu(exon1_control, exon1_treated, alternative='two-sided')
        u_stat2, mw_p2 = mannwhitneyu(exon2_control, exon2_treated, alternative='two-sided')
        
        results['exon1_mannwhitney_p'] = mw_p1
        results['exon2_mannwhitney_p'] = mw_p2
        
    except Exception as e:
        print(f"Warning: Error in statistical tests for {mouse_exon1}-{mouse_exon2}: {e}")
        return None
    
    # Test 2: Correlation differences between exons
    try:
        # Calculate correlation in each group
        if len(exon1_control) > 2 and len(exon1_treated) > 2:
            corr_control = np.corrcoef(exon1_control, exon2_control)[0, 1]
            corr_treated = np.corrcoef(exon1_treated, exon2_treated)[0, 1]
            
            # Fisher's Z transformation for correlation comparison
            from scipy.stats import fisher_exact
            z_control = np.arctanh(corr_control) if not np.isnan(corr_control) else 0
            z_treated = np.arctanh(corr_treated) if not np.isnan(corr_treated) else 0
            
            results['corr_control'] = corr_control
            results['corr_treated'] = corr_treated
            results['corr_difference'] = abs(corr_control - corr_treated)
            
            # Simple correlation difference test (not statistically rigorous but informative)
            results['corr_diff_significant'] = results['corr_difference'] > 0.3
        else:
            results['corr_control'] = np.nan
            results['corr_treated'] = np.nan
            results['corr_difference'] = np.nan
            results['corr_diff_significant'] = False
            
    except Exception as e:
        results['corr_control'] = np.nan
        results['corr_treated'] = np.nan
        results['corr_difference'] = np.nan
        results['corr_diff_significant'] = False
    
    # Test 3: Combined EEI expression (sum or product)
    try:
        # Sum of both exons
        combined_control = exon1_control + exon2_control
        combined_treated = exon1_treated + exon2_treated
        
        t_stat_combined, p_combined = ttest_ind(combined_control, combined_treated)
        results['combined_ttest_p'] = p_combined
        results['combined_ttest_stat'] = t_stat_combined
        
        # Product of both exons (captures interaction effects)
        product_control = exon1_control * exon2_control
        product_treated = exon1_treated * exon2_treated
        
        t_stat_product, p_product = ttest_ind(product_control, product_treated)
        results['product_ttest_p'] = p_product
        results['product_ttest_stat'] = t_stat_product
        
    except Exception as e:
        results['combined_ttest_p'] = np.nan
        results['combined_ttest_stat'] = np.nan
        results['product_ttest_p'] = np.nan
        results['product_ttest_stat'] = np.nan
    
    # Calculate effect sizes and summary statistics
    try:
        # Mean expression differences
        results['exon1_mean_diff'] = exon1_treated.mean() - exon1_control.mean()
        results['exon2_mean_diff'] = exon2_treated.mean() - exon2_control.mean()
        results['combined_mean_diff'] = (exon1_treated + exon2_treated).mean() - (exon1_control + exon2_control).mean()
        
        # Fold changes
        results['exon1_fold_change'] = exon1_treated.mean() / exon1_control.mean() if exon1_control.mean() > 0 else np.nan
        results['exon2_fold_change'] = exon2_treated.mean() / exon2_control.mean() if exon2_control.mean() > 0 else np.nan
        
        # Coefficient of variation (measure of expression stability)
        results['exon1_cv_control'] = exon1_control.std() / exon1_control.mean() if exon1_control.mean() > 0 else np.nan
        results['exon1_cv_treated'] = exon1_treated.std() / exon1_treated.mean() if exon1_treated.mean() > 0 else np.nan
        results['exon2_cv_control'] = exon2_control.std() / exon2_control.mean() if exon2_control.mean() > 0 else np.nan
        results['exon2_cv_treated'] = exon2_treated.std() / exon2_treated.mean() if exon2_treated.mean() > 0 else np.nan
        
    except Exception as e:
        # Set defaults if calculations fail
        results['exon1_mean_diff'] = np.nan
        results['exon2_mean_diff'] = np.nan
        results['combined_mean_diff'] = np.nan
        results['exon1_fold_change'] = np.nan
        results['exon2_fold_change'] = np.nan
        results['exon1_cv_control'] = np.nan
        results['exon1_cv_treated'] = np.nan
        results['exon2_cv_control'] = np.nan
        results['exon2_cv_treated'] = np.nan
    
    # Store basic metadata
    results['mouse_exon1'] = mouse_exon1
    results['mouse_exon2'] = mouse_exon2
    results['human_exon1'] = eei_row.get('human_exon1', '')
    results['human_exon2'] = eei_row.get('human_exon2', '')
    results['avg_identity'] = eei_row.get('avg_identity', 0)
    
    # Calculate overall significance (minimum p-value across all tests)
    p_values = [results.get('exon1_ttest_p', 1), results.get('exon2_ttest_p', 1), 
                results.get('combined_ttest_p', 1), results.get('product_ttest_p', 1)]
    p_values = [p for p in p_values if not np.isnan(p)]
    
    if p_values:
        results['min_p_value'] = min(p_values)
        results['best_test'] = ['exon1_ttest', 'exon2_ttest', 'combined_ttest', 'product_ttest'][p_values.index(min(p_values))]
    else:
        results['min_p_value'] = np.nan
        results['best_test'] = 'none'
    
    return results

def print_results_summary(summary_stats, results_df, significant_eeis):
    """Print a summary of the analysis results."""
    
    print("\n" + "="*50)
    print("TREATMENT ASSOCIATION ANALYSIS RESULTS")
    print("="*50)
    
    print(f"Test variable: {summary_stats['test_variable']}")
    print(f"Total EEIs tested: {summary_stats['total_eeis_tested']}")
    print(f"Significant EEIs (p≤0.05): {summary_stats['significant_eeis']}")
    print(f"Significance rate: {summary_stats['significance_rate']:.1%}")
    print(f"Median p-value: {summary_stats['median_p_value']:.4f}")
    print(f"Minimum p-value: {summary_stats['min_p_value']:.2e}")
    
    if len(significant_eeis) > 0:
        print(f"\nTop 5 most significant EEIs:")
        for idx, row in significant_eeis.head().iterrows():
            print(f"  {row['mouse_exon1']} - {row['mouse_exon2']}: p={row['min_p_value']:.2e}")
            print(f"    Best test: {row.get('best_test', 'unknown')}")
            if 'exon1_fold_change' in row and not pd.isna(row['exon1_fold_change']):
                print(f"    Exon1 fold change: {row['exon1_fold_change']:.2f}")
            if 'exon2_fold_change' in row and not pd.isna(row['exon2_fold_change']):
                print(f"    Exon2 fold change: {row['exon2_fold_change']:.2f}")
            if 'corr_difference' in row and not pd.isna(row['corr_difference']):
                print(f"    Correlation difference: {row['corr_difference']:.3f}")
    
    # Print summary of test types that were most successful
    if len(results_df) > 0 and 'best_test' in results_df.columns:
        test_counts = results_df['best_test'].value_counts()
        print(f"\nMost successful test types:")
        for test, count in test_counts.items():
            print(f"  {test}: {count} EEIs")
    
    print("\n" + "="*50)

def save_results(results_df, significant_eeis, summary_stats, output_dir, test_variable):
    """Save analysis results to files."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save all results
    results_df.to_csv(f"{output_dir}/eei_{test_variable}_association_results.tsv", 
                     sep='\t', index=False)
    print(f"Saved all results to {output_dir}/eei_{test_variable}_association_results.tsv")
    
    # Save significant results
    if len(significant_eeis) > 0:
        significant_eeis.to_csv(f"{output_dir}/significant_eei_{test_variable}_associations.tsv", 
                               sep='\t', index=False)
        print(f"Saved significant results to {output_dir}/significant_eei_{test_variable}_associations.tsv")
    
    # Save summary statistics
    import json
    with open(f"{output_dir}/{test_variable}_association_summary.json", 'w') as f:
        json.dump(summary_stats, f, indent=2)
    print(f"Saved summary to {output_dir}/{test_variable}_association_summary.json")
    
    # Print column information for debugging
    print(f"\nResults columns: {list(results_df.columns)}")
    print(f"Number of columns: {len(results_df.columns)}")

def generate_association_plots(top_eeis, expression_df, metadata_df, cpm_threshold,
                              test_variable, output_dir, coord_to_exon):
    """Generate association plots for top EEIs."""
    
    os.makedirs(f"{output_dir}/association_plots", exist_ok=True)
    
    for idx, eei_row in top_eeis.iterrows():
        try:
            plot_eei_association(eei_row, expression_df, metadata_df, 
                                cpm_threshold, test_variable, output_dir, coord_to_exon)
        except Exception as e:
            print(f"Error plotting EEI {eei_row['mouse_exon1']}-{eei_row['mouse_exon2']}: {e}")

def plot_eei_association(eei_row, expression_df, metadata_df, cpm_threshold,
                        test_variable, output_dir, coord_to_exon):
    """Create box plots showing EEI expression differences across treatment groups."""
    
    mouse_exon1 = eei_row['mouse_exon1']
    mouse_exon2 = eei_row['mouse_exon2']
    
    # Get expression data
    exon1_expr = find_expression_with_mapping(mouse_exon1, expression_df, coord_to_exon)
    exon2_expr = find_expression_with_mapping(mouse_exon2, expression_df, coord_to_exon)
    
    if exon1_expr is None or exon2_expr is None:
        return
    
    # Create DataFrame for plotting
    plot_data = pd.DataFrame({
        'exon1_expression': exon1_expr,
        'exon2_expression': exon2_expr,
        'combined_expression': exon1_expr + exon2_expr,
        test_variable: metadata_df[test_variable]
    })
    
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot 1: Exon 1 expression
    sns.boxplot(data=plot_data, x=test_variable, y='exon1_expression', ax=axes[0])
    axes[0].set_title(f'Exon 1: {mouse_exon1}')
    axes[0].set_ylabel('Expression (CPM)')
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Exon 2 expression
    sns.boxplot(data=plot_data, x=test_variable, y='exon2_expression', ax=axes[1])
    axes[1].set_title(f'Exon 2: {mouse_exon2}')
    axes[1].set_ylabel('Expression (CPM)')
    axes[1].grid(True, alpha=0.3)
    
    # Plot 3: Combined expression
    sns.boxplot(data=plot_data, x=test_variable, y='combined_expression', ax=axes[2])
    axes[2].set_title('Combined Expression (Exon1 + Exon2)')
    axes[2].set_ylabel('Expression (CPM)')
    axes[2].grid(True, alpha=0.3)
    
    # Add p-values to titles if available
    if 'exon1_ttest_p' in eei_row:
        axes[0].set_title(f'Exon 1: {mouse_exon1}\np={eei_row["exon1_ttest_p"]:.2e}')
    if 'exon2_ttest_p' in eei_row:
        axes[1].set_title(f'Exon 2: {mouse_exon2}\np={eei_row["exon2_ttest_p"]:.2e}')
    if 'combined_ttest_p' in eei_row:
        axes[2].set_title(f'Combined Expression\np={eei_row["combined_ttest_p"]:.2e}')
    
    # Rotate x-labels if needed
    for ax in axes:
        if len(plot_data[test_variable].unique()) > 4:
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    
    plt.suptitle(f'EEI: {mouse_exon1} - {mouse_exon2}\nOverall p-value: {eei_row["min_p_value"]:.2e}', 
                 fontsize=14, y=0.95)
    
    plt.tight_layout()
    
    # Save plot
    filename = f"eei_quantitative_{mouse_exon1.replace(':', '_')}_{mouse_exon2.replace(':', '_')}.png"
    filepath = f"{output_dir}/association_plots/{filename}"
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Main function to run the EEI treatment association analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Run EEI treatment association analysis')
    parser.add_argument('--mapped-eeis', required=True, 
                       help='Path to mapped mouse EEIs file')
    parser.add_argument('--expression', required=True,
                       help='Path to exon expression data (CPM)')
    parser.add_argument('--metadata', required=True,
                       help='Path to sample metadata/samplesheet')
    parser.add_argument('--mapping-files-dir', required=True,
                       help='Directory containing exon mapping files')
    parser.add_argument('--output-dir', default='outputs',
                       help='Output directory for results')
    parser.add_argument('--cpm-threshold', type=float, default=0.5,
                       help='CPM threshold for exon expression')
    parser.add_argument('--p-threshold', type=float, default=0.05,
                       help='P-value threshold for significance')
    parser.add_argument('--min-group-size', type=int, default=3,
                       help='Minimum samples per group for testing')
    parser.add_argument('--test-variable', default='treatment',
                       choices=['treatment', 'treatment_binary', 'model', 'timepoint'],
                       help='Which metadata variable to test')
    
    args = parser.parse_args()
    
    # Check if mapping files exist
    mapping_file = os.path.join(args.mapping_files_dir, "exon_coord_mapping_mouse.tsv")
    if not os.path.exists(mapping_file):
        print(f"Error: Mouse exon mapping file not found in {args.mapping_files_dir}")
        return
    
    print("=== EEI TREATMENT ASSOCIATION ANALYSIS ===")
    print(f"Mapped EEIs file: {args.mapped_eeis}")
    print(f"Expression file: {args.expression}")
    print(f"Metadata file: {args.metadata}")
    print(f"Test variable: {args.test_variable}")
    print("=" * 50)
    
    # Run the analysis
    try:
        significant_eeis, all_results, summary_stats = test_eei_treatment_association(
            mapped_eeis_file=args.mapped_eeis,
            expression_file=args.expression,
            metadata_file=args.metadata,
            cpm_threshold=args.cpm_threshold,
            p_value_threshold=args.p_threshold,
            min_group_size=args.min_group_size,
            output_dir=args.output_dir,
            mapping_files_dir=args.mapping_files_dir,
            test_variable=args.test_variable
        )
        
        if significant_eeis is not None:
            print(f"\nAnalysis completed successfully!")
            print(f"Results saved to: {args.output_dir}")
        
    except Exception as e:
        print(f"Error running analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()