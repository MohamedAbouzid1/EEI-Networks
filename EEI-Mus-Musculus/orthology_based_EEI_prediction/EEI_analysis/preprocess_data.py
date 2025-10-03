import pandas as pd
import numpy as np
import os
import re

def fix_column_formats(df, file_type):
    """
    Fix common format issues in the dataframes
    
    Args:
        df: DataFrame to fix
        file_type: Type of file ('human_eei', 'mouse_eei', 'egio', 'predicted')
        
    Returns:
        Fixed DataFrame
    """
    # Make a copy to avoid modifying the original
    fixed_df = df.copy()
    
    # Handle common issues based on file type
    if file_type in ['human_eei', 'mouse_eei']:
        # Ensure coverage columns are numeric
        for col in ['exon1_coverage_percent', 'exon2_coverage_percent', 'jaccard_percent']:
            if col in fixed_df.columns:
                fixed_df[col] = pd.to_numeric(fixed_df[col], errors='coerce')
                
        # Ensure length columns are numeric
        for col in ['exon1_length', 'exon2_length', 'exon1_coverage', 'exon2_coverage']:
            if col in fixed_df.columns:
                fixed_df[col] = pd.to_numeric(fixed_df[col], errors='coerce')
    
    elif file_type == 'egio':
        # Convert identity to numeric
        if 'Iden' in fixed_df.columns:
            fixed_df['Iden'] = pd.to_numeric(fixed_df['Iden'], errors='coerce')
            
        # Handle 'None' values in position columns
        for col in ['hsaPos', 'musPos']:
            if col in fixed_df.columns:
                fixed_df[col] = fixed_df[col].apply(lambda x: None if x == 'None' else x)
    
    elif file_type == 'predicted':
        # Ensure confidence columns are numeric
        for col in ['confidence', 'identity1', 'identity2']:
            if col in fixed_df.columns:
                fixed_df[col] = pd.to_numeric(fixed_df[col], errors='coerce')
    
    return fixed_df

def standardize_exon_ids(df, columns=None):
    """
    Standardize exon IDs to a common format
    
    Args:
        df: DataFrame with exon IDs
        columns: List of column names containing exon IDs (default: ['exon1', 'exon2'])
        
    Returns:
        DataFrame with standardized exon IDs
    """
    if columns is None:
        columns = ['exon1', 'exon2']
    
    # Make a copy to avoid modifying the original
    std_df = df.copy()
    
    # Function to standardize IDs
    def standardize_id(exon_id):
        if exon_id is None or pd.isna(exon_id):
            return exon_id
            
        # Convert string to string
        exon_id = str(exon_id)
        
        # Check if it's already in standard format (ENSE00000000123)
        if re.match(r'^ENS[A-Z]*E\d+$', exon_id):
            return exon_id
            
        # Check if it's a coordinate (chr1:12345:67890:1)
        if re.match(r'^chr[\dXY]+:\d+:\d+:[-+\d]+$', exon_id):
            return exon_id
            
        # Otherwise, return as is
        return exon_id
    
    # Apply to specified columns
    for col in columns:
        if col in std_df.columns:
            std_df[col] = std_df[col].apply(standardize_id)
    
    return std_df

def check_data_consistency(human_df, mouse_df, egio_df, predicted_df=None):
    """
    Check consistency across datasets and report issues
    
    Args:
        human_df: Human EEI DataFrame
        mouse_df: Mouse EEI DataFrame
        egio_df: EGIO orthology DataFrame
        predicted_df: Optional predicted EEIs DataFrame
        
    Returns:
        Dictionary with consistency check results
    """
    results = {
        'warnings': [],
        'info': [],
        'is_consistent': True
    }
    
    # Check human EEI data
    if human_df is not None:
        human_exon_count = len(set(human_df['exon1'].unique()).union(set(human_df['exon2'].unique())))
        results['info'].append(f"Human EEI network: {len(human_df)} interactions between {human_exon_count} unique exons")
        
        # Check for missing required columns
        required_cols = ['exon1', 'exon2']
        missing_cols = [col for col in required_cols if col not in human_df.columns]
        if missing_cols:
            results['warnings'].append(f"Human EEI data is missing required columns: {missing_cols}")
            results['is_consistent'] = False
    else:
        results['warnings'].append("Human EEI data is missing")
        results['is_consistent'] = False
    
    # Check mouse EEI data
    if mouse_df is not None:
        mouse_exon_count = len(set(mouse_df['exon1'].unique()).union(set(mouse_df['exon2'].unique())))
        results['info'].append(f"Mouse EEI network: {len(mouse_df)} interactions between {mouse_exon_count} unique exons")
        
        # Check for missing required columns
        required_cols = ['exon1', 'exon2']
        missing_cols = [col for col in required_cols if col not in mouse_df.columns]
        if missing_cols:
            results['warnings'].append(f"Mouse EEI data is missing required columns: {missing_cols}")
            results['is_consistent'] = False
    else:
        results['warnings'].append("Mouse EEI data is missing")
        results['is_consistent'] = False
    
    # Check EGIO data
    if egio_df is not None:
        results['info'].append(f"EGIO orthology mappings: {len(egio_df)} total mappings")
        
        # Check for missing required columns
        required_cols = ['hsaPos', 'musPos', 'Iden', 'Type']
        missing_cols = [col for col in required_cols if col not in egio_df.columns]
        if missing_cols:
            results['warnings'].append(f"EGIO data is missing required columns: {missing_cols}")
            results['is_consistent'] = False
            
        # Check for valid identity values
        if 'Iden' in egio_df.columns:
            egio_df['Iden'] = pd.to_numeric(egio_df['Iden'], errors='coerce')
            invalid_iden = egio_df['Iden'].isna().sum()
            if invalid_iden > 0:
                results['warnings'].append(f"EGIO data contains {invalid_iden} invalid identity values")
                results['is_consistent'] = False
            
            # Report identity distribution
            min_iden = egio_df['Iden'].min()
            max_iden = egio_df['Iden'].max()
            mean_iden = egio_df['Iden'].mean()
            results['info'].append(f"EGIO identity range: {min_iden:.2f} to {max_iden:.2f}, mean: {mean_iden:.2f}")
    else:
        results['warnings'].append("EGIO orthology data is missing")
        results['is_consistent'] = False
    
    # Check predicted EEI data if provided
    if predicted_df is not None:
        pred_exon_count = len(set(predicted_df['exon1'].unique()).union(set(predicted_df['exon2'].unique())))
        results['info'].append(f"Predicted EEI network: {len(predicted_df)} interactions between {pred_exon_count} unique exons")
        
        # Check for missing required columns
        required_cols = ['exon1', 'exon2']
        missing_cols = [col for col in required_cols if col not in predicted_df.columns]
        if missing_cols:
            results['warnings'].append(f"Predicted EEI data is missing required columns: {missing_cols}")
            results['is_consistent'] = False
            
        # Check for valid confidence values
        if 'confidence' in predicted_df.columns:
            predicted_df['confidence'] = pd.to_numeric(predicted_df['confidence'], errors='coerce')
            invalid_conf = predicted_df['confidence'].isna().sum()
            if invalid_conf > 0:
                results['warnings'].append(f"Predicted EEI data contains {invalid_conf} invalid confidence values")
                results['is_consistent'] = False
            
            # Report confidence distribution
            min_conf = predicted_df['confidence'].min()
            max_conf = predicted_df['confidence'].max()
            mean_conf = predicted_df['confidence'].mean()
            results['info'].append(f"Prediction confidence range: {min_conf:.2f} to {max_conf:.2f}, mean: {mean_conf:.2f}")
    
    # Check for common exons between human EEI and EGIO
    if human_df is not None and egio_df is not None:
        human_exons = set(human_df['exon1'].unique()).union(set(human_df['exon2'].unique()))
        egio_human_exons = set(egio_df['hsaPos'].dropna().unique())
        common_human_exons = human_exons.intersection(egio_human_exons)
        
        if len(common_human_exons) == 0:
            results['warnings'].append("No common exons found between human EEI data and EGIO human exons")
            results['is_consistent'] = False
        else:
            results['info'].append(f"Common human exons: {len(common_human_exons)} ({len(common_human_exons)/len(human_exons)*100:.1f}% of human EEI exons)")
    
    # Check for common exons between mouse EEI and EGIO
    if mouse_df is not None and egio_df is not None:
        mouse_exons = set(mouse_df['exon1'].unique()).union(set(mouse_df['exon2'].unique()))
        egio_mouse_exons = set(egio_df['musPos'].dropna().unique())
        common_mouse_exons = mouse_exons.intersection(egio_mouse_exons)
        
        if len(common_mouse_exons) == 0:
            results['warnings'].append("No common exons found between mouse EEI data and EGIO mouse exons")
            results['is_consistent'] = False
        else:
            results['info'].append(f"Common mouse exons: {len(common_mouse_exons)} ({len(common_mouse_exons)/len(mouse_exons)*100:.1f}% of mouse EEI exons)")
    
    return results

def filter_high_quality_data(human_df, mouse_df, egio_df, min_identity=0.8, min_coverage=10):
    """
    Filter data to include only high quality entries
    
    Args:
        human_df: Human EEI DataFrame
        mouse_df: Mouse EEI DataFrame
        egio_df: EGIO orthology DataFrame
        min_identity: Minimum identity score for orthologous exons
        min_coverage: Minimum interface coverage percentage
        
    Returns:
        Dictionary with filtered DataFrames
    """
    # Filter EGIO data by identity
    filtered_egio = egio_df.copy()
    filtered_egio['Iden'] = pd.to_numeric(filtered_egio['Iden'], errors='coerce')
    filtered_egio = filtered_egio[filtered_egio['Iden'] >= min_identity]
    
    # Filter human EEI data by coverage
    filtered_human = human_df.copy()
    if all(col in filtered_human.columns for col in ['exon1_coverage_percent', 'exon2_coverage_percent']):
        filtered_human['exon1_coverage_percent'] = pd.to_numeric(filtered_human['exon1_coverage_percent'], errors='coerce')
        filtered_human['exon2_coverage_percent'] = pd.to_numeric(filtered_human['exon2_coverage_percent'], errors='coerce')
        filtered_human = filtered_human[
            (filtered_human['exon1_coverage_percent'] >= min_coverage) & 
            (filtered_human['exon2_coverage_percent'] >= min_coverage)
        ]
    
    # Filter mouse EEI data by coverage
    filtered_mouse = mouse_df.copy()
    if all(col in filtered_mouse.columns for col in ['exon1_coverage_percent', 'exon2_coverage_percent']):
        filtered_mouse['exon1_coverage_percent'] = pd.to_numeric(filtered_mouse['exon1_coverage_percent'], errors='coerce')
        filtered_mouse['exon2_coverage_percent'] = pd.to_numeric(filtered_mouse['exon2_coverage_percent'], errors='coerce')
        filtered_mouse = filtered_mouse[
            (filtered_mouse['exon1_coverage_percent'] >= min_coverage) & 
            (filtered_mouse['exon2_coverage_percent'] >= min_coverage)
        ]
    
    filtered = {
        'human_eei': filtered_human,
        'mouse_eei': filtered_mouse,
        'egio': filtered_egio
    }
    
    # Print stats about filtering
    print(f"EGIO: {len(filtered_egio)}/{len(egio_df)} entries retained ({len(filtered_egio)/len(egio_df)*100:.1f}%)")
    print(f"Human EEI: {len(filtered_human)}/{len(human_df)} entries retained ({len(filtered_human)/len(human_df)*100:.1f}%)")
    print(f"Mouse EEI: {len(filtered_mouse)}/{len(mouse_df)} entries retained ({len(filtered_mouse)/len(mouse_df)*100:.1f}%)")
    
    return filtered

def save_processed_data(dataframes, output_dir='processed_data'):
    """
    Save processed dataframes to files
    
    Args:
        dataframes: Dictionary of dataframes to save
        output_dir: Directory to save the files
        
    Returns:
        Dictionary with paths to saved files
    """
    os.makedirs(output_dir, exist_ok=True)
    
    saved_paths = {}
    
    for name, df in dataframes.items():
        filename = f"{name}_processed.tsv"
        filepath = os.path.join(output_dir, filename)
        df.to_csv(filepath, sep='\t', index=False)
        saved_paths[name] = filepath
        print(f"Saved {name} data to {filepath}")
    
    return saved_paths

def preprocess_and_check_data(human_eei_file, mouse_eei_file, egio_file, predicted_file=None, 
                             output_dir='processed_data', min_identity=0.8, min_coverage=10):
    """
    Main function to process and check all data files
    
    Args:
        human_eei_file: Path to human EEI network file
        mouse_eei_file: Path to mouse EEI network file
        egio_file: Path to EGIO orthology mapping file
        predicted_file: Optional path to predicted EEIs file
        output_dir: Directory to save processed data
        min_identity: Minimum identity score for orthologous exons
        min_coverage: Minimum interface coverage percentage
        
    Returns:
        Dictionary with paths to processed files and consistency check results
    """
    print("Loading data files...")
    
    # Load data
    human_df = pd.read_csv(human_eei_file, sep='\t')
    mouse_df = pd.read_csv(mouse_eei_file, sep='\t')
    egio_df = pd.read_csv(egio_file, sep='\t')
    
    predicted_df = None
    if predicted_file:
        predicted_df = pd.read_csv(predicted_file, sep='\t')
    
    print(f"Loaded {len(human_df)} human EEIs, {len(mouse_df)} mouse EEIs, and {len(egio_df)} orthology mappings")
    if predicted_df is not None:
        print(f"Loaded {len(predicted_df)} predicted EEIs")
    
    # Fix column formats
    print("Fixing column formats...")
    human_df = fix_column_formats(human_df, 'human_eei')
    mouse_df = fix_column_formats(mouse_df, 'mouse_eei')
    egio_df = fix_column_formats(egio_df, 'egio')
    
    if predicted_df is not None:
        predicted_df = fix_column_formats(predicted_df, 'predicted')
    
    # Standardize exon IDs
    print("Standardizing exon IDs...")
    human_df = standardize_exon_ids(human_df)
    mouse_df = standardize_exon_ids(mouse_df)
    
    # EGIO has different column names for exon positions
    egio_df = standardize_exon_ids(egio_df, columns=['hsaPos', 'musPos'])
    
    if predicted_df is not None:
        predicted_df = standardize_exon_ids(predicted_df)
    
    # Check data consistency
    print("Checking data consistency...")
    consistency_results = check_data_consistency(human_df, mouse_df, egio_df, predicted_df)
    
    # Print warnings and info
    for warning in consistency_results['warnings']:
        print(f"WARNING: {warning}")
    
    for info in consistency_results['info']:
        print(f"INFO: {info}")
    
    if not consistency_results['is_consistent']:
        print("WARNING: Data consistency check failed. Results may be unreliable.")
    else:
        print("Data consistency check passed.")
    
    # Filter data for high quality
    print(f"Filtering data (min_identity={min_identity}, min_coverage={min_coverage})...")
    filtered_data = filter_high_quality_data(human_df, mouse_df, egio_df, min_identity, min_coverage)
    
    # Prepare dictionary of all dataframes
    all_data = {
        'human_eei': human_df,
        'mouse_eei': mouse_df,
        'egio': egio_df,
        'human_eei_filtered': filtered_data['human_eei'],
        'mouse_eei_filtered': filtered_data['mouse_eei'],
        'egio_filtered': filtered_data['egio']
    }
    
    if predicted_df is not None:
        all_data['predicted'] = predicted_df
    
    # Save processed data
    print(f"Saving processed data to {output_dir}...")
    saved_paths = save_processed_data(all_data, output_dir)
    
    # Return results
    return {
        'paths': saved_paths,
        'consistency': consistency_results,
        'data': all_data
    }

def generate_summary_report(result, output_file='data_summary.txt'):
    """
    Generate a summary report of the data preprocessing
    
    Args:
        result: Result dictionary from preprocess_and_check_data
        output_file: Path to save the summary report
    """
    with open(output_file, 'w') as f:
        f.write("===== EEI Data Preprocessing Summary =====\n\n")
        
        # Write data counts
        f.write("Data Counts:\n")
        f.write(f"  Human EEI network: {len(result['data']['human_eei'])} interactions\n")
        f.write(f"  Mouse EEI network: {len(result['data']['mouse_eei'])} interactions\n")
        f.write(f"  EGIO orthology mappings: {len(result['data']['egio'])} mappings\n")
        if 'predicted' in result['data']:
            f.write(f"  Predicted EEIs: {len(result['data']['predicted'])} interactions\n")
        f.write("\n")
        
        # Write filtering results
        f.write("Filtering Results:\n")
        f.write(f"  Human EEI filtered: {len(result['data']['human_eei_filtered'])}/{len(result['data']['human_eei'])} retained "
                f"({len(result['data']['human_eei_filtered'])/len(result['data']['human_eei'])*100:.1f}%)\n")
        f.write(f"  Mouse EEI filtered: {len(result['data']['mouse_eei_filtered'])}/{len(result['data']['mouse_eei'])} retained "
                f"({len(result['data']['mouse_eei_filtered'])/len(result['data']['mouse_eei'])*100:.1f}%)\n")
        f.write(f"  EGIO filtered: {len(result['data']['egio_filtered'])}/{len(result['data']['egio'])} retained "
                f"({len(result['data']['egio_filtered'])/len(result['data']['egio'])*100:.1f}%)\n")
        f.write("\n")
        
        # Write consistency warnings
        f.write("Consistency Warnings:\n")
        if not result['consistency']['warnings']:
            f.write("  No warnings\n")
        else:
            for i, warning in enumerate(result['consistency']['warnings'], 1):
                f.write(f"  {i}. {warning}\n")
        f.write("\n")
        
        # Write file paths
        f.write("Processed File Paths:\n")
        for name, path in result['paths'].items():
            f.write(f"  {name}: {path}\n")
    
    print(f"Summary report saved to {output_file}")
    return output_file

def load_exon_mappings(human_mapping_file, mouse_mapping_file):
    """
    Load exon ID to coordinate mapping files
    
    Args:
        human_mapping_file: Path to human exon mapping file
        mouse_mapping_file: Path to mouse exon mapping file
        
    Returns:
        Dictionary with mapping dictionaries
    """
    print(f"Loading exon ID mappings from {human_mapping_file} and {mouse_mapping_file}")
    
    # Load human mappings
    human_df = pd.read_csv(human_mapping_file, sep='\t')
    human_id_to_coord = dict(zip(human_df['exon_id'], human_df['coord']))
    human_coord_to_id = dict(zip(human_df['coord'], human_df['exon_id']))
    
    print(f"Loaded {len(human_id_to_coord)} human exon ID mappings")
    
    # Load mouse mappings
    mouse_df = pd.read_csv(mouse_mapping_file, sep='\t')
    mouse_id_to_coord = dict(zip(mouse_df['exon_id'], mouse_df['coord']))
    mouse_coord_to_id = dict(zip(mouse_df['coord'], mouse_df['exon_id']))
    
    print(f"Loaded {len(mouse_id_to_coord)} mouse exon ID mappings")
    
    return {
        'human_id_to_coord': human_id_to_coord,
        'human_coord_to_id': human_coord_to_id,
        'mouse_id_to_coord': mouse_id_to_coord,
        'mouse_coord_to_id': mouse_coord_to_id
    }

def convert_eei_to_coordinates(eei_df, id_to_coord_map, columns=['exon1', 'exon2']):
    """
    Convert EEI dataframe to use coordinate format for exon IDs
    
    Args:
        eei_df: EEI dataframe with exon IDs
        id_to_coord_map: Dictionary mapping exon IDs to coordinates
        columns: Columns containing exon IDs
        
    Returns:
        Converted dataframe with coordinates
    """
    # Create a copy to avoid modifying the original
    converted_df = eei_df.copy()
    
    # Track how many IDs were successfully converted
    conversion_stats = {'total': 0, 'converted': 0}
    
    # Convert each column
    for col in columns:
        if col in converted_df.columns:
            # Count total IDs
            conversion_stats['total'] += len(converted_df[col].dropna())
            
            # Apply conversion
            converted_df[col] = converted_df[col].apply(
                lambda x: id_to_coord_map.get(x, x) if pd.notna(x) else x
            )
            
            # Count successful conversions
            conversion_stats['converted'] += sum(converted_df[col].isin(id_to_coord_map.values()))
    
    # Calculate conversion rate
    if conversion_stats['total'] > 0:
        conversion_rate = conversion_stats['converted'] / conversion_stats['total'] * 100
        print(f"Converted {conversion_stats['converted']}/{conversion_stats['total']} "
              f"exon IDs to coordinates ({conversion_rate:.1f}%)")
    
    return converted_df

if __name__ == "__main__":
    # Example usage
    result = preprocess_and_check_data(
        human_eei_file="path/to/human_eei_network.tsv",
        mouse_eei_file="path/to/mouse_eei_network.tsv",
        egio_file="path/to/egio_output.tsv",
        predicted_file="path/to/predicted_eeis.tsv",
        output_dir="processed_data",
        min_identity=0.8,
        min_coverage=5
    )
    
    generate_summary_report(result, "data_summary.txt")