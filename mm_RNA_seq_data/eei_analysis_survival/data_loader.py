import pandas as pd
import os
import numpy as np

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

