import pandas as pd

def create_coordinate_to_id_map(egio_df, species='mus'):
    """Create mapping from genomic coordinates to exon IDs"""
    coordinate_to_id = {}
    
    if species == 'mus':
        pos_col = 'musPos'
    else:  # Assuming human otherwise
        pos_col = 'hsaPos'
    
    for _, row in egio_df.iterrows():
        if row[pos_col] != 'None' and not pd.isna(row[pos_col]):
            coordinate_to_id[row[pos_col]] = f"Group_{row['Group']}"
    
    return coordinate_to_id

def create_orthology_map(egio_df, identity_threshold=0.8):
    """Create mapping between orthologous exons"""
    print(f"Creating orthology map with identity threshold: {identity_threshold}...")
    
    # Initialize orthology maps
    mus_to_hsa = {}
    hsa_to_mus = {}
    
    # Process each row
    for _, row in egio_df.iterrows():
        # Skip entries with None or low identity
        if (row['musPos'] == 'None' or row['hsaPos'] == 'None' or 
            pd.isna(row['musPos']) or pd.isna(row['hsaPos'])):
            continue
        
        try:
            identity = float(row['Iden'])
            if identity < identity_threshold:
                continue
        except:
            continue
            
        # Store mappings in both directions
        mus_to_hsa[row['musPos']] = {
            'hsaPos': row['hsaPos'],
            'identity': identity,
            'type': row['Type']
        }
        
        hsa_to_mus[row['hsaPos']] = {
            'musPos': row['musPos'],
            'identity': identity,
            'type': row['Type']
        }
    
    return {
        'mus_to_hsa': mus_to_hsa,
        'hsa_to_mus': hsa_to_mus
    }

def create_id_to_coordinate_map(exon_file):
    """Create mapping from exon IDs to genomic coordinates
    
    Args:
        exon_file: Path to file containing exon ID to coordinate mappings
                   (this could be a processed GTF file)
    
    Returns:
        Dictionary mapping exon IDs to coordinates
    """
    print(f"Loading exon ID to coordinate mappings from {exon_file}...")
    id_to_coordinate = {}
    
    # Read the exon file - format will depend on your data
    # Example format: exon_id\tcoordinate
    df = pd.read_csv(exon_file, sep='\t', header=None, names=['exon_id', 'coordinate'])
    
    for _, row in df.iterrows():
        id_to_coordinate[row['exon_id']] = row['coordinate']
    
    print(f"Loaded {len(id_to_coordinate)} exon ID to coordinate mappings")
    return id_to_coordinate