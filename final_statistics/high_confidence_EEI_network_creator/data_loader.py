import pandas as pd

def load_contact_network(filepath):
    """Load Contact-based network and extract EEI pairs."""
    print(f"Loading Contact network from: {filepath}")
    # Read as strings to avoid mixed types (e.g., floats from NaNs)
    df = pd.read_csv(filepath, sep='\t', dtype=str)
    
    # Extract unique EEI pairs (exon1-exon2 combinations)
    # Also keep protein information for reference
    eei_data = df[['exon1', 'exon2', 'protein1', 'protein1.1']].copy()
    eei_data.columns = ['exon1', 'exon2', 'protein1', 'protein2']
    
    # Drop rows with missing exon identifiers and strip whitespace
    eei_data = eei_data.dropna(subset=['exon1', 'exon2'])
    eei_data['exon1'] = eei_data['exon1'].str.strip()
    eei_data['exon2'] = eei_data['exon2'].str.strip()
    
    # Create a set of EEI pairs (sorted to handle bidirectional interactions)
    eei_pairs = set()
    protein_map = {}
    
    for _, row in eei_data.iterrows():
        # Sort exons to handle bidirectional interactions consistently
        exon_pair = tuple(sorted([row['exon1'], row['exon2']]))
        eei_pairs.add(exon_pair)
        
        # Store protein mapping
        protein_map[exon_pair] = (row['protein1'], row['protein2'])
    
    print(f"  Found {len(eei_pairs)} unique EEI pairs in Contact network")
    return eei_pairs, protein_map, df

def load_eppic_network(filepath):
    """Load EPPIC network and extract EEI pairs."""
    print(f"Loading EPPIC network from: {filepath}")
    df = pd.read_csv(filepath, sep='\t', dtype=str)
    
    # Extract unique EEI pairs
    eei_data = df[['exon1', 'exon2', 'protein1', 'protein2']].copy()
    
    # Drop rows with missing exon identifiers and strip whitespace
    eei_data = eei_data.dropna(subset=['exon1', 'exon2'])
    eei_data['exon1'] = eei_data['exon1'].str.strip()
    eei_data['exon2'] = eei_data['exon2'].str.strip()
    
    eei_pairs = set()
    protein_map = {}
    
    for _, row in eei_data.iterrows():
        exon_pair = tuple(sorted([row['exon1'], row['exon2']]))
        eei_pairs.add(exon_pair)
        protein_map[exon_pair] = (row['protein1'], row['protein2'])
    
    print(f"  Found {len(eei_pairs)} unique EEI pairs in EPPIC network")
    return eei_pairs, protein_map, df

def load_pisa_network(filepath):
    """Load PISA network and extract EEI pairs."""
    print(f"Loading PISA network from: {filepath}")
    df = pd.read_csv(filepath, sep='\t', dtype=str)
    
    # Extract unique EEI pairs
    eei_data = df[['exon1', 'exon2', 'protein1', 'protein2']].copy()
    
    # Drop rows with missing exon identifiers and strip whitespace
    eei_data = eei_data.dropna(subset=['exon1', 'exon2'])
    eei_data['exon1'] = eei_data['exon1'].str.strip()
    eei_data['exon2'] = eei_data['exon2'].str.strip()
    
    eei_pairs = set()
    protein_map = {}
    
    for _, row in eei_data.iterrows():
        exon_pair = tuple(sorted([row['exon1'], row['exon2']]))
        eei_pairs.add(exon_pair)
        protein_map[exon_pair] = (row['protein1'], row['protein2'])
    
    print(f"  Found {len(eei_pairs)} unique EEI pairs in PISA network")
    return eei_pairs, protein_map, df