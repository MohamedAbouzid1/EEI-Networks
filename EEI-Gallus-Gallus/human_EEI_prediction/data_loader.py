import pandas as pd
from species_config import get_species_config

def load_egio_data(egio_file, species_code='dme'):
    """Load EGIO orthologous exon mappings"""
    print(f"Loading orthology data from {egio_file}...")
    
    # Get species-specific configuration
    config = get_species_config(species_code)
    print(f"Using configuration for {config['species_name']} ({config['species_code']})")
    
    # Use species-specific column names
    columns = config['egio_columns']
    
    # Load data with low_memory=False to avoid DtypeWarning
    df = pd.read_csv(egio_file, sep='\t', names=columns, low_memory=False)
    
    # Clean up whitespace in string columns
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip()
    
    return df

def load_eei_network(eei_file):
    """Load EEI network from PISA or EPPIC output"""
    print(f"Loading EEI network from {eei_file}...")
    
    # Define column names based on the format you provided
    columns = ['exon1', 'exon2', 'AA1', 'AA2', 'protein1', 'protein2', 
               'FreeEnergy', 'BuriedArea', 'Hydrogen', 'Disulphide', 
               'Saltbridge', 'Covalent', 'BuriedAreaAbs', 'SolAccAreaAbs', 
               'PDBID', 'allAA']
    
    # Load data with low_memory=False to avoid potential warnings
    df = pd.read_csv(eei_file, sep='\t', names=columns, low_memory=False)
    
    # Clean up whitespace in string columns
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip()
    
    return df