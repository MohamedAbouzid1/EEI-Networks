import pandas as pd

def load_egio_data(egio_file):
    """Load EGIO orthologous exon mappings"""
    print(f"Loading orthology data from {egio_file}...")
    
    # Define column names based on the format you provided
    columns = ['Group', 'hsaEnsemblG', 'hsaPos', 'musEnsemblG', 'musPos', 'Iden', 'Type']
    
    # Load data
    df = pd.read_csv(egio_file, sep='\t', names=columns)
    
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
    
    # Load data
    df = pd.read_csv(eei_file, sep='\t', names=columns)
    
    # Clean up whitespace in string columns
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip()
    
    return df