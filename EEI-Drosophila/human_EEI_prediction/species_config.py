# Species configuration for the EEI prediction pipeline
# This file allows easy switching between different species

# Configuration for Drosophila melanogaster (dme)
DROSOPHILA_CONFIG = {
    'species_code': 'dme',
    'species_name': 'Drosophila melanogaster',
    'egio_columns': ['Group', 'hsaEnsemblG', 'hsaPos', 'dmeEnsemblG', 'dmePos', 'Iden', 'Type'],
    'ortholog_pos_col': 'dmePos',
    'output_columns': {
        'ortholog_exon1': 'drosophila_exon1',
        'ortholog_exon2': 'drosophila_exon2'
    }
}

# Configuration for Mus musculus (mus) - for reference
MOUSE_CONFIG = {
    'species_code': 'mus',
    'species_name': 'Mus musculus',
    'egio_columns': ['Group', 'hsaEnsemblG', 'hsaPos', 'musEnsemblG', 'musPos', 'Iden', 'Type'],
    'ortholog_pos_col': 'musPos',
    'output_columns': {
        'ortholog_exon1': 'mouse_exon1',
        'ortholog_exon2': 'mouse_exon2'
    }
}

# Configuration for Rattus norvegicus (rno) - for reference
RAT_CONFIG = {
    'species_code': 'rno',
    'species_name': 'Rattus norvegicus',
    'egio_columns': ['Group', 'hsaEnsemblG', 'hsaPos', 'rnoEnsemblG', 'rnoPos', 'Iden', 'Type'],
    'ortholog_pos_col': 'rnoPos',
    'output_columns': {
        'ortholog_exon1': 'rat_exon1',
        'ortholog_exon2': 'rat_exon2'
    }
}

# Default configuration (currently set to Drosophila)
DEFAULT_CONFIG = DROSOPHILA_CONFIG

def get_species_config(species_code=None):
    """
    Get species configuration based on species code.
    
    Args:
        species_code (str): Species code (e.g., 'dme', 'mus', 'rno')
        
    Returns:
        dict: Species configuration
    """
    if species_code == 'dme':
        return DROSOPHILA_CONFIG
    elif species_code == 'mus':
        return MOUSE_CONFIG
    elif species_code == 'rno':
        return RAT_CONFIG
    else:
        print(f"Warning: Unknown species code '{species_code}'. Using default configuration.")
        return DEFAULT_CONFIG

def print_available_species():
    """Print all available species configurations."""
    print("Available species configurations:")
    print("  dme - Drosophila melanogaster")
    print("  mus - Mus musculus (Mouse)")
    print("  rno - Rattus norvegicus (Rat)")
    print("\nTo use a different species, update the species_config.py file")
    print("or modify the data_loader.py to use the appropriate column names.")
