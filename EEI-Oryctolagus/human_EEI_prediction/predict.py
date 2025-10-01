import pandas as pd
from species_config import get_species_config

def predict_novel_human_eeis_from_ortholog(ortholog_eei_df, egio_df, human_eei_df=None, species_code='dme'):
    """
    Predict novel human EEIs based on orthologous species EEI data and EGIO orthology mapping.
    
    Args:
        ortholog_eei_df: DataFrame containing orthologous species EEI data with genomic coordinate format for exon1/exon2
        egio_df: DataFrame with EGIO orthology mapping between orthologous species and human exons
        human_eei_df: Optional DataFrame of known human EEIs to filter out
        species_code: Species code for the orthologous species (e.g., 'dme', 'mus', 'rno')
    
    Returns:
        DataFrame of predicted novel human EEIs
    """
    print("Predicting novel human EEIs based on orthologous species data...")
    
    # Get species-specific configuration
    config = get_species_config(species_code)
    ortholog_pos_col = config['ortholog_pos_col']
    output_cols = config['output_columns']
    
    # Create a mapping from orthologous species coordinates to human coordinates
    ortholog_to_hsa = {}
    for _, row in egio_df.iterrows():
        if row[ortholog_pos_col] != 'None' and row['hsaPos'] != 'None' and row[ortholog_pos_col] != '' and row['hsaPos'] != '':
            # Map based on genomic coordinates
            ortholog_coord = row[ortholog_pos_col]
            human_coord = row['hsaPos']
            
            try:
                identity = float(row['Iden'])
                if identity < 0.0:  # Skip low identity mappings - use your threshold here
                    continue
            except (ValueError, TypeError):
                # Skip if identity is not a valid float
                continue
                
            # Only include 1-1 mappings
            if row['Type'] != '1-1':
                continue
                
            ortholog_to_hsa[ortholog_coord] = {
                'hsaPos': human_coord,
                'identity': identity,
                'type': row['Type']
            }
    
    # Set of existing human EEI pairs for filtering
    existing_eei_pairs = set()
    if human_eei_df is not None:
        for _, row in human_eei_df.iterrows():
            if pd.notna(row['exon1']) and pd.notna(row['exon2']):
                pair = (row['exon1'], row['exon2'])
                existing_eei_pairs.add(pair)
                existing_eei_pairs.add((row['exon2'], row['exon1']))  # Add reverse pair too
    
    # Process orthologous species EEIs
    predicted_human_eeis = []
    
    for _, eei in ortholog_eei_df.iterrows():
        ortholog_exon1 = eei['exon1']  # This is now a genomic coordinate
        ortholog_exon2 = eei['exon2']  # This is now a genomic coordinate
        
        # Check if both orthologous species exons have human orthologs
        if ortholog_exon1 in ortholog_to_hsa and ortholog_exon2 in ortholog_to_hsa:
            human_exon1 = ortholog_to_hsa[ortholog_exon1]['hsaPos']
            human_exon2 = ortholog_to_hsa[ortholog_exon2]['hsaPos']
            
            # Skip if either human exon is empty
            if not human_exon1 or not human_exon2:
                continue
                
            # Calculate confidence score
            identity1 = ortholog_to_hsa[ortholog_exon1]['identity']
            identity2 = ortholog_to_hsa[ortholog_exon2]['identity']
            confidence = (identity1 * identity2) ** 0.5  # Geometric mean
            
            # Skip if confidence is too low
            #if confidence < 0.8:  # Adjust threshold as needed
            #    continue
                
            # Check if this is a novel human EEI
            if (human_exon1, human_exon2) not in existing_eei_pairs:
                # Create a dictionary to hold the predicted EEI data
                predicted_eei = {
                    'exon1': human_exon1,
                    'exon2': human_exon2,
                    output_cols['ortholog_exon1']: ortholog_exon1,
                    output_cols['ortholog_exon2']: ortholog_exon2,
                    'confidence': confidence,
                    'identity1': identity1,
                    'identity2': identity2,
                    'type1': ortholog_to_hsa[ortholog_exon1]['type'],
                    'type2': ortholog_to_hsa[ortholog_exon2]['type'],
                }
                
                # Copy all available fields from the orthologous species EEI data
                # This allows the code to work with different formats (PISA, CONTACT, etc.)
                for column in ortholog_eei_df.columns:
                    if column not in ['exon1', 'exon2']:  # Skip exon coordinates as they're already handled
                        if column in eei:
                            predicted_eei[column] = eei[column]
                
                predicted_human_eeis.append(predicted_eei)
    
    # Convert to DataFrame and sort by confidence
    result_df = pd.DataFrame(predicted_human_eeis)
    if not result_df.empty:
        result_df = result_df.sort_values('confidence', ascending=False)
    
    print(f"Predicted {len(result_df)} novel human EEIs based on {config['species_name']} data")
    return result_df