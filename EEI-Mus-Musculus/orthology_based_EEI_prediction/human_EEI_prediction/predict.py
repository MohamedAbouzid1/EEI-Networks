import pandas as pd

def predict_novel_human_eeis_from_mouse(mouse_eei_df, egio_df, human_eei_df=None):
    """
    Predict novel human EEIs based on mouse EEI data and EGIO orthology mapping.
    
    Args:
        mouse_eei_df: DataFrame containing mouse EEI data with genomic coordinate format for exon1/exon2
        egio_df: DataFrame with EGIO orthology mapping between mouse and human exons
        human_eei_df: Optional DataFrame of known human EEIs to filter out
    
    Returns:
        DataFrame of predicted novel human EEIs
    """
    print("Predicting novel human EEIs based on mouse data...")
    
    # Create a mapping from mouse coordinates to human coordinates
    mus_to_hsa = {}
    for _, row in egio_df.iterrows():
        if row['musPos'] != 'None' and row['hsaPos'] != 'None' and row['musPos'] != '' and row['hsaPos'] != '':
            # Map based on genomic coordinates
            mouse_coord = row['musPos']
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
                
            mus_to_hsa[mouse_coord] = {
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
    
    # Process mouse EEIs
    predicted_human_eeis = []
    
    for _, eei in mouse_eei_df.iterrows():
        mouse_exon1 = eei['exon1']  # This is now a genomic coordinate
        mouse_exon2 = eei['exon2']  # This is now a genomic coordinate
        
        # Check if both mouse exons have human orthologs
        if mouse_exon1 in mus_to_hsa and mouse_exon2 in mus_to_hsa:
            human_exon1 = mus_to_hsa[mouse_exon1]['hsaPos']
            human_exon2 = mus_to_hsa[mouse_exon2]['hsaPos']
            
            # Skip if either human exon is empty
            if not human_exon1 or not human_exon2:
                continue
                
            # Calculate confidence score
            identity1 = mus_to_hsa[mouse_exon1]['identity']
            identity2 = mus_to_hsa[mouse_exon2]['identity']
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
                    'mouse_exon1': mouse_exon1,
                    'mouse_exon2': mouse_exon2,
                    'confidence': confidence,
                    'identity1': identity1,
                    'identity2': identity2,
                    'type1': mus_to_hsa[mouse_exon1]['type'],
                    'type2': mus_to_hsa[mouse_exon2]['type'],
                }
                
                # Copy all available fields from the mouse EEI data
                # This allows the code to work with different formats (PISA, CONTACT, etc.)
                for column in mouse_eei_df.columns:
                    if column not in ['exon1', 'exon2']:  # Skip exon coordinates as they're already handled
                        if column in eei:
                            predicted_eei[column] = eei[column]
                
                predicted_human_eeis.append(predicted_eei)
    
    # Convert to DataFrame and sort by confidence
    result_df = pd.DataFrame(predicted_human_eeis)
    if not result_df.empty:
        result_df = result_df.sort_values('confidence', ascending=False)
    
    print(f"Predicted {len(result_df)} novel human EEIs based on mouse data")
    return result_df