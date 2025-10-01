import pandas as pd
from utils import parse_arguments
from mapping import create_orthology_map, create_coordinate_to_id_map
from data_loader import load_egio_data, load_eei_network
from predict import predict_novel_human_eeis_from_ortholog
import os

def create_id_to_coordinate_map(exon_map_file):
    """Create mapping from exon IDs to genomic coordinates"""
    print(f"Loading exon ID to coordinate mappings from {exon_map_file}...")
    id_to_coordinate = {}
    
    # Read the exon mapping file
    df = pd.read_csv(exon_map_file, sep='\t')
    
    for _, row in df.iterrows():
        id_to_coordinate[row['exon_id']] = row['coord']
    
    print(f"Loaded {len(id_to_coordinate)} exon ID to coordinate mappings")
    return id_to_coordinate

def create_coordinate_to_id_map(exon_map_file):
    """Create mapping from genomic coordinates to exon IDs"""
    print(f"Loading coordinate to exon ID mappings from {exon_map_file}...")
    coordinate_to_id = {}
    
    # Read the exon mapping file
    df = pd.read_csv(exon_map_file, sep='\t')
    
    for _, row in df.iterrows():
        coordinate_to_id[row['coord']] = row['exon_id']
    
    print(f"Loaded {len(coordinate_to_id)} coordinate to exon ID mappings")
    return coordinate_to_id

def main():
    args = parse_arguments()
    
    print("Starting EEI prediction pipeline...")
    
    # Determine species code from EGIO filename
    species_code = 'dme'  # Default to Drosophila
    if 'hsa_dme' in args.egio:
        species_code = 'dme'
    elif 'hsa_mus' in args.egio:
        species_code = 'mus'
    elif 'hsa_rno' in args.egio:
        species_code = 'rno'
    
    print(f"Detected species code: {species_code}")
    
    # Load data
    egio_df = load_egio_data(args.egio, species_code)
    ortholog_eei_df = pd.read_csv(args.eei, sep='\t')
    
    # Load optional human EEI data if provided
    human_eei_df = None
    if args.human_eei:
        human_eei_df = pd.read_csv(args.human_eei, sep='\t')
        print(f"Loaded {len(human_eei_df)} existing human EEIs for filtering")
    
    # Create ID to coordinate mappings
    ortholog_id_to_coord = create_id_to_coordinate_map(args.ortholog_exon_map)
    human_id_to_coord = create_id_to_coordinate_map(args.human_exon_map)
    
    # Create coordinate to ID mappings
    ortholog_coord_to_id = create_coordinate_to_id_map(args.ortholog_exon_map)
    human_coord_to_id = create_coordinate_to_id_map(args.human_exon_map)
    
    # Create orthology map
    orthology_map = create_orthology_map(egio_df, args.identity_threshold, species_code)
    
    # Create a copy of the ortholog EEI dataframe with coordinates instead of IDs
    ortholog_eei_coords_df = ortholog_eei_df.copy()
    for i, row in ortholog_eei_coords_df.iterrows():
        if row['exon1'] in ortholog_id_to_coord:
            ortholog_eei_coords_df.at[i, 'exon1'] = ortholog_id_to_coord[row['exon1']]
        if row['exon2'] in ortholog_id_to_coord:
            ortholog_eei_coords_df.at[i, 'exon2'] = ortholog_id_to_coord[row['exon2']]
    
    # Similarly for human EEIs if present
    human_eei_coords_df = None
    if human_eei_df is not None:
        human_eei_coords_df = human_eei_df.copy()
        for i, row in human_eei_coords_df.iterrows():
            if row['exon1'] in human_id_to_coord:
                human_eei_coords_df.at[i, 'exon1'] = human_id_to_coord[row['exon1']]
            if row['exon2'] in human_id_to_coord:
                human_eei_coords_df.at[i, 'exon2'] = human_id_to_coord[row['exon2']]
    
    # Predict novel human EEIs
    predicted_eeis = predict_novel_human_eeis_from_ortholog(
        ortholog_eei_df=ortholog_eei_coords_df,
        egio_df=egio_df,
        human_eei_df=human_eei_coords_df,
        species_code=species_code
    )
    
    # Convert human coordinates back to IDs for output if possible
    if not predicted_eeis.empty:
        for i, row in predicted_eeis.iterrows():
            if row['exon1'] in human_coord_to_id:
                predicted_eeis.at[i, 'exon1'] = human_coord_to_id[row['exon1']]
            if row['exon2'] in human_coord_to_id:
                predicted_eeis.at[i, 'exon2'] = human_coord_to_id[row['exon2']]
    
    if not os.path.exists(args.output):
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
    # Save results
    predicted_eeis.to_csv(args.output, sep='\t', index=False)
    print(f"Predicted {len(predicted_eeis)} human EEIs, saved to {args.output}")

if __name__ == "__main__":
    main()