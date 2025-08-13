#!/usr/bin/env python3
"""
Script to demonstrate how to use the mapping files in the main pipeline.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from coord_gene_mapping import (
    coordinate_to_gene_mapping_with_files,
    load_mapping_file,
    load_expression_data,
    load_survival_data
)
import pandas as pd

def process_coordinates_with_mapping_files(coordinates_list, mapping_files_dir, expression_file):
    """
    Process a list of coordinates using pre-created mapping files.
    
    Parameters:
    -----------
    coordinates_list : list
        List of genomic coordinates to process
    mapping_files_dir : str
        Directory containing mapping files
    expression_file : str
        Path to expression data file
        
    Returns:
    --------
    dict: Mapping of coordinates to genes and expression data
    """
    print(f"Processing {len(coordinates_list)} coordinates using mapping files...")
    
    # Load expression data
    expression_df = load_expression_data(expression_file)
    
    results = {}
    
    for coordinate in coordinates_list:
        # Map coordinate to gene
        mapped_gene = coordinate_to_gene_mapping_with_files(coordinate, mapping_files_dir)
        
        if mapped_gene:
            # Get expression data for the mapped gene
            if mapped_gene in expression_df.columns:
                expression_data = expression_df[mapped_gene]
                results[coordinate] = {
                    'mapped_gene': mapped_gene,
                    'expression_data': expression_data,
                    'status': 'success'
                }
            else:
                results[coordinate] = {
                    'mapped_gene': mapped_gene,
                    'expression_data': None,
                    'status': 'gene_not_in_expression'
                }
        else:
            results[coordinate] = {
                'mapped_gene': None,
                'expression_data': None,
                'status': 'mapping_failed'
            }
    
    return results

def create_coordinate_expression_matrix(coordinates_list, mapping_files_dir, expression_file, output_file):
    """
    Create a matrix of coordinate expression data using mapping files.
    
    Parameters:
    -----------
    coordinates_list : list
        List of genomic coordinates
    mapping_files_dir : str
        Directory containing mapping files
    expression_file : str
        Path to expression data file
    output_file : str
        Path to save the expression matrix
    """
    print("Creating coordinate expression matrix...")
    
    # Process coordinates
    results = process_coordinates_with_mapping_files(coordinates_list, mapping_files_dir, expression_file)
    
    # Load expression data
    expression_df = load_expression_data(expression_file)
    
    # Create results matrix
    successful_mappings = {}
    failed_coordinates = []
    
    for coordinate, result in results.items():
        if result['status'] == 'success':
            successful_mappings[coordinate] = result['expression_data']
        else:
            failed_coordinates.append(coordinate)
    
    # Create matrix
    if successful_mappings:
        matrix_df = pd.DataFrame(successful_mappings)
        matrix_df.to_csv(output_file, sep='\t')
        print(f"Expression matrix saved to: {output_file}")
        print(f"Successfully mapped {len(successful_mappings)} coordinates")
        print(f"Failed to map {len(failed_coordinates)} coordinates")
    else:
        print("No successful mappings found!")
    
    return successful_mappings, failed_coordinates

def main():
    """
    Example usage of mapping files in pipeline.
    """
    # Example coordinates (replace with your actual coordinates)
    example_coordinates = [
        "chr4:148568767:148568879:1",
        "chr1:1000000:1000100:1",
        "chr2:2000000:2000100:1"
    ]
    
    # File paths (replace with your actual paths)
    mapping_files_dir = "mapping_files"
    expression_file = "path/to/your/expression_data.csv"
    output_file = "coordinate_expression_matrix.tsv"
    
    print("Example: Using mapping files in pipeline")
    print("=" * 50)
    
    # Check if mapping files exist
    combined_mapping_file = os.path.join(mapping_files_dir, "combined_mapping.json")
    if not os.path.exists(combined_mapping_file):
        print(f"Error: Mapping files not found in {mapping_files_dir}")
        print("Please run the mapping creation script first:")
        print("python coord_gene_mapping.py --gtf <gtf_file> --expression <expression_file> --output-dir mapping_files")
        return
    
    # Process coordinates
    results = process_coordinates_with_mapping_files(example_coordinates, mapping_files_dir, expression_file)
    
    # Print results
    print("\nMapping results:")
    for coordinate, result in results.items():
        status = result['status']
        mapped_gene = result['mapped_gene']
        print(f"  {coordinate} -> {mapped_gene} ({status})")
    
    # Create expression matrix
    print("\nCreating expression matrix...")
    successful, failed = create_coordinate_expression_matrix(
        example_coordinates, 
        mapping_files_dir, 
        expression_file, 
        output_file
    )

if __name__ == "__main__":
    main() 