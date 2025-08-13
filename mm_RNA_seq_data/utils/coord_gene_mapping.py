import pybedtools
from pybedtools import BedTool
import requests
import time
from io import StringIO
import pandas as pd
import json
import os
import pickle

# Option 1: Using GTF/GFF annotation file (Recommended)
def load_gene_annotation(gtf_file_path):
    """
    Load gene annotation from GTF/GFF file and create coordinate->gene mapping.
    
    Parameters:
    -----------
    gtf_file_path : str
        Path to GTF/GFF annotation file (e.g., from GENCODE, Ensembl)
        
    Returns:
    --------
    dict: coordinate->gene mapping and BedTool object for interval queries
    """
    print(f"Loading gene annotation from {gtf_file_path}...")
    
    # Parse GTF file
    annotation_bed = BedTool(gtf_file_path)
    
    # Create coordinate->gene mapping
    coord_to_gene = {}
    gene_intervals = []
    
    for feature in annotation_bed:
        if feature.fields[2] == 'exon':  # Focus on exons
            chrom = feature.chrom
            start = int(feature.start)
            end = int(feature.end)
            
            # Extract gene information from attributes
            attributes = feature.fields[8]
            gene_id = extract_attribute(attributes, 'gene_id')
            gene_name = extract_attribute(attributes, 'gene_name')
            
            if gene_name:
                gene_intervals.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'gene': gene_name,
                    'gene_id': gene_id
                })
    
    print(f"Loaded {len(gene_intervals)} exon intervals")
    return gene_intervals

def extract_attribute(attr_string, attr_name):
    """Extract specific attribute from GTF attribute string."""
    import re
    pattern = f'{attr_name} "([^"]+)"'
    match = re.search(pattern, attr_string)
    return match.group(1) if match else None

def coordinate_to_gene_mapping_gtf(coordinate, gene_intervals, available_genes):
    """
    Map genomic coordinate to gene using GTF annotation.
    
    Parameters:
    -----------
    coordinate : str
        Genomic coordinate in format "chr4:148568767:148568879:1"
    gene_intervals : list
        List of gene interval dictionaries from load_gene_annotation()
    available_genes : list
        Available genes in expression data
        
    Returns:
    --------
    str or None: Mapped gene name
    """
    try:
        # Parse coordinate
        parts = coordinate.replace('chr', '').split(':')
        if len(parts) < 3:
            return None
            
        chrom = f"chr{parts[0]}"
        start = int(parts[1])
        end = int(parts[2])
        
        # Find overlapping genes
        overlapping_genes = []
        for interval in gene_intervals:
            if (interval['chrom'] == chrom and 
                not (end < interval['start'] or start > interval['end'])):
                overlapping_genes.append(interval['gene'])
        
        # Return gene if available in expression data
        for gene in overlapping_genes:
            if gene in available_genes:
                return gene
                
        # Try gene symbols or IDs
        for gene in overlapping_genes:
            # Try exact matches or partial matches
            matches = [g for g in available_genes if gene.upper() in g.upper() or g.upper() in gene.upper()]
            if matches:
                return matches[0]
                
        return None
        
    except Exception as e:
        print(f"Error mapping coordinate {coordinate}: {e}")
        return None

# Option 2: Using biomaRt/Ensembl API (for online queries)
def coordinate_to_gene_biomart(coordinate, species="mmusculus"):
    """
    Map coordinate to gene using Ensembl REST API.
    
    Parameters:
    -----------
    coordinate : str
        Genomic coordinate in format "chr4:148568767:148568879:1"
    species : str
        Species (default: mmusculus for mouse)
        
    Returns:
    --------
    list: Gene names overlapping the coordinate
    """
    try:
        # Parse coordinate
        parts = coordinate.replace('chr', '').split(':')
        if len(parts) < 3:
            return []
            
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        
        # Query Ensembl REST API
        url = f"https://rest.ensembl.org/overlap/region/{species}/{chrom}:{start}-{end}?feature=gene;content-type=application/json"
        
        response = requests.get(url)
        time.sleep(0.1)  # Be nice to the API
        
        if response.status_code == 200:
            genes = response.json()
            gene_names = [gene.get('external_name', gene.get('id', '')) for gene in genes if gene.get('external_name')]
            return gene_names
        else:
            return []
            
    except Exception as e:
        print(f"Error querying biomart for {coordinate}: {e}")
        return []

def coordinate_to_gene_mapping_biomart(coordinate, available_genes, cache={}):
    """
    Map coordinate to gene using biomart with caching.
    """
    if coordinate in cache:
        return cache[coordinate]
    
    gene_candidates = coordinate_to_gene_biomart(coordinate)
    
    # Find best match in available genes
    for gene in gene_candidates:
        if gene in available_genes:
            cache[coordinate] = gene
            return gene
            
        # Try partial matches
        matches = [g for g in available_genes if gene.upper() in g.upper() or g.upper() in gene.upper()]
        if matches:
            cache[coordinate] = matches[0]
            return matches[0]
    
    cache[coordinate] = None
    return None

# Option 3: Pre-built coordinate mapping using gene symbols
def create_gene_coordinate_database(expression_genes, species="mouse"):
    """
    Create a mapping database for genes in your expression data.
    This uses gene symbols to approximate genomic locations.
    """
    gene_mapping = {}
    
    # This is a simplified approach - you'd want to use proper annotation
    # For each gene in expression data, you could query its coordinates
    for gene in expression_genes:
        # Store gene symbol for lookup
        gene_mapping[gene.upper()] = gene
    
    return gene_mapping

def coordinate_to_gene_mapping_symbol(coordinate, gene_mapping, available_genes):
    """
    Map coordinate to gene using gene symbol approximation.
    This is the most basic approach when other methods aren't available.
    """
    try:
        # Extract potential gene information from coordinate
        # This is very basic and may not work for all coordinate formats
        
        # Method 1: If coordinate contains gene-like information
        coord_parts = coordinate.replace('chr', '').split(':')
        
        # Look for gene symbols in available genes that might match
        # This is a fallback method and quite limited
        
        for gene in available_genes:
            # Try various matching strategies
            if any(part.upper() in gene.upper() for part in coord_parts if len(part) > 2):
                return gene
        
        return None
        
    except:
        return None

# Updated main coordinate mapping function
def coordinate_to_gene_mapping(coordinate, available_genes, gene_intervals=None, use_biomart=False):
    """
    Enhanced coordinate to gene mapping with multiple strategies.
    
    Parameters:
    -----------
    coordinate : str
        Genomic coordinate
    available_genes : list
        Available genes in expression data
    gene_intervals : list, optional
        Gene intervals from GTF annotation
    use_biomart : bool
        Whether to use biomart API queries
        
    Returns:
    --------
    str or None: Mapped gene name
    """
    
    # Strategy 1: Direct match (if coordinates are actually gene names)
    if coordinate in available_genes:
        return coordinate
    
    # Strategy 2: Use GTF annotation if available
    if gene_intervals:
        result = coordinate_to_gene_mapping_gtf(coordinate, gene_intervals, available_genes)
        if result:
            return result
    
    # Strategy 3: Use biomart if enabled
    if use_biomart:
        result = coordinate_to_gene_mapping_biomart(coordinate, available_genes)
        if result:
            return result
    
    # Strategy 4: Symbol-based approximation (fallback)
    gene_mapping = create_gene_coordinate_database(available_genes)
    result = coordinate_to_gene_mapping_symbol(coordinate, gene_mapping, available_genes)
    if result:
        return result
    
    # Strategy 5: Try chromosome-based matching (your original approach)
    try:
        chrom = coordinate.split(':')[0].replace('chr', '')
        for gene in available_genes:
            if chrom.lower() in gene.lower():
                return gene
    except:
        pass
    
    return None

# Enhanced find_expression_for_coordinate function
def find_expression_for_coordinate(coordinate, expression_df, gene_intervals=None):
    """
    Enhanced function to find expression data for a genomic coordinate.
    """
    
    # Method 1: Direct coordinate match
    if coordinate in expression_df.columns:
        return expression_df[coordinate]
    
    # Method 2: Use coordinate to gene mapping
    mapped_gene = coordinate_to_gene_mapping(
        coordinate, 
        expression_df.columns.tolist(), 
        gene_intervals=gene_intervals,
        use_biomart=False  # Set to True if you want to use API queries
    )
    
    if mapped_gene and mapped_gene in expression_df.columns:
        return expression_df[mapped_gene]
    
    return None

# Usage example for your main function
def load_annotation_and_test(mapped_eeis_file, expression_file, survival_file, 
                           gtf_file=None, **kwargs):
    """
    Enhanced version of your test function with proper coordinate mapping.
    
    Parameters:
    -----------
    gtf_file : str, optional
        Path to GTF annotation file for coordinate mapping
    """
    
    # Load gene annotation if provided
    gene_intervals = None
    if gtf_file:
        gene_intervals = load_gene_annotation(gtf_file)
    
    # Load your data as before
    mapped_eeis = pd.read_csv(mapped_eeis_file, sep='\t')
    expression_df = load_expression_data(expression_file)
    survival_df = load_survival_data(survival_file)
    
    # Now your coordinate mapping will work much better
    # ... rest of your analysis code
    
    return test_orthologous_eeis_survival(
        mapped_eeis_file, expression_file, survival_file,
        gene_intervals=gene_intervals, **kwargs
    )

def create_coordinate_gene_mapping_file(gtf_file_path, output_file, species="mouse"):
    """
    Create a comprehensive coordinate-to-gene mapping file from GTF annotation.
    
    Parameters:
    -----------
    gtf_file_path : str
        Path to GTF/GFF annotation file
    output_file : str
        Path to save the mapping file (JSON format)
    species : str
        Species name for biomart queries if needed
        
    Returns:
    --------
    dict: The created mapping dictionary
    """
    print(f"Creating coordinate-to-gene mapping from {gtf_file_path}...")
    
    # Load gene annotation
    gene_intervals = load_gene_annotation(gtf_file_path)
    
    # Create comprehensive mapping
    coord_to_gene_mapping = {}
    
    for interval in gene_intervals:
        coord_key = f"{interval['chrom']}:{interval['start']}-{interval['end']}"
        
        if coord_key not in coord_to_gene_mapping:
            coord_to_gene_mapping[coord_key] = {
                'gene_name': interval['gene'],
                'gene_id': interval['gene_id'],
                'chrom': interval['chrom'],
                'start': interval['start'],
                'end': interval['end']
            }
        else:
            # If multiple intervals for same gene, keep the one with longer overlap
            current_length = coord_to_gene_mapping[coord_key]['end'] - coord_to_gene_mapping[coord_key]['start']
            new_length = interval['end'] - interval['start']
            if new_length > current_length:
                coord_to_gene_mapping[coord_key] = {
                    'gene_name': interval['gene'],
                    'gene_id': interval['gene_id'],
                    'chrom': interval['chrom'],
                    'start': interval['start'],
                    'end': interval['end']
                }
    
    # Save mapping to file
    with open(output_file, 'w') as f:
        json.dump(coord_to_gene_mapping, f, indent=2)
    
    print(f"Created mapping file with {len(coord_to_gene_mapping)} coordinate-gene pairs")
    print(f"Mapping saved to: {output_file}")
    
    return coord_to_gene_mapping

def create_expression_gene_mapping_file(expression_file, output_file):
    """
    Create a mapping file for genes in expression data.
    
    Parameters:
    -----------
    expression_file : str
        Path to expression data file
    output_file : str
        Path to save the gene mapping file
        
    Returns:
    --------
    dict: The created gene mapping dictionary
    """
    print(f"Creating gene mapping from expression file: {expression_file}")
    
    # Load expression data to get gene names
    expression_df = load_expression_data(expression_file)
    
    # For transposed files, gene names are in the index (rows), not columns
    # Check if this is a transposed file by looking at the first few column names
    first_cols = expression_df.columns.tolist()[:5]
    if any('mouse' in str(col).lower() or 'BM' in str(col) for col in first_cols):
        # This is a transposed file - genes are in the index
        available_genes = expression_df.index.tolist()
        print(f"Detected transposed file - using {len(available_genes)} genes from index")
    else:
        # This is a regular file - genes are in the columns
        available_genes = expression_df.columns.tolist()
        print(f"Detected regular file - using {len(available_genes)} genes from columns")
    
    # Create gene mapping dictionary
    gene_mapping = {}
    for gene in available_genes:
        gene_mapping[gene.upper()] = gene
        gene_mapping[gene.lower()] = gene
        gene_mapping[gene] = gene
    
    # Save mapping to file
    with open(output_file, 'w') as f:
        json.dump(gene_mapping, f, indent=2)
    
    print(f"Created gene mapping file with {len(available_genes)} genes")
    print(f"Mapping saved to: {output_file}")
    
    return gene_mapping

def load_mapping_file(mapping_file_path):
    """
    Load a mapping file created by the above functions.
    
    Parameters:
    -----------
    mapping_file_path : str
        Path to the mapping file
        
    Returns:
    --------
    dict: The loaded mapping dictionary
    """
    with open(mapping_file_path, 'r') as f:
        return json.load(f)

def create_comprehensive_mapping_files(gtf_file_path, expression_file, output_dir, species="mouse"):
    """
    Create all necessary mapping files for the pipeline.
    
    Parameters:
    -----------
    gtf_file_path : str
        Path to GTF annotation file
    expression_file : str
        Path to expression data file
    output_dir : str
        Directory to save mapping files
    species : str
        Species name
        
    Returns:
    --------
    dict: Paths to created mapping files
    """
    os.makedirs(output_dir, exist_ok=True)
    
    mapping_files = {}
    
    # Create coordinate-to-gene mapping
    coord_mapping_file = os.path.join(output_dir, "coordinate_gene_mapping.json")
    coord_mapping = create_coordinate_gene_mapping_file(gtf_file_path, coord_mapping_file, species)
    mapping_files['coordinate_gene_mapping'] = coord_mapping_file
    
    # Create expression gene mapping
    gene_mapping_file = os.path.join(output_dir, "expression_gene_mapping.json")
    gene_mapping = create_expression_gene_mapping_file(expression_file, gene_mapping_file)
    mapping_files['expression_gene_mapping'] = gene_mapping_file
    
    # Create combined mapping for quick lookup
    combined_mapping_file = os.path.join(output_dir, "combined_mapping.json")
    combined_mapping = {
        'coordinate_to_gene': coord_mapping,
        'gene_symbols': gene_mapping,
        'available_genes': list(gene_mapping.values())
    }
    
    with open(combined_mapping_file, 'w') as f:
        json.dump(combined_mapping, f, indent=2)
    mapping_files['combined_mapping'] = combined_mapping_file
    
    # Save mapping file paths
    paths_file = os.path.join(output_dir, "mapping_file_paths.json")
    with open(paths_file, 'w') as f:
        json.dump(mapping_files, f, indent=2)
    
    print(f"All mapping files created in: {output_dir}")
    print(f"Mapping file paths saved to: {paths_file}")
    
    return mapping_files

def coordinate_to_gene_mapping_with_files(coordinate, mapping_files_dir):
    """
    Map coordinate to gene using pre-created mapping files.
    
    Parameters:
    -----------
    coordinate : str
        Genomic coordinate (can be in EEI format: "chr4:148568767:148568879:1")
    mapping_files_dir : str
        Directory containing mapping files
        
    Returns:
    --------
    str or None: Mapped gene name
    """
    # Handle NaN/empty values
    if pd.isna(coordinate) or coordinate == '' or coordinate is None:
        return None
    
    # Convert to string if it's not already
    coordinate = str(coordinate)
    
    # Load mapping files
    combined_mapping_file = os.path.join(mapping_files_dir, "combined_mapping.json")
    if not os.path.exists(combined_mapping_file):
        raise FileNotFoundError(f"Mapping file not found: {combined_mapping_file}")
    
    with open(combined_mapping_file, 'r') as f:
        mappings = json.load(f)
    
    coord_mapping = mappings['coordinate_to_gene']
    gene_symbols = mappings['gene_symbols']
    available_genes = mappings['available_genes']
    
    # Try direct coordinate match (for mapping file format)
    if coordinate in coord_mapping:
        gene_name = coord_mapping[coordinate]['gene_name']
        if gene_name in available_genes:
            return gene_name
    
    # Try parsing EEI coordinate format: "chr4:148568767:148568879:1"
    try:
        # Handle EEI format: "chr4:148568767:148568879:1"
        if coordinate.startswith('chr') and coordinate.count(':') >= 3:
            parts = coordinate.split(':')
            if len(parts) >= 4:
                chrom = parts[0]  # "chr4"
                start = int(parts[1])  # 148568767
                end = int(parts[2])    # 148568879
                strand = parts[3]       # "1"
                
                # Convert to mapping file format: "4:148568767-148568879"
                chrom_num = chrom.replace('chr', '')
                mapping_key = f"{chrom_num}:{start}-{end}"
                
                if mapping_key in coord_mapping:
                    gene_name = coord_mapping[mapping_key]['gene_name']
                    if gene_name in available_genes:
                        return gene_name
                
                # Also try without strand info: "chr4:148568767:148568879"
                alt_coord = f"{chrom}:{start}:{end}"
                if alt_coord in coord_mapping:
                    gene_name = coord_mapping[alt_coord]['gene_name']
                    if gene_name in available_genes:
                        return gene_name
                
                # Try to find overlapping coordinates
                for coord_key, gene_info in coord_mapping.items():
                    if gene_info['chrom'] == chrom_num:
                        coord_start = gene_info['start']
                        coord_end = gene_info['end']
                        
                        # Check for overlap
                        if not (end < coord_start or start > coord_end):
                            gene_name = gene_info['gene_name']
                            if gene_name in available_genes:
                                return gene_name
    except (ValueError, IndexError):
        pass
    
    # Try parsing coordinate and finding overlaps (original logic)
    try:
        parts = coordinate.replace('chr', '').split(':')
        if len(parts) >= 3:
            chrom = f"chr{parts[0]}"
            start = int(parts[1])
            end = int(parts[2])
            
            # Find overlapping coordinates
            for coord_key, gene_info in coord_mapping.items():
                if gene_info['chrom'] == chrom:
                    coord_start = gene_info['start']
                    coord_end = gene_info['end']
                    
                    # Check for overlap
                    if not (end < coord_start or start > coord_end):
                        gene_name = gene_info['gene_name']
                        if gene_name in available_genes:
                            return gene_name
    except:
        pass
    
    # Try gene symbol matching with proper type checking
    for gene in available_genes:
        try:
            if isinstance(gene, str) and gene:
                if gene.upper() in coordinate.upper() or coordinate.upper() in gene.upper():
                    return gene
        except (AttributeError, TypeError):
            continue
    
    return None

def load_expression_data(expression_file):
    """
    Load expression data from file.
    """
    if expression_file.endswith('.csv'):
        return pd.read_csv(expression_file, index_col=0)
    elif expression_file.endswith('.tsv') or expression_file.endswith('.txt'):
        return pd.read_csv(expression_file, sep='\t', index_col=0)
    else:
        return pd.read_csv(expression_file, index_col=0)

def load_survival_data(survival_file):
    """
    Load survival data from file.
    """
    if survival_file.endswith('.csv'):
        return pd.read_csv(survival_file, index_col=0)
    elif survival_file.endswith('.tsv') or survival_file.endswith('.txt'):
        return pd.read_csv(survival_file, sep='\t', index_col=0)
    else:
        return pd.read_csv(survival_file, index_col=0)

def main():
    """
    Main function to create mapping files for the pipeline.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Create coordinate-to-gene mapping files for the pipeline')
    parser.add_argument('--gtf', required=True, help='Path to GTF annotation file')
    parser.add_argument('--expression', required=True, help='Path to expression data file')
    parser.add_argument('--output-dir', required=True, help='Output directory for mapping files')
    parser.add_argument('--species', default='mouse', help='Species name (default: mouse)')
    
    args = parser.parse_args()
    
    # Create mapping files
    mapping_files = create_comprehensive_mapping_files(
        args.gtf,
        args.expression,
        args.output_dir,
        args.species
    )
    
    print("\nMapping files created successfully!")
    print("Files created:")
    for file_type, file_path in mapping_files.items():
        print(f"  {file_type}: {file_path}")
    
    print(f"\nTo use these mapping files in your pipeline:")
    print(f"1. Use coordinate_to_gene_mapping_with_files(coordinate, '{args.output_dir}')")
    print(f"2. Or load specific mapping files using load_mapping_file()")

def example_usage():
    """
    Example usage of the mapping functions.
    """
    # Example 1: Create mapping files
    print("Example 1: Creating mapping files")
    print("=" * 50)
    
    # This would be your actual file paths
    gtf_file = "path/to/your/annotation.gtf"
    expression_file = "path/to/your/expression_data.csv"
    output_dir = "mapping_files"
    
    # Uncomment to run:
    # mapping_files = create_comprehensive_mapping_files(gtf_file, expression_file, output_dir)
    
    print("\nExample 2: Using mapping files in pipeline")
    print("=" * 50)
    
    # Example coordinate
    coordinate = "chr4:148568767:148568879:1"
    mapping_dir = "mapping_files"
    
    # Uncomment to run:
    # mapped_gene = coordinate_to_gene_mapping_with_files(coordinate, mapping_dir)
    # print(f"Coordinate {coordinate} maps to gene: {mapped_gene}")
    
    print("\nExample 3: Loading specific mapping files")
    print("=" * 50)
    
    # Uncomment to run:
    # coord_mapping = load_mapping_file("mapping_files/coordinate_gene_mapping.json")
    # gene_mapping = load_mapping_file("mapping_files/expression_gene_mapping.json")
    # print(f"Loaded {len(coord_mapping)} coordinate mappings")
    # print(f"Loaded {len(gene_mapping)} gene mappings")

if __name__ == "__main__":
    main()