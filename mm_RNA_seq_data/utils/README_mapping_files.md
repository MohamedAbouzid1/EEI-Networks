# Coordinate-to-Gene Mapping Files

This script has been modified to create persistent mapping files that can be used later in the main pipeline, improving performance and allowing for offline processing.

## Overview

The script now supports creating mapping files that store coordinate-to-gene relationships, which can be reused across multiple pipeline runs without needing to re-process the GTF annotation files each time.

## Key Features

1. **Create Mapping Files**: Generate comprehensive mapping files from GTF annotation and expression data
2. **Fast Lookup**: Use pre-created mapping files for fast coordinate-to-gene mapping
3. **Multiple Formats**: Support for JSON mapping files that are easy to read and modify
4. **Pipeline Integration**: Easy integration with existing pipeline workflows

## Usage

### Step 1: Create Mapping Files

First, create the mapping files using your GTF annotation and expression data:

```bash
python coord_gene_mapping.py \
    --gtf path/to/your/annotation.gtf \
    --expression path/to/your/expression_data.csv \
    --output-dir mapping_files \
    --species mouse
```

This will create the following files in the `mapping_files` directory:
- `coordinate_gene_mapping.json`: Coordinate-to-gene mappings from GTF
- `expression_gene_mapping.json`: Gene symbol mappings from expression data
- `combined_mapping.json`: Combined mappings for fast lookup
- `mapping_file_paths.json`: Paths to all mapping files

### Step 2: Use Mapping Files in Pipeline

In your main pipeline, use the mapping files for fast coordinate-to-gene conversion:

```python
from coord_gene_mapping import coordinate_to_gene_mapping_with_files

# Map a coordinate to a gene
coordinate = "chr4:148568767:148568879:1"
mapping_dir = "mapping_files"
mapped_gene = coordinate_to_gene_mapping_with_files(coordinate, mapping_dir)

if mapped_gene:
    print(f"Coordinate {coordinate} maps to gene: {mapped_gene}")
else:
    print(f"Could not map coordinate {coordinate}")
```

### Step 3: Process Multiple Coordinates

For processing multiple coordinates efficiently:

```python
from use_mapping_files import process_coordinates_with_mapping_files

coordinates_list = [
    "chr4:148568767:148568879:1",
    "chr1:1000000:1000100:1",
    "chr2:2000000:2000100:1"
]

results = process_coordinates_with_mapping_files(
    coordinates_list, 
    "mapping_files", 
    "path/to/expression_data.csv"
)

for coordinate, result in results.items():
    print(f"{coordinate} -> {result['mapped_gene']} ({result['status']})")
```

## File Structure

The mapping files are stored in JSON format for easy reading and modification:

### coordinate_gene_mapping.json
```json
{
  "chr1:1000000-1000100": {
    "gene_name": "GeneA",
    "gene_id": "ENSG00000123456",
    "chrom": "chr1",
    "start": 1000000,
    "end": 1000100
  }
}
```

### expression_gene_mapping.json
```json
{
  "GENE1": "Gene1",
  "gene1": "Gene1",
  "GENE2": "Gene2",
  "gene2": "Gene2"
}
```

### combined_mapping.json
```json
{
  "coordinate_to_gene": { ... },
  "gene_symbols": { ... },
  "available_genes": ["Gene1", "Gene2", ...]
}
```

## Performance Benefits

1. **Faster Processing**: No need to re-parse GTF files for each coordinate
2. **Offline Capability**: Mapping files can be used without internet connection
3. **Reduced API Calls**: No need for online biomart queries
4. **Batch Processing**: Efficient processing of large coordinate lists

## Integration with Existing Pipeline

To integrate with your existing pipeline, replace coordinate mapping calls with:

```python
# Old way (slow, requires GTF parsing each time)
mapped_gene = coordinate_to_gene_mapping(coordinate, available_genes, gene_intervals)

# New way (fast, uses pre-created mapping files)
mapped_gene = coordinate_to_gene_mapping_with_files(coordinate, mapping_files_dir)
```

## Error Handling

The script includes comprehensive error handling:
- Missing mapping files
- Invalid coordinate formats
- Gene not found in expression data
- Mapping failures

## Example Scripts

- `coord_gene_mapping.py`: Main script for creating mapping files
- `use_mapping_files.py`: Example script showing how to use mapping files in pipeline

## Dependencies

Required Python packages:
- pandas
- pybedtools
- requests (for biomart queries, optional)
- json (built-in)
- os (built-in)

## Notes

- Mapping files are created once and can be reused across multiple pipeline runs
- Update mapping files when you get new annotation or expression data
- The script supports multiple coordinate formats and species
- Mapping files are human-readable JSON format for easy inspection and modification 