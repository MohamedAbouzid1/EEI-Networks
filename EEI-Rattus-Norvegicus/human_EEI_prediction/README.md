# Species-Agnostic EEI Prediction Pipeline

This pipeline predicts novel human Exon-Exon Interactions (EEIs) based on orthologous species EEI data and EGIO orthology mapping. The pipeline has been updated to work with different species, not just mouse.

## Features

- **Species-agnostic**: Works with Drosophila, Mouse, Rat, and other species
- **Automatic species detection**: Detects species from EGIO filename
- **Configurable column names**: Automatically adapts to different EGIO file formats
- **Flexible input formats**: Supports PISA, EPPIC, and other EEI network formats

## Supported Species

Currently configured for:
- **dme** - Drosophila melanogaster (default)
- **mus** - Mus musculus (Mouse)
- **rno** - Rattus norvegicus (Rat)

## Usage

### Basic Command

```bash
python main.py \
  --egio ../EGIO/ExonGroup_testpro_hsa_dme.txt \
  --eei ../data/PISA_networks_filtered/PISA_EEIN_0.5.txt \
  --human_exon_map ../../orthology_based_EEI_prediction/results_CONTACT_BASED/exon_coord_mapping_human.tsv \
  --ortholog_exon_map ../exon_coord_map/exon_coord_map_dse.tsv \
  --human_eei ../../EEI-Homo-Sapiens-2024/results/2-PISA_results/PISA_EEIN_0.5_human.txt \
  --identity_threshold 0.0 \
  --output ../data/predicted/predicted_pisa_human_eeis.tsv
```

### Arguments

- `--egio`: EGIO output file with orthologous exon mappings
- `--eei`: Orthologous species EEI network file
- `--human_exon_map`: Human exon ID to coordinate mapping file
- `--ortholog_exon_map`: Orthologous species exon ID to coordinate mapping file
- `--human_eei`: Optional existing human EEI file to filter out known interactions
- `--identity_threshold`: Minimum sequence identity to consider orthology (default: 0.0)
- `--output`: Output file for predicted EEIs

## Adding New Species

To add support for a new species:

1. **Update `species_config.py`**:
   ```python
   NEW_SPECIES_CONFIG = {
       'species_code': 'new',
       'species_name': 'New Species Name',
       'egio_columns': ['Group', 'hsaEnsemblG', 'hsaPos', 'newEnsemblG', 'newPos', 'Iden', 'Type'],
       'ortholog_pos_col': 'newPos',
       'output_columns': {
           'ortholog_exon1': 'new_species_exon1',
           'ortholog_exon2': 'new_species_exon2'
       }
   }
   ```

2. **Add to the configuration dictionary** in `get_species_config()` function

3. **Update the species detection logic** in `main.py` if needed

## File Formats

### EGIO File Format
Expected columns (order matters):
- `Group`: Orthology group identifier
- `hsaEnsemblG`: Human Ensembl gene ID
- `hsaPos`: Human genomic position
- `[species]EnsemblG`: Species-specific Ensembl gene ID
- `[species]Pos`: Species-specific genomic position
- `Iden`: Sequence identity score
- `Type`: Orthology type (1-1, 1-many, etc.)

### EEI Network File Format
Expected columns:
- `exon1`: First exon identifier or coordinate
- `exon2`: Second exon identifier or coordinate
- Additional columns are preserved in the output

### Exon Mapping File Format
Expected columns:
- `exon_id`: Exon identifier
- `coord`: Genomic coordinate

## Output Format

The pipeline outputs a TSV file with the following columns:
- `exon1`, `exon2`: Human exon identifiers
- `[species]_exon1`, `[species]_exon2`: Original orthologous species exon coordinates
- `confidence`: Prediction confidence score (geometric mean of identities)
- `identity1`, `identity2`: Individual identity scores
- `type1`, `type2`: Orthology types
- All additional columns from the input EEI file

## Examples

### Drosophila Example
```bash
# Run with Drosophila data
python main.py \
  --egio ../EGIO/ExonGroup_testpro_hsa_dme.txt \
  --eei ../data/PISA_networks_filtered/PISA_EEIN_0.5.txt \
  --human_exon_map ../../orthology_based_EEI_prediction/results_CONTACT_BASED/exon_coord_mapping_human.tsv \
  --ortholog_exon_map ../exon_coord_map/exon_coord_map_dse.tsv \
  --output ../data/predicted/predicted_pisa_human_eeis.tsv
```

### Mouse Example
```bash
# Run with Mouse data
python main.py \
  --egio ../EGIO/ExonGroup_testpro_hsa_mus.txt \
  --eei ../data/mouse_EEI_network.txt \
  --human_exon_map ../../orthology_based_EEI_prediction/results_CONTACT_BASED/exon_coord_mapping_human.tsv \
  --ortholog_exon_map ../exon_coord_map/exon_coord_map_mouse.tsv \
  --output ../data/predicted/predicted_mouse_human_eeis.tsv
```

## Troubleshooting

### DtypeWarning
If you see warnings about mixed data types, the pipeline automatically handles this with `low_memory=False`. These warnings are not critical.

### Column Mismatch Errors
Ensure your EGIO file has the expected column order. The pipeline automatically detects species from the filename and uses appropriate column names.

### Memory Issues
For very large files, consider splitting the data or using a machine with more RAM.

## Dependencies

- Python 3.6+
- pandas
- numpy (optional, for advanced operations)

## License

This pipeline is part of the EEI-Conservation project.
