# EEI Survival Analysis Plotting Script

This script has been updated to handle multiple data results instead of hardcoded paths. It can now process multiple results directories and automatically find data files.

## Usage Options

### 1. List Available Data
First, see what data is available:
```bash
python plot_survival_results.py --list-results
```

This will show you:
- Available results directories
- Available expression files
- Available survival files

### 2. Process Single Results Directory
Process a specific results directory:
```bash
python plot_survival_results.py --results-dir my_results --expression-file my_expression.tsv --survival-file my_survival.csv
```

### 3. Batch Process All Results
Process all results directories found in the current directory:
```bash
python plot_survival_results.py --batch-mode
```

### 4. Custom Output Directory
Specify a custom output directory:
```bash
python plot_survival_results.py --results-dir my_results --output-dir custom_plots
```

### 5. Custom Coordinate Mapping Directory
Specify a custom coordinate mapping directory:
```bash
python plot_survival_results.py --results-dir my_results --coord-mapping-dir my_mapping_files
```

## Command Line Arguments

- `--results-dir, -r`: Directory containing analysis results (default: results)
- `--expression-file, -e`: Expression data file (TSV format)
- `--survival-file, -s`: Survival data file (CSV format)
- `--output-dir, -o`: Output directory for plots (default: results/plots)
- `--coord-mapping-dir, -c`: Directory for coordinate to gene mapping files (default: coord_gene_mapping_files)
- `--list-results`: List available results directories and exit
- `--batch-mode`: Process all results directories found in current directory

## File Structure Expected

Each results directory should contain:
- `all_eei_survival_results.tsv`
- `significant_eei_survival_results.tsv`
- `survival_analysis_summary.json`

## Data File Formats

### Expression File
- TSV format with genes as rows and samples as columns
- First column should be gene identifiers

### Survival File
- CSV format with samples as rows
- Should contain columns for:
  - Survival time (one of: survival_days, survival_time, time, days)
  - Event status (one of: event_status, status, event, dead)

## Examples

### Example 1: Process Default Results
```bash
python plot_survival_results.py
```

### Example 2: Process Specific Results with Custom Files
```bash
python plot_survival_results.py \
  --results-dir experiment_1_results \
  --expression-file data/expression_matrix.tsv \
  --survival-file data/survival_data.csv \
  --output-dir experiment_1_plots
```

### Example 3: Batch Process All Experiments
```bash
python plot_survival_results.py --batch-mode
```

### Example 4: List Available Data
```bash
python plot_survival_results.py --list-results
```

## Output

The script will create:
- P-value distribution plots
- Significance summary plots
- Survival curves (if data files are provided)
- Analysis overview plots
- Text report with summary statistics

All plots are saved in the specified output directory (or `results/plots` by default).
