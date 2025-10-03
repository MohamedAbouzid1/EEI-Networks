# EEI Survival Analysis Pipeline

## Overview

This pipeline validates Exon-Exon Interactions (EEIs) using treatment response patterns as a survival proxy. It analyzes how EEI expression correlates with treatment outcomes in mouse cancer models, specifically focusing on the differential response between E0771 (good responder) and 4T1 (moderate responder) models treated with cyclophosphamide.

## Background

Based on cancer immunotherapy research, this pipeline implements a survival-like analysis using:
- **E0771 + Cyclophosphamide**: Shows robust tumor regression (good prognosis)
- **4T1 + Cyclophosphamide**: Shows growth stasis (moderate prognosis)  
- **Vehicle controls**: Show continued growth (poor prognosis)

The pipeline leverages early IFN response patterns that predict later treatment efficacy.

## Pipeline Components

### 1. Data Loading (`data_loader.py`)
- **Expression data**: Exon-level CPM (Counts Per Million) expression matrix
- **Metadata**: Sample information including model, treatment, timepoint
- **Exon mappings**: Coordinate-to-exon ID mappings for expression lookup

### 2. Response Group Assignment (`analysis.py`)
Creates survival proxy groups based on experimental conditions:

| Model | Treatment | Response Group | Survival Score |
|-------|-----------|----------------|----------------|
| E0771 | Cyclophosphamide | good_responder | 1.0 |
| 4T1 | Cyclophosphamide | moderate_responder | 0.5 |
| Any | Vehicle | poor_control | 0.0 |

### 3. EEI Analysis (`analysis.py`)

#### Core Analysis Functions:

**`analyze_eei_ifn_correlation()`**
- Calculates EEI expression as average of two exon expressions
- Compares expression between response groups using t-tests
- Computes response scores (good_responder_mean - poor_control_mean)
- Applies Benjamini-Hochberg FDR correction for multiple testing

**`identify_time_dependent_eeis()`**
- Analyzes temporal changes in EEI presence
- Compares early (day1-3) vs late (day6+) timepoints
- Identifies EEIs with significant temporal changes (>30% threshold)

**`create_pseudo_survival_analysis()`**
- Groups EEIs by significance (FDR < 0.05)
- Identifies top positive/negative response-associated EEIs
- Creates risk stratification framework

### 4. Visualization (`generate_plots.py`)
Generates four validation plots:
1. **Response Score Distribution**: Histogram of EEI response scores
2. **P-value Distribution**: FDR-adjusted significance distribution
3. **Model Comparison**: E0771 vs 4T1 expression scatter plot
4. **Top Significant EEIs**: Bar chart of most significant associations

### 5. Main Pipeline (`main.py`)
Orchestrates the complete analysis workflow and generates summary statistics.

## Input Requirements

### Required Files:
1. **Mapped EEIs file** (TSV format):
   ```
   mouse_exon1	mouse_exon2
   ENSMUSE00000217453	ENSMUSE00001217049
   ENSMUSE00000217453	ENSMUSE00001238106
   ```

2. **Expression data** (TSV format):
   - Rows: Exon IDs (ENSMUSE format)
   - Columns: Sample IDs
   - Values: CPM expression values

3. **Metadata** (CSV format):
   - Required columns: `model`, `treatment`, `timepoint`, `sample_id`
   - Models: E0771, 4T1
   - Treatments: Cyclophosphamide variants, Vehicle

4. **Exon mapping files** (directory):
   - Coordinate-to-exon ID mappings
   - Used for expression data lookup

## Output Files

### Primary Results:
- **`eei_response_associations.tsv`**: Complete analysis results
  - Columns: `mouse_exon1`, `mouse_exon2`, `response_score`, `p_value`, `p_adj_fdr_bh`, group means/stds
- **`time_dependent_eeis.tsv`**: Temporal analysis results
  - Columns: `mouse_exon1`, `mouse_exon2`, `early_presence`, `late_presence`, `temporal_change`
- **`eei_response_validation_plots.png`**: Four-panel validation plots

### Summary Statistics:
- Total EEIs analyzed
- EEIs with expression data
- Significant EEIs (FDR < 0.05)
- Time-dependent EEIs count

## Statistical Methods

### Multiple Testing Correction:
- **Method**: Benjamini-Hochberg FDR correction
- **Implementation**: Custom implementation ensuring monotonicity
- **Threshold**: FDR < 0.05 for significance

### Statistical Tests:
- **Primary**: Two-sample t-test comparing good_responder vs poor_control
- **Temporal**: Presence/absence analysis across timepoints
- **Effect Size**: Response score (mean difference between groups)

## Usage

### Basic Usage:
```bash
cd mm_RNA_seq_data/eei_analysis_survival
python main.py
```

### Custom Parameters:
Modify the `main_analysis_pipeline()` call in `main.py`:
```python
results = main_analysis_pipeline(
    mapped_eeis_file="path/to/your/eeis.tsv",
    expression_file="path/to/expression.tsv", 
    metadata_file="path/to/metadata.csv",
    mapping_files_dir="path/to/mappings/",
    output_dir="path/to/output/"
)
```

## Key Features

### Robust Analysis:
- Handles missing expression data gracefully
- Applies appropriate statistical corrections
- Provides comprehensive validation plots

### Flexible Input:
- Works with any EEI network format (converted to required schema)
- Compatible with different expression data formats
- Adaptable to various experimental designs

### Quality Control:
- Reports data loading statistics
- Validates sample group distributions
- Provides detailed analysis summaries

## Interpretation Guidelines

### Response Scores:
- **Positive**: EEI higher in good responders (potential therapeutic targets)
- **Negative**: EEI higher in poor responders (potential resistance markers)
- **Magnitude**: Larger absolute values indicate stronger associations

### Significance:
- **FDR < 0.05**: Statistically significant after multiple testing correction
- **Raw p < 0.05**: Significant before correction (use with caution)

### Temporal Patterns:
- **Gained**: EEI presence increases over time
- **Lost**: EEI presence decreases over time
- **Threshold**: >30% change considered significant

## Troubleshooting

### Common Issues:
1. **No expression data**: Check exon ID format compatibility
2. **Empty results**: Verify sample group assignments
3. **Missing mappings**: Ensure coordinate-to-exon mappings are complete

### Data Requirements:
- Minimum 3 samples per group for statistical tests
- Consistent exon ID formats across files
- Complete metadata for all samples

## Dependencies

- pandas
- numpy  
- scipy
- matplotlib
- seaborn
- sklearn

## Version History

- **v1.0**: Initial implementation with raw p-values
- **v1.1**: Added FDR correction and improved visualization
- **v1.2**: Enhanced documentation and error handling

## Citation

If using this pipeline, please cite the original research that inspired the methodology and acknowledge the EEI conservation analysis framework.
