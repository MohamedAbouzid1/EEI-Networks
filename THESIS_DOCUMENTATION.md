# EEI Conservation Project: Comprehensive Documentation

## Project Overview

This repository contains the complete implementation of a three-pipeline approach for constructing and expanding Human Exon-Exon Interaction (EEI) networks through evolutionary conservation analysis. The project aims to identify and predict human EEIs by leveraging structural data from multiple eukaryote species and orthology relationships.

## Research Context

The work builds upon the research presented in:
**"Prognostic importance of splicing-triggered aberrations of protein complex interfaces in cancer"** by *Khalique Newaz, Jan Baumbach and Dmitrij Frishman*

This thesis project extends the original work by implementing a comprehensive conservation-based approach to expand human EEI networks beyond what can be detected directly from human structural data.

## Project Architecture

The project is organized into three main pipelines that work sequentially:

### Pipeline 1: Human EEI Network Construction
**Location**: `EEI-Homo-Sapiens/`

**Objective**: Construct the baseline human EEI network using three complementary detection methods.

**Methods Implemented**:
1. **Contact-based Detection** (`1-Contact-based/`)
   - Identifies EEIs based on physical contacts between exons in protein structures
   - Uses UniProt data and PDB structures
   - Maps exon coordinates to protein structures via SIFTS

2. **PISA-based Detection** (`2-PISA-based/`)
   - Energy-based approach using PISA (Protein Interfaces, Surfaces and Assemblies)
   - Identifies stable protein interfaces that correspond to EEIs
   - Uses thermodynamic analysis of protein-protein interactions

3. **EPPIC-based Detection** (`3-EPPIC-based/`)
   - Evolution-based approach using EPPIC (Evolutionary Protein-Protein Interface Classifier)
   - Leverages evolutionary conservation patterns
   - Identifies interfaces that are evolutionarily conserved

**Output**: High-confidence human EEI network combining all three methods

### Pipeline 2: Species-Specific EEI Network Construction
**Location**: `EEI-[Species-Name]/`

**Objective**: Construct EEI networks for 7 selected eukaryote species using the same three detection methods.

**Species Selection Criteria**:
The 7 eukaryote species were selected based on their abundance of high-resolution PDB structures:

1. **Mus musculus** (Mouse) - `EEI-Mus-Musculus/`
2. **Bos taurus** (Cattle) - `EEI-Bos-Taurus/`
3. **Drosophila melanogaster** (Fruit fly) - `EEI-Drosophila/`
4. **Gallus gallus** (Chicken) - `EEI-Gallus-Gallus/`
5. **Oryctolagus cuniculus** (Rabbit) - `EEI-Oryctolagus/`
6. **Rattus norvegicus** (Rat) - `EEI-Rattus-Norvegicus/`
7. **Saccharomyces cerevisiae** (Yeast) - `EEI-Saccharomyces/`

**Implementation**: Each species folder contains the same three-method pipeline as human:
- `1-Contact-based/` - Contact-based EEI detection
- `2-PISA-based/` - PISA-based EEI detection  
- `3-EPPIC-based/` - EPPIC-based EEI detection

### Pipeline 3: Orthology Mapping and EEI Prediction
**Location**: `orthology_based_EEI_prediction/`

**Objective**: Expand human EEI networks by leveraging orthology relationships between human and other species.

**Components**:

#### 3.1 EGIO (Exon Group Ideogram) Analysis
**Location**: `EGIO/` and species-specific `EGIO/` folders

**Purpose**: Detect orthologous exons between human and each of the 7 species using EGIO methodology.

**EGIO Implementation**:
- Uses dynamic programming for isoform alignment
- Reciprocal BLASTN results guide the alignment process
- Identifies orthologous exon groups with identity and coverage thresholds
- Generates two main outputs:
  - `ExonGroup.txt`: Orthologous exon groups
  - `OrthoIso.txt`: Orthologous isoforms

**Parameters**:
- Identity threshold: 0.8 (default)
- Coverage threshold: 0.8 (default)
- Match score: 2, Mismatch penalty: -2, Gap penalty: -1

#### 3.2 Orthology-Based EEI Prediction
**Location**: `orthology_based_EEI_prediction/`

**Process**:
1. **Data Integration**: Combine EGIO orthology results with species-specific EEI networks
2. **Mapping**: Map orthologous exons between species and human
3. **Prediction**: Predict human EEIs based on conserved interactions in other species
4. **Validation**: Cross-validate predictions against known human EEIs

**Output Analysis**:
- Coverage analysis of predicted vs. known EEIs
- Confidence scoring based on conservation patterns
- Statistical evaluation of prediction accuracy
- Visualization of conservation patterns

## Technical Implementation

### Data Sources
- **UniProt**: Protein sequence and annotation data
- **PDB**: Protein structure database
- **SIFTS**: Structure Integration with Function, Taxonomy and Sequence
- **Ensembl**: Genomic annotations and orthology data
- **PISA**: Protein interface analysis
- **EPPIC**: Evolutionary interface classification

### Software Dependencies
- **R**: Statistical analysis and data processing
- **Python**: Machine learning and data analysis
- **BLAST+**: Sequence alignment for EGIO
- **GCC**: Compilation of EGIO C components
- **Pandas/NumPy**: Data manipulation
- **Plotnine**: Visualization

### Pipeline Execution Order

1. **Human Baseline** (`EEI-Homo-Sapiens/`)
   ```
   1-Contact-based/ → 2-PISA-based/ → 3-EPPIC-based/ → 4_find_global_EEINs.sh
   ```

2. **Species Networks** (for each of 7 species)
   ```
   EGIO/ → 1-Contact-based/ → 2-PISA-based/ → 3-EPPIC-based/
   ```

3. **Orthology Mapping**
   ```
   EGIO results + Species EEI networks → human_EEI_prediction/ → Analysis
   ```

## Results and Outputs

### Human EEI Network
- **Contact-based**: Physical interaction-based EEIs
- **PISA-based**: Energy-stable interface EEIs  
- **EPPIC-based**: Evolutionarily conserved EEIs
- **Combined**: High-confidence network (NETHIGH)

### Species-Specific Networks
- Individual EEI networks for each of 7 species
- Same three-method approach as human
- Species-specific structural data and annotations

### Conservation-Based Predictions
- **Coverage Analysis**: Percentage of human EEIs that can be predicted from other species
- **Confidence Scoring**: Reliability metrics for predicted EEIs
- **Conservation Patterns**: Evolutionary conservation analysis
- **Validation Results**: Accuracy of predictions against known human EEIs

## File Structure

```
EEI-Conservation-main/
├── EEI-Homo-Sapiens/           # Pipeline 1: Human baseline
│   ├── 1-Contact-based/        # Contact-based detection
│   ├── 2-PISA-based/           # PISA-based detection
│   ├── 3-EPPIC-based/          # EPPIC-based detection
│   └── results/                # Combined human EEI network
├── EEI-[Species]/              # Pipeline 2: Species networks (7 species)
│   ├── EGIO/                   # Orthology detection
│   ├── 1-Contact-based/        # Contact-based detection
│   ├── 2-PISA-based/           # PISA-based detection
│   ├── 3-EPPIC-based/          # EPPIC-based detection
│   └── human_EEI_prediction/   # Human EEI predictions
├── EGIO/                       # Global EGIO results
├── orthology_based_EEI_prediction/  # Pipeline 3: Orthology mapping
│   ├── human_EEI_prediction/   # Prediction algorithms
│   └── results_*/              # Analysis results by method
└── data/                       # Shared data resources
```

## Detailed Methodology

### Pipeline 1: Human EEI Network Construction

#### Contact-Based Detection Method
**Input Data**: UniProt sequences, PDB structures, SIFTS mappings
**Process**:
1. Download and preprocess UniProt data for human proteins
2. Map UniProt IDs to PDB structures via SIFTS
3. Download relevant PDB structures
4. Map exon coordinates to protein structures
5. Identify physical contacts between exons in protein complexes
6. Generate contact-based EEI network

**Key Scripts**:
- `1a_download_uniprot.r` - UniProt data download
- `2_download_sifts_mapping.r` - SIFTS mapping
- `3_download_PDBs.r` - PDB structure download
- `8_exonexon_contact.r` - Contact identification
- `9_save_contact_networks.r` - Network generation

#### PISA-Based Detection Method
**Input Data**: PDB structures, PISA interface analysis
**Process**:
1. Run PISA analysis on protein complexes
2. Parse PISA results for interface information
3. Map interfaces to exon coordinates
4. Identify energy-stable interfaces corresponding to EEIs
5. Generate PISA-based EEI network

**Key Scripts**:
- `7_runPISA.r` - PISA analysis execution
- `8_parse_files.r` - Result parsing
- `9a_map_EEIs.r` - EEI mapping
- `9b_preprocess_EEINs.r` - Network preprocessing

#### EPPIC-Based Detection Method
**Input Data**: PDB structures, evolutionary conservation data
**Process**:
1. Run EPPIC analysis on protein interfaces
2. Identify evolutionarily conserved interfaces
3. Map conserved interfaces to exon coordinates
4. Generate EPPIC-based EEI network

**Key Scripts**:
- `1_run_eppic.r` - EPPIC analysis execution
- `2_map_EEI.r` - EEI mapping
- `3_preprocess_EEINs.r` - Network preprocessing

### Pipeline 2: Species-Specific EEI Networks

#### Species Selection Rationale
The 7 eukaryote species were selected based on:
1. **PDB Structure Abundance**: High-resolution PDB structures available
2. **Evolutionary Distance**: Diverse phylogenetic representation
3. **Genomic Completeness**: Well-annotated genomes
4. **Orthology Data**: Available orthology relationships with human

**Species Details**:
- **Mus musculus**: Model organism, extensive structural data
- **Bos taurus**: Mammalian conservation, agricultural importance
- **Drosophila melanogaster**: Invertebrate model, evolutionary distance
- **Gallus gallus**: Avian model, evolutionary perspective
- **Oryctolagus cuniculus**: Mammalian model, research applications
- **Rattus norvegicus**: Model organism, biomedical research
- **Saccharomyces cerevisiae**: Simple eukaryote, fundamental biology

#### Implementation for Each Species
Each species follows the same three-method approach as human:
1. **Data Preparation**: Species-specific UniProt, PDB, and annotation data
2. **Contact-Based**: Physical interaction detection
3. **PISA-Based**: Energy-stable interface detection
4. **EPPIC-Based**: Evolutionarily conserved interface detection

### Pipeline 3: Orthology Mapping and Prediction

#### EGIO (Exon Group Ideogram) Implementation
**Purpose**: Detect orthologous exons between human and each species

**Technical Details**:
- **Algorithm**: Dynamic programming for isoform alignment
- **Sequence Alignment**: Reciprocal BLASTN for exon sequences
- **Parameters**:
  - Identity threshold: 0.8 (80% sequence identity)
  - Coverage threshold: 0.8 (80% sequence coverage)
  - Match score: 2, Mismatch penalty: -2, Gap penalty: -1
- **Input**: cDNA/CDS sequences, GTF annotations, orthologous gene pairs
- **Output**: Orthologous exon groups and isoforms

**EGIO Workflow**:
1. Prepare species-specific sequence data (cDNA, CDS, GTF)
2. Run reciprocal BLASTN between human and species exons
3. Apply dynamic programming alignment
4. Identify orthologous exon groups
5. Generate orthology mappings

#### Orthology-Based EEI Prediction
**Process**:
1. **Data Integration**: Combine EGIO orthology results with species EEI networks
2. **Mapping**: Map orthologous exons between species and human
3. **Prediction**: Transfer EEIs from species to human based on orthology
4. **Validation**: Cross-validate against known human EEIs
5. **Confidence Scoring**: Assign reliability scores based on conservation patterns

**Prediction Algorithm**:
- For each species EEI, identify orthologous human exons
- Transfer EEI if orthology confidence exceeds threshold
- Combine predictions from multiple species
- Apply confidence scoring based on conservation patterns

## Key Innovations

1. **Multi-Method Integration**: Combines three complementary EEI detection approaches
2. **Evolutionary Conservation**: Leverages orthology to expand human EEI networks
3. **Cross-Species Validation**: Uses multiple species for robust predictions
4. **Confidence Scoring**: Provides reliability metrics for predicted interactions
5. **Comprehensive Analysis**: Statistical evaluation and visualization of results

## Usage Instructions

### Running the Complete Pipeline

1. **Human Baseline Construction**:
   ```bash
   cd EEI-Homo-Sapiens/
   # Run scripts in order: 1-Contact-based/, 2-PISA-based/, 3-EPPIC-based/
   bash 4_find_global_EEINs.sh
   ```

2. **Species Network Construction** (for each species):
   ```bash
   cd EEI-[Species-Name]/
   # Run EGIO first
   cd EGIO/
   ./_RUN_egio.sh [parameters]
   # Then run the three detection methods
   ```

3. **Orthology-Based Prediction**:
   ```bash
   cd orthology_based_EEI_prediction/
   python human_EEI_prediction/main.py
   ```

## Results and Analysis Framework

### Human EEI Network Results
**Location**: `EEI-Homo-Sapiens/results/`

**Output Files**:
- `1-CONTACT_results/human_network_final.txt` - Contact-based EEI network
- `2-PISA_results/PISA_EEIN_0.5_human.txt` - PISA-based EEI network
- `3-EPPIC_results/EPPIC_EEIN_filtered.txt` - EPPIC-based EEI network

**Network Statistics**:
- Total EEIs detected by each method
- Overlap between methods
- High-confidence network (NETHIGH) combining all methods

### Species-Specific Network Results
**Location**: `EEI-[Species]/results/`

**For Each Species**:
- Contact-based EEI networks
- PISA-based EEI networks  
- EPPIC-based EEI networks
- Combined species-specific networks

### Orthology-Based Prediction Results
**Location**: `orthology_based_EEI_prediction/results_*/`

**Analysis Components**:

#### Coverage Analysis
- **Coverage Percentage**: Percentage of human EEIs that can be predicted from other species
- **Species Contribution**: Individual species contribution to predictions
- **Method Comparison**: Performance of Contact, PISA, and EPPIC predictions

#### Confidence Analysis
- **Confidence Scoring**: Reliability metrics for predicted EEIs
- **Conservation Patterns**: Evolutionary conservation analysis
- **Validation Results**: Accuracy against known human EEIs

#### Statistical Evaluation
- **Prediction Accuracy**: True positive, false positive rates
- **Coverage Comparison**: Predicted vs. known EEI coverage
- **Conservation Rates**: Conservation patterns across species
- **Identity Distributions**: Sequence identity patterns

### Visualization and Reporting
**Location**: `orthology_based_EEI_prediction/*/visualizations/`

**Generated Figures**:
- Coverage analysis plots
- Confidence distribution histograms
- Conservation pattern visualizations
- Prediction accuracy comparisons
- Species contribution analysis

## Expected Outcomes

This three-pipeline approach is expected to:
- Significantly expand the human EEI network beyond direct detection
- Provide confidence-scored predictions based on evolutionary conservation
- Enable identification of novel EEIs through cross-species analysis
- Contribute to understanding of protein complex evolution and conservation

### Quantitative Expectations
- **Network Expansion**: 2-5x increase in human EEI network size
- **Prediction Accuracy**: >80% accuracy for high-confidence predictions
- **Species Coverage**: >70% of human EEIs predictable from multiple species
- **Conservation Insights**: Identification of evolutionarily conserved interaction patterns

## Contact and Support

For questions about the implementation or methodology, please refer to the individual README files in each pipeline folder or contact the project maintainers.

---

*This documentation provides a comprehensive overview of the EEI Conservation project. For detailed implementation instructions, refer to the README files in each pipeline folder.*
