# Technical Implementation Guide

## Overview
This document provides detailed technical implementation information for the EEI Conservation project, complementing the main documentation.

## Pipeline 1: Human EEI Network Construction

### Contact-Based Detection

#### Data Flow
```
UniProt Data → SIFTS Mapping → PDB Structures → Exon Mapping → Contact Detection → EEI Network
```

#### Key Scripts and Functions

**1. Data Download and Preprocessing**
- `1a_download_uniprot.r`: Downloads human protein data from UniProt
- `1b_preprocess_uniprot.r`: Processes and cleans UniProt data
- `2_download_sifts_mapping.r`: Downloads SIFTS mappings for PDB-UniProt relationships

**2. Structure Processing**
- `3_download_PDBs.r`: Downloads relevant PDB structures
- `4a_create_UniprotPDB_map.r`: Creates mapping between UniProt and PDB
- `4b_add_CIF_num.r`: Adds CIF numbering information

**3. Exon Mapping**
- `5_Ensembl_exon_map.r`: Maps Ensembl exons to protein coordinates
- `6_all_map.r`: Integrates all mapping information

**4. Network Construction**
- `7_chain2net.r`: Converts protein chains to network format
- `8_exonexon_contact.r`: Identifies exon-exon contacts
- `9_save_contact_networks.r`: Saves final contact-based EEI network

#### Technical Parameters
- **Contact Distance Threshold**: 5Å (default)
- **Minimum Contact Residues**: 3
- **Sequence Identity Threshold**: 95% for redundancy removal

### PISA-Based Detection

#### Data Flow
```
PDB Structures → PISA Analysis → Interface Parsing → Exon Mapping → EEI Network
```

#### Key Scripts and Functions

**1. PISA Analysis**
- `7_runPISA.r`: Executes PISA analysis on protein complexes
- `PisaAuto_file.py`: Automated PISA file processing
- `PisaAuto_id.py`: PISA analysis by PDB ID

**2. Result Processing**
- `8_parse_files.r`: Parses PISA output files
- `Pisa_residue_parser.py`: Extracts residue-level interface information
- `Pisa_table_parser.py`: Processes PISA table outputs

**3. Network Construction**
- `9a_map_EEIs.r`: Maps PISA interfaces to EEIs
- `9b_preprocess_EEINs.r`: Preprocesses EEI networks

#### Technical Parameters
- **Interface Stability Threshold**: 0.5 (default)
- **Minimum Interface Area**: 1000 Å²
- **Energy Threshold**: -1.0 kcal/mol

### EPPIC-Based Detection

#### Data Flow
```
PDB Structures → EPPIC Analysis → Conservation Analysis → Exon Mapping → EEI Network
```

#### Key Scripts and Functions

**1. EPPIC Setup**
- `EPPIC_setup.sh`: Configures EPPIC environment
- `config_file.txt`: EPPIC configuration parameters

**2. Analysis Execution**
- `1_run_eppic.r`: Runs EPPIC analysis
- `2_map_EEI.r`: Maps EPPIC results to EEIs
- `3_preprocess_EEINs.r`: Preprocesses EEI networks

#### Technical Parameters
- **Conservation Threshold**: 0.7 (default)
- **Minimum Interface Size**: 500 Å²
- **Evolutionary Distance**: All available species

## Pipeline 2: Species-Specific Networks

### Implementation for Each Species

#### Data Requirements
For each species, the following data is required:
- **Genomic Sequences**: cDNA and CDS FASTA files
- **Annotations**: GTF files with exon coordinates
- **Protein Data**: UniProt sequences and PDB structures
- **Orthology Data**: Gene orthology relationships with human

#### Species-Specific Adaptations

**Mus musculus (Mouse)**
- Extensive PDB coverage
- Well-annotated genome
- High orthology with human

**Bos taurus (Cattle)**
- Agricultural importance
- Mammalian conservation patterns
- Moderate PDB coverage

**Drosophila melanogaster (Fruit fly)**
- Invertebrate model
- Evolutionary distance from human
- Extensive research data

**Gallus gallus (Chicken)**
- Avian model
- Evolutionary perspective
- Moderate PDB coverage

**Oryctolagus cuniculus (Rabbit)**
- Mammalian model
- Research applications
- Limited PDB coverage

**Rattus norvegicus (Rat)**
- Model organism
- Biomedical research
- High orthology with human

**Saccharomyces cerevisiae (Yeast)**
- Simple eukaryote
- Fundamental biology
- Limited PDB coverage

## Pipeline 3: Orthology Mapping and Prediction

### EGIO Implementation

#### Prerequisites
- **BLAST+**: Sequence alignment software
- **GCC**: C compiler for EGIO components
- **Python**: pandas, numpy, plotnine
- **R**: Statistical analysis

#### EGIO Workflow

**1. Data Preparation**
```bash
# Prepare sequence files
python __prepare_egio_blastn.py --species1 hsa --species2 gge
python __prepare_egio_extra.py --species1 hsa --species2 gge
```

**2. BLASTN Analysis**
```bash
# Run reciprocal BLASTN
blastn -query species1_cds.fa -db species2_cds -outfmt 6 > blast_results.txt
blastn -query species2_cds.fa -db species1_cds -outfmt 6 > blast_results_rev.txt
```

**3. EGIO Execution**
```bash
./_RUN_egio.sh -s hsa -S gge -r hsa.gtf -R gge.gtf \
               -e hsa_cdna.fa -E gge_cdna.fa \
               -o hsa_cds.fa -O gge_cds.fa \
               -h ortholog_pairs.txt -p 6 -i 0.8 -c 0.8
```

#### EGIO Parameters
- `-s, -S`: Species names
- `-r, -R`: GTF annotation files
- `-e, -E`: cDNA sequence files
- `-o, -O`: CDS sequence files
- `-h`: Orthologous gene pairs
- `-p`: Number of processors
- `-i`: Identity threshold (0.8)
- `-c`: Coverage threshold (0.8)
- `-m`: Match score (2)
- `-n`: Mismatch penalty (-2)
- `-g`: Gap penalty (-1)

### Orthology-Based Prediction

#### Prediction Algorithm

**1. Data Integration**
```python
# Load EGIO results
egio_results = load_egio_results('ExonGroup.txt', 'OrthoIso.txt')

# Load species EEI networks
species_eeis = load_species_eeis('species_network.txt')

# Load human exon coordinates
human_exons = load_human_exons('human_exons.txt')
```

**2. Orthology Mapping**
```python
# Map orthologous exons
orthology_map = map_orthologous_exons(egio_results, human_exons)

# Filter by confidence
high_confidence_orthologs = filter_by_confidence(orthology_map, threshold=0.8)
```

**3. EEI Transfer**
```python
# Transfer EEIs based on orthology
predicted_eeis = transfer_eeis(species_eeis, orthology_map)

# Apply confidence scoring
scored_predictions = apply_confidence_scoring(predicted_eeis)
```

#### Confidence Scoring Algorithm

**Factors Considered**:
1. **Sequence Identity**: Higher identity = higher confidence
2. **Coverage**: Higher coverage = higher confidence
3. **Conservation**: Multiple species support = higher confidence
4. **Method Agreement**: Multiple methods = higher confidence

**Scoring Formula**:
```
Confidence = (Identity × 0.4) + (Coverage × 0.3) + (Conservation × 0.2) + (Method_Agreement × 0.1)
```

## Database Integration

### EEINet Database
**Location**: `EEINet/`

**Purpose**: Web-based interface for EEI network exploration and analysis

**Components**:
- **Backend**: Node.js with PostgreSQL database
- **Frontend**: Web interface for network visualization
- **API**: RESTful API for data access
- **Analysis Tools**: Statistical analysis and visualization

### Data Schema
```sql
-- EEI Networks Table
CREATE TABLE eei_networks (
    id SERIAL PRIMARY KEY,
    species VARCHAR(50),
    method VARCHAR(20),
    exon1_id VARCHAR(100),
    exon2_id VARCHAR(100),
    confidence FLOAT,
    evidence_type VARCHAR(50)
);

-- Orthology Mappings Table
CREATE TABLE orthology_mappings (
    id SERIAL PRIMARY KEY,
    human_exon_id VARCHAR(100),
    species_exon_id VARCHAR(100),
    species_name VARCHAR(50),
    identity FLOAT,
    coverage FLOAT,
    confidence FLOAT
);
```

## Performance Optimization

### Computational Requirements

**Hardware Requirements**:
- **CPU**: Multi-core processor (8+ cores recommended)
- **RAM**: 32GB+ for large-scale analysis
- **Storage**: 1TB+ for PDB structures and results
- **Network**: High-speed internet for data downloads

**Software Dependencies**:
- **R**: 4.0+ with required packages
- **Python**: 3.8+ with scientific computing packages
- **BLAST+**: 2.10+ for sequence alignment
- **PostgreSQL**: 12+ for database operations

### Optimization Strategies

**1. Parallel Processing**
- Multi-threaded BLASTN analysis
- Parallel PISA/EPPIC execution
- Distributed EGIO processing

**2. Memory Management**
- Chunked data processing
- Efficient data structures
- Garbage collection optimization

**3. Storage Optimization**
- Compressed data storage
- Efficient indexing
- Data archiving strategies

## Quality Control and Validation

### Data Quality Checks

**1. Sequence Quality**
- Sequence length validation
- Ambiguous base filtering
- Duplicate sequence removal

**2. Structure Quality**
- Resolution filtering (>2.5Å)
- Completeness checks
- Interface validation

**3. Annotation Quality**
- Coordinate validation
- Gene model consistency
- Orthology verification

### Validation Metrics

**1. Prediction Accuracy**
- True positive rate
- False positive rate
- Precision and recall

**2. Coverage Analysis**
- Network coverage percentage
- Species contribution
- Method comparison

**3. Conservation Analysis**
- Evolutionary conservation patterns
- Functional conservation
- Structural conservation

## Troubleshooting Guide

### Common Issues

**1. EGIO Execution Failures**
- Check BLAST+ installation
- Verify sequence file formats
- Ensure sufficient disk space

**2. PISA Analysis Errors**
- Verify PDB file integrity
- Check PISA installation
- Monitor memory usage

**3. Database Connection Issues**
- Verify PostgreSQL service
- Check connection parameters
- Validate database schema

### Performance Issues

**1. Slow BLASTN Analysis**
- Increase number of threads
- Optimize database size
- Use faster hardware

**2. Memory Limitations**
- Reduce batch sizes
- Optimize data structures
- Increase system memory

**3. Storage Space**
- Implement data archiving
- Use compression
- Clean temporary files

## Maintenance and Updates

### Regular Maintenance Tasks

**1. Data Updates**
- UniProt data refresh
- PDB structure updates
- Ensembl annotation updates

**2. Software Updates**
- Package version updates
- Security patches
- Performance improvements

**3. Database Maintenance**
- Index optimization
- Data cleanup
- Backup procedures

### Monitoring and Logging

**1. Process Monitoring**
- Execution time tracking
- Resource usage monitoring
- Error logging

**2. Quality Monitoring**
- Result validation
- Performance metrics
- User feedback

**3. Backup Procedures**
- Regular data backups
- Version control
- Disaster recovery

---

*This technical implementation guide provides detailed information for maintaining and extending the EEI Conservation project. For specific implementation questions, refer to the individual script documentation or contact the development team.*
