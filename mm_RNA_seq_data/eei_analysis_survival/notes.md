we ran the pipeline once on the orthologous EEIs for human BRCA
now we want to try to run it on the high confidence EEIs from mouse


---
### Output for high confidence mouse data not adjusted p-value
/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/mm_RNA_seq_data/eei_analysis_survival $ python main.py 
=== EEI VALIDATION USING TREATMENT RESPONSE ===
Loading expression data from ../crc2022_exon_expression/crc2022_exon_cpm.tsv...
Expression data shape: (588565, 25)
First few exon IDs: ['ENSMUSE00001334242', 'ENSMUSE00001327336', 'ENSMUSE00001328607', 'ENSMUSE00001323371', 'ENSMUSE00001326783']
Sample IDs (first 5): ['SRR17179989', 'SRR17179988', 'SRR17179987', 'SRR17179986', 'SRR17179985']
Loading metadata from ../crc2022_exon_expression/samplesheet.csv...
Metadata shape: (25, 13)
Metadata columns: ['geo_accession', 'series', 'title', 'model', 'treatment', 'timepoint', 'replicate', 'Run', 'LibraryLayout', 'fastq_1', 'fastq_2', 'bam', 'strand']
treatment_binary counts -> Control: 5, Treated: 20

Response group distribution:
response_group
moderate_responder    10
good_responder        10
poor_control           5
Name: count, dtype: int64

Analyzing EEI associations with treatment response...

Identifying time-dependent EEIs...

Performing pseudo-survival analysis...

Generating validation plots...
Validation plots saved to ../results_survival_proxy_high_confidence_mm/eei_response_validation_plots.png
Saved response associations to ../results_survival_proxy_high_confidence_mm/eei_response_associations.tsv
Saved time-dependent EEIs to ../results_survival_proxy_high_confidence_mm/time_dependent_eeis.tsv

=== ANALYSIS SUMMARY ===
Total EEIs analyzed: 524
EEIs with expression data: 521
Significant EEIs (p<0.05): 227
Time-dependent EEIs: 28
Top response-associated EEIs: 227

---

### Output for high confidence mouse data adjusted p-value
(geo_env) bbf3630@mario:/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/mm_RNA_seq_data/eei_analysis_survival $ python main.py 
=== EEI VALIDATION USING TREATMENT RESPONSE ===
Loading expression data from ../crc2022_exon_expression/crc2022_exon_cpm.tsv...
Expression data shape: (588565, 25)
First few exon IDs: ['ENSMUSE00001334242', 'ENSMUSE00001327336', 'ENSMUSE00001328607', 'ENSMUSE00001323371', 'ENSMUSE00001326783']
Sample IDs (first 5): ['SRR17179989', 'SRR17179988', 'SRR17179987', 'SRR17179986', 'SRR17179985']
Loading metadata from ../crc2022_exon_expression/samplesheet.csv...
Metadata shape: (25, 13)
Metadata columns: ['geo_accession', 'series', 'title', 'model', 'treatment', 'timepoint', 'replicate', 'Run', 'LibraryLayout', 'fastq_1', 'fastq_2', 'bam', 'strand']
treatment_binary counts -> Control: 5, Treated: 20

Response group distribution:
response_group
moderate_responder    10
good_responder        10
poor_control           5
Name: count, dtype: int64

Analyzing EEI associations with treatment response...

Identifying time-dependent EEIs...

Performing pseudo-survival analysis...

Generating validation plots...
Validation plots saved to ../results_survival_proxy_high_confidence_mm/eei_response_validation_plots.png
Saved response associations to ../results_survival_proxy_high_confidence_mm/eei_response_associations.tsv
Saved time-dependent EEIs to ../results_survival_proxy_high_confidence_mm/time_dependent_eeis.tsv

=== ANALYSIS SUMMARY ===
Total EEIs analyzed: 524
EEIs with expression data: 521
Significant EEIs (FDR<0.05): 88
Time-dependent EEIs: 28
Top response-associated EEIs: 88

### Output for CRPES mouse EEIs adjusted P-Value

(geo_env) bbf3630@mario:/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/mm_RNA_seq_data/eei_analysis_survival $ python main.py 
=== EEI VALIDATION USING TREATMENT RESPONSE ===
Loading expression data from ../crc2022_exon_expression/crc2022_exon_cpm.tsv...
Expression data shape: (588565, 25)
First few exon IDs: ['ENSMUSE00001334242', 'ENSMUSE00001327336', 'ENSMUSE00001328607', 'ENSMUSE00001323371', 'ENSMUSE00001326783']
Sample IDs (first 5): ['SRR17179989', 'SRR17179988', 'SRR17179987', 'SRR17179986', 'SRR17179985']
Loading metadata from ../crc2022_exon_expression/samplesheet.csv...
Metadata shape: (25, 13)
Metadata columns: ['geo_accession', 'series', 'title', 'model', 'treatment', 'timepoint', 'replicate', 'Run', 'LibraryLayout', 'fastq_1', 'fastq_2', 'bam', 'strand']
treatment_binary counts -> Control: 5, Treated: 20

Response group distribution:
response_group
moderate_responder    10
good_responder        10
poor_control           5
Name: count, dtype: int64

Analyzing EEI associations with treatment response...

Identifying time-dependent EEIs...

Performing pseudo-survival analysis...

Generating validation plots...
Validation plots saved to ../results_survival_proxy_mouse_crpes_brca_corrected_pvalues/eei_response_validation_plots.png
Saved response associations to ../results_survival_proxy_mouse_crpes_brca_corrected_pvalues/eei_response_associations.tsv
Saved time-dependent EEIs to ../results_survival_proxy_mouse_crpes_brca_corrected_pvalues/time_dependent_eeis.tsv

=== ANALYSIS SUMMARY ===
Total EEIs analyzed: 331
EEIs with expression data: 269
Significant EEIs (FDR<0.05): 42
Time-dependent EEIs: 44
Top response-associated EEIs: 42

