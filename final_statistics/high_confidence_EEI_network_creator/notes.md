## Tasks:
1. write scripts to create high confidence EEI network between the 3 EEI network of mouse


output:
(geo_env) bbf3630@mario:/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/high_confidence_EEI_network_creator $ python main.py 
Loading Contact network from: ../EEI-Mus-Musculus/results/Mus_musculus/CONTACT_networks/CONTACT_net_6_1.txt
  Found 2321 unique EEI pairs in Contact network
Loading EPPIC network from: ../EEI-Mus-Musculus/3-EPPIC-based/results/EPPIC_EEIN_filtered.txt
  Found 1952 unique EEI pairs in EPPIC network
Loading PISA network from: ../EEI-Mus-Musculus/PISA-based/PISA_networks_filtered/PISA_EEIN_0.5.txt
  Found 1387 unique EEI pairs in PISA network

Finding high-confidence EEIs (present in all three networks)...
  Found 524 high-confidence EEI pairs

Pairwise overlaps:
  Contact ∩ EPPIC: 1097 EEIs
  Contact ∩ PISA: 829 EEIs
  EPPIC ∩ PISA: 926 EEIs

Creating comprehensive high-confidence network...

High-confidence network saved to: ../EEI-Mus-Musculus/high_confidence_network/high_confidence_eei_network_mm.txt

Summary:
  Total high-confidence EEIs: 524
  Unique protein pairs: 80

First 5 high-confidence EEIs:
                  exon1               exon2       protein1       protein2 contact_jaccard_percent contact_exon1_coverage_percent contact_exon2_coverage_percent eppic_buried_area eppic_cs_score eppic_cr_score eppic_pdbid pisa_free_energy pisa_buried_area pisa_hydrogen_bonds pisa_salt_bridges pisa_pdbid
210  ENSMUSE00000217453  ENSMUSE00001217049  B8ZXI1_7b2i_A  Q9JMA2_7b2i_C                    4.76                           7.69                            2.7            285.32          -1.66           0.34        7b2i      -2.72024794        346.80794                   0                 0       7b2i
421  ENSMUSE00000217453  ENSMUSE00001238106  B8ZXI1_7b2i_A  Q9JMA2_7b2i_C                       5                           3.17                           8.11            540.35          -1.66           0.34        7b2i      -5.57342206        612.03693                   1                 0       7b2i
314  ENSMUSE00000217453  ENSMUSE00001304607  B8ZXI1_7b2i_A  Q9JMA2_7b2i_C                      20                          38.89                          10.81            561.13          -1.66           0.34        7b2i      -5.49762632        609.92201                   0                 0       7b2i
288  ENSMUSE00000217455  ENSMUSE00001465794  B8ZXI1_7b2i_A  Q9JMA2_7b2i_C                   17.35                          13.33                          30.43           1111.67          -1.66           0.34        7b2i      -7.16560874        968.11135                   0                 1       7b2i
342  ENSMUSE00000217456  ENSMUSE00001217049  B8ZXI1_7b2i_A  Q9JMA2_7b2i_C                    5.45                           3.85                            6.9            460.53          -1.66           0.34        7b2i      -3.44850472        520.91476                   0                 0       7b2i
