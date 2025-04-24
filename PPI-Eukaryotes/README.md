### PPI IN EUKARYOTES

In this folder it's contained all the necessary for building PPI-networks for UniProt data referring to available Eukaryotes species, collected in 2024.

There are three folders, each collecting scripts for each of the three methods used in the reference paper: Contact-based, Energy-based, Evolution-based.

The folders are numbered and should be run in ascending order, as should the scripts inside them.

Please notice that:

+ due to the large amount of data, all analysis, in particular PISA- and EPPIC-based, might take a lot.

#### **EXTRA SCRIPTS:**
+ *species_detection.sh* creates, for each species, a txt file listing PDBs in which at least one protein of that species is contained. It uses as input the *processed_complex_file.txt* created in script *1-Contact-based/2_Preprocess_database.r*.
