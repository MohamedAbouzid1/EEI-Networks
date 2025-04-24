### EEI in Homo Sapiens

In this folder it's contained all the necessary for building EEI-networks for Homo Sapiens data collected in 2024.

There are three folders, each collecting scripts for each of the three methods used in the reference paper: **Contact-based**, **Energy-based**, **Evolution-based**. 

The folders are numbered and should be run in ascending order, as should the scripts inside them.

 *4_find_global_EEINs.sh* script generate global EEI networks at different confidence levels, based on the three EEI definition approaches. In particular, the high confidence network (NETHIGH), with the smallest coverage and having edges supported by all three EEI definition approaches. It uses the output files of the scripts in the folders, it should be run after all of them.

*Please notice that:*
+ *the scripts from 1a to 7 in **Contact-based** folder are needed for preparation of data. The resulting file of script 7 is indeed subsequently used by all three methods for building the EEI network.*
+ *the two scripts number 5 in 1-Contact-based folder are mutually exclusive. Only one of the two have to be executed.*
