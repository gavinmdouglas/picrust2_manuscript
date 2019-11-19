## PICRUSt2 manuscript repository
This repository contains the scripts and key datafiles used in the PICRUSt2 manuscript. Note that the actual PICRUSt2 code is not in this repository and can be found [here](https://github.com/picrust/picrust2) instead.

PICRUSt2 is described in this [bioRxiv pre-print](https://www.biorxiv.org/content/10.1101/672295v1).



### Repository outline

The ```scripts``` directory is split into three sub-directories:

* ```analyses``` - scripts used for creating key intermediate files for plotting, summarizing data, and calculating statistics for the paper.

* ```figures``` - scripts for plotting main-text and supplementary figures. Each script plots a different figure, which can be re-made using the intermediate files in this repository.

* ```processing_scripts``` - scripts used for processing raw data, either for validation datasets or reference files for PICRUSt2 database.

* The script ```picrust2_ms_functions.R``` is also found in the ```scripts``` directory and contains key functions re-used across multiple scripts in this repository (such as for reading in all prediction tables for a given dataset).


The ```data``` directory contains a number of intermediate files that are called by scripts in this repository. These include all necessary files for re-generating the manuscript figures.

The two most important sub-directories are:

* ```16S_validation``` - contains all prediction tables for all 16S rRNA gene validation datasets.

* ```mgs_validation``` - contains all HUMAnN2 output tables for paired shotgun metagenomics data matching the same 16S rRNA gene samples as in ```16S_validation```.
