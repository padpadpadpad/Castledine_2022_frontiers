# Analyses and data of:

_Castledine et al. (2022) Greater phage genotypic diversity constrains arms-race coevolution. Frontiers in Cellular and Infection Microbiology. 12: 834406_. 

DOI of paper: [https://doi.org/10.3389%2Ffcimb.2022.834406](https://doi.org/10.3389%2Ffcimb.2022.834406)

ENA Accession number for sequencing: [PRJEB50009](https://www.ebi.ac.uk/ena/browser/view/PRJEB50009)

### Outline

This repository contains the final datasets, analyses and figures of the above-mentioned paper. It can recreate all of the analyses and figures in the main text (Figure 1 and 2) and all the tables in the Supplementary Information

### Feedback

- Please report any problems or bugs in the code in the [Issues](https://github.com/padpadpadpad/Castledine_2022_frontiers) tab of the GitHub repository. Alternatively, please email _d.padfield@exeter.ac.uk_.

### Licensing

This code is licensed under GPL-3.

### Running the scripts and analyses

- The project can be `cloned` or for those not familiar with GitHub, a zip file of this project can be downloaded using the "Clone or download" button at the top right of this page.
- All of the scripts for the analyses can be found in `scripts/`.
- `scripts/map_phage_sequencing.sh` provides the code used to map the re-sequenced phage back to the reference genome and use freebayes to call genetic variants.
- Before running any of the R scripts, open the R project file in the downloaded folder. [R projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects) automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. These are very easy to open using RStudio. 
- `scripts/phenotypic_analysis.R` analyses the time-shift assays of phage infectivity and resistance of different phage treatments. Recreates Figures 1 & 2 and Tables S1 to S3.
- `scripts/check_phage_snps.R` annotates the phage SNPs and indels identified by mapping the sequencing to the phage reference genome. Recreates Table 1.
- `scripts/map_phage.sh` provides the script used for cleaning, mapping, and variant calling of the re-sequencing of phage 1 and phage 2.
