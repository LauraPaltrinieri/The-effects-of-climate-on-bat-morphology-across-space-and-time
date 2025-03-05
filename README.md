This repository contains all the data and scripts needed to reproduce the results of the manuscript: The effects of climate on bat morphology across space and time (https://doi.org/10.1002/ecog.07663).

The main dataset is stored in the file: FAL_raw.csv, which includes over 64,000 forearm length (FAL) measurements for 39 European bat species.

The subset dataset is available in the file: pairwise_FAL_BM_raw.csv, containing more than 46,000 records of both FAL and body mass (BM) for 35 bat species.

Scripts are organized to replicate the analysis presented in the manuscript. The workflow begins with 01_Data_Cleaning.R, where the data is cleaned and processed. The following scripts are divided into spatial and temporal analyses:

Spatial Analysis: These scripts include environmental variables and analyze spatial patterns in bat morphology.
Temporal Analysis: These scripts focus on how bat body size changes over time without incorporating environmental variables.
Each script is divided into two sections:

"FAL dataset": Refers to the analyses performed with the main FAL data.
"pairwise_FAL_BM dataset": Refers to the analyses performed with the pairwise_FAL_BM data, which includes both FAL and BM records.
Correlation coefficients and phylogenetic correlation matrices are calculated, and a phylogenetic meta-analysis is performed. Some scripts include code to reproduce the figures presented in the manuscript.
