# Script: Spatial weighted correlation 
# Date: 11th November, 2024


# This script performs weighted correlation analysis on environmental data 
# for 'FAL' and 'pairwise FAL-BM' datasets. It calculates the weighted Spearman 
# correlations of FAL and BM with time (variable "Year"). 

#Library
library(wCorr)
library(ggplot2) 
library(mapview)
library(ggpubr)
library(dplyr)
library(rlang)

# Clear environment
rm(list=ls())

# From now on, the analysis is divided based on the dataset used
#----------------------
## FAL dataset
#----------------------
## Upload the dataset
data <- read.csv("Data/Temporal/FAL_data.csv")

# Weighted correlation
FAL_wcor <- data %>%
  group_by(Bat_species, cell_id, fal) %>%
  dplyr::summarize(
    r = as.numeric(format(round(weightedCorr(y = wFAL_mean, x = Year,  method = "spearman", weights = log10(wFAL_n + 1)), 3), nsmall = 3)),
    n = n())

write.csv(FAL_wcor, "Results/Temporal/Temporal_FAL_wcor.csv", row.names = F)
rm(list = c("data"))

#--------------------------
## Pairwise FAL-BM dataset 
#--------------------------
## Upload the dataset
data <- read.csv("Data/Temporal/pairwiseFAL_BM_data.csv")

# Weighted correlation
subFAL_wcor <- data %>%
  group_by(Bat_species, cell_id, fal, bm) %>%
  dplyr::summarize(
    r_FAL = as.numeric(format(round(weightedCorr(y = wFAL_mean, x = Year,  method = "spearman", weights = log10(wFAL.wBM_n + 1)), 3), nsmall = 3)),
    r_BM  = as.numeric(format(round(weightedCorr(y = wBM_mean,  x = Year,  method = "spearman", weights = log10(wFAL.wBM_n + 1)), 3), nsmall = 3)),
    n = n())

write.csv(subFAL_wcor, "Results/Temporal/Temporal_pairwiseFAL_BM_wcor.csv", row.names = F)
rm(list = c("data"))

