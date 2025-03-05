# Script: Spatial weighted correlation 
# Date: 11th November, 2024


# This script performs weighted correlation analysis on environmental data 
# for 'FAL' and 'pairwise FAL-BM' datasets. It calculates the weighted Spearman 
# correlations of FAL and BM with climate variables (Tmax, Tmin, Summer precipitation).

#Library
library(wCorr)
library(ggplot2) 
library(mapview)
library(ggpubr)
library(dplyr)
library(rlang)
library(raster)
# Clear environment
rm(list=ls())

# From now on, the analysis is divided based on the dataset used: 
# the "FAL dataset" and the "pairwise FAL-BM dataset"
#----------------------
## FAL dataset
#----------------------
## Upload the dataset
data <- read.csv("Data/Spatial/FAL_env_var.csv")

# Weighted correlation
FAL_wcor <- data %>%
  group_by(Bat_species, Sex, fal) %>%
  dplyr::summarize(
    r_Tmax = as.numeric(format(round(weightedCorr(y = wFAL_mean, x = Tmax,                 method = "spearman", weights = log10(wFAL_n + 1)), 3), nsmall = 3)),
    r_Tmin = as.numeric(format(round(weightedCorr(y = wFAL_mean, x = Tmin,                 method = "spearman", weights = log10(wFAL_n + 1)), 3), nsmall = 3)),
    r_SP   = as.numeric(format(round(weightedCorr(y = wFAL_mean, x = Summer_precipitation, method = "spearman", weights = log10(wFAL_n + 1)), 3), nsmall = 3)),
    n = n())

write.csv(FAL_wcor, "Results/Spatial/Spatial_FAL_wcor.csv", row.names = F)
rm(list = c("data"))

#--------------------------
## Pairwise FAL-BM dataset 
#--------------------------
## Upload the dataset
data_p <- read.csv("Data/Spatial/pairwiseFAL_BM_env_var.csv")

# Weighted correlation
wcor_p <- data_p %>%
  group_by(Bat_species, Sex, fal, bm) %>%
  dplyr::summarize(
    r_Tmax_FAL = as.numeric(format(round(weightedCorr(y = wFAL_mean, x = Tmax,  method                = "spearman", weights = log10(wFAL.wBM_n + 1)), 3), nsmall = 3)),
    r_Tmin_FAL = as.numeric(format(round(weightedCorr(y = wFAL_mean, x = Tmin,  method                = "spearman", weights = log10(wFAL.wBM_n + 1)), 3), nsmall = 3)),
    r_SP_FAL   = as.numeric(format(round(weightedCorr(y = wFAL_mean, x = Summer_precipitation, method = "spearman", weights = log10(wFAL.wBM_n + 1)), 3), nsmall = 3)),
    r_Tmax_BM  = as.numeric(format(round(weightedCorr(y = wBM_mean,  x = Tmax,  method                = "spearman", weights = log10(wFAL.wBM_n + 1)), 3), nsmall = 3)),
    r_Tmin_BM  = as.numeric(format(round(weightedCorr(y = wBM_mean,  x = Tmin,  method                = "spearman", weights = log10(wFAL.wBM_n + 1)), 3), nsmall = 3)),
    r_SP_BM    = as.numeric(format(round(weightedCorr(y = wBM_mean,  x = Summer_precipitation, method = "spearman", weights = log10(wFAL.wBM_n + 1)), 3), nsmall = 3)),
    n = n())

write.csv(wcor_p, "Results/Spatial/Spatial_pairwiseFAL_BM_wcor.csv", row.names = F)
rm(list = c("data_p"))
