# Script: Spatial effect size
# Date: 11th November, 2024


# This script is used to calculate Fisher's r-to-z transformed correlation coefficient

# load libraries
library(metafor)
library(ggplot2)

# Clear environment
rm(list=ls())

# From now on, the analysis is divided based on the dataset used: 
# the "FAL dataset" and the "pairwise FAL-BM dataset"
#----------------------
## FAL dataset
#----------------------

# Importing datasets
# Load data
data <- read.csv("Results/Spatial/Spatial_FAL_wcor.csv") # Load data with weighted corr.
head(data)

# Calculate effect size: 
## Tmax 
data_Tmax <- escalc(measure = "ZCOR", ri = r_Tmax, ni = n, data = data)
head(data_Tmax)

# Add yi_Tmax to the main data
data$yi_Tmax <- data_Tmax$yi

rm(list = c("data_Tmax"))

## Tmin
data_Tmin <- escalc(measure = "ZCOR", ri = r_Tmin, ni = n, data = data)
head(data_Tmin)

# Add yi_Tmin to the main data
data$yi_Tmin <- data_Tmin$yi

rm(list = c("data_Tmin"))

## P.Seasonality
data_SP <- escalc(measure = "ZCOR", ri = r_SP, ni = n, data = data)
head(data_SP)

# Add yi_SP to the main data
data$yi_SP <- data_SP$yi

# Add vi to the main data (it is the same for all the env. variables)
data$vi <- data_SP$vi

write.csv(data, "Results/Spatial/Spatial_FAL_effect_size.csv", row.names = F)
rm("data_SP")

#--------------------------
## Pairwise FAL-BM dataset 
#--------------------------
# Importing datasets
# Load data
p_data <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_wcor.csv") #load data with weighted corr
head(p_data)

# Calculate effect size:
## Tmax 
data_Tmax_FAL <- escalc(measure = "ZCOR", ri = r_Tmax_FAL, ni = n, data = p_data)
data_Tmax_BM  <- escalc(measure = "ZCOR", ri = r_Tmax_BM, ni  = n, data = p_data)

# Add yi_Tmax to the main data
p_data$yi_Tmax_FAL <- data_Tmax_FAL$yi
p_data$yi_Tmax_BM  <- data_Tmax_BM$yi

rm(list = c("data_Tmax_FAL", "data_Tmax_BM"))

## Tmin
data_Tmin_FAL <- escalc(measure = "ZCOR", ri = r_Tmin_FAL, ni = n, data = p_data)
data_Tmin_BM  <- escalc(measure = "ZCOR", ri = r_Tmin_BM, ni  = n, data = p_data)

# Add yi_Tmin to the main data
p_data$yi_Tmin_FAL <- data_Tmin_FAL$yi
p_data$yi_Tmin_BM  <- data_Tmin_BM$yi

rm(list = c("data_Tmin_FAL", "data_Tmin_BM"))

## P.Seasonality
data_SP_FAL <- escalc(measure = "ZCOR", ri = r_SP_FAL, ni = n, data = p_data)
data_SP_BM  <- escalc(measure = "ZCOR", ri = r_SP_BM, ni  = n, data = p_data)

# Add yi_SP to the main data
p_data$yi_SP_FAL <- data_SP_FAL$yi
p_data$yi_SP_BM  <- data_SP_BM$yi

# Add vi (variance) to the main data (it is the same for all the variables)
p_data$vi <- data_SP_FAL$vi


write.csv(p_data, "Results/Spatial/Spatial_pairwiseFAL_BM_effect_size.csv", row.names = F)
rm(list = c("data_SP_FAL", "data_SP_BM"))
