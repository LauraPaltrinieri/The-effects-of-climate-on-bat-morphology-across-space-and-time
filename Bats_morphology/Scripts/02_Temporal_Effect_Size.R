# Script: Spatial effect size
# Date: 11th November, 2024


# This script is used to calculate Fisher's r-to-z transformed correlation coefficient

# load libraries
library(metafor)
library(ggplot2)

# Clear environment
rm(list=ls())

#----------------------
## FAL 
#----------------------
#load data
data <- read.csv("Results/Temporal/Temporal_FAL_wcor.csv") #load data with weighted corr
head(data)

# Calculate effect size 
data_s <- escalc(measure = "ZCOR", ri = r, ni = n, data = data)
head(data_s)

colnames(data_s)[6] <- "yi_s"
colnames(data_s)[7] <- "vi_s"

write.csv(data_s, "Results/Temporal/Temporal_FAL_effect_size.csv", row.names = F)

#-------------------------- 
## Pairwise FAL-BM dataset
#--------------------------
#load data
sub_data <- read.csv("Results/Temporal/Temporal_pairwiseFAL_BM_wcor.csv") #load data with weighted corr
head(sub_data)

## Calculate the effect size
#Calculate the effect sizes and add them to the main dataframe
data_sFAL <- escalc(measure = "ZCOR", ri = r_FAL, ni = n, data = sub_data)
data_sBM  <- escalc(measure = "ZCOR", ri = r_BM, ni = n, data = sub_data)
head(sub_data)

sub_data$yi_s_FAL <- data_sFAL$yi
sub_data$yi_s_BM  <- data_sBM$yi
sub_data$vi_s <- data_sBM$vi

write.csv(sub_data, "Results/Temporal/Temporal_pairwiseFAL_BM_effect_size.csv", row.names = F)
