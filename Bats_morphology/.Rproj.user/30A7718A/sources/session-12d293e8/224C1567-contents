# Script: Spatial Phylogenetic meta-analysis
# Date: 11th November, 2024


# This script is used to run the phylogenetic meta-analysis  
# based on the Spearman correlations between the environmental variable (year)
# and bats body size measurements (FAL and BM) for both "FAL" and "Pairwise FAL-BM" datasets

# load libraries
library(metafor)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readxl)
library(mapview)
library(sf)
library(raster)

# Clear environment
rm(list=ls())
# From now on, the analysis is divided based on the dataset used:
#----------------------
## FAL dataset 
#----------------------
#load data
data <- read.csv("Data/Temporal/FAL_ph.csv", stringsAsFactors = F)
head(data)

# loading phylogenetic matrixes 
load("Results/Temporal/Temporal_FAL_bat_phylo_cor.Rdata") #bat_phylo_cor

# define phylo vcov matrix and random effects
phylocor <- list(Bat_species = bat_phylo_cor)
RE = list( ~1|Bat_species, ~1|SPID)

# Phylogenetic meta-analyses
# Run phylogenetic meta-analysis based on Spearman Correlations
m1 <- rma.mv(yi_s, vi_s, data = data, random = RE, R = phylocor)
summary(m1)
m1Y <- data.frame(ES=as.numeric(m1$b), CIL=as.numeric(m1$ci.lb), CUL=as.numeric(m1$ci.ub), Variable = "Year")

round(transf.ztor(m1$beta) , digits=3) # estimate
round(transf.ztor(m1$ci.lb), digits=3) # ci lower
round(transf.ztor(m1$ci.ub), digits=3) # ci upper
round(m1$pval, digits=3) #pval
round(m1$QE, digits=1) # QE
m1$QEp # pval for QE


write.csv(m1Y, "Results/Temporal/Temporal_FAL_meta-analysis.csv", row.names = FALSE)
rm(list = c("m1", "m1Y"))

#--------------------------
## Pairwise FAL-BM dataset 
#--------------------------
#load data
data <- read.csv("Data/Temporal/pairwiseFAL_BM_ph.csv", stringsAsFactors = F)
head(data)

# loading phylogenetic matrixes 
load("Results/Temporal/Temporal_pairwiseFAL_BM_bat_phylo_cor.Rdata") #bat_phylo_cor

# define phylo vcov matrix and random effects
phylocor <- list(Bat_species = bat_phylo_cor)
RE = list( ~1|Bat_species, ~1|SPID)

# Phylogenetic meta-analyses
# Run phylogenetic meta-analysis based on Spearman Correlations
 
#FAL
m2 <- rma.mv(yi_s_FAL, vi_s, data = data, random = RE, R = phylocor)
summary(m2)
m2Y <- data.frame(ES=as.numeric(m2$b), CIL=as.numeric(m2$ci.lb), CUL=as.numeric(m2$ci.ub), Variable = "Year")

round(transf.ztor(m2$beta) , digits=3) # estimate
round(transf.ztor(m2$ci.lb), digits=3) # ci lower
round(transf.ztor(m2$ci.ub), digits=3) # ci upper
round(m2$pval, digits=3) #pval
round(m2$QE, digits=1) # QE
m2$QEp # pval for QE


write.csv(m2Y, "Results/Temporal/Temporal_pairwiseFAL_BM_meta-analysis_FAL.csv", row.names = FALSE)
rm(list = c("m2", "m2Y"))

#BM
m3 <- rma.mv(yi_s_BM, vi_s, data = data, random = RE, R = phylocor)
summary(m3)
m3Y <- data.frame(ES=as.numeric(m3$b), CIL=as.numeric(m3$ci.lb), CUL=as.numeric(m3$ci.ub), Variable = "Year")

round(transf.ztor(m3$beta) , digits=3) # estimate
round(transf.ztor(m3$ci.lb), digits=3) # ci lower
round(transf.ztor(m3$ci.ub), digits=3) # ci upper
round(m3$pval, digits=3) #pval
round(m3$QE, digits=1) # QE
m3$QEp # pval for QE


write.csv(m3Y, "Results/Temporal/Temporal_pairwiseFAL_BM_meta-analysis_BM.csv", row.names = FALSE)
rm(list = c("m3", "m3Y"))

