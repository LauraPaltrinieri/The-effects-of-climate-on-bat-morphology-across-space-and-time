# Script: Spatial Phylogenetic meta-regression
# Date: 11th November, 2024


# This script is used to run the phylogenetic meta-regression  
# between the Spearman correlations and the species-level body sizes 
# for both "FAL" and "Pairwise FAL-BM" datasets

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

## Tmax
# effect of size
#males
em1 <- rma.mv(yi_s ~ log10(fal), vi_s, data = data, random = RE, R = phylocor)
summary(em1)

temp <- log10(seq(from = min(data$fal), to = max(data$fal), length.out = 1000))
em1_Y <- predict(em1, newmods = cbind(temp), addx=TRUE)
em1_Y <- data.frame(em1_Y)
em1_Y$temp <- em1_Y$X.temp

em1_Y$cor <- transf.ztor(em1_Y$pred)
em1_Y$cor_ci.lb <- transf.ztor(em1_Y$ci.lb)
em1_Y$cor_ci.ub <- transf.ztor(em1_Y$ci.ub)

write.csv(em1_Y, "Results/Temporal/Temporal_FAL_meta-regression.csv", row.names = FALSE)

#--------------------------
## Pairwise FAL-BM dataset
#--------------------------
rm(list=ls())

#load data
data <- read.csv("Data/Temporal/pairwiseFAL_BM_ph.csv", stringsAsFactors = F)
head(data)

# loading phylogenetic matrixes 
load("Results/Temporal/Temporal_pairwiseFAL_BM_bat_phylo_cor.Rdata") #bat_phylo_cor

# define phylo vcov matrix and random effects
phylocor <- list(Bat_species = bat_phylo_cor)
RE = list( ~1|Bat_species, ~1|SPID)


# effect of size
#FAL
em1 <- rma.mv(yi_s_FAL ~ log10(fal), vi_s, data = data, random = RE, R = phylocor)
summary(em1)

temp <- log10(seq(from = min(data$fal), to = max(data$fal), length.out = 1000))
em1_Y_FAL <- predict(em1, newmods = cbind(temp), addx=TRUE)
em1_Y_FAL <- data.frame(em1_Y_FAL)
em1_Y_FAL$temp <- em1_Y_FAL$X.temp

em1_Y_FAL$cor <- transf.ztor(em1_Y_FAL$pred)
em1_Y_FAL$cor_ci.lb <- transf.ztor(em1_Y_FAL$ci.lb)
em1_Y_FAL$cor_ci.ub <- transf.ztor(em1_Y_FAL$ci.ub)

write.csv(em1_Y_FAL, "Results/Temporal/Temporal_pairwiseFAL_BM_meta-regression_FAL.csv", row.names = FALSE)

#BM
em1 <- rma.mv(yi_s_BM ~ log10(bm), vi_s, data = data, random = RE, R = phylocor)
summary(em1)

temp <- log10(seq(from = min(data$bm), to = max(data$bm), length.out = 1000))
em1_Y_BM <- predict(em1, newmods = cbind(temp), addx=TRUE)
em1_Y_BM <- data.frame(em1_Y_BM)
em1_Y_BM$temp <- em1_Y_BM$X.temp

em1_Y_BM$cor <- transf.ztor(em1_Y_BM$pred)
em1_Y_BM$cor_ci.lb <- transf.ztor(em1_Y_BM$ci.lb)
em1_Y_BM$cor_ci.ub <- transf.ztor(em1_Y_BM$ci.ub)

write.csv(em1_Y_BM, "Results/Temporal/Temporal_pairwiseFAL_BM_meta-regression_BM.csv", row.names = FALSE)



