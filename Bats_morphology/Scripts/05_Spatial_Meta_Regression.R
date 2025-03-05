# Script: Spatial Phylogenetic meta-regression
# Date: 11th November, 2024


# This script is used to run the phylogenetic meta-regression  
# between the Spearman correlations and the species-level body sizes 
#for both "FAL" and "Pairwise FAL-BM" datasets

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
# the "FAL dataset" and the "pairwise FAL-BM dataset"
#----------------------
## FAL dataset
#----------------------
#load data
data <- read.csv("Data/Spatial/FAL_ph.csv", stringsAsFactors = F)
head(data)

# separate by sex
data_m <- subset(data, Sex == "male")
data_f <- subset(data, Sex == "female")


# loading phylogenetic matrixes 
load("Results/Spatial/Spatial_FAL_bat_phylo_cor.Rdata") #bat_phylo_cor

# define phylo vcov matrix and random effects
phylocor <- list(Bat_species = bat_phylo_cor)
RE = list( ~1|Bat_species, ~1|SPID)

## Tmax
# effect of size
#males
em1 <- rma.mv(yi_Tmax ~ log10(fal), vi, data = data_m, random = RE, R = phylocor)
summary(em1)

temp <- log10(seq(from = min(data_m$fal), to = max(data_m$fal), length.out = 1000))
em1_Tmax <- predict(em1, newmods = cbind(temp), addx=TRUE)
em1_Tmax <- data.frame(em1_Tmax)
em1_Tmax$temp <- em1_Tmax$X.temp

em1_Tmax$cor <- transf.ztor(em1_Tmax$pred)
em1_Tmax$cor_ci.lb <- transf.ztor(em1_Tmax$ci.lb)
em1_Tmax$cor_ci.ub <- transf.ztor(em1_Tmax$ci.ub)

write.csv(em1_Tmax, "Results/Spatial/Spatial_FAL_meta-regression_m_Tmax.csv", row.names = FALSE)

#females
ef1 <- rma.mv(yi_Tmax ~ log10(fal), vi, data = data_f, random = RE, R = phylocor)
summary(ef1)

temp <- log10(seq(from = min(data_f$fal), to = max(data_f$fal), length.out = 1000))

ef1_Tmax <- predict(ef1, newmods = cbind(temp), addx=TRUE)
ef1_Tmax <- data.frame(ef1_Tmax)
ef1_Tmax$temp <- ef1_Tmax$X.temp

ef1_Tmax$cor <- transf.ztor(ef1_Tmax$pred)
ef1_Tmax$cor_ci.lb <- transf.ztor(ef1_Tmax$ci.lb)
ef1_Tmax$cor_ci.ub <- transf.ztor(ef1_Tmax$ci.ub)

write.csv(ef1_Tmax, "Results/Spatial/Spatial_FAL_meta-regression_f_Tmax.csv", row.names = FALSE)


#--------------------------
## Pairwise FAL-BM dataset
#--------------------------
rm(list=ls())
#load data
data <- read.csv("Data/Spatial/pairwiseFAL_BM_ph.csv", stringsAsFactors = F)
head(data)


# separate by sex
data_m <- subset(data, Sex == "male")
data_f <- subset(data, Sex == "female")


# loading phylogenetic matrixes 
load("Results/Spatial/Spatial_pairwiseFAL_BM_bat_phylo_cor.Rdata") #bat_phylo_cor

# define phylo vcov matrix and random effects
phylocor <- list(Bat_species = bat_phylo_cor)
RE = list( ~1|Bat_species, ~1|SPID)

## Tmax
# effect of size
#males_FAL
em1 <- rma.mv(yi_Tmax_FAL ~ log10(fal), vi, data = data_m, random = RE, R = phylocor)
summary(em1)

temp <- log10(seq(from = min(data_m$fal), to = max(data_m$fal), length.out = 1000))
em1_Tmax_FAL <- predict(em1, newmods = cbind(temp), addx=TRUE)
em1_Tmax_FAL <- data.frame(em1_Tmax_FAL)
em1_Tmax_FAL$temp <- em1_Tmax_FAL$X.temp

em1_Tmax_FAL$cor <- transf.ztor(em1_Tmax_FAL$pred)
em1_Tmax_FAL$cor_ci.lb <- transf.ztor(em1_Tmax_FAL$ci.lb)
em1_Tmax_FAL$cor_ci.ub <- transf.ztor(em1_Tmax_FAL$ci.ub)

write.csv(em1_Tmax_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-regression_m_Tmax_FAL.csv", row.names = FALSE)

#males_BM
em1 <- rma.mv(yi_Tmax_BM ~ log10(bm), vi, data = data_m, random = RE, R = phylocor)
summary(em1)

temp <- log10(seq(from = min(data_m$bm), to = max(data_m$bm), length.out = 1000))
em1_Tmax_BM <- predict(em1, newmods = cbind(temp), addx=TRUE)
em1_Tmax_BM <- data.frame(em1_Tmax_BM)
em1_Tmax_BM$temp <- em1_Tmax_BM$X.temp

em1_Tmax_BM$cor <- transf.ztor(em1_Tmax_BM$pred)
em1_Tmax_BM$cor_ci.lb <- transf.ztor(em1_Tmax_BM$ci.lb)
em1_Tmax_BM$cor_ci.ub <- transf.ztor(em1_Tmax_BM$ci.ub)

write.csv(em1_Tmax_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-regression_m_Tmax_BM.csv", row.names = FALSE)

#females_FAL
ef1 <- rma.mv(yi_Tmax_FAL ~ log10(fal), vi, data = data_f, random = RE, R = phylocor)
summary(ef1)

temp <- log10(seq(from = min(data_f$fal), to = max(data_f$fal), length.out = 1000))

ef1_Tmax_FAL <- predict(ef1, newmods = cbind(temp), addx=TRUE)
ef1_Tmax_FAL <- data.frame(ef1_Tmax_FAL)
ef1_Tmax_FAL$temp <- ef1_Tmax_FAL$X.temp

ef1_Tmax_FAL$cor <- transf.ztor(ef1_Tmax_FAL$pred)
ef1_Tmax_FAL$cor_ci.lb <- transf.ztor(ef1_Tmax_FAL$ci.lb)
ef1_Tmax_FAL$cor_ci.ub <- transf.ztor(ef1_Tmax_FAL$ci.ub)

write.csv(ef1_Tmax_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-regression_f_Tmax_FAL.csv", row.names = FALSE)

#females_BM
ef1 <- rma.mv(yi_Tmax_BM ~ log10(bm), vi, data = data_f, random = RE, R = phylocor)
summary(ef1)

temp <- log10(seq(from = min(data_f$bm), to = max(data_f$bm), length.out = 1000))

ef1_Tmax_BM <- predict(ef1, newmods = cbind(temp), addx=TRUE)
ef1_Tmax_BM <- data.frame(ef1_Tmax_BM)
ef1_Tmax_BM$temp <- ef1_Tmax_BM$X.temp

ef1_Tmax_BM$cor <- transf.ztor(ef1_Tmax_BM$pred)
ef1_Tmax_BM$cor_ci.lb <- transf.ztor(ef1_Tmax_BM$ci.lb)
ef1_Tmax_BM$cor_ci.ub <- transf.ztor(ef1_Tmax_BM$ci.ub)

write.csv(ef1_Tmax_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-regression_f_Tmax_BM.csv", row.names = FALSE)

