# Script: Spatial Phylogenetic meta-analysis
# Date: 11th November, 2024


# This script is used to run the phylogenetic meta-analysis  
# based on the Spearman correlations between the environmental variables (Tmax, Tmin, SP)
# and bats body size measurements (FAL and BM) for both "FAL" and "Pairwise FAL-BM" datasets

# Load libraries
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

# Separate by sex
data_m <- subset(data, Sex == "male")
data_f <- subset(data, Sex == "female")

# Loading phylogenetic matrices 
load("Results/Spatial/Spatial_FAL_bat_phylo_cor.Rdata") # Bat_phylo_cor

# Define phylo vcov matrix and random effects
phylocor <- list(Bat_species = bat_phylo_cor) # Phylogenetic vcov matrix
RE = list( ~1|Bat_species, ~1|SPID) # Random effects

# Phylogenetic meta-analyses

# Run phylogenetic meta-analysis based on Spearman Correlations
# Tmax variable 
# Males
m1 <- rma.mv(yi_Tmax, vi, data = data_m, random = RE, R = phylocor)
summary(m1)
m1Tmax <- data.frame(ES=as.numeric(m1$b), CIL=as.numeric(m1$ci.lb), CUL=as.numeric(m1$ci.ub), Variable = "Tmax")

round(transf.ztor(m1$beta) , digits=3) # estimate
round(transf.ztor(m1$ci.lb), digits=3) # ci lower
round(transf.ztor(m1$ci.ub), digits=3) # ci upper
round(m1$QE, digits=3) # pval
round(m1$QE, digits=1) # QE
m1$QEp # pval for QE

# Create a subset df that includes only Tmax
data_mTmax <- data_m %>%
  dplyr::select(-yi_Tmin, -yi_SP) %>%
  dplyr::rename(yi_s = yi_Tmax)
data_mTmax$Variable <- "Tmax"

# Save the results
write.csv(m1Tmax, "Results/Spatial/Spatial_FAL_meta-analysis_m_Tmax.csv", row.names = FALSE)
write.csv(data_mTmax, "Results/Spatial/Spatial_FAL_data_m_Tmax.csv", row.names = FALSE)
rm(list = c("data", "data_mTmax", "m1", "m1Tmax"))

# Females
f1 <- rma.mv(yi_Tmax, vi, data = data_f, random = RE, R = phylocor)
summary(f1)
f1Tmax <- data.frame(ES=as.numeric(f1$b), CIL=as.numeric(f1$ci.lb), CUL=as.numeric(f1$ci.ub), Variable = "Tmax")

round(transf.ztor(f1$beta) , digits=3) # estimate
round(transf.ztor(f1$ci.lb), digits=3) # ci lower
round(transf.ztor(f1$ci.ub), digits=3) # ci upper
round(f1$QE, digits=3) #pval
round(f1$QE, digits=1) # QE
f1$QEp # pval for QE

# Create a subset df that includes only Tmax
data_fTmax <- data_f %>%
  dplyr::select(-yi_Tmin, -yi_SP) %>%
  dplyr::rename(yi_s = yi_Tmax)
data_fTmax$Variable <- "Tmax"

# Save the results
write.csv(f1Tmax, "Results/Spatial/Spatial_FAL_meta-analysis_f_Tmax.csv", row.names = FALSE)
write.csv(data_fTmax, "Results/Spatial/Spatial_FAL_data_f_Tmax.csv", row.names = FALSE)
rm(list = c("data_fTmax", "f1", "f1Tmax"))


# Tmin  
# Males
m2 <- rma.mv(yi_Tmin, vi, data = data_m, random = RE, R = phylocor)
summary(m2)
m2Tmin <- data.frame(ES=as.numeric(m2$b), CIL=as.numeric(m2$ci.lb), CUL=as.numeric(m2$ci.ub), Variable = "Tmin")

round(transf.ztor(m2$beta) , digits=3) # estimate
round(transf.ztor(m2$ci.lb), digits=3) # ci lower
round(transf.ztor(m2$ci.ub), digits=3) # ci upper
round(m2$QE, digits=3) #pval
round(m2$QE, digits=1) # QE
m2$QEp # p val for QE

data_mTmin <- data_m %>%
  dplyr::select(-yi_Tmax, -yi_SP) %>%
  dplyr::rename(yi_s = yi_Tmin)
data_mTmin$Variable <- "Tmin"

write.csv(m2Tmin, "Results/Spatial/Spatial_FAL_meta-analysis_m_Tmin.csv", row.names = FALSE)
write.csv(data_mTmin, "Results/Spatial/Spatial_FAL_data_m_Tmin.csv", row.names = FALSE)
rm(list = c("data_mTmin", "m2", "m2Tmin"))

#females
f2 <- rma.mv(yi_Tmin, vi, data = data_f, random = RE, R = phylocor)
summary(f2)
f2Tmin <- data.frame(ES=as.numeric(f2$b), CIL=as.numeric(f2$ci.lb), CUL=as.numeric(f2$ci.ub), Variable = "Tmin")

round(transf.ztor(f2$beta) , digits=3) # estimate
round(transf.ztor(f2$ci.lb), digits=3) # ci lower
round(transf.ztor(f2$ci.ub), digits=3) # ci upper
round(f2$QE, digits=3) #pval
round(f2$QE, digits=1) # QE
f2$QEp # p val for QE

data_fTmin <- data_f %>%
  dplyr::select(-yi_Tmax, -yi_SP) %>%
  dplyr::rename(yi_s = yi_Tmin)
data_fTmin$Variable <- "Tmin"

write.csv(f2Tmin, "Results/Spatial/Spatial_FAL_meta-analysis_f_Tmin.csv", row.names = FALSE)
write.csv(data_fTmin, "Results/Spatial/Spatial_FAL_data_f_Tmin.csv", row.names = FALSE)

# SP  
# Males
m3 <- rma.mv(yi_SP, vi, data = data_m, random = RE, R = phylocor)
summary(m3)
m3SP <- data.frame(ES=as.numeric(m3$b), CIL=as.numeric(m3$ci.lb), CUL=as.numeric(m3$ci.ub), Variable = "S. Precipitation")

round(transf.ztor(m3$beta) , digits=3) # estimate
round(transf.ztor(m3$ci.lb), digits=3) # ci lower
round(transf.ztor(m3$ci.ub), digits=3) # ci upper
round(m3$QE, digits=3) #pval
round(m3$QE, digits=1) # QE
m3$QEp # p val for QE

data_mSP <- data_m %>%
  dplyr::select(-yi_Tmax, -yi_Tmin) %>%
  dplyr::rename(yi_s = yi_SP)
data_mSP$Variable <- "S. Precipitation"

write.csv(m3SP, "Results/Spatial/Spatial_FAL_meta-analysis_m_SP.csv", row.names = FALSE)
write.csv(data_mSP, "Results/Spatial/Spatial_FAL_data_m_SP.csv", row.names = FALSE)
rm(list = c("data_mSP", "m3", "m3SP"))

#females
f3 <- rma.mv(yi_SP, vi, data = data_f, random = RE, R = phylocor)
summary(f3)
f3SP <- data.frame(ES=as.numeric(f3$b), CIL=as.numeric(f3$ci.lb), CUL=as.numeric(f3$ci.ub), Variable = "S. Precipitation")

round(transf.ztor(f3$beta) , digits=3) # estimate
round(transf.ztor(f3$ci.lb), digits=3) # ci lower
round(transf.ztor(f3$ci.ub), digits=3) # ci upper
round(f3$QE, digits=3) #pval
round(f3$QE, digits=1) # QE
f3$QEp # p val for QE

data_fSP <- data_f %>%
  dplyr::select(-yi_Tmax, -yi_Tmin) %>%
  dplyr::rename(yi_s = yi_SP)
data_fSP$Variable <- "S. Precipitation"

write.csv(f3SP, "Results/Spatial/Spatial_FAL_meta-analysis_f_SP.csv", row.names = FALSE)
write.csv(data_fSP, "Results/Spatial/Spatial_FAL_data_f_SP.csv", row.names = FALSE)


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

# Phylogenetic meta-analyses

# Run phylogenetic meta-analysis based on Spearman Correlations
# Tmax  
# Males_FAL
m1 <- rma.mv(yi_Tmax_FAL, vi, data = data_m, random = RE, R = phylocor)
summary(m1)
m1Tmax_FAL <- data.frame(ES=as.numeric(m1$b), CIL=as.numeric(m1$ci.lb), CUL=as.numeric(m1$ci.ub), Variable = "Tmax")

round(transf.ztor(m1$beta) , digits=3) # estimate
round(transf.ztor(m1$ci.lb), digits=3) # ci lower
round(transf.ztor(m1$ci.ub), digits=3) # ci upper
round(m1$QE, digits=3) #pval
round(m1$QE, digits=1) # QE
m1$QEp # p val for QE

data_mTmax_FAL <- data_m %>%
  dplyr::select(-yi_Tmin_FAL, -yi_SP_FAL, -yi_Tmin_BM, -yi_SP_BM, -yi_Tmax_BM) %>%
  dplyr::rename(yi_s = yi_Tmax_FAL)
data_mTmax_FAL$Variable <- "Tmax"

write.csv(m1Tmax_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_Tmax_FAL.csv", row.names = FALSE)
write.csv(data_mTmax_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_data_m_Tmax_FAL.csv", row.names = FALSE)

# Males_BM
m1 <- rma.mv(yi_Tmax_BM, vi, data = data_m, random = RE, R = phylocor)
summary(m1)
m1Tmax_BM <- data.frame(ES = as.numeric(m1$b), CIL = as.numeric(m1$ci.lb), CUL = as.numeric(m1$ci.ub), Variable = "Tmax")

round(transf.ztor(m1$beta) , digits=3) # estimate
round(transf.ztor(m1$ci.lb), digits=3) # ci lower
round(transf.ztor(m1$ci.ub), digits=3) # ci upper
round(m1$QE, digits=3) #pval
round(m1$QE, digits=1) # QE
m1$QEp # p val for QE

data_mTmax_BM <- data_m %>%
  dplyr::select(-yi_Tmin_FAL, -yi_SP_FAL, -yi_Tmin_BM, -yi_SP_BM, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_Tmax_BM)
data_mTmax_BM$Variable <- "Tmax"

write.csv(m1Tmax_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_Tmax_BM.csv", row.names = FALSE)
write.csv(data_mTmax_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_data_m_Tmax_BM.csv", row.names = FALSE)
# rm(list = c("data_mTmax_FAL", "m1", "m1T", "m1Tmax_FAL"))

#females_FAL
f1 <- rma.mv(yi_Tmax_FAL, vi, data = data_f, random = RE, R = phylocor)
summary(f1)
f1Tmax_FAL <- data.frame(ES=as.numeric(f1$b), CIL=as.numeric(f1$ci.lb), CUL=as.numeric(f1$ci.ub), Variable = "Tmax")

round(transf.ztor(f1$beta) , digits=3) # estimate
round(transf.ztor(f1$ci.lb), digits=3) # ci lower
round(transf.ztor(f1$ci.ub), digits=3) # ci upper
round(f1$QE, digits=3) #pval
round(f1$QE, digits=1) # QE
f1$QEp # p val for QE

data_fTmax_FAL <- data_f %>%
  dplyr::select(-yi_Tmin_FAL, -yi_SP_FAL, -yi_Tmin_BM, -yi_SP_BM, -yi_Tmax_BM) %>%
  dplyr::rename(yi_s = yi_Tmax_FAL)
data_fTmax_FAL$Variable <- "Tmax"

write.csv(f1Tmax_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_Tmax_FAL.csv", row.names = FALSE)
write.csv(data_fTmax_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_data_f_Tmax_FAL.csv", row.names = FALSE)

#females_BM
f1 <- rma.mv(yi_Tmax_BM, vi, data = data_f, random = RE, R = phylocor)
summary(f1)
f1Tmax_BM <- data.frame(ES=as.numeric(f1$b), CIL=as.numeric(f1$ci.lb), CUL=as.numeric(f1$ci.ub), Variable = "Tmax")

round(transf.ztor(f1$beta) , digits=3) # estimate
round(transf.ztor(f1$ci.lb), digits=3) # ci lower
round(transf.ztor(f1$ci.ub), digits=3) # ci upper
round(f1$QE, digits=3) #pval
round(f1$QE, digits=1) # QE
f1$QEp # p val for QE

data_fTmax_BM <- data_f %>%
  dplyr::select(-yi_Tmin_FAL, -yi_SP_FAL, -yi_Tmin_BM, -yi_SP_BM, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_Tmax_BM)
data_fTmax_BM$Variable <- "Tmax"

write.csv(f1Tmax_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_Tmax_BM.csv", row.names = FALSE)
write.csv(data_fTmax_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_data_f_Tmax_BM.csv", row.names = FALSE)


# Males_FAL
m2 <- rma.mv(yi_Tmin_FAL, vi, data = data_m, random = RE, R = phylocor)
summary(m2)
m2Tmin_FAL <- data.frame(ES=as.numeric(m2$b), CIL=as.numeric(m2$ci.lb), CUL=as.numeric(m2$ci.ub), Variable = "Tmin")

round(transf.ztor(m2$beta) , digits=3) # estimate
round(transf.ztor(m2$ci.lb), digits=3) # ci lower
round(transf.ztor(m2$ci.ub), digits=3) # ci upper
round(m2$QE, digits=3) #pval
round(m2$QE, digits=1) # QE
m2$QEp # p val for QE

data_mTmin <- data_m %>%
  dplyr::select(-yi_Tmax_BM, -yi_SP_FAL, -yi_Tmin_BM, -yi_SP_BM, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_Tmin_FAL)
data_mTmin$Variable <- "Tmin"

write.csv(m2Tmin_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_Tmin_FAL.csv", row.names = FALSE)
write.csv(data_mTmin, "Results/Spatial/Spatial_pairwiseFAL_BM_data_m_Tmin_FAL.csv", row.names = FALSE)
# rm(list = c("data_mTmax", "m2", "m2T", "m2Tmax"))

# Males_BM
m2 <- rma.mv(yi_Tmin_BM, vi, data = data_m, random = RE, R = phylocor)
summary(m2)
m2Tmin_BM <- data.frame(ES=as.numeric(m2$b), CIL=as.numeric(m2$ci.lb), CUL=as.numeric(m2$ci.ub), Variable = "Tmin")

round(transf.ztor(m2$beta) , digits=3) # estimate
round(transf.ztor(m2$ci.lb), digits=3) # ci lower
round(transf.ztor(m2$ci.ub), digits=3) # ci upper
round(m2$QE, digits=3) #pval
round(m2$QE, digits=1) # QE
m2$QEp # p val for QE

data_mTmin <- data_m %>%
  dplyr::select(-yi_Tmin_FAL, -yi_SP_FAL, -yi_Tmax_BM, -yi_SP_BM, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_Tmin_BM)
data_mTmin$Variable <- "Tmin"

write.csv(m2Tmin_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_Tmin_BM.csv", row.names = FALSE)
write.csv(data_mTmin, "Results/Spatial/Spatial_pairwiseFAL_BM_data_m_Tmin_BM.csv", row.names = FALSE)
# rm(list = c("data_mTmax", "m2", "m2T", "m2Tmax"))

#females_FAL
f2 <- rma.mv(yi_Tmin_FAL, vi, data = data_f, random = RE, R = phylocor)
summary(f2)
f1Tmin_FAL <- data.frame(ES=as.numeric(f2$b), CIL=as.numeric(f2$ci.lb), CUL=as.numeric(f2$ci.ub), Variable = "Tmin")

round(transf.ztor(f2$beta) , digits=3) # estimate
round(transf.ztor(f2$ci.lb), digits=3) # ci lower
round(transf.ztor(f2$ci.ub), digits=3) # ci upper
round(f2$QE, digits=3) #pval
round(f2$QE, digits=1) # QE
f2$QEp # p val for QE

data_fTmin <- data_f %>%
  dplyr::select(-yi_Tmax_BM, -yi_SP_FAL, -yi_Tmin_BM, -yi_SP_BM, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_Tmin_FAL)
data_fTmin$Variable <- "Tmin"

write.csv(f1Tmin_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_Tmin_FAL.csv", row.names = FALSE)
write.csv(data_fTmin, "Results/Spatial/Spatial_pairwiseFAL_BM_data_f_Tmin_FAL.csv", row.names = FALSE)

#females_BM
f2 <- rma.mv(yi_Tmin_BM, vi, data = data_f, random = RE, R = phylocor)
summary(f2)
f1Tmin_BM <- data.frame(ES=as.numeric(f2$b), CIL=as.numeric(f2$ci.lb), CUL=as.numeric(f2$ci.ub), Variable = "Tmin")

round(transf.ztor(f2$beta) , digits=3) # estimate
round(transf.ztor(f2$ci.lb), digits=3) # ci lower
round(transf.ztor(f2$ci.ub), digits=3) # ci upper
round(f2$QE, digits=3) #pval
round(f2$QE, digits=1) # QE
f2$QEp # p val for QE

data_fTmin <- data_f %>%
  dplyr::select(-yi_Tmax_BM, -yi_SP_FAL, -yi_Tmin_FAL, -yi_SP_BM, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_Tmin_BM)
data_fTmin$Variable <- "Tmin"

write.csv(f1Tmin_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_Tmin_BM.csv", row.names = FALSE)
write.csv(data_fTmin, "Results/Spatial/Spatial_pairwiseFAL_BM_data_f_Tmin_BM.csv", row.names = FALSE)


# Males_FAL
m3 <- rma.mv(yi_SP_FAL, vi, data = data_m, random = RE, R = phylocor)
summary(m3)
m3SP_FAL <- data.frame(ES=as.numeric(m3$b), CIL=as.numeric(m3$ci.lb), CUL=as.numeric(m3$ci.ub), Variable = "S. Precipitation")

round(transf.ztor(m3$beta) , digits=3) # estimate
round(transf.ztor(m3$ci.lb), digits=3) # ci lower
round(transf.ztor(m3$ci.ub), digits=3) # ci upper
round(m3$QE, digits=3) #pval
round(m3$QE, digits=1) # QE
m3$QEp # p val for QE

data_mSP <- data_m %>%
  dplyr::select(-yi_Tmax_BM, -yi_Tmin_BM, -yi_Tmin_FAL, -yi_SP_BM, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_SP_FAL)
data_mSP$Variable <- "S. Precipitation"

write.csv(m3SP_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_SP_FAL.csv", row.names = FALSE)
write.csv(data_mSP, "Results/Spatial/Spatial_pairwiseFAL_BM_data_m_SP_FAL.csv", row.names = FALSE)
rm(list = c("data_mSP", "m3", "m3SP_FAL"))

# Males_BM
m3 <- rma.mv(yi_SP_BM, vi, data = data_m, random = RE, R = phylocor)
summary(m3)
m3SP_BM <- data.frame(ES=as.numeric(m3$b), CIL=as.numeric(m3$ci.lb), CUL=as.numeric(m3$ci.ub), Variable = "S. Precipitation")

round(transf.ztor(m3$beta) , digits=3) # estimate
round(transf.ztor(m3$ci.lb), digits=3) # ci lower
round(transf.ztor(m3$ci.ub), digits=3) # ci upper
round(m3$QE, digits=3) #pval
round(m3$QE, digits=1) # QE
m3$QEp # p val for QE

data_mSP <- data_m %>%
  dplyr::select(-yi_Tmax_BM, -yi_Tmin_BM, -yi_Tmin_FAL, -yi_SP_FAL, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_SP_BM)
data_mSP$Variable <- "S. Precipitation"

write.csv(m3SP_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_SP_BM.csv", row.names = FALSE)
write.csv(data_mSP, "Results/Spatial/Spatial_pairwiseFAL_BM_data_m_SP_BM.csv", row.names = FALSE)

#females_FAL
f3 <- rma.mv(yi_SP_FAL, vi, data = data_f, random = RE, R = phylocor)
summary(f3)

f3SP_FAL <- data.frame(ES=as.numeric(f3$b), CIL=as.numeric(f3$ci.lb), CUL=as.numeric(f3$ci.ub), Variable = "S. Precipitation")

round(transf.ztor(f3$beta) , digits=3) # estimate
round(transf.ztor(f3$ci.lb), digits=3) # ci lower
round(transf.ztor(f3$ci.ub), digits=3) # ci upper
round(f3$QE, digits=3) #pval
round(f3$QE, digits=1) # QE
f3$QEp # p val for QE

data_fSP <- data_f %>%
  dplyr::select(-yi_Tmax_BM, -yi_Tmin_BM, -yi_Tmin_FAL, -yi_SP_BM, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_SP_FAL)
data_fSP$Variable <- "S. Precipitation"

write.csv(f3SP_FAL, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_SP_FAL.csv", row.names = FALSE)
write.csv(data_fSP, "Results/Spatial/Spatial_pairwiseFAL_BM_data_f_SP_FAL.csv", row.names = FALSE)

#females_BM
f3 <- rma.mv(yi_SP_BM, vi, data = data_f, random = RE, R = phylocor)
summary(f3)
f3SP_BM <- data.frame(ES=as.numeric(f3$b), CIL=as.numeric(f3$ci.lb), CUL=as.numeric(f3$ci.ub), Variable = "S. Precipitation")

round(transf.ztor(f3$beta) , digits=3) # estimate
round(transf.ztor(f3$ci.lb), digits=3) # ci lower
round(transf.ztor(f3$ci.ub), digits=3) # ci upper
round(f3$QE, digits=3) #pval
round(f3$QE, digits=1) # QE
f3$QEp # p val for QE

data_fSP <- data_f %>%
  dplyr::select(-yi_Tmax_BM, -yi_Tmin_BM, -yi_Tmin_FAL, -yi_SP_FAL, -yi_Tmax_FAL) %>%
  dplyr::rename(yi_s = yi_SP_BM)
data_fSP$Variable <- "S. Precipitation"

write.csv(f3SP_BM, "Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_SP_BM.csv", row.names = FALSE)
write.csv(data_fSP, "Results/Spatial/Spatial_pairwiseFAL_BM_data_f_SP_BM.csv", row.names = FALSE)

