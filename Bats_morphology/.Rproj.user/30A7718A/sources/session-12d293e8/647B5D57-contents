# Script: Spatial model figures
# Date: 11th November, 2024


##With this script we create the figure 4

# load libraries
library(metafor)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readxl)
library(readr)

# Clear environment
rm(list=ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Meta regression------------------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------
## FAL dataset
#----------------------
# Clear environment
rm(list=ls())
# This script generates the figures related to the meta-regression analysis. 

# Load the dataset 
data_m <- read.csv("Data/Temporal/FAL_ph.csv")

# Load the meta-regression results for males.
df_m <- read.csv("Results/Temporal/Temporal_FAL_meta-regression.csv")

# Create a sequence of log-transformed forearm lengths for plotting the regression line.
temp <- log10(seq(from = min(data_m$fal), to = max(data_m$fal), length.out = 1000))

# Calculate the size of points based on the inverse of the square root of the sampling variance,
# which adjusts the point size to reflect the precision of each data point.
wi <- 1/sqrt(data_m$vi)
size <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

# Set up the range for the x-axis with a small padding.
x_range <- range(log10(data_m$fal))
x_padding <- 0.02

# Generate the plot for males.
M <- ggplot() +
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray") +
  geom_line(data = df_m, aes(temp, cor), color="grey30", size = .6) +
  geom_point(data_m, mapping = aes(x = log10(fal), y = r), fill = "#ffa617", color='black', 
             shape=21, size = size * 0.7, alpha = .4) +
  scale_size(range=c(2, 10))+
  theme_bw(base_size=18) +
  geom_ribbon(data = df_m, aes(x = temp, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "#ffa617", alpha = .1) +
  xlab("log10(forearm length (mm))") +
  ylab("Spearman r (FAL-Year)") +
  ggtitle("FAL") +
  theme(legend.position='none',
        axis.text = element_text(size = 10),        # Match axis text size
        axis.title = element_text(size = 11),       # Match axis label size
        plot.title = element_text(size = 13)) +   
  scale_x_continuous(limits = c(min(x_range) - x_padding, max(x_range) + x_padding)) +
  scale_y_continuous(breaks = seq(-1.0,1.0,0.2), limits= c(-1.0,1.0))


#--------------------------
## Pairwise FAL-BM dataset 
#--------------------------
# The following sections repeat the same analysis for subsets of the FAL data and for body mass.
# FAL
data_mfa <- read.csv("Data/Temporal/pairwiseFAL_BM_ph.csv")
df_mfa <- read.csv("Results/Temporal/Temporal_pairwiseFAL_BM_meta-regression_FAL.csv")
tempfa <- log10(seq(from = min(data_mfa$fal), to = max(data_mfa$fal), length.out = 1000))

wifa <- 1/sqrt(data_mfa$vi)
sizefa <- 2 + 20.0 * (wifa - min(wifa))/(max(wifa) - min(wifa))

x_range <- range(log10(data_mfa$fal))
x_padding <- 0.02

Mfa <- ggplot() +
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray") +
  geom_line(data = df_mfa, aes(tempfa, cor), color="grey30", size = .6) +
  geom_point(data_mfa, mapping = aes(x = log10(fal), y = r_FAL), fill = "#ffa617", color='black', 
             shape=21, size = sizefa * 0.7, alpha = .4) +
  scale_size(range=c(2, 10)) +
  theme_bw(base_size=18) +
  geom_ribbon(data = df_mfa, aes(x = tempfa, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "#ffa617", alpha = .1) +
  xlab("log10(forearm length (mm))") +
  ylab("Spearman r (FAL-Year)") +
  ggtitle("Pairwise FAL") +
  theme(legend.position='none',
        axis.text = element_text(size = 10),        # Match axis text size
        axis.title = element_text(size = 11),       # Match axis label size
        plot.title = element_text(size = 13)) +   
  scale_x_continuous(limits = c(min(x_range) - x_padding, max(x_range) + x_padding)) +
  scale_y_continuous(breaks = seq(-1.0,1.0,0.2), limits= c(-1.0,1.0))

# BM
data_bm <- read.csv("Data/Temporal/pairwiseFAL_BM_ph.csv")
df_bm <- read.csv("Results/Temporal/Temporal_pairwiseFAL_BM_meta-regression_BM.csv")
tempbm <- log10(seq(from = min(data_bm$bm), to = max(data_bm$bm), length.out = 1000))

wibm <- 1/sqrt(data_bm$vi)
sizebm <- 2 + 20.0 * (wibm - min(wibm))/(max(wibm) - min(wibm))

x_range <- range(log10(data_bm$bm))
x_padding <- 0.02

Mbm <- ggplot() +
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray") +
  geom_line(data = df_bm, aes(tempbm, cor), color="grey30", size = .6) +
  geom_point(data_bm, mapping = aes(x = log10(bm), y = r_BM), fill = "#ffa617", color='black', 
             shape=21, size = sizebm * 0.7, alpha = .4) +
  scale_size(range=c(2, 10)) +
  theme_bw(base_size=18) +
  geom_ribbon(data = df_bm, aes(x = tempbm, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "#ffa617", alpha = .1) +
  xlab("log10(body mass (g))") +
  ylab("Spearman r (BM-Year)") +
  ggtitle("Pairwise BM") +
  theme(legend.position='none',
        axis.text = element_text(size = 10),        # Match axis text size
        axis.title = element_text(size = 11),       # Match axis label size
        plot.title = element_text(size = 13)) +   
  scale_x_continuous(limits = c(min(x_range) - x_padding, max(x_range) + x_padding)) +
  scale_y_continuous(breaks = seq(-1.0,1.0,0.2), limits= c(-1.0,1.0))


## Final image (Figure 4) ----
final <- ggarrange(M, Mfa, Mbm, ncol = 3, nrow = 1, align = "h", labels = c("a", "b", "c"))
final
# ggsave(filename='Figures/Figure_4.png', 
#        width = 430, height = 150, units = 'mm', dpi = 900)
ggsave(filename='Figures/Figure_4.pdf',
       width = 300, height = 100, units = 'mm', dpi = 900)
