# Script: Spatial model figures
# Date: 11th November, 2024


##With this script we create the figures 2 and 3 

# load libraries
library(metafor)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readxl)
library(readr)
library(sf)
library(raster)
# Clear environment
rm(list=ls())
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Meta analysis------------------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------
## FAL dataset 
#----------------------
## Males
#Loading male model parameter data
m1 <- read.csv("Results/Spatial/Spatial_FAL_meta-analysis_m_Tmax.csv") # Load meta-analysis results for Tmax
m2 <- read.csv("Results/Spatial/Spatial_FAL_meta-analysis_m_Tmin.csv") # Load meta-analysis results for Tmin
m3 <- read.csv("Results/Spatial/Spatial_FAL_meta-analysis_m_SP.csv")   # Load meta-analysis results for Summer precipitation
ESMActi <- rbind(m1,m2,m3)

# Load male data for the different variables
d1 <- read.csv("Results/Spatial/Spatial_FAL_data_m_Tmax.csv")
d2 <- read.csv("Results/Spatial/Spatial_FAL_data_m_Tmin.csv")
d3 <- read.csv("Results/Spatial/Spatial_FAL_data_m_SP.csv")
data <- rbind(d1,d2,d3)

# Calculate the correlation coefficient r from Fisher's Z-transformed correlation (yi)
data$r <- tanh(data$yi)

# Calculate the weights (wi) based on the inverse of the square root of the variance (vi)
# These weights will be used to adjust the size of the points in the plot
wi    <- 1/sqrt(data$vi)    
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi)) 

#Plot
# Create a ggplot to visualize the effect size (Spearman r) for the different environmental variables
sm <- ggplot(data, aes(x = Variable, y = r, fill = Variable, size = size))+
  geom_hline(yintercept= 0, size=1, linetype="dashed", col="dark gray")+
  geom_jitter(alpha=0.4, width=0.15, height=0, shape=21, col="black")+
  scale_fill_manual(values = alpha(c("#20a380","#ff5729","#3186d6")))+
  scale_size(range=c(.1, 8))+
  geom_pointrange(data = ESMActi, aes(x = Variable, y = ES), ymin=ESMActi$CIL, 
                  ymax=ESMActi$CUL, shape= 21, size=0.7, lwd=0.8,  fill='#CB2A04FF', col="black")+
  xlab("")+
  ylab("Spearman r (FAL)")+
  theme_bw()+ 
  theme(legend.position='none')+
  theme(axis.text = element_text(size = 10))+
  ggtitle("FAL Males   n=766") 


## females 
#Loading female model parameter data
f1ff <- read.csv("Results/Spatial/Spatial_FAL_meta-analysis_f_Tmax.csv")
f2ff <- read.csv("Results/Spatial/Spatial_FAL_meta-analysis_f_Tmin.csv")
f3ff <- read.csv("Results/Spatial/Spatial_FAL_meta-analysis_f_SP.csv")
ESMActi1ff <- rbind(f1ff,f2ff,f3ff)

# Load female data for the different variables
e1ff <- read.csv("Results/Spatial/Spatial_FAL_data_f_Tmax.csv")
e2ff <- read.csv("Results/Spatial/Spatial_FAL_data_f_Tmin.csv")
e3ff <- read.csv("Results/Spatial/Spatial_FAL_data_f_SP.csv")
data1ff <- rbind(e1ff,e2ff,e3ff)

# Calculate the correlation coefficient r from Fisher's Z-transformed correlation (yi)
data1ff$r <- tanh(data1ff$yi_s)

# Calculate the weights (wi) based on the inverse of the square root of the variance (vi)
# These weights will be used to adjust the size of the points in the plot
wi1ff    <- 1/sqrt(data1ff$vi)   
size1ff  <- 2 + 20.0 * (wi1ff - min(wi1ff))/(max(wi1ff) - min(wi1ff)) 

#Plot
# Create a ggplot to visualize the effect size (Spearman r) for the different environmental variables
sff <- ggplot(data1ff, aes(x = Variable, y = r, fill = Variable, size = size1ff))+
  geom_hline(yintercept= 0, size=1, linetype="dashed", col="dark gray")+
  geom_jitter(alpha=0.4, width=0.15, height=0, shape=21, col="black")+
  scale_fill_manual(values = alpha(c("#20a380","#ff5729","#3186d6")))+
  scale_size(range=c(.1, 8))+
  geom_pointrange(data = ESMActi1ff, aes(x = Variable, y = ES), ymin=ESMActi1ff$CIL, 
                  ymax=ESMActi1ff$CUL, shape= 21, size=0.7, lwd=0.8, fill='#CB2A04FF', col="black")+
  xlab("")+
  ylab("Spearman r (FAL)")+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text = element_text(size = 10))+
  ggtitle("FAL Females   n=740")

# ## Save figure
# fig_spearman <- ggarrange(sf, sm, ncol = 2, nrow = 1, align = "h")
# fig_spearman

#ggsave(filename='Results/3_Spearman_FA2.png', 
    #   width = 250, height = 130, units = 'mm', dpi = 600)

#--------------------------------------
## Pairwise FAL-BM dataset (Allometry) 
#--------------------------------------
## Tmax
# Males
data_mfa <- subset(read.csv("Data/Spatial/pairwiseFAL_BM_ph.csv"), Sex == "male")
tempfa <- seq(from = min(data_mfa$yi_Tmax_BM), to = max(data_mfa$yi_Tmax_BM), length.out = 1000)

wifa <- 1/sqrt(data_mfa$vi)
sizefa <- 2 + 20.0 * (wifa - min(wifa))/(max(wifa) - min(wifa))

x_range <- range(data_mfa$yi_Tmax_BM)
x_padding <- 0.02

m7 <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_Tmax_FAL.csv")
m9 <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_Tmax_BM.csv")
m7 <- m7 %>%
  dplyr::rename(CILFAL = CIL,
         CULFAL = CUL,
         ESFAL = ES,
         VariableFAL = Variable)
a <- cbind(m7,m9)

data_mfa$color_group <- with(data_mfa,
                             ifelse(r_Tmax_BM == r_Tmax_FAL, "Red Border", 
                             ifelse(r_Tmax_BM < 0 & r_Tmax_FAL < r_Tmax_BM, "Color 1", 
                             ifelse(r_Tmax_BM < 0 & r_Tmax_FAL > r_Tmax_BM, "Color 2", 
                             ifelse(r_Tmax_BM > 0 & r_Tmax_FAL > r_Tmax_BM, "Color 3", "Other")))))

Mfa <- ggplot() +
  annotate("polygon", x = c(-Inf, 0, 0), y = c(-Inf, -Inf, 0), 
           fill = "#98C2EA", alpha = 0.2) +  
  annotate("polygon", x = c(0, 0, Inf), y = c(0, Inf, Inf), #FFE23D#F8CE8A
           fill = "#F8CE8A", alpha = 0.2) +
  annotate("polygon", x = c(-Inf, -Inf, 0), y = c(-Inf, 0, 0), #A6DAD8
           fill = "#79C8B3", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = 0, fill = "#79C8B3", alpha = 0.2)+
  geom_hline(yintercept = 0, size = .7, linetype = "dashed", color = "dark gray") +
  geom_vline(xintercept = 0, size = .7, linetype = "dashed", color = "dark gray") +
  geom_abline(slope = 1, intercept = 0, color = "tomato 4", linetype = "dashed", size = .7) +
  # geom_line(data = df_mfa, aes(tempfa, cor), color="grey30", size = .4) +
  # geom_ribbon(data = df_mfa, aes(x = tempfa, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "grey20", alpha = .1) +
  geom_point(data_mfa, mapping = aes(x = r_Tmax_BM, y = r_Tmax_FAL), 
             shape = 21, size = sizefa * 0.5, alpha = 1) +
  scale_size(range = c(2, 10)) +
  #scale_fill_manual(values = c("Color 1" = "#E67B7B", "Color 2" = "#FFC552", "Color 3" = "#315A83", "Other" = "grey90")) +
  theme_bw(base_size = 18) +
  xlab("Spearman r (BM-Tmax)") +
  ylab("Spearman r (FAL-Tmax)") +
  #ggtitle("Males") +
  theme(legend.position='none') +
  scale_x_continuous(breaks = seq(-0.6,0.6,0.2), limits= c(-0.7,0.7), 
                     labels = scales::label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(-0.6,0.6,0.2), limits= c(-0.7,0.7), 
                     labels = scales::label_number(accuracy = 0.1))+
  geom_pointrange(data = a, aes(x = ES, y = ESFAL), 
                  ymin = a$CILFAL, ymax = a$CULFAL,
                  shape = 21, size = 0.7, lwd = 0.8, fill = '#CB2A04FF', col = "black")+
  geom_pointrange(data = a, aes(x = ES, y = ESFAL), 
                  xmin = a$CIL, xmax = a$CUL,
                  shape = 21, size = 0.7, lwd = 0.8, fill = '#CB2A04FF', col = "black")+
  theme(
    axis.text = element_text(size = 10),        # Match axis text size
    axis.title = element_text(size = 11),       # Match axis label size
    plot.title = element_text(size = 13)        # Match title size
  ) +
  ggtitle("Allometry BM-FAL/Tmax Males n=567")+
  theme(panel.border = element_rect(colour = "black", size = 0.5))


# Females
data_ffa <- subset(read.csv("Data/Spatial/pairwiseFAL_BM_ph.csv"), Sex == "female")
tempfaf <- seq(from = min(data_ffa$yi_Tmax_BM), to = max(data_ffa$yi_Tmax_BM), length.out = 1000)

wifaf <- 1/sqrt(data_ffa$vi)
sizefaf <- 2 + 20.0 * (wifaf - min(wifaf))/(max(wifaf) - min(wifaf))

x_range <- range(data_ffa$yi_Tmax_BM)
x_padding <- 0.02

m7fa <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_Tmax_FAL.csv")
m9fa <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_Tmax_BM.csv")
m7fa <- m7fa %>%
  dplyr::rename(CILFAL = CIL,
         CULFAL = CUL,
         ESFAL = ES,
         VariableFAL = Variable)
afa <- cbind(m7fa,m9fa)

data_ffa$color_group <- with(data_ffa,
                             ifelse(r_Tmax_BM == r_Tmax_FAL, "Red Border", 
                             ifelse(r_Tmax_BM < 0 & r_Tmax_FAL < r_Tmax_BM, "Color 1", 
                             ifelse(r_Tmax_BM < 0 & r_Tmax_FAL > r_Tmax_BM, "Color 2", 
                             ifelse(r_Tmax_BM > 0 & r_Tmax_FAL > r_Tmax_BM, "Color 3", "Other")))))
Ffa <- ggplot() +
  annotate("polygon", x = c(-Inf, 0, 0), y = c(-Inf, -Inf, 0), 
           fill = "#98C2EA", alpha = 0.2) +  
  annotate("polygon", x = c(0, 0, Inf), y = c(0, Inf, Inf), #FFE23D#F8CE8A
           fill = "#F8CE8A", alpha = 0.2) +
  annotate("polygon", x = c(-Inf, -Inf, 0), y = c(-Inf, 0, 0), #A6DAD8
           fill = "#79C8B3", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = 0, fill = "#79C8B3", alpha = 0.2)+
  geom_hline(yintercept = 0, size = .7, linetype = "dashed", color = "dark gray") +
  geom_vline(xintercept= 0, size = .7, linetype = "dashed", color = "dark gray") +
  geom_abline(slope = 1, intercept = 0, color = "tomato 4", linetype = "dashed", size = .7) +
  #geom_ribbon(data = df_ffaf, aes(x = tempfaf, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "grey30", alpha = .1) +
  #geom_line(data = df_ffaf, aes(tempfaf, cor), color="grey30", size = .4) +
  geom_point(data_ffa, mapping = aes(x = r_Tmax_BM, y = r_Tmax_FAL),
             shape=21, size = sizefaf * 0.5, alpha = 1) +
  #scale_size(range=c(2, 10)) +
  #scale_fill_manual(values = c("Color 1" = "#E67B7B", "Color 2" = "#FFC552", "Color 3" = "#315A83", "Other" = "grey90")) +
  theme_bw(base_size=18) +
  xlab("Spearman r (BM-Tmax)") +
  ylab("Spearman r (FAL-Tmax)") +
  #ggtitle("Females") +
  theme(legend.position='none') +
  scale_x_continuous(breaks = seq(-0.8,0.8,0.2), limits= c(-0.9,0.9), 
                     labels = scales::label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(-0.8,0.8,0.2), limits= c(-0.9,0.9), 
                     labels = scales::label_number(accuracy = 0.1))+
  geom_pointrange(data = afa, aes(x = ES, y = ESFAL), 
                  ymin = afa$CILFAL, ymax = afa$CULFAL,
                  shape = 21, size = 0.7, lwd = 0.8, fill = '#CB2A04FF', col = "black")+
  geom_pointrange(data = afa, aes(x = ES, y = ESFAL), 
                  xmin = afa$CIL, xmax = afa$CUL,
                  shape = 21, size = 0.7, lwd = 0.8, fill = '#CB2A04FF', col = "black")+
  theme(
    axis.text = element_text(size = 10),        # Match axis text size
    axis.title = element_text(size = 11),       # Match axis label size
    plot.title = element_text(size = 13)        # Match title size
  ) +
  ggtitle("Allometry BM-FAL/Tmax Females n=566")+
  theme(panel.border = element_rect(colour = "black", size = 0.5))


## Tmin
# Males
data_mfa <- subset(read.csv("Data/Spatial/pairwiseFAL_BM_ph.csv"), Sex == "male")
tempfa <- seq(from = min(data_mfa$yi_Tmin_BM), to = max(data_mfa$yi_Tmin_BM), length.out = 1000)

wifa <- 1/sqrt(data_mfa$vi)
sizefa <- 2 + 20.0 * (wifa - min(wifa))/(max(wifa) - min(wifa))

x_range <- range(data_mfa$yi_Tmax_BM)
x_padding <- 0.02

m7 <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_Tmin_FAL.csv")
m9 <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_m_Tmin_BM.csv")
m7 <- m7 %>%
  dplyr::rename(CILFAL = CIL,
         CULFAL = CUL,
         ESFAL = ES,
         VariableFAL = Variable)
a <- cbind(m7,m9)

data_mfa$color_group <- with(data_mfa,
                             ifelse(r_Tmin_BM == r_Tmin_FAL, "Red Border", 
                             ifelse(r_Tmin_BM < 0 & r_Tmin_FAL < r_Tmin_BM, "Color 1", 
                             ifelse(r_Tmin_BM < 0 & r_Tmin_FAL > r_Tmin_BM, "Color 2", 
                             ifelse(r_Tmin_BM > 0 & r_Tmin_FAL > r_Tmin_BM, "Color 3", "Other")))))

Mfa1 <- ggplot() +
  annotate("polygon", x = c(-Inf, 0, 0), y = c(-Inf, -Inf, 0), 
           fill = "#98C2EA", alpha = 0.2) +  
  annotate("polygon", x = c(0, 0, Inf), y = c(0, Inf, Inf), #FFE23D#F8CE8A
           fill = "#F8CE8A", alpha = 0.2) +
  annotate("polygon", x = c(-Inf, -Inf, 0), y = c(-Inf, 0, 0), #A6DAD8
           fill = "#79C8B3", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = 0, fill = "#79C8B3", alpha = 0.2)+
  geom_hline(yintercept = 0, size = .7, linetype = "dashed", color = "dark gray") +
  geom_vline(xintercept = 0, size = .7, linetype = "dashed", color = "dark gray") +
  geom_abline(slope = 1, intercept = 0, color = "tomato 4", linetype = "dashed", size = .7) +
  # geom_line(data = df_mfa, aes(tempfa, cor), color="grey30", size = .4) +
  # geom_ribbon(data = df_mfa, aes(x = tempfa, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "grey20", alpha = .1) +
  geom_point(data_mfa, mapping = aes(x = r_Tmin_BM, y = r_Tmin_FAL), 
             shape = 21, size = sizefa * 0.5, alpha = 1) +
  scale_size(range = c(2, 10)) +
  #scale_fill_manual(values = c("Color 1" = "#E67B7B", "Color 2" = "#FFC552", "Color 3" = "#315A83", "Other" = "grey90")) +
  theme_bw(base_size = 18) +
  xlab("Spearman r (BM-Tmin)") +
  ylab("Spearman r (FAL-Tmin)") +
  #ggtitle("Males") +
  theme(legend.position='none') +
  scale_x_continuous(breaks = seq(-0.6,0.6,0.2), limits= c(-0.7,0.7), 
                     labels = scales::label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(-0.6,0.6,0.2), limits= c(-0.7,0.7), 
                     labels = scales::label_number(accuracy = 0.1))+
  geom_pointrange(data = a, aes(x = ES, y = ESFAL), 
                  ymin = a$CILFAL, ymax = a$CULFAL,
                  shape = 21, size = 0.7, lwd = 0.8, fill = '#CB2A04FF', col = "black")+
  geom_pointrange(data = a, aes(x = ES, y = ESFAL), 
                  xmin = a$CIL, xmax = a$CUL,
                  shape = 21, size = 0.7, lwd = 0.8, fill = '#CB2A04FF', col = "black")+
  theme(
    axis.text = element_text(size = 10),        # Match axis text size
    axis.title = element_text(size = 11),       # Match axis label size
    plot.title = element_text(size = 13)        # Match title size
  ) +
  ggtitle("Allometry BM-FAL/Tmin Males n=567")+
  theme(panel.border = element_rect(colour = "black", size = 0.5))


# Females-Tmin
data_ffa <- subset(read.csv("Data/Spatial/pairwiseFAL_BM_ph.csv"), Sex == "female")
tempfaf <- seq(from = min(data_ffa$yi_Tmin_BM), to = max(data_ffa$yi_Tmin_BM), length.out = 1000)

wifaf <- 1/sqrt(data_ffa$vi)
sizefaf <- 2 + 20.0 * (wifaf - min(wifaf))/(max(wifaf) - min(wifaf))

x_range <- range(data_ffa$yi_Tmin_BM)
x_padding <- 0.02

m7 <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_Tmin_FAL.csv")
m9 <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-analysis_f_Tmin_BM.csv")
m7 <- m7 %>%
  dplyr::rename(CILFAL = CIL,
         CULFAL = CUL,
         ESFAL = ES,
         VariableFAL = Variable)
a <- cbind(m7,m9)

data_ffa$color_group <- with(data_ffa,
                             ifelse(r_Tmin_BM == r_Tmin_FAL, "Red Border", 
                             ifelse(r_Tmin_BM < 0 & r_Tmin_FAL < r_Tmin_BM, "Color 1", 
                             ifelse(r_Tmin_BM < 0 & r_Tmin_FAL > r_Tmin_BM, "Color 2", 
                             ifelse(r_Tmin_BM > 0 & r_Tmin_FAL > r_Tmin_BM, "Color 3", "Other")))))

Ffa1 <- ggplot() +
  annotate("polygon", x = c(-Inf, 0, 0), y = c(-Inf, -Inf, 0), 
           fill = "#98C2EA", alpha = 0.2) +  
  annotate("polygon", x = c(0, 0, Inf), y = c(0, Inf, Inf), #FFE23D#F8CE8A
           fill = "#F8CE8A", alpha = 0.2) +
  annotate("polygon", x = c(-Inf, -Inf, 0), y = c(-Inf, 0, 0), #A6DAD8
           fill = "#79C8B3", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = 0, fill = "#79C8B3", alpha = 0.2)+
  geom_hline(yintercept = 0, size = .7, linetype = "dashed", color = "dark gray") +
  geom_vline(xintercept= 0, size = .7, linetype = "dashed", color = "dark gray") +
  geom_abline(slope = 1, intercept = 0, color = "tomato 4", linetype = "dashed", size = .7) +
  # geom_ribbon(data = df_ffaf, aes(x = tempfaf, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "grey30", alpha = .1) +
  # geom_line(data = df_ffaf, aes(tempfaf, cor), color="grey30", size = .4) +
  geom_point(data_ffa, mapping = aes(x = r_Tmin_BM, y = r_Tmin_FAL), 
             shape = 21, size = sizefaf * 0.5, alpha = 1) +
  scale_size(range = c(2, 10)) +
  #scale_fill_manual(values = c("Color 1" = "#E67B7B", "Color 2" = "#FFC552", "Color 3" = "#315A83", "Other" = "grey90")) +
  theme_bw(base_size = 18) +
  xlab("Spearman r (BM-Tmin)") +
  ylab("Spearman r (FAL-Tmin)") +
  #ggtitle("Females") +
  theme(legend.position='none') +
  scale_x_continuous(breaks = seq(-0.6,0.6,0.2), limits= c(-0.6,0.6), labels = scales::label_number()) +
  scale_y_continuous(breaks = seq(-0.6,0.6,0.2), limits= c(-0.6,0.6), labels = scales::label_number())+
  geom_pointrange(data = a, aes(x = ES, y = ESFAL), 
                  ymin = a$CILFAL, ymax = a$CULFAL,
                  shape = 21, size = 0.7, lwd = 0.8, fill = '#CB2A04FF', col = "black")+
  geom_pointrange(data = a, aes(x = ES, y = ESFAL), 
                  xmin = a$CIL, xmax = a$CUL,
                  shape = 21, size = 0.7, lwd = 0.8, fill = '#CB2A04FF', col = "black")+
  theme(
    axis.text = element_text(size = 10),        # Match axis text size
    axis.title = element_text(size = 11),       # Match axis label size
    plot.title = element_text(size = 13)        # Match title size
  ) +
  ggtitle("Allometry BM-FAL/Tmin Females n=566") +
  theme(panel.border = element_rect(colour = "black", size = 0.5))

## Final image ----
final <- ggarrange(sff, Ffa, Ffa1, sm, Mfa,  Mfa1, labels = c("a", "b", "c", "d", "e", "f"), align = "h")
final 
# ggsave(filename='Figures/Figure_2.png', 
#        width = 300, height = 200, units = 'mm', dpi = 600)
ggsave(filename='Figures/Figure_2.pdf',
       width = 300, height = 200, units = 'mm', dpi = 600)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Meta regression------------------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clear environment
rm(list=ls())
# This script generates the figures related to the meta-regression analysis. 

#----------------------
## FAL 
#----------------------
# Males
# Load the dataset and filter it to include only male bats.
data_m <- subset(read.csv("Data/Spatial/FAL_ph.csv"), Sex == "male")

# Load the meta-regression results for males.
df_m <- read.csv("Results/Spatial/Spatial_FAL_meta-regression_m_Tmax.csv")

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
  geom_point(data_m, mapping = aes(x = log10(fal), y = r_Tmax), fill = "#ff5729", color='black', 
             shape=21, size = size * 0.5, alpha = .4) +
  scale_size(range=c(2, 10))+
  theme_bw(base_size=18) +
  geom_ribbon(data = df_m, aes(x = temp, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "#ff5729", alpha = .1) +
  xlab("log10(forearm length (mm))") +
  ylab("Spearman r (FAL-Tmax)") +
  ggtitle("FAL Males") +
  theme(legend.position='none',
        axis.text = element_text(size = 10),        # Match axis text size
        axis.title = element_text(size = 11),       # Match axis label size
        plot.title = element_text(size = 13)) +        # Match title size) +
  scale_x_continuous(limits = c(min(x_range) - x_padding, max(x_range) + x_padding)) +
  scale_y_continuous(breaks = seq(-1.0,1.0,0.2), limits= c(-1.0,1.0))

# Females
# Load the dataset and filter it to include only female bats.
data_f <- subset(read.csv("Data/Spatial/FAL_ph.csv"), Sex == "female")

# Load the meta-regression results for females.
df_f <- read.csv("Results/Spatial/Spatial_FAL_meta-regression_f_Tmax.csv")

# Create a sequence of log-transformed forearm lengths for plotting the regression line.
tempf <- log10(seq(from = min(data_f$fal), to = max(data_f$fal), length.out = 1000))

# Calculate the size of points based on the inverse of the square root of the sampling variance.
wif <- 1/sqrt(data_f$vi)
sizef <- 2 + 20.0 * (wif - min(wif))/(max(wif) - min(wif))

# Set up the range for the x-axis with a small padding.
x_range <- range(log10(data_f$fal))
x_padding <- 0.02

# Generate the plot for females.
Ff <- ggplot() +
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray") +
  geom_line(data = df_f, aes(tempf, cor), color="grey30", size = .6) +
  geom_point(data_f, mapping = aes(x = log10(fal), y = r_Tmax), fill = "#ff5729", color='black', 
             shape=21, size = sizef * 0.5, alpha = .4) +
  scale_size(range=c(2, 10)) +
  theme_bw(base_size=18) +
  geom_ribbon(data = df_f, aes(x = tempf, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "#ff5729", alpha = .1) +
  xlab("log10(forearm length (mm))") +
  ylab("Spearman r (FAL-Tmax)") +
  ggtitle("FAL Females") +
  theme(legend.position='none',
        axis.text = element_text(size = 10),        # Match axis text size
        axis.title = element_text(size = 11),       # Match axis label size
        plot.title = element_text(size = 13)) +        # Match title size) +
  scale_x_continuous(limits = c(min(x_range) - x_padding, max(x_range) + x_padding)) +
  scale_y_continuous(breaks = seq(-1.0,1.0,0.2), limits= c(-1.0,1.0))

#--------------------------
## Pairwise FAL-BM dataset 
#--------------------------
# The following sections repeat the same analysis for subsets of the FAL data and for body mass (BM).

# Males for subFAL
data_mfa <- subset(read.csv("Data/Spatial/pairwiseFAL_BM_ph.csv"), Sex == "male")
df_mfa <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-regression_m_Tmax_FAL.csv")
tempfa <- log10(seq(from = min(data_mfa$fal), to = max(data_mfa$fal), length.out = 1000))

wifa <- 1/sqrt(data_mfa$vi)
sizefa <- 2 + 20.0 * (wifa - min(wifa))/(max(wifa) - min(wifa))

x_range <- range(log10(data_mfa$fal))
x_padding <- 0.02

Mfa <- ggplot() +
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray") +
  geom_line(data = df_mfa, aes(tempfa, cor), color="grey30", size = .6) +
  geom_point(data_mfa, mapping = aes(x = log10(fal), y = r_Tmax_FAL), fill = "#ff5729", color='black', 
             shape=21, size = sizefa * 0.5, alpha = .4) +
  scale_size(range=c(2, 10)) +
  theme_bw(base_size=18) +
  geom_ribbon(data = df_mfa, aes(x = tempfa, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "#ff5729", alpha = .1) +
  xlab("log10(forearm length (mm))") +
  ylab("Spearman r (FAL-Tmax)") +
  ggtitle("Pairwise FAL Males")+
  theme(legend.position='none',
        axis.text = element_text(size = 10),        # Match axis text size
        axis.title = element_text(size = 11),       # Match axis label size
        plot.title = element_text(size = 13)) +        # Match title size) +
  scale_x_continuous(limits = c(min(x_range) - x_padding, max(x_range) + x_padding)) +
  scale_y_continuous(breaks = seq(-1.0,1.0,0.2), limits= c(-1.0,1.0))

# Females for subFAL
data_ffa <- subset(read.csv("Data/Spatial/pairwiseFAL_BM_ph.csv"), Sex == "female")
df_ffaf <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-regression_f_Tmax_FAL.csv")
tempfaf <- log10(seq(from = min(data_ffa$fal), to = max(data_ffa$fal), length.out = 1000))

wifaf <- 1/sqrt(data_ffa$vi)
sizefaf <- 2 + 20.0 * (wifaf - min(wifaf))/(max(wifaf) - min(wifaf))

x_range <- range(log10(data_ffa$fal))
x_padding <- 0.02

Ffa <- ggplot() +
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray") +
  geom_line(data = df_ffaf, aes(tempfaf, cor), color="grey30", size = .6) +
  geom_point(data_ffa, mapping = aes(x = log10(fal), y = r_Tmax_FAL), fill = "#ff5729", color='black', 
             shape=21, size = sizefaf * 0.5, alpha = .4) +
  scale_size(range=c(2, 10)) +
  theme_bw(base_size=18) +
  geom_ribbon(data = df_ffaf, aes(x = tempfaf, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "#ff5729", alpha = .1) +
  xlab("log10(forearm length (mm))") +
  ylab("Spearman r (FAL-Tmax)") +
  ggtitle("Pairwise FAL Females") +
  theme(legend.position='none',
        axis.text = element_text(size = 10),        # Match axis text size
        axis.title = element_text(size = 11),       # Match axis label size
        plot.title = element_text(size = 13)) +        # Match title size) +
  scale_x_continuous(limits = c(min(x_range) - x_padding, max(x_range) + x_padding)) +
  scale_y_continuous(breaks = seq(-1.0,1.0,0.2), limits= c(-1.0,1.0))


# Males for BM
data_bm <- subset(read.csv("Data/Spatial/pairwiseFAL_BM_ph.csv"), Sex == "male")
df_bm <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-regression_m_Tmax_BM.csv")
tempbm <- log10(seq(from = min(data_bm$bm), to = max(data_bm$bm), length.out = 1000))

wibm <- 1/sqrt(data_bm$vi)
sizebm <- 2 + 20.0 * (wibm - min(wibm))/(max(wibm) - min(wibm))

x_range <- range(log10(data_bm$bm))
x_padding <- 0.02

Mbm <- ggplot() +
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray") +
  geom_line(data = df_bm, aes(tempbm, cor), color="grey30", size = .6) +
  geom_point(data_bm, mapping = aes(x = log10(bm), y = r_Tmax_BM), fill = "#ff5729", color='black', 
             shape=21, size = sizebm * 0.5, alpha = .4) +
  scale_size(range=c(2, 10)) +
  theme_bw(base_size=18) +
  geom_ribbon(data = df_bm, aes(x = tempbm, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "#ff5729", alpha = .1) +
  xlab("log10(body mass (g))") +
  ylab("Spearman r (BM-Tmax)") +
  ggtitle("Pairwise BM Males") +
  theme(legend.position='none',
        axis.text = element_text(size = 10),        # Match axis text size
        axis.title = element_text(size = 11),       # Match axis label size
        plot.title = element_text(size = 13)) +        # Match title size) +
  scale_x_continuous(limits = c(min(x_range) - x_padding, max(x_range) + x_padding)) +
  scale_y_continuous(breaks = seq(-1.0,1.0,0.2), limits= c(-1.0,1.0))


# Females for BM
data_bf <- subset(read.csv("Data/Spatial/pairwiseFAL_BM_ph.csv"), Sex == "female")
df_bf <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_meta-regression_f_Tmax_BM.csv")
tempbf <- log10(seq(from = min(data_bf$bm), to = max(data_bf$bm), length.out = 1000))

wibf <- 1/sqrt(data_bf$vi)
sizebf <- 2 + 20.0 * (wibf - min(wibf))/(max(wibf) - min(wibf))

x_range <- range(log10(data_bf$bm))
x_padding <- 0.02

Fbm <- ggplot() +
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray") +
  geom_line(data = df_bf, aes(tempbf, cor), color="grey30", size = .6) +
  geom_point(data_bf, mapping = aes(x = log10(bm), y = r_Tmax_BM), fill = "#ff5729", color='black', 
             shape=21, size = sizebf * 0.5, alpha = .4) +
  scale_size(range=c(2, 10)) +
  theme_bw(base_size=18) +
  geom_ribbon(data = df_bf, aes(x = tempbf, ymin = cor_ci.lb, ymax = cor_ci.ub), fill = "#ff5729", alpha = .1) +
  xlab("log10(body mass (g))") +
  ylab("Spearman r (BM-Tmax)") +
  ggtitle("Pairwise BM Females") +
  theme(legend.position='none',
        axis.text = element_text(size = 10),        # Match axis text size
        axis.title = element_text(size = 11),       # Match axis label size
        plot.title = element_text(size = 13)) +        # Match title size) +
  scale_x_continuous(limits = c(min(x_range) - x_padding, max(x_range) + x_padding)) +
  scale_y_continuous(breaks = seq(-1.0,1.0,0.2), limits= c(-1.0,1.0))



## Final image ----
final <- ggarrange(Ff, Ffa, Fbm, M, Mfa, Mbm,labels = c("a", "b", "c", "d", "e", "f"))
final
# ggsave(filename='Figures/Figure_3.png',
#        width = 300, height = 205, units = 'mm', dpi = 600)
ggsave(filename='Figures/Figure_3.pdf',
       width = 300, height = 200, units = 'mm', dpi = 600)


