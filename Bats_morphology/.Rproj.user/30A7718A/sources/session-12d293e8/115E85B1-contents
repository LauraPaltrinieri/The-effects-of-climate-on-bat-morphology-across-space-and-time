# Script: Figure 1
# Date: 11th November, 2024


##With this script we create the figure 1

#library
library(sf)
library(ggpubr)
library(sf)
library(here)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(magrittr)
library(paletteer)
library(scales)
library(RColorBrewer)
# Clear environment
rm(list=ls())

#----------------------
## Spatial Data
#----------------------
b_S2 <- read.csv("Data/Spatial/FAL_data.csv")

# Calculate "n", the number of cells where each species is present
aa <-  b_S2 %>%
  group_by(Bat_species) %>%
  dplyr::summarise(n1 = length(unique(cell_id)))

# Join aa with the b_S dataset to include coordinate variables
b_S3 <- left_join(b_S2, aa, by = c("Bat_species"))

b_S3 <-  b_S3 %>%
  group_by(Bat_species, n1, Lat_c, Long_c)%>%
  dplyr::summarise(n = n())

## Get world polygon for background
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot the map
map.color_S <- ggplot(data = world) +
  theme(plot.margin = unit(c(0.2, 2, 0.2, 0), 'lines'))+
  geom_sf(fill = "white") +
  coord_sf(xlim = c(-15, 45),
           ylim = c(27,65 ),
           expand = FALSE) +
  geom_point(data = b_S3, aes(x = Long_c, y = Lat_c,
                              fill = Bat_species),
             size = 2.5, shape = 21, color = "gray50",  stroke = 0.3) +
  scale_fill_manual(values = alpha(c("#6FB7A7", "#ABD2CB", "#D2E6CD", "#C3E1A8", "#B8D28C",
                                     "#F0F0BF", "#EBD89F", "#EEC9AD", "#E4B888", "#E49A91", 
                                     "#EFAFAF", "#D8B2C4", "#BD98BD", "#C6C4D8", "#ADA8D9",
                                     "#A1C3D9", "#B4CDD9", "#B6AFA4", "#D2D2D2"),0.7))+
  xlab("") +
  ylab("") +
  theme_bw() 

## Display map
map.color_S <- map.color_S + 
  theme(legend.position = "none")

map.color_S <-  map.color_S+
  theme(legend.text = element_text(size = 8)
  )+ 
  guides(color = guide_legend(override.aes = list(size = 1.2), nrow = 3))

map.color_S <-  map.color_S + theme(axis.text.y = element_text(size = 8), 
                                    axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8))


# Create the bar plot 
bar_plot_S <- ggplot(aa) +
  geom_bar(aes(x = n1, y = factor(Bat_species, levels = rev(levels(factor(Bat_species)))), 
               fill = Bat_species), colour = "black", 
           stat = "identity", alpha = 0.7, linewidth = 0.35,
           position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c("#6FB7A7", "#ABD2CB", "#D2E6CD", "#C3E1A8", "#B8D28C",
                               "#F0F0BF", "#EBD89F", "#EEC9AD", "#E4B888", "#E49A91", 
                               "#EFAFAF", "#D8B2C4", "#BD98BD", "#C6C4D8", "#ADA8D9",
                               "#A1C3D9", "#B4CDD9", "#B6AFA4", "#D2D2D2"))+
theme_bw() + 
  scale_x_continuous(breaks = seq(0, 700, by = 50))

bar_plot_S <- bar_plot_S + 
  guides(fill = "none") 

bar_plot_S <-  bar_plot_S + theme(axis.text.y = element_text(size = 12, face = "bold.italic"), 
                                  axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = ))+
  theme(plot.margin = unit(c(0.2, 0, 0.2, 2), 'lines'))+
  labs(x ="n cells")+
  ggtitle("Spatial")+ 
  theme(plot.title = element_text(size = 20, face = "bold"))

#----------------------
## Temporal Data 
#----------------------
b_T <- read.csv("Data/Temporal/FAL_data.csv")

# Calculate "n", the number of years where each species is present
aa <-  b_T %>%
  group_by(Bat_species) %>%
  dplyr::summarise(n1 = length(unique(Year)))

# Join aa with the b_T dataset to include coordinate variables
n2_T <- left_join(b_T, aa, by = c("Bat_species"))
n2_T <-  n2_T %>%
  group_by(Bat_species, cell_id, n1, Lat_c, Long_c)%>%
  dplyr::summarise(n = n())

## Get world polygon for background
world <- ne_countries(scale = "medium", returnclass = "sf")

nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

map.color_T <- ggplot(data = world) +
  theme(plot.margin = unit(c(0.2, 4, 0.2, 0), 'lines'))+
  geom_sf(fill = "white") +
  coord_sf(xlim = c(-15, 45),
           ylim = c(27,65 ),
           expand = F) +
  geom_point(data = n2_T, aes(x = Long_c, y = Lat_c,
                              fill = Bat_species),
             size = 2.5, alpha = 0.7, shape = 21, color = "gray50",  stroke = 0.3) +
  scale_fill_manual(values = alpha(c("#6FB7A7", "#ABD2CB", "#B8D28C",
                                     "#EBD89F", "#EEC9AD", "#E4B888", "#E49A91", 
                                     "#EFAFAF", "#D8B2C4", "#BD98BD", "#C6C4D8", "#ADA8D9",
                                     "#8DAAEE"),0.7)) +
  xlab("") +
  ylab("") +
  theme_bw() 

## Display map
map.color_T <- map.color_T + 
  theme(legend.position = "none")
map.color_T <-  map.color_T+
  theme(legend.text = element_text(size = 8)
  )+ 
  guides(color = guide_legend(override.aes = list(size = 1.2), nrow = 3))
map.color_T <-  map.color_T + theme(axis.text.y = element_text(size = 8), 
                                    axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8))

# Create the barplot
bar_plot_T <- ggplot(aa) +
  geom_bar(aes(x = n1, y = factor(Bat_species, levels = rev(levels(factor(Bat_species)))), 
               fill = Bat_species), colour = "black", 
           stat = "identity", alpha = 0.7, linewidth = 0.35,
           position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c("#6FB7A7", "#ABD2CB", "#B8D28C",
                               "#EBD89F", "#EEC9AD", "#E4B888", "#E49A91", 
                               "#EFAFAF", "#D8B2C4", "#BD98BD", "#C6C4D8", "#ADA8D9",
                               "#8DAAEE"))+
  theme(text = element_text(family = "Source Code Pro"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9)) +
  scale_x_discrete(labels = label_wrap(2)) + #,guide = guide_axis(n.dodge = 1.5))+
  theme_bw() + 
  scale_x_continuous(breaks = seq(0, 54, by = 10))
bar_plot_T <- bar_plot_T + theme(legend.position = "none")+
  guides(fill = "none") 
bar_plot_T <-  bar_plot_T + theme(axis.text.y = element_text(size = 12, face = "bold.italic"), 
                                  axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = ))+
  theme(plot.margin = unit(c(0.2, 0, 0.2, 2), 'lines'))+
  labs(x ="n years")+
  ggtitle("Temporal")+ 
  theme(plot.title = element_text(size = 20, face = "bold"))


# Arrange the plots
arrange_ST <- ggarrange(bar_plot_S,   bar_plot_T,map.color_S, map.color_T, 
                        nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"), #common.legend = T, legend = "bottom", 
                        align = "hv", widths = c(0.1,0.1) )+ theme(plot.background = element_rect(color = "black"))
arrange_ST

# ##Save the Figure
# ggsave('Figures/Figure_1.png', arrange_ST,
#        width = 330, height = 230, units = 'mm', dpi = 300, bg="white")
# 
ggsave('Figures/Figure_1.pdf', arrange_ST,
       width = 330, height = 230, units = 'mm', dpi = 300, bg="white")
