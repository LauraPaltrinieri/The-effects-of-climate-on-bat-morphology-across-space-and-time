# Script: Spatial phylogeny OTL
# Date: 11th November, 2024


# This script is used to build the phylogeny for European bat species

#load libraries
library(dplyr)
library(ape)
library(treebase)
library(rotl)
library(diagram)
library(stringr)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("treeio")
# BiocManager::install("ggtree")
library(treeio)
library(ggtree)
library(ggpubr)
library(viridisLite)
library(viridis)

# Clear environment
rm(list=ls())

# From now on, the analysis is divided based on the dataset used: 
# the "FAL dataset" and the "pairwise FAL-BM dataset"
#----------------------
## FAL dataset
#----------------------

# Load database and elton traits data to remove marine species
data <- read.csv("Results/Spatial/Spatial_FAL_effect_size.csv")

# Generating list of species
species <- sort(unique(as.character(data$Bat_species))) #19 species

# Formatting species data
# Obtaining dataframe listing the Open Tree identifiers potentially 
# Matching our list of species.

taxa <- tnrs_match_names(species) 

# According to the `approximate_match` column, there might be 
# 0 typos in the species list 
taxa[taxa$approximate_match == TRUE,] # no species returned

# Exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa[taxa$number_matches != 1,]
ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa, ott_id = ott_id_tocheck[i]))
}

# Check synonyms and change name accordingly
fix_taxa <- taxa[taxa$is_synonym==TRUE,] # 0 species returned

# Retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name")
plot(tree, cex=.5, label.offset =.1, no.margin = TRUE)

# Dealing with polytomies
# We run the following function to check for the existence of polytomies.
# If polytomies exist, the output will be `FALSE`, and vice versa.

is.binary(tree) # there are some polytomies

# To take care of these polytomies, we are going to use a 
# randomization approach
set.seed(111) # Making it replicable, at least for this version of R (i.e. v.4.0.2)
tree_random <- multi2di(tree,random=TRUE)
is.binary(tree_random)


# Final checks
# Exploring whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.

tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
intersect(as.character(tree_random$tip.label), as.character(species))

species[!species %in% as.character(tree_random$tip.label)] # Listed in our database but not in the tree
tree_random$tip.label[!as.character(tree_random$tip.label) %in% species] # Listed in the tree but not in our database

tiff("Results/Spatial/Spatial_FAL_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# We can now save the tree
save(tree_random, file = "Results/Spatial/Spatial_FAL_tree_random.Rdata")

# Computing branch lengths
# We are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# Before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(data$speciesname, as.character(tree_random$tip.label)) # Null
setdiff(as.character(tree_random$tip.label),data$Bat_species) # 0

# Exclude species in the tree that are not in dataset
drops <- tree_random$tip.label[!tree_random$tip.label %in% data$Bat_species]
tree_random.fixed <- drop.tip(tree_random, drops)

# Save the new tree
write.tree(tree_random.fixed, file = "Results/Spatial/Spatial_FAL_tree_random_fixed.tre")

# Compute branch lengths of tree
phylo_branch <- compute.brlen(tree_random.fixed, method = "Grafen", power = 1)

# Check if tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# Phylogenetic matrix

# Matrix to be included in the models
bat_phylo_cor <- vcv(phylo_branch, cor = T)

# Remove rows not in correlation matrix
data_ph <- data[which(data$Bat_species %in% rownames(bat_phylo_cor)),] 

## Create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID <- data.frame(Bat_species = unique(data_ph$Bat_species), SPID = paste0("SP",1:length(unique(data_ph$Bat_species))))
SpID$Bat_species <- as.character(SpID$Bat_species)
data_ph <- inner_join(data_ph,SpID, by = "Bat_species")

# Finally, save matrix for future analyses
save(bat_phylo_cor, file = "Results/Spatial/Spatial_FAL_bat_phylo_cor.Rdata")

# Exporting fixed dataset for analyses
write.csv(data_ph, 
          "Data/Spatial/FAL_ph.csv", row.names = FALSE)
# Exporting species list
write.csv(species, "Results/Spatial/Spatial_FAL_species_list.csv")

#---
#Create the philogenetic tree with the branch colored based on the z-scores
# Convert the tree in treedata to can add new data
tree_data <- treeio::as.treedata(tree_random)

# Merge the z-scores data with the tree data
# Estrai l'albero `phylo` dal `treedata`
phylo_tree <- treeio::get.tree(tree_data)
# Modifica i nomi delle etichette
phylo_tree$tip.label <- gsub("^(\\w)\\w*\\s(\\w+)", "\\1. \\2", phylo_tree$tip.label)

data_m <- subset(data, Sex == "male")
data_m$Bat_species <- sapply(data_m$Bat_species, function(x) {
  parts <- strsplit(x, " ")[[1]]          
  paste0(substr(parts[1], 1, 1), ". ", parts[2])  
})
tree_datam <- full_join(phylo_tree, data_m, by = c("label" = "Bat_species"))


data_f <- subset(data, Sex == "female")
data_f$Bat_species <- sapply(data_f$Bat_species, function(x) {
  parts <- strsplit(x, " ")[[1]]          
  paste0(substr(parts[1], 1, 1), ". ", parts[2])  
})
tree_dataf <- full_join(phylo_tree, data_f, by = c("label" = "Bat_species"))

##Plot
#Tmax
Tmax_m <- ggtree(tree_datam, layout='circular', ladderize = T, size=4) + 
  geom_tree(aes(color = yi_Tmax),  size=3.5) + 
  scale_color_viridis(name = "Tmax 
z-scores")+
  geom_tiplab(aes(color = yi_Tmax), hjust = -.1, size = 6.5) + 
  xlim(0, 15) + 
  theme(legend.position = c(-0.02, .7),  # Aumenta la dimensione del titolo
        legend.text = element_text(size = 15),  # Aumenta la dimensione del testo della legenda
        legend.title = element_text(size = 18),  # Aumenta la dimensione del titolo della legenda
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(size = 25))+ 
  labs(title = "Males") 


Tmax_f <- ggtree(tree_dataf, layout='circular', ladderize = T, size=4) + 
  geom_tree(aes(color = yi_Tmax),  size=3.5) + 
  scale_color_viridis(name = "Tmax 
z-scores")+
  geom_tiplab(aes(color=yi_Tmax), hjust = -.1, size = 6.5) + 
  xlim(0, 15) + 
  theme(legend.position = c(-0.02, .7),  # Aumenta la dimensione del titolo
        legend.text = element_text(size = 15),  # Aumenta la dimensione del testo della legenda
        legend.title = element_text(size = 18),  # Aumenta la dimensione del titolo della legenda
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(size = 25)) + 
  labs(title = "Females") 


#Summer Precipitation
SP_m <- ggtree(tree_datam, layout='circular', ladderize = T, size=4) + 
  geom_tree(aes(color = yi_SP),  size=3.5) + 
  scale_color_viridis(name = "SP  
z-scores")+
  geom_tiplab(aes(color=yi_SP), hjust = -.1, size = 6.5) + 
  xlim(0, 15) + 
  theme(legend.position = c(-0.02, .7),  
        legend.text = element_text(size = 15),  
        legend.title = element_text(size = 18),  
        legend.key.size = unit(1, "cm") ) 

SP_f <- ggtree(tree_dataf, layout='circular', ladderize = T, size=4) + 
  geom_tree(aes(color = yi_SP), size=3.5) + 
  scale_color_viridis(name = "SP 
z-scores")+
  geom_tiplab(aes(color=yi_SP), hjust = -.1, size = 6.5) + 
  xlim(0, 15) + 
  theme(legend.position = c(-0.02, .7),  
        legend.text = element_text(size = 15),  
        legend.title = element_text(size = 18),  
        legend.key.size = unit(1, "cm") )

#Save figures
ggarrange(Tmax_f, Tmax_m, SP_f,SP_m)
# ggsave(filename='Figures/Figure_S3.png',
#        width = 500, height = 380, units = 'mm', dpi = 1000)
ggsave(filename='Figures/Figure_S3.pdf',
       width = 500, height = 380, units = 'mm', dpi = 1000)
#---

rm(list=ls())

#--------------------------
## Pairwise FAL-BM dataset 
#--------------------------
# Load database and elton traits data to remove marine species
data <- read.csv("Results/Spatial/Spatial_pairwiseFAL_BM_effect_size.csv")

# Generating list of species
species <- sort(unique(as.character(data$Bat_species))) # 19 species

# Formatting species data
# Obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.

taxa <- tnrs_match_names(species) 

# According to the `approximate_match` column, there might be 
# 0 typos in the species list 
taxa[taxa$approximate_match ==TRUE,] # no species returned

# Exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa[taxa$number_matches != 1,]
ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa, ott_id = ott_id_tocheck[i]))
}

# Check synonyms and change name accordingly
fix_taxa <- taxa[taxa$is_synonym==TRUE,] # 0 species returned

# Retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name")
plot(tree, cex=.5, label.offset =.1, no.margin = TRUE)

# Dealing with polytomies
# We run the following function to check for the existence of polytomies.
# If polytomies exist, the output will be `FALSE`, and vice versa.

is.binary(tree) # There are some polytomies

# To take care of these polytomies, we are going to use a 
# randomization approach
set.seed(111) # Making it replicable, at least for this version of R (i.e. v.4.0.2)
tree_random <- multi2di(tree,random=TRUE)
is.binary(tree_random)


# Final checks
# Exploring whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.

tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
intersect(as.character(tree_random$tip.label), as.character(species))

species[!species %in% as.character(tree_random$tip.label)] # Listed in our database but not in the tree
tree_random$tip.label[!as.character(tree_random$tip.label) %in% species] # Listed in the tree but not in our database

tiff("Results/Spatial/Spatial_pairwiseFAL_BM_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# We can now save the tree
save(tree_random, file = "Results/Spatial/Spatial_pairwiseFAL_BM_tree_random.Rdata")

# Computing branch lengths
# We are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# Before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(data$speciesname, as.character(tree_random$tip.label)) # Null
setdiff(as.character(tree_random$tip.label),data$Bat_species)# 0

# Exclude species in the tree that are not in dataset
drops <- tree_random$tip.label[!tree_random$tip.label %in% data$Bat_species]
tree_random.fixed <- drop.tip(tree_random, drops)

# Save the new tree
write.tree(tree_random.fixed, file = "Results/Spatial/Spatial_pairwiseFAL_BM_tree_random_fixed.tre")

# Compute branch lengths of tree
phylo_branch <- compute.brlen(tree_random.fixed, method = "Grafen", power = 1)

# Check if tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# Phylogenetic matrix

# Matrix to be included in the models
bat_phylo_cor <- vcv(phylo_branch, cor = T)

# Remove rows not in correlation matrix
data_ph <- data[which(data$Bat_species %in% rownames(bat_phylo_cor)),] 

## Create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID <- data.frame(Bat_species = unique(data_ph$Bat_species), SPID = paste0("SP",1:length(unique(data_ph$Bat_species))))
SpID$Bat_species <- as.character(SpID$Bat_species)
data_ph <- inner_join(data_ph,SpID, by = "Bat_species")

# Finally, save matrix for future analyses
save(bat_phylo_cor, file = "Results/Spatial/Spatial_pairwiseFAL_BM_bat_phylo_cor.Rdata")

# Exporting fixed dataset for analyses
write.csv(data_ph, 
          "Data/Spatial/pairwiseFAL_BM_ph.csv", row.names = FALSE)
# Exporting species list
write.csv(species, "Results/Spatial/Spatial_pairwiseFAL_BM_species_list.csv")
