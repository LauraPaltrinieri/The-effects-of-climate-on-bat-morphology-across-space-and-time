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

# Clear environment
rm(list=ls())

# From now on, the analysis is divided based on the dataset used
#----------------------
## FAL dataset
#----------------------

# load database and elton traits data to remove marine species
data <- read.csv("Results/Temporal/Temporal_FAL_effect_size.csv")

# generating list of species
species <- sort(unique(as.character(data$Bat_species))) #33species

# Formatting species data
# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.

taxa <- tnrs_match_names(species) 

# according to the `approximate_match` column, there might be 
# 0 typos in the species list 
taxa[taxa$approximate_match==TRUE,] # no species returned

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa[taxa$number_matches != 1,]
ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa, ott_id = ott_id_tocheck[i]))
}

# check synonyms and change name accordingly
fix_taxa <- taxa[taxa$is_synonym==TRUE,] # 0 species returned


# retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name")
plot(tree, cex=.5, label.offset =.1, no.margin = TRUE)

# Dealing with polytomies
# we run the following function to check for the existence of polytomies.
# If polytomies exist, the output will be `FALSE`, and vice versa.

is.binary(tree) # there are some polytomies

# to take care of these polytomies, we are going to use a 
# randomization approach
set.seed(111) #making it replicable, at least for this version of R (i.e. v.4.0.2)
tree_random <- multi2di(tree,random=TRUE)
is.binary(tree_random)


# Final checks
# exploring whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.

tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
intersect(as.character(tree_random$tip.label), as.character(species))

species[!species %in% as.character(tree_random$tip.label)] #listed in our database but not in the tree
tree_random$tip.label[!as.character(tree_random$tip.label) %in% species] # listed in the tree but not in our database

tiff("Results/Temporal/Temporal_FAL_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "Results/Temporal/Temporal_FAL_tree_random.Rdata")

# Computing branch lengths
# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(data$speciesname, as.character(tree_random$tip.label)) # Null
setdiff(as.character(tree_random$tip.label),data$Bat_species)# 0

# exclude species in the tree that are not in dataset
drops <- tree_random$tip.label[!tree_random$tip.label %in% data$Bat_species]
tree_random.fixed <- drop.tip(tree_random, drops)

# save the new tree
write.tree(tree_random.fixed, file = "Results/Temporal/Temporal_FAL_tree_random_fixed.tre")

# compute branch lengths of tree
phylo_branch <- compute.brlen(tree_random.fixed, method = "Grafen", power = 1)

# check if tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# Phylogenetic matrix

# matrix to be included in the models
bat_phylo_cor <- vcv(phylo_branch, cor = T)

# remove rows not in correlation matrix
data_ph <- data[which(data$Bat_species %in% rownames(bat_phylo_cor)),] 

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID <- data.frame(Bat_species = unique(data_ph$Bat_species), SPID = paste0("SP",1:length(unique(data_ph$Bat_species))))
SpID$Bat_species <- as.character(SpID$Bat_species)
data_ph <- inner_join(data_ph,SpID, by = "Bat_species")

# finally, save matrix for future analyses
save(bat_phylo_cor, file = "Results/Temporal/Temporal_FAL_bat_phylo_cor.Rdata")

# exporting fixed dataset for analyses
write.csv(data_ph, 
          "Data/Temporal/FAL_ph.csv", row.names = FALSE)
# exporting species list
write.csv(species, "Results/Temporal/Temporal_FAL_species_list.csv")

rm(list=ls())

#--------------------------
## Pairwise FAL-BM dataset
#--------------------------
# load database and elton traits data to remove marine species
data <- read.csv("Results/Temporal/Temporal_pairwiseFAL_BM_effect_size.csv")

# generating list of species
species <- sort(unique(as.character(data$Bat_species))) #33species

# Formatting species data
# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.

taxa <- tnrs_match_names(species) 

# according to the `approximate_match` column, there might be 
# 0 typos in the species list  
taxa[taxa$approximate_match==TRUE,] # no species returned

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa[taxa$number_matches != 1,]
ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa, ott_id = ott_id_tocheck[i]))
}

# check synonyms and change name accordingly
fix_taxa <- taxa[taxa$is_synonym==TRUE,] # 0 species returned

# retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name")
plot(tree, cex=.5, label.offset =.1, no.margin = TRUE)

# Dealing with polytomies
# we run the following function to check for the existence of polytomies.
# If polytomies exist, the output will be `FALSE`, and vice versa.

is.binary(tree) # there are some polytomies

# to take care of these polytomies, we are going to use a 
# randomization approach
set.seed(111) #making it replicable, at least for this version of R (i.e. v.4.0.2)
tree_random <- multi2di(tree,random=TRUE)
is.binary(tree_random)


# Final checks
# exploring whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.

tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
intersect(as.character(tree_random$tip.label), as.character(species))

species[!species %in% as.character(tree_random$tip.label)] #listed in our database but not in the tree
tree_random$tip.label[!as.character(tree_random$tip.label) %in% species] # listed in the tree but not in our database

tiff("Results/Temporal/Temporal_pairwiseFAL_BM_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "Results/Temporal/Temporal_pairwiseFAL_BM_tree_random.Rdata")

# Computing branch lengths
# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(data$speciesname, as.character(tree_random$tip.label)) # Null
setdiff(as.character(tree_random$tip.label),data$Bat_species)# 0

# exclude species in the tree that are not in dataset
drops <- tree_random$tip.label[!tree_random$tip.label %in% data$Bat_species]
tree_random.fixed <- drop.tip(tree_random, drops)

# save the new tree
write.tree(tree_random.fixed, file = "Results/Temporal/Temporal_pairwiseFAL_BM_tree_random_fixed.tre")

# compute branch lengths of tree
phylo_branch <- compute.brlen(tree_random.fixed, method = "Grafen", power = 1)

# check if tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# Phylogenetic matrix

# matrix to be included in the models
bat_phylo_cor <- vcv(phylo_branch, cor = T)

# remove rows not in correlation matrix
data_ph <- data[which(data$Bat_species %in% rownames(bat_phylo_cor)),] 

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID <- data.frame(Bat_species = unique(data_ph$Bat_species), SPID = paste0("SP",1:length(unique(data_ph$Bat_species))))
SpID$Bat_species <- as.character(SpID$Bat_species)
data_ph <- inner_join(data_ph,SpID, by = "Bat_species")

# finally, save matrix for future analyses
save(bat_phylo_cor, file = "Results/Temporal/Temporal_pairwiseFAL_BM_bat_phylo_cor.Rdata")

# exporting fixed dataset for analyses
write.csv(data_ph, 
          "Data/Temporal/pairwiseFAL_BM_ph.csv", row.names = FALSE)
# exporting species list
write.csv(species, "Results/Temporal/Temporal_pairwiseFAL_BM_species_list.csv")

