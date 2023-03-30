# venelin tutorial w cetacean data

library(PCMBase)
library(PCMFit)
library(data.table)
library(ggplot2)
library(ggtree)
library(cowplot)

mySpreadsheet <- read.csv("data/MySpreadsheet.csv", sep = ";",na.strings = c("", NA))

# the tree
 tree.dangerous <- ape::read.tree("data/trees/Cetacea_Dangerous_Median.tre")
 tree.d <- PCMTree(tree.dangerous)

# the data
 traits <- read.csv2("data/traits/Racicotetal2019_data_manually.csv") # Cl and Cw as 2 test traits
 
 #test <- traits[, c(6,8)] # 90 specimen but 719 tips in tree
 
# pruning tree with geiger?
 #geiger::treedata(traits$Taxa.clean, tree, sort = T)
 # or use prune2data(tree, vector e.g species names)
 
# reformat data to make it work (matrix)
 # remove empty rows
 tmp <- traits[-c(32:43),]
 
 # begin with one complete trait 
 test <- tmp[,c(1,6)]
 colnames(test) <- c("taxon", "trait")  # change colnames
 row.names(test) <- seq_len(nrow(test)) # change rownames to seq
 
 matrix <- as.matrix(test)              # as matrix
 