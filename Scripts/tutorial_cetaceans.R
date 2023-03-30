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
 
# use Venelins code
 # but still taxa don't correspond (names are different, amount of taxa
 # and order is still wrong)
 
 ## cross check taxa (are taxa from the data in the tree) from thesis
 idx <- c() #empty vector to put in taxa which are in the tree
 
 for(i in matrix$taxon){  # for each taxon
   
   tmp <- gsub(" ", "_", i) # substitute space for _ in taxon name for every tip label
   
   # is the taxa in the tree
   idx <- c(idx, 
            (ifelse (tmp %in% tree.d$tip.label, TRUE, FALSE)) ) # write TRUE or FALSE
   
 }

 matrix$present <- idx  # new col with TRUE or FALSE
 
 # with merged.csv: remove all taxa absent from tree and all with NA for mean.Cl
 
 
 
 
 # create plots of tree and data
 plTree.d <- PCMTreePlot(tree.d, layout="fan") +
   geom_tiplab2(size = 2) + 
   geom_nodelab(size = 2, color = "black") + 
   geom_treescale(width = max(PCMTreeNodeTimes(tree.d)), x = 0, linesize = .25, fontsize = 2, offset = 79)
 
 plmat <- PCMPlotTraitData2D(
   matrix[, seq_len(PCMTreeNumTips(tree.d))], 
   tree.d, 
   scaleSizeWithTime = FALSE,
   numTimeFacets = 4) +
   geom_text(
     aes(x = x, y = y, label = id, color = regime), 
     size=2, 
     position = position_jitter(.4, .4)) +
   theme_bw() +
   theme(legend.position = "bottom") 
 