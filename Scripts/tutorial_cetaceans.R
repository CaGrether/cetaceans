# venelin tutorial w cetacean data

mySpreadsheet <- read.csv("data/MySpreadsheet.csv", sep = ";",na.strings = c("", NA))

# the tree
 tree.dangerous <- ape::read.tree("data/trees/Cetacea_Dangerous_Median.tre")
 tree <- PCMTree(tree.dangerous)

# the data
 traits <- read.csv2("data/traits/Racicotetal2019_data_manually.csv") # Cl and Cw as 2 test traits
 
 test <- traits[, c(6,8)] # 90 specimen but 719 tips in tree
 
# pruning tree with geiger?
 geiger::treedata(traits$Taxa.clean, test, sort = T)
