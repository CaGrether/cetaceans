setwd("C:/Users/Caro/Documents/Uni_Erlangen/Studium/Thesis")


#read in data
merge <- read.csv("Data/merged_data.csv", sep = ";", dec = ",")
tree.dangerous <- ape::read.tree("Trees/Cetacea_Dangerous_Median.tre")
tree.risky <- ape::read.tree("Trees/Cetacea_Risky_Median.tre")
tree.safe <- ape::read.tree("Trees/Cetacea_Safe_Median.tre")

View(merge)

# start pipeline ----------------------------------------------------------

# extract genus names from tip labels
genus.names.tree <- c()

for (i in tree.safe$tip.label){
  
  genus.names.tree <- c(genus.names.tree, unlist(strsplit(i, "_"))[1])
}


# delete the ones which occur more than once

keepers <- c()


for (i in unique(genus.names.tree)) {
  
  tmp <- grep(i, tree.safe$tip.label)
  if(length(tmp) > 1)
    keepers <- c(keepers, tmp[1])
  #  keepers <- c(keepers, sample(tmp, 1))
  else 
    keepers <- c(keepers, tmp)
  
}
new.tree <- ape::keep.tip(tree.safe, tree.safe$tip.label[keepers])


# rename tip labels as only genera name

tree.genus <- new.tree
genus <- c()

for (i in tree.genus$tip.label) {
  
  genus <- c(genus, unlist(strsplit(i, "_"))[1])
}

tree.genus$tip.label <- genus



# write.tree(tree.genus, file = "Trees/genus_tree.tre")


# cross check which tip labels in tree.genus are also represented in my data

tree.done <- tree.genus
keep <- c()


for (i in unique(tree.done$tip.label)) {
  
  # check whether it's in merge$genus
  row <- which(i == merge$genus)
  
  # if yes, keep it
  if (length(row)>0){
    keep <- c(keep, (which(i==tree.done$tip.label)))
  }
  
}




tree.finished <- ape::keep.tip(tree.done, tree.done$tip.label[keep])
plot(tree.finished)
s


# tree without NA consisting of $mean.Cl

merge.sub <- subset(merge, !is.na(mean.Cl) & safe.genus==TRUE)#49

mean.Cl <- merge.sub$mean.Cl
names(mean.Cl) <- merge.sub$genus

#make a tree with all genera
tree.mean.Cl <- tree.genus

dont.keep <- c()

for (i in unique(tree.mean.Cl$tip.label)) {
  
  # check whether it's in merge$genus
  row <- which(i == merge.sub$genus)
  
  # if no, dont keep it
  if(length(row) == 0){
    dont.keep <- c(dont.keep, (which(i==tree.mean.Cl$tip.label)))
  }
  
}
#tree.mean.Cl <- ape::keep.tip(tree.mean.Cl, tree.mean.Cl$tip.label[keep])
tree.mean.Cl <- ape::drop.tip(tree.mean.Cl, tree.mean.Cl$tip.label[dont.keep])

length(tree.mean.Cl$tip.label)
length(merge.sub$genus)

plot(tree.mean.Cl)



phytools::phenogram(tree.mean.Cl, mean.Cl, ftype = "off")



