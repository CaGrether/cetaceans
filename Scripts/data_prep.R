# cross-check data
# 1. read the input data

#setwd("C:/Users/Caro/Documents/Uni_Erlangen/Studium/Thesis")

mySpreadsheet <- read.csv("data/MySpreadsheet.csv", sep = ";",na.strings = c("", NA))



# 2. read the trees
## set up library
tree.dangerous <- ape::read.tree("data/trees/Cetacea_Dangerous_Median.tre") #very big
tree.risky <- ape::read.tree("data/trees/Cetacea_Risky_Median.tre")         #big tree
tree.safe <- ape::read.tree("data/trees/Cetacea_Safe_Median.tre")           #smaller tree



# a) Dangerous
## names of taxa at tip of tree as vector
tree.dang.labels <- tree.dangerous$tip.label

## cross check taxa (are taxa from the data in the tree)

dangerous <- c() #empty vector to put in taxa which are in the tree

for(i in mySpreadsheet$Taxa.new){  # for each taxon
  
  tmp <- gsub(" ", "_", i) # substitute space for _ in taxon name for every tip label
  
  # is the taxa in the tree
  dangerous <- c(dangerous, 
                 (ifelse (tmp %in% tree.dang.labels, TRUE, FALSE)) ) # write TRUE or FALSE
  
}

mySpreadsheet$dangerous <- dangerous  # new col with TRUE or FALSE



# b) risky
## names of taxa at tip of tree as vector
tree.risk.labels <- tree.risky$tip.label

## cross check taxa (are taxa from the data in the tree)
risky <- c()

for (i in mySpreadsheet$Taxa.new) {   # for each taxon
  
  tmp2 <- gsub(" ","_",i)  # substitute space for "_"
  
  risky <- c(risky,        # is it in the tree
             (ifelse(tmp2 %in% tree.risk.labels, TRUE, FALSE))) #TRUE/FALSE
}

mySpreadsheet$risky <- risky    # new col with TRUE or FALSE



# c) safe
## names of taxa at tip of tree
tree.safe.labels <- tree.safe$tip.label

## cross check taxa (are taxa from the data in the tree)
safe <- c()

for (i in mySpreadsheet$Taxa.new) {    #for each taxon
  
  tmp3 <- gsub(" ","_",i)    # substitue space for "_"
  
  
  safe <- c(safe,            # is it in the tree
            (ifelse(tmp3 %in% tree.safe.labels,TRUE,FALSE))) # TRUE or FALSE
  
}

mySpreadsheet$safe <- safe
View(mySpreadsheet)   # check if everything worked


## check how much of the data is covered
sum(mySpreadsheet$risky==TRUE)    # how many taxa are in risky tree
#length(mySpreadsheet$Specimen.number.from.Racicot.et.al.2019)
nrow(mySpreadsheet)               # how many are in the data


######### genus in tree -----------------------------------------------------------


## check for genus coverage

#my taxa as in mySpreadsheet$taxa_corrected
split.taxa <- unlist(strsplit(mySpreadsheet$Taxa.new, " "))  # split taxon names (genus / species)

#the taxa in trees as in tree.safe/risky/dangerous$tip.label

split.dang.labels <- unlist(strsplit(tree.dang.labels, "_")) # split taxon names in the trees
split.risk.labels <- unlist(strsplit(tree.risk.labels, "_")) # i.e split tip labels
split.safe.labels <- unlist(strsplit(tree.safe.labels, "_")) # and store in container


# get the genus names (without species names)

genus.name <- c()                             # empty container

for (i in 1:length(mySpreadsheet$Taxa.new)) { # for every corrected taxon
  
  genus <- mySpreadsheet$Taxa.new[i]          # put into temp container
  
  if (!is.na(genus)){                         # if there is a taxon present
    genus.name <- c(genus.name, unlist(strsplit(mySpreadsheet$Taxa.new[i], " "))[1])
  }                                           # put only genus name into container
  
}

duplicated(genus.name)                        # how many genera are present multiple times
genus.o <- genus.name[!duplicated(genus.name)]# put only unique genus names into genus.o



## loops for cross-checking how many genus names are in trees
# dangerous tree

dangerous.genus <- c()  # empty container for genera present in tree
genus <- c()            # empty container for storing all genera

# loop
for (i in genus.o) {   # for each genus
  
  #check whether it's in tree
  dangerous.genus <- c(dangerous.genus, (ifelse(i %in% split.dang.labels, TRUE, FALSE)) )
                     # if present in tree, put into container
  genus <- c(genus, print(i))    # put it into genus container
  
}

testing.genus <- data.frame(genus, dangerous.genus)  # make a data frame with all genera
                                                     # to compare with coverage in all trees

# risky tree

risky.genus <- c()  # empty container for genera present in tree

for (i in genus.o) {
  
  risky.genus <- c(risky.genus, (ifelse(i %in% split.risk.labels, TRUE, FALSE)))
                 # if present in tree, put into container
}

testing.genus$risky.genus <- risky.genus  # add to others for comparison

# safe tree

safe.genus <- c()     # empty container for genera present in tree

for (i in genus.o) {
  
  safe.genus <- c(safe.genus, (ifelse(i %in% split.safe.labels, TRUE, FALSE)))
                # if present in tree, put into container
}

testing.genus$safe.genus <- safe.genus     # add to others for comparison


# how many taxa are represented in tree?
sum(testing.genus$dangerous.genus==TRUE)   # number of taxa in dangerous tree
length(genus)                              # 

# inner ear data ----------------------------------------------------------

#read in data

inner.ear <- read.csv("Data/Racicotetal2019_data_manually.csv", sep = ";")

names(inner.ear)
names(inner.ear)[1] <- "Taxon.clean"


# get NA into empty rows

inner.ear[inner.ear==""] <- NA


#  create new column for genus names

inner.ear$genus <- NA

# put genus names in inner.ear

for (i in 1:length(inner.ear$Taxon.clean)) {
  
  name <- inner.ear$Taxon.clean[i]
  
  if (!is.na(name)){
    genus <- unlist(strsplit(name, " "))[1]
    
    inner.ear$genus[i] <- testing.genus$genus[which(testing.genus$genus == genus)]
  }
  else {
    inner.ear$genus[i] <- NA
  }
}



# take the mean -----------------------------------------------------------


# create columns for the trait data

testing.genus$mean.Cl <- NA
testing.genus$mean.SBL <- NA
testing.genus$mean.Cw <- NA
testing.genus$mean.Ch <- NA
testing.genus$mean.W2 <- NA
testing.genus$mean.ITD <- NA
testing.genus$mean.GAN <- NA
testing.genus$mean.FC <- NA
testing.genus$"mean.T" <- NA



# Cl
for (i in 1:length(testing.genus$genus)) {
  
  genus <- testing.genus$genus[i]
  
  row <- which(inner.ear$genus == genus)
  
  if (length(row)>0) { 
    #mean or median
    testing.genus$mean.Cl[i] <- mean(as.numeric(gsub(",", ".", inner.ear$Cl[row], fixed = T)), na.rm = T) 
  }
  else {
    testing.genus$mean.Cl[i] <- NA
  }
  
}


# SBL

for (i in 1:length(testing.genus$genus)) {
  
  genus <- testing.genus$genus[i]
  
  row <- which(inner.ear$genus == genus)
  
  if (length(row)>0) { 
    #mean or median
    testing.genus$mean.SBL[i] <- mean(as.numeric(gsub(",", ".", inner.ear$SBL[row], fixed = T)), na.rm = T) 
  }
  else {
    testing.genus$mean.SBL[i] <- NA
  }
  
}


# Cw

for (i in 1:length(testing.genus$genus)) {
  
  genus <- testing.genus$genus[i]
  
  row <- which(inner.ear$genus == genus)
  
  if (length(row)>0) { 
    #mean or median
    testing.genus$mean.Cw[i] <- mean(as.numeric(gsub(",", ".", inner.ear$Cw[row], fixed = T)), na.rm = T) 
  }
  else {
    testing.genus$mean.Cw[i] <- NA
  }
  
}

# Ch

for (i in 1:length(testing.genus$genus)) {
  
  genus <- testing.genus$genus[i]
  
  row <- which(inner.ear$genus == genus)
  
  if (length(row)>0) { 
    #mean or median
    testing.genus$mean.Ch[i] <- mean(as.numeric(gsub(",", ".", inner.ear$Ch[row], fixed = T)), na.rm = T) 
  }
  else {
    testing.genus$mean.Ch[i] <- NA
  }
  
}

# W2

for (i in 1:length(testing.genus$genus)) {
  
  genus <- testing.genus$genus[i]
  
  row <- which(inner.ear$genus == genus)
  
  if (length(row)>0) { 
    #mean or median
    testing.genus$mean.W2[i] <- mean(as.numeric(gsub(",", ".", inner.ear$W2[row], fixed = T)), na.rm = T) 
  }
  else {
    testing.genus$mean.W2[i] <- NA
  }
  
}

# ITD

for (i in 1:length(testing.genus$genus)) {
  
  genus <- testing.genus$genus[i]
  
  row <- which(inner.ear$genus == genus)
  
  if (length(row)>0) { 
    #mean or median
    testing.genus$mean.ITD[i] <- mean(as.numeric(gsub(",", ".", inner.ear$ITD[row], fixed = T)), na.rm = T) 
  }
  else {
    testing.genus$mean.ITD[i] <- NA
  }
  
}

# GAN 

for (i in 1:length(testing.genus$genus)) {
  
  genus <- testing.genus$genus[i]
  
  row <- which(inner.ear$genus == genus)
  
  if (length(row)>0) { 
    #mean or median
    testing.genus$mean.GAN[i] <- mean(as.numeric(gsub(",", ".", inner.ear$GAN[row], fixed = T)), na.rm = T) 
  }
  else {
    testing.genus$mean.GAN[i] <- NA
  }
  
}

# FC

for (i in 1:length(testing.genus$genus)) {
  
  genus <- testing.genus$genus[i]
  
  row <- which(inner.ear$genus == genus)
  
  if (length(row)>0) { 
    #mean or median
    testing.genus$mean.FC[i] <- mean(as.numeric(gsub(",", ".", inner.ear$FC[row], fixed = T)), na.rm = T) 
  }
  else {
    testing.genus$mean.FC[i] <- NA
  }
  
}

# T

for (i in 1:length(testing.genus$genus)) {
  
  genus <- testing.genus$genus[i]
  
  row <- which(inner.ear$genus == genus)
  
  if (length(row)>0) { 
    #mean or median
    testing.genus$mean.T[i] <- mean(as.numeric(gsub(",", ".", inner.ear$"T"[row], fixed = T)), na.rm = T) 
  }
  else {
    testing.genus$mean.T[i] <- NA
  }
  
}


# result: if inner ear measurements don't exist, output is NA
# if measurements for specific trait are NA, output is NaN



# merge with Leon ---------------------------------------------------------


# read in Leons data

Kasuya <- read.csv("Data/sheet_merge.csv", sep = ",")
View(Kasuya)

names(Kasuya)
names(Kasuya)[1] <- "genus"

merge <- merge(Kasuya, testing.genus, by = "genus" , all= T )
View(merge)


# write.csv2(merge, "Data/merged_data.csv", row.names = T)


