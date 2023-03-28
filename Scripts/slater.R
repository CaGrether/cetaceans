library(ape)
library(geiger)
library(phytools)
library(mvtnorm)
library(ggplot2)
library(OUwie)
library(RPANDA) # use github version for fit_t_env() as of January 2017

# measurements:
#CL: Length of the cochlear canal
#SBL: Length of the secondary bony lamina
#CW: Cochlear width
#CH: Cochlear length
#W2: Cochlear width perpendicular to CW
#ITD: Inter-turn distance
#GAN: Spiral ganglion canal width at the first quarter-turn
#FC: Surface area of the fenestra cochleae
#T: Number of turns


setwd("C:/Users/Caro/Documents/Uni_Erlangen/Studium/Thesis")

source("Scripts/slater_functions.R")


## read in tree and data

phy <- ladderize(ape::read.tree("Trees/genus_tree.tre")) 

d <- read.csv2("Data/merged_data.csv", stringsAsFactors = F, row.names = 1)


# testing phylogenetic signal ---------------------------------------------


# empty vectors for phylogenetic signal
psk <- c()
psk.p <- c()
psl <- c()
n.trait <- c()

# test phylogenetic signal
for (i in c(24:32)){
        
        # take one trait
        a <- data.frame(d[i])
        # convert zeros to NA
        a[which(a == 0),] = NA 
        
        # take the log
        logx1 <- log(a[1])
        a <- data.frame(logx1)
        row.names(a) <- d$genus
        x1 <- na.omit(a)
        
        td <- treedata(phy, x1, sort=T)
        
        if(length(x1[,1]) == length(td$data)) row.names(td$data)<-row.names(x1)
        
        phy.tmp <- td$phy
        traitx1 <-setNames(td$data[,1], rownames(td$data))
        
        psk.tmp <- phylosig(phy.tmp, traitx1, method="K", test=TRUE)
        psl.tmp <- phylosig(phy.tmp, traitx1, method="lambda", test=TRUE)
        
        psk <- c(psk, psk.tmp$K)
        psk.p <- c(psk.p, psk.tmp$P)
        psl <- c(psl, psl.tmp$lambda)
        n.trait <- c(n.trait, i)
}

kappa <- round(psk, 3)
p.kappa <- round(psk.p, 3)
lambda <- round(psl, 3)

signal <- data.frame(kappa, p.kappa, lambda, n.trait)
# write.csv2(signal, "Data/phyl_sig.csv")



# calculate models --------------------------------------------------------


res.alltraits <- data.frame()

for (i in 24:32) {
        # take one trait
        a <- data.frame(d[i])
        # convert zeros to NA
        a[which(a == 0),] = NA 
        
        # take the log
        logx1 <- log(a[1])
        a <- data.frame(logx1)
        row.names(a) <- d$genus
        x1 <- na.omit(a)
        
        td <- treedata(phy, x1, sort=T)
        
        if(length(x1[,1]) == length(td$data)) row.names(td$data)<-row.names(x1)
        
        phy.tmp <- td$phy
        traitx1 <-setNames(td$data[,1], rownames(td$data))
        
        ## model fitting 
        ## basic models ##
        bm <- fitContinuous(phy.tmp, traitx1, model="BM")
        acdc <- fitContinuous(phy.tmp, traitx1 , model= "EB", bounds=list(a=c(-1, 1)))
        (ou <- fitContinuous(phy.tmp, traitx1, model="OU")) ##
        trend <- fitContinuous(phy.tmp, traitx1, model="drift")
        
        
        # temp-dependent rates
        
        #zd <- read.csv("2008CompilationData.csv") #d18o data
        #zd <- zd[-which(zd$age>40), ] 
        #zd <- zd[-which(is.na(zd$oxygen)),]
        
        #spl <- smooth.spline(-zd$age, zd$oxygen, df = 15) 
        #xx <- cbind(-spl$x, spl$y)
        #temp.rate <- fit_t_env(phy, logTL, xx, error=sem, model= "EnvLin", 
        #method="Nelder-Mead", control=list(maxit=20000))
        
        ## trend - shift model ##
        ts <- trend.shift(phy.tmp, traitx1, model="trend.shift", n.iter=40)
        
        #####  use OUwie to test for a shift in rate through time  #####
        
        #ouwie.states <- traitx1
        #ouwie.data <- data.frame(species=names(traitx1), trait=traitx1)
        
        # slice1bms <- OUwie.slice(phy, ouwie.data,model=c("BMS"), root.age = max(diag(vcv(phy))), root.station=FALSE, timeslices=c(NA),mserr="known")
        
        ## compute AICc weights
        
        # all.model.fit <- aicw(setNames(c(bm$opt$aicc, acdc$opt$aicc, ou$opt$aicc, trend$opt$aicc, ts$aicc), c("BM", "ACDC","SSP", "Trend",   "BM to \nTrend")))
        # temp.rate$aicc, slice1bms$AICc, "Temp-dep \nRates", "Rate Shift",
        
        
        # all.model.fit[rev(order(all.model.fit$w)),]
        
        #### model parameters and fits ####
        
        sigsq <- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$sigsq))
        root <- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$z0))
        aicc <- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$aicc))
        lnL<- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$lnL))
        param1 <- c(NA, acdc$opt$a, trend$opt$drift, ou$opt$alpha)
        shift.time <- rep(NA, length(param1))
        k <- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$k))
        
        res <- cbind(sigsq, root, param1, shift.time, lnL, k, aicc)
        #rownames(res) <- c("BM", "ACDC", "trend", "SSP")
        #res <- rbind(res, c(temp.rate$param[1], temp.rate$root, temp.rate$param[2], NA, temp.rate$LH,temp.rate$free.parameters, temp.rate$aicc))
        
        #rownames(res) <- c("BM", "ACDC", "trend", "SSP", "TTR")
        
        #bms <- c(slice1bms$solution[2,1], slice1bms$theta[1,1], slice1bms$solution[2,2], slice1bms$timeslices[2], slice1bms$loglik, 4, slice1bms$AICc)
        tshift <- c(ts$rate1, ts$root, ts$mu, ts$shift.time, ts$lnl, 4, ts$aicc)
        
        #res <- rbind(res, bms, tshift) 
        res <- rbind(as.data.frame(res),tshift, make.row.names = F)
        models <- c("BM", "ACDC", "trend", "OU", "tshift")
        res <- cbind(models,res)
        
        all.model.fit <- aicw(setNames(c(bm$opt$aicc, acdc$opt$aicc, trend$opt$aicc, ou$opt$aicc, ts$aicc), c("BM", "ACDC", "trend", "OU", "tshift")))
        
        res$deltaAIC <- round(all.model.fit$delta, 3)
        res$weights <- round(all.model.fit$w, 2)
        
        res <- res[rev(order(res$aicc)), ]
        
        
        res$trait = i
        res.alltraits = rbind(res.alltraits, res)
}

res.alltraits

# mean of shift time
## all shift time
shift.time.all <- na.omit(res.alltraits$shift.time)
time.all <- mean(shift.time.all)

## only for deltaAIC of tshift <2
shift.time.best <- res.alltraits$shift.time[which(res.alltraits$deltaAIC<2 & res.alltraits$models=="tshift")]
time.best <- mean(shift.time.best)

# output
#mean.times <- data.frame(time.all, time.best)
#write.csv2(mean.times, file = "Data/mean_times.csv", row.names = F)


# compare results ---------------------------------------------------------


# read in data
res.alltraits <- read.csv2("Data/res_alltraits.csv", )

## unique tables for each trait

for (k in unique(res.alltraits$trait)) {
        # take trait x
        trait1 <-  res.alltraits[which(res.alltraits$trait == k ),]
        trait <- unique(trait1$trait)
        
        # paste the models which have deltaAIC less than 2
        model <- paste0(trait1$models[which(trait1$deltaAIC <2)], collapse = " or ")
        
        # create dataframe with name
        assign(paste0("df", trait), data.frame(trait, model))
}

df24$description <- paste("Length of the cochlear canal")
df25$description <- paste("Length of the secondary bony lamina")
df26$description <- paste("Cochlear width")
df27$description <- paste("Cochlear length")
df28$description <- paste("Cochlear width perpendicular to CW")
df29$description <- paste("Inter-turn distance")
df30$description <- paste("Spiral ganglion canal width at the first quarter-turn")
df31$description <- paste("Surface area of the fenestra cochleae")
df32$description <- paste("Number of turns")

# output
# write.csv2(rbind(df24, df25, df26, df27, df28, df29, df30, df31, df32), file = "Data/traitsdf.csv")


