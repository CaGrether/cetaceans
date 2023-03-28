library(ape)
library(geiger)
library(phytools)
library(mvtnorm)
library(ggplot2)
library(OUwie)
library(RPANDA) # use github version for fit_t_env() as of January 2017
setwd("C:/Users/Caro/Documents/Uni_Erlangen/Studium/Thesis/Scripts/Slater Scripts")
source("slater_functions.R")

## read in tree and data

phy <- ladderize(read.nexus("slater_mcct.tre")) 

d <- read.csv("slater_length_data.csv", stringsAsFactors = F, row.names=1)

td <- treedata(phy, d, sort=T)

phy <- td$phy

logTL <-setNames(td$data[,1], rownames(td$data))

sem=setNames(d[match(names(logTL), rownames(d)),2], names(logTL))

rm(d)

### if NOT truncating tips to create Pleistocene fossil record, comment out the following section ###

 phy$edge.length[match(grep("robustus", phy$tip.label), phy$edge[,2])] <-phy$edge.length[match(grep("robustus", phy$tip.label), phy$edge[,2])] - 0.2

 phy$edge.length[match(grep("Megaptera_novaeangliae", phy$tip.label), phy$edge[,2])] <-phy$edge.length[match(grep("Megaptera_novaeangliae", phy$tip.label), phy$edge[,2])] - 0.125

###############################################################################################

 ### to add Horopeta and Whakakai (Tsai and Fordyce, 2015, 2016), uncomment the following block
 
 # phy<- bind.tip(phy, tip.label="Horopeta_umarere", edge.length = 2, where = 82, position = 3.5)
 # 
 # phy<- bind.tip(phy, tip.label="Whakakai_waipata", edge.length = 1.5, where = which(phy$tip.label=="Horopeta_umarere"), position = 1.5)
 # 
 # plot(phy, cex=0.5) # check they're in the correct place, sister to Mauicetus + Aglaocetus + crown mysticetes and one node up from Eomysticetus
 # 
 # logTL <- c(logTL, setNames(rep(log(1000, base=10),2), c("Horopeta_umarere", "Whakakai_waipata")))
 # 
 # sem <-c(sem, setNames(rep(0.068110790, 2), c("Horopeta_umarere", "Whakakai_waipata")))
 # 
 ###############################################################################################
 
 

## model fitting 

## basic models ##

bm <- fitContinuous(phy, logTL, model="BM", SE = sem)
acdc <- fitContinuous(phy, logTL, model= "EB", bounds=list(a=c(-1, 1)), SE = sem)
(ou <- fitContinuous(phy, logTL, model="OU",SE = sem)) ##
trend <- fitContinuous(phy, logTL, model="drift",SE = sem)

# temp-dependent rates

zd <- read.csv("2008CompilationData.csv") #d18o data
zd <- zd[-which(zd$age>40), ] 
zd <- zd[-which(is.na(zd$oxygen)),]

spl <- smooth.spline(-zd$age, zd$oxygen, df = 15) 
xx <- cbind(-spl$x, spl$y)
temp.rate <- fit_t_env(phy, logTL, xx, error=sem, model= "EnvLin", 
          method="Nelder-Mead", control=list(maxit=20000))

## trend - shift model ##
ts <- trend.shift(phy, logTL, sem = sem, model="trend.shift", n.iter=40)

#####  use OUwie to test for a shift in rate through time  #####

ouwie.states <- logTL
ouwie.data <- data.frame(species=names(logTL), trait=logTL, sem = sem)

slice1bms <- OUwie.slice(phy, ouwie.data,model=c("BMS"), root.age = max(diag(vcv(phy))), root.station=FALSE, timeslices=c(NA),mserr="known")

## compute AICc weights

all.model.fit <- aicw(setNames(c(bm$opt$aicc, acdc$opt$aicc, ou$opt$aicc, trend$opt$aicc,  ts$aicc), c("BM", "ACDC","SSP", "Trend",   "BM to \nTrend")))
#slice1bms$AICc,"Rate Shift", temp.rate$aicc,"Temp-dep \nRates",
all.model.fit[rev(order(all.model.fit$w)),]

#### model parameters and fits ####

sigsq <- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$sigsq))
root <- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$z0))
aicc <- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$aicc))
lnL<- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$lnL))
param1 <- c(NA, acdc$opt$a, trend$opt$drift, ou$opt$alpha)
shift.time <- rep(NA, length(param1))
k <- unlist(lapply(list(bm$opt, acdc$opt, trend$opt, ou$opt), foo <- function(x) x$k))

res <- cbind(sigsq, root, param1, shift.time, lnL, k, aicc)
rownames(res) <- c("BM", "ACDC", "trend", "SSP")
res <- rbind(res, c(temp.rate$param[1], temp.rate$root, temp.rate$param[2], NA, temp.rate$LH,temp.rate$free.parameters, temp.rate$aicc))

rownames(res) <- c("BM", "ACDC", "trend", "SSP", "TTR")

bms <- c(slice1bms$solution[2,1], slice1bms$theta[1,1], slice1bms$solution[2,2], slice1bms$timeslices[2], slice1bms$loglik, 4, slice1bms$AICc)
tshift <- c(ts$rate1, ts$root, ts$mu, ts$shift.time, ts$lnl, 4, ts$aicc)

res <- rbind(res, bms, tshift) 

res # <- final results


#### tests of phylogenetic signal #####

psk <- phylosig(phy, logTL, method="K", test=TRUE)
psl <- phylosig(phy, logTL, method="lambda", test=TRUE)

dtt_fossil <- fossilDTT(phy, logTL, nsim = 10000)
geiger:::getMDIp(dtt_fossil)
cbind(dtt_fossil$times, dtt_fossil$dtt)
root.age <- max(diag(vcv(phy)))

quartz(width=7, height=3)
par(mar=c(4,5,1,1), mfrow=c(1,2))
plot((root.age* dtt_fossil$times), dtt_fossil$dtt, type="n", axes=FALSE, xlab="", ylab=paste("mean relative \nsubclade disparity"))


axis(1, at=(root.age - c(35, 30, 25, 20, 15, 10, 5,0)), labels = c(35, 30, 25, 20, 15, 10, 5,0))
mtext(text="millions of years ago", side=1, line=2.5, at=root.age/2)
axis(2, las=2)
poly <- geiger:::.dtt.polygon(dtt_fossil$sim, (root.age* dtt_fossil$times), alpha = 1 - 0.95)
polygon(poly[, "x"], poly[, "y"], col = geiger:::.transparency("lightgray", 0.5), border = NA)
lines((root.age* dtt_fossil$times), apply(dtt_fossil$sim, 1, median), lty = 2, col="black")
lines((root.age* dtt_fossil$times), dtt_fossil$dtt, col="blue", lwd=1.5)

par(xpd=NA)
text(root.age-50, 1.05, label="(a)", font=1, cex=1.5)
text(root.age+5, 1.05, label="(b)", font=1, cex=1.5)
par(xpd=F)

hist(dtt_fossil$sim[68,], main="", xlab="mean relative subclade disparity\n at 5.01 Ma", col='gray', breaks=100, border="gray", xlim=c(0,2.5))
arrows(dtt_fossil$dtt[68],y1=500, y0=0, col="red", angle=20, code=1, lwd=2)

quantile(dtt_fossil$sim[68,], c(0.95))

(root.age - (root.age*dtt_fossil$times))[68]

text(x=dtt_fossil$dtt[68], y= 500, paste("MDI = ", round(dtt_fossil$dtt[68], 3),"\n", "P = ", 
                                         round(length(which(dtt_fossil$sim[68,] > dtt_fossil$dtt[68])) / ncol(dtt_fossil$sim), 3), sep=""), pos=4)




