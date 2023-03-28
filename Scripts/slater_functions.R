##### Functions used in Slater et al. Mysticete body size evolution Paper #####
##### 						compiled January 6 2017 					  #####

trend.shift <- function(phy, d, sem = NULL, model = c("trend.shift"), lb=NULL, ub=NULL, n.iter=50) {
  
  td <- treedata(phy, d, sort=T)
  phy <- td$phy
  d <- td$data[,1]
  
  ntax <- length(d)
  
  if(is.null(sem)) {
    sem <- setNames(rep(0, ntax), names(d))
  } else {
    sem <-treedata(phy, sem, sort=T)$data[,1]
  }
  
  v <- vcv(phy)
  tip.times <- diag(v)
  root.age <- max(tip.times)
  
  if(model=="trend.shift") {
    if(is.null(lb)) {
      lb = c(log(10^-8), -1, min(d),  10^-2)
    } 
    if(is.null(ub)){
      ub = c(log(10), 1, max(d), root.age)
    }
    bounds <- cbind(lb, ub)
    foo.unif <- function(x) runif(1, x[1], x[2])
    
    trend.shift <- function(x) {
      st <- root.age-x[4]
      trend.tips <- which(tip.times >= st)
      Ex <- rep(x[3], length(tip.times))
      Ex[trend.tips] <- x[3] + (x[2] *(tip.times[trend.tips] - st))
      V<-exp(x[1])*v
      diag(V) <- diag(V) + sem^2 
      -dmvnorm(x = d, mean = Ex, sigma =V , log=T)  
    }
    tmp <- list()
    for(i in 1:n.iter){
    start <- apply(bounds, 1, foo.unif)    
    start[3] <- phylogMean(v, d)
    tmp[[i]] <- optim(par=start, fn=trend.shift, method ="L-BFGS-B", lower=lb, upper=ub, hessian=T)
    }  
    o <- tmp[[which.min(unlist(lapply(tmp, foo <- function(x) return(x$value))))]]
    aic <- 2 * 4 - 2 * -o$value
    aicc <- 2 * 4 * (ntax - 1)/(ntax - 4 - 2) - 2 * -o$value
    hessian <- o$hessian
    rownames(hessian) <- colnames(hessian) <-    c("rate", "mu", "root", "shift")
      results <- list(lnl = -o$value, aic=aic, aicc=aicc, root = o$par[3], rate1=exp(o$par[1]),  mu=o$par[2], shift.time=o$par[4], convergence=o$convergence, hessian = hessian)
    
  }
 
  results
  
}



BranchingTimesFossil <- function (phy, tol = .Machine$double.eps^0.5) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")
    phy2 <- phy
    phy <- new2old.phylo(phy)
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    xx <- as.numeric(rep(NA, nb.tip + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(as.numeric(phy$edge[, 2]) == nod)
        base <- phy$edge[ind, 1]
        xx[i] <- xx[base] + phy$edge.length[ind]
    }
    bt <- abs(xx - max(xx));
	
	for(i in 1:length(bt)) {
		
		if(bt[i]<.Machine$double.eps^0.5) bt[i] <- 0; 	}
	
	names(bt) <- c(seq(nb.tip+1, nb.tip+nb.node), phy$tip.label)
	
	
	return(bt);
}

fossilDTT <- function (phy, data, index = c("avg.sq", "avg.manhattan", "num.states"), 
          mdi.range = c(0, 1), nsim = 0, CI = 0.95, plot = TRUE, calculateMDIp = F) 
{
  disp = match.arg(index, c("avg.sq", "avg.manhattan", "num.states"))
  td <- treedata(phy, data)
  dtt.data <- geiger:::.dtt(td$phy, td$data, disp = disp)
  ltt <- sort(BranchingTimesFossil(td$phy)[1:(td$phy$Nnode)], decreasing = TRUE)
  ltt <- c(0, (max(ltt) - ltt)/max(ltt))
  s <- ratematrix(td$phy, td$data)
  dtt.sims = NULL
  MDI = NULL
  ylim = c(range(pretty(dtt.data)))
  if (is.numeric(nsim)) {
    if (nsim > 0) {
      sims <- sim.char(td$phy, s, nsim)
      dtt.sims <- geiger:::.dtt(td$phy, sims)
      mean.sims <- apply(dtt.sims, 1, mean)
      median.sims <- apply(dtt.sims, 1, median)
      MDI <- unname(geiger:::.area.between.curves(ltt, apply(dtt.sims, 
                                                    1, median), dtt.data, sort(mdi.range)))
      names(MDI) = disp
      colnames(dtt.sims) = NULL
      yy = range(dtt.sims)
      ylim = range(c(ylim, yy))
    }
  }
  if (plot) {
    plot(ltt, dtt.data, xlab = "relative time", ylab = "disparity", 
         ylim = ylim, bty = "n", type = "n")
    if (!is.null(dtt.sims)) {
      poly = geiger:::.dtt.polygon(dtt.sims, ltt, alpha = 1 - CI)
      polygon(poly[, "x"], poly[, "y"], col = geiger:::.transparency("lightgray", 
                                                            0.5), border = NA)
      lines(ltt, median.sims, lty = 2)
    }
    lines(ltt, dtt.data, type = "l", lwd = 2)
  }
  res = list(dtt = dtt.data, times = ltt, sim = dtt.sims, MDI = MDI)
  drp = sapply(res, function(x) is.null(x))
  if (any(drp)) 
    res = res[-which(drp)]
  if (calculateMDIp) {
    pVal <- getMDIp(res)
    res <- c(res, MDIpVal = pVal)
  }
  return(res)
}


phylogMean <- function (phyvcv, data) 
{
    o <- rep(1, length(data))
    ci <- solve(phyvcv)
    m1 <- solve(t(o) %*% ci %*% o)
    m2 <- t(o) %*% ci %*% data
    return(m1 %*% m2)
}


biasedSamplingSim <- function(sig2, iter) {
##  sig2 is rate of morphological evolution, iter is number of starting positions used for ML inference of trend shift
  while(1) {
    phy <- pbtree(b = 0.2, d = 0.15, n=15) # condition on 15 extant taxa 
    if(length(phy$tip.label) > 100 & length(which(round(max(diag(vcv(phy))) - diag(vcv(phy)),2)==0)) ==15 ) break # require more than 60 fossils
  }
  #write.tree(phy, "~/Desktop/simtree.tre")
  phy <- rescale(phy, "depth", 35) # rescale to 35 Ma for age of mysticetes
  
  d <- fastBM(phy, sig2=sig2)
  fossils <- setdiff(phy$tip.label, drop.fossil(phy)$tip.label)
  extant <- drop.fossil(phy)$tip.label 
  
  ## get params needed for logistic sampling probs ##
  z <-  (d-min(d)) / (max(d)-min(d))

  scale <- c(0, 2.5, 5, 10, 25, 100)

  sprob <- function(x,c) {
    1-  (1 / (1+exp(c*(0.5-x))))
  }

  sampling.prob <- matrix(data=NA, nrow=length(d), ncol=length(scale))
  for(i in 1:length(scale)) {
      sampling.prob[,i] <- sprob(z, scale[i])
  }
  sampling.prob<- sampling.prob[match(fossils, phy$tip.label), ]
  #plot(seq(0, 1, 0.1), apply(as.matrix(seq(0,1,0.1)), 1, sprob, c=25 ))
  
  ## generate sampled datasets and store in list ##
  sampled.data <- list()
  for(s in 1: length(scale)) {
    sampled <- character()
    for(f in 1:length(fossils)) {
      p <- runif(1,0,1)
      if( p <= sampling.prob[f, s] ) {
        sampled <- c(sampled, fossils[f])
      } 
    }
    sampled.data[[s]] <- c(d[sampled], d[extant])    
  }
  
  bm <- fitContinuous(phy, d, model="BM")$opt$aicc
  ts <- trend.shift(phy, d, model="trend.shift", n.iter = iter)
  res <- c(aicw(c(bm, ts$aic))[,3], ts$mu, ts$shift.time)

  for(m in 1: length(sampled.data)) {
    bm <- suppressWarnings(fitContinuous(phy, sampled.data[[m]], model="BM")$opt$aicc)
    ts <- suppressWarnings(trend.shift(phy, sampled.data[[m]], model="trend.shift", n.iter = iter))
    res <- c(res, c(aicw(c(bm, ts$aicc))[,3], ts$mu, ts$shift.time))
    
  }
  x <- c("BM", "TS", "MU", "shift")
  
  names(res) <- c(paste(x, "all",sep="_"),paste(x, 0,sep="_"), paste(x, 2.5,sep="_"),paste(x, 5,sep="_"),paste(x, 10,sep="_"),paste(x, 25,sep="_"), paste(x, 100,sep="_"))
  
  res
}



###### trend shift simulation ####


trend.shift.sim <- function(vcv, root.state, sigmasq, mu, shift, n) {
  
  root.age <- max(vcv)
  tip.times <- root.age - diag(vcv)
  v <- vcv*sigmasq
  ex <- rep(root.state, length(tip.times))
  ex[which(tip.times<=shift)] <- root.state + (mu*(shift-tip.times[which(tip.times<=shift)]) )
  simd <- t(mvrnorm(n=n, mu=ex, Sigma=v))
  res <- array(simd, dim=c(1,length(ex), n))
  res <- aperm(res, c(2,1,3))
  rownames(res) <- phy$tip.label
  return(res)
}

