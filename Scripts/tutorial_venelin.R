# fancy package: venelin Tutorial

library(PCMBase)
library(PCMFit)
library(data.table)
library(ggplot2)
library(ggtree)
library(cowplot)

# make results reproducible
set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")

tree <- PCMTree(PCMFitDemoObjects$dtSimulated$tree[[1]])
tree

X <- PCMFitDemoObjects$dtSimulated$X[[1]][, seq_len(PCMTreeNumTips(tree))]
dim(X)

plTree <- PCMTreePlot(tree, layout="fan") +
  geom_tiplab2(size = 2) + 
  geom_nodelab(size = 2, color = "black") + 
  geom_treescale(width = max(PCMTreeNodeTimes(tree)), x = 0, linesize = .25, fontsize = 2, offset = 79)

plX <- PCMPlotTraitData2D(
  X[, seq_len(PCMTreeNumTips(tree))], 
  tree, 
  scaleSizeWithTime = FALSE,
  numTimeFacets = 4) +
  geom_text(
    aes(x = x, y = y, label = id, color = regime), 
    size=2, 
    position = position_jitter(.4, .4)) +
  theme_bw() +
  theme(legend.position = "bottom")

cowplot::plot_grid(plTree, plX, labels = LETTERS[1:2], nrow=2, rel_heights = c(2,1))

PCMDefaultModelTypes()

modelBM <- PCM(
  PCMDefaultModelTypes()["B"], modelTypes = PCMDefaultModelTypes(), k = 2)

modelBM

modelOU <- PCM(
  PCMDefaultModelTypes()["F"], modelTypes = PCMDefaultModelTypes(), k = 2)

modelOU

modelTrueTypeMapping <- MixedGaussian(
  k = 2,
  modelTypes = MGPMDefaultModelTypes(),
  mapping = c(4, 3, 2),
  X0 = structure(
    0, class = c("VectorParameter",
                 "_Global"), description = "trait values at the root"), 
  Sigmae_x = structure(
    0, class = c("MatrixParameter", "_Omitted",
                 "_Global"),
    description =
      "Upper triangular Choleski factor of the non-phylogenetic variance-covariance"))

treeWithTrueShifts <- PCMTree(PCMFitDemoObjects$dtSimulated$tree[[1]])
PCMTreeSetPartRegimes(
  treeWithTrueShifts, 
  part.regime = c(`81` = 1, `105` = 2, `125` = 3), 
  setPartition = TRUE)

modelTrue <- PCMFitDemoObjects$dtSimulated$model[[1]]
modelTrue

# calculate parameter count likelihood and AIC for it:
attr(modelTrue, "tree") <- treeWithTrueShifts
attr(modelTrue, "X") <- X
attr(modelTrue, "SE") <- X * 0.0


plTree <- PCMTreePlot(treeWithTrueShifts, layout="fan") %<+% 
  data.table(
    node = c(12, 77, 45), 
    part.model = c(" 1.D ", " 2.C ", " 3.B "),
    offset = 5) + 
  geom_tiplab2(size = 2) + 
  geom_tiplab2(aes(label = part.model), offset = 16) + 
  geom_nodelab(size = 2, color = "black") + 
  geom_treescale(
    width = max(PCMTreeNodeTimes(treeWithTrueShifts)), x = 0, 
    linesize = .25, fontsize = 2, offset = 79)

plX <- PCMPlotTraitData2D(
  X[, seq_len(PCMTreeNumTips(treeWithTrueShifts))], 
  treeWithTrueShifts, 
  scaleSizeWithTime = FALSE,
  numTimeFacets = 4) +
  geom_text(
    aes(x = x, y = y, label = id, color = regime), 
    size=2, 
    position = position_jitter(.4, .4)) +
  theme_bw() +
  theme(legend.position = "bottom")

cowplot::plot_grid(plTree, plX, labels = LETTERS[1:2], nrow = 2, rel_heights = c(2,1))


## File: DefineParameterLimits.R
# lower limits for models A and B
PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$Sigma_x)) {
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -.0
    }
  } else {
    if(!is.Diagonal(o$Sigma_x)) {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 2, r] <- -.0
      }
    }
  }
  o
}

# upper limits for models A and B
PCMParamUpperLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 1.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- 1.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 1.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- 1.0
      }
    }
  }
  o
}

# lower limits for models C, ..., F.
PCMParamLowerLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$Theta)) {
    o$Theta[1] <- 0.0
    o$Theta[2] <- -1.2
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 0.0
      o$Theta[2, r] <- -1.2
    }
  }
  if(is.Global(o$Sigma_x)) {
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -.0
    }
  } else {
    if(!is.Diagonal(o$Sigma_x)) {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 2, r] <- -.0
      }
    }
  }
  o
}

# upper limits for models C, ..., F.
PCMParamUpperLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$Theta)) {
    o$Theta[1] <- 7.8
    o$Theta[2] <- 4.2
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 7.8
      o$Theta[2, r] <- 4.2
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 1.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- 1.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 1.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- 1.0
      }
    }
  }
  o
}

###################################################
###################################################
# warnings
fitBM <- PCMFit(model = modelBM, tree = tree, X = X, 
                metaI = PCMBaseCpp::PCMInfoCpp)
fitOU <- PCMFit(model = modelOU, tree = tree, X = X, 
                metaI = PCMBaseCpp::PCMInfoCpp) 

# This can take about 5 minutes to finish
fitMGPMTrueTypeMapping <- PCMFit(
  model = modelTrueTypeMapping, 
  tree = treeWithTrueShifts, 
  X = X,
  metaI = PCMBaseCpp::PCMInfoCpp)

fitMGPMTrueTypeMappingCheat <- PCMFit(
  model = modelTrueTypeMapping, 
  tree = treeWithTrueShifts, 
  X = X,
  matParInit = matrix(PCMParamGetShortVector(modelTrue), 1L),
  numRunifInitVecParams = 1000L,
  numGuessInitVecParams = 100L,
  metaI = PCMBaseCpp::PCMInfoCpp)

listModels <- list(
  RetrieveBestModel(fitBM), 
  RetrieveBestModel(fitOU),
  RetrieveBestModel(fitMGPMTrueTypeMapping), 
  RetrieveBestModel(fitMGPMTrueTypeMappingCheat), 
  modelTrue)

dtSummary <- data.table(
  model = c(
    "Global BM", 
    "Global OU", 
    "True MGPM, unknown parameters", 
    "True MGPM, known true parameters", 
    "True MGPM, true parameters"),
  p = sapply(listModels, PCMParamCount),
  logLik = sapply(listModels, logLik), 
  AIC = sapply(listModels, AIC))
knitr::kable(dtSummary)
################ very diff vals for True MGPM 1 and 2