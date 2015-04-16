
##' @title          myCols
##' @description    Function that defines nice colors
##'
##' @return
##' A vector of colors
##'
##' @author Dorothee Nickles
myCols <- function(){
  myColsSet1 <- c("antiquewhite3", "aquamarine3", "azure3", "lightpink1", "cadetblue4", "coral2",
                  "firebrick", "lightblue", "olivedrab", "palevioletred3")
  myColsSet2 <- c("darkgoldenrod", "darkgoldenrod2", "cornflowerblue", "darkkhaki", "coral3",
                  "chartreuse3", "dodgerblue4", "blueviolet", "darkturquoise", "thistle4") 
  myColsSet3 <- c("darkgreen", "darkmagenta", "cyan4", "darkorange", "darkred", "gold",
                  "darkblue", "deeppink2", "yellowgreen", "orangered3")
  c(myColsSet3, myColsSet2, myColsSet1)
}

##' @title          dCoVar
##' @description    Determine standard deviation at a certain percent CV
##'
##' @param  x       numeric vector 
##' @param  d       scalar numeric, percent CV
##' @param  ...     additional parameters passed on to mean()
##'
##' @return standard deviation of x at d percent of CV
##'
##' @author Dorothee Nickles
##'
dCoVar <- function(x, d, ...) {
  d * mean(x, ...)
}

##' @title          cutoffByVar
##' @description    Determine cutoffs of x (for outlier detection) based on a certain percent CV
##'
##' @param  x       numeric vector 
##' @param  d       scalar numeric, percent CV; passed on to dCoVar
##' @param  ...     additional parameters passed on to mean()
##'
##' @return
##' A list of length 2, with a scalar numeric in each slot, one giving the lower threshold (mean(x) - CV% based cutoff),
##' the other giving the upper threshold (mean(x) + percent CV based cutoff) for outlier definition.
##'
##' @author Dorothee Nickles
##'
cutoffByVar <- function(x, d, ...) {
  thresholds <- list(mean(x, ...) - dCoVar(x, d, ...), mean(x, ...) + dCoVar(x, d, ...))
  return(thresholds)
}

##' @title          cutoffByMMAD
##' @description    Determine cutoffs of x (for outlier detection) based on median
##'
##' @param  x       numeric vector 
##' @param  d       scalar numeric, factor by which to multiply MAD of x
##' @param  ...     additional parameters passed on to median()
##'
##' @return
##' A list of length 2, with a scalar numeric in each slot, one giving the lower threshold (median(x) - d * mad(x)),
##' the other giving the upper threshold (median(x) + d * mad(x)) for outlier definition.
##'
##' @author Dorothee Nickles
##'
cutoffByMMAD <- function(x, d, ...) {
  thresholds <- list(median(x, ...) - d * mad(x, ...), median(x, ...) + d * mad(x, ...))
  return(thresholds)
}

##' @title          colByFun
##' @description    Color x based on upper and lower thresholds
##'
##' @param  x           Numeric vector 
##' @param  thresholds  List of length 2, with a scalar numeric in each slot, one giving the lower the upper threshold (for outlier definition)
##'
##' @return
##' A vector of colors, with "red" for all values of x exceeding thresholds and "black" for all other values
##'
##' @author Dorothee Nickles
##'
colByFun <- function(x, thresholds) {
  ifelse(x < thresholds[[1]], "red", 
      ifelse(x >  thresholds[[2]], "red", "black"))
}

##' @title          colByCovar
##' @description    Define colors based on a covariate of a NanoString ExpressionSet object
##'
##' @param  pdata   pData() of a NanoString ExpressionSet object
##' @param  covar   character, colname in pData(rccSet); used to stratify (color) data
##'
##' @return
##' A list of length 2, with [["color"]] being a character vector of colors (one color for each level of covar) of length=number of observations
##' and [["legend"]] providing the levels of covar to map colors to covar
##'
##' @author Dorothee Nickles
##'
colByCovar <- function(pdata, covar) {
  if(!covar %in% colnames(pdata)) {
    stop("The covariate you chose is not defined.")
  } else {
    diffTypes <- list()
    diffTypes[["color"]] <- myCols()[as.numeric(as.vector(factor(pdata[, covar], 
            labels=c(1:length(unique(pdata[, covar]))))))]
    diffTypes[["legend"]] <- unique(pdata[, covar])
    return(diffTypes)
  }         
}

##' @title fovPlot
##'
##' @description
##' Plot the fraction of successfully imaged fields of view (FOV) in the given NanoString ExpressionSet.
##' The ExpressionSet's pData should have 'FovCount' and 'FovCounted' columns populated with the total
##' and successfully imaged FOV counts, respectively. Samples with a FOV counted/FOV count of less than
##' 80% are marked in red (dashed red line indicates threshold).
##'
##' @param  rccSet  NanoString ExpressionSet object
##'
##' @return
##' A plot
##'
##' @author Dorothee Nickles
##'
fovPlot <- function(rccSet) {
  stopifnot(is(rccSet, "ExpressionSet"))
  fov <- as.numeric(as.vector(pData(rccSet)$FovCounted))/
         as.numeric(as.vector(pData(rccSet)$FovCount))
  plot(fov, pch=16, main="FOV Counted versus FOV Count", ylab="FOV Counted/FOV Count",
       xlab="Sample Index",
       col=ifelse(fov < 0.8, "red", "black"),
       ylim=c(min(fov, 0.8), 1))
  abline(h=0.8, col="red", lty=2)	
}

##' @title bdPlot
##'
##' @description
##' Plot the binding density of each sample in a NanoString ExpressionSet object.
##' Samples with a binding density < 0.05 or > 2.25 (thresholds defined by NanoString) 
##' are marked in red (dashed red line indicates threshold).
##'
##' @param  rccSet  NanoString ExpressionSet object
##'
##' @return
##' A plot
##'
##' @author Dorothee Nickles
##'
bdPlot <- function(rccSet) {
  stopifnot(is(rccSet, "ExpressionSet"))
  plot(as.numeric(as.vector(pData(rccSet)$BindingDensity)), pch=16, 
       ylim=c(0, max(c(as.numeric(as.vector(pData(rccSet)$BindingDensity)), 2.25))),
       main="Binding Density", ylab="Binding density",
       xlab="Sample Index",
       col=ifelse(as.numeric(as.vector(pData(rccSet)$BindingDensity)) > 2.25, "red", 
            ifelse(as.numeric(as.vector(pData(rccSet)$BindingDensity)) < 0.05, "red", "black")))
  abline(h=c(0.05, 2.25), col="red", lty=2)
}

##' @title          flagSamplesTech
##' @description    Flag samples based on their technical performance, i.e. field of vision (FOV) counted and binding density
##'
##' @details
##' Samples with a FOV counted/FOV count of less than 80% will be flagged. Also samples with a binding density
##'
##' @param rccSet NanoString ExpressionSet object
##'
##' @return
##' A numeric vector giving the indices of samples with outlier values in FOV counted and
##' binding density < 0.05 or > 2.25 (thresholds defined by NanoString) will be flagged.
##'
##' @author Dorothee Nickles
##'
flagSamplesTech <- function(rccSet) {
  return(unique(c(which(as.numeric(as.vector(pData(rccSet)$FovCounted))/as.numeric(as.vector(pData(rccSet)$FovCount)) < 0.8),
                  which(as.numeric(as.vector(pData(rccSet)$BindingDensity)) > 2.25 | as.vector(pData(rccSet)$BindingDensity) < 0.05))))
}

##' @title          ctrlsOverviewPlot
##' @description    Plot individual negative and positive controls across all samples in a NanoString ExpressionSet object.
##'
##' @param  rccSet  NanoString ExpressionSet object
##'
##' @return
##' A plot with two panels, one for the negative controls, one for the positive controls.
##'
##' @author Dorothee Nickles
##'
ctrlsOverviewPlot <- function(rccSet) {
  stopifnot(is(rccSet, "ExpressionSet"))
  par(mfrow=c(1,2))

  M <- assayData(rccSet)$exprs

  negCtrls <- (fData(rccSet)$CodeClass == "Negative")
  M.negCtrls <- M[ negCtrls, ]
  rownames(M.negCtrls) <- fData(rccSet)$GeneName[negCtrls]

  boxplot(t(M.negCtrls),
          las=2,
          cex.axis=0.6,
          ylab="counts (raw)",
          log="y",
          main="Negative Controls")

  posCtrls <- which(fData(rccSet)$CodeClass == "Positive")
  M.posCtrls <- M[ posCtrls, ]
  rownames(M.posCtrls) <- fData(rccSet)$GeneName[posCtrls]
  spikein <- fData(rccSet)$SpikeInInput[ posCtrls ]
  M.posCtrls.ordered <- M.posCtrls[ order(spikein), ]

  boxplot(t(M.posCtrls.ordered),
          las=2,
          cex.axis=0.6, 
          ylab="counts (raw)",
          log="y",
          main="Positive Controls")
}

##' @title          negCtrlsPlot
##' @description    Plot negative controls across all samples in a NanoString ExpressionSet object
##'
##' @details
##' In the second panel, boxplots are colored by lane (as specified in the pData slot). Bars on top of
##' the panel indicate the stage position for each cartridge/sample (as specified in the pData slot).
##'
##' @param  rccSet  NanoString ExpressionSet object
##'
##' @return
##' A plot with two panels: one showing boxplots for the individual negative controls across all samples, and 
##' one showing boxplots for the negative control counts for each individual sample (lane-specific background).
##'
##' @author Dorothee Nickles
##'
negCtrlsPlot <- function(rccSet) {
  stopifnot(is(rccSet, "ExpressionSet"))
  
  negCtrls <- (fData(rccSet)$CodeClass == "Negative")
  M <- assayData(rccSet)$exprs

  par(mfrow=c(1,2))

  boxplot(t(M[ negCtrls, ]),
          las=2,
          cex.axis=0.6,
          ylab="counts (raw)",
          log="y",
          main="Negative Controls")

  boxplot(M[ negCtrls, ],
          log="y",
          xaxt="n",
          col=myCols()[ pData(rccSet)$LaneID ],
          ylab="counts",
          xlab="lane")

  axis(1,
       at = 1:ncol(M),
       labels = as.integer(pData(rccSet)$LaneID),
       cex.axis=0.6)

  n <- 1
  for (i in 1:length(unique(pData(rccSet)$StagePosition))) {
    tmp <- which(pData(rccSet)$StagePosition == pData(rccSet)$StagePosition[i])
    lines(y = rep(max(M[ negCtrls, ]) + 5, length(tmp)),
          x = n:(n+length(tmp)-1),
          col=i,
          lwd=4,
          xpd=TRUE)
    n <- n + length(tmp)
  }
}

##' @title          negCtrlsByLane
##' @description    Plot negative controls per lane in a NanoString ExpressionSet object
##'
##' @details
##' Boxplots are colored by lane (as specified in the pData slot). Bars on top of
##' the panel indicate the stage position for each cartridge/sample (as specified in the pData slot).
##'
##' @param rccSet NanoString ExpressionSet object
##'
##' @return
##' A plot with boxplots for the negative control counts for each individual sample (lane-specific background)
##'
##' @author Dorothee Nickles
##'
negCtrlsByLane <- function(rccSet) {
  stopifnot(is(rccSet, "ExpressionSet"))

  negCtrls <- (fData(rccSet)$CodeClass == "Negative")

  M <- assayData(rccSet)$exprs
  laneid <- as.numeric(pData(rccSet)$LaneID)

  boxplot(
     M[ negCtrls, ]
    ,log="y"
    ,xaxt="n"
    ,col=myCols()[laneid]
    ,ylab="counts"
    ,xlab="lane"
    ,main="Negative Controls by lane -\nlane-specific background"
  )
  axis(1,
       at = 1:ncol(M),
       labels = laneid,
       cex.axis = min(0.2 + 1/log2(length(laneid)), 0.8), las=2)
  ## disable addition of StagePosition indicator for now
  #n <- 1
  #for (i in 1:length(unique(pData(rccSet)$StagePosition))) {
  #  tmp <- which(pData(rccSet)$StagePosition == pData(rccSet)$StagePosition[i])
  #  lines(y = rep(max(M[ negCtrls, ]) + 5, length(tmp)),
  #        x = n:(n+length(tmp)-1),
  #        col = i,
  #        lwd = 4,
  #        xpd = TRUE)
  #  n <- n + length(tmp)
  #}
  legend(x=(length(laneid)+3), y=max(M[ negCtrls, ]), cex=0.5,
          xpd=TRUE, bty="n", fill=myCols()[sort(unique(laneid))],
          legend=paste(rep("Lane", length(unique(laneid))), sort(unique(laneid))))
}

##' @title          scatterPair
##' @description    Helper function for a scatter plot inside a pairs plots
##'
##' @param  x   integer, x positions
##' @param  y   integer, y positions
##'
##' @return
##' A scatter plot x versus y.
##'
scatterPair <- function(x,y) {points(x,y, pch=20)}

##' @title          panelCor
##' @description    Helper function for printing correlation coefficients inside a pairs plots
##'
##' @param  x       integer
##' @param  y       integer, same length as x
##' @param  digits  scalar integer, indicating the number of decimal positions for displaying the correlation coefficient
##' @param  cex.cor scalar numeric to specify relative font size for priting the correlation coefficient
##' @param  doTest  boolean, whether a results of cor.test should be displayed as well
##'
##' @return
##' Prints correlation coefficients (and p-values if doTest = TRUE) within a pairs plot.
##'
panelCor <- function(x, y, digits=2, cex.cor=0.75, doTest=FALSE)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  if (doTest) {
    test <- cor.test(x,y)
    Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
    text(0.5, 0.25, paste("r=",txt))
    text(.5, .75, Signif)
  } else {
    text(0.5, 0.5, paste("r=",txt), cex=cex.cor)
  }
}

##' @title          negCtrlsPairs
##' @description    Pairs plot of negative controls across all samples in a NanoString ExpressionSet object
##'
##' @param  rccSet          NanoString ExpressionSet object
##' @param  log.transform   boolean, whether data needs to be log2 transformed
##'
##' @return
##' Pairs plot of the negative controls with a scatter plot in the lower panel and
##' correlatation coefficients printed in the upper panel.
##'
##' @author Dorothee Nickles
##'
negCtrlsPairs <- function(rccSet, log.transform=FALSE) {

  M <- assayData(rccSet)$exprs
  negCtrls <- fData(rccSet)$CodeClass == "Negative"

  tmp <- M[ negCtrls, ]
  if (log.transform) {
    tmp <- log2(tmp)
  }

  pairs(t(tmp),
        upper.panel=panelCor,
        lower.panel=scatterPair)
}

##' @title          posNormFactPlot
##' @description    Plot positive control scaling factor for each sample in a NanoString ExpressionSet object
##'
##' @param  rccSet  NanoString ExpressionSet object
##'
##' @return
##' A plot of the positive control scaling factor for each sample in a NanoString ExpressionSet object.
##' Samples with a positive control scaling factor < 0.3 or > 3 (thresholds defined by NanoString) are
##' marked in red (dashed red line indicates threshold).
##'
##' @author Dorothee Nickles
##'
posNormFactPlot <- function(rccSet) {
  stopifnot(is(rccSet, "ExpressionSet"))

  M <- assayData(rccSet)$exprs
  positives <- (fData(rccSet)$CodeClass == "Positive")
  
  posSignal <- colSums( M[positives, ] )
  signalMean <- mean(posSignal)
  posFactor <- signalMean/posSignal
  plot(posFactor,
       las = 2,
       cex.axis = 0.6,
       ylab = "pos ctrl scaling factor",
       pch = 16,
       xlab = "Sample index",
       main = "Positive control scaling factors across samples", 
       col = ifelse(posFactor < 0.3, "red", ifelse(posFactor > 3, "red", "black")),
       ylim = c(min(c(posFactor, 0.3)), max(c(posFactor, 3))))
  abline(h = c(0.3, 3),
         col = "red",
         lty = 2)
}

##' @title          zfacFun
##' @description    Calculate Z' Factor
##'
##' @param  p   numeric vector: measurements for the positive controls (or actual measurement)
##' @param  n   numeric vector: measurements for the negative controls
##'
##' @return     Scalar numeric: the Z' Factor
##'
##' @author Dorothee Nickles
##'
zfacFun <- function(p, n) {
  return( 1 - ((3* (mad(p) + mad(n))) / abs(median(p) - median(n))) )
}

##' @title          ctrlsZprimePlot
##' @description    Plot distribution of counts and the Z' Factors comparing the nagative controls and
##'                 the three highest input positive controls of a NanoString ExpressionSet object.
##'
##' @param  rccSet  NanoString ExpressionSet object
##'
##' @return         A plot
##'
##' @author Dorothee Nickles
##'
ctrlsZprimePlot <- function(rccSet) {
  stopifnot(is(rccSet, "ExpressionSet"))

  M <- assayData(rccSet)$exprs
  negatives    <- (fData(rccSet)$CodeClass == "Negative")
  positive_128 <- (fData(rccSet)$SpikeInInput == 128)
  positive_32  <- (fData(rccSet)$SpikeInInput == 32)
  positive_8   <- (fData(rccSet)$SpikeInInput == 8)

  posN <- log2(apply( M[negatives, ], 2, median ))
  
  zfac1 <- round(zfacFun(log2(M[ positive_128, ]), posN), digits=2)
  zfac2 <- round(zfacFun(log2(M[ positive_32,  ]), posN), digits=2)
  zfac3 <- round(zfacFun(log2(M[ positive_8,   ]), posN), digits=2)
  plot(density(posN),
       col = myCols()[1],
       xlim = c(1, log2(max(M[ positive_128, ]))),
       main = "Frequency distributions of log2 transformed counts\nfor negative and positive controls")
  lines(density(log2(M[ positive_8,   ])), col=myCols()[2])        		
  lines(density(log2(M[ positive_32,  ])), col=myCols()[3])
  lines(density(log2(M[ positive_128, ])), col=myCols()[4])
  legend("top",
         lty = 1,
         col = myCols()[1:4],
         legend = c("negCtrls",
                    paste0("pos8: Z'=", zfac3),
                    paste0("pos32: Z'=", zfac2),
                    paste0("pos128: Z'=", zfac1)))
}

##' @title          iqrPlot
##' @description    Plot the interquartile range (IQR) for a certain code class of probes in a NanoString ExpressionSet object.
##'
##' @details
##' IQR of the specified code class for each sample in the ExpressionSet are plotted and outliers (as determined
##' by method) are marked in red (thresholds for outlier definition are plotted as red dashed lines).
##'
##' @param  rccSet      NanoString ExpressionSet object
##' @param  codeClass   Scalar character, specifying the code class (as annotated in the fData(rccSet)$CodeClass column)
##'                     for which the IQR shall be determined
##' @param  method      Scalar character, specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff      Scalar numeric, cutoff in method to determine outliers
##'
##' @return             A plot
##'
##' @seealso \code{\link{cutoffByMMAD}}, \code{\link{cutoffByVar}}
##'
##' @author Dorothee Nickles
##'
iqrPlot <- function(rccSet,
                    codeClass = c("Negative", "Positive", "Endogenous", "Housekeeping"),
                    method = c("cutoffByMMAD", "cutoffByVar"),
                    cutoff) {
  stopifnot(is(rccSet, "ExpressionSet"))
  codeClass <- match.arg(codeClass)
  method <- match.arg(method)
  fun <- match.fun(method)

  M <- assayData(rccSet)$exprs

  iqr <- apply(log2(M[ fData(rccSet)$CodeClass == codeClass, ]), 2, IQR)
  thresholds <- fun(iqr, cutoff)

  plot(iqr,
       las = 2,
       cex.axis = 0.6,
       ylab = "counts [log2]",
       pch = 16,
       xlab = "Sample index",
       main = sprintf("IQR of code class %s across samples", codeClass), 
       col = colByFun(iqr, thresholds), 
       ylim = c(min(c(iqr, unlist(thresholds))), max(c(iqr, unlist(thresholds)))))

  abline(h = unlist(thresholds),
         col = "red",
         lty = 2)
}

##' @title          posR2Plot
##' @description    Plot the R squared of linear fit of counts versus input for positive controls in a NanoString ExpressionSet object.
##'
##' @details
##' R squared for each sample in the ExpressionSet are plotted and samples with R squared < 0.95 are marked in red
##' (threshold indicated by dashed red line).
##'
##' @param  rccSet  NanoString ExpressionSet object
##'
##' @return         A plot
##'
##' @author Dorothee Nickles
##'
posR2Plot <- function(rccSet) {
  stopifnot(is(rccSet, "ExpressionSet"))
  par(mfrow=c(1,1))

  M <- assayData(rccSet)$exprs
  positives <- (fData(rccSet)$CodeClass == "Positive")

  tmp <- apply(log2(M[positives, ]),
               2,
               function(x) { summary(lm(x ~ log2(fData(rccSet)$SpikeInInput[ positives ])))$r.squared })
  plot(tmp,
       pch = 16,
       col = ifelse(tmp < 0.95, "red", "black"),
       main = "R2 for linearity in positive controls")

  abline(h=0.95, col="red", lty=2)
}

##' @title          posSlopePlot
##' @description    Plot the slope of linear fit of counts versus input for positive controls in a NanoString ExpressionSet object.
##'
##' @param  rccSet  NanoString ExpressionSet object
##' @param  method  Scalar character specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff  Scalar numeric: cutoff in method to determine outliers
##'
##' @return         A plot
##'
##' @details
##' The slope for each sample in the ExpressionSet are plotted and and outliers (as determined
##' by method) are marked in red (thresholds for outlier definition are plotted as red dashed lines).
##'
##' @author Dorothee Nickles
##'
posSlopePlot <- function(rccSet, method=c("cutoffByMMAD", "cutoffByVar"), cutoff) {
  stopifnot(is(rccSet, "ExpressionSet"))
  method <- match.arg(method)
  fun <- match.fun(method)

  M <- assayData(rccSet)$exprs
  positives <- (fData(rccSet)$CodeClass == "Positive")

  tmp <- apply(log2(M[positives, ]),
               2,
               function(x) { coefficients(lm(x ~ log2(fData(rccSet)$SpikeInInput[ positives ])))[2] })
  thresholds <- fun(tmp, cutoff)
  plot(tmp,
       pch = 16,
       col = colByFun(tmp, thresholds),
       main = "Slope in positive controls",
       ylab = "slope",
       xlab = "Sample index",
       ylim = c(min(c(tmp, unlist(thresholds))), max(c(tmp, unlist(thresholds)))))

  abline(h=unlist(thresholds), col="red", lty=2)
}

##' @title          posCtrlsPlot
##' @description    Wrapper to plot the individual positive controls across all samples as well as the slope of 
##'                 a linear fit of counts versus input for positive controls in a NanoString ExpressionSet object
##'
##' @param  rccSet  NanoString ExpressionSet object
##' @param  method  Scalar character specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff  Scalar numeric: cutoff in method to determine outliers
##'
##' @return         A plot with two panels
##'
##' @seealso \code{\link{posSlopePlot}}
##'
##' @author Dorothee Nickles
##'
posCtrlsPlot <- function(rccSet, method=c("cutoffByMMAD", "cutoffByVar"), cutoff) {
  stopifnot(is(rccSet, "ExpressionSet"))

  M <- assayData(rccSet)$exprs
  positives <- (fData(rccSet)$CodeClass == "Positive")

  M.pos <- M[positives, ]
  spikein <- fData(rccSet)$SpikeInInput[ positives ]
  M.pos.ordered <- M.pos[ order(spikein), ]

  par(mfrow=c(1,2))

  boxplot(t(M.pos.ordered),
          las = 2,
          cex.axis = 0.6,
          ylab = "counts (raw)",
          log = "y",
          main = "Positive Controls")

  posSlopePlot(rccSet, method, cutoff)

  par(mfrow=c(1,1))
}

##' @title          posMaxPlot
##' @description    Plot the maximum count of the positive controls for each sample in a NanoString ExpressionSet object.
##'
##' @details
##' The maximum count of the positive controls for each sample in the ExpressionSet are plotted and and outliers (as determined
##' by method) are marked in red (thresholds for outlier definition are plotted as red dashed lines).
##'
##' @param  rccSet  NanoString ExpressionSet object
##' @param  method  Scalar character, specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff  Scalar numeric, cutoff in method to determine outliers
##'
##' @return         A plot
##'
##' @author Dorothee Nickles
##'
posMaxPlot <- function(rccSet, method=c("cutoffByMMAD", "cutoffByVar"), cutoff) {
  stopifnot(is(rccSet, "ExpressionSet"))
  method <- match.arg(method)
  fun <- match.fun(method)

  M <- assayData(rccSet)$exprs
  positive_128 <- (fData(rccSet)$SpikeInInput == 128)

  thresholds <- fun(log2(M[positive_128, ]), cutoff)

  plot(y = M[positive_128, ],
       x = 1:ncol(M),
       las = 2,
       cex.axis = 0.6,
       ylab = "counts",
       main = "Max Positive Control across samples",
       pch = 16,
       xlab = "samples",
       col = colByFun(log2(M[positive_128, ]), thresholds),
       log = "y",
       ylim = c(min(c(M[positive_128, ], 2^unlist(thresholds))), 
                max(c(M[positive_128, ], 2^unlist(thresholds)))))

  abline(h=2^unlist(thresholds), col="red", lty=2)
}

##' @title          posRatioPlot
##' @description    Plot the ratio of the mean of positive control counts for each sample and the overall mean of positive control counts
##'                 in a NanoString ExpressionSet object.
##'
##' @details
##' The ratio for each sample in the ExpressionSet are plotted and and outliers (as determined
##' by method) are marked in red (thresholds for outlier definition are plotted as red dashed lines).
##'
##' @param  rccSet  NanoString ExpressionSet object
##' @param  method  Scalar character, specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff  Scalar numeric, cutoff in method to determine outliers
##'
##' @return         A plot
##'
##' @author Dorothee Nickles
##'
posRatioPlot <- function(rccSet, method=c("cutoffByMMAD", "cutoffByVar"), cutoff) {
  stopifnot(is(rccSet, "ExpressionSet"))
  method <- match.arg(method)
  fun <- match.fun(method)

  M <- assayData(rccSet)$exprs
  positives <- (fData(rccSet)$CodeClass == "Positive")

  meanPos <- log2(colMeans(M[positives, ]))
  meanPosMean <- mean(meanPos)
  meanPosRatio <- meanPos/meanPosMean
  thresholds <- fun(meanPosRatio, cutoff)

  plot(y = meanPosRatio,
       x = 1:ncol(M),
       las = 2,
       cex.axis = 0.6, 
       ylab = "mean pos ctrl counts/overall mean of pos ctrl counts",
       xlab="Sample index",
       main = "", 
       pch = 16,
       xlab = "samples",
       col = colByFun(meanPosRatio, thresholds),
       log = "y",
       ylim = c(min(c(meanPosRatio, unlist(thresholds))), max(c(meanPosRatio, unlist(thresholds)))))

  abline(h=unlist(thresholds), col="red", lty=2)
}

##' @title          allSumPlot
##' @description    Plot the sum of all counts (endogenous and housekeeping genes only) for each sample in a NanoString ExpressionSet object.
##'
##' @details
##' The sum of counts for each sample in the ExpressionSet are plotted and and outliers (as determined
##' by method) are marked in red (thresholds for outlier definition are plotted as red dashed lines).
##'
##' @param  rccSet      NanoString ExpressionSet object
##' @param  method      scalar character, specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff      scalar numeric, cutoff in method to determine outliers
##'
##' @return         A plot
##'
##' @author Dorothee Nickles
##'
allSumPlot <- function(rccSet,
                       method=c("cutoffByMMAD", "cutoffByVar"),
                       cutoff)
{
  stopifnot(is(rccSet, "ExpressionSet"))

  method <- match.arg(method)
  fun <- match.fun(method)

  M <- assayData(rccSet)$exprs
  nonctrls <- (fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping"))

  allSum <- colSums(M[nonctrls, ])
  thresholds <- fun(log2(allSum), cutoff)
  blankLabel <- getBlankLabel(rccSet)

  plot(y = allSum,
       x = 1:ncol(M),
       log = "y",
       las = 2,
       cex.axis = 0.6,
       ylab = "sum of raw counts", 
       xlab = "Sample index",
       main = "Sum of counts across samples", 
       col = colByFun(log2(allSum), thresholds),
       pch = ifelse( (pData(rccSet)$SampleType %in% blankLabel), 17, 16),
       ylim = c(min(c(allSum, 2^unlist(thresholds))), max(c(allSum, 2^unlist(thresholds)))))

  abline(h=2^unlist(thresholds), col="red", lty=2)
}

##' @title          sumRatioPlot
##' @description    Plot the ratio of counts for tho different code classes for all samples in a NanoString ExpressionSet object.
##'
##' @details
##' The ratio of counts is determined the following way: counts for code class == pair[[1]] / counts for code class == pair[[2]]. 
##' The ratio for each sample in the ExpressionSet is plotted and and outliers (as determined
##' by method) are marked in red (thresholds for outlier definition are plotted as red dashed lines).
##'
##' @param  rccSet      NanoString ExpressionSet object
##' @param  pair        list of length two, specifying the two code classes that are compared to each other
##' @param  method      scalar character, specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff      scalar numeric, cutoff in method to determine outliers
##'
##' @return         A plot
##'
##' @author Dorothee Nickles
##'
sumRatioPlot <- function(rccSet,
                         pair,
                         method=c("cutoffByMMAD", "cutoffByVar"),
                         cutoff)
{
  stopifnot(is(rccSet, "ExpressionSet"))

  method <- match.arg(method)
  fun <- match.fun(method)

  M <- assayData(rccSet)$exprs

  iSum <- colSums(M[fData(rccSet)$CodeClass %in% pair[[1]],])
  nSum <- colSums(M[fData(rccSet)$CodeClass %in% pair[[2]],])
  if (median(iSum/nSum) < 1) {
    cutoff <- cutoff*3					
  }
  thresholds <- fun(iSum/nSum, cutoff)
  blankLabel <- getBlankLabel(rccSet)

  plot(y = iSum/nSum,
       x = 1:ncol(M),
       las = 2,
       cex.axis = 0.6,
       ylab = "ratio",
       xlab = "Sample index",
       main = sprintf("Ratio of Sums - %s/%s", names(pair)[1], names(pair)[2]),
       col = colByFun(iSum/nSum, thresholds),
       pch = ifelse( (pData(rccSet)$SampleType %in% blankLabel), 17, 16),
       ylim = c(min(c(iSum/nSum, unlist(thresholds))), max(c(iSum/nSum, unlist(thresholds)))))

  abline(h=unlist(thresholds), col="red", lty=2)
}

##' @title          ctrlsPlot
##' @description    Wrapper function for plotting QC metrics involving NanoString controls
##'
##' @param  rccSet  NanoString ExpressionSet object
##' @param  method  scalar character, specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff  scalar numeric, cutoff in method to determine outliers
##'
##' @return
##' Plots generated by ctrlsOverviewPlot, iqrPlot, posSlopePlot, posMaxPlot, allSumPlot and sumRatioPlot.
##'
##' @seealso \code{\link{ctrlsOverviewPlot}}, \code{\link{iqrPlot}}, \code{\link{posSlopePlot}}, \code{\link{posMaxPlot}}, 
##' \code{\link{allSumPlot}}, \code{\link{sumRatioPlot}}
##'
##' @author Dorothee Nickles
##'
ctrlsPlot <- function(rccSet, method=c("cutoffByMMAD", "cutoffByVar"), cutoff) {
  par(mfrow=c(1,1))
  ctrlsOverviewPlot(rccSet)
  #negToPosPlot(rccSet)
  
  ## plot more on negatives - sample wise
  par(mfrow=c(2,1))
  iqrPlot(rccSet, "Negative", method, cutoff)
  iqrPlot(rccSet, "Positive", method, cutoff)
  
  ## plot more on positives - sample wise
  par(mfrow=c(2,1))
  posSlopePlot(rccSet, method, cutoff)
  posMaxPlot(rccSet, method, cutoff)
  
  ## plots on overall counts 
  par(mfrow=c(3,1))
  allSumPlot(rccSet, method, cutoff)
  sumRatioPlot(rccSet, pair=list(allSum=c("Endogenous", "Housekeeping"), posSum="Positive"), method, cutoff)
  sumRatioPlot(rccSet, pair=list(posSum="Positive", allSum=c("Endogenous", "Housekeeping")), method, cutoff)
}

##' @title          flagSamplesCtrl
##' @description    Flag samples based on the performance of their controls.
##'
##' @details
##' Outliers are determined based on method and cutoff. Samples defined as outliers in the following categories 
##' will be flagged: interquartile range of negative and positive controls, slope of the linear fit of count versus input 
##' of positive controls, a positive control scaling factor (based on sum of positive controls) below 0.3 or greater than 3.
##'
##' @param  rccSet  NanoString ExpressionSet object
##' @param  method  scalar character, specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff  scalar numeric, cutoff in method to determine outliers
##'
##' @return         A numeric vector giving the indices of samples with outlier values in control metrics
##'
##' @author Dorothee Nickles
##'
flagSamplesCtrl <- function(rccSet,
                            method = c("cutoffByMMAD", "cutoffByVar"),
                            cutoff) {
  stopifnot(is(rccSet, "ExpressionSet"))
  method <- match.arg(method)
  fun <- match.fun(method)

  M <- assayData(rccSet)$exprs

  negatives    <- (fData(rccSet)$CodeClass == "Negative")
  positives    <- (fData(rccSet)$CodeClass == "Positive")
  positive_128 <- (fData(rccSet)$SpikeInInput == 128)

  a <- apply(log2(M[negatives, ]), 2, IQR)
  ctrls <- which(colByFun(a, fun(a,cutoff)) == "red")

  a <- apply(log2(M[positives, ]), 2, IQR)
  ctrls <- c(ctrls, which(colByFun(a, fun(a,cutoff)) == "red"))

  #a <- log2(M[positive_128, ])
  #ctrls <- c(ctrls, which(colByFun(a, fun(a,cutoff)) == "red"))

  a <- apply(log2(M[positives, ]),
             2,
             function(x) { coefficients(lm(x ~ log2(fData(rccSet)$SpikeInInput[ positives ])))[2] })
  ctrls <- c(ctrls, which(colByFun(a, fun(a,cutoff)) == "red"))

  posSignal <- colSums(M[positives, ])
  signalMean <- mean(posSignal)
  posFactor <- signalMean/posSignal

  ctrls2 <- c(which(posFactor < 0.3), which(posFactor > 3))

  return(unique(ctrls, ctrls2))
}

##' @title          pcaPlot
##' @description    Wrapper function to perform a PCA analysis on the exprs slot of an ExpressionSet object
##'                 and plot some results
##'
##' @param  exx     exprs() of a NanoString ExpressionSet object (exprs(rccSet))
##' @param  ...     additional parameters passed on to the plotting functions
##'
##' @return
##' PCA screeplot and a plot with two panels, one plotting PC1 versus PC2, the other plotting
##' PC1 versus PC3.
##'
##' @author Dorothee Nickles
##'
pcaPlot <- function(exx, ...) {		## takes exprs of an expression set object and title for plots
  if (!is.matrix(exx)) {
    stop("Input is not a matrix")   # was: "stop("This function takes a matrix. Are you sure you put in exprs(rccSet)?")" -RZ 2015-03-24
  }
  pca_test <- prcomp(exx, scale=TRUE, tol=0)
  par(mfrow=c(1,1))
  screeplot(pca_test, ...)
  par(mfrow=c(1,2))
  plot(pca_test$rotation[,1], pca_test$rotation[,2], ylab="PC2", xlab="PC1", pch=16, ...)
  #text(x=pca_test$rotation[,1], y=pca_test$rotation[,2], label=rownames(pca_test$rotation),
  #  xpd=TRUE)
  plot(pca_test$rotation[,1], pca_test$rotation[,3], pch=16, ylab="PC3", xlab="PC1", ...)
  #text(x=pca_test$rotation[,1], y=pca_test$rotation[,3], label=rownames(pca_test$rotation),
  #  xpd=TRUE)
  #plot(pca_test$rotation[,2], pca_test$rotation[,3], pch=16, ylab="PC3", xlab="PC2", ...)
  #text(x=pca_test$rotation[,2], y=pca_test$rotation[,3], label=rownames(pca_test$rotation),
  #  xpd=TRUE)
  #plot(pca_test$rotation[,2], pca_test$rotation[,4], pch=16, ylab="PC4", xlab="PC2", ...)
}

##' @title          lodAssess
##' @description    Assess how many genes per sample in a NanoString ExpressionSet object are below the limit of detection.
##'
##' @param  rccSet      NanoString ExpressionSet object
##' @param  reference   scalar character, reference to determine limit of detection; either "blank" or "negatives"
##' @param  stringency  factor by which deviation (MAD) of the median of the reference counts is multiplied
##'
##' @return
##' A numeric vector giving the number of missing genes (endogenous and housekeeping genes) for each sample 
##' in a NanoString ExpressionSet.
##'
##' @author Dorothee Nickles
##'
lodAssess <- function(rccSet,
                      reference=c("blank", "negatives"),
                      stringency=2.5)
{
  stopifnot(is(rccSet, "ExpressionSet"))
  reference <- match.arg(reference)

  M <- assayData(rccSet)$exprs

  negatives <- (fData(rccSet)$CodeClass == "Negative")
  nonctrls  <- (fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping"))

  if (reference == "negatives") {
    negMed <- apply(M[negatives, ], 2, median)
    negMad <- apply(M[negatives, ], 2, mad)
    threshold <- negMed + negMad * stringency
    misValues <- vapply(1:ncol(M),
                        function(i) { sum(M[nonctrls, i] <= threshold[i]) },
                        integer(1))
    names(misValues) <- colnames(rccSet)
  } else {
    blankLabel <- getBlankLabel(rccSet)
    blanks <- (pData(rccSet)$SampleType %in% blankLabel)
    blankMed <- apply(as.matrix(M[, blanks]), 1, median, na.rm=TRUE)
    blankMad <- apply(as.matrix(M[, blanks]), 1, mad, na.rm=TRUE)
    threshold <- blankMed + blankMad * stringency
    threshold <- threshold[ nonctrls ]
    misValues <- apply(M[nonctrls, ], 
                       2,
                       function(x) { sum(x <= threshold) })
  }
  return(misValues)
}

##' @title          flagSamplesCount
##' @description    Flag samples based on overall counts
##'
##' @details
##' Outliers are determined based on method and cutoff. Samples defined as outliers in the following categories 
##' will be flagged: the sum of counts of endogeneous genes, the ratio of sums for positive controls to the sums of
##' endogenous genes, the number of genes below the lower limit of detection
##'
##' @param  rccSet      NanoString ExpressionSet object
##' @param  method      scalar character, specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar
##' @param  cutoff      scalar numeric, cutoff in method to determine outliers
##' @param  reference   scalar character, reference to determine limit of detection, passed on to \code{\link{lodAssess}}
##' @param  stringency  passed on to \code{\link{lodAssess}}
##' @param  maxMiss     scalar numeric, indicating the fraction of genes accepted to be missed 
##'
##' @return             A numeric vector giving the indices of samples with outlier values based on overall counts
##'
##' @author Dorothee Nickles
##'
flagSamplesCount <- function(rccSet,
                             method = c("cutoffByMMAD", "cutoffByVar"),
                             cutoff,
                             reference = c("blank", "negatives"),
                             stringency=2.5,
                             maxMiss) {
  stopifnot(is(rccSet, "ExpressionSet"))
  method <- match.arg(method)
  reference <- match.arg(reference)
  fun <- match.fun(method)

  M <- assayData(rccSet)$exprs

  positives  <- (fData(rccSet)$CodeClass == "Positive")
  endogenous <- (fData(rccSet)$CodeClass == "Endogenous")

  posSum <- colSums(M[positives, ])
  allSum <- colSums(M[endogenous, ])
  misValues <- lodAssess(rccSet, reference, stringency) 
  sampC <- unique(c(
    which(colByFun(allSum, fun(allSum, cutoff)) == "red"),
    which(colByFun(posSum/allSum, fun(posSum/allSum, cutoff*3)) == "red"),
    which(misValues/nrow(M) >= maxMiss)
  ))

  return(sampC)
}

##' @title          lodPlot
##' @description    Function to plot the number of missing genes per sample in a NanoString ExpressionSet object
##'
##' @details
##' Samples with more than 50% missingness are marked in red. If blank measurements are present, they are represented as triangles.
##'
##' @param  rccSet          NanoString ExpressionSet object
##' @param  reference       scalar character, reference to determine limit of detection; either "blank" or "negatives";
##'                         passed on to \code{\link{lodAssess}}
##' @param  stringency      factor by which deviation (MAD) of the median of the reference counts is multiplied;
##'                         passed on to \code{\link{lodAssess}}
##' @param  maxMiss         scalar numeric, indicating the fraction of genes accepted to be missed
##'
##' @return
##' A plot
##'
##' @author Dorothee Nickles
##'
lodPlot <- function(rccSet,
                    reference=c("blank", "negatives"),
                    stringency=2.5,
                    maxMiss)
{
  stopifnot(is(rccSet, "ExpressionSet"))

  blankLabel <- getBlankLabel(rccSet)
  nonctrls <- which(fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping"))

  misValues <- lodAssess(rccSet, reference, stringency)

  plot(misValues/length(nonctrls), 
       col = ifelse(misValues/length(nonctrls) >= maxMiss, "red", "black"),
       pch = ifelse( (pData(rccSet)$SampleType %in% blankLabel), 17, 16),
       main = "Sample Missingness",
       ylab = "Fraction of genes below background estimate",
       xlab="Sample index")

  abline(h=maxMiss, col="red", lty=2)   # was: abline(h=length(nonctrls)*maxMiss, col="red", lty=2) -RZ 2015-03-27
}

##' @title          sampleClustering
##' @description    Clustering by sample correlation
##'
##' @param  rccSet          ExpressionSet
##' @param  annCol          See \code{\link{aheatmap}}
##' @param  log.transform   (currently unused)
##' @param  covar           Covariate (e.g. "SampleType")
##' @param  ...             Extra arguments passed to aheatmap() (optional)
##'
##' @return
##' Nothing, a plot is generated on the active device.
##'
##' @importFrom NMF aheatmap
##' @importFrom RColorBrewer brewer.pal
##'
sampleClustering <- function(rccSet,                ## takes exprSet with raw and normalized counts, respectively
                             annCol=NULL,
                             log.transform=FALSE,   # TODO: check if this is still necessary -RZ 2015-03-10
                             covar,
                             ...)
{
  stopifnot(is(rccSet, "ExpressionSet"))

  M <- assayData(rccSet)$exprs

  zeroVarSamples <- which(apply(M, 2, function(x){ length(unique(x)) }) == 1)

  if (length(zeroVarSamples) > 0)
  {
    M.without_zeroVarSamples <- M[, -zeroVarSamples]
    pdata <- pData(rccSet)[ -zeroVarSamples ]

    warning(sprintf("%i samples were excluded from correlation plot because of zero variance:\n %s", 
                    length(zeroVarSamples),
                    paste(rownames(pData(rccSet))[zeroVarSamples], collapse=",")))
  } else {
    M.without_zeroVarSamples <- M
    pdata <- pData(rccSet)
  }

  cols <- colByCovar(pdata, covar)
  if (missing(annCol)) {
    annCol <- data.frame( SampleSubgroup = pdata[, covar])
  }
  par(mfrow=c(1,1))
  aheatmap(cor(M.without_zeroVarSamples),
           annCol = annCol,
           color = brewer.pal(8, "Blues"),
           ...)
}

##' @title          sampleFlag
##' @description    Flag samples in an experiment based on the number of non-detectable genes and sample-by-sample correlations - obsolete
##'
##' @details
##' Deprecated! - Samples are flagged if their median Pearson correlation coefficient with all other samples in the experiment is below 0.2 
##' or if the number of detected genes (as determined by lodAssess) is less than half the genes assayed.
##'
##' @param  rccSet      NanoString ExpressionSet object
##' @param  reference   scalar character, reference to determine limit of detection; either "blank" or "negatives"; passed on to lodAssess
##' @param  stringency  factor by which deviation (MAD) of the median of the reference counts is multiplied; passed on to lodAssess 
##'
##' @return             A numeric vector giving the indices of samples to be flagged
##'
##' @seealso \code{\link{lodAssess}} 
##'
##' @author Dorothee Nickles
##'
sampleFlag <- function(rccSet,
                       reference = c("blank", "negatives"),
                       stringency = 2.5) {
  stopifnot(is(rccSet, "ExpressionSet"))
  reference <- match.arg(reference)

  M <- assayData(rccSet)$exprs

  misValues <- lodAssess(rccSet, reference, stringency)
  toFlag <- unique(c(which(apply(cor(M), 2, median) < 0.2), 
                     which(misValues/sum(fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping")) > 0.5)))
  return(toFlag)
}

##' @title          geneClustering
##' @description    Clustering of genes by correlation across an experiment
##'
##' @param  M               One of the matrices from assayData(rccSet) where rccSet is a NanoString ExpressionSet object
##' @param  main            Plot title
##' @param  labels          Character vector for the axis labels
##' @param  mainImage       Filename for the main PNG file that will be generated
##' @param  previewImage    Filename for the preview PNG file that will be generated
##' @param  covar           Character; colname in fdata that can be used to label genes by a category of interest
##' @param  fdata           fData(rccSet)
##' @param  log.transform   Scalar boolean
##' @param  gsubset         Either a boolean vector of length nrow(exprs(rccSet)) or a numeric vector of indices, indicating
##'                         to which genes in exprs(rccSet) to subset
##' @param  ...             Extra arguments passed to aheatmap() (optional)
##'
##' @return
##' A NanoString ExpressionSet object that has count data adjusted by positive control counts.
##'
##' @importFrom NMF aheatmap
##' @importFrom png readPNG
##'
##' @author Dorothee Nickles, Robert Ziman
##'
geneClusteringPNG <- function(M,
                              main,
                              labels,
                              mainImage,
                              previewImage,
                              log.transform=FALSE,
                              covar=NULL,
                              fdata=NULL,
                              gsubset=NULL,
                              ...)
{
  stopifnot(is(M, "matrix"))

  if (nrow(M) > 1000)
    stop("This function is designed to plot a heatmap of up to 1000 genes")

  if (log.transform) {
    M <- log2(M)
  }

  #if (any(apply(exprs(rccTmp), 2, function(x) { length(unique(x))}) == 1)) {
  #    noVar <- which(apply(exprs(rccTmp), 2, function(x) { length(unique(x))}) == 1)
  #    tmp <- tmp[,-noVar]
  #    warning(sprintf("%i samples were excluded from correlation plot because of zero variance:\n %s", 
  #      length(noVar), paste(rownames(pData(rccTmp))[noVar], collapse=",")))
  #}

  rownames(M) <- labels

  if (!missing(covar) && !missing(fdata)) {
    if (!(covar %in% colnames(fdata)))
      stop("covar not found in fdata")
    plabel <- fdata[, covar]
    if (!missing(gsubset)) {
      M <- M[gsubset,]
      plabel <- plabel[gsubset]
    }
  } else {
    plabel <- NULL
  }

  if (nrow(M) < 50) {

    width = (20*nrow(M) + 250) / 72
    height = (20*nrow(M) + 250) / 72

    png( filename = mainImage, width = width, height = height, units="in", res=72 )

    aheatmap( cor(t(M)),
              annCol      = plabel,
              color       = "-RdBu",
              treeheight  = 0,
              legend      = TRUE,
              fontsize    = 10,
              cexRow      = 1,
              cexCol      = 1,
              cellwidth   = 20,
              cellheight  = 20,
             #gp=grid::gpar
              main        = main)

    dev.off()

  } else {

    width <- (10 * nrow(M) + 500) / 300
    height <- (10 * nrow(M) + 200) / 300

    png( filename = mainImage, width = width, height = height, units = "in", res = 300 )

    aheatmap( cor(t(M)),
              annCol      = plabel,
              color       = "-RdBu",
              treeheight  = 0,
              legend      = TRUE,
              fontsize    = 2,
             #cexRow      = 0.5,       # These won't work when the input matrix is large (>200 genes or so);
             #cexCol      = 0.5,       # use cellwidth/cellheight instead. See aheatmap() source.
              cellwidth   = 2,
              cellheight  = 2,
              gp          = grid::gpar(cex=2),  # Required to compensate for hardcoded values in the aheatmap() source
              main        = main )

    dev.off()

  }

  raster <- readPNG( mainImage )
  png( filename = previewImage, width = 500/72, nrow(raster)/ncol(raster)*500/72, units = "in", res = 72 )
  op <- par( mar = rep(0,4) )
  nPoints <- ncol(raster)
  plot( c(1:nrow(raster), rep(nrow(raster), nPoints - nrow(raster))),
    1:ncol(raster), type = "n", xlab = "", ylab = "", axes = FALSE )
  rasterImage( raster, 1, 1, nrow(raster), ncol(raster) )
  dev.off()
  par(op)

}

##' @title          densityPlot
##' @description    Plot the density of counts for all endogenous and housekeeping genes for each sample in a NanoString ExpressionSet object
##'
##' @param  M               One of the matrices from the assayData() of a NanoString ExpressionSet object (make sure
##'                         to set the log.transform parameter accordingly)
##' @param  log.transform   Scalar boolean
##' @param  pdata           pData() of the NanoString ExpressionSet object
##' @param  covar           character; colname in pData() that can be used to label genes by a category of interest (passed on to colByCovar)
##' @param  ...             additional plotting parameters
##'
##' @return
##' A density plot
##'
##' @author Dorothee Nickles
##'
densityPlot <- function(M, log.transform=FALSE, pdata, covar, ...) {	## takes exprSet
  
  stopifnot(is(M, "matrix"))

  if (log.transform) {
    M <- log2(M)  
  }
  
  cols <- colByCovar(pdata, covar)
  plot(density(M[,1]),
       ylim=c(0,0.5), col="white", ...)
  sapply(1:ncol(M),
         function(x) { lines(density(M[,x]),
                       col=cols[["color"]][x]) })
  legend("topright", fill=unique(cols[["color"]]), legend=cols[["legend"]])
}

##' @title          assessHousekeeping
##' @description    Assess correlation and variance/variability of housekeeping genes
##'
##' @details
##' Pairwise correlations of all defined housekeeping genes will be assessed and pairwise scatterplots will be generated. This function does not only
##' output pairwise correlation coefficients, but also - for each housekeeping gene - the variance, the interquartile range (IQR) and median expression
##' level across all samples in the experiment.
##'
##' @param  rccSet      NanoString ExpressionSet object
##' @param  hk          Either a boolean vector of length nrow(exprs(rccSet)) or a numeric vector of indices
##'                     which genes in exprs(rccSet) are housekeeping genes
##' @param  covar       character; colname in fData(rccSet) that can be used to label genes by a category of interest
##' @param  annotate    Scalar boolean; if TRUE (default), probes will be "annotated" using the "GeneName" column in the fData(rccSet) slot
##' @param  digits      Scalar integer, the number of decimal places
##' @param  plot        Scalar boolean, plot pairwise relationships ?
##'
##' @return
##' A dataframe with one row per housekeeping genes and several columns with metrics suggested to
##' assess performance of defined housekeeping genes.
##'
##' @author Dorothee Nickles
##'
assessHousekeeping <- function(rccSet, hk, covar, annotate=TRUE, plot=TRUE, digits=2) {
  stopifnot(is(rccSet, "ExpressionSet"))

  if (!is.logical(hk)) {
    stop("Please provide logical vector assigning whether each gene serves as a housekeeping gene or not.")
  }
  if (sum(hk) == 1) {
    stop("You need more than one housekeeping gene to use this function.")
  }

  M <- assayData(rccSet)$exprs

  values <- rownames(M[hk, ])
  if (annotate) {
    values <- fData(rccSet)$GeneName[hk]
  }
  if (any(duplicated(values))) {
    n <- 1
    for (r in which(duplicated(values))) {
      values[r] <- paste(values[r], n, sep="v")
      n <- n+1
    }
  }

  if(plot == TRUE){
    pairs(t(log2(M[hk, ])), 
          main = "Correlation of housekeeping genes",
          col = colByCovar(pData(rccSet), covar)[["color"]],
          labels = values,
          pch = 20)
  }

  output <- data.frame(matrix(nrow=sum(hk), ncol=0))
  output[,1] <- values
  output[,2] <- var(t(log2(M[hk, ])))[,1]
  output[,3] <- apply(log2(M[hk, ]), 1, IQR)
  output[,4] <- apply(log2(M[hk, ]), 1, median)
  output[, 5:(sum(hk)+4)] <- cor(t(log2(M[hk, ])))
  output[,c(2:ncol(output))] <- round(output[,c(2:ncol(output))], digits=digits)
  colnames(output) <- c("Gene", "Variance", "IQR", "median", paste("cor", values, sep="_"))

  return(output)
}

