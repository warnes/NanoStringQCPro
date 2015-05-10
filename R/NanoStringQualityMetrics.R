
##' @title          myCols
##' @description    Function that defines nice colors
##'
##' @return
##' A vector of colors
##'
##' @author Dorothee Nickles
myCols <- function() {
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
      ifelse(x > thresholds[[2]], "red", "black"))
}

##' @title          colByCovar
##' @description    Define colors based on a covariate of an RccSet object
##'
##' @param  pdata   pData() of an RccSet object
##' @param  covar   character, colname in the pdata used to stratify (color) data
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

setGeneric( "fovPlot", function( rccSet ) standardGeneric( "fovPlot" ) )
##' @rdname fovPlot
##' @aliases fovPlot
##'
##' @title Fields of view (FOV) plot
##'
##' @description
##' Plot the fraction of successfully imaged fields of view (FOV) in the given RccSet.
##' The RccSet's phenoData should have 'FovCount' and 'FovCounted' columns populated with the total
##' and successfully imaged FOV counts, respectively. Samples with a FOV counted/FOV count of less than
##' 80% are marked in red (dashed red line indicates threshold).
##'
##' @param  rccSet  An RccSet object
##'
##' @return
##' A plot
##'
##' @author Dorothee Nickles
##'
setMethod(
  "fovPlot",
  "RccSet",
  function(rccSet) {
    fov <- as.numeric(as.vector(pData(rccSet)$FovCounted))/
           as.numeric(as.vector(pData(rccSet)$FovCount))
    plot(fov,
         pch=16,
         main="FOV Counted versus FOV Count",
        #ylab="FOV Counted/FOV Count",
         ylab="ratio",
         xlab="Sample Index",
         las=1,
         col=ifelse(fov < 0.8, "red", "black"),
         ylim=c(min(fov, 0.8), 1))
    abline(h=0.8, col="red", lty=2)	
  }
  )

setGeneric( "bdPlot", function( rccSet ) standardGeneric( "bdPlot" ) )
##' @rdname bdPlot
##' @aliases bdPlot
##'
##' @title Binding density plot
##'
##' @description
##' Plot the binding density of each sample in an RccSet object.
##' Samples with a binding density < 0.05 or > 2.25 (thresholds defined by NanoString) 
##' are marked in red (dashed red line indicates threshold).
##'
##' @param  rccSet  An RccSet object
##'
##' @return
##' A plot
##'
##' @author Dorothee Nickles
##'
setMethod(
  "bdPlot",
  "RccSet",
  function(rccSet) {
    plot(as.numeric(as.vector(pData(rccSet)$BindingDensity)), pch=16, 
         ylim=c(0, max(c(as.numeric(as.vector(pData(rccSet)$BindingDensity)), 2.25))),
         main="Binding Density", ylab="Binding density",
         xlab="Sample Index",
         las=1,
         col=ifelse(as.numeric(as.vector(pData(rccSet)$BindingDensity)) > 2.25, "red", 
              ifelse(as.numeric(as.vector(pData(rccSet)$BindingDensity)) < 0.05, "red", "black")))
    abline(h=c(0.05, 2.25), col="red", lty=2)
  }
  )

setGeneric( "flagSamplesTech", function( rccSet ) standardGeneric( "flagSamplesTech" ) )
##' @rdname flagSamplesTech
##' @aliases flagSamplesTech
##'
##' @title          flagSamplesTech
##' @description    Flag samples based on their technical performance, i.e. field of vision (FOV) counted and binding density
##'
##' @details
##' Samples with a FOV counted/FOV count of less than 80% will be flagged. Also samples with a binding density
##'
##' @param  rccSet  An RccSet object
##'
##' @return
##' A numeric vector giving the indices of samples with outlier values in FOV counted and
##' binding density < 0.05 or > 2.25 (thresholds defined by NanoString) will be flagged.
##'
##' @author Dorothee Nickles
##'
setMethod(
  "flagSamplesTech",
  "RccSet",
  function(rccSet) {
    return(unique(c(which(as.numeric(as.vector(pData(rccSet)$FovCounted))/as.numeric(as.vector(pData(rccSet)$FovCount)) < 0.8),
                    which(as.numeric(as.vector(pData(rccSet)$BindingDensity)) > 2.25 | as.vector(pData(rccSet)$BindingDensity) < 0.05))))
  }
  )

setGeneric( "ctrlsOverviewPlot", function( rccSet ) standardGeneric( "ctrlsOverviewPlot" ) )
##' @rdname ctrlsOverviewPlot
##' @aliases ctrlsOverviewPlot
##'
##' @title          ctrlsOverviewPlot
##' @description    Plot individual negative and positive controls across all samples in an RccSet object.
##'
##' @param  rccSet  An RccSet object
##'
##' @return
##' A plot with two panels, one for the negative controls, one for the positive controls.
##'
##' @author Dorothee Nickles
##'
setMethod(
  "ctrlsOverviewPlot",
  "RccSet",
  function(rccSet) {
    par(mfrow=c(1,2))
  
    M <- assayData(rccSet)$exprs
  
    negCtrls <- (fData(rccSet)$CodeClass == "Negative")
    M.negCtrls <- M[ negCtrls, ]
    rownames(M.negCtrls) <- fData(rccSet)$GeneName[negCtrls]
    M.negCtrls.ordered <- M.negCtrls[ order(rownames(M.negCtrls)), ]
  
    boxplot(t(M.negCtrls.ordered),
            las=2,
            cex.axis=0.6,
            ylab="counts (raw)",
            log="y",
            main="Negative Controls")
  
    posCtrls <- which(fData(rccSet)$CodeClass == "Positive")
    M.posCtrls <- M[ posCtrls, ]
    rownames(M.posCtrls) <- fData(rccSet)$GeneName[posCtrls]
    M.posCtrls.ordered <- M.posCtrls[ order(rownames(M.posCtrls)), ]
  
    boxplot(t(M.posCtrls.ordered),
            las=2,
            cex.axis=0.6, 
            ylab="counts (raw)",
            log="y",
            main="Positive Controls")
  }
  )

setGeneric( "negCtrlsPlot", function( rccSet ) standardGeneric( "negCtrlsPlot" ) )
##' @rdname negCtrlsPlot
##' @aliases negCtrlsPlot
##'
##' @title          negCtrlsPlot
##' @description    Plot negative controls across all samples in an RccSet object
##'
##' @details
##' In the second panel, boxplots are colored by lane (as specified in the pData slot). Bars on top of
##' the panel indicate the stage position for each cartridge/sample (as specified in the pData slot).
##'
##' @param  rccSet  An RccSet object
##'
##' @return
##' A plot with two panels: one showing boxplots for the individual negative controls across all samples, and 
##' one showing boxplots for the negative control counts for each individual sample (lane-specific background).
##'
##' @author Dorothee Nickles
##'
setMethod(
  "negCtrlsPlot",
  "RccSet",
  function(rccSet) {
    
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
  )

setGeneric( "negCtrlsByLane", function( rccSet ) standardGeneric( "negCtrlsByLane" ) )
##' @rdname negCtrlsByLane
##' @aliases negCtrlsByLane
##'
##' @title          negCtrlsByLane
##' @description    Plot negative controls per lane in an RccSet object
##'
##' @details
##' Boxplots are colored by lane (as specified in the pData slot). Bars on top of
##' the panel indicate the stage position for each cartridge/sample (as specified in the pData slot).
##'
##' @param  rccSet  An RccSet object
##'
##' @return
##' A plot with boxplots for the negative control counts for each individual sample (lane-specific background)
##'
##' @author Dorothee Nickles
##'
setMethod(
  "negCtrlsByLane",
  "RccSet",
  function(rccSet) {
  
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
  )

##' Plot of negative controls by lane (vertical orientation)
##'
##' @param rccSet
##' An RccSet
##'
##' @param outputFile
##' Output PNG filename
##'
##' @return
##' A PNG file containing a boxplot of the counts for negative controls by lane
##' in the input. The PNG is set to a fixed resolution of 300 pixels per inch and
##' a fixed width of 2250 pixels (i.e. 7.5" at 300ppi), but the height varies
##' with the size of the input. The font size is also fixed so that the labels
##' will be legible even for large datasets.
##'
negCtrlsByLane_verticalPlot <- function(rccSet, outputFile)
{
  rccSet.sorted <- rccSet[ , rev(order(rownames(pData(rccSet))))]
  
  negCtrls <- (fData(rccSet.sorted)$CodeClass == "Negative")
  
  M <- assayData(rccSet.sorted)$exprs
  laneid <- as.numeric(pData(rccSet.sorted)$LaneID)

  png(filename=outputFile,
      width=2250,
      height=200 + 30*ncol(M),
      res=300,
      units="px",
      pointsize=7
  )
  
  par(mai=c(0.5,3,0.5,0.75))
  
  par(yaxs="i")
  
  boxplot(M[negCtrls, ],
          horizontal=TRUE,
          log="x",
          xlab="Counts (raw)",
          ylab="",
          xlim=c(0,(ncol(M)+1)),  # This gets assigned to ylim by boxplot() when horizontal=TRUE.
          las=1,
          col=myCols()[laneid],
          main="Negative controls by lane -\nlane-specific background"
  )

  ## disable addition of StagePosition indicator for now
  #n <- 1
  #for (i in 1:length(unique(pData(rccSet)$StagePosition))) {
  #  tmp <- which(pData(rccSet)$StagePosition == pData(rccSet)$StagePosition[i])
  #  lines(y = rep(max(M[ negCtrls, ]) + 5, length(tmp)),
  #        x = n:(n+length(tmp)-1),
  #        col = i,
  #        lwd = 4,
  #        ypd = TRUE)
  #  n <- n + length(tmp)
  #}

  legend(x=max(M[ negCtrls, ])*1.2,
         y=(length(laneid)+3),
        #cex=0.5,
         xpd=TRUE,
         bty="n",
         fill=myCols()[sort(unique(laneid))],
         legend=paste(rep("Lane", length(unique(laneid))), sort(unique(laneid))))
  
  dev.off()
  
  return(TRUE)
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

setGeneric( "negCtrlsPairs", function( rccSet, ... ) standardGeneric( "negCtrlsPairs" ) )
##' @rdname negCtrlsPairs
##' @aliases negCtrlsPairs
##'
##' @title          negCtrlsPairs
##' @description    Pairs plot of negative controls across all samples in an RccSet object
##'
##' @param  rccSet          An RccSet object
##' @param  log.transform   boolean, whether data needs to be log2 transformed
##'
##' @return
##' Pairs plot of the negative controls with a scatter plot in the lower panel and
##' correlatation coefficients printed in the upper panel.
##'
##' @author Dorothee Nickles
##'
setMethod(
  "negCtrlsPairs",
  "RccSet",
  function(rccSet, log.transform=FALSE) {
  
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
  )

setGeneric( "posNormFactPlot", function( rccSet ) standardGeneric( "posNormFactPlot" ) )
##' @rdname posNormFactPlot
##' @aliases posNormFactPlot
##'
##' @title          posNormFactPlot
##' @description    Plot positive control scaling factor for each sample in an RccSet object
##'
##' @param  rccSet  An RccSet object
##'
##' @return
##' A plot of the positive control scaling factor for each sample in an RccSet object.
##' Samples with a positive control scaling factor < 0.3 or > 3 (thresholds defined by NanoString) are
##' marked in red (dashed red line indicates threshold).
##'
##' @author Dorothee Nickles
##'
setMethod(
  "posNormFactPlot",
  "RccSet",
  function(rccSet)
  {
    if ("posCtrlData" %in% ls(assayData(rccSet))) {
        M <- assayData(rccSet)$exprs
        positives <- (fData(rccSet)$CodeClass == "Positive")

        posFactor <- pData(rccSet)$PosFactor
      
        plot(posFactor,
             las = 1,
             cex.axis = 0.6,
             ylab = "Positive control scaling factor",
             pch = 16,
             xlab = "Sample index",
             main = "Positive control scaling factors across samples", 
             col = ifelse(posFactor < 0.3, "red", ifelse(posFactor > 3, "red", "black")),
             ylim = c(min(c(posFactor, 0.3)), max(c(posFactor, 3))))
        abline(h = c(0.3, 3),
               col = "red",
               lty = 2)
    } else {
        print("This plot requires the RccSet to have positive control normalized data")
    }
  }
  )

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

setGeneric( "ctrlsZprimePlot", function( rccSet ) standardGeneric( "ctrlsZprimePlot" ) )
##' @rdname ctrlsZprimePlot
##' @aliases ctrlsZprimePlot
##'
##' @title          ctrlsZprimePlot
##' @description    Plot distribution of counts and the Z' Factors comparing the nagative controls and
##'                 the three highest input positive controls of an RccSet object.
##'
##' @param  rccSet  An RccSet object
##'
##' @return         A plot
##'
##' @author Dorothee Nickles
##'
setMethod(
  "ctrlsZprimePlot",
  "RccSet",
  function(rccSet) {
  
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
  )

setGeneric( "iqrPlot", function( rccSet, ... ) standardGeneric( "iqrPlot" ) )
##' @rdname iqrPlot
##' @aliases iqrPlot
##'
##' @title
##' iqrPlot
##'
##' @description
##' Plot the interquartile range (IQR) for a certain code class of probes in an
##' RccSet object.
##'
##' @details
##' IQR of the specified code class for each sample in the RccSet are plotted
##' and outliers (as determined by the function specified in the method
##' argument) are marked in red (thresholds for outlier definition are plotted
##' as red dashed lines).
##'
##' @param rccSet
##' An RccSet object
##'
##' @param codeClass
##' Character string specifying the code class (as annotated in the
##' fData(rccSet)$CodeClass column) for which the IQR shall be determined.
##'
##' @param method
##' Character string specifying the method for outlier detection: either
##' "cutoffByMMAD" or "cutoffByVar".
##'
##' @param stringency
##' Numeric value passed to the cutoff function specified by the method
##' argument (see the 'd' argument of cutoffByMMAD and cutoffByVar).
##'
##' @return
##' A plot
##'
##' @seealso \code{\link{cutoffByMMAD}}, \code{\link{cutoffByVar}}
##'
##' @author Dorothee Nickles
##'
setMethod(
  "iqrPlot",
  "RccSet",
  function(rccSet,
           codeClass = c("Negative", "Positive", "Endogenous", "Housekeeping"),
           method = c("cutoffByMMAD", "cutoffByVar"),
           stringency = 4)
  {
    codeClass <- match.arg(codeClass)
    fun <- match.fun(method)
  
    M <- assayData(rccSet)$exprs
  
    iqr <- apply(log2(M[ fData(rccSet)$CodeClass == codeClass, ]), 2, IQR)
    thresholds <- fun(iqr, stringency)
  
    plot(iqr,
         las = 1,
         cex.axis = 0.6,
         ylab = "counts (log2)",
         pch = 16,
         xlab = "Sample index",
         main = sprintf("IQR of code class %s across samples", codeClass), 
         col = colByFun(iqr, thresholds), 
         ylim = c(min(c(iqr, unlist(thresholds))), max(c(iqr, unlist(thresholds)))))
  
    abline(h = unlist(thresholds),
           col = "red",
           lty = 2)
  }
  )

setGeneric( "posR2Plot", function( rccSet ) standardGeneric( "posR2Plot" ) )
##' @rdname posR2Plot
##' @aliases posR2Plot
##'
##' @title          posR2Plot
##' @description    Plot the R squared of linear fit of counts versus input for positive controls in an RccSet object.
##'
##' @details
##' R squared for each sample in the RccSet are plotted and samples with R squared < 0.95 are marked in red
##' (threshold indicated by dashed red line).
##'
##' @param  rccSet  RccSet object
##'
##' @return         A plot
##'
##' @author Dorothee Nickles
##'
setMethod(
  "posR2Plot",
  "RccSet",
  function(rccSet)
  {
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
  )

setGeneric( "posSlopePlot", function( rccSet, ... ) standardGeneric( "posSlopePlot" ) )
##' @rdname posSlopePlot
##' @aliases posSlopePlot
##'
##' @title
##' posSlopePlot
##'
##' @description
##' Plot the slope of linear fit of counts versus input for positive controls
##' in an RccSet object.
##'
##' @param rccSet
##' An RccSet object
##'
##' @param method
##' Character string specifying the method for outlier detection: either
##' "cutoffByMMAD" or "cutoffByVar".
##'
##' @param stringency
##' Numeric value passed to the cutoff function specified by the method
##' argument (see the 'd' argument of cutoffByMMAD and cutoffByVar).
##'
##' @return
##' A plot
##'
##' @details
##' The slope for each sample in the RccSet are plotted and and outliers (as
##' determined by the function specified by the method argument) are marked in
##' red (thresholds for outlier definition are plotted as red dashed lines).
##'
##' @author Dorothee Nickles
##'
setMethod(
  "posSlopePlot",
  "RccSet",
  function(rccSet,
           method = c("cutoffByMMAD", "cutoffByVar"),
           stringency = 4)
  {
    fun <- match.fun(method)
  
    M <- assayData(rccSet)$exprs
    positives <- (fData(rccSet)$CodeClass == "Positive")
  
    tmp <- apply(log2(M[positives, ]),
                 2,
                 function(x) { coefficients(lm(x ~ log2(fData(rccSet)$SpikeInInput[ positives ])))[2] })
    thresholds <- fun(tmp, stringency)
    plot(tmp,
         pch = 16,
         col = colByFun(tmp, thresholds),
         main = "Slope in positive controls",
         ylab = "slope",
         xlab = "Sample index",
         ylim = c(min(c(tmp, unlist(thresholds))), max(c(tmp, unlist(thresholds)))),
         las = 1)
  
    abline(h=unlist(thresholds), col="red", lty=2)
  }
  )


setGeneric( "posRatioPlot", function( rccSet, ... ) standardGeneric( "posRatioPlot" ) )
##' @rdname posRatioPlot
##' @aliases posRatioPlot
##'
##' @title
##' posRatioPlot
##'
##' @description
##' Plot the ratio of the mean of positive control counts for each sample and
##' the overall mean of positive control counts in an RccSet object.
##'
##' @details
##' The ratio for each sample in the RccSet is plotted and and outliers (as
##' determined the cutoff function specified by the method argument) are marked
##' in red (thresholds for outlier definition are plotted as red dashed lines).
##'
##' @param rccSet
##' An RccSet object
##'
##' @param method
##' Character string specifying the method for outlier detection: either
##' "cutoffByMMAD" or "cutoffByVar".
##'
##' @param stringency
##' Numeric value passed to the cutoff function specified by the method
##' argument (see the 'd' argument of cutoffByMMAD and cutoffByVar).
##'
##' @return
##' A plot
##'
##' @author Dorothee Nickles
##'
setMethod(
  "posRatioPlot",
  "RccSet",
  function(rccSet,
           method = c("cutoffByMMAD", "cutoffByVar"),
           stringency = 4)
  {
    fun <- match.fun(method)
  
    M <- assayData(rccSet)$exprs
    positives <- (fData(rccSet)$CodeClass == "Positive")
  
    meanPos <- log2(colMeans(M[positives, ]))
    meanPosMean <- mean(meanPos)
    meanPosRatio <- meanPos/meanPosMean
    thresholds <- fun(meanPosRatio, stringency)
  
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
  )

setGeneric( "allSumPlot", function( rccSet, ... ) standardGeneric( "allSumPlot" ) )
##' @rdname allSumPlot
##' @aliases allSumPlot
##'
##' @title
##' allSumPlot
##'
##' @description
##' Plot the sum of all counts (endogenous and housekeeping genes only) for each
##' sample in an RccSet object.
##'
##' @details
##' The sum of counts for each sample in the RccSet are plotted and and outliers
##' (as determined the cutoff function specified by the method argument) are
##' marked in red (thresholds for outlier definition are plotted as red dashed
##' lines).
##'
##' @param rccSet
##' An RccSet object
##'
##' @param method
##' Character string specifying the method for outlier detection: either
##' "cutoffByMMAD" or "cutoffByVar".
##'
##' @param stringency
##' Numeric value passed to the cutoff function specified by the method
##' argument (see the 'd' argument of cutoffByMMAD and cutoffByVar).
##'
##' @return
##' A plot
##'
##' @author Dorothee Nickles
##'
setMethod(
  "allSumPlot",
  "RccSet",
  function(rccSet,
           method = c("cutoffByMMAD", "cutoffByVar"),
           stringency = 4)
  {
    fun <- match.fun(method)
  
    M <- assayData(rccSet)$exprs
    nonctrls <- (fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping"))
  
    allSum <- colSums(M[nonctrls, ])
    thresholds <- fun(log2(allSum), stringency)
    blankLabel <- getBlankLabel(rccSet)
  
    plot(y = allSum,
         x = 1:ncol(M),
         log = "y",
         las = 1,
         cex.axis = 0.6,
         ylab = "sum of counts",
         xlab = "Sample index",
         main = "Sum of raw counts for all probes", 
         col = colByFun(log2(allSum), thresholds),
         pch = ifelse( (pData(rccSet)$SampleType %in% blankLabel), 17, 16),
         ylim = as.integer(c(min(c(allSum, 2^unlist(thresholds))), max(c(allSum, 2^unlist(thresholds)))))
         )
  
    abline(h=2^unlist(thresholds), col="red", lty=2)
  }
  )

setGeneric( "posSumVsAllSumPlot", function( rccSet, ... ) standardGeneric( "posSumVsAllSumPlot" ) )
##' @rdname posSumVsAllSumPlot
##' @aliases posSumVsAllSumPlot
##'
##' @title
##' posSumVsAllSumPlot
##'
##' @description
##' Plot the ratio of sums of positive control counts to all counts for all
##' samples in an RccSet object.
##'
##' @details
##' The ratio for each sample in the RccSet is plotted and and outliers (as
##' determined by the cutoff function specified by the method argument) are
##' marked in red (thresholds for outlier definition are plotted as red
##' dashed lines).
##'
##' @param rccSet
##' An RccSet object
##'
##' @param method
##' Character string specifying the method for outlier detection: either
##' "cutoffByMMAD" or "cutoffByVar".
##'
##' @param stringency
##' Numeric value passed to the cutoff function specified by the method
##' argument (see the 'd' argument of cutoffByMMAD and cutoffByVar).
##' (If the median ratio is less than 1, three times this value will be
##' used.)
##'
##' @return
##' A plot
##'
##' @author Dorothee Nickles
##'
setMethod(
  "posSumVsAllSumPlot",
  "RccSet",
  function(rccSet,
           method = c("cutoffByMMAD", "cutoffByVar"),
           stringency = 4)
  {
    fun <- match.fun(method)
  
    M <- assayData(rccSet)$exprs
  
    iSum <- colSums(M[fData(rccSet)$CodeClass %in% "Positive",])
    nSum <- colSums(M[fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping"),])

    if (median(iSum/nSum) < 1) {
      stringency <- 3*stringency
    }

    thresholds <- fun(iSum/nSum, stringency)

    blankLabel <- getBlankLabel(rccSet)
  
    plot(y = iSum/nSum,
         x = 1:ncol(M),
         las = 1,
         cex.axis = 0.6,
         ylab = "ratio",
         xlab = "Sample index",
         main = sprintf("Sum of positive control counts /\nsum of endogenous and housekeeping gene counts"),
         col = colByFun(iSum/nSum, thresholds),
         pch = ifelse( (pData(rccSet)$SampleType %in% blankLabel), 17, 16),
         ylim = c(min(c(iSum/nSum, unlist(thresholds))), max(c(iSum/nSum, unlist(thresholds)))))

    abline(h=unlist(thresholds), col="red", lty=2)
  }
  )


setGeneric( "flagSamplesCtrl", function( rccSet, ... ) standardGeneric( "flagSamplesCtrl" ) )
##' @rdname flagSamplesCtrl
##' @aliases flagSamplesCtrl
##'
##' @title
##' flagSamplesCtrl
##'
##' @description
##' Flag samples based on the performance of their controls.
##'
##' @details
##' The method and stringency arguments determine a cutoff value used to
##' flag outlier samples based on the interquartile range of negative and
##' positive controls and in the slope of the linear fit of count versus input
##' of positive controls. Outliers will also be flagged if their positive
##' control scaling factor (see posCtrlNorm) is outside the NanoString
##' recommended range (i.e. below 0.3 or greater than 3).
##'
##' @param rccSet
##' An RccSet object
##'
##' @param method
##' Character string specifying the method for outlier detection: either
##' "cutoffByMMAD" or "cutoffByVar".
##'
##' @param stringency
##' Numeric value passed to the cutoff function specified by the method
##' argument (see the 'd' argument of cutoffByMMAD and cutoffByVar).
##'
##' @return
##' A numeric vector giving the indices of samples with outlier values as
##' described above.
##'
##' @author Dorothee Nickles
##'
setMethod(
  "flagSamplesCtrl",
  "RccSet",
  function(rccSet,
           method = c("cutoffByMMAD", "cutoffByVar"),
           stringency = 4)
  {
    method <- match.arg(method)
    fun <- match.fun(method)

    M <- assayData(rccSet)$exprs
  
    negatives    <- (fData(rccSet)$CodeClass == "Negative")
    positives    <- (fData(rccSet)$CodeClass == "Positive")
    positive_128 <- (fData(rccSet)$SpikeInInput == 128)
  
    a <- apply(log2(M[negatives, ]), 2, IQR)
    ctrls <- which(colByFun(a, fun(a, stringency)) == "red")
  
    a <- apply(log2(M[positives, ]), 2, IQR)
    ctrls <- c(ctrls, which(colByFun(a, fun(a, stringency)) == "red"))
  
    #a <- log2(M[positive_128, ])
    #ctrls <- c(ctrls, which(colByFun(a, fun(a, stringency)) == "red"))
  
    a <- apply(log2(M[positives, ]),
               2,
               function(x) { coefficients(lm(x ~ log2(fData(rccSet)$SpikeInInput[ positives ])))[2] })
    ctrls <- c(ctrls, which(colByFun(a, fun(a, stringency)) == "red"))

    if ("PosFactor" %in% rownames(pData(rccSet))) {
        posFactor <- pData(rccSet)$PosFactor
        ctrls2 <- c(which(posFactor < 0.3), which(posFactor > 3))
    } else {
        ctrls2 <- integer(0)
    }

    return(unique(ctrls, ctrls2))
  }
  )

##' @title          pcaPlot
##' @description    Wrapper function to perform a PCA analysis on the exprs slot of an RccSet object
##'                 and plot some results
##'
##' @param  exx     exprs() of an RccSet object
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

setGeneric( "lodAssess", function( rccSet, ... ) standardGeneric( "lodAssess" ) )
##' @rdname lodAssess
##' @aliases lodAssess
##'
##' @title
##' lodAssess
##'
##' @description
##' Assess how many genes in each sample in an RccSet object are below the
##' limit of detection. (The current implementation does a straightforward
##' column sum on the presence/absence matrix (paData) in assayData.)
##'
##' @param rccSet
##' An RccSet object
##'
##' @return
##' A numeric vector giving the number of missing genes (endogenous and
##' housekeeping genes) for each sample in an RccSet. If paData is not
##' found in the input's assayData, NULL is returned.
##'
##' @author Dorothee Nickles
##'
setMethod(
  "lodAssess",
  "RccSet",
  function(rccSet)
  {
    if ("paData" %in% ls(assayData(rccSet))) {
      nonctrls <- (fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping"))
      misValues <- colSums( assayData(rccSet)$paData[nonctrls, ] == FALSE )
      return(misValues)
    } else {
      return(NULL)
    }
  }
  )

setGeneric( "flagSamplesCount", function( rccSet, ... ) standardGeneric( "flagSamplesCount" ) )
##' @rdname flagSamplesCount
##' @aliases flagSamplesCount
##'
##' @title
##' flagSamplesCount
##'
##' @description
##' Flag samples based on overall counts
##'
##' @details
##' The method and stringency arguments deterine a cutoff value used to flag
##' samples as outliers: samples will be flagged if the sum of counts of their
##' endogeneous genes exceeds the cutoff or if the ratio of the sums of their
##' positive controls to the sums of their endogenous genes exceeds three
##' times the cutoff. Samples will also be flagged if the fraction of
##' genes below the lower limit of detection exceeds the maxMiss value.
##'
##' @param rccSet
##' An RccSet object
##'
##' @param method
##' Character string specifying the method for outlier detection: either
##' "cutoffByMMAD" or "cutoffByVar".
##'
##' @param stringency
##' Numeric value passed to the cutoff function specified by the method
##' argument (see the 'd' argument of cutoffByMMAD and cutoffByVar).
##'
##' @param maxMiss
##' Numeric specifying the allowable fraction of genes below the lower limit
##' of detection in a sample.
##'
##' @return
##' A numeric vector giving the indices of samples with outlier values
##' according to the criteria described above.
##'
##' @author Dorothee Nickles
##'
setMethod(
  "flagSamplesCount",
  "RccSet",
  function(rccSet,
           method = c("cutoffByMMAD", "cutoffByVar"),
           stringency = 4,
           maxMiss = 0.2)
  {
    fun <- match.fun(method)
  
    M <- assayData(rccSet)$exprs
  
    positives  <- (fData(rccSet)$CodeClass == "Positive")
    endogenous <- (fData(rccSet)$CodeClass == "Endogenous")
  
    posSum <- colSums(M[positives, ])
    allSum <- colSums(M[endogenous, ])

    if ("paData" %in% ls(assayData(rccSet))) {
      misValues <- lodAssess(rccSet)
      flagged_maxMiss <- which(misValues/nrow(M) > maxMiss)
    } else {
      flagged_maxMiss <- integer(0)
    }

    sampC <- unique(c(
      which(colByFun(allSum, fun(allSum, stringency)) == "red"),
      which(colByFun(posSum/allSum, fun(posSum/allSum, stringency*3)) == "red"),
      flagged_maxMiss
    ))
  
    return(sampC)
  }
  )

setGeneric( "lodPlot", function( rccSet, ... ) standardGeneric( "lodPlot" ) )
##' @rdname lodPlot
##' @aliases lodPlot
##'
##' @title
##' lodPlot
##'
##' @description
##' Function to plot the number of missing genes per sample in an RccSet object.
##'
##' @details
##' Samples with more than 50% missingness are marked in red. If blank
##' measurements are present, they are represented as triangles.
##'
##' @param rccSet
##' An RccSet object
##'
##' @param maxMiss
##' Numeric specifying the allowable fraction of genes below the lower limit
##' of detection in a sample.
##'
##' @return
##' A plot
##'
##' @author Dorothee Nickles
##'
setMethod(
  "lodPlot",
  "RccSet",
  function(rccSet,
           maxMiss = 0.2)
  {
    if ("paData" %in% ls(assayData(rccSet))) {

        blankLabel <- getBlankLabel(rccSet)
        nonctrls <- (fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping"))
        n_nonctrls <- sum(nonctrls)
      
        misValues <- lodAssess(rccSet)
      
        plot(misValues/n_nonctrls, 
             col = ifelse(misValues/n_nonctrls >= maxMiss, "red", "black"),
             pch = ifelse(pData(rccSet)$SampleType %in% blankLabel, 17, 16),
             main = "Sample Missingness",
             ylab = "Fraction of genes below background estimate",
             xlab = "Sample index",
             las = 1)
      
        abline(h=maxMiss, col="red", lty=2)   # was: abline(h=length(nonctrls)*maxMiss, col="red", lty=2) -RZ 2015-03-27

    } else {
        print("This plot requires the RccSet to have the presence/absence matrix")
    }
  }
  )

setGeneric( "sampleClustering", function( rccSet, ... ) standardGeneric( "sampleClustering" ) )
##' @rdname sampleClustering
##' @aliases sampleClustering
##'
##' @title
##' Clustering by sample correlation
##' 
##' @description
##' Clustering by sample correlation
##'
##' @param rccSet
##' An RccSet
##'
##' @param outputFile
##' Output PNG filename
##'
##' @param main
##' Plot title
##' 
##' @param annCol
##' See \code{\link{aheatmap}}
##'
##' @param covar
##' Covariate (e.g. "SampleType")
##'
##' @return
##' A PNG file showing clustering of samples by correlation. Any zero-variance
##' samples are omitted from the heatmap.
##' The width and height of the PNG file are set to vary with the size of the input.
##'
##' @importFrom NMF aheatmap
##' @importFrom RColorBrewer brewer.pal
##'
##' @author Dorothee Nickles, Robert Ziman
##'
setMethod(
  "sampleClustering",
  "RccSet",
  function(rccSet,
           outputFile,
           main = "Sample correlations in raw data",
           annCol = NULL,
           covar = "SampleType")
  {
    M <- assayData(rccSet)$exprs
  
    # Remove any zero-variance rows. These rows are problematic
    # in the heatmap-plotting code that follows (they cause NA values in the output of cor() which
    # in turn causes hclust() to error out).

    zero_variance <- apply(M, 2, function(x){ var(x) == 0 })
    if (any(zero_variance)) {
      warning(
        paste(
          c(
            "The following samples were excluded because of zero variance:",
            paste0("    ", colnames(M)[zero_variance])
          ),
          collapse="\n"
        )
      )
      M <- M[ , !zero_variance]
      pdata <- pData(rccSet)[ !zero_variance ] 
    } else {
     #M <- M
      pdata <- pData(rccSet)
    }

    # Assign annCol (specifications of column annotation tracks displayed as
    # coloured rows on top of the heatmaps)

    cols <- colByCovar(pdata, covar)
    if (missing(annCol)) {
      annCol <- data.frame(SampleSubgroup = pdata[ , covar])
    }

    # Check if the input is too large to be comfortably plotted.

    if (ncol(M) > 1000)
      stop("This function is designed to plot a heatmap of up to 1000 samples")

    # Generate the heatmap. Given the finicky aheatmap() behavior, need different
    # settings when the input is small (<50 samples) or large (>50). 

    if (ncol(M) < 50) {

      width = (20*ncol(M) + 250) / 72
      height = (20*ncol(M) + 250) / 72

      png( filename = outputFile, width = width, height = height, units="in", res=72 )

      aheatmap( cor(M),
                annCol      = annCol,
                color       = brewer.pal(8, "Blues"),
                treeheight  = 0,
                legend      = TRUE,
                fontsize    = 10,
                cexRow      = 1,
                cexCol      = 1,
                cellwidth   = 20,
                cellheight  = 20,
                main        = main)

      dev.off()

    } else {

      width <- (10 * ncol(M) + 600) / 300
      height <- (10 * ncol(M) + 200) / 300

      png( filename = outputFile, width = width, height = height, units = "in", res = 300 )

      aheatmap( cor(M),
                annCol      = annCol,
                color       = brewer.pal(8, "Blues"),
                treeheight  = 0,
                legend      = TRUE,
                fontsize    = 2,
               #cexRow      = 0.5,       # These won't work when the input matrix is large (>200 columns);
               #cexCol      = 0.5,       # use cellwidth/cellheight instead. See aheatmap() source.
                cellwidth   = 2,
                cellheight  = 2,
                gp          = grid::gpar(cex=2),  # Required to compensate for hardcoded values in the aheatmap() source
                main        = main )

      dev.off()

    }

  }
  )

##' Plot counts in blank samples (vertical orientation)
##'
##' @param rccSet
##' An RccSet
##'
##' @param outputFile
##' Output PNG filename
##'
##' @return
##' A PNG file containing a boxplot of the gene-wise counts for the blank samples
##' in the input. The PNG is set to a fixed resolution of 300 pixels per inch and
##' a fixed width of 2250 pixels (i.e. 7.5" at 300ppi), but the height varies
##' with the size of the input. The font size is also fixed so that the labels
##' will be legible even for large datasets.
##'
countsInBlankSamples_verticalPlot <- function(rccSet, outputFile)
{
  if (hasBlanks(rccSet))
  {
    blankLabel <- getBlankLabel(rccSet)
    blanks <- (pData(rccSet)$SampleType %in% blankLabel)
    endo   <- which(fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping"))
    
    M <- assayData(rccSet)$exprs
    rownames(M) <- fData(rccSet)$GeneName   
    
    png(filename=outputFile,
        width=2250,
        height=200 + 30*nrow(M),
        res=300,
        units="px",
        pointsize=7
        )
    
    par(mai=c(0.5,1.5,0.5,0.25))

    par(yaxs="i")
    
    boxplot(t(M[endo, blanks]),
            horizontal=TRUE,
            log="x",
            xlab="Counts (raw)",
            ylab="",
            las=1,
            pch=".",
            main="Gene-wise counts in blank measurements"
            )
    sapply(endo,
           function(n) {
             points(x = median(M[n, blanks]) + (2.5*mad(M[n, blanks])),
                      #
                      # ***** TODO: replace with actual background estimate *****
                      #
                    y=n,
                    col="red",
                    pch=20,
                    cex=0.7)
           })
    
    dev.off()
  }
  else {
    print("No blank samples are present")
  }
  
  return(TRUE)
}


##' @title          densityPlot
##' @description    Plot the density of counts for all endogenous and housekeeping genes for each sample in an RccSet object
##'
##' @param  M               One of the matrices from the assayData() of an RccSet object (make sure
##'                         to set the log.transform parameter accordingly)
##' @param  log.transform   Scalar boolean
##' @param  pdata           pData() of the RccSet object
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
       ylim=c(0,0.5),
       col="white",
       las=1,
       ...)
  sapply(1:ncol(M),
         function(x) { lines(density(M[,x]),
                       col=cols[["color"]][x]) })
  legend("topright", fill=unique(cols[["color"]]), legend=cols[["legend"]])
}

setGeneric( "assessHousekeeping", function( rccSet, ... ) standardGeneric( "assessHousekeeping" ) )
##' @rdname assessHousekeeping
##' @aliases assessHousekeeping
##'
##' @title          assessHousekeeping
##' @description    Assess correlation and variance/variability of housekeeping genes
##'
##' @details
##' Pairwise correlations of all defined housekeeping genes will be assessed and pairwise scatterplots will be generated. This function does not only
##' output pairwise correlation coefficients, but also - for each housekeeping gene - the variance, the interquartile range (IQR) and median expression
##' level across all samples in the experiment.
##'
##' @param  rccSet      An RccSet object
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
setMethod(
  "assessHousekeeping",
  "RccSet",
  function(rccSet, hk, covar, annotate=TRUE, plot=TRUE, digits=2) {
  
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
            pch = 20,
            las = 1)
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
  )


##' Gene clustering heatmap
##'
##' @param rccSet
##' An RccSet
##'
##' @param outputFile
##' Output PNG filename
##'
##' @param main
##' Plot title
##'
##' @param covar
##' Colname in the rccSet's fData that can be used to label genes by a category of interest
##'
##' @return
##' A PNG file showing clustering of genes by correlation across an experiment
##' Positive and negative control probes and any zero-variance genes (typically
##' housekeeping genes) are omitted from the heatmap.
##' The width and height of the PNG file are set to vary with the size of the input.
##'
##' @importFrom NMF aheatmap
##'
##' @author Dorothee Nickles, Robert Ziman
##'
geneClustering <- function(rccSet,
                           outputFile,
                           main = "Gene clustering",
                           covar = NULL
                          #gsubset = NULL
                               ## @param gsubset
                               ## Either a boolean vector of length nrow(exprs(rccSet)) or a numeric vector of
                               ## indices indicating a subset of genes in exprs(rccSet) with which to generate
                               ## the heatmap instead of the full set.
                           )
{
  stopifnot(is(rccSet, "RccSet"))

  M <- assayData(rccSet)$normData
  rownames(M) <- fData(rccSet)$GeneName

  # Remove positive and negative control probes from the matrix.

  nonctrls <- (fData(rccSet)$CodeClass %in% c("Endogenous", "Housekeeping"))
  M <- M[nonctrls, ]

  # Remove any zero-variance rows -- typically the housekeeping genes. These rows are problematic
  # in the heatmap-plotting code that follows (they cause NA values in the output of cor() which
  # in turn causes hclust() to error out), and the housekeeping genes are not informative anyway.

  zero_variance <- apply(M, 1, function(x) { var(x) == 0 })
  if (any(zero_variance)) {
    warning(
      paste(
        c(
          "The following features were excluded because of zero variance:",
          paste0("    ", rownames(M)[zero_variance])
        ),
        collapse="\n"
      )
    )
    M <- M[ !zero_variance, ]
  }

  # Check if the input is too large to be comfortably plotted.

  if (nrow(M) > 1000)
    stop("This function is designed to plot a heatmap of up to 1000 genes")

  # Generate the heatmap. Given the finicky aheatmap() behavior, need different
  # settings when the input is small (<50 genes) or large (>50).

  if (nrow(M) < 50) {

    width = (20*nrow(M) + 250) / 72
    height = (20*nrow(M) + 250) / 72

    png( filename = outputFile, width = width, height = height, units="in", res=72 )

    aheatmap( cor(t(M)),
             #annCol      = plabel,
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

    png( filename = outputFile, width = width, height = height, units = "in", res = 300 )

    aheatmap( cor(t(M)),
             #annCol      = plabel,
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
}


centeredSampleClustering <- function(gmrccSet,
                                     outputFile,
                                     main,
                                     flags,
                                     excludeFlagged)
{
    if (missing(gmrccSet) || missing(outputFile) || missing(main) || missing(flags) || missing(excludeFlagged))
        stop("missing args")

    if (excludeFlagged)
    {
        flagged <- unique(unlist(apply(flags, 2, which)))

        temp_rccSet <- copyRccSet(gmrccSet)[, -flagged]
        
        gmrccSet <- contentNorm(temp_rccSet,
                                method = "global",
                                summaryFunction = "median",
                                inputMatrix = preproc(gmrccSet)$normData_inputMatrix,
                                quietly = TRUE)
    }

    nonctrls <- (fData(gmrccSet)$CodeClass %in% c("Endogenous", "Housekeeping"))

    M <- t(
            apply( assayData(gmrccSet)$normData[nonctrls, ],
                  1,
                  scale, center = TRUE, scale = FALSE )
          )
    colnames(M) <- colnames( assayData(gmrccSet)$normData[nonctrls, ] )
        #
        # TODO: The apply(... 1, scale, ...) combination above transposes the matrix
        # and strips the rownames. See if the same output -- un-transposed and with names
        # intact -- can be generated with a direct call to scale(). -RZ 2015-04-24
        #
    rownames(M) <- fData(gmrccSet)$GeneName[nonctrls]

    # Remove any zero-variance rows or cols. These are problematic
    # in the heatmap-plotting code that follows (they cause NA values in the output of cor() which
    # in turn causes hclust() to error out).

    zero_variance <- apply(M, 1, function(x) { var(x) == 0 })
    if (any(zero_variance)) {
      warning(
        paste(
          c(
            "The following genes were excluded because of zero variance:",
            paste0("    ", rownames(M)[zero_variance])
          ),
          collapse="\n"
        )
      )
      M <- M[!zero_variance, ]
    }

    zero_variance <- apply(M, 2, function(x) { var(x) == 0 })
    if (any(zero_variance)) {
      warning(
        paste(
          c(
            "The following samples were excluded because of zero variance:",
            paste0("    ", colnames(M)[zero_variance])
          ),
          collapse="\n"
        )
      )
      M <- M[ , !zero_variance]
    }

    # Check if the input is too large to be comfortably plotted.

    if (nrow(M) > 1000 || ncol(M) > 1000)
        stop("This function is designed to plot a heatmap of up to 1000 genes and 1000 samples")

    # Generate the heatmap. Given the finicky aheatmap() behavior, need different
    # settings when the input is small (<50 samples) or large (>50).

    row.cl <- as.dendrogram( hclust( as.dist( 1-cor(t(M)) ), method="ward.D" ) )
    col.cl <- as.dendrogram( hclust( as.dist( 1-cor(M) ), method="ward.D" ) )

    #annColors=lapply( seq_along( colnames( flags )), function(n){
    #  setNames( c("grey90", brewer.pal(ncol(flags), "Set1")[n]), c(FALSE, TRUE))
    #  })
    #names(annColors) <- colnames(flags)

    if (excludeFlagged)
        annCol <- data.frame(Counts = colSums(exprs(gmrccSet)[nonctrls, ]))
    else
        annCol <- flags

    if (nrow(M) < 50 && ncol(M) < 50) {

      width = (20*ncol(M) + 350) / 72
      height = (20*nrow(M) + 350) / 72

      png( filename = outputFile, width = width, height = height, units="in", res=72 )

      aheatmap( M, 
               #annColors   = annColors,
                annCol      = annCol,
                color       = "-RdBu",
                Colv        = col.cl,
                Rowv        = row.cl,
#                labRow      = NA,
#                labCol      = NA,
                scale       = "none",
                breaks      = 0,
                treeheight  = 0,
                legend      = TRUE,
                fontsize    = 10,
                cexRow      = 1,
                cexCol      = 1,
                cellwidth   = 20,
                cellheight  = 20,
                main        = main )

      dev.off()

    } else {

      width <- (10 * ncol(M) + 500) / 300
      height <- (10 * nrow(M) + 200) / 300

      png( filename = outputFile, width = width, height = height, units = "in", res = 300 )

      aheatmap( M, 
               #annColors   = annColors,
                annCol      = annCol,
                color       = "-RdBu",
                Colv        = col.cl,
                Rowv        = row.cl,
#                labRow      = NA,
#                labCol      = NA,
                scale       = "none",
                breaks      = 0,
                treeheight  = 0,
                legend      = TRUE,
                fontsize    = 2,
               #cexRow      = 0.5,       # These won't work when the input matrix is large (>200 columns);
               #cexCol      = 0.5,       # use cellwidth/cellheight instead. See aheatmap() source.
                cellwidth   = 2,
                cellheight  = 2,
                gp          = grid::gpar(cex=2),  # Required to compensate for hardcoded values in the aheatmap() source
                main        = main )

      dev.off()
      
    }
}

