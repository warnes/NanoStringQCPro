setGeneric( "nSolverBackground", function( rccSet, ... ) standardGeneric( "nSolverBackground" ) )

##' @rdname nSolverBackground
##' @aliases nSolverBackground
##' 
##' @title
##'
##' nSolver Analysis Software background estimation
##' 
##' @description
##'
##' Calculates initial probe- and lane-specific background estimates using an
##' algorithm that mimics the implementation in NanoString's nSolver Analysis
##' Software (see details below for the exact algorithm).
##'
##' @param rccSet
##'
##' NanoString RccSet object
##' 
##' @param stringency
##'
##' Multiplier with which to adjust final values.
##' 
##' @param shrink
##'
##' Boolean specifying if probe-specific estimates should be shrunken towards
##' their global mean.
##' 
##' @param w1
##'
##' Shrink weight "w1".
##'
##' @return
##'
##' A matrix containing lane- and probe-specific background estimates.
##'
##' @details
##'
##' The mean values for each blank lane (not including positive control probes)
##' are computed from the original data, and a vector of probe-specific
##' background is established by taking the rowMeans of the blank measurements
##' for each probe after subtracting out these values. If shrink=TRUE, the
##' vector is adjusted via the following formula (where 'probe.bg' represents
##' the vector):
##'
##' \preformatted{w2 <- 1/length(blanks)
##' probe.bg <- (w1*probe.bg + w2*mean(probe.bg)) / (w1 + w2)
##' }
##'
##' This probe-specific background is further adjusted by subtracting the mean
##' of its values for the negative control probes. A lane-specific "affinity" is
##' calculated for all lanes in the original data by taking the colMeans of the
##' negative control probe values in the original data, and background estimates
##' for each probe and lane in the original data are computed by summing the
##' corresponding probe-specific background and lane-specific affinity. Any
##' resulting values less than zero are set to zero, and the last step before
##' returning these values is to multiply them by the given stringency.
##'
##' @seealso \code{\link{getBackground}}, \code{\link{subtractBackground}}
##'
##' @author Dorothee Nickles, Thomas Sandmann
##'
setMethod(
    "nSolverBackground",
    "RccSet",
    function(rccSet,
             stringency=1,
             shrink=TRUE,
             w1=2.18)
    {
        stopifnot(hasBlanks(rccSet))

        blankLabel <- getBlankLabel(rccSet)

        state <- preproc(rccSet)$state
        if (state == "newRccSet") {
            M <- assayData(rccSet)$exprs
        } else if (state == "posCtrlNorm") {
            M <- assayData(rccSet)$posCtrlData
        } else {
            #warning("Input already contains background subtracted data")
            M <- assayData(rccSet)$posCtrlData
        }

	neg.ctrl.probes <- featureNames(rccSet)[ fData(rccSet)$CodeClass == "Negative" ]
  	not.positive.ctrl.probes <- featureNames(rccSet)[ fData(rccSet)$CodeClass != "Positive" ]
  	blanks <- (pData(rccSet)$SampleType %in% blankLabel)

	blank.bg <- colMeans( M[not.positive.ctrl.probes, blanks, drop=FALSE], na.rm=TRUE )
	probe.bg.diff <- sweep( M[, blanks, drop=FALSE], 2, blank.bg )
	probe.bg <- rowMeans(probe.bg.diff, na.rm=TRUE)
	
	if (shrink) {
            w2 <- 1/length(blanks)
            probe.bg <- (w1 * probe.bg + w2 * mean(probe.bg))/(w1 + w2) 
	}

	mean.neg.probe.bg <- mean( probe.bg[neg.ctrl.probes] )
  	probe.bg.model <- probe.bg - mean.neg.probe.bg

	laneAffinity <- colMeans( M[neg.ctrl.probes, ] )
	bg <- sapply(laneAffinity, function(x) probe.bg.model + x)
	bg[ bg < 0 ] <- 0
	bgEstimates <- bg * stringency

	return(bgEstimates)
    }
    )




setGeneric( "getBackground", function( rccSet, ... ) standardGeneric( "getBackground" ) )

##' @rdname getBackground
##' @aliases getBackground
##' 
##' @title
##' 
##' Get background estimates for a NanoString RccSet
##'
##' @description
##'
##' Returns background estimates for a NanoString RccSet object. The function
##' depends upon correct annotation in the RccSet: if the \code{reference}
##' argument is set to \code{"blank"}, it expects blank measurements (i.e.,
##' water runs) to have their phenoData SampleType set to the value indicating
##' blanks (see getBlankLabel(); normally this value would have been set using
##' an argument to newRccSet()). If \code{reference} is set to
##' \code{"Negative"}, then it expects to find the negative control probes via
##' CodeClass == "Negative". If set to "both", it expects both of the above and
##' will calculate initial background estimates using an algorithm that mimics
##' the implementation in NanoString's nSolver Analysis Software; see the
##' nSolverBackground() man page for details on the algorithm.
##'
##' @param  rccSet
##'
##' NanoString RccSet object.
##' 
##' @param  reference
##'
##' Measurements to use for background estimates: either "blank" (for blank
##' samples), "Negative" (for negative control probes), or "both". Blanks are
##' assumed to be indicated as in the description above.
##' 
##' @param  metric
##'
##' Summarization metric for background measurements: either "mean" or "median".
##' 
##' @param  stringency
##'
##' Factor by which deviation (SD or MAD) of the summarization metric will be
##' multiplied to obtain final background estimates.
##' 
##' @param  nSolverBackground.shrink
##'
##' Value to use for the 'shrink' argument to nSolverBackground().
##' 
##' @param  nSolverBackground.w1
##'
##' Value to use for the 'w1' argument to nSolverBackground().
##'
##' @return
##'
##' A matrix containing background estimates for a NanoString RccSet object.
##'
##' @export
##'
##' @examples
##'
##' data(example_rccSet)
##'
##' ## Calculate probe-specific background based on negative control probes
##' bg <- getBackground(example_rccSet, reference="Negative", metric="mean")
##'
##' ## Calculate sample-specific background based on blanks
##' bg <- getBackground(example_rccSet, reference="blank")
##'
##' ## Calculate background that is both sample- and probe-specific
##' bg <- getBackground(example_rccSet, reference="both", stringency=1)
##'
##' @seealso \code{\link{subtractBackground}}
##'
##' @author Dorothee Nickles
##'
setMethod(
    "getBackground",
    "RccSet",
    function(rccSet,
             reference=c("blank", "Negative", "both"),
             metric=c("median", "mean"),
             stringency=0,
             nSolverBackground.shrink=TRUE,
             nSolverBackground.w1=2.18)
    {
        reference <- match.arg(reference)
        metric <- match.arg(metric)
        stopifnot(stringency >= 0)

        state <- preproc(rccSet)$state
        if (state == "newRccSet") {
            M <- assayData(rccSet)$exprs
        } else if (state == "posCtrlNorm") {
            M <- assayData(rccSet)$posCtrlData
        } else {
            #warning("Input already contains background subtracted data")
            M <- assayData(rccSet)$posCtrlData
        }

        metric1 <- match.fun(metric)
        if (metric == "median") {
            metric2 <- mad
        } else {
            metric2 <- sd
        }

        blankLabel <- getBlankLabel(rccSet)

        if(reference == "both") {

            bgEstimates <- nSolverBackground(rccSet, stringency, shrink=nSolverBackground.shrink, w1=nSolverBackground.w1)

        } else if (reference == "blank") {

            stopifnot(hasBlanks(rccSet))

            blanks <- (pData(rccSet)$SampleType %in% blankLabel)
            refMed <- apply( M[, blanks], 1, metric1 )
            refMAD <- apply( M[, blanks], 1, metric2 )
            row_estimates <- refMed + refMAD * stringency
            bgEstimates <- matrix(nrow(M), ncol(M), byrow=FALSE, data=rep(row_estimates, nrow(pData(rccSet))))

        } else {                        # reference == "Negative"

            refMed <- apply( M[fData(rccSet)$CodeClass == "Negative", ], 2, metric1 )
            refMAD <- apply( M[fData(rccSet)$CodeClass == "Negative", ], 2, metric2 )
            col_estimates <- refMed + refMAD * stringency
            bgEstimates <- matrix(nrow(M), ncol(M), byrow=TRUE, data=rep(col_estimates, nrow(fData(rccSet))))

        }

        rownames(bgEstimates) <- rownames(fData(rccSet))
        colnames(bgEstimates) <- rownames(pData(rccSet))

        return(bgEstimates)
    }
    )




setGeneric( "subtractBackground", function( rccSet, ... ) standardGeneric( "subtractBackground" ) )

##' @rdname subtractBackground
##' @aliases subtractBackground
##'
##' @title
##'
##' Subtract background estimates for a NanoString RccSet
##'
##' @description
##'
##' Returns a NanoString \code{\linkS4class{RccSet}} with background-corrected
##' count data. During subtraction, any counts below 1 will be truncated to 1 to
##' enable subsequent log transformation of the data.
##'
##' @param  rccSet       NanoString RccSet object
##' @param  bgEstimates  Matrix containing the background estimates.
##' @param  description  Description for the background estimates.
##'
##' @return
##' 
##' A NanoString \code{linkS4class{RccSet}} object with background estimates
##' subtracted from the count data.
##'
##' @export
##'
##' @examples
##'
##' data(example_rccSet)
##'
##' bg1 <- getBackground(example_rccSet, reference="Negative", metric="mean")
##' bg2 <- getBackground(example_rccSet, reference="blank")
##' bg3 <- getBackground(example_rccSet, reference="both", stringency=1)
##'
##' bgCor1 <- subtractBackground(
##'     example_rccSet,
##'     bgEstimates=bg1,
##'     description="Background estimates using only negative control probes")
##'
##' bgCor2 <- subtractBackground(
##'     example_rccSet,
##'     bgEstimates=bg2,
##'     description="Background estimates using only blank samples")
##'
##' bgCor3 <- subtractBackground(
##'     example_rccSet,
##'     bgEstimates=bg3,
##'     description="Background estimates using both blank samples and negative control probes")
##'
##' @seealso \code{\link{getBackground}}
##'
##' @author Dorothee Nickles
##'
setMethod(
    "subtractBackground",
    "RccSet",
    function( rccSet, bgEstimates, description="" )
    {
        state <- preproc(rccSet)$state
        if (state == "newRccSet") {

            warning("Input does not contain positive control normalized data")

            M <- assayData(rccSet)$exprs
            bgCorrData_preprocList_value <- "Background corrected data"
            bgEstimates_preprocList_value <- description

        } else if (state == "posCtrlNorm") {

            M <- assayData(rccSet)$posCtrlData
            bgCorrData_preprocList_value <- "Positive control normalized and background corrected data"
            bgEstimates_preprocList_value <- description

        } else {

            #warning("Input already contains background subtracted data")

            M <- assayData(rccSet)$posCtrlData
            bgCorrData_preprocList_value <- "Positive control normalized and background corrected data"
            bgEstimates_preprocList_value <- description
        }

        M <- M - bgEstimates

        # Truncate to 1 any counts below 1 (cf. pseudocount in pdata_fdata_adata.to.rccSet())
        M[ M < 1 ] <- 1

        srccSet <- copyRccSet(rccSet) # Using copyRccSet() to be *sure* the original doesn't get affected in later code!

        # Update assayData and the preproc list
        assayData(srccSet)$bgCorrData <- M
        assayData(srccSet)$bgEstimates <- bgEstimates
        preproc(srccSet) <- c(preproc(srccSet),
                              assayData_bgCorrData=bgCorrData_preprocList_value,
                              assayData_bgEstimates=description)
        preproc(srccSet)$state <- "subtractBackground"

        return(srccSet)
    }
    )

