
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
##' @param inputMatrix
##'
##' Name of the matrix in the RccSet's assayData to use as input for calculating
##' background estimates (one of "exprs" or "posCtrlData"). If posCtrlData is
##' specified but not present in the assayData, an error will be generated.
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
             stringency = 1,
             shrink = TRUE,
             w1 = 2.18,
             inputMatrix = c("posCtrlData", "exprs"))
    {
        stopifnot(hasBlanks(rccSet))
        blankLabel <- getBlankLabel(rccSet)

        inputMatrix <- match.arg(inputMatrix)

        M <- assayData(rccSet)[[ inputMatrix ]]
        if (is.null(M))
            stop(sprintf("Specified input matrix ('%s') is not present in the RccSet's assayData", inputMatrix))

        neg.ctrl.probes <- featureNames(rccSet)[ fData(rccSet)$CodeClass == "Negative" ]
        not.positive.ctrl.probes <- featureNames(rccSet)[ fData(rccSet)$CodeClass != "Positive" ]
        blanks <- (pData(rccSet)$SampleType %in% blankLabel)

        blank.bg <- colMeans( M[not.positive.ctrl.probes, blanks, drop=FALSE], na.rm=TRUE )
        probe.bg.diff <- sweep( M[, blanks, drop=FALSE], 2, blank.bg )
        probe.bg <- rowMeans(probe.bg.diff, na.rm=TRUE)
        
        if (shrink) {
            w2 <- 1/sum(blanks)
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
##' depends upon correct annotation in the RccSet: if the \code{bgReference}
##' argument is set to \code{"blanks"}, it expects blank measurements (i.e.,
##' water runs) to have their phenoData SampleType set to the value indicating
##' blanks (see getBlankLabel(); normally this value would have been set using
##' an argument to newRccSet()). If \code{bgReference} is set to
##' \code{"negatives"}, then it expects to find the negative control probes via
##' CodeClass == "Negative". If set to \code{"both"}, it expects both of the above and
##' will calculate initial background estimates using an algorithm that mimics
##' the implementation in NanoString's nSolver Analysis Software (see the
##' nSolverBackground() man page for details on the algorithm).
##'
##' @param rccSet
##'
##' NanoString RccSet object.
##'
##' @param bgReference
##'
##' Measurements to use for background estimates: one of "blanks" (for blank
##' samples), "negatives" (for negative control probes), or "both". Blanks are
##' assumed to be indicated as in the description above.
##' 
##' @param summaryFunction
##'
##' Summary function for background measurements (e.g. "mean" or "median").
##' User-defined functions similar to these can be specified here as well.
##' 
##' @param stringency
##'
##' Factor by which deviation (SD or MAD) of the summarization output will be
##' multiplied to obtain final background estimates.
##' 
##' @param nSolverBackground.shrink
##'
##' Value to use for the 'shrink' argument to nSolverBackground().
##' 
##' @param nSolverBackground.w1
##'
##' Value to use for the 'w1' argument to nSolverBackground().
##'
##' @param inputMatrix
##'
##' Name of the matrix in the RccSet's assayData to use as input for calculating
##' background estimates (one of "exprs" or "posCtrlData"). If posCtrlData is
##' specified but not present in the assayData, an error will be generated.
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
##' bg <- getBackground(example_rccSet, bgReference="negatives", summaryFunction="mean",
##'     inputMatrix="exprs")
##'
##' ## Calculate sample-specific background based on blanks
##' bg <- getBackground(example_rccSet, bgReference="blanks", inputMatrix="exprs")
##'
##' ## Calculate background that is both sample- and probe-specific
##' bg <- getBackground(example_rccSet, bgReference="both", stringency=1,
##'     inputMatrix="exprs")
##'
##' @seealso \code{\link{subtractBackground}}
##'
##' @author Dorothee Nickles
##'
setMethod(
    "getBackground",
    "RccSet",
    function(rccSet,
             bgReference                = c("both", "blanks", "negatives"),
             summaryFunction            = "median",
             stringency                 = 0,
             nSolverBackground.shrink   = TRUE,
             nSolverBackground.w1       = 2.18,
             inputMatrix                = c("posCtrlData", "exprs"))
    {
        bgReference <- match.arg(bgReference)
        if (class(summaryFunction) == "character")
            sfun <- get(summaryFunction)
        else
            sfun <- summaryFunction
        stopifnot(stringency >= 0)
        inputMatrix <- match.arg(inputMatrix)

        M <- assayData(rccSet)[[ inputMatrix ]]
        if (is.null(M))
            stop(sprintf("Specified input matrix ('%s') is not present in the RccSet's assayData", inputMatrix))

        metric1 <- sfun
        if (as.character(quote(sfun)) == "median") {
            metric2 <- mad
        } else {
            metric2 <- sd
        }

        blankLabel <- getBlankLabel(rccSet)

        if (bgReference == "both") {

            bgEstimates <- nSolverBackground(rccSet,
                                             stringency,
                                             shrink = nSolverBackground.shrink,
                                             w1 = nSolverBackground.w1,
                                             inputMatrix)

        } else if (bgReference == "blanks") {

            stopifnot(hasBlanks(rccSet))

            blanks <- (pData(rccSet)$SampleType %in% blankLabel)
            refMed <- apply( M[, blanks], 1, metric1 )
            refMAD <- apply( M[, blanks], 1, metric2 )
            row_estimates <- refMed + refMAD * stringency
            bgEstimates <- matrix(nrow(M), ncol(M), byrow=FALSE, data=rep(row_estimates, nrow(pData(rccSet))))

        } else { # bgReference == "negatives"

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
##' @param rccSet
##'
##' NanoString RccSet object
##'
##' @param inputMatrix
##'
##' Name of the matrix in the RccSet's assayData to use as input for subtracting
##' background estimates (one of "exprs" or "posCtrlData"). If posCtrlData is
##' specified but not found in the assayData, an error will be generated.
##'
##' @param bgEstimates
##'
##' Matrix containing the background estimates to subtract.
##'
##' @param bgEstimatesParams
##'
##' A list with the parameters that were used to generate the background
##' estimates (see getBackground()):
##'
##' \itemize{
##'   \item bgReference
##'   \item summaryFunction
##'   \item stringency
##'   \item nSolverBackground.w1
##'   \item nSolverBackground.shrink
##'   \item inputMatrix
##' }
##'
##' The values of these list elements will be assigned to corresponding elements
##' in the output's experimentData@@preprocessing list. If any element is NULL,
##' the corresponding element in the output's preprocessing list will be NA.
##'
##' @param quietly
##'
##' Boolean specifying whether or not messages and warnings should be omitted.
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
##' pcnorm_rccSet <- posCtrlNorm(example_rccSet)
##'
##' bg1 <- getBackground(pcnorm_rccSet, bgReference="negatives", summaryFunction="mean")
##' bg2 <- getBackground(pcnorm_rccSet, bgReference="blanks")
##' bg3 <- getBackground(pcnorm_rccSet, bgReference="both", stringency=1)
##'
##' bgCor1 <- subtractBackground(pcnorm_rccSet, bgEstimates=bg1)
##' bgCor2 <- subtractBackground(pcnorm_rccSet, bgEstimates=bg2)
##' bgCor3 <- subtractBackground(pcnorm_rccSet, bgEstimates=bg3)
##'
##' @seealso \code{\link{getBackground}}
##'
##' @author Dorothee Nickles
##'
setMethod(
    "subtractBackground",
    "RccSet",
    function( rccSet,
              bgEstimates,
              bgEstimatesParams = list(),
              inputMatrix = c("posCtrlData", "exprs"),
              quietly = FALSE )
    {
        if (missing(rccSet) || missing(bgEstimates))
            stop("subtractBackground: missing some arguments")

        inputMatrix <- match.arg(inputMatrix)

        M <- assayData(rccSet)[[ inputMatrix ]]
        if (is.null(M))
            stop(sprintf("Specified input matrix ('%s') is not present in the RccSet's assayData", inputMatrix)) 

        if (!quietly && ("bgEstimates" %in% ls(assayData(rccSet))))
            warning("Input already contains background estimates")
        if (!quietly && ("bgCorrData" %in% ls(assayData(rccSet))))
            warning("Input already contains background corrected data")

        M <- M - bgEstimates

        M[ M < 1 ] <- 1     # Truncate to 1 any counts below 1 (cf. pseudocount in pdata_fdata_adata.to.rccSet())

        srccSet <- copyRccSet(rccSet)

        assayData(srccSet)$bgCorrData <- M
        assayData(srccSet)$bgEstimates <- bgEstimates

        preproc(srccSet)$bgEstimates_bgReference              <- NAifNULL( bgEstimatesParams$bgReference )
        preproc(srccSet)$bgEstimates_summaryFunction          <- NAifNULL( bgEstimatesParams$summaryFunction )
        preproc(srccSet)$bgEstimates_stringency               <- NAifNULL( bgEstimatesParams$stringency )
        preproc(srccSet)$bgEstimates_nSolverBackground.w1     <- NAifNULL( bgEstimatesParams$nSolverBackground.w1 )
        preproc(srccSet)$bgEstimates_nSolverBackground.shrink <- NAifNULL( bgEstimatesParams$nSolverBackground.shrink )
        preproc(srccSet)$bgEstimates_inputMatrix              <- NAifNULL( bgEstimatesParams$inputMatrix )

        return(srccSet)
    }
    )

