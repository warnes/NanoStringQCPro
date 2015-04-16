
##' @title
##' Deep-copy a NanoString RccSet
##'
##' @description
##' 
##' Returns a copy of the input RccSet where the copy's assayData has been
##' produced via copyEnv() rather than a simple assignment -- hence deep-copying
##' the environment pointed to by assayData rather than just copying the
##' pointer. This guarantees that if the copy's assayData is affected later in
##' the code, assayData for the original won't be affected.
##'
##' @param rccSet A NanoString RccSet to be copied.
##'
##' @return
##'
##' A new RccSet that is a deep copy of the original.
##'
##' @export
##'
##' @examples
##' data(example_rccSet)
##' example_rccSet_2 <- copyRccSet(example_rccSet)
##' assayData(example_rccSet)
##' assayData(example_rccSet_2) # Should be different
##'
##' @author Robert Ziman
##'
copyRccSet <- function(rccSet)
{
    stopifnot( is( rccSet, "RccSet") )
    RccSet( rccSet ) # Constructor always deep copies
}

##' @title
##' Validate a NanoString ExpressionSet
##'
##' @description
##' Checks to see if the given ExpressionSet contains all the columns and other
##' member elements that should be present if the ExpressionSet contains
##' NanoString mRNA assay data.
##'
##' @param  rccSet          A NanoString ExpressionSet to be validated.
##' @param  stopOnError     Logical. If TRUE and validity fails, the function generates
##'                         an error message. If FALSE and validity fails, the function
##'                         returns a character vector with the errors found.
##' @param  reportWarnings  Logical. If TRUE, warnings are reported.
##' @param  showMessages    Logical. If TRUE, notes are shown indicating any optional missing
##'                         columns and the like.
##'
##' @return                 Returns TRUE if the input object is an ExpressionSet with NanoString
##'                         data. If problems are found and stopOnError=TRUE, the function
##'                         will generate an error message. If problems are found and
##'                         stopOnError=FALSE, the function will return a character vector with
##'                         the errors found.
##'
##' @examples
##' data(example_rccSet)
##' validRccSet(example_rccSet)
##' fData(example_rccSet)$CodeClass <- NULL
##' #validRccSet(example_rccSet)    # returns an error
##'
##' @export
##'
##' @author Robert Ziman
##'
validRccSet <- function(rccSet, stopOnError=TRUE, reportWarnings=TRUE, showMessages=FALSE)
{
    result.messages <- character()
    result.warnings <- character()
    result.errors <- character()

    if (is(rccSet, "ExpressionSet"))
    {
        fData.requiredCols <- c(
             "CodeClass"
            ,"GeneName"
            ,"Accession"
            ,"BarCode"
            ,"ProbeID"
        )
        fData.optionalCols <- c(
             "TargetSeq"
            ,"Species"
           #,"RefSeqMrna"
           #,"GeneID"
           #,"Hgnc_Symbol"
            ,"Comments"
        )

        pData.requiredCols <- c(
             "FileName"
            ,"SampleID"
            ,"LaneID"
            ,"FovCount"
            ,"FovCounted"
            ,"StagePosition"
            ,"BindingDensity"
            ,"CartridgeID"
            ,"SampleType"
        )
        #pData.recommendedCols <- c(
        #)
        pData.optionalCols <- c(
             "Comments"
            ,"Date"
        )
        pData.qcCols <- c(
             "TechnicalFlags"
            ,"ControlFlags"
            ,"CountFlags"
            ,"pData_qc"
        )

        check_cols <- function(expectedCols, actualCols, msg)
        {
            notfound <- is.na(match(expectedCols, actualCols))
            if (any(notfound)) {
                return(
                    paste0(msg, ": ", paste(collapse=", ", expectedCols[notfound]))
                )
            }
            else {
                return(character(0))
            }
        }

        result.messages <- c(result.messages, check_cols(fData.optionalCols, colnames(fData(rccSet)), "Optional featureData cols not found"))
        result.messages <- c(result.messages, check_cols(pData.optionalCols, colnames(pData(rccSet)), "Optional phenoData cols not found"))

        #result.warnings <- c(result.warnings, check_cols(pData.recommendedCols, colnames(pData(rccSet)), "Recommended pData cols not found"))

        result.errors <- c(result.errors, check_cols(fData.requiredCols, colnames(fData(rccSet)), "Required featureData cols not found"))
        result.errors <- c(result.errors, check_cols(pData.requiredCols, colnames(pData(rccSet)), "Required phenoData cols not found"))

        blankLabel <- getBlankLabel(rccSet, showWarnings=FALSE)
        if (is.na(blankLabel)) {
            result.warnings <- c(result.warnings, "blankLabel not recorded in varMetadata for phenoData SampleType column")
        }

        if (!is.null(pData(rccSet)$FovCount)) {
            if (length(unique(pData(rccSet)$FovCount)) > 1) {
                result.errors <- c(result.errors, "The data seems to have been generated using two different settings (more than one value appears in the FOV Count field)")
            }
        }

        if (!is.null(fData(rccSet)$CodeClass))
        {
            CodeClass.values <- unique(fData(rccSet)$CodeClass)
            CodeClass.nonstandard_values <- grep(invert=TRUE, value=TRUE, "Endogenous", CodeClass.values)
            CodeClass.nonstandard_values <- grep(invert=TRUE, value=TRUE, "Housekeeping", CodeClass.nonstandard_values)
            CodeClass.nonstandard_values <- grep(invert=TRUE, value=TRUE, "Positive", CodeClass.nonstandard_values)
            CodeClass.nonstandard_values <- grep(invert=TRUE, value=TRUE, "Negative", CodeClass.nonstandard_values)
            if (length(CodeClass.nonstandard_values) > 0) {
                result.warnings <- c(result.warnings,
                    paste0("Non-standard CodeClass values found in featureData (values currently recognized are \"Endogenous\", \"Housekeeping\", \"Positive\", and \"Negative\"): ", paste(collapse=", ", CodeClass.nonstandard_values), ")")
                )
            }

            if ("Housekeeping" %in% fData(rccSet)$CodeClass)
            {
                result.messages <- c(result.messages,
                    sprintf( "The following panel housekeeping genes were found: %s",
                        paste(collapse=",", fData(rccSet)$GeneName[ fData(rccSet)$CodeClass == "Housekeeping" ]) )
                )
            } else {
                result.messages <- c(result.messages, "No panel housekeeping genes were found")
            }
        }

        if (!is.null(fData(rccSet)$Accession)) {
            if (any(duplicated(fData(rccSet)$Accession))) {
                result.messages <- c(result.messages, "Duplicate accessions were found")
            }
        }

        #@assayData
        #@annotation

    } else {
        result.errors <- "Input is not an ExpressionSet"
    }

    if (showMessages && length(result.messages) > 0) {
        message(paste(collapse="\n", c(
            "validRccSet() messages:",
            paste0("  ", result.messages)
        )))
    }

    if (reportWarnings && length(result.warnings) > 0) {
        warning(paste(collapse="\n", c(
            "",     # force wordwrap of the first line
            paste0("  ", result.warnings)
        )))
    }

    if (length(result.errors) > 0)
    {
        if (stopOnError) {
            indent <- c("", rep("  ", length(result.errors)-1))
            stop(paste(collapse="\n",
                paste0(indent, result.errors)
            ))
        } else {
            return(result.errors)
        }
    }
    else {
        return(TRUE)
    }
}

##' @title
##' Get the SampleType value that indicates blank samples
##'
##' @description
##' Returns the phenoData SampleType value that indicates blank samples (i.e. water runs).
##' This value is parsed from the single-quoted string enclosed by "blankLabel='...'"
##' in the varMetadata for SampleType.
##'
##' @param rccSet           A NanoString Expressionset
##' @param showWarnings     Logical. If FALSE, no warnings will be generated (if any).
##'
##' @return
##' NULL if the SampleType column is missing altogether, NA if the varMetadata doesn't
##' have blankLabel recorded, or the blankLabel value otherwise.
##'
##' @examples
##' data(example_rccSet)
##' blankLabel <- getBlankLabel(example_rccSet)
##'
##' @export
##'
##' @author Robert Ziman
##'
getBlankLabel <- function(rccSet, showWarnings=TRUE)
{
    stopifnot(is(rccSet, "ExpressionSet"))

    SampleType.rownum <- which(rownames(varMetadata(phenoData(rccSet))) == "SampleType")
    if (length(SampleType.rownum) == 0) {   # "SampleType" was not found amongst the rownames
        if (showWarnings) {
            warning("phenoData doesn't have a 'SampleType' column")
        }
        return(NULL)
    } else {
        SampleType.metadata <- varMetadata(phenoData(rccSet))$labelDescription[ SampleType.rownum ]
        if (!grepl("blankLabel='.*'", SampleType.metadata)) {
            if (showWarnings) {
                warning("phenoData 'SampleType' varMetadata doesn't contain blankLabel string")
            }
            return(NA)
        } else {
            blankLabel <- sub("'.*", "", sub(".*blankLabel='", "", SampleType.metadata))
            return (blankLabel)
        }
    }
}

hasBlanks <- function(rccSet) {
    stopifnot(is(rccSet, "ExpressionSet"))
    blankLabel <- getBlankLabel(rccSet)
    any(pData(rccSet)$SampleType %in% blankLabel)
}

hasPanelHousekeeping <- function(rccSet) {
    stopifnot(is(rccSet, "ExpressionSet"))
    any(fData(rccSet)$CodeClass == "Housekeeping")
}

preprocState <- function(rccSet)
{
    stopifnot(is(rccSet, "ExpressionSet"))

    preprocList <- preproc(rccSet)

    if (is.null(preprocList$exprs)) {
        stop("Malformed preproc list (preproc(rccSet)$exprs is NULL)")
    } else if (preprocList$exprs == "Raw data") {                                                   state <- "newRccSet"
    } else if (preprocList$exprs == "Positive control normalized data") {                           state <- "posCtrlNorm"
    } else if (preprocList$exprs == "Background corrected data") {                                  state <- "subtractBackground"
    } else if (preprocList$exprs == "Positive control normalized and background corrected data") {  state <- "subtractBackground"
    } else if (preprocList$exprs == "Preprocessed and median-normalized data") {                    state <- "preprocRccSet"
    } else if (preprocList$exprs == "Preprocessed and mean-normalized data") {                      state <- "preprocRccSet"
    } else if (preprocList$exprs == "Preprocessed and housekeeping-normalized data") {              state <- "preprocRccSet"
    } else {
        stop("Unrecognized preproc state")
    }

    return(state)
}

