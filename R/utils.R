
setGeneric( "copyRccSet", function( rccSet ) standardGeneric( "copyRccSet" ) )

##' @rdname copyRccSet
##' @aliases copyRccSet
##'
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
setMethod(
    "copyRccSet",
    "RccSet",
    function(rccSet)
    {
        RccSet( rccSet ) # Constructor always deep copies
    }
    )



setGeneric( "checkRccSet", function( rccSet, ... ) standardGeneric( "checkRccSet" ) )

##' @rdname checkRccSet
##' @aliases checkRccSet
##'
##' @title
##' Check an RccSet
##'
##' @description
##' Provides additional checks and generates warnings for unexpected or unusual conditions
##' which, though permitted by the RccSet class, may indicate data import errors.
##'
##' @param  rccSet          An RccSet to be checked.
##' @param  reportWarnings  Logical. If TRUE, warnings are reported.
##' @param  showMessages    Logical. If TRUE, notes are shown indicating any optional missing
##'                         columns and the like.
##'
##' @return                 Returns TRUE if no warnings were generated and FALSE otherwise.
##'
##' @examples
##' data(example_rccSet)
##' checkRccSet(example_rccSet)
##'
##' @export
##'
##' @author Robert Ziman
##'
setMethod(
    "checkRccSet",
    "RccSet",
    function(rccSet, reportWarnings=TRUE, showMessages=FALSE)
    {
        result.messages <- character()
        result.warnings <- character()
        result.errors <- character()

        if (is(rccSet, "RccSet"))
        {
            fData.optionalCols <- c(
                 "BarCode"
                ,"ProbeID"
                ,"TargetSeq"
                ,"Species"
               #,"RefSeqMrna"
               #,"GeneID"
               #,"Hgnc_Symbol"
                ,"Comments"
            )
            result.messages <- c(result.messages, check_cols(fData.optionalCols, colnames(fData(rccSet)), "Optional featureData cols not found"))

            pData.optionalCols <- c(
                 "Comments"
                ,"Date"
            )
            result.messages <- c(result.messages, check_cols(pData.optionalCols, colnames(pData(rccSet)), "Optional phenoData cols not found"))

            if ("normData" %in% ls(assayData(rccSet))) {
                pData.qcCols <- c(
                     "TechnicalFlags"
                    ,"ControlFlags"
                    ,"CountFlags"
                    #,"pData_qc"
                )
                result.messages <- c(result.messages, check_cols(pData.qcCols, colnames(pData(rccSet)), "QC-related phenoData cols not found"))
            }

            blankLabel <- getBlankLabel(rccSet, showWarnings=FALSE)
            if (is.na(blankLabel)) {
                result.warnings <- c(result.warnings, "blankLabel not recorded in varMetadata for phenoData SampleType column")
            }

            if (!is.null(pData(rccSet)$FovCount)) {
                if (length(unique(pData(rccSet)$FovCount)) > 1) {
                    result.warnings <- c(result.warnings, "The data seems to have been generated using two different settings (more than one value appears in the FOV Count field)")
                }
            }

            if (!is.null(fData(rccSet)$CodeClass))
            {
                CodeClass.values <- unique(fData(rccSet)$CodeClass)
                CodeClass.nonstandard_values <- setdiff( CodeClass.values, c("Endogenous", "Housekeeping", "Positive", "Negative") )
                if (length(CodeClass.nonstandard_values) > 0) {
                    result.warnings <- c(result.warnings,
                        paste0("Non-standard CodeClass values found in featureData (values currently recognized are \"Endogenous\", \"Housekeeping\", \"Positive\", and \"Negative\"): ", paste(collapse=", ", CodeClass.nonstandard_values), ")")
                    )
                }

                if ("Housekeeping" %in% fData(rccSet)$CodeClass)
                {
                    result.messages <- c(result.messages,
                        sprintf( "The following panel housekeeping genes were found: %s",
                            paste(fData(rccSet)$GeneName[ fData(rccSet)$CodeClass == "Housekeeping" ],
                                  collapse=", "))
                    )
                } else {
                    result.messages <- c(result.messages, "No panel housekeeping genes were found")
                }
            }

            # Check for duplicate accessions. (This can happen when the GeneName field contains
            # something which is not a legitimate gene symbol. A message or warning is appropriate
            # in this case since downstream analysis is unlikely to anticipate multiple probe pairs
            # for the same annotated transcript and will need to handle this specially.)

            if (!is.null(fData(rccSet)$Accession)) {
                if (any(duplicated(fData(rccSet)$Accession))) {
                    result.messages <- c(result.messages, "Duplicate accessions were found")
                }
            }

            #@assayData
            #@annotation

        } else {
            stop("Input is not an RccSet")
        }

        if (showMessages && length(result.messages) > 0) {
            message(paste(collapse="\n", c(
                "checkRccSet() messages:",
                paste0("  ", result.messages)
            )))
        }

        if (reportWarnings && length(result.warnings) > 0) {
            warning(paste(collapse="\n", c(
                "",     # force wordwrap of the first line
                paste0("  ", result.warnings)
            )))
        }

        if (length(result.warnings) > 0)
            return(FALSE)
        else
            return(TRUE)
    }
    )



setGeneric( "getBlankLabel", function( rccSet, ... ) standardGeneric( "getBlankLabel" ) )

##' @rdname getBlankLabel
##' @aliases getBlankLabel
##'
##' @title
##' Get the SampleType value that indicates blank samples
##'
##' @description
##' Returns the phenoData SampleType value that indicates blank samples (i.e. water runs).
##' This value is parsed from the single-quoted string enclosed by "blankLabel='...'"
##' in the varMetadata for SampleType.
##'
##' @param rccSet           An RccSet
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
setMethod(
    "getBlankLabel",
    "RccSet",
    function(rccSet, showWarnings=TRUE)
    {
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
    )

setGeneric( "hasBlanks", function( rccSet ) standardGeneric( "hasBlanks" ) )
setMethod(
    "hasBlanks",
    "RccSet",
    function(rccSet) {
        blankLabel <- getBlankLabel(rccSet)
        any(pData(rccSet)$SampleType %in% blankLabel)
    }
    )

setGeneric( "hasPanelHousekeeping", function( rccSet ) standardGeneric( "hasPanelHousekeeping" ) )
setMethod(
    "hasPanelHousekeeping",
    "RccSet",
    function(rccSet) {
        any(fData(rccSet)$CodeClass == "Housekeeping")
    }
    )

setGeneric( "paStringency", function( rccSet ) standardGeneric( "paStringency" ) )
setMethod(
    "paStringency",
    "RccSet",
    function(rccSet)
    {
        paData_preprocList_value <- preproc(rccSet)$assayData_paData

        if (!is.null(paData_preprocList_value)) {
            m <- regexec(".*stringency factor: ([0-9\\.]*)", paData_preprocList_value)
            unlist(regmatches(paData_preprocList_value, m))[2]
        } else {
            NULL
        }
    }
    )


NAifNULL <- function(x)
{
    if (is.null(x))
        NA
    else
        x
}

