
check_cols <- function(expectedCols, actualCols, msg) {
    notfound <- is.na(match(expectedCols, actualCols))
    if (any(notfound)) 
        return(paste0(msg, ": ", paste(collapse=", ", expectedCols[notfound])))
    else return(character(0))
}


##' @title
##'
##' RccSet class, derived from ExpressionSet
##'
##' @description
##'
##' The RccSet class is a trivial extension of
##' \code{\linkS4class{ExpressionSet}}, but with additional validation
##' criteria. \code{RccSet} is a class generator function.
##'
##' @details
##'
##' A valid \code{RccSet} object must have the following columns in
##' \code{featureData}: \code{"CodeClass"}, \code{"GeneName"}, and
##' \code{"Accession"}. It must also have the following
##' \code{phenoData} columns: \code{"FileName"}, \code{"SampleID"},
##' \code{"LaneID"}, \code{"FovCount"}, \code{"FovCounted"},
##' \code{"StagePosition"}, \code{"BindingDensity"}, \code{"CartridgeID"}, and
##' \code{"SampleType"}. A final requirement is that the \code{"FovCount"}
##' column of \code{phenoData} have at most one distinct value.
##'
##' @seealso
##'
##' See \code{\link{checkRccSet}}, which provides additional checks and
##' generates warnings for unexpected or unusual conditions which, though
##' permitted by the class, may indicate data import errors.
##' 
##' @export
##'
##' @import methods
##'
##' @examples
##' 
##' data("example_rccSet")
##' e <- example_rccSet
##' 
##' # "ExpressionSet" constructor makes a new assayData environment
##' r1 <- RccSet(e)
##' validObject(r1)
##' assayData(e)
##' assayData(r1)
##' head(pData(r1))
##' head(fData(r1))
##' 
##' # For other constructors, if not explicitly supplied, blank phenoData and
##' # featureData objects are populated with mandatory columns (and NA values).
##' r2 <- RccSet(assayData(e))
##' validObject(r2)
##' head(pData(r2))
##' head(fData(r2))
##' 
##' r3 <- RccSet(assayData(e), phenoData(e), featureData(e))
##' identical(pData(r1), pData(r3))
##' identical(fData(r1), fData(r3))
##' identical(annotation(r1), annotation(r3)) # We forgot it!
##' annotation(e)
##' r3 <- RccSet(assayData(e), phenoData(e), featureData(e), annotation = annotation(e))
##' identical(annotation(r1), annotation(r3)) # Better
##' identical(r1, r3) # False, due to assayData environments
##' assayData(r1)
##' assayData(r3)
##' 
##' # Matrix contructor is similar
##' r4 <- RccSet(exprs(e), phenoData(e), featureData(e), annotation = annotation(e))
##' identical(exprs(r1), exprs(r4))
##' 
##' # Blank object constructor
##' r0 <- RccSet()
##' dim(r0)
##' pData(r0)
##' fData(r0)

.RccSet <- setClass(
    "RccSet",
    contains = "ExpressionSet",
    validity = function( object ) {
        fData.requiredCols <- c(
            "CodeClass"
            ,"GeneName"
            ,"Accession"
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
        result.errors <- check_cols(fData.requiredCols, colnames(fData(object)), "Required featureData cols not found")
        result.errors <- c(result.errors, check_cols(pData.requiredCols, colnames(pData(object)), "Required phenoData cols not found"))        
        if (length(unique(pData(object)$FovCount)) > 1)
            result.errors <- c(result.errors, "The data seems to have been generated using two different settings (more than one value appears in the 'FovCount' field)")
        if (length(result.errors) > 0) {
            indent <- c("", rep("  ", length(result.errors)-1))
            stop(paste(collapse="\n", paste0(indent, result.errors)))
        }
        else return(TRUE)            
    }
    )

##' @title
##'
##' RccSet constructor methods
##'
##' @description
##'
##' Constructor methods for making new \code{\linkS4class{RccSet}} objects.
##'
##' @param obj An object of appropriate class
##' @param ... Passed to methods.
##'
##' @return
##' A new \code{\linkS4class{RccSet}} object.
##'
##' @details
##'
##' Arguments accepted by constructors are identical to those for the
##' \code{\link{ExpressionSet}} constructors.
##'
##' See \code{\linkS4class{RccSet}} class documentation for examples of
##' constructor use.
##'
##' Constructor calls for which mandatory \code{phenoData} or \code{featureData}
##' columns are missing will successfully create objects that include mandatory
##' columns, but with \code{NA} values. See \code{\linkS4class{RccSet}}
##' documentation for a list of mandatory columns.
##'
##' @export

setGeneric( "RccSet", function( obj, ... ) standardGeneric( "RccSet" ) )

##' @rdname RccSet
##' @export

setMethod(
    "RccSet",
    "ExpressionSet",
    function( obj, ... ) {
        .RccSet(
            assayData = copyEnv( assayData( obj ) ),
            phenoData = phenoData( obj ),
            featureData = featureData( obj ),
            experimentData = experimentData( obj ),
            annotation = annotation( obj ),
            protocolData = protocolData( obj ),
            ...
            )
    }
    )

##' @rdname RccSet
##' @export

setMethod(
    "RccSet",
    "environment",
    function( obj, ... ) {

        e <- ExpressionSet( obj, ... )

        fData.char <- c( "CodeClass", "GeneName", "Accession" )
        pData.char <- c( "FileName", "SampleID", "CartridgeID", "SampleType" )
        pData.numeric <- c( "LaneID", "FovCount", "FovCounted", "StagePosition", "BindingDensity" )

        for ( f in fData.char )
            if ( !( f %in% colnames( fData( e ) ) ) )
                fData( e )[[ f ]] <- rep( NA_character_, nrow( e ) ) # Be sure to handle length 0!

        for ( p in pData.char )
            if ( !( p %in% colnames( pData( e ) ) ) )
                pData( e )[[ p ]] <- rep( NA_character_, ncol( e ) )

        for ( p in pData.numeric )
            if ( !( p %in% colnames( pData( e ) ) ) )
                pData( e )[[ p ]] <- rep( NA_real_, ncol( e ) )

        RccSet( e )

    }
)

##' @rdname RccSet
##' @export

setMethod(
    "RccSet",
    "matrix",
    function( obj, ... ) RccSet( assayDataNew( exprs = obj ), ... )
    )

##' @rdname RccSet
##' @export

setMethod(
    "RccSet",
    "missing",
    function( obj, ... ) RccSet( matrix( nrow = 0, ncol = 0 ), ... )
    )
