
##' @title
##' Make NanoString QC report
##'
##' @description
##' Creates an html QC report for a NanoString ExpressionSet object. Since the report depends upon
##' correct NanoString-specific annotation in the input ExpressionSet, this function calls
##' validRccSet() to check for any missing columns or other elements. If any problems are found,
##' the function will abort with an error that describes them. Otherwise, the function will return
##' a copy of the input ExpressionSet with additional columns in phenoData showing the QC flags for
##' each sample.
##' 
##' @param rccSet            NanoString ExpressionSet object for which to generate the QC report
##' @param outputBaseName    Scalar character, the base (without extension) to use for various output filenames
##' @param outputDir         Scalar character, the path to the output directory for the QC report and associated files
##' @param experimentTitle   Scalar character, an easy to read identifier of the experiment
##' @param covar             Scalar character, covariate for stratifying samples
##' @param reference         Scalar character, reference to determine limit of detection: "blank" specifies that blank samples
##'                             should be used; "negatives" specifies that negative control
##'                             probes should be used; "auto" specifies that the decision should be automatic (blank samples will be
##'                             used if present, otherwise negative control probes will be used).
##' @param method            Scalar character, method to determine outliers; either "cutoffByVar" or "cutoffByMMAD"
##' @param cutoff            Scalar numeric, cutoff in method to determine outliers
##' @param maxMiss           Scalar numeric, indicating the fraction of genes accepted to be missed
##' @param SampleName        Scalar character, name of the phenoData column holding the sample names
##' @param cleanMarkdown     Should markdown files used to produce QC report be deleted upon completion?
##'
##' @return An html report is written to disk, and a copy of the input ExpressionSet is invisibly
##' returned, with columns added to phenoData that show the QC flags for each sample.
##'
##' @examples
##' data(example_rccSet)
##'
##' pcnorm_example_rccSet <- posCtrlNorm(example_rccSet)
##' bg <- getBackground(pcnorm_example_rccSet)
##' bgcorr_example_rccSet <- subtractBackground(pcnorm_example_rccSet, bg)
##'
##' qc_example_rccSet <- makeQCReport(bgcorr_example_rccSet, "QC_report")
##'
##' @import knitr
##'
##' @export
##'
##' @author Dorothee Nickles, Thomas Sandmann, Robert Ziman
##'
makeQCReport <- function(rccSet,
                         outputBaseName    = "NanoStringQCPro_QC_report",
                         outputDir         = getwd(),
                         experimentTitle   = expinfo(experimentData(rccSet))["title"],
                         covar             = "SampleType",
                         reference         = c("auto", "blank", "negatives"),
                         method            = c("cutoffByMMAD", "cutoffByVar"),
                         cutoff            = 4,
                         maxMiss           = 0.2,
                         SampleName        = "SampleID",
                         cleanMarkdown     = TRUE)
{
  reference       <- match.arg(reference)
  method          <- match.arg(method)

  # The 'rccSet' argument in previous versions of this function was a string specifying the name of the variable containing the object
  # rather than the object itself. The following check is for old scripts that may still be calling the function as such. -RZ 2014-12-20

  if (class(rccSet) == "character") {
    stop("the 'rccSet' argument that the current version of this function now requires is *the actual NanoString ExpressionSet object* rather than the character vector specifying the name of the variable containing it!")
  }

  # Check that the output html filename is properly specified. A directory of matching
  # name will be created next to it to hold extra output files (SampleFlags.txt,
  # the large versions of the QC report images, etc). An .Rmd filename of matching name
  # will be specified as the brew() output that will then be passed to render() the
  # final html.

  outputExtrasDir <- outputBaseName

  fileNames <- sapply(
    c( "Rmd", "md", "html" ),
    function( ext ) sprintf( "%s.%s", outputBaseName, ext ),
    simplify = FALSE
    )
  
  # Create the output dirs if they don't exist.

  if (!file.exists(outputDir))
    if ( !dir.create(outputDir, recursive=TRUE) )
      stop( "outputDir doesn't exist and couldn't be created" )

  oldDir <- getwd()
  on.exit( setwd( oldDir ) )
  setwd( outputDir ) # Fails automatically if not successful
  
  if (!file.exists(outputExtrasDir))
    if ( !dir.create(outputExtrasDir) )
      stop( "outputExtrasDir doesn't exist in outputDir and couldn't be created" )

  # Validate the input rccSet and arguments.

  validRccSet(rccSet, stopOnError=TRUE)

  if (!is.numeric(cutoff) || length(cutoff) != 1)
    stop("cutoff needs to be a scalar numeric value")

  if (!is.numeric(maxMiss) || length(maxMiss) != 1)
    stop("maxMiss needs to be a scalar numeric value")

  if (!(covar %in% colnames(pData(rccSet))))
    stop("covar not found amongst the phenoData columns")

  if (!(SampleName %in% colnames(pData(rccSet))))
    stop("SampleName not found amongst the phenoData columns")

  # Assign any additional variables required by the template.

  blankLabel <- getBlankLabel(rccSet)

  if (reference == "auto") {
    if (hasBlanks(rccSet)) {
      reference <- "blank"
    } else {
      reference <- "negatives"
    }
  }

  report.date <- date()

  gmMainImage    <- file.path( outputExtrasDir, "GeneClustering_median.png" )
  gmPreviewImage <- file.path( outputExtrasDir, "GeneClustering_median_preview.png" )

  # Generate the report

  message("Generating QC report...")

  qcReportHeader <- system.file("misc", "qcReportHeader.html", package="NanoStringQCPro")
  rmdTemplate <- system.file("templates", "QualityControlReportTemplate.Rmd", package="NanoStringQCPro")

  file.copy( rmdTemplate, fileNames$Rmd )

  # Some other potentially relevant arguments that can be passed, which go to markdownToHTML: stylesheet, title.
  knit2html(
    fileNames$Rmd,
    fileNames$html,
    options = c( "use_xhtml", "smartypants", "toc" ), # passed to markdownToHTML. Skipping base64_images, mathjax, and highlight_code
    header = qcReportHeader,
    quiet = TRUE
    )

  if ( cleanMarkdown ) {
    unlink( fileNames$Rmd )
    unlink( fileNames$md )
  }

  message( sprintf( "Report file has been generated: %s", file.path( outputDir, fileNames$html ) ) )

  # Add sample flags and return.

  sampleFlagFile <- file.path( outputExtrasDir, "SampleFlags.txt" )

  if ( !file.exists(sampleFlagFile) )
    stop(paste0("File not found: [", sampleFlagFile, "]"))

  invisible( addQCFlags(rccSet, sampleFlagFile ) )

  # on.exit() above returns us to previous working directory
  
}
