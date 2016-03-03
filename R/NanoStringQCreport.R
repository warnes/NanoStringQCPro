
setGeneric( "makeQCReport", function( rccSet, ... ) standardGeneric( "makeQCReport" ) )

##' @rdname makeQCReport
##' @aliases makeQCReport
##'
##' @title
##' Make NanoString QC report
##'
##' @description
##' Creates an html QC report for an RccSet object. Alongside the html file, a
##' directory with matching filename is produced that contains additional files
##' as well as high resolution versions of the various plots in the report.
##' In addition to generating the QC report, the function returns a copy of the
##' input RccSet with columns added to phenoData that show the QC flags for each
##' sample.
##'
##' @description
##' The various plots in the report depend upon correct annotation and
##' preprocessing in the input object. If the input is missing any elements
##' required for a given plot, the plot will be replaced with a message
##' indicating the missing elements. If the preprocOverride argument is set to
##' TRUE, the input's preprocessing will be ignored and a default configuration
##' will be used so that all applicable plots will be rendered.
##' 
##' @param rccSet
##' RccSet object for which to generate the QC report.
##'
##' @param outputBaseName
##' Character string specifying the base filename (without extension) to use
##' for the output file.
##'
##' @param outputDir
##' Character string specifying the path to the output directory for the QC
##' report and associated files.
##'
##' @param preprocOverride
##' Logical. If TRUE, the input's preprocessing will be ignored, and a default
##' preprocessing configuration (specifically, the defaults for
##' preprocRccSet()) will be applied so that all applicable plots can be
##' rendered in the report.
##'
##' @param experimentTitle
##' Character string specifying an easy to read identifier of the experiment.
##'
##' @param covar
##' Character string specifying a covariate for stratifying samples (e.g.
##' "SampleType").
##'
##' @param method
##' Method to determine outlier samples: either "cutoffByVar" or
##' "cutoffByMMAD".
##'
##' @param stringency
##' Multiplier with which to adjust cutoff values for determining outlier
##' samples.
##'
##' @param maxMiss
##' Numeric specifying the allowable fraction of genes below the lower limit
##' of detection in a sample.
##'
##' @param sampleNameCol
##' Character string specifying the name of the phenoData column holding the
##' sample names.
##'
##' @param heatmaps
##' Logical: render and show heatmaps?
##'
##' @param cleanMarkdown
##' Logical: upon completion, delete markdown files used to produce QC report?
##'
##' @param verbose
##' Logical: print progress messages?
##'
##' @return
##' An html report is written to disk and a copy of the input RccSet is
##' invisibly returned with columns added to phenoData that show the QC flags
##' for each sample.
##'
##' @examples
##' data(example_rccSet)
##' norm_example_rccSet <- preprocRccSet(example_rccSet)
##' qc_example_rccSet <- makeQCReport(norm_example_rccSet, "example_QC_report")
##'
##' @import knitr
##'
##' @export
##'
##' @author Dorothee Nickles, Thomas Sandmann, Robert Ziman, Richard Bourgon
##'
setMethod(
  "makeQCReport",
  "RccSet",
  function(rccSet,
           outputBaseName    = "NanoStringQCPro_QC_report",
           outputDir         = getwd(),
           preprocOverride   = FALSE,
           experimentTitle   = expinfo(experimentData(rccSet))["title"],
           covar             = "SampleType",
           method            = c("cutoffByMMAD", "cutoffByVar"),
           stringency        = 4,
           maxMiss           = 0.2,
           sampleNameCol     = "SampleID",
           heatmaps          = FALSE,
           cleanMarkdown     = TRUE,
           verbose           = FALSE)
  {
    method <- match.arg(method)

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
  
    # Validate arguments and check the input rccSet
  
    checkRccSet(rccSet, reportWarnings=TRUE, showMessages=FALSE)
  
    if (!is.numeric(stringency) || length(stringency) != 1)
      stop("stringency needs to be a scalar numeric value")
  
    if (!is.numeric(maxMiss) || length(maxMiss) != 1)
      stop("maxMiss needs to be a scalar numeric value")
  
    if (!(covar %in% colnames(pData(rccSet))))
      stop("covar not found amongst the phenoData columns")
  
    if (!(sampleNameCol %in% colnames(pData(rccSet))))
      stop(sprintf("sampleNameCol ('%s') not found amongst the phenoData columns", sampleNameCol))
  
    # Assign any additional variables required by the template.
  
    blankLabel <- getBlankLabel(rccSet)

    scMainImage         <- file.path( outputExtrasDir, "SampleClustering.png" )
    scPreviewImage      <- file.path( outputExtrasDir, "SampleClustering_preview.png" ) 
    gcibsMainImage      <- file.path( outputExtrasDir, "GenewiseCountsInBlankSamples.png" )
    gcibsPreviewImage   <- file.path( outputExtrasDir, "GenewiseCountsInBlankSamples_preview.png" )
    ncblMainImage       <- file.path( outputExtrasDir, "NegativeControlsByLane.png" )
    ncblPreviewImage    <- file.path( outputExtrasDir, "NegativeControlsByLane_preview.png" )
    gcgMainImage        <- file.path( outputExtrasDir, "GeneClusteringGlobal.png" )
    gcgPreviewImage     <- file.path( outputExtrasDir, "GeneClusteringGlobal_preview.png" )
    gchMainImage        <- file.path( outputExtrasDir, "GeneClusteringHousekeeping.png" )
    gchPreviewImage     <- file.path( outputExtrasDir, "GeneClusteringHousekeeping_preview.png" )
    cscMainImage        <- file.path( outputExtrasDir, "CenteredSampleClustering.png" )
    cscPreviewImage     <- file.path( outputExtrasDir, "CenteredSampleClustering_preview.png" )
    cscxfsMainImage     <- file.path( outputExtrasDir, "CenteredSampleClusteringExcludingFlaggedSamples.png" )
    cscxfsPreviewImage  <- file.path( outputExtrasDir, "CenteredSampleClusteringExcludingFlaggedSamples_preview.png" )

    # If specified, override the RccSet's preprocessing with a default
    # configuration. This ensures the availability of assayData matrices
    # and other elements required so that all applicable QC plots will be shown
    # in the report.

    if (preprocOverride) {
        rccSet <- preprocRccSet(rccSet)
    }

    # Generate the report
  
    message("Generating QC report...")

    if (verbose) message("  Gene-wise counts in blank samples: main image...")
    countsInBlankSamples_verticalPlot(rccSet, gcibsMainImage)

    if (verbose) message("  Gene-wise counts in blank samples: preview image...")
    previewPNG(gcibsMainImage, gcibsPreviewImage, width=1500, cropHeight=600, res=300)

    if (verbose) message("  Negative controls by lane: main image...")
    negCtrlsByLane_verticalPlot(rccSet, ncblMainImage)

    if (verbose) message("  Negative controls by lane: preview image...")
    previewPNG(ncblMainImage, ncblPreviewImage, width=1500, cropHeight=600, res=300)

    if (verbose) message("  Raw and content normalized data for density plots and gene clustering...")

    rawData <- assayData(rccSet)$exprs

    if ("bgCorrData" %in% ls(assayData(rccSet)))
      normInputMatrix <- "bgCorrData"
    else if ("posCtrlData" %in% ls(assayData(rccSet)))
      normInputMatrix <- "posCtrlData"
    else
      normInputMatrix <- "exprs"

    if (preproc(rccSet)$normData_method %in% "global") {    # Don't use == here (normData_method can be NA)
      gmrccSet <- copyRccSet(rccSet)
      gmnormData <- assayData(gmrccSet)$normData
    } else {
      gmrccSet <- contentNorm(rccSet,
                              method = "global",
                              summaryFunction = "median",
                              inputMatrix = normInputMatrix,
                              quietly = TRUE)
      gmnormData <- assayData(gmrccSet)$normData
    }

    if (preproc(rccSet)$normData_method %in% "housekeeping") {  # Don't use == here (normData_method can be NA)
      hkrccSet <- copyRccSet(rccSet)
      hknormData <- assayData(hkrccSet)$normData
    } else if (hasPanelHousekeeping(rccSet)) {
      hkrccSet <- contentNorm(rccSet,
                              method="housekeeping",
                              summaryFunction = "median",
                              inputMatrix = normInputMatrix,
                              quietly=TRUE)
      hknormData <- assayData(hkrccSet)$normData
    } else {
      hkrccSet <- NULL
      hknormData <- NULL
    }

    flags <- rep(FALSE, ncol(rccSet)*3)
    dim(flags) <- c(ncol(rccSet), 3)
    rownames(flags) <- colnames(rccSet)
    colnames(flags) <- c("TechnicalFlags", "ControlFlags", "CountFlags")
    flags[,"TechnicalFlags"][flagSamplesTech(rccSet)] <- TRUE
    flags[,"ControlFlags"][flagSamplesCtrl(rccSet, method, stringency)] <- TRUE
    flags[,"CountFlags"][flagSamplesCount(rccSet, method, stringency, maxMiss)] <- TRUE
    flags2 <- cbind(SampleIdentifier=rownames(flags),
                    SampleName=as.character(pData(rccSet)[, sampleNameCol]),
                    flags)
    flags2 <- cbind(flags2, SampleType=pData(rccSet)$SampleType)
    write.table(flags2, file=file.path(outputExtrasDir, "SampleFlags.txt"), sep="\t", row.names=FALSE, quote=FALSE)

    if (heatmaps) {

        if (verbose) message("  Heatmaps...")

        if (verbose) message("    Sample clustering: main image...")
        sampleClustering(rccSet     = rccSet,
                         outputFile = scMainImage,
                         covar      = covar)
        
        if (verbose) message("    Sample clustering: preview image...")
        previewPNG(scMainImage, scPreviewImage, width=1500, cropHeight=1500, res=300)

        if (verbose) message("    Gene clustering (global normalization): main image...")
        geneClustering(rccSet     = gmrccSet,
                       outputFile = gcgMainImage,
                       main       = "Gene clustering (global normalization)")

        if (verbose) message("    Gene clustering (global normalization): preview image...")
        previewPNG(gcgMainImage, gcgPreviewImage, width=1500, cropHeight=1500, res=300)

        if (!is.null(hkrccSet)) {

            if (verbose) message("    Gene clustering (housekeeping normalization): main image...")
            geneClustering(rccSet     = hkrccSet,
                           outputFile = gchMainImage,
                           main       = "Gene clustering (housekeeping normalization)")

            if (verbose) message("    Gene clustering (housekeeping normalization): preview image...")
            previewPNG(gchMainImage, gchPreviewImage, width=1500, cropHeight=1500, res=300)

        } else {
            if (verbose) message("    Gene clustering (housekeeping normalization): main image... skipped")
            if (verbose) message("    Gene clustering (housekeeping normalization): preview image... skipped")
        }

        if (verbose) message("    Centered sample clustering: main image...")
        centeredSampleClustering(gmrccSet,
                                 outputFile = cscMainImage,
                                 main = "Sample clustering of all samples run",
                                 flags = flags,
                                 excludeFlagged = FALSE)

        if (verbose) message("    Centered sample clustering: preview image...")
        previewPNG(cscMainImage, cscPreviewImage, width=1500, cropHeight=1500, res=300)

        if (verbose) message("    Centered sample clustering (excluding flagged samples): main image...")
        centeredSampleClustering(gmrccSet,
                                 outputFile = cscxfsMainImage,
                                 main = "Sample clustering after exclusion of outlier samples",
                                 flags = flags,
                                 excludeFlagged = TRUE)

        if (verbose) message("    Centered sample clustering (excluding flagged samples): preview image...")
        previewPNG(cscxfsMainImage, cscxfsPreviewImage, width=1500, cropHeight=1500, res=300)

    } else {
        if (verbose) message("  Heatmaps... skipped")
    }
 
    report.date <- date()
    qcReportHeader <- system.file("misc", "qcReportHeader.html", package="NanoStringQCPro")
    rmdTemplate <- system.file("templates", "QualityControlReportTemplate.Rmd", package="NanoStringQCPro")
 
    file.copy( rmdTemplate, fileNames$Rmd, overwrite = TRUE )
  
    # Some other potentially relevant arguments that can be passed, which go to markdownToHTML: stylesheet, title.
    knit2html(
      fileNames$Rmd,
      fileNames$html,
      options = c( "use_xhtml", "smartypants", "toc" ), # passed to markdownToHTML. Skipping base64_images, mathjax, and highlight_code
      header = qcReportHeader,
      quiet = TRUE,
      force_v1 = TRUE
        #
        # Fixes the "It seems you should call rmarkdown::render() instead of knitr::knit2html() because NanoStringQCPro_QC_report.Rmd
        # appears to be an R Markdown v2 document" error that started coming up with knitr_1.12 (no such error was coming up with knitr_1.11).
        # -RZ 2016-03-03
        #
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
    rccSet <- addQCFlags(rccSet, sampleFlagFile)
  
    invisible(rccSet)
  
    # on.exit() above returns us to previous working directory
    
  }
  )

