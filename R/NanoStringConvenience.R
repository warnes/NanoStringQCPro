
## Data import and normalization functions -- and some small convenience functions

##' @title          Read an .RCC file
##' @description    Parse an .RCC file into a list with each part of the file (Header,
##'                 Sample_Attributes, Lane_Attributes, Code_Summary, etc) stored as a vector
##'                 or data frame. (Note: for reading a set of .RCC files together, use
##'                 readRccBatch().)
##'
##' @param  rccPath             Path to the .RCC file.
##' @param  removeSpikeInLabels Logical. If TRUE (the default), RNA ``spike-in'' input labels (if any)
##'                             in the GeneName for positive and negative control probes will be removed.
##'
##' @return
##' A list where each element holds the contents of one part of the .RCC file
##' (Header, Sample_Attributes, Lane_Attributes, Code_Summary, etc) as a
##' vector or data frame.
##'
##' @export
##'
##' @examples
##' rccPath <- system.file("extdata", "RCC", "20140604_C1-unstim_C1-unstim_01.RCC", package="NanoStringQCPro")
##' rcc <- readRcc(rccPath)
##'
##' @author Robert Ziman
##'
readRcc <- function(rccPath, removeSpikeInLabels=TRUE)
{
	suppressWarnings(
		lines <- readLines(rccPath)
	)
	lines.df <-
	    data.frame(stringsAsFactors=FALSE,
	         line   = lines
	        ,tag    = rep("", length(lines))
	    )

	tags <- c("Header", "Sample_Attributes", "Lane_Attributes", "Code_Summary", "Messages")
	for(i in 1:length(tags))
	{
	    tag_start <- which(lines.df$line == paste0("<", tags[i], ">"))
	    tag_end <- which(lines.df$line == paste0("</", tags[i], ">"))
	    lines.df$tag[ tag_start:tag_end ] <- tags[i]
	}

	unrecognized_lines <- lines.df$line[ lines.df$tag == "" ]
	unrecognized_tags <- grep(value=TRUE, "^<", unrecognized_lines)
	unrecognized_tags <- grep(invert=TRUE, value=TRUE, "^</", unrecognized_tags)
	if (length(unrecognized_tags) > 0) {
	    warning(paste0("Unrecognized tags: ", paste(collapse=" ", unrecognized_tags)))
	}

	lines.df <- lines.df[ (lines.df$tag != ""), ]
	lines.df <- lines.df[ !grepl("<", lines.df$line), ]
	lines.df <- lines.df[ !grepl("</", lines.df$line), ]

	ID_lnum <- grep("^ID,", lines.df$line)
	lines.df$line[ID_lnum] <- paste0( sub("_Attributes", "", lines.df$tag[ID_lnum]), lines.df$line[ID_lnum] )

	rcc <- list()

	rcc$File_path <- rccPath

	for (attr_tag in c("Header", "Sample_Attributes", "Lane_Attributes")) {
		attr_lnum <- which(lines.df$tag == attr_tag)
		attr.strsplit <- strsplit(lines.df$line[attr_lnum], ",")
		attr.v <- vapply(attr.strsplit, function(x) {if (!is.na(x[2])) {x[2]} else {""}}, character(1))
		names(attr.v) <- vapply(attr.strsplit, function(x) {x[1]}, character(1))
		rcc[[attr_tag]] <- attr.v
	}

	Code_Summary_lnum.all <- which(lines.df$tag == "Code_Summary")
	Code_Summary_lnum.header <- Code_Summary_lnum.all[1]
	Code_Summary_lnum.body <- Code_Summary_lnum.all[2:length(Code_Summary_lnum.all)]

	if (lines.df$line[Code_Summary_lnum.header] != "CodeClass,Name,Accession,Count") {
	    stop(
	        paste0("Unrecognized header line for Code_Summary section of \"", rccPath, "\"",
	            " (expected \"CodeClass,Name,Accession,Count\" but got \"", lines.df$line[Code_Summary_lnum.header], "\")")
	    )
	}

	Code_Summary_fields <- strsplit(lines.df$line[Code_Summary_lnum.body], ",")
	rcc$Code_Summary <- data.frame(stringsAsFactors=FALSE
		,CodeClass   = vapply(Code_Summary_fields, function(x) {x[1]}, character(1))
		,Name        = vapply(Code_Summary_fields, function(x) {x[2]}, character(1))
		,Accession   = vapply(Code_Summary_fields, function(x) {x[3]}, character(1))
		,Count       = vapply(Code_Summary_fields, function(x) {x[4]}, character(1))
	)

    if(removeSpikeInLabels)
    {
        spikein <- getSpikeInInput(CodeClass=rcc$Code_Summary$CodeClass, GeneName=rcc$Code_Summary$Name)
        rcc$Code_Summary$Name <- spikein$GeneName
    }

	return(rcc)
}

##' @title          rccDir.to.pdata_fdata_adata
##' @description    First stage of readRccBatch(): produces a list containing matrices (for pdata and adata) and a data frame (for fdata) that
##'                 pdata_fdata_adata.to.rccSet() then transforms into a full ExpressionSet (after some further checks and adjustments).
##'                 See also nSolverCsv.to.pdata_fdata_adata().
##'
##' @param  rccDir  Directory containing .RCC files
##'
##' @return
##' A list containing matrices (for pdata and adata) and a data frame (for fdata)
##' that pdata_fdata_adata.to.rccSet() then transforms into a full ExpessionSet.
##'
##' @author Robert Ziman
##'
rccDir.to.pdata_fdata_adata <- function(rccDir)
{
	rccFiles <- list.files(path=rccDir, pattern="*\\.RCC$")
	if (length(rccFiles) == 0) {
		stop( paste0("No RCC files found in \"", rccDir, "\"") )
	}

	rcc_1 <- readRcc( file.path(rccDir, rccFiles[1]) )

    filename <- rccFiles[1]
    names(filename) <- "FileName"
	rcc_1.pdata <- c(filename, rcc_1$Header, rcc_1$Sample_Attributes, rcc_1$Lane_Attributes)
	rcc_1.fdata <- data.frame(stringsAsFactors=FALSE
		,CodeClass   = rcc_1$Code_Summary$CodeClass
		,GeneName    = rcc_1$Code_Summary$Name
		,Accession   = rcc_1$Code_Summary$Accession
	)
	rcc_1.adata <- as.integer(rcc_1$Code_Summary$Count)

	pdata.m <- matrix(NA, nrow=length(rccFiles), ncol=length(rcc_1.pdata))
	colnames(pdata.m) <- names(rcc_1.pdata)
	pdata.m[1, ] <- rcc_1.pdata

	fdata.df <- rcc_1.fdata

	adata.m.t <- matrix(NA, nrow=length(rccFiles), ncol=nrow(rcc_1.fdata))
	adata.m.t[1, ] <- rcc_1.adata

	for (i in 2:length(rccFiles))
	{
		rcc <- readRcc( file.path(rccDir, rccFiles[i]) )

        filename <- rccFiles[i]
        names(filename) <- "FileName"
		rcc.pdata <- c(filename, rcc$Header, rcc$Sample_Attributes, rcc$Lane_Attributes)
		rcc.fdata <- data.frame(stringsAsFactors=FALSE
			,CodeClass   = rcc$Code_Summary$CodeClass
			,GeneName    = rcc$Code_Summary$Name
			,Accession   = rcc$Code_Summary$Accession
		)
		rcc.adata <- as.numeric(rcc$Code_Summary$Count)

		if (!identical(names(rcc.pdata), names(rcc_1.pdata))) {
			stop(
				paste0("The attribute labels in the Header, Sample_Attributes, and Lane_Attributes sections of \"", rccFiles[i],
					"\" don't exactly match those in the other RCC files parsed so far in the given directory (\"", rccDir, "\"")
			)
		}
		if (!identical(rcc.fdata, rcc_1.fdata)) {
			stop(
				paste0("Values in the key cols (i.e. CodeClass, Name, and Accession) in the Code_Summary section of \"", rccFiles[i],
					"\" don't exactly match those in the other RCC files parsed so far in the given directory (\"", rccDir, "\"")
			)
		}

		pdata.m[i, ] <- rcc.pdata
		adata.m.t[i, ] <- rcc.adata
	}

    return( list(pdata.m=pdata.m, fdata.df=fdata.df, adata.m=t(adata.m.t)) )
}

##' @title          nSolverCsv.to.pdata_fdata_fdata
##' @description    First stage of readRccCollectorToolExport(): produces a list containing
##'                 matrices (for pdata and adata) and a data frame (for fdata) that
##'                 pdata_fdata_adata.to.rccSet then transforms into a full ExpressionSet
##'                 (after some further checks and adjustments). Not intended for
##'                 external use; see also rccDir.to.pdata_fdata_adata().
##'
##' @param  rccCollectorToolExport  Path to the nSolver RCC Collector Tool .csv export.
##'
##' @return
##' A list containing matrices (for pdata and adata) and a data frame (for fdata)
##' that pdata_fdata_adata.to.rccSet() then transforms into a full ExpessionSet.
##'
##' @author Dorothee Nickles, Thomas Sandmann, Robert Ziman
##'
nSolverCsv.to.pdata_fdata_adata <- function(rccCollectorToolExport)
{
	csv <- read.csv(rccCollectorToolExport, header=FALSE, stringsAsFactors=FALSE)

    #
    # pdata
    #

    csv.fdata_linenums <- which(csv[, 1] %in% c("Endogenous", "Housekeeping"))
    if (length(csv.fdata_linenums) == 0)
    {
        stop("Expected but didn't find \"Endogenous\" or \"Housekeeping\" in the CodeClass column (the first column) of the input .CSV; check that the file contains the correct kind of data and is in the correct format")
    }
	pdata.m.t <- as.matrix( csv[1:(min(csv.fdata_linenums)-1), ] )
	rownames(pdata.m.t) <- pdata.m.t[, 1]

    # Strip the first col ("File Name" "Description" "Sample ID" ...) and the
    # next two cols (which should be blank in the pdata rows since they're
    # placeholders over the fdata GeneName and Accession cols so that the
    # pdata entries line up with the adata values).
    if ( !(all(pdata.m.t[, 2] == "") && all(pdata.m.t[, 3] == "")) ) {
        stop("Expected 2nd and 3rd cols of the .CSV's phenoData header (the first ~20 lines) to be empty")
    }
    pdata.m.t <- pdata.m.t[, -c(1,2,3)]

    # Check that the RLF field looks correct
    if ( !("Gene RLF" %in% rownames(pdata.m.t)) ) {
        stop("Expected but didn't find \"Gene RLF\" in the .CSV header section")
    }
    Gene_RLF.nonempty.unique <- unique(grep(value=TRUE, invert=TRUE, pattern="^ *$", x=pdata.m.t["Gene RLF", ]))
    if (length(Gene_RLF.nonempty.unique) == 0) {
        stop("Missing or empty values in \"Gene RLF\"")
    }
    if (length(Gene_RLF.nonempty.unique) > 1) {
        stop("The data seems to have been generated using more than one codeset (more than one unique value appears in \"Gene RLF\")")
    }

    # Remove header rows that have no info
    pdata_row_is_empty_throughout <- apply(pdata.m.t, 1, function(x) all(x == ""))
    if (any(pdata_row_is_empty_throughout)) {
        droppedRows <- rownames(pdata.m.t)[ pdata_row_is_empty_throughout ]
	    pdata.m.t <- pdata.m.t[ !pdata_row_is_empty_throughout, ]
		warning(
            paste0("The following rows in the .CSV's phenoData section (the first ~20 lines) were dropped since they have blank values throughout: \"",
                paste(collapse="\", \"", droppedRows), "\"")
        )
	}

    # Adjust "FOV Count" to "FovCount" and "FOV Counted" to "FovCounted"
    # (n.b. further error checking is in pdata_fdata_adata.to.rccSet())
    if ("FOV Count" %in% rownames(pdata.m.t)) {
        rownames(pdata.m.t)[ rownames(pdata.m.t) == "FOV Count" ] <- "FovCount"
    }
    if ("FOV Counted" %in% rownames(pdata.m.t)) {
        rownames(pdata.m.t)[ rownames(pdata.m.t) == "FOV Counted" ] <- "FovCounted"
    }

    # Remove any spaces in the labels and ensure that the first, which is the FileName column,
    # is actually labeled "FileName" (sometimes it may come in as "File name")
    rownames(pdata.m.t) <- gsub(" ", "", rownames(pdata.m.t))
    rownames(pdata.m.t)[1] <- "FileName"

    colnames(pdata.m.t) <- NULL

    #
    # fdata
    #

    fdata.df <- csv[min(csv.fdata_linenums):nrow(csv), 1:3]
    colnames(fdata.df) <- c("CodeClass", "GeneName", "Accession")

	#if (any(duplicated( fdata[, 2] )))
        #
        # GeneName duplication check: appeared in original readRccCollectorToolExport but no
        # longer needed since the feature names are now the concatenation of
        # CodeClass, GeneName, and Accession. -RZ 2015-01
        #

    #
    # adata
    #

    csv.adata <- csv[min(csv.fdata_linenums):nrow(csv), 4:ncol(csv)]
    adata.m <- data.matrix(csv.adata)
    dimnames(adata.m) <- NULL

    #
    # Final output
    #

    return( list(pdata.m=t(pdata.m.t), fdata.df=fdata.df, adata.m=adata.m) )
}

##' @title
##' pdata_fdata_adata.to.rccSet
##'
##' @description
##' Second stage of readRccBatch()/readRccCollectorToolExport() -- not intended for external use.
##'
##' @param  pdata_fdata_adata   List containing the pdata, fdata, and adata returned by
##'                             rccDir.to.pdata_fdata_adata() or nSolverCsv.to.pdata_fdata_adata().
##'
##' @return
##'
##' An \code{\linkS4class{RccSet}} whose contents reflect the input data.
##'
##' @details
##' 
##' Note that a pseudo-count of 1 is always added to all measurements, to enable
##' subsequent log transformation of the data in cases where zero-counts are
##' present.
##'
##' N.B. The function currently expects certain columns to be present in
##' \code{pdata_fdata_a_data$pdata.m}, and it converts these to numerics. These
##' expectations should be incorporated into the class definition, and
##' conversion should only take place with a warning. Future updates will
##' address this.
##' 
##' @author Robert Ziman
##'
pdata_fdata_adata.to.rccSet <- function(pdata_fdata_adata)
{
    pdata.m  <- pdata_fdata_adata$pdata.m
    fdata.df <- pdata_fdata_adata$fdata.df
    adata.m  <- pdata_fdata_adata$adata.m

    # phenoData

    rownames(pdata.m) <- pdata.m[,"FileName"]

    pdata.adf <- new("AnnotatedDataFrame", data=as.data.frame(pdata.m, stringsAsFactors = FALSE))

    pdata.adf$BindingDensity <- as.numeric(pdata.adf$BindingDensity)
    pdata.adf$LaneID <- as.numeric(pdata.adf$LaneID)
    pdata.adf$FovCount <- as.numeric(pdata.adf$FovCount)
    pdata.adf$FovCounted <- as.numeric(pdata.adf$FovCounted)
    pdata.adf$StagePosition <- as.numeric(pdata.adf$StagePosition)

    # featureData

    rownames(fdata.df) <- with(fdata.df, paste0(CodeClass, "_", GeneName, "_", Accession)) # N.B. the same string is used as the merge key for adding annotations (see addCodesetAnnotation())
    fdata.adf <- new("AnnotatedDataFrame", data=fdata.df)

    # assayData

    dimnames(adata.m) <- list(rownames(fdata.df), rownames(pdata.m))

    # Output RccSet

    rccSet <- RccSet( adata.m, phenoData = pdata.adf, featureData = fdata.adf )
    
    # Pseudocount

    exprs(rccSet) <- exprs(rccSet) + 1

    return(rccSet)
}

##' @title
##' Read a directory of .RCC files to produce an RccSet
##'
##' @description
##'
##' Reads the contents of all .RCC files from a given directory into a new
##' RccSet object. Note: this function is not intended for external use. For
##' that, see newRccSet().
##'
##' @param  rccDir  Directory containing the .RCC files.
##'
##' @return
##'
##' An RccSet object that has raw counts in assayData, probe information
##' in fData, and sample annotation in pData.
##'
##' @export
##'
##' @examples
##' rccDir <- system.file("extdata", "RCC", package="NanoStringQCPro")
##' rccSet <- readRccBatch(rccDir)
##'
##' @author Robert Ziman
##'
readRccBatch <- function(rccDir)
{
    pdata_fdata_adata <- rccDir.to.pdata_fdata_adata(rccDir)
    rccSet <- pdata_fdata_adata.to.rccSet(pdata_fdata_adata)
    return(rccSet)
}

##' @title
##' Read RCC Collector Tool Export
##'
##' @description
##' Reads the contents of a .CSV file generated from the RCC Collector Tool Export feature
##' of NanoString's nSolver Analysis software into a new NanoString ExpressionSet object.
##' (Note: this function is not intended for external use. For that, see newRccSet().)
##'
##' @details
##' See 'details' in the readRccBatch() help page.
##'
##' @param  file    Path to the NSolver .csv file to be read.
##'
##' @return
##' A NanoString ExpressionSet object that has count data in exprs, probe
##' information in fData and sample annotation in pData.
##'
##' @export
##'
##' @examples
##' csvPath <- system.file("extdata", "nSolver", "RCC_collector_tool_export.csv", package="NanoStringQCPro")
##' rccSet <- readRccCollectorToolExport(csvPath)
##'
##' @author Dorothee Nickles, Thomas Sandmann
##'
readRccCollectorToolExport <- function(file)
{
    #stringsAsFactors.prev <- getOption("stringsAsFactors")

    pdata_fdata_adata <- nSolverCsv.to.pdata_fdata_adata(file)
    rccSet <- pdata_fdata_adata.to.rccSet(pdata_fdata_adata)

    #options(stringsAsFactors = stringsAsFactors.prev)

    if(!validObject(rccSet)) {
        stop( paste0("validObject() returns FALSE on rccSet generated from NSolver csv file \"", file, "\"") )
    }

    return(rccSet)
}

##' @title
##' getSpikeInInput
##'
##' @description
##' Gets the RNA ``spike-in'' input levels for positive and negative control probes from the label in their GeneName.
##' Note that this is a helper function for readRlf() and elsewhere and is not intended for external use.
##'
##' @param CodeClass        Character vector with code classes for each probe.
##' @param GeneName         Character vector with gene names/symbols for each probe.
##' @param nonCtrlProbeVal  Value to assign as the spike-in input for the non-control probes.
##'
##' @return
##' A data frame with the input CodeClass and GeneName but where the latter has been split into
##' two columns: one showing the GeneName for each probe with spike-in input labels removed
##' -- and another with the spike-in input levels.
##'
##' @note
##' "concn" abbreviation taken from http://www.cas.org/content/cas-standard-abbreviations.
##'
##' @author Robert Ziman
##
getSpikeInInput <- function(CodeClass, GeneName, nonCtrlProbeVal=NA)
{
    stopifnot(!is.null(CodeClass))
    stopifnot(!is.null(GeneName))
    stopifnot(length(CodeClass) == length(GeneName))

    GeneName.split <- strsplit(GeneName, '\\(')
    GeneName.without_parentheses <- vapply(GeneName.split, function(x) {x[1]}, character(1))
    SpikeInInput.char <- sub('\\)$', '', vapply(GeneName.split, function(x) {x[2]}, character(1)))

    ctrls <- (CodeClass %in% c("Positive", "Negative"))

    GeneName.adjusted <- GeneName
    GeneName.adjusted[ ctrls ] <- GeneName.without_parentheses[ ctrls ]     # Only strip the labels from the control probes (in some cases there will be GeneNames with parentheses!)

    SpikeInInput <- rep(nonCtrlProbeVal, length(GeneName))
    SpikeInInput[ ctrls ] <- as.numeric(SpikeInInput.char[ ctrls ])

    return(
        data.frame(stringsAsFactors=FALSE,
            CodeClass   = CodeClass,
            GeneName    = GeneName.adjusted,
            SpikeInInput = SpikeInInput
        )
    )
}

##' @title
##' Read RLF file
##'
##' @description
##' Reads the contents of an .RLF file into a data frame. RNA ``spike-in'' concentrations recorded in the
##' GeneName for positive and negative control probes are stripped and stored in a separate
##' column in the output. An error will be generated for any recognized deviations from the
##' expected file format.
##'
##' @param  rlfPath     Path to the .RLF file
##'
##' @return
##' A data frame containing the contents of the .RLF file.
##'
##' @export
##'
##' @examples
##' rlfPath <- system.file("extdata", "RLF", "NQCP_example.rlf", package="NanoStringQCPro")
##' rlf <- readRlf(rlfPath)
##'
##' @author Robert Ziman
##'
readRlf <- function(rlfPath)
{
    if (!file.exists(rlfPath)) {
        stop(
            paste0("rlfPath \"", rlfPath, "\" not found")
        )
    }

	suppressWarnings(
		all.lines <- readLines(rlfPath)
                        #
                        # N.B. readLines() appears to automatically strip any trailing carriage returns (^M, 0xD)
                        # from files coming from Windows, etc.
                        #
	)

	i <- which(all.lines == "[Content]")

	if (!(    grepl("^ColumnCount=", all.lines[i+1])
           && grepl("^RecordCount=", all.lines[i+2])
           && grepl("^Columns=",     all.lines[i+3]) ))
	{
		stop("expected but did not find \"ColumnCount=\", \"RecordCount=\", and \"Columns=\" lines after \"[Content]\"")
	}

	nrec <- as.integer(sub("RecordCount=", "", all.lines[i+2]))

	content.lines.header <- all.lines[i+3]
	content.lines.header <- sub("Columns=", "", content.lines.header)

    i_content.lines.records_start <- i+4
    i_content.lines.records_end <- i+4+nrec-1
    if (i_content.lines.records_end != length(all.lines)) {
        stop("mismatch between RecordCount and the number of records actually in the file")
    }
	content.lines.records <- all.lines[ (i+4):(i+4+nrec-1) ]
	content.lines.records <- sub("Record[0-9]*=", "", content.lines.records)

	content.lines <- c(content.lines.header, content.lines.records)
	
	tc <- textConnection(content.lines)
	rlf <- read.csv(tc, header=TRUE, as.is=TRUE)

    classnames.lines <- all.lines[ grepl("^ClassName", all.lines) ]
    classnames.lines.strsplit <- strsplit(classnames.lines, "=")
    classnames        <- vapply(classnames.lines.strsplit, function (x) {x[2]}, character(1))
    names(classnames) <- vapply(classnames.lines.strsplit, function (x) {x[1]}, character(1))

    rlf$ClassNameID <- paste0("ClassName", rlf$Classification)
    rlf$ClassName <- as.character( classnames[rlf$ClassNameID] )
    rlf$Classification <- NULL
    rlf$ClassNameID <- NULL

    rlf <- rlf[!(rlf$ClassName == "Reserved"),]     # TODO: check if this will ever be an issue
    rlf <- rlf[!(rlf$ClassName == "Binding"),]
    rlf <- rlf[!(rlf$ClassName == "Purification"),]

    spikein <- getSpikeInInput(CodeClass=rlf$ClassName, GeneName=rlf$GeneName, nonCtrlProbeVal=0)
    rlf$GeneName <- spikein$GeneName
    rlf$SpikeInInput <- spikein$SpikeInInput

    ctrls <- (rlf$ClassName %in% c("Positive", "Negative"))
    if (any(is.na(rlf$SpikeInInput[ ctrls ]))) {
        stop("Missing or malformed spike-in labels in GeneName for some control probes")
    }

    if (any(is.na(rlf$Accession))) {
        stop("Missing accessions for some entries")
    }
    #if (any(duplicated(rlf$GeneName)))
        #
        # GeneName duplication check: appeared in original readRccCollectorToolExport but no
        # longer needed since the feature names are now the concatenation of
        # CodeClass, GeneName, and Accession. -RZ 2015-01
        #

    #colnames(rlf)[ colnames(rlf) == "Accession" ] <- "Accession_RLF"
    colnames(rlf)[ colnames(rlf) == "ClassName" ] <- "CodeClass"
    #colnames(rlf)[ colnames(rlf) == "Comments" ] <- "Comments_RLF"
    #colnames(rlf)[ colnames(rlf) == "GeneName" ] <- "GeneName_RLF"

    return(rlf)
}

##' @title          Read .csv containing CDR 'Design Data' extract
##' @description    Return a data frame containing the contents of the 'Design Data' tab
##'                 extracted from a CDR spreadsheet. The extract, a .csv file, must be
##'                 manually prepared in advance (see 'details' section in the
##'                 buildCodesetAnnotation() help page for more info).
##'
##' @param  cdrDesignData    Path to the .csv file containing the content extracted from the CDR's 'Design Data' tab
##'
##' @return
##' A data frame containing the contents of the CDR 'Design Data' tab.
##'
##' @export
##'
##' @examples
##' path <- system.file("extdata", "CDR", "CDR-DesignData.csv", package="NanoStringQCPro")
##' cdr <- readCdrDesignData(path)
##'
##' @author Robert Ziman
##'
readCdrDesignData <- function(cdrDesignData)
{
    if (!file.exists(cdrDesignData)) {
        stop(
            paste0("CDR Design Data csv \"", cdrDesignData, "\" not found")
        )
    }

    cdr <- read.csv(cdrDesignData, header=TRUE, as.is=TRUE)

    if (!"Customer.Identifier" %in% colnames(cdr)) {
        stop(
            paste0("'Customer Identifier' not found in CDR Design Data \"", cdrDesignData,
                "\"; check that the file has the right format (see details in buildCodesetAnnotation() help page)")
        )
    }

    if (which(colnames(cdr) == "Customer.Identifier") == 2) {
        cdr <- cdr[, -1]     # TODO: Confirm exactly what the column is that's being dropped here -RZ 2015-02-12
    }

    tmp <- which(colnames(cdr) == "Comments")
    cdr <- cdr[, 1:tmp]

    if (any(is.na(cdr$Accession))) {
        cdr <- cdr[!is.na(cdr$Accession), ]
        warning(
            paste0("Dropped some rows in CDR Design Data \"", cdrDesignData, "\" since their Accession field was empty")
        )
    }

    return(cdr)
}

##' @title          Build NanoString codeset annotation
##' @description    This function returns a data frame whose content is the combination of
##'                 the NanoString-provided codeset annotation (.RLF file and the
##'                 "Design Data" tab of the CDR spreadsheet) with gene annotation in the
##'                 org.Hs.eg.db package. 
##'
##' @param  rlfPath                 Path to the RLF file
##' @param  cdrDesignData    Path to a manually prepared csv export of the
##'                                 "Design Data" tab of the CDR file (optional;
##'                                 see 'details' sectio below for how the export
##'                                 should be prepared)
##' @param  removeRedundantCols     Logical. If TRUE, cols in the CDR that are redundant
##'                                 with those in the RLF will be omitted from the output.
##' @param  addEgAnnotations        Logical indicating whether or not to add
##'                                 EntrezGene IDs and HGNC symbols from the
##'                                 org.Hs.eg.db package.
##'
##' @return
##' A data frame whose content is the combination of the NanoString-provided
##' codeset annotation with gene annotation in the org.Hs.eg.db package.
##'
##' @details
##' The original NanoString provided .RLF file is expected as input. This
##' file is the master (i.e. only probes listed here will be annotated;
##' any extra ones in the CDR export will be dropped). If
##' the CDR "Design Data" .csv is specified, the function expects this .csv
##' file to be generated from the "Design Data" tab of the original NanoString
##' provided Excel CDR file. This tab needs to be trimmed by skipping the
##' NanoString header and first column containing only integers;
##' the resulting .csv should contain the actual table (including its header
##' -- beginning with "Customer Identifier"). The function will match and
##' join the .RLF and CDR .csv using their "ProbeID" and "NSID" fields, and
##' then it will add gene annotation (EntrezGene ID, HGNC symbol, and
##' chromosomal position) by doing lookups in the org.Hs.eg.db package using
##' the RefSeq accessions from the RLF.
##'
##' @import
##' org.Hs.eg.db
##' AnnotationDbi
##'
##' @export
##'
##' @examples
##'	rlfPath <- system.file("extdata", "RLF", "NQCP_example.rlf", package="NanoStringQCPro")
##'	cdrDesignData <- system.file("extdata", "CDR", "CDR-DesignData.csv", package="NanoStringQCPro")
##'	annot <- buildCodesetAnnotation(rlfPath, cdrDesignData)
##'
##' @author Dorothee Nickles, Robert Ziman
##'
buildCodesetAnnotation <- function(rlfPath=NULL,
                                   cdrDesignData=NULL,
                                   removeRedundantCols=TRUE,
                                   addEgAnnotations=FALSE)
{
    if (is.null(rlfPath)) {
        stop("RLF not specified")
    }
    rlf <- readRlf(rlfPath)

    if (!is.null(cdrDesignData))
    {
        cdr <- readCdrDesignData(cdrDesignData)

        misProbes <- summary(cdr$NSID %in% rlf$ProbeID)["FALSE"]
        if (!is.na(misProbes)) {
            warning(sprintf("%s probe(s) in the CDR have no entry in the RLF. These will be excluded from annotation.", misProbes))
        }
    
        annot <- merge(rlf, cdr, by.x="ProbeID", by.y="NSID", all.x=TRUE, suffixes=c("", "_CDR"))

        tmp <- is.na(annot$PN.CP.RP.)
        annot$PN.CP.RP.[tmp] <- annot$ProbeID[tmp]

        if (removeRedundantCols)
        {
            annot$Accession_CDR <- NULL     # redundant with Accession
            annot$Gene <- NULL              # redundant with GeneName
            if ("TargetSeq" %in% colnames(annot)) {
                annot$Target.Sequence <- NULL
            }
        }
    }
    else {
        annot <- rlf
    }

    if (addEgAnnotations == TRUE)
    {
        #
        # N.B. Before renaming any of the cols below, check for references to them.
        # Other parts of the code, QC template, or old import scripts may depend on
        # these columns having the precise names they have here. In turn, see
        # createExpressionplotAnnot() for renaming of some of the columns
        # to the precise names ExpressionPlot expects...
        # -RZ 2015-01-19
        #

        annot$RefSeqMrna <- vapply(as.character(annot$Accession), function(x) {unlist(strsplit(x, "[.]"))[1]}, character(1))

        #
        # EntrezGene ID
        #
        gene_ids <- AnnotationDbi::mget(
            annot$RefSeqMrna,
            org.Hs.egREFSEQ2EG,
            ifnotfound = NA_character_
        )
        if ( !all( vapply(gene_ids, length, integer(1)) == 1 ) ) {
            stop( "Unexpected length for some org.Hs.egREFSEQ2EG lookup results" )
        }
        annot$GeneID <- unlist( gene_ids )

        nonNA_geneID <- !is.na(annot$GeneID)

        #
        # HGNC symbol
        #
        hgnc_symbols <- AnnotationDbi::mget(
            annot$GeneID[ nonNA_geneID ],
            org.Hs.egSYMBOL,
            ifnotfound = NA_character_
        )
        if (!all( vapply(hgnc_symbols, length, integer(1)) == 1 )) {
            stop("Unexpected length for some org.Hs.egSYMBOL lookup results")
        }
        annot$Hgnc_Symbol[ nonNA_geneID ] <- unlist(hgnc_symbols)
    }

    #
    # Output
    #

    annot <- annot[ order(annot$Accession), ]

    #if (!is.null(annotPath)) {
    #    write.table(annot, file=annotPath, sep=",", row.names=FALSE)           # Removed as per decision
    #    message(paste0("Wrote codeset annotation to \"", annotPath, "\""))     # following the code review.
    #}                                                                          # -RZ 2015-03-25

    return(annot)
}

##' @title          Add NanoString codeset annotation to a NanoString ExpressionSet object
##' @description    Returns a copy of the input ExpressionSet where the codeset annotation
##'                 has been merged into its fData slot. The merge key for each is a string
##'                 formed from the concatenation of their CodeClass, GeneName, and
##'                 Accession columns ("<CodeClass>_<GeneName>_<Accession>"). For creating
##'                 the codeset annotation object, see buildCodesetAnnotation().
##'
##' @param  rccSet          NanoString ExpressionSet object.
##' @param  annot           Data frame containing the codeset annotation.
##' @param  reorder         Logical indicating whether the probes should be reordered
##'                         according to their barcodes (this can help in identifying
##'                         barcode-specific artifacts -- i.e. background noise).
##' @param  showWarnings    Logical indicating whether or not warnings should be shown, if any.
##'
##' @return
##' A copy of the input ExpressionSet where the codeset annotation has been
##' merged into its fData slot.
##'
##' @export
##'
##' @examples
##' rccDir <- system.file("extdata", "RCC", package="NanoStringQCPro")
##' rccSet <- readRccBatch(rccDir)
##' rlfPath <- system.file("extdata", "RLF", "NQCP_example.rlf", package="NanoStringQCPro")
##' annot <- buildCodesetAnnotation(rlfPath)
##' rccSet.annotated <- addCodesetAnnotation(rccSet, annot)
##'
##' @author Dorothee Nickles, Robert Ziman
##'
addCodesetAnnotation <- function(rccSet, annot, reorder=TRUE, showWarnings=TRUE)
{
    fdata <- fData(rccSet)

    CodeClass <- NULL   # prevent BiocCheck warning
    GeneName <- NULL
    Accession <- NULL
    fdata$mergekey <- with(fdata, paste0(CodeClass, "_", GeneName, "_", Accession))
    annot$mergekey <- with(annot, paste0(CodeClass, "_", GeneName, "_", Accession))

    annot$mergeflag <- rep(TRUE, nrow(annot))

    fdata.annot <- merge(
        fdata
        ,annot
        ,by.x="mergekey"
        ,by.y="mergekey"
        ,all.x=TRUE
        ,suffixes=c("", "_codesetAnnot")
	)

    fdata.annot$CodeClass_codesetAnnot <- NULL
    fdata.annot$Accession_codesetAnnot <- NULL
    fdata.annot$GeneName_codesetAnnot <- NULL

    if (showWarnings && any(grepl('_codesetAnnot$', colnames(fdata.annot))))
        {
            warning("The names of some columns in the codeset annotation were identical to those already in the rccSet's fData. To distinguish them from the fData columns in the input, the output's fData has these columns suffixed with \"_codesetAnnot\".")
        }

    # merge() messes up the order of the rows -- put them back as they were!
    fdata.annot <- fdata.annot[ match(fdata$mergekey, fdata.annot$mergekey), ]

    # Error out if any of the features didn't have a matching entry in annot. This should
    # hopefully take care of anything incomplete, mistaken, or unusual in the RLF.
    mergeflag_NA <- which(is.na(fdata.annot$mergeflag))
    if (length(mergeflag_NA) > 0)
        {
            stop(
                paste0("could not find annotations for some features ",
                       ifelse(length(mergeflag_NA) <= 3, "", "(only the first few shown here)"),
                       ": ", paste(collapse=", ", fdata.annot$mergekey[ head(mergeflag_NA, n=3) ]))
                )
        }
    fdata.annot$mergeflag <- NULL

    fdata.annot$mergekey <- NULL
    rownames(fdata.annot) <- rownames(fdata)

    rccSet.annot <- copyRccSet(rccSet) # Using copyRccSet() to be *sure* the original doesn't get affected in later code!
    fData(rccSet.annot) <- fdata.annot
    if (reorder == TRUE)
        {
	    rccSet.annot <- rccSet.annot[order(fData(rccSet.annot)$BarCode),]
	    rccSet.annot <- rccSet.annot[order(fData(rccSet.annot)$CodeClass),]
	}

    return(rccSet.annot)
}

##' @title
##' Create NanoString RccSet object
##'
##' @description
##' 
##' This is the main wrapper function for generating an ExpressionSet object
##' from NanoString data. The function takes as input either a directory
##' containing NanoString .RCC files (with the raw data) or a .CSV file
##' generated via the RCC Collector Tool Export feature of NanoString's nSolver
##' Analysis Software; a path to the .RLF file describing the codeset used;
##' optional paths to additional annotation about the features and samples; and
##' details about the experiment. It returns an RccSet object.
##'
##' @param  rccDir                  Directory containing the NanoString .RCC files
##'                                 with the raw count data. All files with the
##'                                 upper case .RCC extension will be used.
##' 
##' @param  rlfPath                 Path to the NanoString .RLF file describing
##'                                 the codeset used in generating the .RCCs.
##' 
##' @param  rccCollectorToolExport  Path to a .CSV file generated via the RCC Collector
##'                                 Tool Export feature of NanoString's nSolver
##'                                 Analysis Software. (Note that this is an alternative
##'                                 to rccDir, and if both arguments are specified at the
##'                                 same time, the function will throw an error.)
##' 
##' @param  cdrDesignData           Path to a .CSV extract of the "Design Data"
##'                                 tab of a CDR spreadsheet corresponding to the
##'                                 rest of the input files. See 'Details'
##'                                 section of the buildCodesetAnnotation() help page
##'                                 for more info on how this extract should be
##'                                 prepared.
##' 
##' @param  extraPdata              Vector of paths to files containing additional
##'                                 annotation about the samples which will be added
##'                                 to the phenoData of the output ExpressionSet. All
##'                                 files should be tab-separated and should contain
##'                                 a column labelled "FileName" whose values correspond
##'                                 exactly to the .RCC filenames in the
##'                                 rccDir or listed in the RCC Collector Tool
##'                                 Export. More than one such file may be
##'                                 used. A SampleType column should be present
##'                                 in at most one file.
##' 
##' @param  blankLabel              Value for the output's phenoData SampleType column
##'                                 that will indicate blank samples. This will be
##'                                 recorded in the varMetadata for
##'                                 SampleType. Blank samples, if available,
##'                                 play an important role in preprocessing.
##' 
##' @param  addEgAnnotations        Logical indicating whether or not to add
##'                                 EntrezGene annotations from the org.Hs.eg.db
##'                                 package.
##' 
##' @param  experimentData.name     String passed to the 'name' slot of the
##'                                 output ExpressionSet's experimentData.
##' 
##' @param  experimentData.lab      String passed to the 'lab' slot of the
##'                                 output ExpressionSet's experimentData.
##' 
##' @param  experimentData.contact  String passed to the 'contact' slot of the
##'                                 output ExpressionSet's experimentData.
##' 
##' @param  experimentData.title    String passed to the 'title' slot of the
##'                                 output ExpressionSet's experimentData.
##' 
##' @param  experimentData.abstract String passed to the 'abstract' slot of the
##'                                 output ExpressionSet's experimentData.
##' 
##' @param  experimentData.url      String passed to the 'url' slot of the
##'                                 output ExpressionSet's experimentData.
##' 
##' @param  experimentData.other    List passed to the 'other' slot of the
##'                                 output ExpressionSet's experimentData.
##' 
##' @param  validate                Logical. If TRUE, the output ExpressionSet
##'                                 will be checked via validRccSet(stopOnError=TRUE).
##'
##' @return
##' 
##' An \code{\linkS4class{RccSet}} containing the raw NanoString data and annotations.
##'
##' @details
##' 
##' In the .RLF (and sometimes in the .RCC files), the GeneName field for
##' positive and negative control probes contains a parenthesized label
##' indicating the RNA \dQuote{spike-in} levels for each probe. These labels are
##' removed from the control probe GeneNames in the output and recorded instead
##' in SpikeInInput in the output's featureData.
##'
##' A pseudocount of 1 is added to all measurements to enable subsequent
##' log transformation of the data.
##'
##' If the phenoData SampleType column is not specified via an annotation file
##' passed in through extraPdata, it will be created and assigned NA for all
##' samples.
##'
##' @import Biobase
##'
##' @export
##'
##' @examples
##' rccSet <- newRccSet(
##'      rccDir = system.file("extdata", "RCC", package="NanoStringQCPro")
##'     ,rlfPath = system.file("extdata", "RLF", "NQCP_example.rlf", package="NanoStringQCPro")
##'     ,extraPdata = system.file("extdata", "extraPdata", "SampleType.txt", package="NanoStringQCPro")
##'     ,blankLabel = "blank"
##'     ,experimentData.name = "Robert Ziman"
##'     ,experimentData.lab = "Richard Bourgon"
##'     ,experimentData.contact = "ziman.robert@@gene.com"
##'     ,experimentData.title = "NanoStringQCPro example dataset"
##'     ,experimentData.abstract = "Example data for the NanoStringQCPro package"
##' )
##'
##' @author Robert Ziman
##'
newRccSet <- function(rccDir,
                      rlfPath,
                      rccCollectorToolExport    = NULL,
                      cdrDesignData             = NULL,
                      extraPdata                = NULL,
                      blankLabel                = "blank",
                      addEgAnnotations          = FALSE,
                      experimentData.name       = "",
                      experimentData.lab        = "",
                      experimentData.contact    = "",
                      experimentData.title      = "",
                      experimentData.abstract   = "",
                      experimentData.url        = "",
                      experimentData.other      = list(),
                      validate                  = TRUE)
{
    if (!missing(rccDir) && !is.null(rccDir))
        {
            if(!missing(rccCollectorToolExport) && !is.null(rccCollectorToolExport)) {
                stop("Both rccDir and rccCollectorToolExport are specified")
            }

            if (grepl('\\.csv$', tolower(rccDir))) {
                stop("A .CSV file appears to have been specified as the rccDir argument (if this is the RCC Collector Tool Export, specify its path in the rccCollectorToolExport argument instead)")
            }

            message("Reading RCC files...")
            rccSet <- readRccBatch(rccDir)

        } else if(!missing(rccCollectorToolExport) && !is.null(rccCollectorToolExport)) {

            message("Reading RCC Collector Tool Export...")
            rccSet <- readRccCollectorToolExport(rccCollectorToolExport)

        } else {
            stop("Either rccDir or rccCollectorToolExport must be specified")
        }

    experimentData(rccSet) <- new(
        "MIAME"
        ,name       = experimentData.name
        ,lab        = experimentData.lab
        ,contact    = experimentData.contact
        ,title      = experimentData.title
        ,abstract   = experimentData.abstract
        ,url        = experimentData.url
        ,other      = experimentData.other
        )

    preproc(rccSet) <- list(state="newRccSet")

    GeneRLF <- unique(as.character(pData(rccSet)$GeneRLF))
    if (length(GeneRLF) != 1) {
        stop("Expecting exactly one RLF file to be referenced by all RCC files; check RCCs or RCC Collector Tool Export")
    }
    rlfFile <- basename(rlfPath)
    if (grepl('\\.rlf$', rlfFile)) {
        rlfFile.noext <- sub('\\.rlf$', '', rlfFile)
    } else if (grepl('\\.RLF$', rlfFile)) {
        rlfFile.noext <- sub('\\.RLF$', '', rlfFile)
    }
    if (GeneRLF != rlfFile.noext) {
        stop(paste0("Specified RLF (\"", rlfPath, "\") doesn't match that in the input RCCs or RCC Collector Tool Export (\"", GeneRLF, "\")"))
    }
    annotName <- GeneRLF
    annotation(rccSet) <- GeneRLF
    
    pData(rccSet)$GeneRLF <- NULL # Removing it from pData since it's now in annotation(rccSet)

    if ((!is.null(extraPdata)) && (length(extraPdata) > 0))
        {
            new_pdata <- pData(rccSet)

            for(i in 1:length(extraPdata))
                {
                    extraPdataPath <- extraPdata[i]

                    if (!file.exists(extraPdataPath))
                        stop( paste0("extraPdata file does not exist: \"", extraPdataPath, "\"") )

                    xpd <- read.table(file=extraPdataPath, header=TRUE, sep="\t", stringsAsFactors=FALSE)

                    colnames(xpd)[ tolower(colnames(xpd)) == "filename" ] <- "FileName"
                    colnames(xpd)[ tolower(colnames(xpd)) == "sampletype" ] <- "SampleType"

                    if (!("FileName" %in% colnames(xpd)))
                        stop( "extraPdata file \"", extraPdataPath, "\" is missing the \"FileName\" column" )

                    missing_samples <- setdiff( pData(rccSet)$FileName, xpd$FileName )
                    if ( length( missing_samples ) > 0 )
                        stop( "Some samples were not found in extraPdata file \"", extraPdataPath, "\": ", paste( missing_samples, collapse = ", " ) )

                    if (any(duplicated(xpd$FileName)))
                        stop( "Duplicate entries found (check FileName column) in extraPdata file \"", extraPdataPath, "\"" )

                    if (nrow(xpd) > ncol(rccSet))
                        warning( "Ignoring extraneous entries found in extraPdata file \"", extraPdataPath, "\"" )

                    if ( "SampleType" %in% colnames( xpd ) ) {
                        if ( !all( is.na( new_pdata$SampleType ) ) )
                            stop( "SampleType column being supplied by more than one extraPdata file" )
                        new_pdata$SampleType <- xpd$SampleType
                        xpd$SampleType <- NULL
                    }

                    new_pdata <- merge(
                        new_pdata
                        ,xpd
                        ,by="FileName"
                        ,all.x=TRUE
                        ,all.y=FALSE
                        ,sort=FALSE
                        ,suffixes=c("", paste0("_", basename(extraPdata[i])))
                        )

                }

            # Ensure that the order of new_pdata exactly matches that in the
            # original (even with sort=FALSE, merge() doesn't necessarily ensure that).
            order_in_orig_pdata <- match(new_pdata$FileName, pData(rccSet)$FileName)
            new_pdata <- new_pdata[order_in_orig_pdata, ]
            stopifnot(identical(new_pdata$FileName, pData(rccSet)$FileName))

            # Restore the rownames (merge() mucks that up too).
            rownames(new_pdata) <- rownames(pData(rccSet))

            # Now, finally, the original pData can be updated.
            pData(rccSet) <- new_pdata
        }

    varMetadata( phenoData( rccSet ) )[ "SampleType", "labelDescription" ] <- paste0("blankLabel='", blankLabel, "'")

    # The codeset annotation should be rebuilt from input RLF/CDR/etc *every* time; the csv output to annotPath will thus
    # be just for reference. -RZ 2014-12-14
    annot <- buildCodesetAnnotation(
        rlfPath                 = rlfPath,
        cdrDesignData           = cdrDesignData,
        addEgAnnotations        = addEgAnnotations,
        removeRedundantCols     = TRUE
        )
    rccSet <- addCodesetAnnotation(
        rccSet          = rccSet,
        annot           = annot,
        reorder         = TRUE,
        showWarnings    = FALSE
        )

    #
    # Record the org.Hs.eg.db version
    #
    if (addEgAnnotations == TRUE)
        {
            if ("org.Hs.eg.db" %in% names(sessionInfo()$otherPkgs)) {
                org.Hs.eg.db_version <- sessionInfo()$otherPkgs$org.Hs.eg.db$Version

            } else if ("org.Hs.eg.db" %in% names(sessionInfo()$loadedOnly)) {
                org.Hs.eg.db_version <- sessionInfo()$loadedOnly$org.Hs.eg.db$Version

            } else {
                warning("Attempted to record org.Hs.eg.db version but did not find it in the expected locations in sessionInfo()")

            }

            preproc(rccSet) <- c(preproc(rccSet), fData=paste0("Added EntrezGene annotations using org.Hs.eg.db version ", org.Hs.eg.db_version))
        }

    # Remove superfluous columns; see also additional code in createOrUploadEpProject() (in NanoStringQCProGNE)
    pData(rccSet)$FileVersion         <- NULL
    pData(rccSet)$SoftwareVersion     <- NULL
    pData(rccSet)$Owner               <- NULL
    pData(rccSet)$SystemAPF           <- NULL
    #pData(rccSet)$FovCount                     # Required in the QC report
    #pData(rccSet)$FovCounted                   # Required in the QC report
    pData(rccSet)$ScannerID           <- NULL
    #pData(rccSet)$StagePosition                # Required in the QC report
    #pData(rccSet)$BindingDensity               # Required in the QC report
    #pData(rccSet)$CartridgeID                  # See comment by corresponding line in createOrUpdateEpProject() (in NanoStringQCProGNE)
    pData(rccSet)$CartridgeBarcode    <- NULL
    fData(rccSet)$CodeClass_codesetAnnot <- NULL
    fData(rccSet)$Accession_codesetAnnot <- NULL
    fData(rccSet)$GeneName_codesetAnnot  <- NULL
    fData(rccSet)$Accession_CDR       <- NULL

    # Initialize preproc
    preproc(rccSet) <- c(preproc(rccSet), exprs="Raw data")

    if (validate)
        validRccSet(rccSet, stopOnError=TRUE, reportWarnings=TRUE, showMessages=TRUE)
    
    return(rccSet)

}

##' @title          Positive control normalization
##' @description    Applies NanoString's recommended positive control normalization to data in a NanoString ExpressionSet object.
##'
##' @param  rccSet          NanoString ExpressionSet object
##' @param  metric          Metric of positive controls used for normalization; one of "mean", "median", or "sum" (the default)
##' @param  recordPosFactor Logical. If TRUE, 'PosFactor' will be added to the output's phenoData to record the positive control
##'                         scaling factor that was computed for each sample.
##'
##' @return
##' A NanoString ExpressionSet object that has count data adjusted by positive control counts.
##' The positive control scaling factor is recorded in PosFactor in the output's phenoData
##' (an error is generated if this column already exists in the input).
##'
##' @export
##'
##' @examples
##' data(example_rccSet)
##' pcnorm_example_rccSet <- posCtrlNorm(example_rccSet)
##'
##' @author Dorothee Nickles
##'
posCtrlNorm <- function(rccSet,
                                 metric=c("sum", "mean", "median"),
                                 recordPosFactor=TRUE)
{
	stopifnot( is(rccSet, "ExpressionSet") )

	metric <- match.arg(metric)

    #if (preproc(rccSet)$state != "newRccSet") {
    #    warning("Input already contains positive control normalized data")
    #}

	prccSet <- copyRccSet(rccSet)    # Using copyRccSet() to be *sure* the original doesn't get affected in later code!

	fun <- match.fun(metric)
	posSignal <- apply(exprs(prccSet)[fData(prccSet)$CodeClass == "Positive", ], 2, fun)
	signalMean <- mean(posSignal)
	posFactor <- signalMean/posSignal

    assayData(prccSet)$posCtrlData <- t(apply(assayData(prccSet)$exprs, 1, function(x) { x * posFactor }))
    preproc(prccSet) <- c(preproc(prccSet), assayData_posCtrlData="Positive control normalized data")
    preproc(prccSet)$state <- "posCtrlNorm"

    if (recordPosFactor) {
        #if ("PosFactor" %in% colnames(pData(prccSet))) {
        #    stop("PosFactor already exists in the input's phenoData")
        #}
        pData(prccSet)$PosFactor <- posFactor
        PosFactor.metadata_rownum <- which(rownames(varMetadata(phenoData(prccSet))) == "PosFactor")
        varMetadata(phenoData(prccSet))$labelDescription[ PosFactor.metadata_rownum ] <- "Positive control scaling factor"
    }

	return(prccSet)
}

##' @title
##' NanoString count data normalization
##'
##' @description
##' Function to normalize count data given a NanoString ExpressionSet object.
##' Note that count data is log2-transformed before normalization and *remains*
##' log2-transformed returned ExpressionSet. If housekeeping normalization is
##' specified, the housekeeping features must be specified as well (see below)
##' and this will be recorded in the returned ExpressionSet in a new fData column
##' named "is.housekeeping".
##'
##' @param  rccSet          NanoString ExpressionSet object
##' @param  method          Normalization method (one of "median", "mean", or
##'                         "housekeeping")
##' @param  hk              Logical (TRUE/FALSE) vector defining, for each
##'                         feature, whether or not it shall be used for
##'                         housekeeping normalization if method="housekeeping"
##'
##' @return
##' A NanoString ExpressionSet object that has log2-transformed and
##' normalized count data.
##'
##' @export
##'
##' @examples
##' data(example_rccSet)
##'
##' pcnorm_example_rccSet <- posCtrlNorm(example_rccSet)
##' bg <- getBackground(pcnorm_example_rccSet)
##' bgcorr_example_rccSet <- subtractBackground(pcnorm_example_rccSet, bg)
##'
##' gmnorm_example_rccSet <- contentNorm(bgcorr_example_rccSet, method="median")
##' hknorm_example_rccSet <- contentNorm(bgcorr_example_rccSet, method="housekeeping",
##'     hk=(tolower(fData(bgcorr_example_rccSet)$CodeClass) == "housekeeping"))
##'
##' @author Dorothee Nickles
##'
contentNorm <- function(rccSet,
                        method=c("median", "mean", "housekeeping"),
                        hk=NULL)
{
	stopifnot( is(rccSet, "ExpressionSet") )

	method <- match.arg(method)

	nrccSet <- copyRccSet(rccSet)    # Using copyRccSet() to be *sure* the original doesn't get affected in later code!

    state <- preproc(rccSet)$state
    if (state == "newRccSet") {
        stop("Positive control normalization and background correction should be applied before calling this function")
    } else if (state == "posCtrlNorm") {
        stop("Background correction should be applied before calling this function")
    } else if (state == "subtractBackground") {
        M <- assayData(nrccSet)$bgCorrData
    } else if (state == "preprocRccSet") {
        #warning("Input already contains content normalized data")
        M <- assayData(nrccSet)$bgCorrData
    }

    # Apply the log2 transformation. Note that all the code below is thus
    # operating on values that are on a log2 scale.
    M <- log2(M)

	if (method == "median") {

        message("The data will be normalized by the median of all features.")
        normData_preprocList_value <- "Preprocessed and median-normalized data (log2 scale)"

        nFactors <- apply( M[ fData(nrccSet)$CodeClass == "Endogenous", ], 2, median )
        nFact <-  nFactors - mean(nFactors)

	} else if (method == "mean") {

        message("The data will be normalized by the mean of all features.")
        normData_preprocList_value <- "Preprocessed and mean-normalized data (log2 scale)"

        nFactors <- apply( M[ fData(nrccSet)$CodeClass == "Endogenous", ], 2, mean )
        nFact <-  nFactors - mean(nFactors)

    } else { # method == "housekeeping"

        if (is.null(hk)) {
            hk <- (fData(nrccSet)$CodeClass == "Housekeeping")
            if (all(hk == FALSE)) {
                stop("no housekeeping features defined ('hk' arg is missing or NULL and no features have CodeClass == \"Housekeeping\")")
            }
        }
        if (!is.logical(hk)) {
            stop("'hk' arg must be a logical (TRUE/FALSE) vector")
        }
        if (all(hk == FALSE)) {
            stop("no housekeeping features defined (all entries in 'hk' are FALSE)")
        }
        if (sum(hk) < 3) {
            warning("Less than three houskeeping features are defined")
        }

        if (sum(hk) > 1) {
            message(
                sprintf("The data will be normalized by the median of the following housekeeping features:\n%s",
                    paste(collapse="\n", rownames(fData(nrccSet))[hk]))
            )
            nFact <- apply( M[hk, ], 2, median )
        } else {
            nFact <- M[hk, ]
        }

        fData(nrccSet)$is.housekeeping <- hk

        normData_preprocList_value <- "Preprocessed and housekeeping-normalized data (log2 scale)"
    }

    M <- sweep(M, 2, nFact, "-")      # N.B. "-" is valid here since the opration is being performed
                                      # with log2-based values, so the matrix is effectively being
                                      # divided by nFact.

    assayData(nrccSet)$normData <- M
    preproc(nrccSet) <- c(preproc(nrccSet), assayData_normData=normData_preprocList_value)
    preproc(nrccSet)$state <- "preprocRccSet"

	return(nrccSet)
}

##' @title
##' Add background correction and normalization to a NanoString ExpressionSet
##'
##' @description
##' This function performs the positive control normalization, background
##' correction, and content normalization steps recommended for NanoString
##' datasets. For each step, a matrix is added to the assayData of the
##' resulting ExpressionSet object, and a string is appended to the
##' experimentData@@preprocessing list (accessible through preproc(rccSet)
##' where rccSet is an ExpressionSet output by this function). Positive
##' control scaling factors are recorded in the output's phenoData in a
##' column named 'PosCtrl', and the matrix used for background subtraction
##' is stored in the assayData as 'bgEstimates'. If housekeeping normalization
##' is performed, a column labeled 'is.housekeeping' is added to the fData that
##' indicates which features were used for it.
##'
##' @details
##' For more information on the rationale behind the recommended
##' preprocessing and normalization steps, see the vignette.
##'
##' @param  rccSet          NanoString ExpressionSet object to be used as input.
##' @param  method     Normalization method to use: "median" or "housekeeping".
##' @param  hkgenes    Optional character vector with gene symbols
##'                         to be used for normalization if method="housekeeping"
##'                         instead of the panel housekeeping features
##'                         defined in the input (i.e. those features with
##'                         featureData CodeClass == "Housekeeping"). If specified,
##'                         all features that match any of the specified symbols will
##'                         be used. (To specify specific features, use the hkfeatures
##'                         argument instead; see below.)
##' @param  hkfeatures Optional character vector with full feature names
##'                         ("<CodeClass>_<GeneName>_<Accession>", e.g.
##'                         "Endogenous_ACTG1_NM_001614.1") to be used for
##'                         normalization if method="housekeeping" instead of the
##'                         panel housekeeping features defined in the input.
##'                         (Note: if this argument is specified at the same time as
##'                         hkgenes, an error will be thrown.)
##'
##' @return
##' A copy of the input ExpressionSet with additional matrices in the
##' assayData for each successive preprocessing step (positive control normalization,
##' background estimation, background correction, and content normalization).
##' The final, fully preprocessed data is held in assayData(rccSet)$normData.
##' Details for each step are stored in the experimentData@@preprocessing
##' list (accessible through preproc(rccSet)).
##'
##' @export
##'
##' @examples
##' data(example_rccSet)
##' hknorm_example_rccSet <- preprocRccSet(rccSet = example_rccSet,
##'     method = "housekeeping")
##'
##' @references
##' \href{http://www.nanostring.com/media/pdf/MAN_nCounter_Gene_Expression_Data_Analysis_Guidelines.pdf}{NanoString nCounter(R) Expression Data Analysis Guide (2012)}
##'
##' @author Dorothee Nickles, Robert Ziman
##'
preprocRccSet <- function(rccSet,
                                     method        = c("median", "housekeeping"),
                                     hkgenes       = NULL,
                                     hkfeatures    = NULL)
{
    # Input validation
    method <- match.arg(method)
    hkgenes     <- hkgenes
    hkfeatures  <- hkfeatures
    if ( (!is.null(hkgenes) || !is.null(hkfeatures)) && (method != "housekeeping") ) {
        stop("hkgenes or hkfeatures is specified but method != \"housekeeping\"")
    }

    state <- preproc(rccSet)$state
    #if (state != "newRccSet") {
    #    warning("Input already contains preprocessed (or partially preprocessed) data")
    #}
    
    # Postive control normalization
    prccSet <- posCtrlNorm(rccSet)

    # Background correction
    if (hasBlanks(rccSet))
    {
        bgEstimates <- getBackground(rccSet = prccSet,
                                        reference = "both",
                                        stringency = 1)
        srccSet <- subtractBackground(rccSet = prccSet,
                                          bgEstimates = bgEstimates,
                                          description = "Background estimates using both blank samples and negative control probes")
    } else {

        bgEstimates <- getBackground(rccSet = prccSet,
                                        reference = "Negative",
                                        stringency = 0)
        srccSet <- subtractBackground(rccSet = prccSet,
                                          bgEstimates = bgEstimates,
                                          description <- "Background estimates using only negative control probes")
    }

    # Content normalization

    if (method == "housekeeping")
        #
        # Note that some of the code in the QC report template may be depend
        # on the exact strings (e.g. "Median", "Housekeeping") in the
        # preproc list, so be careful when changing them below. -RZ 2015-01-18
        #
    {
        if (length(hkgenes) > 0)
        {
            fdataGenes <- fData(rccSet)$GeneName

            hkgenesNotFound <- setdiff(hkgenes, fdataGenes)
            if (length(hkgenesNotFound) > 0) {
                stop(
                    paste0("the following genes specified in hkgenes were not found in the input's featureData: ",
                        paste(collapse=" ", hkgenesNotFound))
                )
            }

            hk <- rep(FALSE, nrow(fData(rccSet)))
            GeneName.match <- match(hkgenes, fdataGenes)
            hk[ GeneName.match ] <- TRUE

        } else if (length(hkfeatures) > 0) {

            fdataRownames <- rownames(fData(rccSet))

            hkfeaturesNotFound <- setdiff(hkfeatures, fdataRownames)
            if (length(hkfeaturesNotFound) > 0) {
                stop(
                    paste0("the following features specified in hkfeatures were not found in the input's featureData: \"",
                        paste(collapse="\", \"", hkfeaturesNotFound), "\"")
                )
            }

            hk <- rep(FALSE, nrow(fData(rccSet)))
            rowname.match <- match(hkfeatures, rownames(fData(rccSet)))
            hk[ rowname.match ] <- TRUE

        } else {
            hk <- fData(rccSet)$CodeClass == "Housekeeping"
        }

        if (all(hk == FALSE)) {
            stop("Housekeeping normalization selected but no housekeeping features are defined")
        }

        nrccSet <- contentNorm(srccSet, method="housekeeping", hk=hk)

    } else {    # method == "median"

        nrccSet <- contentNorm(srccSet, method="median")
    }

    return(nrccSet)
}

##' @title          Add sample QC flags to an rccSet
##' @description    Returns a copy of the input ExpressionSet with columns added to pData
##'                 from the provided sample QC flag annotation file. (That file is
##'                 produced by makeQCReport(); see its help page for more
##'                 details.)
##'
##' @param  rccSet      NanoString ExpressionSet object to be used as input
##' @param  flagFile    Path to a sample QC flag file as generated by the
##'                     NanoStringQCPro QC report (see makeQCReport())
##'
##' @return
##' A copy of the input ExpressionSet with columns added to pData from the QC
##' flag file.
##'
##' @author Dorothee Nickles
##'
addQCFlags <- function(rccSet, flagFile)
{
    stopifnot( is(rccSet, "ExpressionSet") )
    if (!file.exists(flagFile)) {
        stop("flagFile not found")
    }

    flags <- read.table(flagFile, header=TRUE, as.is=TRUE, sep="\t")
    stopifnot(all(flags$SampleIdentifier %in% pData(rccSet)$FileName))

	rccSet.withflags <- copyRccSet(rccSet)    # Using copyRccSet() to be *sure* the original doesn't get affected in later code!

    has.QCflag <-
        apply(
            flags,
            1,
            function(x) {
                return(ifelse(length(grep("TRUE", x[c("TechnicalFlags", "ControlFlags", "CountFlags")])) != 0, TRUE, FALSE))
            }
        )
    flags$has.QCflag <- has.QCflag

    mt <- match(pData(rccSet.withflags)$FileName, flags$SampleIdentifier)
    pData(rccSet.withflags) <- cbind(pData(rccSet.withflags),
        flags[mt,c("TechnicalFlags", "ControlFlags", "CountFlags", "has.QCflag")])

    return(rccSet.withflags)
}

