
## Data import and normalization functions -- and some small convenience functions

##' @title          Read an .RCC file
##' @description    Parse an .RCC file into a list with each part of the file (Header,
##'                 Sample_Attributes, Lane_Attributes, Code_Summary, etc) stored as a vector
##'                 or data frame.
##'
##' @param  rcc                 Path to the .RCC file.
##' @param  removeSpikeInLabels Logical. If TRUE (the default), RNA \dQuote{spike-in} input labels (if any)
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
##' rcc <- system.file("extdata", "RCC", "20140604_C1-unstim_C1-unstim_01.RCC", package="NanoStringQCPro")
##' rcc.ls <- readRcc(rcc)
##'
##' @author Robert Ziman
##'
readRcc <- function(rcc, removeSpikeInLabels=TRUE)
{
	suppressWarnings(
		lines <- readLines(rcc)
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

	rcc.ls <- list()

	rcc.ls$File_path <- rcc

	for (attr_tag in c("Header", "Sample_Attributes", "Lane_Attributes")) {
		attr_lnum <- which(lines.df$tag == attr_tag)
		attr.strsplit <- strsplit(lines.df$line[attr_lnum], ",")
		attr.v <- vapply(attr.strsplit, function(x) {if (!is.na(x[2])) {x[2]} else {""}}, character(1))
		names(attr.v) <- vapply(attr.strsplit, function(x) {x[1]}, character(1))
		rcc.ls[[attr_tag]] <- attr.v
	}

	Code_Summary_lnum.all <- which(lines.df$tag == "Code_Summary")
	Code_Summary_lnum.header <- Code_Summary_lnum.all[1]
	Code_Summary_lnum.body <- Code_Summary_lnum.all[2:length(Code_Summary_lnum.all)]

	if (lines.df$line[Code_Summary_lnum.header] != "CodeClass,Name,Accession,Count") {
	    stop(
	        paste0("Unrecognized header line for Code_Summary section of \"", rcc.ls, "\"",
	            " (expected \"CodeClass,Name,Accession,Count\" but got \"", lines.df$line[Code_Summary_lnum.header], "\")")
	    )
	}

	Code_Summary_fields <- strsplit(lines.df$line[Code_Summary_lnum.body], ",")
	rcc.ls$Code_Summary <- data.frame(stringsAsFactors=FALSE
		,CodeClass   = vapply(Code_Summary_fields, function(x) {x[1]}, character(1))
		,Name        = vapply(Code_Summary_fields, function(x) {x[2]}, character(1))
		,Accession   = vapply(Code_Summary_fields, function(x) {x[3]}, character(1))
		,Count       = vapply(Code_Summary_fields, function(x) {x[4]}, character(1))
	)

    if(removeSpikeInLabels)
    {
        spikein <- getSpikeInInput(CodeClass=rcc.ls$Code_Summary$CodeClass, GeneName=rcc.ls$Code_Summary$Name)
        rcc.ls$Code_Summary$Name <- spikein$GeneName
    }

	return(rcc.ls)
}

##' @title          rccFiles.to.pdata_fdata_adata
##' @description    First stage of readRccBatch(): produces a list containing matrices (for pdata and adata) and a data frame (for fdata) that
##'                 pdata_fdata_adata.to.rccSet() then transforms into a full RccSet (after some further checks and adjustments).
##'                 See also nSolverCsv.to.pdata_fdata_adata().
##'
##' @param  rccFiles  Vector of .RCC paths
##'
##' @return
##' A list containing matrices (for pdata and adata) and a data frame (for fdata)
##' that pdata_fdata_adata.to.rccSet() then transforms into a full ExpessionSet.
##'
##' @author Robert Ziman
##'
rccFiles.to.pdata_fdata_adata <- function(rccFiles)
{
    if (missing(rccFiles) || is.null(rccFiles) || length(rccFiles) == 0)
        stop("rccFiles argument is missing, NULL, or empty")

	rcc_1 <- readRcc( rccFiles[1] )

    filename <- basename( rccFiles[1] )
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
		rcc <- readRcc( rccFiles[i] )

        filename <- basename( rccFiles[i] )
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
				paste0("The attribute labels in the Header, Sample_Attributes, and Lane_Attributes sections of \"", filename,
					"\" don't exactly match those in the other RCC files parsed so far")
			)
		}
		if (!identical(rcc.fdata, rcc_1.fdata)) {
			stop(
				paste0("Values in the key cols (i.e. CodeClass, Name, and Accession) in the Code_Summary section of \"", filename,
					"\" don't exactly match those in the other RCC files parsed so far")
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
##'                 pdata_fdata_adata.to.rccSet then transforms into a full RccSet
##'                 (after some further checks and adjustments). Not intended for
##'                 external use; see also rccFiles.to.pdata_fdata_adata().
##'
##' @param  rccCollectorToolExport  Path to the nSolver RCC Collector Tool .CSV export.
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
##'                             rccFiles.to.pdata_fdata_adata() or nSolverCsv.to.pdata_fdata_adata().
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
##' Read RCC files
##'
##' @description
##'
##' Reads the contents of all .RCC files from a given directory into a new
##' RccSet object. Note: this function is not intended for external use. For
##' that, see newRccSet().
##'
##' @param  rccFiles  Vector of .RCC file paths
##'
##' @return
##'
##' An RccSet object that has raw counts in assayData, probe information
##' in fData, and sample annotation in pData.
##'
##' @author Robert Ziman
##'
readRccBatch <- function(rccFiles)
{
    pdata_fdata_adata <- rccFiles.to.pdata_fdata_adata(rccFiles)
    rccSet <- pdata_fdata_adata.to.rccSet(pdata_fdata_adata)
    return(rccSet)
}

##' @title
##' Read RCC Collector Tool Export
##'
##' @description
##' Reads the contents of a .CSV file generated from the RCC Collector Tool Export feature
##' of NanoString's nSolver Analysis software into a new RccSet object.
##' (Note: this function is not intended for external use. For that, see newRccSet().)
##'
##' @details
##' See 'details' in the readRccBatch() help page.
##'
##' @param  file    Path to the NSolver .CSV file to be read.
##'
##' @return
##' An RccSet object that has count data in exprs, probe
##' information in fData and sample annotation in pData.
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
        stop( paste0("validObject() returns FALSE on rccSet generated from NSolver .CSV file \"", file, "\"") )
    }

    return(rccSet)
}

##' @title
##' getSpikeInInput
##'
##' @description
##' Gets the RNA \dQuote{spike-in} input levels for positive and negative control probes from the label in their GeneName.
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
##' Reads the contents of an .RLF file into a data frame. RNA \dQuote{spike-in} concentrations recorded in the
##' GeneName for positive and negative control probes are stripped and stored in a separate
##' column in the output. An error will be generated for any recognized deviations from the
##' expected file format.
##'
##' @param  rlf     Path to the .RLF file
##'
##' @return
##' A data frame containing the contents of the .RLF file.
##'
##' @export
##'
##' @examples
##' rlf <- system.file("extdata", "RLF", "NQCP_example.rlf", package="NanoStringQCPro")
##' rlf.df <- readRlf(rlf)
##'
##' @author Robert Ziman
##'
readRlf <- function(rlf)
{
    if (!file.exists(rlf)) {
        stop(sprintf("rlf \"%s\" not found", rlf))
    }

	suppressWarnings(
		all.lines <- readLines(rlf)
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
	rlf.df <- read.csv(tc, header=TRUE, as.is=TRUE)

    classnames.lines <- all.lines[ grepl("^ClassName", all.lines) ]
    classnames.lines.strsplit <- strsplit(classnames.lines, "=")
    classnames        <- vapply(classnames.lines.strsplit, function (x) {x[2]}, character(1))
    names(classnames) <- vapply(classnames.lines.strsplit, function (x) {x[1]}, character(1))

    rlf.df$ClassNameID <- paste0("ClassName", rlf.df$Classification)
    rlf.df$ClassName <- as.character( classnames[rlf.df$ClassNameID] )
    rlf.df$Classification <- NULL
    rlf.df$ClassNameID <- NULL

    #rlf.df <- rlf.df[ !(rlf.df$ClassName == "Reserved"), ]     # TODO: check if this will ever be an issue
    #rlf.df <- rlf.df[ !(rlf.df$ClassName == "Binding"), ]
    #rlf.df <- rlf.df[ !(rlf.df$ClassName == "Purification"), ]
    rlf.df <- rlf.df[ !(rlf.df$ClassName %in% c("Reserved", "Binding", "Purification")), ]

    spikein <- getSpikeInInput(CodeClass=rlf.df$ClassName, GeneName=rlf.df$GeneName, nonCtrlProbeVal=0)
    rlf.df$GeneName <- spikein$GeneName
    rlf.df$SpikeInInput <- spikein$SpikeInInput

    ctrls <- (rlf.df$ClassName %in% c("Positive", "Negative"))
    if (any(is.na(rlf.df$SpikeInInput[ ctrls ]))) {
        stop("Missing or malformed spike-in labels in GeneName for some control probes")
    }

    if (any(is.na(rlf.df$Accession))) {
        stop("Missing accessions for some entries")
    }

    colnames(rlf.df)[ colnames(rlf.df) == "ClassName" ] <- "CodeClass"

    return(rlf.df)
}

##' @title          Read .CSV containing CDR 'Design Data' extract
##' @description    Return a data frame containing the contents of the 'Design Data' tab
##'                 extracted from a CDR spreadsheet. The extract, a .CSV file, must be
##'                 manually prepared in advance (see 'details' section in the
##'                 buildCodesetAnnotation() help page for more info).
##'
##' @param  cdrDesignData    Path to the .CSV file containing the content extracted from the CDR's 'Design Data' tab
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
            paste0("CDR Design Data .CSV \"", cdrDesignData, "\" not found")
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
##' @param  rlf                     Path to the RLF file
##' @param  cdrDesignData           Path to a manually prepared .CSV export of the
##'                                 "Design Data" tab of the CDR file (optional;
##'                                 see 'details' section below for how the export
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
##' the CDR "Design Data" .CSV is specified, the function expects this .CSV
##' file to be generated from the "Design Data" tab of the original NanoString
##' provided Excel CDR file. This tab needs to be trimmed by skipping the
##' NanoString header and first column containing only integers;
##' the resulting .CSV should contain the actual table (including its header
##' -- beginning with "Customer Identifier"). The function will match and
##' join the .RLF and CDR .CSV using their "ProbeID" and "NSID" fields, and
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
##'	rlf <- system.file("extdata", "RLF", "NQCP_example.rlf", package="NanoStringQCPro")
##'	cdrDesignData <- system.file("extdata", "CDR", "CDR-DesignData.csv", package="NanoStringQCPro")
##'	annot <- buildCodesetAnnotation(rlf, cdrDesignData)
##'
##' @author Dorothee Nickles, Robert Ziman
##'
buildCodesetAnnotation <- function(rlf=NULL,
                                   cdrDesignData=NULL,
                                   removeRedundantCols=TRUE,
                                   addEgAnnotations=FALSE)
{
    if (is.null(rlf)) {
        stop("RLF not specified")
    }
    rlf.df <- readRlf(rlf)

    if (!is.null(cdrDesignData))
    {
        cdr <- readCdrDesignData(cdrDesignData)

        misProbes <- summary(cdr$NSID %in% rlf.df$ProbeID)["FALSE"]
        if (!is.na(misProbes)) {
            warning(sprintf("%s probe(s) in the CDR have no entry in the RLF. These will be excluded from annotation.", misProbes))
        }
    
        annot <- merge(rlf.df, cdr, by.x="ProbeID", by.y="NSID", all.x=TRUE, suffixes=c("", "_CDR"))

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
        annot <- rlf.df
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


setGeneric( "addCodesetAnnotation", function( rccSet, ... ) standardGeneric( "addCodesetAnnotation" ) )

##' @rdname addCodesetAnnotation
##' @aliases addCodesetAnnotation
##'
##' @title          Add NanoString codeset annotation to an RccSet
##' @description    Returns a copy of the input RccSet where the codeset annotation
##'                 has been merged into its fData slot. The merge key for each is a string
##'                 formed from the concatenation of their CodeClass, GeneName, and
##'                 Accession columns ("<CodeClass>_<GeneName>_<Accession>"). For creating
##'                 the codeset annotation object, see buildCodesetAnnotation().
##'
##' @param  rccSet          An RccSet object.
##' @param  annot           Data frame containing the codeset annotation.
##' @param  reorder         Logical indicating whether the probes should be reordered
##'                         according to their barcodes (this can help in identifying
##'                         barcode-specific artifacts -- i.e. background noise).
##' @param  showWarnings    Logical indicating whether or not warnings should be shown, if any.
##'
##' @return
##' A copy of the input RccSet where the codeset annotation has been
##' merged into its fData slot.
##'
##' @export
##'
##' @examples
##' rccDir <- system.file("extdata", "RCC", package="NanoStringQCPro")
##' rccSet <- newRccSet(rccFiles = dir(rccDir, full.names=TRUE))
##' rlf <- system.file("extdata", "RLF", "NQCP_example.rlf", package="NanoStringQCPro")
##' annot <- buildCodesetAnnotation(rlf)
##' rccSet.annotated <- addCodesetAnnotation(rccSet, annot)
##'
##' @author Dorothee Nickles, Robert Ziman
##'
setMethod(
    "addCodesetAnnotation",
    "RccSet",
    function( rccSet, annot, reorder=TRUE, showWarnings=TRUE )
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

        if (reorder == TRUE) {
            if (is.null(fData(rccSet.annot)$BarCode)) {
                stop("reorder == TRUE but featureData BarCode column is missing")
            }
            rccSet.annot <- rccSet.annot[order(fData(rccSet.annot)$BarCode),]
            rccSet.annot <- rccSet.annot[order(fData(rccSet.annot)$CodeClass),]
        }

        return(rccSet.annot)
    }
    )



##' @title
##' Create a new RccSet object
##'
##' @description
##' 
##' This is the main wrapper function for generating an RccSet from NanoString
##' data. The function takes as input a vector of NanoString .RCC files with the
##' raw data or a .CSV file generated via the RCC Collector Tool Export feature
##' of NanoString's nSolver Analysis Software, an optional path to the .RLF file
##' describing the codeset used, optional paths to additional annotation about
##' the features and samples, and details about the experiment. It returns an
##' RccSet object.
##'
##' @param  rccFiles                Vector of paths to .RCC files with the raw count data.
##' 
##' @param  rccCollectorToolExport  Path to a .CSV file generated via the RCC Collector
##'                                 Tool Export feature of NanoString's nSolver
##'                                 Analysis Software. (Note that this is an alternative
##'                                 to rccFiles, and if both arguments are specified at the
##'                                 same time, the function will throw an error.)
##'
##' @param  rlf                     Path to the NanoString .RLF file describing
##'                                 the codeset used in generating the .RCCs. 
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
##'                                 to the phenoData of the output RccSet. All
##'                                 files should be tab-separated and should contain
##'                                 a column labelled "FileName" whose values correspond
##'                                 exactly to the basenames (including .RCC extension)
##'                                 of the files specified in rccFiles or listed in
##'                                 the RCC Collector Tool Export. More than one such file
##'                                 may be used. A SampleType column should be present
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
##' @param  dropPdataCols           Character vector specifying phenoData columns
##'                                 to be dropped from the output object (if empty or
##'                                 NULL, no columns will be dropped).
##'
##' @param  dropFdataCols           Character vector specifying featureData columns
##'                                 to be dropped from the output object (if empty or
##'                                 NULL, no columns will be dropped).
##'
##' @param  experimentData.name     String passed to the 'name' slot of the
##'                                 output RccSet's experimentData.
##' 
##' @param  experimentData.lab      String passed to the 'lab' slot of the
##'                                 output RccSet's experimentData.
##' 
##' @param  experimentData.contact  String passed to the 'contact' slot of the
##'                                 output RccSet's experimentData.
##' 
##' @param  experimentData.title    String passed to the 'title' slot of the
##'                                 output RccSet's experimentData.
##' 
##' @param  experimentData.abstract String passed to the 'abstract' slot of the
##'                                 output RccSet's experimentData.
##' 
##' @param  experimentData.url      String passed to the 'url' slot of the
##'                                 output RccSet's experimentData.
##' 
##' @param  experimentData.other    List passed to the 'other' slot of the
##'                                 output RccSet's experimentData.
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
##' rccDir <- system.file("extdata", "RCC", package="NanoStringQCPro")
##' rccSet <- newRccSet(
##'     rccFiles = dir(rccDir, full.names=TRUE),
##'     rlf = system.file("extdata", "RLF", "NQCP_example.rlf", package="NanoStringQCPro"),
##'     extraPdata = system.file("extdata", "extraPdata", "SampleType.txt", package="NanoStringQCPro"),
##'     blankLabel = "blank",
##'     experimentData.name = "Robert Ziman",
##'     experimentData.lab = "Richard Bourgon",
##'     experimentData.contact = "ziman.robert@@gene.com",
##'     experimentData.title = "NanoStringQCPro example dataset",
##'     experimentData.abstract = "Example data for the NanoStringQCPro package"
##' )
##'
##' @author Robert Ziman
##'
newRccSet <- function(rccFiles,
                      rccCollectorToolExport    = NULL,
                      rlf                       = NULL,
                      cdrDesignData             = NULL,
                      extraPdata                = NULL,
                      blankLabel                = "blank",
                      addEgAnnotations          = FALSE,
                      dropPdataCols             = c("FileVersion",
                                                    "SoftwareVersion",
                                                    "Owner",
                                                    "SystemAPF",
                                                   #"FovCount",             # Required in the QC report
                                                   #"FovCounted",           # Required in the QC report
                                                    "ScannerID",
                                                   #"StagePosition",        # Required in the QC report
                                                   #"BindingDensity",       # Required in the QC report
                                                   #"CartridgeID",          # See comment by corresponding line in createOrUpdateEpProject() (in NanoStringQCProGNE)
                                                    "CartridgeBarcode"
                                                    ),
                      dropFdataCols             = c("CodeClass_codesetAnnot",
                                                    "Accession_codesetAnnot",
                                                    "GeneName_codesetAnnot",
                                                    "Accession_CDR"
                                                    ),
                      experimentData.name       = "",
                      experimentData.lab        = "",
                      experimentData.contact    = "",
                      experimentData.title      = "",
                      experimentData.abstract   = "",
                      experimentData.url        = "",
                      experimentData.other      = list())
{
    if (!missing(rccFiles) && !is.null(rccFiles))
    {
        if(!missing(rccCollectorToolExport) && !is.null(rccCollectorToolExport)) {
            stop("Both rccFiles and rccCollectorToolExport are specified")
        }

        if (grepl('\\.csv$', tolower(rccFiles[1]))) {
            stop("A .CSV file appears to have been specified as the rccFiles argument (if this is the RCC Collector Tool Export, specify its path in the rccCollectorToolExport argument instead)")
        }

        message("Reading RCC files...")
        rccSet <- readRccBatch(rccFiles)

    } else if(!missing(rccCollectorToolExport) && !is.null(rccCollectorToolExport)) {

        message("Reading RCC Collector Tool Export...")
        rccSet <- readRccCollectorToolExport(rccCollectorToolExport)

    } else {
        stop("Either rccFiles or rccCollectorToolExport must be specified")
    }
                                  
    experimentData(rccSet) <- new("MIAME",
                                  name       = experimentData.name,
                                  lab        = experimentData.lab,
                                  contact    = experimentData.contact,
                                  title      = experimentData.title,
                                  abstract   = experimentData.abstract,
                                  url        = experimentData.url,
                                  other      = experimentData.other)

    GeneRLF <- unique(as.character(pData(rccSet)$GeneRLF))
    if (length(GeneRLF) != 1) {
        stop("Expecting exactly one RLF file to be referenced by all RCC files; check RCCs or RCC Collector Tool Export")
    }
    annotation(rccSet) <- GeneRLF
    pData(rccSet)$GeneRLF <- NULL # Removing it from pData since it's now in annotation(rccSet)

    if (!is.null(rlf))
    {
        rlfFile <- basename(rlf)
        if (grepl('\\.rlf$', rlfFile)) {
            rlfFile.noext <- sub('\\.rlf$', '', rlfFile)
        } else if (grepl('\\.RLF$', rlfFile)) {
            rlfFile.noext <- sub('\\.RLF$', '', rlfFile)
        }
        if (GeneRLF != rlfFile.noext) {
            stop(paste0("Specified RLF (\"", rlf, "\") doesn't match that in the input RCCs or RCC Collector Tool Export (\"", GeneRLF, "\")"))
        }
    }

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
                        new_pdata$SampleType <- NULL
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

    if (!is.null(rlf))
    {
        # The codeset annotation should be rebuilt from input RLF/CDR/etc *every* time; the .CSV output to annotPath will thus
        # be just for reference. -RZ 2014-12-14
        annot <- buildCodesetAnnotation(
            rlf                 = rlf,
            cdrDesignData       = cdrDesignData,
            addEgAnnotations    = addEgAnnotations,
            removeRedundantCols = TRUE
            )
        rccSet <- addCodesetAnnotation(
            rccSet          = rccSet,
            annot           = annot,
            reorder         = TRUE,
            showWarnings    = FALSE
            )
    }

    #
    # Record the org.Hs.eg.db version
    #
    if (addEgAnnotations == TRUE) {

        if ("org.Hs.eg.db" %in% names(sessionInfo()$otherPkgs)) {
            org.Hs.eg.db_version <- sessionInfo()$otherPkgs$org.Hs.eg.db$Version

        } else if ("org.Hs.eg.db" %in% names(sessionInfo()$loadedOnly)) {
            org.Hs.eg.db_version <- sessionInfo()$loadedOnly$org.Hs.eg.db$Version

        } else {
            warning("Attempted to record org.Hs.eg.db version but did not find it in the expected locations in sessionInfo()")
        }

        preproc(rccSet)$org.Hs.eg.db_version <- org.Hs.eg.db_version

    } else {
        preproc(rccSet)$org.Hs.eg.db_version <- NA
    }

    # Remove superfluous columns (see also additional code in createOrUploadEpProject() (in NanoStringQCProGNE)).
    if (!is.null(dropPdataCols) && (length(dropPdataCols) > 0)) {
        for (x in dropPdataCols) {
            pData(rccSet)[[ x ]] <- NULL
        }
    }
    if (!is.null(dropFdataCols) && (length(dropFdataCols) > 0)) {
        for (x in dropFdataCols) {
            fData(rccSet)[[ x ]] <- NULL
        }
    }

    checkRccSet(rccSet, reportWarnings=TRUE, showMessages=TRUE)
    
    return(rccSet)

}



setGeneric( "posCtrlNorm", function( rccSet, ... ) standardGeneric( "posCtrlNorm" ) )

##' @rdname posCtrlNorm
##' @aliases posCtrlNorm
##'
##' @title
##' Positive control normalization
##'
##' @description
##' Applies positive control normalization to the data in an RccSet object.
##'
##' @param rccSet
##' An RccSet object.
##'
##' @param summaryFunction
##' Function to be used for the normalization (e.g. "mean", "median", or
##' "sum"). User-defined functions similar to these can be specified here as
##' well.
##'
##' @param quietly
##' Logical. If TRUE, messages and warnings will not be shown.
##'
##' @return
##' A copy of the input RccSet that has count data adjusted by positive control
##' counts. The positive control scaling factor is recorded in PosFactor in the
##' output's phenoData (if this column already exists in the input, it will
##' be overwritten in the output copy).
##'
##' @export
##'
##' @examples
##' data(example_rccSet)
##' pcnorm_example_rccSet <- posCtrlNorm(example_rccSet)
##'
##' @author Dorothee Nickles
##'
setMethod(
    "posCtrlNorm",
    "RccSet",
    function(rccSet,
             summaryFunction = "sum",
             quietly = FALSE)
    {
        if (!quietly && "posCtrlData" %in% ls(assayData(rccSet)))
            warning("Input already contains positive control normalized data")

        prccSet <- copyRccSet(rccSet)

        if (class(summaryFunction) == "character") {
            fun <- get(summaryFunction)
            summaryFunction.char <- summaryFunction
        }
        else {
            fun <- summaryFunction
            summaryFunction.char <- as.character(quote(summaryFunction))
        }

        posSignal <- apply(exprs(prccSet)[fData(prccSet)$CodeClass == "Positive", ], 2, fun)
        signalMean <- mean(posSignal)
        posFactor <- signalMean/posSignal

        if (!quietly && "posCtrlData" %in% ls(assayData(prccSet)))
            warning("Input already contains positive-control normalized data")

        assayData(prccSet)$posCtrlData <- t(apply(assayData(prccSet)$exprs, 1, function(x) { x * posFactor }))
        pData(prccSet)$PosFactor <- posFactor
        PosFactor.metadata_rownum <- which(rownames(varMetadata(phenoData(prccSet))) == "PosFactor")
        varMetadata(phenoData(prccSet))$labelDescription[ PosFactor.metadata_rownum ] <- "Positive control scaling factor"

        preproc(prccSet)$posCtrlData_summaryFunction <- summaryFunction.char

        return(prccSet)
    }
    )



setGeneric( "presAbsCall", function( rccSet, ... ) standardGeneric( "presAbsCall" ) )

##' @rdname presAbsCall
##' @aliases presAbsCall
##'
##' @title
##' Presence/absence call
##'
##' @description
##' Adds a matrix to assayData (`paData') which indicates the presence/absence
##' call for each gene in each sample using the background estimates and
##' a stringency value. A gene is considered present in a sample if
##' its count in that sample exceeds the corresponding background estimate
##' times the stringency. The count values can be taken from either the
##' positive control normalized data or the raw data (see
##' the inputMatrix agrument). If the input doesn't contain background-corrected
##' data, an error will be generated.
##'
##' @param rccSet
##' An RccSet with background-corrected data.
##'
##' @param stringency
##' Multiplier to use in establishing the presence/absence call as
##' mentioned in the description.
##'
##' @param inputMatrix
##' Name of the matrix in the RccSet's assayData on which to apply the
##' presence/absence call (either "posCtrlData" or "exprs").
##'
##' @param quietly
##' Logical. If TRUE, messages and warnings will not be shown.
##' 
##' @return
##' A copy of the input is returned with a new matrix named `paData' added to
##' the assayData that contains the presence/absence calls.
##'
##' @export
##'
##' @examples
##' data(example_rccSet)
##' pcnorm_rccSet <- posCtrlNorm(example_rccSet)
##' bgEst <- getBackground(pcnorm_rccSet)
##' bgcorr_rccSet <- subtractBackground(pcnorm_rccSet, bgEst)
##' pa_rccset <- presAbsCall(bgcorr_rccSet)
##'
setMethod(
    "presAbsCall",
    "RccSet",
    function(rccSet,
             stringency = 2,
             inputMatrix = c("posCtrlData", "exprs"),
             quietly = FALSE
            )
    {
        inputMatrix <- match.arg(inputMatrix)

        if (!("bgEstimates" %in% ls(assayData(rccSet))))
            stop("Background-corrected data required but not found in the input RccSet")

        bgEstimates <- assayData(rccSet)$bgEstimates

        M <- assayData(rccSet)[[ inputMatrix ]]
        if (is.null(M)) 
            stop(sprintf("Specified input matrix ('%s') is not present in the RccSet's assayData", inputMatrix)) 

        if (!quietly && "paData" %in% ls(assayData(rccSet)))
            warning("Input already contains a presence/absence matrix")

        parccSet <- copyRccSet(rccSet)
        assayData(parccSet)$paData <- (M > bgEstimates * stringency)
        preproc(parccSet)$paData_stringency <- stringency
        preproc(parccSet)$paData_inputMatrix <- inputMatrix

        return(parccSet)
    }
    )



setGeneric( "contentNorm", function( rccSet, ... ) standardGeneric( "contentNorm" ) )

##' @rdname contentNorm
##' @aliases contentNorm
##'
##' @title
##' Content normalization
##'
##' @description
##' Performs content normalization on the given RccSet.
##'
##' @param rccSet
##' An RccSet.
##'
##' @param method
##' Specifies the features to be used for normalization. "global" indicates that all
##' features should be used and "housekeeping" indicates that only housekeeping
##' features should be used. If "housekeeping" is specified and the `hk' argument
##' (below) is also specified, then the features indicated by `hk' will be used.
##' If "housekeeping" is specified and `hk' is left NULL, then the default
##' housekeeping features (i.e. those with CodeClass == "Housekeeping") will be used.
##'
##' @param summaryFunction
##' Character specifying the summary function to apply to the selected features
##' (e.g. "mean" or "median"). User-defined functions similar to these can be
##' specified here as well.
##'
##' @param hk
##' Logical vector defining, for each feature, whether or not it shall
##' be used for housekeeping normalization if housekeeping is specified as the
##' normalization method.
##'
##' @param inputMatrix
##' Name of the matrix in the RccSet's assayData to use as input for performing
##' content normalization (one of "exprs", "posCtrlData", or "bgCorrData"). If
##' posCtrlData or bgCorrData are specified but not found in the assayData, an
##' error will be generated.
##'
##' @param quietly
##' Boolean specifying whether or not messages and warnings should be omitted.
##'
##' @return
##' A copy of the input is returned with a new matrix named `normData' added to
##' the assayData that contains the content-normalized counts. (\bold{NOTE}: normData
##' contains values on a log2 scale while all other matrices in assayData are
##' on a linear scale.) If housekeeping is specified as the normalization method,
##' then the housekeeping features used will be recorded in the returned RccSet in
##' a new featureData column named `Housekeeping'. Parameters specified in the
##' function call are also recorded in the output's experimentData@@preprocessing
##' list.
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
##' gmnorm_example_rccSet <- contentNorm(bgcorr_example_rccSet, method="global",
##'     inputMatrix="exprs")
##' hknorm_example_rccSet <- contentNorm(bgcorr_example_rccSet, method="housekeeping",
##'     summaryFunction="mean")
##'
##' @author Dorothee Nickles
##'
setMethod(
    "contentNorm",
    "RccSet",
    function(rccSet,
             method = c("global", "housekeeping"),
             summaryFunction = "median",
             hk = NULL,
             inputMatrix = c("bgCorrData", "posCtrlData", "exprs"),
             quietly = FALSE)
    {
        if (class(summaryFunction) == "character") {
            sfun <- get(summaryFunction)
            summaryFunction.char <- summaryFunction
        } else {
            sfun <- summaryFunction
            summaryFunction.char <- as.character(quote(summaryFunction))
        }

        inputMatrix <- match.arg(inputMatrix)

        M <- assayData(rccSet)[[ inputMatrix ]]
        if (is.null(M))
            stop(sprintf("Specified input matrix ('%s') is not present in the RccSet's assayData", inputMatrix)) 

        if (!quietly && "normData" %in% ls(assayData(rccSet)))
            warning("Input already contains content-normalized data")

        # Apply the log2 transformation. Note that all the code below is thus
        # operating on values that are on a log2 scale.
        M <- log2(M)

        if (method == "global") {

            nFactors <- apply( M[ fData(rccSet)$CodeClass == "Endogenous", ], 2, sfun )
            nFact <-  nFactors - mean(nFactors)

        } else { # method == "housekeeping"

            if (is.null(hk)) {
                hk <- (fData(rccSet)$CodeClass == "Housekeeping")
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
            if (!quietly && (sum(hk) < 3)) {
                warning("Less than three housekeeping features are defined")
            }

            fData(rccSet)$Housekeeping <- hk

            if (sum(hk) > 1)
                nFact <- apply( M[hk, ], 2, sfun )
            else
                nFact <- sfun( M[hk, ] )
        }

        M <- sweep(M, 2, nFact, "-")      # N.B. "-" is valid here since the opration is being performed
                                          # with log2-based values, so the matrix is effectively being
                                          # divided by nFact.

        nrccSet <- copyRccSet(rccSet)

        assayData(nrccSet)$normData <- M

        preproc(nrccSet)$normData_method          <- method
        preproc(nrccSet)$normData_summaryFunction <- summaryFunction.char
        if (method == "housekeeping") {
            preproc(nrccSet)$normData_hkgenes     <- fData(nrccSet)$GeneName[ hk ]
            preproc(nrccSet)$normData_hkfeatures  <- featureNames(nrccSet)[ hk ]
        } else {
            preproc(nrccSet)$normData_hkgenes     <- NA
            preproc(nrccSet)$normData_hkfeatures  <- NA
        }
        preproc(nrccSet)$normData_inputMatrix     <- inputMatrix

        return(nrccSet)
    }
    )



setGeneric( "preprocRccSet", function( rccSet, ... ) standardGeneric( "preprocRccSet" ) )

##' @rdname preprocRccSet
##' @aliases preprocRccSet
##'
##' @title
##' Preprocess an RccSet
##'
##' @description
##' This function is a wrapper to perform any combination of positive control
##' normalization, background correction, and content normalization on the
##' input RccSet. For each completed preprocessing step, a matrix is added to
##' the assayData of the resulting RccSet object:
##'
##' \itemize{
##'   \item posCtrlData: expression data after positive control normalization
##'   \item bgEstimates: background estimates
##'   \item bgCorrData: expression data after positive control normalization and
##'         background correction
##'   \item normData: expression data after positive control normalization,
##'         background correction, and content normalization
##' }
##'
##' (\bold{NOTE}: normData is on a log2 scale while all the other matrices are
##' on a linear scale.)
##'
##' If any step is omitted, the corresponding matrix will not be present in
##' the output's assayData. The parameters for all steps are recorded in the
##' output's experimentData@@preprocessing list (accessible through
##' preproc(rccSet) where rccSet is an RccSet output by this function). In
##' addition:
##'
##' \itemize{
##'   \item If positive control normalization is performed, a column named
##'         'PosCtrl' is added to the output's phenoData to record the
##'         positive control scaling factors.
##'   \item If the presence/absence call is performed, a matrix named `paData'
##'         is added to the output's assayData to indicate the
##'         presence/absence of each feature in each sample. See the `pa'
##'         argument for details.
##'   \item If housekeeping normalization is performed, a column labeled
##'         `Housekeeping' is added to the featureData to indicate which
##'         features were used for it.
##' }
##'
##' @details
##' For more information on the rationale behind the recommended
##' preprocessing and normalization steps, please see the vignette.
##'
##' @param rccSet
##' An RccSet.
##'
##' @param doPosCtrlNorm
##' Boolean specifying whether or not to perform positive control normalization.
##' (`pcd' is short for `posCtrlData', the matrix which gets added to assayData
##' when this step is performed.)
##'
##' @param doBackground
##' Boolean specifying whether or not to perform background correction.
##'
##' @param doPresAbs
##' Boolean specifying whether or not the presence/absence call should be
##' performed. For details, see presAbsCall().
##'
##' @param doContentNorm
##' Boolean specifying whether or not content normalization should be performed.
##'
##' @param pcnSummaryFunction
##' Function to be used for the positive control normalization (e.g. "mean",
##' "median", or "sum"). User-defined functions similar to these can be
##' specified here as well.
##'
##' @param bgReference
##' Measurements to use for background estimates: either "blank" (for blank
##' samples), "negatives" (for negative control probes), or "both". For
##' details on exactly how the background estimates are computed in each
##' case, see getBackground().
##'
##' @param bgSummaryFunction
##' Summary function for background measurements (e.g. "mean" or "median").
##' User-defined functions similar to these can be specified here as well.
##'
##' @param bgStringency
##' Factor by which deviation (SD or MAD) of the summarization output will be
##' multiplied to obtain final background estimates.
##'
##' @param nSolverBackground.w1
##' Value to use for the 'w1' argument to nSolverBackground(). (Only takes
##' effect if bgReference == "both"; see getBackground().)
##'
##' @param nSolverBackground.shrink
##' Value to use for the 'shrink' argument to nSolverBackground(). (Only takes
##' effect if bgReference == "both"; see getBackground().)
##'
##' @param paStringency
##' Multiplier to use in establishing the presence/absence call. For details,
##' see presAbsCall().
##'
##' @param normMethod
##' Specifies the features to be used for content normalization. "global" indicates that all
##' features should be used and "housekeeping" indicates that only housekeeping
##' features should be used. If "housekeeping" is specified and the `hk' argument
##' (below) is also specified, then the features indicated by `hk' will be used.
##' If "housekeeping" is specified and `hk' is left NULL, then the default
##' housekeeping features (i.e. those with CodeClass == "Housekeeping") will be used.
##'
##' @param normSummaryFunction
##' Character specifying the summary function to apply to the selected features
##' (e.g. "mean" or "median") during the content normalization step. User-defined
##' functions similar to these can be specified here as well.
##'
##' @param hkgenes
##' Character vector with gene symbols to be used for content normalization if
##' housekeeping is specified as the normalization method. If specified, all
##' features that match any of the specified symbols will be used. (To specify
##' specific features, use the `hkfeatures' argument instead; see below.)
##'
##' @param hkfeatures
##' Character vector with full feature names
##' ("<CodeClass>_<GeneName>_<Accession>", e.g. "Endogenous_ACTG1_NM_001614.1")
##' to be used for content normalization if housekeeping is specified as the normalization
##' method. (Note: if this argument is specified at the same time as `hkgenes',
##' an error will be thrown.)
##'
##' @param quietly
##' Boolean specifying whether or not messages and warnings should be omitted.
##'
##' @return
##' A copy of the input RccSet with additional matrices in the assayData for each
##' successive preprocessing step along with parameters for each step recorded in the
##' experimentData@@preprocessing list.
##'
##' @export
##'
##' @examples
##' data(example_rccSet)
##' hknorm_example_rccSet <- preprocRccSet(example_rccSet)
##'
##' @references
##' \href{http://www.nanostring.com/media/pdf/MAN_nCounter_Gene_Expression_Data_Analysis_Guidelines.pdf}{NanoString nCounter(R) Expression Data Analysis Guide (2012)}
##'
##' @author Dorothee Nickles, Robert Ziman
##'
setMethod(
    "preprocRccSet",
    "RccSet",
    function(rccSet,
             doPosCtrlNorm = TRUE,      # pcd = TRUE,
             doBackground = TRUE,       # bg = TRUE,
             doPresAbs = TRUE,          # pa = TRUE,
             doContentNorm = TRUE,      # cn = TRUE,
             pcnSummaryFunction = "sum",
             bgReference = c("both", "blanks", "negatives"),
             bgSummaryFunction = "median",
             bgStringency = 1,
             nSolverBackground.w1 = 2.18,
             nSolverBackground.shrink = TRUE,
             paStringency = 2,
             normMethod = c("global", "housekeeping"),
             normSummaryFunction = "median",
             hkgenes = NULL,
             hkfeatures = NULL,
             quietly = FALSE)
    {
        bgReference <- match.arg(bgReference)
        normMethod <- match.arg(normMethod)
       
        # Positive control normalization

        if (doPosCtrlNorm) {

            temp_rccSet_1 <- posCtrlNorm(rccSet,
                                         summaryFunction = pcnSummaryFunction)

            bgInputMatrix <- "posCtrlData"

        } else {

            temp_rccSet_1 <- copyRccSet(rccSet)
            preproc(temp_rccSet_1)$posCtrlData_summaryFunction <- NA

            bgInputMatrix <- "exprs"
        }

        # Background correction

        if (doBackground) {

            bgEstimates <- getBackground(rccSet                   = temp_rccSet_1,
                                         bgReference              = bgReference,
                                         summaryFunction          = bgSummaryFunction,
                                         stringency               = bgStringency,
                                         nSolverBackground.w1     = nSolverBackground.w1,
                                         nSolverBackground.shrink = nSolverBackground.shrink,
                                         inputMatrix              = bgInputMatrix)

            temp_rccSet_2 <- subtractBackground(rccSet            = temp_rccSet_1,
                                                bgEstimates       = bgEstimates,
                                                bgEstimatesParams = list(bgReference              = bgReference,
                                                                         summaryFunction          = bgSummaryFunction,
                                                                         stringency               = bgStringency,
                                                                         nSolverBackground.w1     = nSolverBackground.w1,
                                                                         nSolverBackground.shrink = nSolverBackground.shrink,
                                                                         inputMatrix              = bgInputMatrix),
                                                inputMatrix       = bgInputMatrix)

            normInputMatrix <- "bgCorrData"

        } else {

            temp_rccSet_2 <- copyRccSet(temp_rccSet_1)
            preproc(temp_rccSet_2)$bgEstimatesParams <- list(bgReference              = NA,
                                                             summaryFunction          = NA,
                                                             stringency               = NA,
                                                             nSolverBackground.w1     = NA,
                                                             nSolverBackground.shrink = NA,
                                                             inputMatrix              = NA)

            if (doPosCtrlNorm)
                normInputMatrix <- "posCtrlData"
            else
                normInputMatrix <- "exprs"

        }

        # Presence/absence matrix

        if (doPresAbs) {
            if (doBackground) {
                if (doPosCtrlNorm)
                    temp_rccSet_3 <- presAbsCall(temp_rccSet_2, stringency=paStringency, inputMatrix="posCtrlData")
                else
                    temp_rccSet_3 <- presAbsCall(temp_rccSet_2, stringency=paStringency, inputMatrix="exprs")
            }
        } else {
            temp_rccSet_3 <- copyRccSet(temp_rccSet_2)
            preproc(temp_rccSet_3)$paData_stringency <- NA
            preproc(temp_rccSet_3)$paData_inputMatrix <- NA
        }

        # Content normalization

        if (doContentNorm) {

            if (normMethod == "housekeeping")
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

            } else {

                hk <- NULL

            }

            temp_rccSet_4 <- contentNorm(rccSet          = temp_rccSet_3,
                                         method          = normMethod,
                                         summaryFunction = normSummaryFunction,
                                         hk              = hk,
                                         inputMatrix     = normInputMatrix)

        } else {

            temp_rccSet_4 <- copyRccSet(temp_rccSet_3)

            preproc(temp_rccSet_4)$normData_method          <- NA
            preproc(temp_rccSet_4)$normData_summaryFunction <- NA
            preproc(temp_rccSet_4)$normData_hkgenes         <- NA
            preproc(temp_rccSet_4)$normData_hkfeatures      <- NA
            preproc(temp_rccSet_4)$normData_inputMatrix     <- NA

        }

        return(temp_rccSet_4)
    }
    )



setGeneric( "addQCFlags", function( rccSet, ... ) standardGeneric( "addQCFlags" ) )

##' @rdname addQCFlags
##' @aliases addQCFlags
##'
##' @title          Add sample QC flags to an rccSet
##' @description    Returns a copy of the input RccSet with columns added to pData
##'                 from the provided sample QC flag annotation file. (That file is
##'                 produced by makeQCReport(); see its help page for more
##'                 details.)
##'
##' @param  rccSet      An RccSet object
##' @param  flagFile    Path to a sample QC flag file as generated by the
##'                     NanoStringQCPro QC report (see makeQCReport())
##'
##' @return
##' A copy of the input RccSet with columns added to pData from the QC
##' flag file.
##'
##' @author Dorothee Nickles
##'
setMethod(
    "addQCFlags",
    "RccSet",
    function(rccSet, flagFile)
    {
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
    )

