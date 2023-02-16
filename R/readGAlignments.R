### =========================================================================
### readGAlignments() and related functions
### -------------------------------------------------------------------------


### A "flag filter" is represented as a 'flag' vector of length 2 with names
### keep0 and keep1. The .combineBamFlagFilters() function performs a logical
### AND between 2 "flag filters". It returns a "flag filter".
.combineBamFlagFilters <- function(flagfilterA, flagfilterB)
{
    if (!identical(names(flagfilterA), c("keep0", "keep1"))
     || !identical(names(flagfilterB), c("keep0", "keep1")))
        stop("input must be BAM flag filters")
    ans <- bamFlagAND(flagfilterA, flagfilterB)
    if (!all(bamFlagAsBitMatrix(ans[["keep0"]]) |
             bamFlagAsBitMatrix(ans[["keep1"]])))
        stop("BAM flag filters to combine are incompatible")
    ans
}

.normargParam <- function(param, flag0, what0)
{
    if (is.null(param))
        param <- ScanBamParam()
    bamFlag(param) <-
        .combineBamFlagFilters(bamFlag(param, asInteger=TRUE), flag0)
    bamWhat(param) <- union(bamWhat(param), what0)
    param
}

### 'x' must be a GAlignments object.
.bindExtraData <- function(x, use.names, param, bamcols,
                              with.which_label=FALSE)
{
    if (use.names)
        names(x) <- bamcols$qname
    if (is.null(param))
        return(x)
    colnames <- c(bamWhat(param), bamTag(param))
    if (with.which_label)
        colnames <- c(colnames, "which_label")
    if (length(colnames) != 0L) {
        df <- do.call(DataFrame, bamcols[colnames])
        ## Sadly, the DataFrame() constructor is mangling the duplicated
        ## colnames to make them unique. Since we of course don't want this,
        ## we need to fix them.
        colnames(df) <- colnames
        mcols(x) <- df
    }
    x
}

.load_bamcols_from_BamFile <- function(file, param, what0, 
                                       with.which_label=FALSE)
{
    flag0 <- scanBamFlag(isUnmappedQuery=FALSE)
    param <- .normargParam(param, flag0, what0)
    res <- scanBam(file, param=param)
    if (length(res) == 0L)  # should never happen
        stop("scanBam() returned a list of length zero")
    Rsamtools:::.load_bamcols_from_scanBam_res(res, param,
                    with.which_label=with.which_label)
}

.load_seqlengths_from_BamFile <- function(file, seqlevels)
{
    seqlengths <- seqlengths(file)
    if (is.null(seqlengths))
        return(NULL)
    bad <- setdiff(seqlevels, names(seqlengths))
    if (length(bad) == 0L)
        return(seqlengths)
    bad <- paste(bad, collapse="' '")
    msg <- sprintf("'rname' lengths not in BamFile header; seqlengths not used\n  file: %s\n  missing rname(s): '%s'",
                   path(file), bad)
    warning(msg)
    NULL
}

.open_BamFile <- function(file, index=file, asMates=FALSE, param=NULL)
{
    if ((missing(index) || identical(index, file))
     && (is.null(param) || length(bamWhich(param)) == 0L))
        index <- character(0)
    open(BamFile(file, index=index, asMates=asMates), "rb")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readGAlignments()
###

setGeneric("readGAlignments", signature="file",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
        standardGeneric("readGAlignments")
)

.readGAlignments.BamFile <- function(file, index=file,
                                     use.names=FALSE, param=NULL,
                                     with.which_label=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (is.null(param))
        param <- ScanBamParam()
    if (!asMates(file))
        bamWhat(param) <- setdiff(bamWhat(param), 
                                  c("groupid", "mate_status"))
    what0 <- c("rname", "strand", "pos", "cigar")
    if (use.names)
        what0 <- c(what0, "qname")
    bamcols <- .load_bamcols_from_BamFile(file, param, what0,
                                          with.which_label=with.which_label)
    seqlengths <- .load_seqlengths_from_BamFile(file,
                                                levels(bamcols[["rname"]]))
    ans <- GAlignments(seqnames=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       seqlengths=seqlengths)
    .bindExtraData(ans, use.names, param, bamcols,
                        with.which_label=with.which_label)
}

setMethod("readGAlignments", "BamFile", .readGAlignments.BamFile)

setMethod("readGAlignments", "character",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
    {
        bam <- .open_BamFile(file, index=index, param=param)
        on.exit(close(bam))
        readGAlignments(bam, character(0),
                        use.names=use.names, param=param,
                        with.which_label=with.which_label)
    }
)

setMethod("readGAlignments", "BamViews",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
    {
        if (missing(index))
            index <- bamIndicies(file)
        if (is.null(param)) {
            param <- ScanBamParam(which=bamRanges(file))
        } else if (!identical(bamRanges(file), bamWhich(param))) {
            warning("'bamRanges(file)' and 'bamWhich(param)' differ; using 'bamRanges(file)'")
            bamWhich(param) <- bamRanges(file)
        }
        fun <- function(i, bamViews, verbose)
            readGAlignments(file=bamPaths(bamViews)[i],
                            index=bamIndicies(bamViews)[i],
                            use.names=use.names,
                            param=param,
                            with.which_label=with.which_label)
        ### Rsamtools:::.BamViews_delegate requires the ShortRead package!
        Rsamtools:::.BamViews_delegate("readGAlignments", file, fun)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readGAlignmentPairs()
###

setGeneric("readGAlignmentPairs", signature="file",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE, strandMode=1)
        standardGeneric("readGAlignmentPairs")
)

### 'use.mcols' can be TRUE, FALSE, or a character vector specifying the
### *inner* metadata columns to return, i.e., the metadata columns to set on
### the 2 halves of the returned GAlignmentPairs object.
.make_GAlignmentPairs_from_GAlignments <- function(gal, strandMode=1L,
                                                        use.mcols=FALSE)
{
    groupid <- Rle(mcols(gal, use.names=FALSE)[ , "groupid"])
    stopifnot(isStrictlySorted(runValue(groupid)))
    mate_status <- Rle(mcols(gal, use.names=FALSE)[ , "mate_status"])

    ## Dump alignments with "ambiguous" mate status.
    flushDumpedAlignments()
    dumped_gal <- gal[mate_status == "ambiguous"]
    dumped_count <- length(dumped_gal)
    if (dumped_count != 0L) {
        dumpAlignments(dumped_gal)
        warning("  ", dumped_count, " alignments with ambiguous pairing ",
                "were dumped.\n    Use 'getDumpedAlignments()' to retrieve ",
                "them from the dump environment.")
    }

    ## Keep alignments that (1) have a mate status set to "mated" and
    ## (2) belong to a group of size 2.
    ok <- mate_status == "mated" &
          Rle(runLength(groupid) == 2L, runLength(groupid))
    if (!all(ok))
        gal <- gal[ok]

    ## Check flag bits 0x40 and 0x80.
    flag <- mcols(gal, use.names=FALSE)[ , "flag"]
    is_first_mate <- as.logical(bamFlagAsBitMatrix(flag,
                                                   bitnames="isFirstMateRead"))
    is_last_mate <- as.logical(bamFlagAsBitMatrix(flag,
                                                  bitnames="isSecondMateRead"))
    bits_0x40_0x80_are_ok <- is_first_mate != is_last_mate
    stopifnot(all(bits_0x40_0x80_are_ok))

    ## Split and order the pairs by ascending start position of the first mate.
    idx1 <- which(is_first_mate)
    idx2 <- which(is_last_mate)
    oo1 <- orderIntegerPairs(as.integer(gal@seqnames)[idx1], gal@start[idx1])
    idx1 <- idx1[oo1]
    idx2 <- idx2[oo1]
    ans_first <- gal[idx1]
    ans_last <- gal[idx2]
    groupid1 <- mcols(ans_first, use.names=FALSE)[ , "groupid"]
    groupid2 <- mcols(ans_last, use.names=FALSE)[ , "groupid"]
    stopifnot(identical(groupid1, groupid2))

    ## Drop the names.
    ans_names <- names(ans_first)
    names(ans_first) <- names(ans_last) <- NULL

    ## Check isProperPair (0x2) and isSecondaryAlignment (0x100) flag bits.
    flag1 <- mcols(ans_first, use.names=FALSE)[ , "flag"]
    flag2 <- mcols(ans_last, use.names=FALSE)[ , "flag"]
    is_proper1 <- bamFlagAsBitMatrix(flag1, bitnames="isProperPair")
    is_proper2 <- bamFlagAsBitMatrix(flag2, bitnames="isProperPair")
    stopifnot(identical(is_proper1, is_proper2))
    is_secondary1 <- bamFlagAsBitMatrix(flag1, bitnames="isSecondaryAlignment")
    is_secondary2 <- bamFlagAsBitMatrix(flag2, bitnames="isSecondaryAlignment")
    stopifnot(identical(is_secondary1, is_secondary2))

    ## Make the GAlignmentPairs object and return it.
    if (is.character(use.mcols)) {
        mcols(ans_first) <- mcols(ans_first, use.names=FALSE)[use.mcols]
        mcols(ans_last) <- mcols(ans_last, use.names=FALSE)[use.mcols]
    } else if (!use.mcols) {
        mcols(ans_first) <- mcols(ans_last) <- NULL
    }
    GAlignmentPairs(ans_first, ans_last, strandMode=strandMode,
                    isProperPair=as.logical(is_proper1), names=ans_names)
}

.readGAlignmentPairs.BamFile <- function(file, index=file,
                                         use.names=FALSE, param=NULL,
                                         with.which_label=FALSE,
                                         strandMode=1)
{
    if (!asMates(file)) {
        asMates(file) <- TRUE
        ## This is required because BamFile objects have a pass-by-address
        ## semantic.
        on.exit(asMates(file) <- FALSE)
    }
    if (is.null(param))
        param <- ScanBamParam()
    flag0 <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE)
    what0 <- c("flag", "groupid", "mate_status")
    param2 <- .normargParam(param, flag0, what0)
    gal <- readGAlignments(file, use.names=use.names, param=param2,
                                 with.which_label=with.which_label)
    use.mcols <- c(bamWhat(param), bamTag(param))
    if (with.which_label)
        use.mcols <- c(use.mcols, "which_label")
    .make_GAlignmentPairs_from_GAlignments(gal, strandMode=strandMode,
                                                use.mcols=use.mcols)
}

setMethod("readGAlignmentPairs", "BamFile", .readGAlignmentPairs.BamFile)

setMethod("readGAlignmentPairs", "character",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE, strandMode=1)
    {
        bam <- .open_BamFile(file, index=index, asMates=TRUE, param=param)
        on.exit(close(bam))
        readGAlignmentPairs(bam, character(0),
                            use.names=use.names, param=param,
                            with.which_label=with.which_label,
                            strandMode=strandMode)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readGAlignmentsList()
###

setGeneric("readGAlignmentsList", signature="file",
    function(file, index=file, use.names=FALSE, param=ScanBamParam(),
                   with.which_label=FALSE, strandMode=1)
        standardGeneric("readGAlignmentsList")
)

.setRealStrand <- function(gal, param, strandMode) {
    if (strandMode == 0L)
        strand(gal) <- Rle(strand("*"), length(gal))
    else {
        gal_mcols <- mcols(gal, use.names=FALSE)
        if (is.null(gal_mcols$flag))
            warning("Flag information missing in GAlignmentsList object. Strand information might not be accurate.")
        else {
            mask_first_mate <- bamFlagTest(gal_mcols$flag, "isFirstMateRead")
            if (strandMode == 1L)
                strand(gal[!mask_first_mate]) <-
                    invertStrand(strand(gal[!mask_first_mate]))
            else if (strandMode == 2L)
                strand(gal[mask_first_mate]) <-
                    invertStrand(strand(gal[mask_first_mate]))
            else
                stop("strandMode should be either 0, 1 or 2.")
        }
    }
    ## if the user didn't request the 'flag' info
    ## then remove it to reduce memory footprint
    if (!"flag" %in% bamWhat(param))
        mcols(gal)$flag <- NULL
    gal
}
.matesFromBam <- function(file, use.names, param, what0, with.which_label,
                          strandMode)
{
    bamcols <- .load_bamcols_from_BamFile(file, param, what0,
                                          with.which_label=with.which_label)
    seqlengths <- .load_seqlengths_from_BamFile(file, levels(bamcols$rname))
    gal <- GAlignments(seqnames=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       seqlengths=seqlengths)
    flag0 <- scanBamFlag()
    what0 <- "flag"
    param2 <- .normargParam(param, flag0, what0)
    gal <- .bindExtraData(gal, use.names=FALSE, param2, bamcols,
                          with.which_label=with.which_label)
    gal <- .setRealStrand(gal, param, strandMode)
    if (asMates(file)) {
        f <- factor(bamcols$groupid)
        gal <- unname(split(gal, f))
        mcols(gal)$mate_status <- 
            bamcols$mate_status[match(levels(f), bamcols$groupid)]
    } else {
        ## groupid=NULL when asMates=FALSE
        gal <- unname(split(gal, seq_along(gal)))
    }
    if (use.names)
        names(gal) <- unique(splitAsList(bamcols$qname, bamcols$groupid))
    gal
}

.readGAlignmentsList.BamFile <- function(file, index=file,
                                         use.names=FALSE, param=ScanBamParam(),
                                         with.which_label=FALSE,
                                         strandMode=1L)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!asMates(file))
        bamWhat(param) <- setdiff(bamWhat(param), 
                                  c("groupid", "mate_status"))
    what0 <- c("rname", "strand", "pos", "cigar", "groupid", "mate_status",
               "flag")
    if (use.names)
        what0 <- c(what0, "qname")
    .matesFromBam(file, use.names, param, what0, with.which_label, strandMode)
}

setMethod("readGAlignmentsList", "BamFile", .readGAlignmentsList.BamFile)

setMethod("readGAlignmentsList", "character",
    function(file, index=file, use.names=FALSE, param=ScanBamParam(),
                   with.which_label=FALSE, strandMode=1)
    {
        bam <- .open_BamFile(file, index=index, asMates=TRUE, param=param)
        on.exit(close(bam))
        readGAlignmentsList(bam, character(0),
                            use.names=use.names, param=param,
                            with.which_label=with.which_label,
                            strandMode=strandMode)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readGappedReads()
###

setGeneric("readGappedReads", signature="file",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
        standardGeneric("readGappedReads")
)

.readGappedReads.BamFile <- function(file, index=file,
                                     use.names=FALSE, param=NULL,
                                     with.which_label=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (is.null(param))
        param <- ScanBamParam()
    if (!asMates(file))
        bamWhat(param) <- setdiff(bamWhat(param), 
                                  c("groupid", "mate_status"))
    what0 <- c("rname", "strand", "pos", "cigar", "seq")
    if (use.names)
        what0 <- c(what0, "qname")
    bamcols <- .load_bamcols_from_BamFile(file, param, what0,
                                          with.which_label=with.which_label)
    seqlengths <- .load_seqlengths_from_BamFile(file,
                                                levels(bamcols[["rname"]]))
    ans <- GappedReads(seqnames=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       qseq=bamcols$seq, seqlengths=seqlengths)
    .bindExtraData(ans, use.names, param, bamcols,
                   with.which_label=with.which_label)
}

setMethod("readGappedReads", "BamFile", .readGappedReads.BamFile)

setMethod("readGappedReads", "character",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
    {
        bam <- .open_BamFile(file, index=index, param=param)
        on.exit(close(bam))
        readGappedReads(bam, character(0),
                        use.names=use.names, param=param,
                        with.which_label=with.which_label)
    }
)

