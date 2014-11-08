### =========================================================================
### readGAlignments() and related functions
### -------------------------------------------------------------------------


setGeneric("readGAlignmentsFromBam", signature="file",
    function(file, index=file, ..., use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
        standardGeneric("readGAlignmentsFromBam")
)

setGeneric("readGAlignmentPairsFromBam", signature="file",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
        standardGeneric("readGAlignmentPairsFromBam")
)

setGeneric("readGAlignmentsListFromBam", signature="file",
    function(file, index=file, ..., use.names=FALSE, param=ScanBamParam(),
                   with.which_label=FALSE)
        standardGeneric("readGAlignmentsListFromBam")
)

setGeneric("readGappedReadsFromBam", signature="file",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
        standardGeneric("readGappedReadsFromBam")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Methods for BamFile objects.
###

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

.load_bamcols_from_bamfile <- function(bamfile, param, what0, ...,
                                       with.which_label=FALSE)
{
    flag0 <- scanBamFlag(isUnmappedQuery=FALSE)
    param <- .normargParam(param, flag0, what0)
    res <- scanBam(bamfile, ..., param=param)
    if (length(res) == 0L)  # should never happen
        stop("scanBam() returned a list of length zero")
    Rsamtools:::.load_bamcols_from_scanBam_res(res, param,
                    with.which_label=with.which_label)
}

.load_seqlengths_from_bamfile <- function(file, seqlevels)
{
    seqlengths <- scanBamHeader(file)[["targets"]]
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

.matesFromBam <- function(file, use.names, param, what0, with.which_label)
{
    bamcols <- .load_bamcols_from_bamfile(file, param, what0,
                                          with.which_label=with.which_label)
    seqlengths <- .load_seqlengths_from_bamfile(file, levels(bamcols$rname))
    gal <- GAlignments(seqnames=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       seqlengths=seqlengths)
    gal <- .bindExtraData(gal, use.names=FALSE, param, bamcols,
                          with.which_label=with.which_label)
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

setMethod("readGAlignmentsFromBam", "BamFile",
    function(file, index=file, ..., use.names=FALSE, param=NULL,
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
        bamcols <- .load_bamcols_from_bamfile(file, param, what0, ...,
                                              with.which_label=with.which_label)
        seqlengths <- .load_seqlengths_from_bamfile(file,
                                                    levels(bamcols[["rname"]]))
        ans <- GAlignments(seqnames=bamcols$rname, pos=bamcols$pos,
                           cigar=bamcols$cigar, strand=bamcols$strand,
                           seqlengths=seqlengths)
        .bindExtraData(ans, use.names, param, bamcols,
                       with.which_label=with.which_label)
    }
)

### 'use.mcols' can be TRUE, FALSE, or a character vector specifying the
### *inner* metadata columns to return, i.e., the metadata columns to set on
### the 2 halves of the returned GAlignmentPairs object.
.make_GAlignmentPairs_from_GAlignments <- function(gal, use.mcols=FALSE)
{
    mate_status <- mcols(gal)[ , "mate_status"]

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

    ## Keep alignments with "mated" mate status only.
    is_mated <- mate_status == "mated"
    if (!all(is_mated)) {
        keep_idx <- which(is_mated)
        gal <- gal[keep_idx]
    }

    ## Check flag bits 0x40 and 0x80.
    flag <- mcols(gal)[ , "flag"]
    is_first_mate <- bamFlagAsBitMatrix(flag, bitnames="isFirstMateRead")
    is_last_mate <- bamFlagAsBitMatrix(flag, bitnames="isSecondMateRead")
    bits_0x40_0x80_are_ok <- is_first_mate != is_last_mate
    if (!all(bits_0x40_0x80_are_ok)) {
        keep_idx <- which(bits_0x40_0x80_are_ok)
        gal <- gal[keep_idx]
        is_first_mate <- is_first_mate[keep_idx]
        is_last_mate <- is_last_mate[keep_idx]
    }

    ## Split and order the pairs by ascending start position of the first mate.
    idx1 <- which(as.logical(is_first_mate))
    oo1 <- S4Vectors:::orderIntegerPairs(as.integer(gal@seqnames)[idx1],
                                         gal@start[idx1])
    idx1 <- idx1[oo1]
    idx2 <- which(as.logical(is_last_mate))[oo1]
    ans_first <- gal[idx1]
    ans_last <- gal[idx2]
    groupid1 <- mcols(ans_first)[ , "groupid"]
    groupid2 <- mcols(ans_last)[ , "groupid"]
    stopifnot(identical(groupid1, groupid2))

    ## Drop the names.
    ans_names <- names(ans_first)
    names(ans_first) <- names(ans_last) <- NULL

    ## Check isProperPair (0x2) and isNotPrimaryRead (0x100) flag bits.
    flag1 <- mcols(ans_first)[ , "flag"]
    flag2 <- mcols(ans_last)[ , "flag"]
    is_proper1 <- bamFlagAsBitMatrix(flag1, bitnames="isProperPair")
    is_proper2 <- bamFlagAsBitMatrix(flag2, bitnames="isProperPair")
    stopifnot(identical(is_proper1, is_proper2))
    is_secondary1 <- bamFlagAsBitMatrix(flag1, bitnames="isNotPrimaryRead")
    is_secondary2 <- bamFlagAsBitMatrix(flag2, bitnames="isNotPrimaryRead")
    stopifnot(identical(is_secondary1, is_secondary2))

    ## Drop discordant pairs. 
    is_discordant <- (seqnames(ans_first) != seqnames(ans_last)) |
                     (strand(ans_first) == strand(ans_last))
    discordant_idx <- which(is_discordant)
    if (length(discordant_idx) != 0L) {
        nb_discordant_proper <- sum(is_proper1[discordant_idx])
        nb_discordant_not_proper <- length(discordant_idx) -
                                    nb_discordant_proper
        warning(length(discordant_idx), " pairs (", nb_discordant_proper,
                " proper, ", nb_discordant_not_proper, " not proper) were ",
                "dropped because the seqname\n  or strand of the alignments ",
                "in the pair were not concordant.\n",
                "  Note that a GAlignmentPairs object can only hold ",
                "concordant pairs at the\n  moment, that is, pairs where ",
                "the 2 alignments are on the opposite strands\n  of the same ",
                "chromosome.")
        keep_idx <- which(!is_discordant)
        ans_first <- ans_first[keep_idx]
        ans_last <- ans_last[keep_idx]
        is_proper1 <- is_proper1[keep_idx]
        ans_names <- ans_names[keep_idx]
    }

    ## Make the GAlignmentPairs object and return it.
    if (is.character(use.mcols)) {
        mcols(ans_first) <- mcols(ans_first)[use.mcols]
        mcols(ans_last) <- mcols(ans_last)[use.mcols]
    } else if (!use.mcols) {
        mcols(ans_first) <- mcols(ans_last) <- NULL
    }
    GAlignmentPairs(ans_first, ans_last, as.logical(is_proper1),
                    names=ans_names)
}

setMethod("readGAlignmentPairsFromBam", "BamFile",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
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
        gal <- readGAlignmentsFromBam(file, use.names=use.names,
                                      param=param2,
                                      with.which_label=with.which_label)
        use.mcols <- c(bamWhat(param), bamTag(param))
        if (with.which_label)
            use.mcols <- c(use.mcols, "which_label")
        .make_GAlignmentPairs_from_GAlignments(gal, use.mcols=use.mcols)
    }
)

setMethod("readGAlignmentsListFromBam", "BamFile",
    function(file, index=file, ..., use.names=FALSE, param=ScanBamParam(),
                   with.which_label=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!asMates(file))
            bamWhat(param) <- setdiff(bamWhat(param), 
                                      c("groupid", "mate_status"))
        what0 <- c("rname", "strand", "pos", "cigar", "groupid", "mate_status")
        if (use.names)
            what0 <- c(what0, "qname")
        .matesFromBam(file, use.names, param, what0, with.which_label)
    }
)

setMethod("readGappedReadsFromBam", "BamFile",
    function(file, index=file, use.names=FALSE, param=NULL,
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
        bamcols <- .load_bamcols_from_bamfile(file, param, what0,
                                              with.which_label=with.which_label)
        seqlengths <- .load_seqlengths_from_bamfile(file,
                                                    levels(bamcols[["rname"]]))
        ans <- GappedReads(seqnames=bamcols$rname, pos=bamcols$pos,
                           cigar=bamcols$cigar, strand=bamcols$strand,
                           qseq=bamcols$seq, seqlengths=seqlengths)
        .bindExtraData(ans, use.names, param, bamcols,
                       with.which_label=with.which_label)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Methods for character vectors.
###

setMethod("readGAlignmentsFromBam", "character",
    function(file, index=file, ..., use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
    {
        if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
            index <- character(0)
        bam <- open(BamFile(file, index), "rb")
        on.exit(close(bam))
        readGAlignmentsFromBam(bam, character(), ..., use.names=use.names,
                               param=param,
                               with.which_label=with.which_label)
    }
)

setMethod("readGappedReadsFromBam", "character",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
    {
        if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
            index <- character(0)
        bam <- open(BamFile(file, index), "rb")
        on.exit(close(bam))
        readGappedReadsFromBam(bam, character(), use.names=use.names,
                               param=param,
                               with.which_label=with.which_label)
    }
)

setMethod("readGAlignmentPairsFromBam", "character",
    function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
    {
        if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
            index <- character(0)
        bam <- open(BamFile(file, index, asMates=TRUE), "rb")
        on.exit(close(bam))
        readGAlignmentPairsFromBam(bam, character(), use.names=use.names,
                                   param=param,
                                   with.which_label=with.which_label)
    }
)

setMethod("readGAlignmentsListFromBam", "character",
    function(file, index=file, ..., use.names=FALSE,
                   param=ScanBamParam(), with.which_label=FALSE)
    {
        if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
            index <- character(0)
        bam <- open(BamFile(file, index, asMates=TRUE), "rb")
        on.exit(close(bam))
        readGAlignmentsListFromBam(bam, character(), use.names=use.names,
                                   param=param,
                                   with.which_label=with.which_label)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Methods for BamViews objects.
###

setMethod("readGAlignmentsFromBam", "BamViews",
    function(file, index=file, ..., use.names=FALSE, param=NULL,
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
            readGAlignmentsFromBam(file=bamPaths(bamViews)[i],
                                   index=bamIndicies(bamViews)[i],
                                   use.names=use.names,
                                   param=param,
                                   with.which_label=with.which_label)
        ### Rsamtools:::.BamViews_delegate requires the ShortRead package!
        Rsamtools:::.BamViews_delegate("readGAlignmentsFromBam", file, fun)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Front-ends
###

readGAlignments <- function(file, format="BAM", use.names=FALSE, ...)
{
    if (!isSingleString(format))
        stop("'format' must be a single string")
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (format == "BAM") {
        ans <- readGAlignmentsFromBam(file=file, use.names=use.names, ...)
        return(ans)
    }
    stop("only BAM format is supported at the moment")
}

readGAlignmentPairs <- function(file, format="BAM", use.names=FALSE, ...)
{
    if (!isSingleString(format))
        stop("'format' must be a single string")
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (format == "BAM") {
        ans <- readGAlignmentPairsFromBam(file=file, use.names=use.names, ...)
        return(ans)
    }
    stop("only BAM format is supported at the moment")
}

readGAlignmentsList <- function(file, format="BAM", use.names=FALSE, ...)
{
    if (!isSingleString(format))
        stop("'format' must be a single string")
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (format == "BAM") {
        ans <- readGAlignmentsListFromBam(file=file, use.names=use.names, ...)
        return(ans)
    }
    stop("only BAM format is supported at the moment")
}

readGappedReads <- function(file, format="BAM", use.names=FALSE, ...)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if (!isSingleString(format))
        stop("'format' must be a single string")
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (format == "BAM") {
        ans <- readGappedReadsFromBam(file=file, use.names=use.names, ...)
        return(ans)
    }
    stop("only BAM format is supported at the moment")
}

