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

setMethod("readGAlignmentPairsFromBam", "BamFile",
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
        if (!is.na(yieldSize(file))) {
            warning("'yieldSize' set to 'NA'", immediate.=TRUE)
            yieldSize(file) <- NA_integer_
        }
        flag0 <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE)
        what0 <- c("flag", "mrnm", "mpos")
        param2 <- .normargParam(param, flag0, what0)
        galn <- readGAlignmentsFromBam(file, use.names=TRUE, param=param2,
                                       with.which_label=with.which_label)
        if (is.null(param)) {
            use.mcols <- FALSE
        } else {
            use.mcols <- c(bamWhat(param), bamTag(param))
            if (with.which_label)
                use.mcols <- c(use.mcols, "which_label")
        }
        makeGAlignmentPairs(galn, use.names=use.names, use.mcols=use.mcols)
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
        bam <- open(BamFile(file, index), "rb")
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Old stuff.
###

readBamGappedAlignments <- function(...)
    .Defunct("readGAlignmentsFromBam")

readBamGappedReads <- function(...)
    .Defunct("readGappedReadsFromBam")

readBamGappedAlignmentPairs <- function(...)
    .Defunct("readGAlignmentPairsFromBam")

readBamGAlignmentsList <- function(...)
    .Defunct("readGAlignmentsListFromBam")

