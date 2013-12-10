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
.make_GAlignmentPairs_from_GAlignmentsList <- function(galist, use.mcols=FALSE)
{
    ## Dump "ambiguous" groups.
    flushDumpedAlignments()
    dumped_gal <- unlist(galist[mcols(galist)$mate_status == "ambiguous"])
    dumped_count <- length(dumped_gal)
    if (dumped_count != 0L) {
        dumpAlignments(dumped_gal)
        warning("  ", dumped_count, " alignments with ambiguous pairing ",
                "were dumped.\n    Use 'getDumpedAlignments()' to retrieve ",
                "them from the dump environment.")
    }

    ## Only keep "mated" groups.
    mate_status <- mcols(galist)[ , "mate_status"]
    galist <- galist[which(mate_status %in% "mated")]

    unlisted_galist <- unlist(galist, use.names=FALSE)

    ## Extract indices into 'unlisted_galist' of the 1st and 2nd elements
    ## of each pair.
    galist_partitioning <- PartitioningByEnd(galist)
    is_a_pair <- width(galist_partitioning) == 2L
    mate2_idx <- end(galist_partitioning)[is_a_pair]
    mate1_idx <- mate2_idx - 1L

    ## Check "flag" metadata column.
    flag <- mcols(unlisted_galist)[ , "flag"]
    flag1 <- flag[mate1_idx]
    flag2 <- flag[mate2_idx]

    is_first_mate1 <- bamFlagAsBitMatrix(flag1, bitnames="isFirstMateRead")
    is_last_mate1 <- bamFlagAsBitMatrix(flag1, bitnames="isSecondMateRead")
    is_first_mate2 <- bamFlagAsBitMatrix(flag2, bitnames="isFirstMateRead")
    is_last_mate2 <- bamFlagAsBitMatrix(flag2, bitnames="isSecondMateRead")
    stopifnot(all(is_first_mate1))
    stopifnot(all(is_first_mate1 != is_last_mate1))
    stopifnot(all(is_first_mate2 != is_last_mate2))
    stopifnot(all(is_first_mate1 == is_last_mate2))
    #switch_mates_idx <- which(is_last_mate1 != 0L)
    #mate1_idx[switch_mates_idx] <- mate1_idx[switch_mates_idx] + 1L
    #mate2_idx[switch_mates_idx] <- mate2_idx[switch_mates_idx] - 1L

    is_proper1 <- bamFlagAsBitMatrix(flag1, bitnames="isProperPair")
    is_proper2 <- bamFlagAsBitMatrix(flag2, bitnames="isProperPair")
    stopifnot(identical(is_proper1, is_proper2))

    is_secondary1 <- bamFlagAsBitMatrix(flag1, bitnames="isNotPrimaryRead")
    is_secondary2 <- bamFlagAsBitMatrix(flag2, bitnames="isNotPrimaryRead")
    stopifnot(identical(is_secondary1, is_secondary2))

    ## Split 'unlisted_galist' in 2 parallel GAlignments objects: 'ans_first'
    ## and 'ans_last'.
    ans_first <- unlisted_galist[mate1_idx]
    ans_last <- unlisted_galist[mate2_idx]
    ans_names <- names(galist)[is_a_pair]

    ## Drop pairs with discordant seqnames or strand. Right now the pairs in
    ## a GAlignmentPairs object are assumed to be concordant but maybe this
    ## should be revisited.
    is_discordant <- (as.character(seqnames(ans_first)) !=
                      as.character(seqnames(ans_last))) |
                     (as.character(strand(ans_first)) ==
                      as.character(strand(ans_last)))
    discordant_idx <- which(is_discordant)
    if (length(discordant_idx) != 0L) {
        nb_discordant_proper <- sum(is_proper1[discordant_idx])
        if (nb_discordant_proper != 0L) {
            ratio <- 100.0 * nb_discordant_proper / length(discordant_idx)
            warning(ratio, "% of the pairs with discordant seqnames or ",
                    "strand were flagged\n",
                    "  as proper pairs by the aligner. Dropping them anyway.")
        }
        concordant_idx <- which(!is_discordant)
        ans_first <- ans_first[concordant_idx]
        ans_last <- ans_last[concordant_idx]
        ans_names <- ans_names[concordant_idx]
        is_proper1 <- is_proper1[concordant_idx]
    }

    ## Order the pairs by ascending start position of the first mate.
    oo <- IRanges:::orderIntegerPairs(as.integer(ans_first@seqnames),
                                      ans_first@start)
    ans_first <- ans_first[oo]
    ans_last <- ans_last[oo]
    ans_names <- ans_names[oo]
    is_proper1 <- is_proper1[oo]

    ## Make the GAlignmentPairs object and return it.
    names(ans_first) <- names(ans_last) <- NULL
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
        what0 <- "flag"
        param2 <- .normargParam(param, flag0, what0)
        galist <- readGAlignmentsListFromBam(file, use.names=use.names,
                                             param=param2,
                                             with.which_label=with.which_label)
        use.mcols <- c(bamWhat(param), bamTag(param))
        if (with.which_label)
            use.mcols <- c(use.mcols, "which_label")
        .make_GAlignmentPairs_from_GAlignmentsList(galist, use.mcols=use.mcols)
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
