### =========================================================================
### coordinate mapping methods
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics 
###

setGeneric("mapToAlignments", signature=c("x", "alignments"),
    function(x, alignments, ...) 
        standardGeneric("mapToAlignments")
)

setGeneric("pmapToAlignments", signature=c("x", "alignments"),
    function(x, alignments, ...) 
        standardGeneric("pmapToAlignments")
)

setGeneric("mapFromAlignments", signature=c("x", "alignments"),
    function(x, alignments, ...) 
        standardGeneric("mapFromAlignments")
)

setGeneric("pmapFromAlignments", signature=c("x", "alignments"),
    function(x, alignments, ...) 
        standardGeneric("pmapFromAlignments")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToAlignments() and mapFromAlignments() methods
###

.mapFromAlignments <- function(x, alignments)
{
    if (!length(x) && !length(alignments))
        return(GRanges(xHits=integer(), transcriptsHits=integer()))
    if (is.null(xNames <- names(x)) || 
        is.null(alignmentsNames <- names(alignments)))
        stop ("both 'x' and 'alignments' must have names")

    ## name matching determines pairs
    match0 <- match(alignmentsNames, alignmentsNames)
    match1 <- match(xNames, alignmentsNames)
    group0 <- splitAsList(seq_along(alignmentsNames), match0)
    group1 <- group0[match(na.omit(match1), names(group0))]
    xHits <- rep(which(!is.na(match1)), elementNROWS(group1))
    alignmentsHits <- unlist(group1, use.names=FALSE)
    if (!length(xHits <- na.omit(xHits)))
        stop ("none of 'names(x)' are in 'names(alignments)'")

    x <- x[xHits]
    alignments <- alignments[alignmentsHits]
    s <- .Call("query_locs_to_ref_locs",
               start(x), cigar(alignments), 
               start(alignments), FALSE)
    e <- .Call("query_locs_to_ref_locs", 
               end(x), cigar(alignments), 
               start(alignments), TRUE)
    e <- pmax(e, s - 1L)

    ## remove non-hits
    keep <- !is.na(s) & !is.na(e)
    seqname <- as.character(seqnames(alignments))
    GRanges(Rle(seqname[keep]), 
            IRanges(s[keep], e[keep], names=names(x)[keep]), 
            xHits=xHits[keep], alignmentsHits=alignmentsHits[keep])
}

.mapToAlignments <- function(x, alignments) 
{
    if (!length(x) && !length(alignments))
        return(GRanges(xHits=integer(), transcriptsHits=integer()))
    if (is.null(names(alignments)))
        stop ("'alignments' must have names")

    ## map all possible pairs; returns hits only
    map <- .Call("map_ref_locs_to_query_locs", 
                 start(x), end(x), cigar(alignments), 
                 start(alignments))
    xHits <- map[[3]]
    alignmentsHits <-  map[[4]]
    if (length(xHits))
        GRanges(Rle(names(alignments)[alignmentsHits]), 
                IRanges(map[[1]], pmax(map[[2]], map[[1]] - 1L), 
                        names=names(x)[xHits]), 
                strand="*", xHits, alignmentsHits)
    else
        GRanges(xHits=integer(), transcriptsHits=integer())
}

setMethod("mapToAlignments", c("Ranges", "GAlignments"),
    function(x, alignments, ...)
        ranges(.mapToAlignments(x, alignments))
)

setMethod("mapToAlignments", c("GenomicRanges", "GAlignments"),
    function(x, alignments, ...)
        .mapToAlignments(x, alignments)
)

setMethod("mapFromAlignments", c("Ranges", "GAlignments"),
    function(x, alignments, ...)
        ranges(.mapFromAlignments(x, alignments))
)

setMethod("mapFromAlignments", c("GenomicRanges", "GAlignments"),
    function(x, alignments, ...) 
        .mapFromAlignments(x, alignments)
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmapToAlignments() and pmapFromAlignments() methods
###

.pmapAlignments <- function(x, alignments, reverse)
{
    if (length(x) && length(alignments)) {
        if (length(x) != length(alignments))
            stop("'x' and 'alignments' must have the same length")

        if (reverse) {
            FUN <- "query_locs_to_ref_locs"
            seqname <- as.character(seqnames(alignments))
        } else {
            if (is.null(names(alignments)))
                stop ("'alignments' must have names")
            FUN <- "ref_locs_to_query_locs"
            seqname <- names(alignments)
        }
        s <- .Call(FUN, start(x), cigar(alignments), start(alignments), FALSE)
        e <- .Call(FUN, end(x), cigar(alignments), start(alignments), TRUE)
        e <- pmax(e, s - 1L)

        ## non-hits
        if (any(skip <- is.na(s) | is.na(e))) {
            s[skip] <- 0L 
            e[skip] <- -1L
            seqname[skip] <- "UNMAPPED"
        }
        GRanges(Rle(seqname), IRanges(s, e, names=names(x)))
    } else {
        GRanges()
    }
}

setMethod("pmapToAlignments", c("Ranges", "GAlignments"),
    function(x, alignments, ...)
        ranges(.pmapAlignments(x, alignments, FALSE))
)

setMethod("pmapToAlignments", c("GenomicRanges", "GAlignments"), 
    function(x, alignments, ...) 
        .pmapAlignments(ranges(x), alignments, FALSE)
)

setMethod("pmapFromAlignments", c("Ranges", "GAlignments"),
    function(x, alignments, ...)
        ranges(.pmapAlignments(x, alignments, TRUE))
)

setMethod("pmapFromAlignments", c("GenomicRanges", "GAlignments"), 
    function(x, alignments, ...) 
        .pmapAlignments(ranges(x), alignments, TRUE)
)
