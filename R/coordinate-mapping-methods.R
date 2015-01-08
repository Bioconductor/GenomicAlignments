### =========================================================================
### mapToAlignments() and pmapToAlignments()
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics 
###

setGeneric("mapToAlignments", signature=c("x", "alignments"),
    function(x, alignments, reverse=FALSE, ...) 
        standardGeneric("mapToAlignments")
)

setGeneric("pmapToAlignments", signature=c("x", "alignments"),
    function(x, alignments, reverse=FALSE, ...) 
        standardGeneric("pmapToAlignments")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToAlignments() methods
###

.mapToAlignments <- function(x, alignments, reverse) 
{
    if (length(x) && length(alignments)) {
        if (reverse)
            FUN <- "map_query_locs_to_ref_locs"
        else
            FUN <- "map_ref_locs_to_query_locs"
        map <- .Call(FUN, start(x), end(x), cigar(alignments), start(alignments))
        starts <- map[[1]]
        if (!all(length(starts))) {
            ans <- GRanges()
            mcols(ans) <- DataFrame(xHits=integer(), alignmentsHits=integer())
            return(ans) 
        }
        ends <- pmax(map[[2]], starts - 1L)
        xHits <- map[[3]]
        alignmentsHits <- map[[4]]
        seqname <- as.character(seqnames(alignments)[alignmentsHits])
        if (any(skip <- is.na(starts) | is.na(ends))) {
            starts[skip] <- 1L 
            ends[skip] <- 0L
            seqname[skip] <- "unmapped"
        }
        ## result is "*" strand
        GRanges(Rle(seqname), IRanges(starts, ends), strand="*", 
                DataFrame(xHits, alignmentsHits))
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), alignmentsHits=integer())
        ans
    }
}

setMethod("mapToAlignments", c("Ranges", "GAlignments"),
    function(x, alignments, reverse=FALSE, ...)
    {
        if (reverse)
            ranges(.mapToAlignments(x, alignments, TRUE))
        else
            ranges(.mapToAlignments(x, alignments, FALSE))
    }
)

setMethod("mapToAlignments", c("GenomicRanges", "GAlignments"),
    function(x, alignments, reverse=FALSE, ...)
    {
        if (reverse)
            .mapToAlignments(x, alignments, TRUE)
        else
            .mapToAlignments(x, alignments, FALSE)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmapToAlignments() methods
###

.pmapToAlignments <- function(x, alignments, reverse)
{
    if (length(x) && length(alignments)) {
        if (length(x) != length(alignments))
            stop("'x' and 'alignments' must have the same length")

        if (reverse)
            FUN <- "query_locs_to_ref_locs"
        else
            FUN <- "ref_locs_to_query_locs"
        s <- .Call(FUN, start(x), cigar(alignments), start(alignments), FALSE)
        e <- .Call(FUN, end(x), cigar(alignments), start(alignments), TRUE)
        e <- pmax(e, s - 1L)

        ## non-hits are zero width with seqname "unmmapped"
        seqname <- as.character(seqnames(alignments))
        if (any(skip <- is.na(s) | is.na(e))) {
            s[skip] <- 1L 
            e[skip] <- 0L
            seqname[skip] <- "unmapped"
        }
        ## result is "*" strand
        GRanges(Rle(seqname), IRanges(s, e))
    } else {
        GRanges()
    }
}

setMethod("pmapToAlignments", c("Ranges", "GAlignments"),
    function(x, alignments, reverse=FALSE, ...)
    {
        if (reverse)
            ranges(.pmapToAlignments(x, alignments, TRUE))
        else
            ranges(.pmapToAlignments(x, alignments, FALSE))
    }
)

setMethod("pmapToAlignments", c("GenomicRanges", "GAlignments"), 
    function(x, alignments, reverse=FALSE, ...) 
    {
        if (reverse)
            .pmapToAlignments(ranges(x), alignments, TRUE)
        else
            .pmapToAlignments(ranges(x), alignments, FALSE)

    }
)

