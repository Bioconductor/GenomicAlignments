### =========================================================================
### mapToAlignment() and pmapToAlignment()
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics 
###

setGeneric("mapToAlignment", signature=c("x", "alignment"),
    function(x, alignment, reverse=FALSE, ...) 
        standardGeneric("mapToAlignment")
)

setGeneric("pmapToAlignment", signature=c("x", "alignment"),
    function(x, alignment, reverse=FALSE, ...) 
        standardGeneric("pmapToAlignment")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToAlignment() methods
###

.mapToAlignment <- function(x, alignment, reverse) 
{
    if (length(x) && length(alignment)) {
        if (reverse)
            FUN <- "map_to_transcript"
        else
            FUN <- "map_to_genome"
        map <- .Call(FUN, start(x), end(x), cigar(alignment), start(alignment))
        starts <- map[[1]]
        if (!all(length(starts))) {
            ans <- GRanges()
            mcols(ans) <- DataFrame(xHits=integer(), alignmentHits=integer())
            return(ans) 
        }
        ends <- pmax(map[[2]], starts - 1L)
        xHits <- map[[3]]
        alignmentHits <- map[[4]]
        seqname <- as.character(seqnames(alignment)[alignmentHits])
        if (any(skip <- is.na(starts) | is.na(ends))) {
            starts[skip] <- 1L 
            ends[skip] <- 0L
            seqname[skip] <- "unmapped"
        }
        ## result is "*" strand
        GRanges(Rle(seqname), IRanges(starts, ends), strand="*", 
                DataFrame(xHits, alignmentHits))
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), alignmentHits=integer())
        ans
    }
}

setMethod("mapToAlignment", c("Ranges", "GAlignments"),
    function(x, alignment, reverse=FALSE, ...)
    {
        if (reverse)
            ranges(.mapToAlignment(x, alignment, TRUE))
        else
            ranges(.mapToAlignment(x, alignment, FALSE))
    }
)

setMethod("mapToAlignment", c("GenomicRanges", "GAlignments"),
    function(x, alignment, reverse=FALSE, ...)
    {
        if (reverse)
            .mapToAlignment(x, alignment, TRUE)
        else
            .mapToAlignment(x, alignment, FALSE)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmapToAlignment() methods
###

.pmapToAlignment <- function(x, alignment, reverse)
{
    if (length(x) && length(alignment)) {
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")

        if (reverse)
            FUN <- "ref_locs_to_query_locs"
        else
            FUN <- "query_locs_to_ref_locs"
        starts <- .Call(FUN, start(x), cigar(alignment), 
                        start(alignment), FALSE)
        ends <- .Call(FUN, end(x), cigar(alignment), 
                      start(alignment), TRUE)
        ends <- pmax(ends, starts - 1L)
        seqname <- as.character(seqnames(alignment))
        if (any(skip <- is.na(starts) | is.na(ends))) {
            starts[skip] <- 1L 
            ends[skip] <- 0L
            seqname[skip] <- "unmapped"
        }
        ## result is "*" strand
        GRanges(Rle(seqname), IRanges(starts, ends))
    } else {
        GRanges()
    }
}

setMethod("pmapToAlignment", c("Ranges", "GAlignments"),
    function(x, alignment, reverse=FALSE, ...)
    {
        if (reverse)
            ranges(.pmapToAlignment(x, alignment, TRUE))
        else
            ranges(.pmapToAlignment(x, alignment, FALSE))
    }
)

setMethod("pmapToAlignment", c("GenomicRanges", "GAlignments"), 
    function(x, alignment, reverse=FALSE, ...) 
    {
        if (reverse)
            .pmapToAlignment(ranges(x), alignment, TRUE)
        else
            .pmapToAlignment(ranges(x), alignment, FALSE)
    }
)

