### =========================================================================
### 'mapCoords' and 'pmapCoords' methods
### -------------------------------------------------------------------------
###

### Generics are in IRanges.

### mapCoords:

setMethod("mapCoords", c("GenomicRanges", "GAlignments"), 
    function(from, to, ...) 
    {
        msg <- c("mapCoords() is deprecated. ",
                 "Use 'mapToTranscripts' from the GenomicFeatures package ",
                 "or 'mapToAlignments' from the GenomicAlignments package ",
                 "instead.")
        .Deprecated(msg=wmsg(msg))
        to_grl <- grglist(to, drop.D.ranges=TRUE)
        from_ol <- findOverlaps(from, to_grl, ignore.strand=TRUE, type="within")
        to_hits <- to[subjectHits(from_ol)]
        from_hits <- ranges(from)[queryHits(from_ol)]
        ranges <- pmapCoords(from_hits, to_hits)
        space <- names(to_hits)
        if (is.null(space))
          space <- as.character(seq_len(length(to))[subjectHits(from_ol)])
 
        GRanges(Rle(space), ranges, fromHits=queryHits(from_ol),
                toHits=subjectHits(from_ol))
    }
)

### pmapCoords:

setMethod("pmapCoords", c("Ranges", "GAlignments"), 
    function(from, to, ...) 
    {
        msg <- c("'pmapCoords' is deprecated. ",
                 "Use 'pmapToTranscripts' from the GenomicFeatures package ",
                 "or 'pmapToAlignments' from the GenomicAlignments package ",
                 "instead.")
        .Deprecated(msg=wmsg(msg))
        starts <- .Call("ref_locs_to_query_locs", start(from), cigar(to), 
                        start(to), FALSE, PACKAGE="GenomicAlignments")
        ends <- .Call("ref_locs_to_query_locs", end(from), cigar(to), 
                      start(to), TRUE, PACKAGE="GenomicAlignments")
        ends <- pmax(ends, starts - 1L)
        IRanges(starts, ends)
    }
)

### prmap (not exported):

setGeneric("prmap", function(from, to) standardGeneric("prmap"))

setMethod("prmap", c("Ranges", "GAlignments"),
    function(from, to) 
    {
        .Deprecated("pmapFromTranscripts", old="prmap")
        starts <- .Call("query_locs_to_ref_locs", start(from), cigar(to), 
                        start(to), FALSE, PACKAGE="GenomicAlignments")
        ends <- .Call("query_locs_to_ref_locs", end(from), cigar(to), 
                      start(to), TRUE, PACKAGE="GenomicAlignments")
        ends <- pmax(ends, starts - 1L)
        IRanges(starts, ends)
    }
)
