### =========================================================================
### "mapCoords" methods
### -------------------------------------------------------------------------
###

setMethod("map", c("GenomicRanges", "GAlignments"), function(from, to) {
  .Defunct(msg="map() is defunct. Use mapCoords() instead.")
})

setMethod("pmap", c("Ranges", "GAlignments"), function(from, to) {
  .Defunct(msg="pmap() is defunct. Use pmapCoords() instead.")
})

setMethod("mapCoords", c("GenomicRanges", "GAlignments"), function(x, to, ...) {
  to_grl <- grglist(to, drop.D.ranges=TRUE)
  from_ol <- findOverlaps(x, to_grl, ignore.strand=TRUE, type="within")
  to_hits <- to[subjectHits(from_ol)]
  from_hits <- ranges(x)[queryHits(from_ol)]
  ranges <- pmapCoords(from_hits, to_hits)
  space <- names(to_hits)
  if (is.null(space))
    space <- as.character(seq_len(length(to))[subjectHits(from_ol)])
  GRanges(Rle(space), ranges, queryHits=queryHits(from_ol),
          subjectHits=subjectHits(from_ol))
})

setMethod("pmapCoords", c("Ranges", "GAlignments"), function(x, to, ...) {
  starts <- .Call("ref_locs_to_query_locs", start(x),
                  cigar(to), start(to), FALSE, PACKAGE="GenomicAlignments")
  ends <- .Call("ref_locs_to_query_locs", end(x),
                cigar(to), start(to), TRUE, PACKAGE="GenomicAlignments")
  ends <- pmax(ends, starts - 1L)
  IRanges(starts, ends)
})

setGeneric("prmap", function(from, to) standardGeneric("prmap")) # not exported

setMethod("prmap", c("Ranges", "GAlignments"), # not exported
  function(from, to) {
  starts <- .Call("query_locs_to_ref_locs", start(from),
                  cigar(to), start(to), FALSE, PACKAGE="GenomicAlignments")
  ends <- .Call("query_locs_to_ref_locs", end(from),
                cigar(to), start(to), TRUE, PACKAGE="GenomicAlignments")
  ends <- pmax(ends, starts - 1L)
  IRanges(starts, ends)
})

