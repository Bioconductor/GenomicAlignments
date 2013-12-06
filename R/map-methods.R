### =========================================================================
### "map" methods
### -------------------------------------------------------------------------
###


setMethod("map", c("GenomicRanges", "GAlignments"), function(from, to) {
  to_grl <- grglist(to, drop.D.ranges=TRUE)
  from_ol <- findOverlaps(from, to_grl, ignore.strand=TRUE, type="within")
  to_hits <- to[subjectHits(from_ol)]
  from_hits <- ranges(from)[queryHits(from_ol)]
  ranges <- pmap(from_hits, to_hits)
  space <- names(to_hits)
  if (is.null(space))
    space <- as.character(seq_len(length(to))[subjectHits(from_ol)])
  new("RangesMapping", hits = from_ol, space = Rle(space),
      ranges = ranges)
})

setMethod("pmap", c("Ranges", "GAlignments"), function(from, to) {
  starts <- .Call("ref_locs_to_query_locs", start(from),
                  cigar(to), start(to), FALSE, PACKAGE="GenomicRanges")
  ends <- .Call("ref_locs_to_query_locs", end(from),
                cigar(to), start(to), TRUE, PACKAGE="GenomicRanges")
  ends <- pmax(ends, starts - 1L)
  IRanges(starts, ends)
})

setGeneric("prmap", function(from, to) standardGeneric("prmap")) # not exported

setMethod("prmap", c("Ranges", "GAlignments"), # not exported
  function(from, to) {
  starts <- .Call("query_locs_to_ref_locs", start(from),
                  cigar(to), start(to), FALSE, PACKAGE="GenomicRanges")
  ends <- .Call("query_locs_to_ref_locs", end(from),
                cigar(to), start(to), TRUE, PACKAGE="GenomicRanges")
  ends <- pmax(ends, starts - 1L)
  IRanges(starts, ends)
})

