### =========================================================================
### [p]mapToGenome() and [p]mapToTranscript() methods
### -------------------------------------------------------------------------
###

### Generics are in GenomicRanges.

### Non-parallel: mapToGenome() and mapToTranscript()
### - All elements in 'from' are mapped to 'to'
### - Result length varies like a Hits object
### - Result only contains mapped records; non-hits are not returned
### FIXME: currently not strand-aware

### Parallel: pmapToGenome() and pmapToTranscript
### - i-th element of 'from' is mapped to i-th element of 'to'
### - Result is the same length as 'from'
### - Ranges with strand mismatch or no hit are returned as zero-width ranges.
###
###   When mapping transcript -> genome the range in 'from' can never
###   be before the range in 'to' and therefore the start of the zero-width
###   range begins after the width of 'to'.
###
###   When mapping genome -> transcript the range in 'from' can be either
###   before or after the range in 'to'. Start of the zero-width range
###   begins at 0 if if falls before the range or after the width of 'to' if 
###   it falls after.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToGenome and pmapToGenome 
###

.mapToGenome <- function(from, to)
{
    map <- .Call("map_to_genome", 
                 start(from), end(from), cigar(to), start(to),
                 PACKAGE="GenomicAlignments")
    starts <- map[[1]]
    ends <- pmax(map[[2]], starts - 1L)

    ans <- IRanges(starts, ends)
    mcols(ans) <- DataFrame(from_hits=map[[3]], to_hits=map[[4]])
    ans
}

setMethod("mapToGenome", c("Ranges", "GAlignments"), .mapToGenome)

setMethod("mapToGenome", c("GRanges", "GAlignments"),
    function(from, to, ...)
    {
        map <- .mapToGenome(from, to)
        index <- mcols(map)$from_hits
        ans <- suppressWarnings(GRanges(seqnames(from)[index], map, 
                                strand(from)[index]))
        mcols(ans) <- mcols(map)
        ans
    }
)

.pmapToGenome <- function(from, to, strand_mismatch)
{
    starts <- .Call("query_locs_to_ref_locs", 
                    start(from), cigar(to), start(to), 
                    FALSE, PACKAGE="GenomicAlignments")
    ends <- .Call("query_locs_to_ref_locs", 
                  end(from), cigar(to), start(to), 
                  TRUE, PACKAGE="GenomicAlignments")
    ends <- pmax(ends, starts - 1L)

    skip <- is.na(starts) | is.na(ends) | strand_mismatch
    if (any(skip)) {
        starts[skip] <- width(to[skip]) + 1L 
        ends[skip] <- starts[skip] - 1L
    }
    IRanges(starts, ends)
}

setMethod("pmapToGenome", c("Ranges", "GAlignments"),
    function(from, to, ...)
    {
        if (length(from) != length(to))
            stop("'from' and 'to' must have the same length")
        .pmapToGenome(from, to, logical(length(from)))
    }
)

setMethod("pmapToGenome", c("GRanges", "GAlignments"), 
    function(from, to, ignore.strand = TRUE, ...) 
    {
        if (length(from) != length(to))
            stop("'from' and 'to' must have the same length")
        if (ignore.strand)
            strand_mismatch <- logical(length(from)) 
        else
            strand_mismatch <- as.logical(strand(from) != strand(to))

        map <- .pmapToGenome(ranges(from), to, strand_mismatch)
        GRanges(seqnames(from), map, strand(from))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToTranscript and pmapToTranscript 
###

.mapToTranscript <- function(from, to)
{
    map <- .Call("map_to_transcript", 
                 start(from), end(from), cigar(to), start(to),
                 PACKAGE="GenomicAlignments")
    starts <- map[[1]]
    ends <- pmax(map[[2]], starts - 1L)

    ans <- IRanges(starts, ends)
    mcols(ans) <- DataFrame(from_hits=map[[3]], to_hits=map[[4]])
    ans
}

setMethod("mapToTranscript", c("Ranges", "GAlignments"), .mapToTranscript)

setMethod("mapToTranscript", c("GRanges", "GAlignments"),
    function(from, to, ...)
    {
        map <- .mapToTranscript(from, to)
        index <- mcols(map)$from_hits
        ans <- suppressWarnings(GRanges(seqnames(from)[index], map, 
                                strand(from)[index]))
        mcols(ans) <- mcols(map)
        ans
    }
)

.pmapToTranscript <- function(from, to, strand_mismatch)
{
    starts <- .Call("ref_locs_to_query_locs", 
                    start(from), cigar(to), start(to), 
                    FALSE, PACKAGE="GenomicAlignments")
    ends <- .Call("ref_locs_to_query_locs", 
                  end(from), cigar(to), start(to), 
                  TRUE, PACKAGE="GenomicAlignments")
    ends <- pmax(ends, starts - 1L)

    before <- start(from) < start(to)
    after <- end(from) > end(to)
    skip <- before | after | strand_mismatch
    if (any(skip)) {
        starts[before | strand_mismatch] <- 1L 
        ends[before | strand_mismatch] <- 0L 
        starts[after] <- width(to[after]) + 1L 
        ends[after] <- starts[after] - 1L
    }
    IRanges(starts, ends)
}

setMethod("pmapToTranscript", c("Ranges", "GAlignments"),
    function(from, to, ...)
    {
        if (length(from) != length(to))
            stop("'from' and 'to' must have the same length")
        .pmapToTranscript(from, to, logical(length(from)))
    }
)

setMethod("pmapToTranscript", c("GRanges", "GAlignments"), 
    function(from, to, ignore.strand = TRUE, ...) 
    {
        if (length(from) != length(to))
            stop("'from' and 'to' must have the same length")
        if (ignore.strand)
            strand_mismatch <- logical(length(from))
        else
            strand_mismatch <- as.logical(strand(from) != strand(to))

        map <- .pmapToTranscript(ranges(from), to, strand_mismatch)
        GRanges(seqnames(from), map, strand(from))
    }
)
