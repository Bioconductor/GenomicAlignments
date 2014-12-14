### =========================================================================
### [p]mapToGenome() and [p]mapToTranscript() methods
### -------------------------------------------------------------------------
###

### Generics are in GenomicRanges.

### Non-parallel: mapToGenome() and mapToTranscript()
### - All elements in 'x' are mapped to 'alignment'
### - Result length varies like a Hits object
### - Result only contains mapped records; non-hits are not returned
### FIXME: currently not strand-aware

### Parallel: pmapToGenome() and pmapToTranscript
### - i-th element of 'x' is mapped to i-th element of 'alignment'
### - Result is the same length as 'x'
### - Ranges with strand mismatch or no hit are returned as zero-width ranges.
###
###   When mapping transcript -> genome the range in 'x' can never
###   be before the range in 'alignment' and therefore the start of the 
###   zero-width range begins after the width of 'alignment'.
###
###   When mapping genome -> transcript the range in 'x' can be either
###   before or after the range in 'alignment'. Start of the zero-width range
###   begins at 0 if if falls before the range or after the width of 
###   'alignment' if it falls after. Strand mismatch zero-width ranges start 
###   at 0.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToGenome and pmapToGenome 
###

.mapToGenome <- function(x, alignment)
{
    map <- .Call("map_to_genome", 
                 start(x), end(x), cigar(alignment), start(alignment),
                 PACKAGE="GenomicAlignments")
    starts <- map[[1]]
    ends <- pmax(map[[2]], starts - 1L)

    ans <- IRanges(starts, ends)
    mcols(ans) <- DataFrame(x_hits=map[[3]], alignment_hits=map[[4]])
    ans
}

setMethod("mapToGenome", c("Ranges", "GAlignments"), .mapToGenome)

setMethod("mapToGenome", c("GRanges", "GAlignments"),
    function(x, alignment, ...)
    {
        map <- .mapToGenome(x, alignment)
        index <- mcols(map)$x_hits
        ans <- suppressWarnings(GRanges(seqnames(x)[index], map, 
                                strand(x)[index]))
        mcols(ans) <- mcols(map)
        ans
    }
)

.pmapToGenome <- function(x, alignment, strand_mismatch)
{
    starts <- .Call("query_locs_to_ref_locs", 
                    start(x), cigar(alignment), start(alignment), 
                    FALSE, PACKAGE="GenomicAlignments")
    ends <- .Call("query_locs_to_ref_locs", 
                  end(x), cigar(alignment), start(alignment), 
                  TRUE, PACKAGE="GenomicAlignments")
    ends <- pmax(ends, starts - 1L)

    skip <- is.na(starts) | is.na(ends) | strand_mismatch
    if (any(skip)) {
        starts[skip] <- width(alignment[skip]) + 1L 
        ends[skip] <- starts[skip] - 1L
    }
    IRanges(starts, ends)
}

setMethod("pmapToGenome", c("Ranges", "GAlignments"),
    function(x, alignment, ...)
    {
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")
        .pmapToGenome(x, alignment, logical(length(x)))
    }
)

setMethod("pmapToGenome", c("GRanges", "GAlignments"), 
    function(x, alignment, ignore.strand = TRUE, ...) 
    {
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")
        if (ignore.strand)
            strand_mismatch <- logical(length(x)) 
        else
            strand_mismatch <- as.logical(strand(x) != strand(alignment))

        map <- .pmapToGenome(ranges(x), alignment, strand_mismatch)
        GRanges(seqnames(x), map, strand(x))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToTranscript and pmapToTranscript 
###

.mapToTranscript <- function(x, alignment)
{
    map <- .Call("map_to_transcript", 
                 start(x), end(x), cigar(alignment), start(alignment),
                 PACKAGE="GenomicAlignments")
    starts <- map[[1]]
    ends <- pmax(map[[2]], starts - 1L)

    ans <- IRanges(starts, ends)
    mcols(ans) <- DataFrame(x_hits=map[[3]], alignment_hits=map[[4]])
    ans
}

setMethod("mapToTranscript", c("Ranges", "GAlignments"), .mapToTranscript)

setMethod("mapToTranscript", c("GRanges", "GAlignments"),
    function(x, alignment, ...)
    {
        map <- .mapToTranscript(x, alignment)
        index <- mcols(map)$x_hits
        ans <- suppressWarnings(GRanges(seqnames(x)[index], map, 
                                strand(x)[index]))
        mcols(ans) <- mcols(map)
        ans
    }
)

.pmapToTranscript <- function(x, alignment, strand_mismatch)
{
    starts <- .Call("ref_locs_to_query_locs", 
                    start(x), cigar(alignment), start(alignment), 
                    FALSE, PACKAGE="GenomicAlignments")
    ends <- .Call("ref_locs_to_query_locs", 
                  end(x), cigar(alignment), start(alignment), 
                  TRUE, PACKAGE="GenomicAlignments")
    ends <- pmax(ends, starts - 1L)

    before <- start(x) < start(alignment)
    after <- end(x) > end(alignment)
    skip <- before | after | strand_mismatch
    if (any(skip)) {
        starts[before | strand_mismatch] <- 1L 
        ends[before | strand_mismatch] <- 0L 
        starts[after] <- width(alignment[after]) + 1L 
        ends[after] <- starts[after] - 1L
    }
    IRanges(starts, ends)
}

setMethod("pmapToTranscript", c("Ranges", "GAlignments"),
    function(x, alignment, ...)
    {
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")
        .pmapToTranscript(x, alignment, logical(length(x)))
    }
)

setMethod("pmapToTranscript", c("GRanges", "GAlignments"), 
    function(x, alignment, ignore.strand = TRUE, ...) 
    {
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")
        if (ignore.strand)
            strand_mismatch <- logical(length(x))
        else
            strand_mismatch <- as.logical(strand(x) != strand(alignment))

        map <- .pmapToTranscript(ranges(x), alignment, strand_mismatch)
        GRanges(seqnames(x), map, strand(x))
    }
)
