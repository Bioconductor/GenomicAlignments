## CIGAR ops M, =, X
GR_1 <- GRanges("chr1", IRanges(c(5, 10, 20, 25), width=2), strand="+")
cigar <- c("11M", "11=", "11X")
GA_1 <- GAlignments(rep("chr1", 3), rep(10L, 3), cigar, strand(rep("+", 3)))

## CIGAR ops S, N, D, I, H, P
GR_2 <- GRanges("chr1", IRanges(c(1, 20), width=6), strand="+")
cigar <- c("1S6M1S", "3M2N3M", "3M2D3M", "3M2I3M", "1H6M1H", "1P6M1P")
GA_2 <- GAlignments(rep("chr1", 6), rep(10L, 6), cigar, strand(rep("+", 6)))
GA_3 <- GAlignments(rep("chr1", 6), rep(20L, 6), cigar, strand(rep("+", 6)))

test_mapToGenome <- function() {
    ans <- mapToGenome(GR_1, GA_1)
    checkIdentical(start(ans), c(14L, 14L, 14L, 19L, 19L, 19L))
    checkIdentical(mcols(ans)$from_hits, c(rep(1L, 3), rep(2L, 3)))
    checkIdentical(mcols(ans)$to_hits, rep(1:3, 2))

    ans <- mapToGenome(GR_2, GA_2)
    checkIdentical(start(ans), rep(10L, 6))
    checkIdentical(end(ans), c(15L, 17L, 17L, 13L, 15L, 15L)) 
    checkIdentical(mcols(ans)$to_hits, as.integer(1:6))

    ## TODO: strand
}

test_pmapToGenome <- function() {
    GA_1P <- rep(GA_1[1], length(GR_1))
    ans <- pmapToGenome(GR_1, GA_1P)
    checkTrue(ncol(mcols(ans)) == 0L)
    checkTrue(length(ans) == length(GR_1))
    checkIdentical(width(ans), c(2L, 2L, 0L, 0L))
    checkIdentical(start(ans)[3:4], width(GA_1P)[3:4] + 1L)

    GR_2P <- rep(GR_2[1], length(GA_2))
    ans <- pmapToGenome(GR_2P, GA_2)
    checkIdentical(width(ans), c(6L, 8L, 8L, 4L, 6L, 6L))
    checkIdentical(end(ans), c(15L, 17L, 17L, 13L, 15L, 15L)) 

    ## strand
    strand(GR_2P) <- "-"
    ans <- pmapToGenome(GR_2P, GA_2, ignore.strand=FALSE)
    checkTrue(length(ans) == length(GR_2P))
    checkTrue(all(width(ans)) == 0L)
}

test_mapToTranscript <- function() {
    ans <- mapToTranscript(GR_1, GA_1)
    checkIdentical(start(ans), rep(1L, 3))
    checkIdentical(end(ans), rep(2L, 3))
    checkIdentical(mcols(ans)$from_hits, c(rep(2L, 3)))
    checkIdentical(mcols(ans)$to_hits, c(1L, 2L, 3L))

    ans <- mapToTranscript(GR_2, GA_3)
    checkIdentical(end(ans), c(7L, 4L, 4L, 8L, 6L, 6L)) 
    checkIdentical(mcols(ans)$to_hits, as.integer(1:6))

    ## TODO: strand
}

test_pmapToTranscript <- function() {
    GA_1P <- rep(GA_1[1], length(GR_1))
    ans <- pmapToTranscript(GR_1, GA_1P)
    checkIdentical(length(ans), length(GR_1))
    checkIdentical(width(ans), c(0L, 2L, 0L, 0L))
    checkIdentical(start(ans), c(1L, 1L, 12L, 12L))
    checkIdentical(end(ans), c(0L, 2L, 11L, 11L))

    GR_2P <- rep(GR_2[2], length(GA_2))
    ans <- pmapToTranscript(GR_2P, GA_3)
    checkIdentical(width(ans), c(6L, 4L, 4L, 8L, 6L, 6L))
    checkIdentical(start(ans), c(2L, rep(1L, 5))) 
    checkIdentical(end(ans), c(7L, 4L, 4L, 8L, 6L, 6L)) 

    ## strand
    strand(GR_2P) <- "-"
    ans <- pmapToTranscript(GR_2P, GA_2, ignore.strand=FALSE)
    checkTrue(length(ans) == length(GR_2P))
    checkTrue(all(width(ans)) == 0L)
    checkIdentical(end(ans), width(GA_2))
}

test_ref_locs_to_query_locs <- function() {
    cigar <- "66S42M2I20M8I18D15M43243N5M1D38M1D85M1D115M139S"
    pos <- 525842L
    ref <- 43425L + pos - 1L
    query <- 238L
    ans <- .Call("ref_locs_to_query_locs", ref, cigar, pos, 
                 FALSE, PACKAGE="GenomicAlignments")
    checkIdentical(ans, query)

    ## out of bounds
    ans_s <- .Call("ref_locs_to_query_locs", 
                   start(GR_1[1]), cigar(GA_1[1]), 
                   start(GA_1[1]), FALSE, 
                   PACKAGE="GenomicAlignments")
    ans_e <- .Call("ref_locs_to_query_locs", 
                   end(GR_1[1]), cigar(GA_1[1]), 
                   start(GA_1[1]), TRUE, 
                   PACKAGE="GenomicAlignments")

    checkIdentical(ans_s, NA_integer_)
    checkIdentical(ans_e, NA_integer_)
}

test_query_locs_to_ref_locs <- function() {
    ## out of bounds
    ans_s <- .Call("query_locs_to_ref_locs", 
                   start(GR_1[4]), cigar(GA_1[1]), 
                   start(GA_1[1]), FALSE, 
                   PACKAGE="GenomicAlignments")
    ans_e <- .Call("query_locs_to_ref_locs", 
                   end(GR_1[4]), cigar(GA_1[1]), 
                   start(GA_1[1]), TRUE, 
                   PACKAGE="GenomicAlignments")

    checkIdentical(ans_s, NA_integer_)
    checkIdentical(ans_e, NA_integer_)
}
