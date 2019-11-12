## CIGAR ops M, =, X
x1 <- GRanges("chr1", IRanges(c(5, 10, 20, 25), width=2, names=LETTERS[1:4]))
align1 <- GAlignments(rep("chr1", 3), rep(10L, 3), c("11M", "11=", "11X"), 
                      strand(rep("+", 3)), names=letters[1:3])

## CIGAR ops S, N, D, I, H, P
x2 <- GRanges("chr1", IRanges(c(1, 20), width=6, names=LETTERS[1:2]))
cigar <- c("1S6M1S", "3M2N3M", "3M2D3M", "3M2I3M", "1H6M1H", "1P6M1P")
align2 <- GAlignments(rep("chr1", 6), rep(10L, 6), cigar, strand(rep("+", 6)))
align3 <- GAlignments(rep("chr1", 6), rep(20L, 6), cigar, strand(rep("+", 6)))
names(align2) <- names(align3) <- letters[1:6]

test_mapToAlignments <- function() {
    ans <- mapToAlignments(x1, align1)
    checkIdentical(start(ans), rep(1L, 3))
    checkIdentical(end(ans), rep(2L, 3))
    checkIdentical(mcols(ans)$xHits, c(rep(2L, 3)))
    checkIdentical(mcols(ans)$alignmentsHits, c(1L, 2L, 3L))
    checkIdentical(names(ans), rep("B", 3))
    checkIdentical(seqlevels(ans), letters[1:3])

    ans <- mapToAlignments(x2, align3)
    checkIdentical(end(ans), c(7L, 4L, 4L, 8L, 6L, 6L)) 
    checkIdentical(mcols(ans)$alignmentsHits, as.integer(1:6))
}

test_mapFromAlignments <- function() {
    x <- x1
    names(x) <- rep("all", length(x)) 
    align <- align1
    names(align) <- rep("all", length(align)) 
    ans <- mapFromAlignments(x, align)
    checkIdentical(start(ans), c(14L, 14L, 14L, 19L, 19L, 19L))
    checkIdentical(mcols(ans)$xHits, c(rep(1L, 3), rep(2L, 3)))
    checkIdentical(mcols(ans)$alignmentsHits, rep(1:3, 2))
    checkIdentical(seqlevels(ans), "chr1")
    checkIdentical(names(ans), rep("all", 6))

    names(x) <- c("hit", "hit", "blank", "blank") 
    names(align) <- c("BLANK", "hit", "BLANK") 
    ans <- mapFromAlignments(x, align)
    checkIdentical(names(ans), c("hit", "hit"))
    checkIdentical(seqlevels(ans), "chr1")

    x <- x2
    names(x) <- rep("all", length(x)) 
    align <- align2
    names(align) <- rep("all", length(align)) 
    ans <- mapFromAlignments(x, align)
    checkIdentical(start(ans), rep(10L, 6))
    checkIdentical(end(ans), c(15L, 17L, 17L, 13L, 15L, 15L)) 
    checkIdentical(mcols(ans)$alignmentsHits, as.integer(1:6))
}

test_pmapToAlignments <- function() {
    x <- x1
    align <- rep(align1[1], length(x1))
    ans <- pmapToAlignments(x, align)
    checkIdentical(length(ans), length(x))
    checkIdentical(width(ans), c(0L, 2L, 0L, 0L))
    checkIdentical(start(ans), c(0L, 1L, 0L, 0L))
    checkIdentical(end(ans), c(-1L, 2L, -1L, -1L))
    checkIdentical(names(ans), names(x))
    checkTrue(all(seqlevels(ans) %in% c("a", "UNMAPPED"))) 

    x <- rep(x2[2], length(align3))
    align <- align3
    ans <- pmapToAlignments(x, align)
    checkIdentical(width(ans), c(6L, 4L, 4L, 8L, 6L, 6L))
    checkIdentical(start(ans), c(2L, rep(1L, 5))) 
    checkIdentical(end(ans), c(7L, 4L, 4L, 8L, 6L, 6L)) 
}

test_pmapFromAlignments <- function() {
    x <- x1
    names(x) <- rep("all", length(x))
    align <- rep(align1[1], length(x1))
    names(align) <- rep("all", length(align))
    ans <- pmapFromAlignments(x, align)
    checkTrue(ncol(mcols(ans)) == 0L)
    checkTrue(length(ans) == length(x1))
    checkIdentical(width(ans), c(2L, 2L, 0L, 0L))
    checkIdentical(start(ans), c(14L, 19L, 0L, 0L))
    checkIdentical(names(ans), names(x))
    checkTrue(all(seqlevels(ans) %in% c("chr1", "UNMAPPED"))) 

    x <- rep(x2[1], length(align2))
    names(x) <- LETTERS[seq_along(x)]
    align <- align2
    ans <- pmapFromAlignments(x, align)
    checkIdentical(width(ans), c(6L, 8L, 8L, 4L, 6L, 6L))
    checkIdentical(end(ans), c(15L, 17L, 17L, 13L, 15L, 15L)) 
    checkIdentical(names(ans), names(x))
    checkTrue(all(seqlevels(ans) %in% "chr1")) 
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
                   start(x1[1]), cigar(align1[1]), 
                   start(align1[1]), FALSE, 
                   PACKAGE="GenomicAlignments")
    ans_e <- .Call("ref_locs_to_query_locs", 
                   end(x1[1]), cigar(align1[1]), 
                   start(align1[1]), TRUE, 
                   PACKAGE="GenomicAlignments")

    checkIdentical(ans_s, NA_integer_)
    checkIdentical(ans_e, NA_integer_)
}

test_query_locs_to_ref_locs <- function() {
    ## out of bounds
    ans_s <- .Call("query_locs_to_ref_locs", 
                   start(x1[4]), cigar(align1[1]), 
                   start(align1[1]), FALSE, 
                   PACKAGE="GenomicAlignments")
    ans_e <- .Call("query_locs_to_ref_locs", 
                   end(x1[4]), cigar(align1[1]), 
                   start(align1[1]), TRUE, 
                   PACKAGE="GenomicAlignments")

    checkIdentical(ans_s, NA_integer_)
    checkIdentical(ans_e, NA_integer_)
}

test_map_ref_locs_to_query_locs <- function() {
    ## hit
    map <- .Call("map_ref_locs_to_query_locs", 
                 12L, 16L, "11M", 10L)
    checkIdentical(unlist(map), c(3L, 7L, 1L, 1L))

    ## first record out of bounds
    map <- .Call("map_ref_locs_to_query_locs", 
                 c(5L, 12L), c(16L, 16L), "11M", 10L)
    checkIdentical(unlist(map), c(3L, 7L, 2L, 1L))

    ## second record out of bounds
    map <- .Call("map_ref_locs_to_query_locs", 
                 c(12L, 5L), c(16L, 16L), "11M", 10L)
    checkIdentical(unlist(map), c(3L, 7L, 1L, 1L))

    ## first alignment out of bounds
    map <- .Call("map_ref_locs_to_query_locs", 
                 12L, 16L, c("11M", "11M"), c(20L, 10L))
    checkIdentical(unlist(map), c(3L, 7L, 1L, 2L))

    ## second alignment out of bounds
    map <- .Call("map_ref_locs_to_query_locs", 
                 12L, 16L, c("11M", "11M"), c(10L, 20L))
    checkIdentical(unlist(map), c(3L, 7L, 1L, 1L))
}
