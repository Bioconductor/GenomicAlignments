### =========================================================================
### Low-level CIGAR utilities
### -------------------------------------------------------------------------


### See p. 4 of the SAM Spec v1.4 at http://samtools.sourceforge.net/ for the
### list of CIGAR operations and their meanings.
CIGAR_OPS <- c("M", "I", "D", "N", "S", "H", "P", "=", "X")

.normarg_cigar <- function(cigar)
{
    if (is.factor(cigar))
        cigar <- as.character(cigar)
    if (!is.character(cigar))
        stop("'cigar' must be a character vector or factor")
    cigar
}

.normarg_flag <- function(flag, cigar)
{
    if (!is.null(flag)) {
        if (!is.numeric(flag))
            stop("'flag' must be NULL or a vector of integers")
        if (!is.integer(flag))
            flag <- as.integer(flag)
        if (length(cigar) != length(flag))
            stop("'cigar' and 'flag' must have the same length")
    }
    flag
}

.normarg_pos <- function(pos, cigar)
{
    if (!is.numeric(pos))
        stop("'pos' must be a vector of integers")
    if (!is.integer(pos))
        pos <- as.integer(pos)
    if (length(pos) != 1L && length(pos) != length(cigar))
        stop("'pos' must have length 1 or the same length as 'cigar'")
    pos
}

.select_reference_space <- function(N.regions.removed)
{
    if (!isTRUEorFALSE(N.regions.removed))
        stop("'N.regions.removed' must be TRUE or FALSE")
    if (N.regions.removed) {
        space <- 2L  # REFERENCE_N_REGIONS_REMOVED
    } else {
        space <- 1L  # REFERENCE
    }
    space
}

.select_query_space <- function(before.hard.clipping, after.soft.clipping)
{
    if (!isTRUEorFALSE(before.hard.clipping))
        stop("'before.hard.clipping' must be TRUE or FALSE")
    if (!isTRUEorFALSE(after.soft.clipping))
        stop("'after.soft.clipping' must be TRUE or FALSE")
    if (before.hard.clipping) {
        if (after.soft.clipping)
            stop("'before.hard.clipping' and 'after.soft.clipping' ",
                 "cannot both be TRUE")
        space <- 4L  # QUERY_BEFORE_HARD_CLIPPING
    } else if (after.soft.clipping) {
        space <- 5L  # QUERY_AFTER_SOFT_CLIPPING
    } else {
        space <- 3L  # QUERY
    }
    space
}

.select_pairwise_space <- function(N.regions.removed, dense)
{
    if (!isTRUEorFALSE(N.regions.removed))
        stop("'N.regions.removed' must be TRUE or FALSE")
    if (!isTRUEorFALSE(dense))
        stop("'dense' must be TRUE or FALSE")
    if (N.regions.removed) {
        if (dense)
           stop("'N.regions.removed' and 'dense' ",
                "cannot both be TRUE")
        space <- 7L  # PAIRWISE_N_REGIONS_REMOVED
    } else if (dense) {
        space <- 8L  # PAIRWISE_DENSE
    } else {
        space <- 6L  # PAIRWISE
    }
    space
}

.normarg_ops <- function(ops)
{
    if (is.null(ops))
        return(ops)
    if (!is.character(ops))
        stop("'ops' must be a character vector")
    if (any(is.na(ops)))
        stop("'ops' cannot contain NAs")
    if (length(ops) == 1L) {
        ops <- strsplit(ops, NULL, fixed=TRUE)[[1L]]
    } else if (any(nchar(ops) != 1L)) {
        stop("when 'length(ops) != 1', all its elements ",
             "must be single letters")
    }
    if (anyDuplicated(ops))
        stop("'ops' cannot contain duplicated letters")
    if (!all(ops %in% CIGAR_OPS))
        stop("'ops' contains invalid CIGAR operations")
    ops
}

validCigar <- function(cigar)
{
    cigar <- .normarg_cigar(cigar)
    .Call2("valid_cigar", cigar, 0L, PACKAGE="GenomicAlignments")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transform CIGARs into other useful representations
###

explodeCigarOps <- function(cigar, ops=CIGAR_OPS)
{
    cigar <- .normarg_cigar(cigar)
    ops <- .normarg_ops(ops)
    .Call2("explode_cigar_ops", cigar, ops, PACKAGE="GenomicAlignments")
}

explodeCigarOpLengths <- function(cigar, ops=CIGAR_OPS)
{
    cigar <- .normarg_cigar(cigar)
    ops <- .normarg_ops(ops)
    .Call2("explode_cigar_op_lengths", cigar, ops, PACKAGE="GenomicAlignments")
}

cigarToRleList <- function(cigar)
{
    cigar_ops <- explodeCigarOps(cigar)
    cigar_op_lengths <- explodeCigarOpLengths(cigar)
    if (length(cigar) == 0L) {
        unlisted_cigar_ops <- character(0)
        unlisted_cigar_op_lengths <- integer(0)
    } else {
        unlisted_cigar_ops <- unlist(cigar_ops, use.names=FALSE)
        unlisted_cigar_op_lengths <- unlist(cigar_op_lengths, use.names=FALSE)
    }

    ## Prepare 'ans_flesh'.
    ans_flesh <- Rle(unlisted_cigar_ops, unlisted_cigar_op_lengths)

    ## Prepare 'ans_skeleton'.
    nops_per_cigar <- elementNROWS(cigar_op_lengths)
    ans_breakpoints <- cumsum(unlisted_cigar_op_lengths)[cumsum(nops_per_cigar)]
    ans_skeleton <- PartitioningByEnd(ans_breakpoints)

    ## Relist.
    relist(ans_flesh, ans_skeleton)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Summarize CIGARs
###

cigarOpTable <- function(cigar)
{
    cigar <- .normarg_cigar(cigar)
    ans <- .Call2("cigar_op_table", cigar, PACKAGE="GenomicAlignments")
    stopifnot(identical(CIGAR_OPS, colnames(ans)))  # sanity check
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### From CIGARs to ranges
###

.cigar_ranges <- function(cigar, flag, space, pos, f,
                          ops, drop.empty.ranges, reduce.ranges, with.ops)
{
    cigar <- .normarg_cigar(cigar)
    flag <- .normarg_flag(flag, cigar)
    if (!isSingleNumber(space))
        stop("'space' must be a single integer")
    if (!is.integer(space))
        space <- as.integer(space)
    pos <- .normarg_pos(pos, cigar)
    if (!is.null(f)) {
        if (!is.factor(f))
            stop("'f' must be NULL or a factor")
        if (length(f) != length(cigar))
            stop("'f' must have the same length as 'cigar'")
    }
    ops <- .normarg_ops(ops)
    if (!isTRUEorFALSE(drop.empty.ranges))
        stop("'drop.empty.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(reduce.ranges))
        stop("'reduce.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.ops))
        stop("'with.ops' must be TRUE or FALSE")
    .Call2("cigar_ranges",
           cigar, flag, space, pos, f,
           ops, drop.empty.ranges, reduce.ranges, with.ops,
           PACKAGE="GenomicAlignments")
}

cigarRangesAlongReferenceSpace <- function(cigar, flag=NULL,
                                    N.regions.removed=FALSE, pos=1L, f=NULL,
                                    ops=CIGAR_OPS,
                                    drop.empty.ranges=FALSE,
                                    reduce.ranges=FALSE,
                                    with.ops=FALSE)
{
    space <- .select_reference_space(N.regions.removed)
    C_ans <- .cigar_ranges(cigar, flag, space, pos, f,
                           ops, drop.empty.ranges, reduce.ranges, with.ops)
    if (is.null(f))
        return(C_ans)
    compress <- length(C_ans) >= 200L
    IRangesList(C_ans, compress=compress)
}

cigarRangesAlongQuerySpace <- function(cigar, flag=NULL,
                                    before.hard.clipping=FALSE,
                                    after.soft.clipping=FALSE,
                                    ops=CIGAR_OPS,
                                    drop.empty.ranges=FALSE,
                                    reduce.ranges=FALSE,
                                    with.ops=FALSE)
{
    space <- .select_query_space(before.hard.clipping, after.soft.clipping)
    .cigar_ranges(cigar, flag, space, 1L, NULL,
                  ops, drop.empty.ranges, reduce.ranges, with.ops)
}

cigarRangesAlongPairwiseSpace <- function(cigar, flag=NULL,
                                    N.regions.removed=FALSE, dense=FALSE,
                                    ops=CIGAR_OPS,
                                    drop.empty.ranges=FALSE,
                                    reduce.ranges=FALSE,
                                    with.ops=FALSE)
{
    space <- .select_pairwise_space(N.regions.removed, dense)
    .cigar_ranges(cigar, flag, space, 1L, NULL,
                  ops, drop.empty.ranges, reduce.ranges, with.ops)
}

### A convenience wrapper to cigarRangesAlongReferenceSpace().
extractAlignmentRangesOnReference <- function(cigar, pos=1L,
                                              drop.D.ranges=FALSE, f=NULL)
{
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    ## Not sure why we include "I" operations here since they don't generate
    ## coverage on the reference (they always produce zero-width ranges on the
    ## reference).
    if (drop.D.ranges) {
        ops <- c("M", "=", "X", "I")
    } else {
        ops <- c("M", "=", "X", "I", "D")
    }
    cigarRangesAlongReferenceSpace(cigar, flag=NULL, pos=pos, f=f,
                                   ops=ops,
                                   drop.empty.ranges=FALSE, reduce.ranges=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### From CIGARs to sequence lengths
###

.cigar_width <- function(cigar, flag, space)
{
    cigar <- .normarg_cigar(cigar)
    flag <- .normarg_flag(flag, cigar)
    if (!isSingleNumber(space))
        stop("'space' must be a single integer")
    if (!is.integer(space))
        space <- as.integer(space)
    .Call2("cigar_width", cigar, flag, space, PACKAGE="GenomicAlignments")
}

cigarWidthAlongReferenceSpace <- function(cigar, flag=NULL,
                                          N.regions.removed=FALSE)
{
    space <- .select_reference_space(N.regions.removed)
    .cigar_width(cigar, flag, space)
}

cigarWidthAlongQuerySpace <- function(cigar, flag=NULL,
                                      before.hard.clipping=FALSE,
                                      after.soft.clipping=FALSE)
{
    space <- .select_query_space(before.hard.clipping, after.soft.clipping)
    .cigar_width(cigar, flag, space)
}

cigarWidthAlongPairwiseSpace <- function(cigar, flag=NULL,
                                         N.regions.removed=FALSE, dense=FALSE)
{
    space <- .select_pairwise_space(N.regions.removed, dense)
    .cigar_width(cigar, flag, space)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Narrow CIGARs
###

cigarNarrow <- function(cigar, start=NA, end=NA, width=NA)
{
    cigar_width <- cigarWidthAlongReferenceSpace(cigar)
    cigar_ranges <- IRanges(start=rep.int(1L, length(cigar_width)),
                            width=cigar_width)
    threeranges <- threebands(cigar_ranges, start=start, end=end, width=width)
    C_ans <- .Call2("cigar_narrow",
                   cigar, width(threeranges$left), width(threeranges$right),
                   PACKAGE="GenomicAlignments")
    ans <- C_ans[[1L]]
    attr(ans, "rshift") <- C_ans[[2L]]
    ans
}

cigarQNarrow <- function(cigar, start=NA, end=NA, width=NA)
{
    cigar_qwidth <- cigarWidthAlongQuerySpace(cigar)
    cigar_qranges <- IRanges(start=rep.int(1L, length(cigar_qwidth)),
                             width=cigar_qwidth)
    threeranges <- threebands(cigar_qranges, start=start, end=end, width=width)
    C_ans <- .Call2("cigar_qnarrow",
                   cigar, width(threeranges$left), width(threeranges$right),
                   PACKAGE="GenomicAlignments")
    ans <- C_ans[[1L]]
    attr(ans, "rshift") <- C_ans[[2L]]
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Translate coordinates between query-based and reference-based
###

queryLoc2refLoc <- function(qloc, cigar, pos=1L)
{
    stop(wmsg("Not implemented yet, sorry! Please use mapToAlignments() ",
              "and family if you need to map coordinates between query reads ",
              "and reference genome on the basis of the reported CIGAR. ",
              "See '?mapToAlignments'"))
}

queryLocs2refLocs <- function(qlocs, cigar, pos=1L, flag=NULL)
{
    stop(wmsg("Not implemented yet, sorry! Please use mapToAlignments() ",
              "and family if you need to map coordinates between query reads ",
              "and reference genome on the basis of the reported CIGAR. ",
              "See '?mapToAlignments'"))
}

