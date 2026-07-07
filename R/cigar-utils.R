### =========================================================================
### Low-level CIGAR utilities
### -------------------------------------------------------------------------


### A convenience wrapper to cigars_as_ranges_along_ref().
### This is the only thing in this file that is not deprecated (we keep
### it in GenomicAlignments for now).
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
    cigars_as_ranges_along_ref(cigar, flags=NULL, lmmpos=pos, f=f, ops=ops,
                               drop.empty.ranges=FALSE, reduce.ranges=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Everything below is deprecated!
###

call_new_fun_in_cigarillo <- function(old_fun, new_fun, ...)
{
    msg <- c(old_fun, "() is defunct in GenomicAlignments >= 1.49.1 ",
             "and replaced with the ", new_fun, "() function from the ",
             "new cigarillo package")
    .Defunct(msg=wmsg(msg))
}

validCigar <- function(cigar)
{
    call_new_fun_in_cigarillo("validCigar", "validate_cigars", cigars=cigar)
}

explodeCigarOps <- function(cigar, ops=CIGAR_OPS)
{
    call_new_fun_in_cigarillo("explodeCigarOps", "explode_cigar_ops",
                              cigars=cigar, ops=ops)
}

explodeCigarOpLengths <- function(cigar, ops=CIGAR_OPS)
{
    call_new_fun_in_cigarillo("explodeCigarOpLengths", "explode_cigar_oplens",
                              cigars=cigar, ops=ops)
}

cigarToRleList <- function(cigar)
{
    call_new_fun_in_cigarillo("cigarToRleList", "cigars_as_RleList",
                              cigars=cigar)
}

cigarOpTable <- function(cigar)
{
    call_new_fun_in_cigarillo("cigarOpTable", "tabulate_cigar_ops",
                              cigars=cigar)
}

cigarRangesAlongReferenceSpace <- function(cigar, flag=NULL,
                                    N.regions.removed=FALSE, pos=1L, f=NULL,
                                    ops=CIGAR_OPS,
                                    drop.empty.ranges=FALSE,
                                    reduce.ranges=FALSE,
                                    with.ops=FALSE)
{
    call_new_fun_in_cigarillo("cigarRangesAlongReferenceSpace",
                              "cigars_as_ranges_along_ref",
                              cigars=cigar,
                              N.regions.removed=N.regions.removed,
                              flags=flag,
                              lmmpos=pos,
                              f=f,
                              ops=ops,
                              drop.empty.ranges=drop.empty.ranges,
                              reduce.ranges=reduce.ranges,
                              with.ops=with.ops)
}

cigarRangesAlongQuerySpace <- function(cigar, flag=NULL,
                                    before.hard.clipping=FALSE,
                                    after.soft.clipping=FALSE,
                                    ops=CIGAR_OPS,
                                    drop.empty.ranges=FALSE,
                                    reduce.ranges=FALSE,
                                    with.ops=FALSE)
{
    call_new_fun_in_cigarillo("cigarRangesAlongQuerySpace",
                              "cigars_as_ranges_along_query",
                              cigars=cigar,
                              before.hard.clipping=before.hard.clipping,
                              after.soft.clipping=after.soft.clipping,
                              flags=flag,
                              ops=ops,
                              drop.empty.ranges=drop.empty.ranges,
                              reduce.ranges=reduce.ranges,
                              with.ops=with.ops)
}

cigarRangesAlongPairwiseSpace <- function(cigar, flag=NULL,
                                    N.regions.removed=FALSE, dense=FALSE,
                                    ops=CIGAR_OPS,
                                    drop.empty.ranges=FALSE,
                                    reduce.ranges=FALSE,
                                    with.ops=FALSE)
{
    call_new_fun_in_cigarillo("cigarRangesAlongPairwiseSpace",
                              "cigars_as_ranges_along_pwa",
                              cigars=cigar,
                              N.regions.removed=N.regions.removed,
                              dense=dense,
                              flags=flag,
                              ops=ops,
                              drop.empty.ranges=drop.empty.ranges,
                              reduce.ranges=reduce.ranges,
                              with.ops=with.ops)
}

cigarWidthAlongReferenceSpace <- function(cigar, flag=NULL,
                                          N.regions.removed=FALSE)
{
    call_new_fun_in_cigarillo("cigarWidthAlongReferenceSpace",
                              "cigar_extent_along_ref",
                              cigars=cigar,
                              N.regions.removed=N.regions.removed,
                              flags=flag)
}

cigarWidthAlongQuerySpace <- function(cigar, flag=NULL,
                                      before.hard.clipping=FALSE,
                                      after.soft.clipping=FALSE)
{
    call_new_fun_in_cigarillo("cigarWidthAlongQuerySpace",
                              "cigar_extent_along_query",
                              cigars=cigar,
                              before.hard.clipping=before.hard.clipping,
                              after.soft.clipping=after.soft.clipping,
                              flags=flag)
}

cigarWidthAlongPairwiseSpace <- function(cigar, flag=NULL,
                                         N.regions.removed=FALSE, dense=FALSE)
{
    call_new_fun_in_cigarillo("cigarWidthAlongPairwiseSpace",
                              "cigar_extent_along_pwa",
                              cigars=cigar,
                              N.regions.removed=N.regions.removed,
                              dense=dense,
                              flags=flag)
}

cigarNarrow <- function(cigar, start=NA, end=NA, width=NA)
{
    call_new_fun_in_cigarillo("cigarNarrow",
                              "narrow_cigars_along_ref",
                              cigars=cigar, start=start, end=end, width=width)
}

cigarQNarrow <- function(cigar, start=NA, end=NA, width=NA)
{
    call_new_fun_in_cigarillo("cigarQNarrow",
                              "narrow_cigars_along_query",
                              cigars=cigar, start=start, end=end, width=width)
}

