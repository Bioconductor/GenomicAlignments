### =========================================================================
### Intra-range methods
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateCigarAndStart() -- NOT exported
###

setGeneric("updateCigarAndStart",
    function(x, cigar=NULL, start=NULL) standardGeneric("updateCigarAndStart")
)

setMethod("updateCigarAndStart", "GAlignments",
    function(x, cigar=NULL, start=NULL)
    {
        if (is.null(cigar)) {
            cigar <- cigar(x)
        } else {
            if (!is.character(cigar) || length(cigar) != length(x))
                stop("when not NULL, 'cigar' must be a character vector ",
                     "of the same length as 'x'")
            ## There might be an "rshift" attribute on 'cigar', typically.
            ## We want to get rid of it as well as any other potential
            ## attribute like names, dim, dimnames etc...
            attributes(cigar) <- NULL
        }
        if (is.null(start))
            start <- start(x)
        else if (!is.integer(start) || length(start) != length(x))
            stop("when not NULL, 'start' must be an integer vector ",
                 "of the same length as 'x'")
        x@cigar <- cigar
        x@start <- start
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### narrow()
###

setMethod("narrow", "GAlignments",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
        .narrowGAlignments(x, cigarNarrow, start, end, width)
)

setMethod("narrow", "GAlignmentsList",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        gal <- narrow(x@unlistData, start=start, end=end, width=width,
                      use.names=use.names)
        relist(gal, x@partitioning)
    }
)

setMethod("narrow", "GappedReads",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        stop("coming soon")
        ## ans_cigar <- cigarNarrow(cigar(x),
        ##                          start=start, end=end, width=width)
        ## ans_start <- start(x) + attr(ans_cigar, "rshift")
        ## updateCigarAndStart(x, cigar=ans_cigar, start=ans_start)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### qnarrow()
###

setGeneric("qnarrow", signature="x",
    function(x, start=NA, end=NA, width=NA) standardGeneric("qnarrow")
)

.narrowGAlignments <- function(x, CIGAR_CUTTER, start, end, width)
{
    ans_cigar <- CIGAR_CUTTER(cigar(x), start=start, end=end, width=width)
    ans_start <- start(x) + attr(ans_cigar, "rshift")
    updateCigarAndStart(x, cigar=ans_cigar, start=ans_start)
}

setMethod("qnarrow", "GAlignments",
    function(x, start=NA, end=NA, width=NA)
        .narrowGAlignments(x, cigarQNarrow, start, end, width)
)

setMethod("qnarrow", "GAlignmentsList",
    function(x, start=NA, end=NA, width=NA)
    {
        gal <- qnarrow(x@unlistData, start=start, end=end, width=width)
        relist(gal, x@partitioning)
    }
)

setMethod("qnarrow", "GappedReads",
    function(x, start=NA, end=NA, width=NA)
    {
        stop("coming soon")
        ## ans_cigar <- cigarQNarrow(cigar(x),
        ##                           start=start, end=end, width=width)
        ## ans_start <- start(x) + attr(ans_cigar, "rshift")
        ## updateCigarAndStart(x, cigar=ans_cigar, start=ans_start)
    }
)

