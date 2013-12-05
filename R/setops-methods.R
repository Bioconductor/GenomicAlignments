### =========================================================================
### Set operations
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pintersect()
###

### TODO: Revisit this method (seems to do strange things).
setMethod("pintersect", c("GAlignments", "GRanges"),
    function(x, y, ...)
    {
        bounds <- try(callGeneric(granges(x), y), silent=TRUE)
        if (inherits(bounds, "try-error"))
            stop("CIGAR is empty after intersection")
        start <- start(bounds) - start(x) + 1L
        start[which(start < 1L)] <- 1L
        end <- end(bounds) - end(x) - 1L
        end[which(end > -1L)] <- -1L
        narrow(x, start=start, end=end)
    }
)

setMethod("pintersect", c("GRanges", "GAlignments"),
    function(x, y, ...)
    {
        callGeneric(y, x)
    }
)

