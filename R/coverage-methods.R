### =========================================================================
### "coverage" methods
### -------------------------------------------------------------------------


setMethod("coverage", "GAlignments",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash"), drop.D.ranges=FALSE)
        coverage(grglist(x, drop.D.ranges=drop.D.ranges),
                 shift=shift, width=width, weight=weight, method=method)
)

setMethod("coverage", "GAlignmentPairs",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash"), drop.D.ranges=FALSE)
        coverage(grglist(x, drop.D.ranges=drop.D.ranges),
                 shift=shift, width=width, weight=weight, method=method)
)

setMethod("coverage", "BamFile",
    function(x, shift=0L, width=NULL, weight=1L, ..., param=ScanBamParam())
{
    if (!isOpen(x)) {
        open(x)
        on.exit(close(x))
    }

    cvg <- NULL
    repeat {
        aln <- readGAlignments(x, param=param)
        if (length(aln) == 0L) {
            if (is.null(cvg))
                cvg <- coverage(aln, shift=shift, width=width, 
                                weight=weight, ...)
            break
        }
        cvg0 <- coverage(aln, shift=shift, width=width, weight=weight, ...)
        if (is.null(cvg))
            cvg <- cvg0
        else
            cvg <- cvg + cvg0
    }
    cvg
})

setMethod("coverage", "character",
    function(x, shift=0L, width=NULL, weight=1L, ..., yieldSize=2500000L)
{
    if (!isSingleString(x))
        stop("'x' must be character(1) for coverage,character-method")
    bf <- BamFile(x, yieldSize=yieldSize)
    coverage(bf, shift=shift, width=width, weight=weight, ...)
})
