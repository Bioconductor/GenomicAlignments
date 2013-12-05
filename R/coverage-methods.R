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
        coverage(readGAlignmentsFromBam(x, param=param),
                 shift=shift, width=width, weight=weight, ...)
)

