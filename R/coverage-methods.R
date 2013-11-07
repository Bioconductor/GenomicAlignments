### =========================================================================
### "coverage" methods
### -------------------------------------------------------------------------


setMethod("coverage", "BamFile",
    function(x, shift=0L, width=NULL, weight=1L, ..., param=ScanBamParam())
        coverage(readGAlignmentsFromBam(x, param=param),
                 shift=shift, width=width, weight=weight, ...)
)

