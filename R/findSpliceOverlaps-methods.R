### =========================================================================
### "findSpliceOverlaps" methods
### -------------------------------------------------------------------------


setMethod("findSpliceOverlaps", c("character", "ANY"),
          function(query, subject, ignore.strand=FALSE, ...,
                   param=ScanBamParam(), singleEnd=TRUE)
{
    findSpliceOverlaps(BamFile(query), subject, ignore.strand, ...,
                       param=param, singleEnd=singleEnd)
})

setMethod("findSpliceOverlaps", c("BamFile", "ANY"),
    function(query, subject, ignore.strand=FALSE, ...,
             param=ScanBamParam(), singleEnd=TRUE)
{
    findSpliceOverlaps(.readRanges(query, param, singleEnd), subject,
                       ignore.strand, ...)
})

.readRanges <- function(bam, param, singleEnd)
{
    if (!"XS" %in% bamTag(param))
        bamTag(param) <- c(bamTag(param), "XS")
    if (singleEnd)
        reads <- readGAlignmentsFromBam(bam, param=param)
    else {
        reads <- readGAlignmentPairsFromBam(path(bam), param=param)
        first_xs <- mcols(first(reads))$XS
        last_xs <- mcols(last(reads))$XS
        if (!is.null(first_xs) && !is.null(last_xs)) {
            xs <- first_xs
            xs[is.na(xs)] <- last_xs[is.na(xs)]
            mcols(reads)$XS <- xs
        }
    }

    metadata(reads)$bamfile <- bam

    reads
}

