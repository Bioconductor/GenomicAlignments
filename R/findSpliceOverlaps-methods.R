### =========================================================================
### "findSpliceOverlaps" methods
### -------------------------------------------------------------------------


setMethod("findSpliceOverlaps", c("GAlignments", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    callGeneric(grglist(query, order.as.in.query=TRUE), subject,
                ignore.strand, ..., cds=cds)
})

setMethod("findSpliceOverlaps", c("GAlignmentPairs", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
### FIXME: order.as.in.query = FALSE needed for insertGaps(). If we
### really want to use insertGaps(), we need to make it robust to
### different orderings.

### FIXME:
### instead of relying on query.break column, maybe we should add a
### 'splice = .gaps(query)' argument to .findSpliceOverlaps that we
### set to introns(query) here. The downside is that a GRangesList
### derived from GAlignmentPairs will no longer work.
    callGeneric(grglist(query, order.as.in.query=FALSE), subject,
                ignore.strand, ..., cds=cds)
})

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

