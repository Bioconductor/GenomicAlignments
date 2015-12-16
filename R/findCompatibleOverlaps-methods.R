### =========================================================================
### findCompatibleOverlaps
### -------------------------------------------------------------------------
###


setGeneric("findCompatibleOverlaps",
    function(query, subject) standardGeneric("findCompatibleOverlaps")
)

.GAlignmentsORGAlignmentPairs.findCompatibleOverlaps <-
    function(query, subject)
{
    ## Starting with BioC 3.2, the 'order.as.in.query' argument is not
    ## supported anymore for GAlignmentPairs objects.
    if (is(query, "GAlignmentPairs"))
        grl <- grglist(query)
    else
        grl <- grglist(query, order.as.in.query=TRUE)
    ## TODO: Use 'type="within"' when it's supported for circular
    ## sequences like the mitochondrial chromosome.
    ov <- findOverlaps(grl, subject, ignore.strand=TRUE)
    ovenc <- encodeOverlaps(grl, subject, hits=ov,
                            flip.query.if.wrong.strand=TRUE)
    ov_is_compat <- isCompatibleWithSplicing(ovenc)
    ov[ov_is_compat]
}

setMethod("findCompatibleOverlaps", c("GAlignments", "GRangesList"),
    .GAlignmentsORGAlignmentPairs.findCompatibleOverlaps
)

setMethod("findCompatibleOverlaps", c("GAlignmentPairs", "GRangesList"),
    .GAlignmentsORGAlignmentPairs.findCompatibleOverlaps
)

countCompatibleOverlaps <- function(query, subject)
{
    compatov <- findCompatibleOverlaps(query, subject)
    tabulate(queryHits(compatov), nbins=queryLength(compatov))
}

