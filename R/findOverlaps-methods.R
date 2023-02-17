### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------


setMethod("findOverlaps", c("GAlignments", "Vector"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("Vector", "GAlignments"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             ignore.strand = FALSE)
    {
        findOverlaps(query, grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

### Not strictly needed! Defining the above 2 methods covers that case but
### with the following note:
###   > findOverlaps(al1, al0)
###   Note: Method with signature "GAlignments#ANY" chosen for
###    function "findOverlaps", target signature
###    "GAlignments#GAlignments".
###    "ANY#GAlignments" would also be valid
setMethod("findOverlaps", c("GAlignments", "GAlignments"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignmentPairs", "Vector"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("Vector", "GAlignmentPairs"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             ignore.strand = FALSE)
    {
        findOverlaps(query, grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignmentPairs", "GAlignmentPairs"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignmentsList", "Vector"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             ignore.strand = FALSE)
    {
        ## to take into account the right (real) strand, when ignore.strand=FALSE,
        ## we duplicate GAligbnmentsList object in memory and update the original
        ## strand with the real one, according to the strandMode parameter
        ## TODO: it may be worth investigating alternative ways to avoid duplicating
        ## a presumably large GAlignmentsList object
        query2 <- query
        if (!ignore.strand)
          strand(query2) <- strand(query2) ## overwrite original strand w/ real one
        hits <- findOverlaps(grglist(unlist(query2, use.names=FALSE)),
                             subject, maxgap = maxgap, minoverlap = minoverlap,
                             type = match.arg(type), select = match.arg(select),
                             ignore.strand = ignore.strand)
        ## TODO: Replace 'factor(togroup(PartitioningByWidth(...)))' with
        ## 'as.factor(PartitioningByWidth(...))' when "as.factor" method for
        ## ManyToOneGrouping objects becomes available.
        query_map <- factor(togroup(PartitioningByWidth(query)))
        remapHits(hits, Lnodes.remapping=query_map)
    }
)

setMethod("findOverlaps", c("Vector", "GAlignmentsList"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             ignore.strand = FALSE)
    {
        ## to take into account the right (real) strand, when ignore.strand=FALSE,
        ## we duplicate GAligbnmentsList object in memory and update the original
        ## strand with the real one, according to the strandMode parameter
        ## TODO: it may be worth investigating alternative ways to avoid duplicating
        ## a presumably large GAlignmentsList object
        subject2 <- subject
        if (!ignore.strand)
          strand(subject2) <- strand(subject2) ## overwrite original strand w/ real one
        hits <- findOverlaps(query, grglist(unlist(subject2, use.names = FALSE)),
                             maxgap = maxgap, minoverlap = minoverlap,
                             type = match.arg(type), select = match.arg(select),
                             ignore.strand = ignore.strand)
        ## TODO: Replace 'factor(togroup(PartitioningByWidth(...)))' with
        ## 'as.factor(PartitioningByWidth(...))' when "as.factor" method for
        ## ManyToOneGrouping objects becomes available.
        subject_map <- factor(togroup(PartitioningByWidth(subject)))
        remapHits(hits, Rnodes.remapping=subject_map)
    }
)

setMethod("findOverlaps", c("GAlignmentsList", "GAlignmentsList"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             ignore.strand = FALSE)
    {
        ## to take into account the right (real) strand, when ignore.strand=FALSE,
        ## we duplicate GAligbnmentsList object in memory and update the original
        ## strand with the real one, according to the strandMode parameter
        ## TODO: it may be worth investigating alternative ways to avoid duplicating
        ## a presumably large GAlignmentsList object
        query2 <- query
        subject2 <- subject
        if (!ignore.strand) {
          strand(query2) <- strand(query2) ## overwrite original strand w/ real one
          strand(subject2) <- strand(subject2) ## overwrite original strand w/ real one
        }
        hits <- findOverlaps(grglist(unlist(query2, use.names = FALSE)), 
                             grglist(unlist(subject2, use.names = FALSE)),
                             maxgap = maxgap, minoverlap = minoverlap,
                             type = match.arg(type), 
                             select = match.arg(select),
                             ignore.strand = ignore.strand)
        ## TODO: Replace 'factor(togroup(PartitioningByWidth(...)))' with
        ## 'as.factor(PartitioningByWidth(...))' when "as.factor" method for
        ## ManyToOneGrouping objects becomes available.
        query_map <- factor(togroup(PartitioningByWidth(query)))
        subject_map <- factor(togroup(PartitioningByWidth(subject)))
        remapHits(hits, Lnodes.remapping=query_map,
                        Rnodes.remapping=subject_map)
    }
)

