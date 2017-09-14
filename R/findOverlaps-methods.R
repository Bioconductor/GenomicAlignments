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
        hits <- findOverlaps(grglist(unlist(query, use.names = FALSE)),
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
        hits <- findOverlaps(query, grglist(unlist(subject, use.names = FALSE)),
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
        hits <- findOverlaps(grglist(unlist(query, use.names = FALSE)), 
                             grglist(unlist(subject, use.names = FALSE)),
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

