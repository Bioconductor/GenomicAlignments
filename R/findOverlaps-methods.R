### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------


setMethod("findOverlaps", c("GAlignments", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     algorithm = match.arg(algorithm),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("Vector", "GAlignments"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand = FALSE)
    {
        findOverlaps(query, grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     algorithm = match.arg(algorithm),
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
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     algorithm = match.arg(algorithm),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignmentPairs", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     algorithm = match.arg(algorithm),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("Vector", "GAlignmentPairs"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand = FALSE)
    {
        findOverlaps(query, grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     algorithm = match.arg(algorithm),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignmentPairs", "GAlignmentPairs"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     algorithm = match.arg(algorithm),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignmentsList", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand = FALSE)
    {
        hits <- findOverlaps(grglist(unlist(query, use.names = FALSE)),
                             subject, maxgap = maxgap, minoverlap = minoverlap,
                             type = match.arg(type), select = match.arg(select),
                             algorithm = match.arg(algorithm),
                             ignore.strand = ignore.strand)
        remapHits(hits, query.map=factor(togroup(query)))
    }
)

setMethod("findOverlaps", c("Vector", "GAlignmentsList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand = FALSE)
    {
        hits <- findOverlaps(query, grglist(unlist(subject, use.names = FALSE)),
                             maxgap = maxgap, minoverlap = minoverlap,
                             type = match.arg(type), select = match.arg(select),
                             algorithm = match.arg(algorithm),
                             ignore.strand = ignore.strand)
        remapHits(hits, subject.map=factor(togroup(subject)))
    }
)

setMethod("findOverlaps", c("GAlignmentsList", "GAlignmentsList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand = FALSE)
    {
        hits <- findOverlaps(grglist(unlist(query, use.names = FALSE)), 
                             grglist(unlist(subject, use.names = FALSE)),
                             maxgap = maxgap, minoverlap = minoverlap,
                             type = match.arg(type), 
                             select = match.arg(select),
                             algorithm = match.arg(algorithm),
                             ignore.strand = ignore.strand)
        remapHits(hits, subject.map=factor(togroup(subject)),
                  query.map=factor(togroup(query)))
    }
)


### =========================================================================
### findOverlaps-based methods
### -------------------------------------------------------------------------

.signatures1 <- list(
    c("GAlignments", "Vector"),
    c("Vector", "GAlignments"),
    c("GAlignments", "GAlignments"),

    c("GAlignmentPairs", "Vector"),
    c("Vector", "GAlignmentPairs"),
    c("GAlignmentPairs", "GAlignmentPairs"),

    c("GAlignmentsList", "Vector"),
    c("Vector", "GAlignmentsList"),
    c("GAlignmentsList", "GAlignmentsList")
)

.signatures2 <- list(
    c("GAlignments", "GenomicRanges"),
    c("GenomicRanges", "GAlignments"),
    c("GAlignments", "GRangesList"),
    c("GRangesList", "GAlignments")
)

setMethods("countOverlaps", c(.signatures1, .signatures2),
    GenomicRanges:::countOverlaps.definition
)

setMethods("overlapsAny", .signatures1,
    GenomicRanges:::overlapsAny.definition
)

setMethods("subsetByOverlaps", .signatures1,
    GenomicRanges:::subsetByOverlaps.definition1
)

