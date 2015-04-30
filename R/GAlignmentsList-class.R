### =========================================================================
### GAlignmentsList objects
### -------------------------------------------------------------------------
###

setClass("GAlignmentsList",
    contains="CompressedList",
    representation(
        unlistData="GAlignments",
        elementMetadata="DataFrame"
    ),
    prototype(
        elementType="GAlignments"
    )
)

### Formal API:
###   names(x)    - NULL or character vector.
###   length(x)   - single integer. Nb of alignments in 'x'.
###   seqnames(x) - 'factor' Rle of the same length as 'x'.
###   rname(x)    - same as 'seqnames(x)'.
###   seqnames(x) <- value - replacement form of 'seqnames(x)'.
###   rname(x) <- value - same as 'seqnames(x) <- value'.
###   cigar(x)    - character vector of the same length as 'x'.
###   strand(x)   - 'factor' Rle of the same length as 'x' (levels: +, -, *).
###   qwidth(x)   - integer vector of the same length as 'x'.
###   start(x), end(x), width(x) - integer vectors of the same length as 'x'.
###   njunc(x)    - integer vector of the same length as 'x'.

###   grglist(x)  - GRangesList object of the same length as 'x'.
###   granges(x)  - GRanges object of the same length as 'x'.
###   rglist(x)   - CompressedIRangesList object of the same length as 'x'.
###   ranges(x)   - IRanges object of the same length as 'x'.

###   show(x)     - compact display in a data.frame-like fashion.
###   GAlignmentsList(x, ...) - constructor.
###   x[i]        - GAlignmentsList object of the same class as 'x'
###                 (endomorphism).
###
###   findOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GAlignments objects. Just a convenient wrapper for
###                 'findOverlaps(grglist(query), subject, ...)', etc...
###
###   countOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GAlignments objects. Just a convenient wrapper for
###                 'countOverlaps(grglist(query), subject, ...)', etc...
###
###   subsetByOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GAlignments objects.
###

###   qnarrow(x, start=NA, end=NA, width=NA) - GAlignmentsList object of the
###                 same length and class as 'x' (endomorphism).
###
###   narrow(x, start=NA, end=NA, width=NA) - GAlignmentsList object of the
###                 same length and class as 'x' (endomorphism).
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("seqnames", "GAlignmentsList", 
    function(x) 
        new2("CompressedRleList",
             unlistData=x@unlistData@seqnames, 
             partitioning=x@partitioning, check=FALSE)
)

setMethod("rname", "GAlignmentsList", 
    function(x) 
        new2("CompressedRleList",
             unlistData=x@unlistData@seqnames, 
             partitioning=x@partitioning, check=FALSE)
)

setMethod("cigar", "GAlignmentsList", 
    function(x) 
        new2("CompressedCharacterList",
             unlistData=x@unlistData@cigar, 
             partitioning=x@partitioning, check=FALSE)
)

setMethod("strand", "GAlignmentsList",
    function(x)
        new2("CompressedRleList",
             unlistData=x@unlistData@strand, 
             partitioning=x@partitioning, check=FALSE)
)

setMethod("qwidth", "GAlignmentsList",
    function(x)
        new2("CompressedIntegerList",
             unlistData=cigarWidthAlongQuerySpace(x@unlistData@cigar),
             partitioning=x@partitioning, check=FALSE)
)

setMethod("njunc", "GAlignmentsList",
    function(x)
        new2("CompressedIntegerList",
             unlistData=unname(elementLengths(rglist(x@unlistData))) - 1L,
             partitioning=x@partitioning, check=FALSE)
)

setMethod("start", "GAlignmentsList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData=x@unlistData@start,
             partitioning=x@partitioning, check=FALSE)
)

setMethod("end", "GAlignmentsList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData=end(x@unlistData),
             partitioning=x@partitioning, check=FALSE)
)

setMethod("width", "GAlignmentsList",
    function(x)
        new2("CompressedIntegerList",
             unlistData=cigarWidthAlongReferenceSpace(x@unlistData@cigar),
             partitioning=x@partitioning, check=FALSE)
)

setMethod("seqinfo", "GAlignmentsList", function(x) seqinfo(x@unlistData))

setMethod("elementMetadata", "GAlignmentsList",
    GenomicRanges:::getElementMetadataList
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

setReplaceMethod("rname", "GAlignmentsList",
    function(x, value) `seqnames<-`(x, value)
)

setReplaceMethod("elementMetadata", "GAlignmentsList", 
    GenomicRanges:::replaceElementMetadataList
)

setReplaceMethod("strand", "GAlignmentsList",
    GenomicRanges:::replaceStrandList
)

setReplaceMethod("seqinfo", "GAlignmentsList",
    GenomicRanges:::replaceSeqinfoList
)

setReplaceMethod("seqnames", "GAlignmentsList",
    GenomicRanges:::replaceSeqnamesList
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GAlignmentsList <- function(x)
{
   ## TDB: Currently known pitfalls are caught by
   ## GAlignments validity. 
}

setValidity2("GAlignmentsList", .valid.GAlignmentsList)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

GAlignmentsList <- function(...)
{
    listData <- list(...)
    if (length(listData) == 0L) {
        unlistData <- GAlignments()
    } else {
        if (length(listData) == 1L && is.list(listData[[1L]]))
            listData <- listData[[1L]]
        if (!all(sapply(listData, is, "GAlignments")))
            stop("all elements in '...' must be GAlignments objects")
        unlistData <- suppressWarnings(do.call("c", unname(listData)))
    }
    relist(unlistData, PartitioningByEnd(listData))
}

setMethod("updateObject", "GAlignmentsList",
    function(object, ..., verbose=FALSE)
    {
        if (verbose)
            message("updateObject(object = 'GAlignmentsList')")
        if (is(try(validObject(object@unlistData, complete=TRUE), silent=TRUE),
               "try-error")) {
            object@unlistData <- updateObject(object@unlistData)
            return(object)
        }
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setMethod("grglist", "GAlignmentsList",
    function(x, use.mcols=FALSE, order.as.in.query=FALSE, 
             drop.D.ranges=FALSE, ignore.strand=FALSE) 
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        if (!identical(order.as.in.query, FALSE)) {
            msg <- c("Starting with BioC 3.2, the \"grglist\" method for ",
                     "GAlignmentsList objects *always* returns the ranges ",
                     "\"ordered as in query\". Therefore the ",
                     "'order.as.in.query' argument is now ignored (and ",
                     "deprecated).")
            .Deprecated(msg=wmsg(msg))
        }
        if (!identical(drop.D.ranges, FALSE)) {
            msg <- c("Starting with BioC 3.2, the 'drop.D.ranges' ",
                     "argument is ignored and deprecated in the ",
                     "\"grglist\" method for GAlignmentsList objects.")
            .Deprecated(msg=wmsg(msg))
        }
        if (ignore.strand)
            strand(x@unlistData) <- "*"
        gr <- granges(x@unlistData, use.mcols=use.mcols)
        ans <- relist(gr, x@partitioning)
        names(ans) <- names(x)
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)
 
setMethod("granges", "GAlignmentsList",
    function(x, use.mcols=FALSE, ignore.strand=FALSE) 
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        if (ignore.strand)
            strand(x@unlistData) <- "*"
        msg <- paste0("For some list elements in 'x', the ranges are ",
                      "not aligned to the same chromosome and strand. ",
                      "Cannot extract a single range for them. ",
                      "As a consequence, the returned GRanges object ",
                      "is not parallel to 'x'.")
        rg <- range(grglist(x, ignore.strand=ignore.strand))
        is_one_to_one <- all(elementLengths(rg) == 1L)
        if (!is_one_to_one && all(width(x@partitioning) > 0)) {
            if (ignore.strand)
                warning(msg)
            else
                warning(paste0(msg, " Consider using 'ignore.strand=TRUE'."))
        }
        ans <- unlist(rg)
        if (is_one_to_one && use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)

setMethod("rglist", "GAlignmentsList",
    function(x, use.mcols=FALSE, order.as.in.query=FALSE, drop.D.ranges=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- relist(ranges(x@unlistData), x@partitioning)
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)

setMethod("ranges", "GAlignmentsList",
    function(x) 
        unlist(range(rglist(x)), use.names=FALSE)
)

setAs("GAlignmentsList", "GRangesList", 
    function(from) grglist(from, use.mcols=TRUE)
)
setAs("GAlignmentsList", "GRanges", 
    function(from) granges(from, use.mcols=TRUE)
)
setAs("GAlignmentsList", "RangesList", 
    function(from) rglist(from, use.mcols=TRUE)
)
setAs("GAlignmentsList", "Ranges", 
    function(from) ranges(from)
)

setAs("GAlignmentPairs", "GAlignmentsList", 
    function(from) 
    {
        if (length(from) == 0L)
            pbe <- PartitioningByEnd()
        else
            pbe <- PartitioningByEnd(seq(2, 2*length(from), 2), names=names(from)) 
        new("GAlignmentsList",
            unlistData=unlist(from, use.names=FALSE),
            partitioning=pbe)
        }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

## "[", "[<-" and "[[", "[[<-" from CompressedList


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going from GAlignments to GAlignmentsList with extractList() and family.
###

setMethod("relistToClass", "GAlignments",
    function(x) "GAlignmentsList"
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

setMethod("show", "GAlignmentsList",
    function(object)
        GenomicRanges:::showList(object, showGAlignments, FALSE)
)

