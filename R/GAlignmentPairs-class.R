### =========================================================================
### GAlignmentPairs objects
### -------------------------------------------------------------------------
###

### TODO: Implement a GAlignmentsList class (CompressedList subclass)
### and derive GAlignmentPairs from it.

### "first" and "last" GAlignments must have identical seqinfo.
setClass("GAlignmentPairs",
    contains="List",
    representation(
        strandMode="integer",         # single integer (0L, 1L, or 2L)
        NAMES="characterORNULL",      # R doesn't like @names !!
        first="GAlignments",          # of length N, no names, no elt metadata
        last="GAlignments",           # of length N, no names, no elt metadata
        isProperPair="logical",       # of length N
        elementMetadata="DataFrame"   # N rows
    ),
    prototype(
        strandMode=1L,
        elementType="GAlignments"
    )
)

### Formal API:
###   strandMode(x) - indicates how to infer the strand of a pair from the
###                 strand of the first and last alignments in the pair:
###                   0: strand of the pair is always *;
###                   1: strand of the pair is strand of its first alignment;
###                   2: strand of the pair is strand of its last alignment.
###                 These modes are equivalent to 'strandSpecific' equal 0, 1,
###                 and 2, respectively, for the featureCounts() function
###                 defined in the Rsubread package.
###   length(x)   - single integer N. Nb of pairs in 'x'.
###   names(x)    - NULL or character vector.
###   first(x)    - returns "first" slot.
###   last(x)     - returns "last" slot.
###   seqnames(x) - same as 'seqnames(first(x))' or 'seqnames(last(x))'.
###   strand(x)   - obeys strandMode(x) (see above).
###   njunc(x)    - same as 'njunc(first(x)) + njunc(last(x))'.
###   isProperPair(x) - returns "isProperPair" slot.
###   seqinfo(x)  - returns 'seqinfo(first(x))' (same as 'seqinfo(last(x))').
###   granges(x)  - GRanges object of the same length as 'x'.
###   grglist(x)  - GRangesList object of the same length as 'x'.
###   show(x)     - compact display in a data.frame-like fashion.
###   GAlignmentPairs(x) - constructor.
###   x[i]        - GAlignmentPairs object of the same class as 'x'
###                 (endomorphism).
###

setGeneric("strandMode", function(x) standardGeneric("strandMode"))
setGeneric("strandMode<-", signature="x",
    function(x, value) standardGeneric("strandMode<-")
)

setGeneric("isProperPair", function(x) standardGeneric("isProperPair"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("strandMode", "GAlignmentPairs",
    function(x) x@strandMode
)

setMethod("length", "GAlignmentPairs",
    function(x) length(x@first)
)

setMethod("names", "GAlignmentPairs",
    function(x) x@NAMES
)

setMethod("first", "GAlignmentPairs",
    function(x, real.strand=FALSE)
    {
        if (!isTRUEorFALSE(real.strand))
            stop("'real.strand' must be TRUE or FALSE")
        ans <- setNames(x@first, names(x))
        if (real.strand) {
            if (strandMode(x) == 0L) {
                strand(ans) <- "*"
            } else if (strandMode(x) == 2L) {
                ans <- invertStrand(ans)
            }
        }
        ans
    }
)

setGeneric("last", function(x, ...) standardGeneric("last"))

setMethod("last", "GAlignmentPairs", function(x, ...) second(x, ...))

setMethod("second", "GAlignmentPairs",
    function(x, real.strand=FALSE)
    {
        if (!isTRUEorFALSE(real.strand))
            stop("'real.strand' must be TRUE or FALSE")
        ans <- setNames(x@last, names(x))
        if (real.strand) {
            if (strandMode(x) == 0L) {
                strand(ans) <- "*"
            } else if (strandMode(x) == 1L) {
                ans <- invertStrand(ans)
            }
        }
        ans
    }
)

setMethod("seqnames", "GAlignmentPairs",
    function(x)
    {
        ans <- seqnames(x@first)
        if (any(ans != seqnames(x@last)))
            stop(wmsg("For some pairs in 'x', the 2 alignments are not on ",
                      "the same chromosome. Cannot associate a sequence name ",
                      "to them. ",
                      "Note that the GAlignmentPairs container only supports ",
                      "pairs where the 2 alignments are on the same ",
                      "chromosome at the moment."))
        ans
    }
)

setMethod("strand", "GAlignmentPairs",
    function(x)
    {
        if (strandMode(x) == 0L)
            return(strand(Rle("*", length(x))))
        x_first_strand <- strand(x@first)
        x_last_strand <- strand(x@last)
        if (strandMode(x) == 1L) {
            ans <- x_first_strand
        } else {
            ans <- x_last_strand
        }
        ans[x_first_strand == x_last_strand] <- "*"
        ans
    }
)

setMethod("njunc", "GAlignmentPairs",
    function(x) {njunc(x@first) + njunc(x@last)}
)

setMethod("isProperPair", "GAlignmentPairs",
    function(x) x@isProperPair
)

setMethod("seqinfo", "GAlignmentPairs",
    function(x) seqinfo(x@first)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

.normarg_strandMode_replace_value <- function(value)
{
    if (!isSingleNumber(value))
        stop("invalid strand mode (must be 0, 1, or 2)")
    if (!is.integer(value))
        value <- as.integer(value)
    if (!(value %in% 0:2))
        stop("invalid strand mode (must be 0, 1, or 2)")
    value
}

setReplaceMethod("strandMode", "GAlignmentPairs",
    function(x, value)
    {
        x@strandMode <- .normarg_strandMode_replace_value(value)
        x
    }
)

setReplaceMethod("names", "GAlignmentPairs",
    function(x, value)
    {
        if (!is.null(value))
            value <- as.character(value)
        x@NAMES <- value
        validObject(x)
        x
    }
)

setMethod("seqlevelsInUse", "GAlignmentPairs",
    function(x)
    {
        in_use1 <- seqlevelsInUse(x@first)
        in_use2 <- seqlevelsInUse(x@last)
        ## We cannot just do union() because we want the returned levels
        ## to be in the order they appear in 'seqlevels(x)'.
        intersect(seqlevels(x), union(in_use1, in_use2))
    }
)

setReplaceMethod("seqinfo", "GAlignmentPairs",
    function(x, new2old=NULL, force=FALSE, value)
    {
        if (!is(value, "Seqinfo"))
            stop("the supplied 'seqinfo' must be a Seqinfo object")
        dangling_seqlevels <- GenomeInfoDb:::getDanglingSeqlevels(x,
                                  new2old=new2old, force=force,
                                  seqlevels(value))
        if (length(dangling_seqlevels) != 0L) {
            dropme_in_first <- seqnames(x@first) %in% dangling_seqlevels
            dropme_in_last <- seqnames(x@last) %in% dangling_seqlevels
            dropme <- dropme_in_first | dropme_in_last
            x <- x[!dropme]
        }
        seqinfo(x@first, new2old=new2old) <- value
        seqinfo(x@last, new2old=new2old) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GAlignmentPairs.strandMode <- function(x)
{
    if (!(isSingleInteger(x@strandMode) && x@strandMode %in% 0:2))
        return("'x@strandMode' must be 0L, 1L, or 2L")
    NULL
}

.valid.GAlignmentPairs.names <- function(x)
{
    x_names <- names(x)
    if (is.null(x_names))
        return(NULL)
    if (!is.character(x_names) || !is.null(attributes(x_names))) {
        msg <- c("'names(x)' must be NULL or a character vector ",
                 "with no attributes")
        return(paste(msg, collapse=""))
    }
    if (length(x_names) != length(x))
        return("'names(x)' and 'x' must have the same length")
    NULL
}

.valid.GAlignmentPairs.first <- function(x)
{
    x_first <- x@first
    if (class(x_first) != "GAlignments")
        return("'x@first' must be a GAlignments instance")
    NULL
}

.valid.GAlignmentPairs.last <- function(x)
{
    x_last <- x@last
    if (class(x_last) != "GAlignments")
        return("'x@last' must be a GAlignments instance")
    x_first <- x@first
    if (length(x_last) != length(x_first))
        return("'x@last' and 'x@first' must have the same length")
    if (!identical(seqinfo(x_last), seqinfo(x_first)))
        return("'seqinfo(x@last)' and 'seqinfo(x@first)' must be identical")
    NULL
}

.valid.GAlignmentPairs.isProperPair <- function(x)
{
    x_isProperPair <- x@isProperPair
    if (!is.logical(x_isProperPair) || !is.null(attributes(x_isProperPair))) {
        msg <- c("'x@isProperPair' must be a logical vector ",
                 "with no attributes")
        return(paste(msg, collapse=""))
    }
    if (length(x_isProperPair) != length(x))
        return("'x@isProperPair' and 'x' must have the same length")
    if (S4Vectors:::anyMissing(x_isProperPair))
        return("'x@isProperPair' cannot contain NAs")
    NULL
}

.valid.GAlignmentPairs <- function(x)
{
    c(.valid.GAlignmentPairs.strandMode(x),
      .valid.GAlignmentPairs.names(x),
      .valid.GAlignmentPairs.first(x),
      .valid.GAlignmentPairs.last(x),
      .valid.GAlignmentPairs.isProperPair(x))
}

setValidity2("GAlignmentPairs", .valid.GAlignmentPairs,
             where=asNamespace("GenomicAlignments"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GAlignmentPairs <- function(first, last,
                            strandMode=1L, isProperPair=TRUE, names=NULL)
{
    if (!(is(first, "GAlignments") && is(last, "GAlignments")))
        stop("'first' and 'last' must be GAlignments objects")
    if (length(first) != length(last))
        stop("'first' and 'last' must have the same length")
    strandMode <- .normarg_strandMode_replace_value(strandMode)
    if (identical(isProperPair, TRUE))
        isProperPair <- rep.int(isProperPair, length(first))
    new2("GAlignmentPairs",
         strandMode=strandMode,
         NAMES=names,
         first=first, last=last,
         isProperPair=isProperPair,
         elementMetadata=S4Vectors:::make_zero_col_DataFrame(length(first)),
         check=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Vector methods.
###

setMethod("extractROWS", "GAlignmentPairs",
    function(x, i)
    {
        i <- normalizeSingleBracketSubscript(i, x, as.NSBS=TRUE)
        ans_NAMES <- extractROWS(x@NAMES, i)
        ans_first <- extractROWS(x@first, i)
        ans_last <- extractROWS(x@last, i)
        ans_isProperPair <- extractROWS(x@isProperPair, i)
        ans_elementMetadata <- extractROWS(x@elementMetadata, i)
        BiocGenerics:::replaceSlots(x,
            NAMES=ans_NAMES,
            first=ans_first,
            last=ans_last,
            isProperPair=ans_isProperPair,
            elementMetadata=ans_elementMetadata)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### List methods.
###

### TODO: Remove the "[[" method below after the definition of the
### GAlignmentPairs class is changed to derive from CompressedList.
### (The "[[" method for CompressedList objects should do just fine i.e. it
### should do something like x@unlistData[x@partitioning[[i]]] and that
### should be optimal.)
.GAlignmentPairs.getElement <- function(x, i)
{
    c(x@first[i], x@last[i])
}

setMethod("[[", "GAlignmentPairs",
    function(x, i, j, ... , drop=TRUE)
    {
        if (missing(i) || !missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        i <- normalizeDoubleBracketSubscript(i, x)
        .GAlignmentPairs.getElement(x, i)
    }
)

### TODO: Remove this method after the definition of the GAlignmentPairs
### class is changed to derive from CompressedList.
setMethod("unlist", "GAlignmentPairs",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        x_first <- x@first
        x_last <- x@last
        collate_subscript <-
            S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(length(x))
        ans <- c(x_first, x_last)[collate_subscript]
        if (use.names)
            names(ans) <- rep(names(x), each=2L)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setMethod("ranges", "GAlignmentPairs",
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        x_first_ranges <- ranges(x@first, use.names=FALSE)
        x_last_ranges <- ranges(x@last, use.names=FALSE)
        ans <- punion(x_first_ranges, x_last_ranges, fill.gap=TRUE)
        if (use.names)
            names(ans) <- names(x)
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)

setMethod("granges", "GAlignmentPairs",
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- GRanges(seqnames(x),
                       ranges(x, use.names=use.names),
                       strand(x),
                       seqinfo=seqinfo(x))
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)

### Shrink CompressedList 'x' (typically a GRangesList) by half by combining
### pairs of consecutive top-level elements.
shrinkByHalf <- function(x)
{
    if (length(x) %% 2L != 0L)
        stop("'x' must have an even length")
    x_eltNROWS <- elementNROWS(x)
    if (length(x_eltNROWS) == 0L) {
        ans_nelt1 <- ans_nelt2 <- integer(0)
    } else {
        ans_nelt1 <- x_eltNROWS[c(TRUE, FALSE)]
        ans_nelt2 <- x_eltNROWS[c(FALSE, TRUE)]
    }
    ans_eltNROWS <- ans_nelt1 + ans_nelt2
    ans_partitioning <- PartitioningByEnd(cumsum(ans_eltNROWS))
    ans <- relist(x@unlistData, ans_partitioning)
    mcols(ans) <- DataFrame(nelt1=ans_nelt1, nelt2=ans_nelt2)
    ans
}

### FIXME: Behavior is currently undefined (and undocumented) when
### strandMode(x) is 0. Fix this!
setMethod("grglist", "GAlignmentPairs",
    function(x, use.names=TRUE, use.mcols=FALSE, drop.D.ranges=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        x_mcols <- mcols(x)
        if (use.mcols && "query.break" %in% colnames(x_mcols))
            stop("'mcols(x)' cannot have reserved column \"query.break\"")
        x_first <- x@first
        x_last <- x@last
        if (strandMode(x) == 1L) {
            x_last <- invertStrand(x@last)
            x_unlisted <- c(x_first, x_last)
        } else if (strandMode(x) == 2L) {
            x_first <- invertStrand(x@first)
            x_unlisted <- c(x_last, x_first)
        }
        ## Not the same as doing 'unlist(x, use.names=FALSE)'.
        collate_subscript <-
            S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(length(x))
        x_unlisted <- x_unlisted[collate_subscript]
        grl <- grglist(x_unlisted,
                       order.as.in.query=TRUE,
                       drop.D.ranges=drop.D.ranges)
        ans <- shrinkByHalf(grl)
        if (use.names)
            names(ans) <- names(x)
        ans_mcols <- DataFrame(query.break=mcols(ans)$nelt1)
        if (use.mcols)
            ans_mcols <- cbind(ans_mcols, x_mcols)
        mcols(ans) <- ans_mcols
        ans
    }
)

setAs("GAlignmentPairs", "Ranges",
    function(from) ranges(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignmentPairs", "GRanges",
    function(from) granges(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignmentPairs", "GenomicRanges",
    function(from) as(from, "GRanges")
)
setAs("GAlignmentPairs", "GRangesList",
    function(from) grglist(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignmentPairs", "GAlignments",
    function(from) unlist(from, use.names=TRUE)
)

setAs("GAlignmentPairs", "DataFrame", function(from) {
          firstDF <- as(first(from), "DataFrame")
          colnames(firstDF) <- paste0(colnames(firstDF), ".first")
          lastDF <- as(last(from), "DataFrame")
          colnames(lastDF) <- paste0(colnames(lastDF), ".last")
          DataFrame(firstDF, lastDF, mcols(from))
      })

setMethod("as.data.frame", "GAlignmentPairs",
          function(x, row.names = NULL, optional = FALSE) {
              as.data.frame(as(x, "DataFrame"), row.names=row.names,
                            optional=optional)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fillJunctionGaps()
###
### Not exported. Used in the SplicingGraphs package.
###

fillJunctionGaps <- function(x)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    query.breaks <- mcols(x)$query.break
    if (is.null(query.breaks))
        stop("'x' must be a GRangesList object with a \"query.breaks\" ",
             "metadata column")
    offsets <- end(x@partitioning)
    if (length(x) != 0L) 
        offsets <- c(0L, offsets[-length(offsets)])
    idx <- S4Vectors:::fancy_mseq(query.breaks, offsets)
    half1_partitioning <- PartitioningByEnd(cumsum(query.breaks))
    half1 <- relist(x@unlistData[idx], half1_partitioning)
    half1 <- range(half1)@unlistData
    half2_eltNROWS <- elementNROWS(x) - query.breaks
    half2_partitioning <- PartitioningByEnd(cumsum(half2_eltNROWS))
    half2 <- relist(x@unlistData[-idx], half2_partitioning)
    half2 <- range(half2)@unlistData
    collate_subscript <- S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(length(x))
    ans_unlistData <- c(half1, half2)[collate_subscript]
    ans_partitioning <- PartitioningByEnd(2L * seq_along(x),
                                          names=names(x))
    relist(ans_unlistData, ans_partitioning)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

.makeNakedMatFromGAlignmentPairs <- function(x)
{
    lx <- length(x)
    nc <- ncol(mcols(x))
    pair_cols <- cbind(seqnames=as.character(seqnames(x)),
                       strand=as.character(strand(x)))
    x_first <- x@first
    first_cols <- cbind(ranges=showAsCell(ranges(x_first)))
    x_last <- x@last
    last_cols <- cbind(ranges=showAsCell(ranges(x_last)))
    ans <- cbind(pair_cols,
                 `:`=rep.int(":", lx),
                 first_cols,
                 `--`=rep.int("--", lx),
                 last_cols)
    if (nc > 0L) {
        tmp <- do.call(data.frame, lapply(mcols(x), showAsCell))
        ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
    }
    ans
}

showGAlignmentPairs <- function(x, margin="",
                                   print.classinfo=FALSE,
                                   print.seqinfo=FALSE)
{
    lx <- length(x)
    nc <- ncol(mcols(x))
    cat(class(x), " object with ",
        lx, " ", ifelse(lx == 1L, "pair", "pairs"),
        ", strandMode=", strandMode(x),
        ", and ",
        nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
        ":\n", sep="")
    out <- S4Vectors:::makePrettyMatrixForCompactPrinting(x,
               .makeNakedMatFromGAlignmentPairs)
    if (print.classinfo) {
        .PAIR_COL2CLASS <- c(
            seqnames="Rle",
            strand="Rle"
        )
        .HALVES_COL2CLASS <- c(
            ranges="IRanges"
        )
        .COL2CLASS <- c(.PAIR_COL2CLASS,
                        ":",
                        .HALVES_COL2CLASS,
                        "--",
                        .HALVES_COL2CLASS)
        classinfo <-
            S4Vectors:::makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        ## A sanity check, but this should never happen!
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste0(margin, rownames(out))
    ## We set 'max' to 'length(out)' to avoid the getOption("max.print")
    ## limit that would typically be reached when 'showHeadLines' global
    ## option is set to Inf.
    print(out, quote=FALSE, right=TRUE, max=length(out))
    if (print.seqinfo) {
        cat(margin, "-------\n", sep="")
        cat(margin, "seqinfo: ", summary(seqinfo(x)), "\n", sep="")
    }
}

setMethod("show", "GAlignmentPairs",
    function(object)
        showGAlignmentPairs(object, margin="  ",
                            print.classinfo=TRUE, print.seqinfo=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining.
###

### 'Class' must be "GAlignmentPairs" or the name of a concrete subclass of
### GAlignmentPairs.
### 'objects' must be a list of GAlignmentPairs objects.
### Returns an instance of class 'Class'.
combine_GAlignmentPairs_objects <- function(Class, objects,
                                            use.names=TRUE, ignore.mcols=FALSE)
{
    if (!isSingleString(Class))
        stop("'Class' must be a single character string")
    if (!extends(Class, "GAlignmentPairs"))
        stop("'Class' must be the name of a class that extends GAlignmentPairs")
    if (!is.list(objects))
        stop("'objects' must be a list")
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    ### TODO: Support 'use.names=TRUE'.
    if (use.names)
        stop("'use.names=TRUE' is not supported yet")
    if (!isTRUEorFALSE(ignore.mcols))
        stop("'ignore.mcols' must be TRUE or FALSE")

    if (length(objects) != 0L) {
        ## TODO: Implement (in C) fast 'elementIsNull(objects)' in IRanges,
        ## that does 'sapply(objects, is.null, USE.NAMES=FALSE)', and use it
        ## here.
        null_idx <- which(sapply(objects, is.null, USE.NAMES=FALSE))
        if (length(null_idx) != 0L)
            objects <- objects[-null_idx]
    }   
    if (length(objects) == 0L)
        return(new(Class))

    ## TODO: Implement (in C) fast 'elementIs(objects, class)' in IRanges, that
    ## does 'sapply(objects, is, class, USE.NAMES=FALSE)', and use it here.
    ## 'elementIs(objects, "NULL")' should work and be equivalent to
    ## 'elementIsNull(objects)'.
    if (!all(sapply(objects, is, Class, USE.NAMES=FALSE)))
        stop("the objects to combine must be ", Class, " objects (or NULLs)")
    objects_names <- names(objects)
    names(objects) <- NULL  # so lapply(objects, ...) below returns an
                            # unnamed list

    ## Combine "NAMES" slots.
    NAMES_slots <- lapply(objects, function(x) x@NAMES)
    ## TODO: Use elementIsNull() here when it becomes available.
    has_no_names <- sapply(NAMES_slots, is.null, USE.NAMES=FALSE)
    if (all(has_no_names)) {
        ans_NAMES <- NULL
    } else {
        noname_idx <- which(has_no_names)
        if (length(noname_idx) != 0L)
            NAMES_slots[noname_idx] <-
                lapply(elementNROWS(objects[noname_idx]), character)
        ans_NAMES <- unlist(NAMES_slots, use.names=FALSE)
    }

    ## Combine "first" slots.
    first_slots <- lapply(objects, function(x) x@first)
    ans_first <- combine_GAlignments_objects("GAlignments", first_slots,
                                             use.names=FALSE,
                                             ignore.mcols=ignore.mcols)

    ## Combine "last" slots.
    last_slots <- lapply(objects, function(x) x@last)
    ans_last <- combine_GAlignments_objects("GAlignments", last_slots,
                                            use.names=FALSE,
                                            ignore.mcols=ignore.mcols)

    ## Combine "isProperPair" slots.
    isProperPair_slots <- lapply(objects, function(x) x@isProperPair)
    ans_isProperPair <- unlist(isProperPair_slots, use.names=FALSE)

    ## Combine "mcols" slots. We don't need to use fancy
    ## IRanges:::rbind.mcols() for this because the "mcols" slot of a
    ## GAlignmentPairs object is guaranteed to be a DataFrame.
    if (ignore.mcols) {
        ans_mcols <- S4Vectors:::make_zero_col_DataFrame(length(ans_first))
    } else  {
        mcols_slots <- lapply(objects, function(x) x@elementMetadata)
        ## Will fail if not all the GAlignmentPairs objects in 'objects' have
        ## exactly the same metadata cols.
        ans_mcols <- do.call(rbind, mcols_slots)
    }

    ## Make 'ans' and return it.
    new(Class, NAMES=ans_NAMES,
               first=ans_first,
               last=ans_last,
               isProperPair=ans_isProperPair,
               elementMetadata=ans_mcols)
}

setMethod("c", "GAlignmentPairs",
    function(x, ..., ignore.mcols=FALSE, recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop("\"c\" method for GAlignmentPairs objects ",
                 "does not support the 'recursive' argument")
        if (missing(x)) {
            objects <- list(...)
            x <- objects[[1L]]
        } else {
            objects <- list(x, ...)
        }
        combine_GAlignmentPairs_objects(class(x), objects,
                                        use.names=FALSE, 
                                        ignore.mcols=ignore.mcols)
    }
)

