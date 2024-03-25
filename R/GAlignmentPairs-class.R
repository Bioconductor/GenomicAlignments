### =========================================================================
### GAlignmentPairs objects
### -------------------------------------------------------------------------
###

### "first" and "last" GAlignments must have identical seqinfo.
setClass("GAlignmentPairs",
    contains="List",
    representation(
        strandMode="integer",         # single integer (0L, 1L, or 2L)
        first="GAlignments",          # of length N, no names, no elt metadata
        last="GAlignments",           # of length N, no names, no elt metadata
        isProperPair="logical",       # of length N
        NAMES="character_OR_NULL",    # R doesn't like @names !!
        elementMetadata="DataFrame"   # N rows
    ),
    prototype(
        strandMode=1L,
        elementType="GAlignments"
    )
)

### Combine the new "parallel slots" with those of the parent class. Make
### sure to put the new parallel slots **first**. See R/Vector-class.R file
### in the S4Vectors package for what slots should or should not be considered
### "parallel".
setMethod("parallel_slot_names", "GAlignmentPairs",
    function(x) c("first", "last", "isProperPair", "NAMES", callNextMethod())
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
### updateObject()
###

setMethod("updateObject", "GAlignmentPairs",
    function(object, ..., verbose=FALSE)
    {
        object@first <- updateObject(object@first, ..., verbose=verbose)
        object@last  <- updateObject(object@last,  ..., verbose=verbose)

        callNextMethod()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("strandMode", "GAlignmentPairs",
    function(x) x@strandMode
)

setMethod("names", "GAlignmentPairs",
    function(x) x@NAMES
)

.make_it_real <- function(strand, strand_mode)
{
    if (strand_mode == 0L) 
        return(Rle(strand("*"), length(strand)))
    invertStrand(strand)
}

.first_strand <- function(x, real.strand=FALSE)
{
    if (!isTRUEorFALSE(real.strand))
        stop("'real.strand' must be TRUE or FALSE")
    ans <- strand(x@first)
    if (real.strand && strandMode(x) != 1L)
        ans <- .make_it_real(ans, strandMode(x))
    ans
}

.last_strand <- function(x, real.strand=FALSE)
{
    if (!isTRUEorFALSE(real.strand))
        stop("'real.strand' must be TRUE or FALSE")
    ans <- strand(x@last)
    if (real.strand && strandMode(x) != 2L)
        ans <- .make_it_real(ans, strandMode(x))
    ans
}

setMethod("first", "GAlignmentPairs",
    function(x, real.strand=FALSE)
    {
        if (!isTRUEorFALSE(real.strand))
            stop("'real.strand' must be TRUE or FALSE")
        ans <- x@first
        x_names <- names(x)
        if (!is.null(x_names))
            ans <- setNames(ans, x_names)
        if (real.strand && strandMode(x) != 1L)
            strand(ans) <- .first_strand(x, real.strand=TRUE)
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
        ans <- x@last
        x_names <- names(x)
        if (!is.null(x_names))
            ans <- setNames(ans, x_names)
        if (real.strand && strandMode(x) != 2L)
            strand(ans) <- .last_strand(x, real.strand=TRUE)
        ans
    }
)

setMethod("seqnames", "GAlignmentPairs",
    function(x)
    {
        ans <- seqnames(x@first)
        ans[seqnames(x@last) != ans] <- NA
        ans
    }
)

setMethod("strand", "GAlignmentPairs",
    function(x)
    {
        ans <- .first_strand(x, real.strand=TRUE)
        x_last_strand <- .last_strand(x, real.strand=TRUE)
        ans[ans != x_last_strand] <- "*"
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

setMethod("invertStrand", "GAlignmentPairs",
    function(x)
    {
        strand_mode <- strandMode(x)
        if (strand_mode != 0L)
            strandMode(x) <- 3L - strand_mode
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
    function(x, new2old=NULL,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
    {
        if (!is(value, "Seqinfo"))
            stop("the supplied 'seqinfo' must be a Seqinfo object")
        dangling_seqlevels <- GenomeInfoDb:::getDanglingSeqlevels(x,
                                  new2old=new2old,
                                  pruning.mode=pruning.mode,
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

.valid.GAlignmentPairs.isProperPair <- function(x)
{
    x_isProperPair <- x@isProperPair
    if (!is.logical(x_isProperPair) || !is.null(attributes(x_isProperPair))) {
        msg <- c("'x@isProperPair' must be a logical vector ",
                 "with no attributes")
        return(paste(msg, collapse=""))
    }
    if (S4Vectors:::anyMissing(x_isProperPair))
        return("'x@isProperPair' cannot contain NAs")
    NULL
}

.valid.GAlignmentPairs.seqinfo <- function(x)
{
    if (!identical(seqinfo(x@first), seqinfo(x@last)))
        return("'seqinfo(x@first)' and 'seqinfo(x@last)' must be identical")
    NULL
}

.valid.GAlignmentPairs <- function(x)
{
    c(.valid.GAlignmentPairs.strandMode(x),
      .valid.GAlignmentPairs.isProperPair(x),
      .valid.GAlignmentPairs.seqinfo(x))
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

### 'caller' must be "ranges" or "granges".
.error_on_discordant_seqnames <- function(caller)
{
    if (caller == "ranges") {
        range_type <- "range"
        returned_object <- "IRanges"
        a_returned_object <- "an IRanges"
        alternate_caller <- "rglist"
        alternate_returned_object <- "an IRangesList"
    } else {
        range_type <- "genomic range"
        returned_object <- "GRanges"
        a_returned_object <- "a GRanges"
        alternate_caller <- "grglist"
        alternate_returned_object <- "a GRangesList"
    }
    wmsg(
        "For some pairs in 'x', the 2 alignments are not on the same ",
        "chromosome. Cannot associate a unique ", range_type, " to ",
        "such pairs. Please call ", caller, "() with ",
        "'on.discordant.seqnames=\"drop\"' to drop these pairs, ",
        "or with 'on.discordant.seqnames=\"split\"' to represent ",
        "each of them with 2 ", range_type, "s in the returned ",
        returned_object, " object. Note that in both cases the returned ",
        "object won't be parallel to 'x'. Alternatively, please ",
        "consider using ", alternate_caller, "() instead of ", caller, "() ",
        "to turn 'x' into ", alternate_returned_object, " object instead of ",
        a_returned_object, " object. See ?GAlignmentPairs for more ",
        "information."
    )
}

.make_split_IRanges_from_GAlignmentPairs <- function(x, use.names=TRUE,
                                                        use.mcols=FALSE)
{
    x_first_ranges <- ranges(x@first, use.names=FALSE)
    x_last_ranges <- ranges(x@last, use.names=FALSE)

    x_len <- length(x)
    collate_subscript <- S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(x_len)
    partitioning <- PartitioningByEnd(2L * seq_len(x_len))

    range_pairs <- relist(c(x_first_ranges, x_last_ranges)[collate_subscript],
                          partitioning)

    x_first_seqnames <- seqnames(x@first)
    x_last_seqnames <- seqnames(x@last)
    merge_idx <- x_first_seqnames == x_last_seqnames

    range_pairs[merge_idx] <- range(range_pairs[merge_idx])
    if (use.names)
        names(range_pairs) <- names(x)

    ans <- unlist(range_pairs, use.names=use.names)
    if (use.mcols) {
        i <- rep.int(seq_len(x_len), elementNROWS(range_pairs))
        mcols(ans) <- extractROWS(mcols(x, use.names=FALSE), i)
    }
    ans
}

setMethod("ranges", "GAlignmentPairs",
    function(x, use.names=TRUE, use.mcols=FALSE,
                on.discordant.seqnames=c("error", "drop", "split"))
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        on.discordant.seqnames <- match.arg(on.discordant.seqnames)

        ans_seqnames <- seqnames(x)
        is_discordant <- is.na(ans_seqnames)
        if (any(is_discordant)) {
            if (on.discordant.seqnames == "error")
                stop(.error_on_discordant_seqnames("ranges"))
            if (on.discordant.seqnames == "split")
                return(.make_split_IRanges_from_GAlignmentPairs(x,
                           use.names=use.names, use.mcols=use.mcols))
            ## on.discordant.seqnames == "drop"
            x <- x[!is_discordant]
        }

        x_first_ranges <- ranges(x@first, use.names=FALSE)
        x_last_ranges <- ranges(x@last, use.names=FALSE)
        ans <- punion(x_first_ranges, x_last_ranges, fill.gap=TRUE)
        if (use.names)
            names(ans) <- names(x)
        if (use.mcols)
            mcols(ans) <- mcols(x, use.names=FALSE)
        ans
    }
)

.make_split_GRanges_from_GAlignmentPairs <- function(x, use.names=TRUE,
                                                        use.mcols=FALSE)
{
    x_first_seqnames <- seqnames(x@first)
    x_last_seqnames <- seqnames(x@last)
    is_discordant <- x_first_seqnames != x_last_seqnames
    ndiscordant <- sum(is_discordant)
    collate_subscript <-
        S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(ndiscordant)
    partitioning <- PartitioningByEnd(2L * seq_len(ndiscordant))

    seqnames1 <- x_first_seqnames[is_discordant]
    seqnames2 <- x_last_seqnames[is_discordant]
    discordant_seqnames <- relist(c(seqnames1, seqnames2)[collate_subscript],
                                  partitioning)

    ans_seqnames <- as(x_first_seqnames, "List")
    ans_seqnames[is_discordant] <- discordant_seqnames

    ans_ranges <- .make_split_IRanges_from_GAlignmentPairs(x,
                      use.names=use.names, use.mcols=use.mcols)

    strand1 <- .first_strand(x, real.strand=TRUE)[is_discordant]
    strand2 <- .last_strand(x, real.strand=TRUE)[is_discordant]
    discordant_strand <- relist(c(strand1, strand2)[collate_subscript],
                                partitioning)
    ans_strand <- as(strand(x), "List")
    ans_strand[is_discordant] <- discordant_strand

    GRanges(unlist(ans_seqnames, use.names=FALSE),
            ans_ranges,
            unlist(ans_strand, use.names=FALSE),
            seqinfo=seqinfo(x))
}

setMethod("granges", "GAlignmentPairs",
    function(x, use.names=TRUE, use.mcols=FALSE,
                on.discordant.seqnames=c("error", "drop", "split"))
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        on.discordant.seqnames <- match.arg(on.discordant.seqnames)

        ans_seqnames <- seqnames(x)
        is_discordant <- is.na(ans_seqnames)
        if (any(is_discordant)) {
            if (on.discordant.seqnames == "error")
                stop(.error_on_discordant_seqnames("granges"))
            if (on.discordant.seqnames == "split")
                return(.make_split_GRanges_from_GAlignmentPairs(x,
                           use.names=use.names, use.mcols=use.mcols))
            ## on.discordant.seqnames == "drop"
            x <- x[!is_discordant]
            ans_seqnames <- seqnames(x)
        }

        ans <- GRanges(ans_seqnames,
                       ranges(x, use.names=use.names),
                       strand(x),
                       seqinfo=seqinfo(x))
        if (use.mcols)
            mcols(ans) <- mcols(x, use.names=FALSE)
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

setMethod("grglist", "GAlignmentPairs",
    function(x, use.names=TRUE, use.mcols=FALSE, drop.D.ranges=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        x_mcols <- mcols(x, use.names=FALSE)
        if (use.mcols && "query.break" %in% colnames(x_mcols))
            stop("'mcols(x)' cannot have reserved column \"query.break\"")
        x_first <- x@first
        x_last <- x@last
        if (strandMode(x) == 0L) {
            x_first <- unstrand(x_first)
            x_last <- unstrand(x_last)
            x_unlisted <- c(x_first, x_last)
        } else if (strandMode(x) == 1L) {
            x_last <- invertStrand(x@last)
            x_unlisted <- c(x_first, x_last)
        } else if (strandMode(x) == 2L) {
            x_first <- invertStrand(x@first)
            x_unlisted <- c(x_last, x_first)
        } else {
            ## Should never happen.
            stop("unsupported strandMode: ", strandMode(x))
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
        ans_mcols <- DataFrame(query.break=mcols(ans, use.names=FALSE)$nelt1)
        if (use.mcols)
            ans_mcols <- cbind(ans_mcols, x_mcols)
        mcols(ans) <- ans_mcols
        ans
    }
)

setAs("GAlignmentPairs", "IntegerRanges",
    function(from) ranges(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignmentPairs", "GRanges",
    function(from) granges(from, use.names=TRUE, use.mcols=TRUE)
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
          DataFrame(firstDF, lastDF, mcols(from, use.names=FALSE))
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
    query.breaks <- mcols(x, use.names=FALSE)$query.break
    if (is.null(query.breaks))
        stop("'x' must be a GRangesList object with a \"query.breaks\" ",
             "metadata column")
    idx <- sequence(query.breaks, from=start(x@partitioning))
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
    nc <- ncol(mcols(x, use.names=FALSE))
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
        tmp <- do.call(data.frame,
                       lapply(mcols(x, use.names=FALSE), showAsCell))
        ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
    }
    ans
}

showGAlignmentPairs <- function(x, margin="",
                                   print.classinfo=FALSE,
                                   print.seqinfo=FALSE)
{
    lx <- length(x)
    nc <- ncol(mcols(x, use.names=FALSE))
    cat(class(x), " object with ",
        lx, " ", ifelse(lx == 1L, "pair", "pairs"),
        ", strandMode=", strandMode(x),
        ", and ",
        nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
        ":\n", sep="")
    out <- makePrettyMatrixForCompactPrinting(x,
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
        classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
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
### Concatenation
###

setMethod("bindROWS", "GAlignmentPairs",
    GenomicRanges:::concatenate_GenomicRanges_objects
)

