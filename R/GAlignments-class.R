### =========================================================================
### GAlignments objects
### -------------------------------------------------------------------------
###

setClass("GAlignments",
    contains="Ranges",
    representation(
        seqnames="Rle",               # 'factor' Rle
        start="integer",              # POS field in SAM
        cigar="character",            # extended CIGAR (see SAM format specs)
        strand="Rle",                 # 'factor' Rle
        #mismatches="character_OR_NULL", # see MD optional field in SAM format
                                         #specs
        NAMES="character_OR_NULL",    # R doesn't like @names !!
        elementMetadata="DataFrame",
        seqinfo="Seqinfo"
    ),
    prototype(
        seqnames=Rle(factor()),
        strand=Rle(strand())
    )
)

### Combine the new parallel slots with those of the parent class. Make sure
### to put the new parallel slots *first*.
setMethod("parallelSlotNames", "GAlignments",
    function(x) c("seqnames", "start", "cigar", "strand", "NAMES",
                  callNextMethod())
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
###   as.data.frame(x) - data.frame with 1 row per alignment in 'x'.
###   show(x)     - compact display in a data.frame-like fashion.
###   GAlignments(x) - constructor.
###   x[i]        - GAlignments object of the same class as 'x' (endomorphism).
###
###   qnarrow(x, start=NA, end=NA, width=NA) - GAlignments object of the
###                 same length and class as 'x' (endomorphism).
###
###   narrow(x, start=NA, end=NA, width=NA) - GAlignments object of the
###                 same length and class as 'x' (endomorphism).
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

setGeneric("rname", function(x) standardGeneric("rname"))
setGeneric("rname<-", function(x, value) standardGeneric("rname<-"))

setGeneric("cigar", function(x) standardGeneric("cigar"))

setGeneric("qwidth", function(x) standardGeneric("qwidth"))

setGeneric("njunc", function(x) standardGeneric("njunc"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###
### Internal representation of GAlignments objects has changed in
### GenomicAlignments 1.15.11 (Bioc 3.7).
###

.get_GAlignments_version <- function(object)
{
    if (.hasSlot(object, "elementType")) "current" else "< 1.15.11"
}

setMethod("updateObject", "GAlignments",
    function(object, ..., verbose=FALSE)
    {
        ## elementType slot.
        version <- .get_GAlignments_version(object)
        if (version == "current") {
            if (verbose)
                message("[updateObject] Internal representation of ",
                        class(object), " object is current.\n",
                        "[updateObject] Nothing to update.")
        } else {
            if (verbose)
                message("[updateObject] ", class(object), " object uses ",
                        "internal representation from\n",
                        "[updateObject] GenomicAlignments ", version, ". ",
                        "Updating it ...")
            object@elementType <- new(class(object))@elementType
        }
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("names", "GAlignments", function(x) x@NAMES)

setMethod("seqnames", "GAlignments", function(x) x@seqnames)
setMethod("rname", "GAlignments", function(x) seqnames(x))

setMethod("cigar", "GAlignments", function(x) x@cigar)

setMethod("width", "GAlignments",
    function(x) cigarWidthAlongReferenceSpace(x@cigar)
)

setMethod("start", "GAlignments", function(x, ...) x@start)

setMethod("strand", "GAlignments", function(x) x@strand)

setMethod("qwidth", "GAlignments",
    function(x) cigarWidthAlongQuerySpace(x@cigar)
)

setMethod("njunc", "GAlignments",
    function(x) {unname(elementNROWS(rglist(x))) - 1L}
)

setMethod("seqinfo", "GAlignments", function(x) x@seqinfo)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

setReplaceMethod("names", "GAlignments",
    function(x, value)
    {
        if (!is.null(value))
            value <- as.character(value)
        x@NAMES <- value
        validObject(x)
        x
    }
)

setReplaceMethod("seqnames", "GAlignments",
    function(x, value)
    {
        value <- GenomicRanges:::.normalize_seqnames_replacement_value(value, x)
        BiocGenerics:::replaceSlots(x, seqnames=value)
    }
)

setReplaceMethod("rname", "GAlignments",
    function(x, value) `seqnames<-`(x, value)
)

setReplaceMethod("strand", "GAlignments",
    function(x, value) 
    {
        x@strand <- GenomicRanges:::normalize_strand_replacement_value(value, x)
        x
    }
)

### Does NOT suppoprt pruning mode "fine". Pruning modes "coarse" and "tidy"
### are equivalent on a GAlignments object.
### FIXME: This repeats most of the code in
###        GenomicRanges:::set_GenomicRanges_seqinfo!
set_GAlignments_seqinfo <-
    function(x, new2old=NULL,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
{
    pruning.mode <- match.arg(pruning.mode)
    if (pruning.mode == "fine")
        stop(wmsg("\"fine\" pruning mode not supported on ",
                  class(x), " objects"))
    if (!is(value, "Seqinfo"))
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    dangling_seqlevels <- GenomeInfoDb:::getDanglingSeqlevels(x,
                              new2old=new2old,
                              pruning.mode=pruning.mode,
                              seqlevels(value))
    if (length(dangling_seqlevels) != 0L) {
        ## Prune 'x'.
        non_dangling_range <- !(seqnames(x) %in% dangling_seqlevels)
        x <- x[non_dangling_range]
    }
    old_seqinfo <- seqinfo(x)
    x@seqnames <- GenomeInfoDb:::makeNewSeqnames(x,
                              new2old, seqlevels(value))
    x@seqinfo <- value
    geom_has_changed <- GenomeInfoDb:::sequenceGeometryHasChanged(
                              seqinfo(x), old_seqinfo, new2old=new2old)
    if (any(geom_has_changed, na.rm=TRUE))
        validObject(x)
    x
}
setReplaceMethod("seqinfo", "GAlignments", set_GAlignments_seqinfo)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GAlignments.start <- function(x)
{
    x_start <- start(x)
    if (!is.integer(x_start) || !is.null(names(x_start)) || S4Vectors:::anyMissing(x_start))
        return("'start(x)' must be an unnamed integer vector with no NAs")
    NULL
}

.valid.GAlignments.cigar <- function(x)
{
    x_cigar <- cigar(x)
    if (!is.character(x_cigar) || !is.null(names(x_cigar)) || any(is.na(x_cigar)))
        return("'cigar(x)' must be an unnamed character vector with no NAs")
    tmp <- validCigar(x_cigar)
    if (!is.null(tmp))
        return(paste("in 'cigar(x)':", tmp))
    NULL
}

.valid.GAlignments.strand <- function(x)
{
    x_strand <- strand(x)
    if (!is(x_strand, "Rle") || !is.factor(runValue(x_strand))
     || !identical(levels(runValue(x_strand)), levels(strand()))
     || !is.null(names(x_strand)) || any(is.na(x_strand)))
        return("'strand(x)' must be an unnamed 'factor' Rle with no NAs (and with levels +, - and *)")
    NULL
}

.valid.GAlignments <- function(x)
{
    c(GenomicRanges:::.valid.GenomicRanges.seqnames(x),
      .valid.GAlignments.start(x),
      .valid.GAlignments.cigar(x),
      .valid.GAlignments.strand(x),
      GenomicRanges:::valid.GenomicRanges.seqinfo(x))
}

setValidity2("GAlignments", .valid.GAlignments,
             where=asNamespace("GenomicAlignments"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.asFactorRle <- function(x)
{
    if (is.character(x)) {
        x <- Rle(as.factor(x))
    } else if (is.factor(x)) {
        x <- Rle(x)
    } else if (is(x, "Rle") && is.character(runValue(x))) {
        runValue(x) <- as.factor(runValue(x))
    } else if (!is(x, "Rle") || !is.factor(runValue(x))) {
        stop("'x' must be a character vector, a factor, ",
             "a 'character' Rle, or a 'factor' Rle")
    }
    x
}

GAlignments <- function(seqnames=Rle(factor()), pos=integer(0),
                        cigar=character(0), strand=NULL,
                        names=NULL, seqlengths=NULL, ...)
{
    ## Prepare the 'seqnames' slot.
    seqnames <- .asFactorRle(seqnames)
    if (any(is.na(seqnames)))
        stop("'seqnames' cannot have NAs")
    ## Prepare the 'pos' slot.
    if (!is.integer(pos) || any(is.na(pos)))
        stop("'pos' must be an integer vector with no NAs")
    ## Prepare the 'cigar' slot.
    if (!is.character(cigar) || any(is.na(cigar)))
        stop("'cigar' must be a character vector with no NAs")
    ## Prepare the 'strand' slot.
    if (is.null(strand)) {
        if (length(seqnames) != 0L)
            stop("'strand' must be specified when 'seqnames' is not empty")
        strand <- Rle(strand())
    } else if (is.factor(strand)) {
        strand <- Rle(strand)
    }
    ## Prepare the 'elementMetadata' slot.
    varlist <- list(...)
    elementMetadata <- 
        if (0L == length(varlist))
            S4Vectors:::make_zero_col_DataFrame(length(seqnames))
        else
            do.call(DataFrame, varlist)
    ## Prepare the 'seqinfo' slot.
    if (is.null(seqlengths)) {
        seqlengths <- rep(NA_integer_, length(levels(seqnames)))
        names(seqlengths) <- levels(seqnames)
    } else if (!is.numeric(seqlengths)
            || is.null(names(seqlengths))
            || any(duplicated(names(seqlengths)))) {
        stop("'seqlengths' must be an integer vector with unique names")
    } else if (!setequal(names(seqlengths), levels(seqnames))) {
        stop("'names(seqlengths)' incompatible with 'levels(seqnames)'")
    } else if (!is.integer(seqlengths)) { 
        storage.mode(seqlengths) <- "integer"
    }
    seqinfo <- Seqinfo(seqnames=names(seqlengths), seqlengths=seqlengths)
    ## Create and return the GAlignments instance.
    new("GAlignments", NAMES=names,
                       seqnames=seqnames, start=pos, cigar=cigar,
                       strand=strand,
                       elementMetadata=elementMetadata,
                       seqinfo=seqinfo)
}

setMethod("update", "GAlignments",
          function(object, ...)
          {
            BiocGenerics:::replaceSlots(object, ...)
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Internal helper function used by higher level coercion functions.
###
### NOT exported.
###

### Names are propagated via 'x@partitioning' ('x' is a CompressedIRangesList).
make_GRangesList_from_CompressedIRangesList <- function(x, seqnames, strand,
                                                           seqinfo)
{
    x_eltNROWS <- elementNROWS(x)
    seqnames <- rep.int(seqnames, x_eltNROWS)
    strand <- rep.int(strand, x_eltNROWS)
    unlisted_ans <- GRanges(seqnames=seqnames, ranges=x@unlistData,
                            strand=strand)
    seqinfo(unlisted_ans) <- seqinfo
    ans <- relist(unlisted_ans, x)
    mcols(ans) <- mcols(x, use.names=FALSE)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setMethod("ranges", "GAlignments",
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        if (use.names) {
            ans_names <- names(x)
        } else {
            ans_names <- NULL
        }
        ans <- IRanges(start=start(x), width=width(x), names=ans_names)
        if (use.mcols)
            mcols(ans) <- mcols(x, use.names=FALSE)
        ans
    }
)

setMethod("granges", "GAlignments",
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- GRanges(seqnames(x),
                       ranges(x, use.names=use.names),
                       strand(x),
                       seqinfo=seqinfo(x))
        if (use.mcols)
            mcols(ans) <- mcols(x, use.names=FALSE)
        ans
    }
)

setMethod("grglist", "GAlignments",
    function(x, use.names=TRUE, use.mcols=FALSE,
                order.as.in.query=FALSE, drop.D.ranges=FALSE)
    {
        rgl <- rglist(x, use.names=use.names,
                         use.mcols=use.mcols,
                         order.as.in.query=order.as.in.query,
                         drop.D.ranges=drop.D.ranges)
        make_GRangesList_from_CompressedIRangesList(rgl,
                         seqnames(x), strand(x), seqinfo(x))
    }
)

setMethod("rglist", "GAlignments",
    function(x, use.names=TRUE, use.mcols=FALSE,
                order.as.in.query=FALSE, drop.D.ranges=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        if (!isTRUEorFALSE(order.as.in.query))
            stop("'order.as.in.query' must be TRUE or FALSE")
        ans <- extractAlignmentRangesOnReference(x@cigar, x@start,
                                                 drop.D.ranges=drop.D.ranges)
        if (order.as.in.query)
            ans <- revElements(ans, strand(x) == "-")
        if (use.names)
            names(ans) <- names(x)
        if (use.mcols)
            mcols(ans) <- mcols(x, use.names=FALSE)
        ans
    }
)

setAs("GAlignments", "IntegerRanges",
    function(from) ranges(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignments", "GRanges",
    function(from) granges(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignments", "GRangesList",
    function(from) grglist(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignments", "IntegerRangesList",
    function(from) rglist(from, use.names=TRUE, use.mcols=TRUE)
)

setAs("GAlignments", "DataFrame", function(from) {
          DataFrame(seqnames=seqnames(from),
                    strand=strand(from),
                    cigar=cigar(from),
                    qwidth=qwidth(from),
                    start=start(from),
                    end=end(from),
                    width=width(from),
                    njunc=njunc(from),
                    mcols(from, use.names=FALSE),
                    row.names=names(from),
                    check.names=FALSE)
      })

setMethod("as.data.frame", "GAlignments",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        as.data.frame(as(x, "DataFrame"), row.names=row.names,
                      optional=optional)
    }
)

setAs("GenomicRanges", "GAlignments",
    function(from)
    {
        from_mcols <- mcols(from, use.names=FALSE)
        ans_cigar <- from_mcols[["cigar"]]
        if (is.null(ans_cigar))
            ans_cigar <- paste0(width(from), "M")
        ans <- GAlignments(seqnames(from), start(from), ans_cigar, strand(from),
                  if (!is.null(names(from))) names(from) else from_mcols$name,
                  seqlengths(from),
                  from_mcols[setdiff(colnames(from_mcols), c("cigar", "name"))]
               )
        metadata(ans) <- metadata(from)
        seqinfo(ans) <- seqinfo(from)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

.makeNakedMatFromGAlignments <- function(x)
{
    lx <- length(x)
    nc <- ncol(mcols(x, use.names=FALSE))
    ans <- cbind(seqnames=as.character(seqnames(x)),
                 strand=as.character(strand(x)),
                 cigar=S4Vectors:::sketchStr(cigar(x), 23L),
                 qwidth=qwidth(x),
                 start=start(x),
                 end=end(x),
                 width=width(x),
                 njunc=njunc(x))
    if (nc > 0L) {
        tmp <- do.call(data.frame,
                       lapply(mcols(x, use.names=FALSE), showAsCell))
        ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
    }
    ans
}

showGAlignments <- function(x, margin="",
                            print.classinfo=FALSE, print.seqinfo=FALSE)
{
    lx <- length(x)
    nc <- ncol(mcols(x, use.names=FALSE))
    cat(class(x), " object with ",
        lx, " ", ifelse(lx == 1L, "alignment", "alignments"),
        " and ",
        nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
        ":\n", sep="")
    out <- S4Vectors:::makePrettyMatrixForCompactPrinting(x,
               .makeNakedMatFromGAlignments)
    if (print.classinfo) {
        .COL2CLASS <- c(
            seqnames="Rle",
            strand="Rle",
            cigar="character",
            qwidth="integer",
            start="integer",
            end="integer",
            width="integer",
            njunc="integer"
        )
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

setMethod("show", "GAlignments",
    function(object)
        showGAlignments(object, margin="  ",
                        print.classinfo=TRUE, print.seqinfo=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Concatenation
###

setMethod("bindROWS", "GAlignments",
    GenomicRanges:::concatenate_GenomicRanges_objects
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparing/ordering.
###

setMethod("pcompare", c("GAlignments", "GAlignments"),
    function(x, y)
    {
        x <- granges(x, use.names=FALSE)
        y <- granges(y, use.names=FALSE)
        callGeneric(x, y)
    }
)

setMethod("is.unsorted", "GAlignments",
    function(x, na.rm=FALSE, strictly=FALSE, ...) 
    {
        x <- granges(x, use.names=FALSE)
        callGeneric()
    }
)

setMethod("order", "GAlignments",
    function(..., na.last=TRUE, decreasing=FALSE,
                  method=c("auto", "shell", "radix"))
    {
        args <- list(...)
        order_args <- c(lapply(args, granges, use.names=FALSE),
                        list(na.last=na.last, decreasing=decreasing,
                             method=method))
        do.call(order, order_args)
    }
)

setMethod("sort", "GAlignments",
    function(x, decreasing=FALSE, ...)
    {
        oo <- GenomicRanges:::order_GenomicRanges(x, decreasing=decreasing, ...)
        extractROWS(x, oo)
    }
)

setMethod("rank", "GAlignments",
    function(x, na.last=TRUE,
             ties.method=c("average", "first", "last", "random", "max", "min"),
             ...)
    {
        x <- granges(x, use.names=FALSE)
        callGeneric()
    }
)

