### =========================================================================
### GAlignments objects
### -------------------------------------------------------------------------
###

setClass("GAlignments",
    contains="Vector",
    representation(
        NAMES="character_OR_NULL",    # R doesn't like @names !!
        seqnames="Rle",               # 'factor' Rle
        start="integer",              # POS field in SAM
        cigar="character",            # extended CIGAR (see SAM format specs)
        strand="Rle",                 # 'factor' Rle
        #mismatches="character_OR_NULL", # see MD optional field in SAM format
                                         #specs
        elementMetadata="DataFrame",
        seqinfo="Seqinfo"
    ),
    prototype(
        seqnames=Rle(factor()),
        strand=Rle(strand())
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
### Getters.
###

setMethod("length", "GAlignments", function(x) length(x@cigar))

setMethod("names", "GAlignments", function(x) x@NAMES)

setMethod("seqnames", "GAlignments", function(x) x@seqnames)
setMethod("rname", "GAlignments", function(x) seqnames(x))

setMethod("cigar", "GAlignments", function(x) x@cigar)

setMethod("width", "GAlignments",
    function(x) cigarWidthAlongReferenceSpace(x@cigar)
)

setMethod("start", "GAlignments", function(x, ...) x@start)
setMethod("end", "GAlignments", function(x, ...) {x@start + width(x) - 1L})

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
        x@strand <-
            GenomicRanges:::normargGenomicRangesStrand(value, length(x))
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

.valid.GAlignments.names <- function(x)
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

.valid.GAlignments.seqnames <- function(x)
{
    if (length(seqnames(x)) != length(cigar(x)))
        return("'seqnames(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.GAlignments.start <- function(x)
{
    x_start <- start(x)
    if (!is.integer(x_start) || !is.null(names(x_start)) || S4Vectors:::anyMissing(x_start))
        return("'start(x)' must be an unnamed integer vector with no NAs")
    if (length(x_start) != length(cigar(x)))
        return("'start(x)' and 'cigar(x)' must have the same length")
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
    if (length(x_strand) != length(cigar(x)))
        return("'strand(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.GAlignments <- function(x)
{
    c(.valid.GAlignments.names(x),
      GenomicRanges:::.valid.GenomicRanges.seqnames(x),
      .valid.GAlignments.seqnames(x),
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

setMethod("updateObject", "GAlignments",
    function(object, ..., verbose=FALSE)
    {
        if (verbose)
            message("updateObject(object = 'GAlignments')")
        if (is(try(object@NAMES, silent=TRUE), "try-error")) {
            object@NAMES <- NULL
            return(object)
        }
        if (is(try(validObject(object@seqinfo, complete=TRUE), silent=TRUE),
               "try-error")) {
            object@seqinfo <- updateObject(object@seqinfo)
            return(object)
        }
        object
    }
)

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
    mcols(ans) <- mcols(x)
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
            mcols(ans) <- mcols(x)
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
            mcols(ans) <- mcols(x)
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
            mcols(ans) <- mcols(x)
        ans
    }
)

setAs("GAlignments", "Ranges",
    function(from) ranges(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignments", "GRanges",
    function(from) granges(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignments", "GRangesList",
    function(from) grglist(from, use.names=TRUE, use.mcols=TRUE)
)
setAs("GAlignments", "RangesList",
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
                    mcols(from),
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
      function(from) {
        ga <- GAlignments(seqnames(from), start(from),
                          if (!is.null(mcols(from)[["cigar"]]))
                            mcols(from)[["cigar"]]
                          else paste0(width(from), "M"),
                          strand(from),
                          if (!is.null(names(from))) names(from)
                          else mcols(from)$name,
                          seqlengths(from),
                          mcols(from)[setdiff(colnames(mcols(from)),
                                              c("cigar", "name"))])
        metadata(ga) <- metadata(from)
        seqinfo(ga) <- seqinfo(from)
        ga
      })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("extractROWS", "GAlignments",
    function(x, i)
    {
        i <- normalizeSingleBracketSubscript(i, x, as.NSBS=TRUE)
        ans_NAMES <- extractROWS(x@NAMES, i)
        ans_seqnames <- extractROWS(x@seqnames, i)
        ans_start <- extractROWS(x@start, i)
        ans_cigar <- extractROWS(x@cigar, i)
        ans_strand <- extractROWS(x@strand, i)
        ans_elementMetadata <- extractROWS(x@elementMetadata, i)
        BiocGenerics:::replaceSlots(x,
            NAMES=ans_NAMES,
            seqnames=ans_seqnames,
            start=ans_start,
            cigar=ans_cigar,
            strand=ans_strand,
            elementMetadata=ans_elementMetadata)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

.makeNakedMatFromGAlignments <- function(x)
{
    lx <- length(x)
    nc <- ncol(mcols(x))
    ans <- cbind(seqnames=as.character(seqnames(x)),
                 strand=as.character(strand(x)),
                 cigar=S4Vectors:::sketchStr(cigar(x), 23L),
                 qwidth=qwidth(x),
                 start=start(x),
                 end=end(x),
                 width=width(x),
                 njunc=njunc(x))
    if (nc > 0L) {
        tmp <- do.call(data.frame, lapply(mcols(x), showAsCell))
        ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
    }
    ans
}

showGAlignments <- function(x, margin="",
                            print.classinfo=FALSE, print.seqinfo=FALSE)
{
    lx <- length(x)
    nc <- ncol(mcols(x))
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
### Combining and splitting.
###

### 'Class' must be "GAlignments" or the name of a concrete subclass of
### GAlignments.
### 'objects' must be a list of GAlignments objects.
### Returns an instance of class 'Class'.
combine_GAlignments_objects <- function(Class, objects,
                                        use.names=TRUE, ignore.mcols=FALSE)
{
    if (!isSingleString(Class))
        stop("'Class' must be a single character string")
    if (!extends(Class, "GAlignments"))
        stop("'Class' must be the name of a class that extends GAlignments")
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

    ## Combine "seqnames" slots.
    seqnames_slots <- lapply(objects, function(x) x@seqnames)
    ## TODO: Implement unlist_list_of_Rle() in IRanges and use it here.
    ans_seqnames <- do.call(c, seqnames_slots)

    ## Combine "start" slots.
    start_slots <- lapply(objects, function(x) x@start)
    ans_start <- unlist(start_slots, use.names=FALSE)

    ## Combine "cigar" slots.
    cigar_slots <- lapply(objects, function(x) x@cigar)
    ans_cigar <- unlist(cigar_slots, use.names=FALSE)

    ## Combine "strand" slots.
    strand_slots <- lapply(objects, function(x) x@strand)
    ## TODO: Implement unlist_list_of_Rle() in IRanges and use it here.
    ans_strand <- do.call(c, strand_slots)

    ## Combine "mcols" slots. We don't need to use fancy
    ## IRanges:::rbind.mcols() for this because the "mcols" slot of a
    ## GAlignments object is guaranteed to be a DataFrame.
    if (ignore.mcols) {
        ans_mcols <- S4Vectors:::make_zero_col_DataFrame(length(ans_start))
    } else  {
        mcols_slots <- lapply(objects, function(x) x@elementMetadata)
        ## Will fail if not all the GAlignments objects in 'objects' have
        ## exactly the same metadata cols.
        ans_mcols <- do.call(rbind, mcols_slots)
    }

    ## Combine "seqinfo" slots.
    seqinfo_slots <- lapply(objects, function(x) x@seqinfo)
    ans_seqinfo <- do.call(merge, seqinfo_slots)

    ## Make 'ans' and return it.
    new(Class, NAMES=ans_NAMES,
               seqnames=ans_seqnames,
               start=ans_start,
               cigar=ans_cigar,
               strand=ans_strand,
               elementMetadata=ans_mcols,
               seqinfo=ans_seqinfo)
}

setMethod("c", "GAlignments",
    function (x, ..., ignore.mcols=FALSE, recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop("\"c\" method for GAlignments objects ",
                 "does not support the 'recursive' argument")
        if (missing(x)) {
            objects <- list(...)
            x <- objects[[1L]]
        } else {
            objects <- list(x, ...)
        }
        combine_GAlignments_objects(class(x), objects,
                                    use.names=FALSE,
                                    ignore.mcols=ignore.mcols)
    }
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
        oo <- GenomicRanges:::order_GenomicRanges(x, decreasing=decreasing,
                                                      ...)
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

