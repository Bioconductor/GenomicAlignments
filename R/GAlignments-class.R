### =========================================================================
### GAlignments objects
### -------------------------------------------------------------------------
###

setClass("GAlignments",
    contains="Vector",
    representation(
        NAMES="characterORNULL",      # R doesn't like @names !!
        seqnames="Rle",               # 'factor' Rle
        start="integer",              # POS field in SAM
        cigar="character",            # extended CIGAR (see SAM format specs)
        strand="Rle",                 # 'factor' Rle
        #mismatches="characterORNULL", # see MD optional field in SAM format specs
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
    function(x) {unname(elementLengths(rglist(x))) - 1L}
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

.normargSeqnamesReplaceValue <- function(x, value, ans.type=c("factor", "Rle"))
{
    ans.type <- match.arg(ans.type)
    if (!is.factor(value)
     && !is.character(value)
     && (!is(value, "Rle") || !is.character(runValue(value))
                              && !is.factor(runValue(value))))
        stop("'seqnames' value must be a character factor/vector, ",
             "or a 'character' Rle, or a 'factor' Rle")
    if (ans.type == "factor") {
        if (!is.factor(value))
            value <- as.factor(value)
    } else if (ans.type == "Rle") {
        ## We want to return a 'factor' Rle.
        if (!is(value, "Rle")) {
            if (!is.factor(value))
                value <- as.factor(value)
            value <- Rle(value)
        } else if (!is.factor(runValue(value))) {
            runValue(value) <- as.factor(runValue(value))
        }
    }
    if (length(value) != length(x))
        stop("'seqnames' value must be the same length as the object")
    value
}

### 'old_seqnames' and 'new_seqnames' must be 'factor' Rle.
.getSeqnamesTranslationTable <- function(old_seqnames, new_seqnames)
{
    old <- runValue(old_seqnames)
    new <- runValue(new_seqnames)
    tmp <- unique(data.frame(old=old, new=new))
    if (!identical(runLength(old_seqnames), runLength(new_seqnames)) ||
        anyDuplicated(tmp$old) || anyDuplicated(tmp$new))
        stop("mapping between old an new 'seqnames' values is not one-to-one")
    if (isTRUE(all.equal(as.integer(tmp$old), as.integer(tmp$new)))) {
        tr_table <- levels(new)
        names(tr_table) <- levels(old)
    } else {
        tr_table <- tmp$new
        names(tr_table) <- tmp$old
    }
    tr_table
}

setReplaceMethod("seqnames", "GAlignments",
    function(x, value)
    {
        value <- .normargSeqnamesReplaceValue(x, value, ans.type="Rle")
        tr_table <- .getSeqnamesTranslationTable(seqnames(x), value)
        x@seqnames <- value
        seqnames(x@seqinfo) <- tr_table[seqlevels(x)]
        x
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

setReplaceMethod("elementMetadata", "GAlignments",
    function(x, ..., value)
    {
        value <-
            GenomicRanges:::normalizeMetadataColumnsReplacementValue(value, x)
        x@elementMetadata <- value
        x
    }
)

setReplaceMethod("seqinfo", "GAlignments",
    function(x, new2old=NULL, force=FALSE, value)
    {
        if (!is(value, "Seqinfo"))
            stop("the supplied 'seqinfo' must be a Seqinfo object")
        dangling_seqlevels <- GenomeInfoDb:::getDanglingSeqlevels(x,
                                  new2old=new2old, force=force,
                                  seqlevels(value))
        if (length(dangling_seqlevels) != 0L)
            x <- x[!(seqnames(x) %in% dangling_seqlevels)]
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
)


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
    x_seqnames <- seqnames(x)
    if (!is(x_seqnames, "Rle") || !is.factor(runValue(x_seqnames))
     || !is.null(names(x_seqnames)) || any(is.na(x_seqnames)))
        return("'seqnames(x)' must be an unnamed 'factor' Rle with no NAs")
    if (length(x_seqnames) != length(cigar(x)))
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
            new("DataFrame", nrows=length(seqnames))
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
            BiocGenerics:::updateS4(object, ...)
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Three helper functions used by higher level coercion functions.
###
### Note that their arguments are the different components of a
### GAlignments object instead of just the GAlignments object
### itself (arg 'x'). This allows them to be used in many different contexts
### e.g. when 'x' doesn't exist yet but is in the process of being constructed.
###

.GAlignmentsToGRanges <- function(seqnames, start, width, strand, seqinfo,
                                  names=NULL)
{
    ranges <- IRanges(start=start, width=width, names=names)
    ans <- GRanges(seqnames=seqnames, ranges=ranges, strand=strand)
    seqinfo(ans) <- seqinfo
    ans
}

### Names are propagated via 'x@partitioning' ('x' is a CompressedIRangesList).
.CompressedIRangesListToGRangesList <- function(x, seqnames, strand, seqinfo)
{
    elt_lens <- elementLengths(x)
    seqnames <- rep.int(seqnames, elt_lens)
    strand <- rep.int(strand, elt_lens)
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

setMethod("grglist", "GAlignments",
    function(x, use.mcols=FALSE, order.as.in.query=FALSE, drop.D.ranges=FALSE)
    {
        rgl <- rglist(x, use.mcols=use.mcols,
                         order.as.in.query=order.as.in.query,
                         drop.D.ranges=drop.D.ranges)
        .CompressedIRangesListToGRangesList(rgl, seqnames(x), strand(x),
                                            seqinfo(x))
    }
)

setMethod("granges", "GAlignments",
    function(x, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- .GAlignmentsToGRanges(seqnames(x), start(x), width(x),
                                     strand(x), seqinfo(x), names(x))
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)

setMethod("rglist", "GAlignments",
    function(x, use.mcols=FALSE, order.as.in.query=FALSE, drop.D.ranges=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        if (!isTRUEorFALSE(order.as.in.query))
            stop("'reorder.ranges.from5to3' must be TRUE or FALSE")
        ans <- extractAlignmentRangesOnReference(x@cigar, x@start,
                                                 drop.D.ranges=drop.D.ranges)
        if (order.as.in.query)
            ans <- revElements(ans, strand(x) == "-")
        names(ans) <- names(x)
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)

setMethod("ranges", "GAlignments",
    function(x) IRanges(start=start(x), width=width(x), names=names(x))
)

setAs("GAlignments", "GRangesList",
    function(from) grglist(from, use.mcols=TRUE)
)
setAs("GAlignments", "GRanges", function(from) granges(from, use.mcols=TRUE))
setAs("GAlignments", "RangesList", function(from) rglist(from, use.mcols=TRUE))
setAs("GAlignments", "Ranges", function(from) ranges(from))

setMethod("as.data.frame", "GAlignments",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (is.null(row.names))
            row.names <- names(x)
        else if (!is.character(row.names))
            stop("'row.names' must be NULL or a character vector")
        ans <- data.frame(seqnames=as.character(seqnames(x)),
                          strand=as.character(strand(x)),
                          cigar=cigar(x),
                          qwidth=qwidth(x),
                          start=start(x),
                          end=end(x),
                          width=width(x),
                          njunc=njunc(x),
                          row.names=row.names,
                          check.rows=TRUE,
                          check.names=FALSE,
                          stringsAsFactors=FALSE)
        if (ncol(mcols(x)))
            ans <- cbind(ans, as.data.frame(mcols(x)))
        return(ans)
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
        ans_names <- extractROWS(names(x), i)
        ans_seqnames <- extractROWS(seqnames(x), i)
        ans_start <- extractROWS(start(x), i)
        ans_cigar <- extractROWS(cigar(x), i)
        ans_strand <- extractROWS(strand(x), i)
        ans_mcols <- extractROWS(mcols(x), i)
        GenomicRanges:::clone(x, NAMES=ans_names,
                                 seqnames=ans_seqnames,
                                 start=ans_start,
                                 cigar=ans_cigar,
                                 strand=ans_strand,
                                 elementMetadata=ans_mcols)
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
                 cigar=cigar(x),
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
                            with.classinfo=FALSE, print.seqlengths=FALSE)
{
    lx <- length(x)
    nc <- ncol(mcols(x))
    cat(class(x), " with ",
        lx, " ", ifelse(lx == 1L, "alignment", "alignments"),
        " and ",
        nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
        ":\n", sep="")
    out <- S4Vectors:::makePrettyMatrixForCompactPrinting(x,
               .makeNakedMatFromGAlignments)
    if (with.classinfo) {
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
            GenomicRanges:::makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        ## A sanity check, but this should never happen!
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste0(margin, rownames(out))
    print(out, quote=FALSE, right=TRUE)
    if (print.seqlengths) {
        cat(margin, "---\n", sep="")
        GenomicRanges:::showSeqlengths(x, margin=margin)
    }
}

setMethod("show", "GAlignments",
    function(object)
        showGAlignments(object, margin="  ",
                        with.classinfo=TRUE, print.seqlengths=TRUE)
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
                lapply(elementLengths(objects[noname_idx]), character)
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
        ans_mcols <- new("DataFrame", nrows=length(ans_start))
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
### Old stuff (deprecated & defunct)
###

setMethod("ngap", "GAlignments",
    function(x)
    {
        .Defunct("njunc")
        njunc(x)
    }
)

