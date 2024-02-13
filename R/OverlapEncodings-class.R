### =========================================================================
### OverlapEncodings objects
### -------------------------------------------------------------------------
###


setClass("OverlapEncodings",
    contains="Vector",
    representation(
        Loffset="integer",      # no NAs, >= 0
        Roffset="integer",      # no NAs, >= 0
        encoding="factor",      # no NAs
        flippedQuery="logical"  # no NAs
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### parallel_slot_names()
###

### Combine the new "parallel slots" with those of the parent class. Make
### sure to put the new parallel slots **first**. See R/Vector-class.R file
### in the S4Vectors package for what slots should or should not be considered
### "parallel".
setMethod("parallel_slot_names", "OverlapEncodings",
    function(x) c("Loffset", "Roffset", "encoding", "flippedQuery",
                  callNextMethod())
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setGeneric("Loffset", function(x) standardGeneric("Loffset"))
setMethod("Loffset", "OverlapEncodings", function(x) x@Loffset)

setGeneric("Roffset", function(x) standardGeneric("Roffset"))
setMethod("Roffset", "OverlapEncodings", function(x) x@Roffset)

### encoding() generic is defined in Biostrings.
setMethod("encoding", "OverlapEncodings", function(x) x@encoding)

### S3/S4 combo for levels.OverlapEncodings
levels.OverlapEncodings <- function(x) levels(encoding(x))
setMethod("levels", "OverlapEncodings", levels.OverlapEncodings)

setGeneric("flippedQuery", function(x) standardGeneric("flippedQuery"))
setMethod("flippedQuery", "OverlapEncodings", function(x) x@flippedQuery)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The encodingHalves(), Lencoding() and Rencoding() low-level utilities.
###

.split_encoding_halves <- function(x, single.end.on.left=FALSE,
                                      single.end.on.right=FALSE,
                                      as.factors=FALSE)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    if (!isTRUEorFALSE(single.end.on.left))
        stop("'single.end.on.left' must be TRUE or FALSE")
    if (!isTRUEorFALSE(single.end.on.right))
        stop("'single.end.on.right' must be TRUE or FALSE")
    if (!isTRUEorFALSE(as.factors))
        stop("'as.factors' must be TRUE or FALSE")

    encoding_blocks <- CharacterList(strsplit(x, ":", fixed=TRUE))
    unlisted_blocks <- unlist(encoding_blocks, use.names=FALSE)
    block_halves <- CharacterList(strsplit(unlisted_blocks, "--", fixed=TRUE))

    ## Check that the blocks in any given encoding are either all single-end
    ## or all paired-end.
    nhalves <- unique(relist(elementNROWS(block_halves), encoding_blocks))
    if (any(elementNROWS(nhalves) != 1L))
        stop("some encodings are ill-formed")
    nhalves <- as.integer(nhalves)
    if (any(nhalves > 2L))
        stop("some encodings are ill-formed")

    halves2encoding <- function(halves, skeleton) {
        blocks <- relist(as.character(halves), skeleton)
        encoding <- unstrsplit(blocks, sep=":")
        if (length(encoding) != 0L)
            encoding <- setNames(paste0(encoding, ":"), names(encoding))
        encoding
    }

    Lencoding <- halves2encoding(heads(block_halves, n=1L), encoding_blocks)
    Rencoding <- halves2encoding(tails(block_halves, n=1L), encoding_blocks)
    if (!(single.end.on.left && single.end.on.right)) {
        idx <- which(nhalves == 1L)
        if (!single.end.on.left)
            Lencoding[idx] <- NA_character_
        if (!single.end.on.right)
            Rencoding[idx] <- NA_character_
    }
    if (as.factors) {
        Lencoding <- factor(Lencoding, levels=unique(Lencoding))
        Rencoding <- factor(Rencoding, levels=unique(Rencoding))
    }
    list(Lencoding, Rencoding)
}

setGeneric("encodingHalves", signature="x",
    function(x, single.end.on.left=FALSE, single.end.on.right=FALSE,
                as.factors=FALSE)
        standardGeneric("encodingHalves")
)

setMethod("encodingHalves", "character", .split_encoding_halves)

setMethod("encodingHalves", "factor",
    function(x, single.end.on.left=FALSE, single.end.on.right=FALSE,
                as.factors=FALSE)
    {
        levels_halves <- encodingHalves(levels(x),
                             single.end.on.left=single.end.on.left,
                             single.end.on.right=single.end.on.right,
                             as.factors=as.factors)
        x <- as.integer(x)
        list(levels_halves[[1L]][x], levels_halves[[2L]][x])
    }
)
setMethod("encodingHalves", "OverlapEncodings",
    function(x, single.end.on.left=FALSE, single.end.on.right=FALSE,
                as.factors=FALSE)
    {
        encodingHalves(encoding(x),
                       single.end.on.left=single.end.on.left,
                       single.end.on.right=single.end.on.right,
                       as.factors=as.factors)
    }
)

Lencoding <- function(x, ...) encodingHalves(x, ...)[[1L]]
Rencoding <- function(x, ...) encodingHalves(x, ...)[[2L]]


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The njunc(), Lnjunc(), and Rnjunc() low-level utilities.
###

### 'x' must be a character vector of single-end encodings.
.njunc_single_end_encodings <- function(x)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    old_warn <- getOption("warn")
    options(warn=2L)
    on.exit(options(warn=old_warn))
    M <- try(as.integer(sub(":.*", "", x)), silent=TRUE)
    if (inherits(M, "try-error"))
        stop("some encodings are ill-formed")
    options(warn=old_warn)
    if (any(M < 1L, na.rm=TRUE))
        warning(wmsg("some encodings start with a value < 1 and that is ",
                     "interpreted as a negative number of junctions)"))
    setNames(M - 1L, names(x))
}

Lnjunc <- function(x, single.end.on.left=FALSE)
{
    Lencoding <- Lencoding(x, single.end.on.left=single.end.on.left,
                              as.factors=TRUE)
    .njunc_single_end_encodings(levels(Lencoding))[as.integer(Lencoding)]
}

Rnjunc <- function(x, single.end.on.right=FALSE)
{
    Rencoding <- Rencoding(x, single.end.on.right=single.end.on.right,
                              as.factors=TRUE)
    .njunc_single_end_encodings(levels(Rencoding))[as.integer(Rencoding)]
}

### We make this the default "njunc" method although it will only work on
### objects supported by encodingHalves() (which is called behind the scene),
### that is, for character, factor, and OverlapEncodings objects. So instead
### of defining 3 "njunc" methods (one for each type of object supported by
### encodingHalves()), we define a single one for expediency. Another advantage
### of this approach is that if, in the future, encodingHalves() is extended to
### support more types of objects, then njunc() will work out-of-the-box on
### them i.e. with no need to define additional "njunc" methods.
setMethod("njunc", "ANY",
    function(x)
    {
        Lnjunc <- Lnjunc(x, single.end.on.left=TRUE)
        Rnjunc <- Rnjunc(x)
        Rnjunc[is.na(Rnjunc)] <- 0L
        Lnjunc + Rnjunc
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

### S3/S4 combo for as.data.frame.OverlapEncodings
as.data.frame.OverlapEncodings <- function(x, row.names=NULL,
                                           optional=FALSE, ...)
{
    if (!(is.null(row.names) || is.character(row.names)))
        stop("'row.names' must be NULL or a character vector")
    data.frame(Loffset=Loffset(x),
               Roffset=Roffset(x),
               encoding=encoding(x),
               flippedQuery=flippedQuery(x),
               row.names=row.names,
               check.rows=TRUE,
               check.names=FALSE,
               stringsAsFactors=FALSE)
}
setMethod("as.data.frame", "OverlapEncodings", as.data.frame.OverlapEncodings)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Displaying
###

.make_naked_matrix_from_OverlapEncodings <- function(x)
{
    x_len <- length(x)
    x_mcols <- mcols(x, use.names=FALSE)
    x_nmc <- if (is.null(x_mcols)) 0L else ncol(x_mcols)
    ans <- cbind(Loffset=as.character(Loffset(x)),
                 Roffset=as.character(Roffset(x)),
                 encoding=as.character(encoding(x)),
                 flippedQuery=as.character(flippedQuery(x)))
    if (x_nmc > 0L) {
        tmp <- do.call(data.frame, c(lapply(x_mcols, showAsCell),
                                     list(check.names=FALSE)))
        ans <- cbind(ans, `|`=rep.int("|", x_len), as.matrix(tmp))
    }
    ans
}

showOverlapEncodings <- function(x, margin="", print.classinfo=FALSE)
{
    x_class <- class(x)
    x_len <- length(x)
    x_mcols <- mcols(x, use.names=FALSE)
    x_nmc <- if (is.null(x_mcols)) 0L else ncol(x_mcols)
    cat(classNameForDisplay(x), " object of length ", x_len,
        " with ",
        x_nmc, " metadata column", ifelse(x_nmc == 1L, "", "s"),
        ":\n", sep="")
    out <- makePrettyMatrixForCompactPrinting(x,
                           .make_naked_matrix_from_OverlapEncodings)
    if (print.classinfo) {
        .COL2CLASS <- c(
            Loffset="integer",
            Roffset="integer",
            encoding="factor",
            flippedQuery="logical"
        )
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
}

setMethod("show", "OverlapEncodings",
    function(object)
        showOverlapEncodings(object, margin="  ", print.classinfo=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Detection of "splice compatible" overlaps.
###

build_compatible_encoding_subpatterns <- function(njunc)
{
    ## Each "atom" must match exactly 1 code in the encoding.
    ATOM0 <- "[fgij]"
    if (njunc == 0L)
        return(ATOM0)
    #ntimes <- function(atom, n) rep.int(atom, n)
    ntimes <- function(atom, n) {
        if (n == 1L) atom else c(atom, "{", n, "}")
    }
    LEFT_ATOM <- "[jg]"
    MIDDLE_ATOM <- "g"
    RIGHT_ATOM <- "[gf]"
    WILDCARD_ATOM <- "[^:-]"
    sapply(seq_len(njunc + 1L),
           function(i) {
               if (i == 1L) {
                   atoms <- c(LEFT_ATOM, ntimes(WILDCARD_ATOM, njunc))
               } else if (i == njunc + 1L) {
                   atoms <- c(ntimes(WILDCARD_ATOM, njunc), RIGHT_ATOM)
               } else {
                   atoms <- c(ntimes(WILDCARD_ATOM, i-1L),
                              MIDDLE_ATOM,
                              ntimes(WILDCARD_ATOM, njunc-i+1L))
               }
               paste0(atoms, collapse="")
           })
}

.build_compatible_encoding_pattern <- function(max.njunc)
{
    subpatterns <- sapply(0:max.njunc,
        function(njunc)
            paste0(build_compatible_encoding_subpatterns(njunc), collapse=":")
    )
    paste0(":(", paste0(subpatterns, collapse="|"), "):")
}

### 'x' must be a character vector of single-end encodings.
.is_compatible_with_splicing <- function(x)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    max.njunc <- max(c(0L, .njunc_single_end_encodings(x)))
    pattern <- .build_compatible_encoding_pattern(max.njunc)
    setNames(grepl(pattern, x), names(x))
}

setGeneric("isCompatibleWithSplicing",
    function(x) standardGeneric("isCompatibleWithSplicing")
)

setMethod("isCompatibleWithSplicing", "character",
    function(x)
    {
        halves <- encodingHalves(x, single.end.on.left=TRUE,
                                    single.end.on.right=TRUE)
        ans1 <- .is_compatible_with_splicing(halves[[1L]])
        ans2 <- .is_compatible_with_splicing(halves[[2L]])
        ans1 & ans2
    }
)

setMethod("isCompatibleWithSplicing", "factor",
    function(x) isCompatibleWithSplicing(levels(x))[as.integer(x)]
)

setMethod("isCompatibleWithSplicing", "OverlapEncodings",
    function(x) isCompatibleWithSplicing(encoding(x))
)

