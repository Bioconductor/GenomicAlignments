### =========================================================================
### findMateAlignment()
### -------------------------------------------------------------------------
###
### For each element in GAlignments object 'x', finds its mate in GAlignments
### object 'y'.
###
### Alignments 'x[i1]' and 'y[i2]' are considered mates iff they pass all the
### following tests:
###
###   (A) names(x[i1]) == names(y[i2])
###
###   (B) mcols(x[i1])$mrnm == seqnames(y[i2]) &
###       mcols(y[i2])$mrnm == seqnames(x[i1])
###
###   (C) mcols(x[i1])$mpos == start(y[i2]) &
###       mcols(y[i2])$mpos == start(x[i1])
###
###   (D) isMateMinusStrand(x[i1]) == isMinusStrand(y[i2]) &
###       isMateMinusStrand(y[i2]) == isMinusStrand(x[i1])
###
###   (E) isFirstSegment(x[i1]) & isLastSegment(y[i2]) |
###       isFirstSegment(y[i2]) & isLastSegment(x[i1])
###
###   (F) isProperPair(x[i1]) == isProperPair(y[i2])
###
###   (G) isSecondaryAlignment(x[i1]) == isSecondaryAlignment(y[i2])

.checkMetadatacols <- function(arg, argname)
{
    if (!is(arg, "GAlignments"))
        stop("'", argname, "' must be a GAlignments object")
    if (is.null(names(arg)))
        stop("'", argname, "' must have names")
    arg_mcols <- mcols(arg)
    REQUIRED_COLNAMES <- c("flag", "mrnm", "mpos")
    if (!all(REQUIRED_COLNAMES %in% colnames(arg_mcols))) {
        colnames_in1string <-
            paste0("\"", REQUIRED_COLNAMES, "\"", collapse=", ")
        stop("required columns in 'mcols(", argname, ")': ",
             colnames_in1string)
    }
    if (!is.integer(arg_mcols$flag))
        stop("'mcols(", argname, ")$flag' must be an integer vector")
    if (!is.factor(arg_mcols$mrnm))
        stop("'mcols(", argname, ")$mrnm' must be a factor")
    if (!identical(levels(arg_mcols$mrnm), levels(seqnames(arg))))
        stop("'mcols(", argname, ")$mrnm' and 'seqnames(", argname, ")' ",
             "must have exactly the same levels in the same order")
    if (!is.integer(arg_mcols$mpos))
        stop("'mcols(", argname, ")$mpos' must be an integer vector")
    arg_mcols
}

### 'names', 'flagbits', 'mrnm', and 'mpos', must all come from the same
###     GAlignments object x.
### 'names': names(x).
### 'flagbits': integer matrix (of 0's and 1's) obtained with
###     bamFlagAsBitMatrix(mcols(x)$flag, bitnames=.MATING_FLAG_BITNAMES)
### 'mrnm': factor obtained with mcols(x)$mrnm
### 'mpos': integer vector obtained with mcols(x)$mpos
### Returns 'names' with NAs injected at positions corresponding to alignments
### that satisfy at least one of following conditions:
###     1. Bit 0x1 (isPaired) is 0
###     2. Read is neither first or last mate
###     3. Bit 0x8 (hasUnmappedMate) is 1
###     4. 'mrnm' is NA (i.e. RNEXT = '*')
###     5. 'mpos' is NA (i.e. PNEXT = 0)
### My understanding of the SAM Spec is that 3., 4. and 5. should happen
### simultaneously even though the Spec don't clearly state this.

.MATING_FLAG_BITNAMES <- c("isPaired", "hasUnmappedMate",
                           "isFirstMateRead", "isSecondMateRead")

.makeGAlignmentsGNames <- function(names, flagbits, mrnm, mpos)
{
    is_paired <- flagbits[ , "isPaired"]
    is_first <- flagbits[ , "isFirstMateRead"]
    is_last <- flagbits[ , "isSecondMateRead"]
    has_unmappedmate <- flagbits[ , "hasUnmappedMate"]
    alter_idx <- which(!is_paired |
                       is_first == is_last |
                       has_unmappedmate |
                       is.na(mrnm) |
                       is.na(mpos))
    names[alter_idx] <- NA_integer_
    names
}

### Puts NAs last.
.getCharacterOrderAndGroupSizes <- function(x)
{
    x2 <- match(x, x, nomatch=.Machine$integer.max,
                      incomparables=NA_character_)
    xo <- base::order(x2)
    ox2 <- Rle(x2[xo])
    group.sizes <- runLength(ox2)
    ngroup <- length(group.sizes)
    if (ngroup != 0L && runValue(ox2)[ngroup] == .Machine$integer.max)
        group.sizes <- group.sizes[-ngroup]
    list(xo=xo, group.sizes=group.sizes)
}

### Should return the same as:
###   args <- as.list(setNames(rep(TRUE, length(bitnames)), bitnames))
###   tmp <- do.call(scanBamFlag, args)
###   tmp[[2L]] - tmp[[1L]]
.makeFlagBitmask <- function(bitnames)
{
    bitpos <- match(bitnames, FLAG_BITNAMES)
    sum(as.integer(2L ^ (bitpos-1L)))
}

### 'x_hits' and 'y_hits' must be 2 integer vectors of the same length N
### representing the N edges of a bipartite graph between the [1, x_len] and
### [1, y_len] intervals (the i-th edge being represented by (x[i], y[i])).
### Returns an integer vector F of length 'x_len' where F[k] is defined by:
###   - If there is no occurence of k in 'x', then F[k] = NA.
###   - If there is more than 1 occurence of k in 'x', then F[k] = 0.
###   - If there is exactly 1 occurence of k in 'x', at index i_k, then
###     F[k] = y[i_k].
### In addition, if more than 1 value of index k is associated to F[k], then
### F[k] is replaced by -F[k].
.makeMateIdx2 <- function(x_hits, y_hits, x_len)
{
    idx1 <- which(has_duplicates(y_hits))
    y_hits[idx1] <- - y_hits[idx1]
    idx2 <- which(has_duplicates(x_hits))
    y_hits[idx2] <- 0L
    ans <- rep.int(NA_integer_, x_len)
    ans[x_hits] <- y_hits
    ans
}

.showGAlignmentsEltsWithMoreThan1Mate <- function(x, idx)
{
    if (length(idx) == 0L)
        return()
    cat("\n!! Found more than 1 mate for the following elements in 'x': ",
        paste(idx, collapse=", "),
        ".\n!! Details:\n!! ", sep="")
    showGAlignments(x[idx], margin="!! ",
                            print.classinfo=TRUE,
                            print.seqinfo=FALSE)
    cat("!! ==> won't assign a mate to them!\n")
}

.dump_envir <- new.env(hash=TRUE, parent=emptyenv())
.dumpEnvir <- function() .dump_envir

flushDumpedAlignments <- function()
{
    objnames <- ls(envir=.dumpEnvir())
    rm(list=objnames, envir=.dumpEnvir())
}

dumpAlignments <- function(gal)
{
    objnames <- ls(envir=.dumpEnvir())
    nobj <- length(objnames)
    if (nobj == 0L) {
        new_objname <- 1L
    } else {
        new_objname <- as.integer(objnames[nobj]) + 1L
    }
    new_objname <- sprintf("%08d", new_objname)
    assign(new_objname, gal, envir=.dumpEnvir())
}

countDumpedAlignments <- function()
{
    sum(unlist(eapply(.dumpEnvir(), length, USE.NAMES=FALSE)))
}

getDumpedAlignments <- function()
{
    objnames <- ls(envir=.dumpEnvir())
    args <- unname(mget(objnames, envir=.dumpEnvir()))
    do.call(c, args)
}

### Takes about 2.3 s and 170MB of RAM to mate 1 million alignments,
### and about 13 s and 909MB of RAM to mate 5 million alignments.
findMateAlignment <- function(x)
{
    x_names <- names(x)
    if (is.null(x_names))
        stop("'x' must have names")
    x_mcols <- .checkMetadatacols(x, "x")
    ## flushDumpedAlignments() must be placed *after* the first reference to
    ## 'x', otherwise, when doing 'findMateAlignment(getDumpedAlignments())',
    ## the flushing would happen before 'x' is evaluated, causing 'x' to be
    ## evaluated to NULL.
    flushDumpedAlignments()
    x_flag <- x_mcols$flag
    bitnames <- c(.MATING_FLAG_BITNAMES, "isMinusStrand", "isMateMinusStrand")
    x_flagbits <- bamFlagAsBitMatrix(x_flag, bitnames=bitnames)
    x_mrnm <- x_mcols$mrnm
    x_mpos <- x_mcols$mpos
    x_gnames <- .makeGAlignmentsGNames(x_names, x_flagbits, x_mrnm, x_mpos)
    x_seqnames <- as.factor(seqnames(x))
    x_start <- start(x)

    xo_and_GS <- .getCharacterOrderAndGroupSizes(x_gnames)
    xo <- xo_and_GS$xo
    group.sizes <- xo_and_GS$group.sizes
    ans <- Rsamtools:::.findMateWithinGroups(group.sizes,
                           x_flag[xo], x_seqnames[xo],
                           x_start[xo], x_mrnm[xo], x_mpos[xo])
    dumpme_idx <- which(ans <= 0L)
    if (length(dumpme_idx) != 0L) {
        dumpAlignments(x[xo[dumpme_idx]])
        ans[dumpme_idx] <- NA_integer_
    }
    ans[xo] <- xo[ans]  # isn't that cute!
    dump_count <- countDumpedAlignments()
    if (dump_count != 0L)
        warning("  ", dump_count, " alignments with ambiguous pairing ",
                "were dumped.\n    Use 'getDumpedAlignments()' to retrieve ",
                "them from the dump environment.")
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeGAlignmentPairs().
###

### TODO: Make isFirstSegment() an S4 generic function with methods for
### matrices, integer vectors, and GAlignments objects. Put this with the
### flag utils in Rsamtools.
.isFirstSegment.matrix <- function(x)
{
    is_paired <- as.logical(x[ , "isPaired"])
    is_first0 <- as.logical(x[ , "isFirstMateRead"])
    is_last0 <- as.logical(x[ , "isSecondMateRead"])
    ## According to SAM Spec, bits 0x40 (isFirstMateRead) and 0x80
    ## (isSecondMateRead) can both be set or unset, even when bit 0x1
    ## (isPaired) is set. However we are not interested in those situations
    ## (which have a special meaning).
    is_paired & is_first0 & (!is_last0)
}

.isFirstSegment.integer <- function(flag)
{
    bitnames <- c("isPaired", "isFirstMateRead", "isSecondMateRead")
    .isFirstSegment.matrix(bamFlagAsBitMatrix(flag, bitnames=bitnames))
}

.isFirstSegment.GAlignments <- function(x)
    .isFirstSegment.integer(mcols(x, use.names=FALSE)$flag)

### TODO: Make isLastSegment() an S4 generic function with methods for
### matrices, integer vectors, and GAlignments objects. Put this with the
### flag utils in Rsamtools.
.isLastSegment.matrix <- function(x)
{
    is_paired <- as.logical(x[ , "isPaired"])
    is_first0 <- as.logical(x[ , "isFirstMateRead"])
    is_last0 <- as.logical(x[ , "isSecondMateRead"])
    ## According to SAM Spec, bits 0x40 (isFirstMateRead) and 0x80
    ## (isSecondMateRead) can both be set or unset, even when bit 0x1
    ## (isPaired) is set. However we are not interested in those situations
    ## (which have a special meaning).
    is_paired & is_last0 & (!is_first0)
}

.isLastSegment.integer <- function(flag)
{
    bitnames <- c("isPaired", "isFirstMateRead", "isSecondMateRead")
    .isLastSegment.matrix(bamFlagAsBitMatrix(flag, bitnames=bitnames))
}

.isLastSegment.GAlignments <- function(x)
    .isLastSegment.integer(mcols(x, use.names=FALSE)$flag)

### 'x' must be a GAlignments objects.
makeGAlignmentPairs <- function(x, use.names=FALSE, use.mcols=FALSE,
                                   strandMode=1)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!isTRUEorFALSE(use.mcols)) {
        if (!is.character(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE or a character vector ",
                 "specifying the metadata columns to propagate")
        if (!all(use.mcols %in% colnames(mcols(x, use.names=FALSE))))
            stop("'use.mcols' must be a subset of 'colnames(mcols(x))'")
    }
    mate <- findMateAlignment(x)
    x_is_first <- .isFirstSegment.GAlignments(x)
    x_is_last <- .isLastSegment.GAlignments(x)
    first_idx <- which(!is.na(mate) & x_is_first)
    last_idx <- mate[first_idx]

    ## Fundamental property of the 'mate' vector: it's a permutation of order
    ## 2 and with no fixed point on the set of indices for which 'mate' is
    ## not NA.
    ## Check there are no fixed points.
    if (!all(first_idx != last_idx))
        stop("findMateAlignment() returned an invalid 'mate' vector")
    ## Check order 2 (i.e. permuting a 2nd time brings back the original
    ## set of indices).
    if (!identical(mate[last_idx], first_idx))
        stop("findMateAlignment() returned an invalid 'mate' vector")
    ## One more sanity check.
    if (!all(x_is_last[last_idx]))
        stop("findMateAlignment() returned an invalid 'mate' vector")

    ## Check the 0x2 bit (isProperPair).
    x_is_proper <- as.logical(bamFlagAsBitMatrix(mcols(x, use.names=FALSE)$flag,
                                                 bitnames="isProperPair"))
    ans_is_proper <- x_is_proper[first_idx]

    ## The big split!
    ans_first <- x[first_idx]
    ans_last <- x[last_idx]
    ans_names <- NULL
    if (use.names)
        ans_names <- names(ans_first)
    names(ans_first) <- names(ans_last) <- NULL
    if (is.character(use.mcols)) {
        mcols(ans_first) <- mcols(ans_first, use.names=FALSE)[use.mcols]
        mcols(ans_last) <- mcols(ans_last, use.names=FALSE)[use.mcols]
    } else if (!use.mcols) {
        mcols(ans_first) <- mcols(ans_last) <- NULL
    }
    GAlignmentPairs(ans_first, ans_last, strandMode=strandMode,
                    isProperPair=ans_is_proper, names=ans_names)
}

