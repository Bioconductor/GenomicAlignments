### =========================================================================
### stackStringsFromGAlignments() & related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### stackStringsFromGAlignments()
###

### All the alignments in GAlignments object 'x' must be on the **same**
### chromosome. This is NOT checked!
.stack_reads <- function(x, from, to, what="seq",
                         D.letter="-", N.letter=".",
                         Lpadding.letter="+", Rpadding.letter="+")
{
    x_mcols <- mcols(x, use.names=FALSE)
    what_col_idx <- match(what, colnames(x_mcols))
    if (is.na(what_col_idx))
        stop(wmsg("'x' does not have a \"", what, "\" metadata column"))
    what_col <- x_mcols[[what_col_idx]]
    if (what == "qual")
        what_col <- BStringSet(what_col)
    layed_seq <- sequenceLayer(what_col, cigar(x),
                               D.letter=D.letter, N.letter=N.letter)
    ans <- stackStrings(layed_seq, from, to,
                        shift=start(x)-1L,
                        Lpadding.letter=Lpadding.letter,
                        Rpadding.letter=Rpadding.letter)
    names(ans) <- names(x)
    mcols(ans) <- x_mcols
    ans
}

stackStringsFromGAlignments <- function(x, region, what="seq",
                                        D.letter="-", N.letter=".",
                                        Lpadding.letter="+",
                                        Rpadding.letter="+")
{
    if (!is(x, "GAlignments"))
        stop(wmsg("'x' must be a GAlignments object"))
    if (!is(region, "GRanges"))
        region <- as(region, "GRanges")
    if (length(region) != 1L)
        stop(wmsg("'region' must contain a single genomic range"))
    region_seqname <- seqlevelsInUse(region)
    if (!(region_seqname %in% seqlevels(x)))
        stop(wmsg("seqlevel not in 'x': ", region_seqname))
    what <- match.arg(what, c("seq", "qual"))
    x <- x[overlapsAny(granges(x), region)]
    .stack_reads(x, start(region), end(region), what=what,
                 D.letter=D.letter, N.letter=N.letter,
                 Lpadding.letter=Lpadding.letter,
                 Rpadding.letter=Rpadding.letter)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### stackStringsFromBam()
###

### Should always return a ScanBamParam object containing exactly 1 genomic
### region.
.normarg_param <- function(param)
{
    if (isSingleString(param)) {
        tmp1 <- strsplit(param, ":", fixed=TRUE)[[1L]]
        if (length(tmp1) != 2L)
            stop(wmsg("when a character string, 'param' must be ",
                      "of the form \"chr14:5201-5300\""))
        tmp2 <- as.integer(strsplit(tmp1[2L], "-", fixed=TRUE)[[1L]])
        if (length(tmp2) != 2L || any(is.na(tmp2)))
            stop(wmsg("when a character string, 'param' must be ",
                      "of the form \"chr14:5201-5300\""))
        param <- GRanges(tmp1[1L], IRanges(tmp2[1L], tmp2[2L]))
    }
    if (is(param, "GenomicRanges")) {
        if (length(param) != 1L)
            stop(wmsg("when a GRanges object, 'param' must have length 1"))
        seqlevels(param) <- seqlevelsInUse(param)
        param <- ScanBamParam(which=param)
        return(param)
    }
    if (is(param, "IntegerRangesList")) {
        ## We support IntegerRangesList just because ScanBamParam() supports
        ## it too and also because that's what's returned by bamWhich().
        param <- param[elementNROWS(param) != 0L]
        if (length(unlist(param, use.names=FALSE)) != 1L)
            stop(wmsg("when an IntegerRangesList object, 'param' must contain ",
                      "exactly 1 genomic region (i.e. 'unlist(param)' must ",
                      "have length 1)"))
        param <- ScanBamParam(which=param)
        return(param)
    }
    if (!is(param, "ScanBamParam"))
        stop(wmsg("'param' must be either a ScanBamParam or IntegerRangesList ",
                  "object containing exactly 1 genomic region, or a GRanges ",
                  "object of length 1, or a character string specifying a ",
                  "singe genomic region (in the \"chr14:5201-5300\" format)"))
    param_which <- bamWhich(param)
    param_which <- param_which[elementNROWS(param_which) != 0L]
    if (length(unlist(param_which, use.names=FALSE)) != 1L)
        stop(wmsg("when a ScanBamParam object, 'param' must contain exactly ",
                  "1 genomic region (i.e. 'unlist(bamWhich(param))' must ",
                  "have length 1)"))
    bamWhich(param) <- param_which
    param
}

stackStringsFromBam <- function(file, index=file, param,
                                what="seq", use.names=FALSE,
                                D.letter="-", N.letter=".",
                                Lpadding.letter="+", Rpadding.letter="+")
{
    param <- .normarg_param(param)
    region_range <- unlist(bamWhich(param), use.names=FALSE)
    what <- match.arg(what, c("seq", "qual"))
    param_what <- bamWhat(param)
    if (!(what %in% param_what))
        bamWhat(param) <- c(param_what, what)
    gal <- readGAlignments(file, index=index,
                           use.names=use.names, param=param)
    ans <- .stack_reads(gal, start(region_range), end(region_range), what=what,
                        D.letter=D.letter, N.letter=N.letter,
                        Lpadding.letter=Lpadding.letter,
                        Rpadding.letter=Rpadding.letter)
    if (!(what %in% param_what)) {
        ## Remove the what metadata column.
        ans_mcols <- mcols(ans, use.names=FALSE)
        what_col_idx <- match(what, colnames(ans_mcols))
        ans_mcols <- ans_mcols[ , -what_col_idx, drop=FALSE]
        ## Sadly, subsetting a DataFrame will mangle the colnames of the
        ## returned DataFrame if it has duplicated colnames. Since we of
        ## course don't want this, we fix them.
        colnames(ans_mcols) <- param_what
        mcols(ans) <- ans_mcols
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### alphabetFrequencyFromBam()
###

alphabetFrequencyFromBam <- function(file, index=file, param,
                                     what="seq", ...)
{
    param <- .normarg_param(param)
    region_range <- unlist(bamWhich(param), use.names=FALSE)
    region_seqname <- names(bamWhich(param))
    what <- match.arg(what, c("seq", "qual"))
    bamWhat(param) <- what
    gal <- readGAlignments(file, index=index, param=param)
    seqlevels(gal) <- region_seqname
    what_col <- mcols(gal, use.names=FALSE)[ , what]
    if (what == "qual")
        what_col <- BStringSet(what_col)
    at <- start(region_range) - 1L + seq_len(width(region_range))
    at <- GRanges(region_seqname, IRanges(at, width=1L))
    piles <- pileLettersAt(what_col, seqnames(gal), start(gal), cigar(gal),
                           at)
    alphabetFrequency(piles, ...)
}

