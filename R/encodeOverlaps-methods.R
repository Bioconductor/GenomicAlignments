### =========================================================================
### encodeOverlaps() and related utilities
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### encodeOverlaps1() - A low-level utility.
###
###   > query <- IRanges(start=c(7, 15, 22), end=c(9, 19, 23))
###   > subject <- IRanges(start=c(1, 4, 15, 22, 1, 30, 25),
###                        end=c(2, 9, 19, 25, 10, 38, 25))
###   > encodeOverlaps1(query, subject, as.matrix=TRUE)
###       [,1] [,2] [,3] [,4] [,5] [,6] [,7]
###   [1,] "m"  "j"  "a"  "a"  "i"  "a"  "a" 
###   [2,] "m"  "m"  "g"  "a"  "m"  "a"  "a" 
###   [3,] "m"  "m"  "m"  "f"  "m"  "a"  "a" 
###   > encodeOverlaps1(query, subject)
###   $Loffset
###   [1] 1
###   
###   $Roffset
###   [1] 2
###   
###   $encoding
###   [1] "3:jmm:agm:aaf:imm:"
###
###   > query.space <- c(0, 1, 0)
###   > encodeOverlaps1(query, subject, query.space=query.space)$encoding
###   [1] "3:mXm:jXm:aXm:aXf:iXm:aXa:aXa:"
###   > query.space <- rep(-1, length(query))
###   > subject.space <- rep(-1, length(subject))
###   > encodeOverlaps1(rev(query), rev(subject),
###                     query.space=query.space, subject.space=subject.space)
###   $Loffset
###   [1] 2
###
###   $Roffset
###   [1] 1
###
###   $encoding
###   [1] "3:aai:jmm:agm:aaf:"
###
###   > encodeOverlaps1(query, subject, query.break=2)$encoding
###   [1] "2--1:jm--m:ag--m:aa--f:im--m:"
###   > encodeOverlaps1(rev(query), rev(subject),
###                     query.space=query.space, subject.space=subject.space,
###                     query.break=1)$encoding
###   [1] "1--2:a--ai:j--mm:a--gm:a--af:"

### 'query.space' must be either an integer vector of the same length as
### 'query', or NULL. If NULL, then it's interpreted as
### 'integer(length(query))' i.e. all the ranges in 'query' are considered to
### be on space 0.
encodeOverlaps1 <- function(query, subject,
                            query.space=NULL, subject.space=NULL,
                            query.break=0L, flip.query=FALSE,
                            as.matrix=FALSE, as.raw=FALSE)
{
    if (!is(query, "IntegerRanges"))
        stop("'query' must be an IntegerRanges object")
    if (!is(subject, "IntegerRanges"))
        stop("'subject' must be an IntegerRanges object")
    if (is.numeric(query.space) && !is.integer(query.space))
        query.space <- as.integer(query.space)
    if (is.numeric(subject.space) && !is.integer(subject.space))
        subject.space <- as.integer(subject.space)
    if (!isSingleNumber(query.break))
        stop("'query.break' must be a single integer value")
    if (!is.integer(query.break))
        query.break <- as.integer(query.break)
    if (!isTRUEorFALSE(flip.query))
        stop("'flip.query' must be TRUE or FALSE")
    if (!isTRUEorFALSE(as.matrix))
        stop("'as.matrix' must be TRUE or FALSE")
    if (!isTRUEorFALSE(as.raw))
        stop("'as.raw' must be TRUE or FALSE")
    .Call2("encode_overlaps1",
           start(query), width(query), query.space,
           query.break, flip.query,
           start(subject), width(subject), subject.space,
           as.matrix, as.raw,
           PACKAGE="GenomicAlignments")
}

### TODO: Put this in the (upcoming) man page for encodeOverlaps().
### A simple (but inefficient) implementation of the "findOverlaps" method for
### IntegerRanges objects. Complexity and memory usage is M x N where M and N
### are the lengths of 'query' and 'subject', respectively.
findRangesOverlaps <- function(query, subject)
{
    ovenc <- encodeOverlaps1(query, subject, as.matrix=TRUE, as.raw=TRUE)
    offsets <- which(charToRaw("c") <= ovenc & ovenc <= charToRaw("k")) - 1L
    q_hits <- offsets %% nrow(ovenc) + 1L
    s_hits <- offsets %/% nrow(ovenc) + 1L
    cbind(queryHits=q_hits, subjectHits=s_hits)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The .RangesList_encodeOverlaps() helper.
###
### This is the power horse behind all the "encodeOverlaps" methods.
###

.RangesList_encode_overlaps <- function(query.starts, query.widths,
                                        query.spaces, query.breaks,
                                        subject.starts, subject.widths,
                                        subject.spaces)
{
    .Call2("RangesList_encode_overlaps",
           query.starts, query.widths, query.spaces, query.breaks,
           subject.starts, subject.widths, subject.spaces,
           PACKAGE="GenomicAlignments")
}

.Hits_encode_overlaps <- function(query.starts, query.widths,
                                  query.spaces, query.breaks,
                                  subject.starts, subject.widths,
                                  subject.spaces,
                                  hits, flip.query)
{
    if (queryLength(hits) != length(query.starts) ||
        subjectLength(hits) != length(subject.starts))
        stop("'hits' is not compatible with 'query' and 'subject'")
    .Call2("Hits_encode_overlaps",
           query.starts, query.widths, query.spaces, query.breaks,
           subject.starts, subject.widths, subject.spaces,
           queryHits(hits), subjectHits(hits), flip.query,
           PACKAGE="GenomicAlignments")
}

.RangesList_encodeOverlaps <- function(query.starts, query.widths,
                                       subject.starts, subject.widths,
                                       hits, flip.query=NULL,
                                       query.spaces=NULL, subject.spaces=NULL,
                                       query.breaks=NULL)
{
    if (is.null(hits)) {
        C_ans <- .RangesList_encode_overlaps(query.starts, query.widths,
                                             query.spaces, query.breaks,
                                             subject.starts, subject.widths,
                                             subject.spaces)
        flip.query <- logical(length(C_ans$encoding))
    } else {
        if (!is(hits, "Hits"))
            stop("'hits' must be a Hits object")
        if (is.null(flip.query)) {
            flip.query <- logical(length(hits))
        } else {
            if (!is.logical(flip.query))
                stop("'flip.query' must be a logical vector")
            if (length(flip.query) != length(hits))
                stop("'flip.query' must have the same length as 'hits'")
        }
        C_ans <- .Hits_encode_overlaps(query.starts, query.widths,
                                       query.spaces, query.breaks,
                                       subject.starts, subject.widths,
                                       subject.spaces,
                                       hits, flip.query)
    }
    encoding <- as.factor(C_ans$encoding)
    new2("OverlapEncodings", Loffset=C_ans$Loffset, Roffset=C_ans$Roffset,
                             encoding=encoding, flippedQuery=flip.query,
                             check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### encodeOverlaps() generic and methods for IntegerRangesList objects.
###

setGeneric("encodeOverlaps", signature=c("query", "subject"),
    function(query, subject, hits=NULL, ...) standardGeneric("encodeOverlaps")
)

setMethods("encodeOverlaps", list(c("IntegerRangesList", "IntegerRangesList"),
                                  c("IntegerRangesList", "IntegerRanges"),
                                  c("IntegerRanges", "IntegerRangesList")),
    function(query, subject, hits=NULL, ...)
    {
        .RangesList_encodeOverlaps(as.list(start(query)),
                                   as.list(width(query)),
                                   as.list(start(subject)),
                                   as.list(width(subject)),
                                   hits)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .isWrongStrand() internal helper.
###

.oneValPerTopLevelElt <- function(x, errmsg)
{
    if (!is(x, "RleList"))
        stop("'x' must be an RleList object")
    vals <- runValue(x)
    vals_eltNROWS <- elementNROWS(vals)
    if (!all(vals_eltNROWS == 1L))
        stop(errmsg)
    unlist(vals, use.names=FALSE)
}

.isWrongStrand <- function(query, subject, hits)
{
    if (!is(query, "GRangesList") || !is(subject, "GRangesList"))
        stop("'query' and 'subject' must be GRangesList objects")

    ## Extract the top-level strand and seqnames of the query.
    errmsg <- c("some alignments in 'query' have ranges on ",
                "more than 1 reference sequence (fusion reads?)")
    query_seqnames <- .oneValPerTopLevelElt(seqnames(query), errmsg)
    errmsg <- c("some alignments in 'query' have ranges on ",
                "both strands")
    query_strand <- .oneValPerTopLevelElt(strand(query), errmsg)

    ## Extract the top-level strand and seqnames of the subject.
    errmsg <- c("some transcripts in 'subject' mix exons from ",
                "different chromosomes (trans-splicing?)")
    subject_seqnames <- .oneValPerTopLevelElt(seqnames(subject), errmsg)
    errmsg <- c("some transcripts in 'subject' mix exons from ",
                "both strands (trans-splicing?)")
    subject_strand <- .oneValPerTopLevelElt(strand(subject), errmsg)

    ## Expand the top-level strand and seqnames of the query and subject.
    if (!is.null(hits)) {
        if (!is(hits, "Hits"))
            stop("'hits' must be NULL or a Hits object")
        if (queryLength(hits) != length(query) ||
            subjectLength(hits) != length(subject))
            stop("'hits' is not compatible with 'query' and 'subject' ",
                 "('queryLength(hits)' and 'subjectLength(hits)' don't ",
                 "match the lengths of 'query' and 'subject')")
        query_seqnames <- query_seqnames[queryHits(hits)]
        query_strand <- query_strand[queryHits(hits)]
        subject_seqnames <- subject_seqnames[subjectHits(hits)]
        subject_strand <- subject_strand[subjectHits(hits)]
    }

    ## Should never happen if 'encodeOverlaps(query, subject, hits)'
    ## was called with 'hits' being the result of a call to
    ## 'findOverlaps(query, subject)'.
    if (!all(query_seqnames == subject_seqnames))
        stop("cannot use 'flip.query.if.wrong.strand=TRUE' to ",
             "encode overlaps across chromosomes")

    query_strand != subject_strand
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### flipQuery()
###

flipQuery <- function(x, i)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    i <- normalizeSingleBracketSubscript(i, x, as.NSBS=TRUE)
    xi <- extractROWS(x, i)
    x <- replaceROWS(x, i, invertStrand(revElements(xi)))
    xi_query.break <- mcols(xi)$query.break
    if (!is.null(xi_query.break)) {
        revxi_query.break <- elementNROWS(xi) - xi_query.break
        mcols(x)$query.break <- replaceROWS(mcols(x)$query.break, i,
                                            revxi_query.break)
    }
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Should we use generic + methods for this?
###

.get_GRanges_spaces <- function(x)
{
        ans <- as.integer(seqnames(x))
        x_strand <- as.integer(strand(x))
        is_minus <- which(x_strand == as.integer(strand("-")))
        ans[is_minus] <- - ans[is_minus]
        ans
}

.get_GRangesList_spaces <- function(x)
{
        unlisted_ans <- .get_GRanges_spaces(x@unlistData)
        as.list(relist(unlisted_ans, x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "encodeOverlaps" method for GRangesList objects.
###

.GRangesList_encodeOverlaps <- function(query, subject, hits,
                                        flip.query.if.wrong.strand)
{
    if (!isTRUEorFALSE(flip.query.if.wrong.strand))
        stop("'flip.query.if.wrong.strand' must be TRUE or FALSE")
    seqinfo <- merge(seqinfo(query), seqinfo(subject))
    seqlevels(query) <- seqlevels(subject) <- seqlevels(seqinfo)
    if (flip.query.if.wrong.strand) {
        flip.query <- .isWrongStrand(query, subject, hits)
    } else {
        flip.query <- NULL
    }
    query.breaks <- mcols(query)$query.break
    .RangesList_encodeOverlaps(as.list(start(query)),
                               as.list(width(query)),
                               as.list(start(subject)),
                               as.list(width(subject)),
                               hits, flip.query,
                               query.spaces=.get_GRangesList_spaces(query),
                               subject.spaces=.get_GRangesList_spaces(subject),
                               query.breaks=query.breaks)
}

setMethod("encodeOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, hits=NULL, flip.query.if.wrong.strand=FALSE)
        .GRangesList_encodeOverlaps(query, subject, hits,
                                    flip.query.if.wrong.strand)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### selectEncodingWithCompatibleStrand().
###

selectEncodingWithCompatibleStrand <- function(ovencA, ovencB,
                                               query.strand, subject.strand,
                                               hits=NULL)
{
    if (!is(ovencA, "OverlapEncodings"))
        stop("'ovencA' must be an OverlapEncodings object")
    if (!is(ovencB, "OverlapEncodings"))
        stop("'ovencB' must be an OverlapEncodings object")
    if (!is.null(hits)) {
        if (!is(hits, "Hits"))
            stop("'hits' must be a Hits object or NULL")
        query.strand <- query.strand[queryHits(hits)]
        subject.strand <- subject.strand[subjectHits(hits)]
    }
    ans <- ovencA
    names(ans) <- NULL
    mcols(ans) <- NULL
    is_wrong_strand <- query.strand != subject.strand
    idx <- which(is_wrong_strand)
    ans@Loffset[idx] <- ovencB@Loffset[idx]
    ans@Roffset[idx] <- ovencB@Roffset[idx]
    ans_encoding <- as.character(ans@encoding)
    ans_encoding[idx] <- as.character(ovencB@encoding[idx])
    ans@encoding <- as.factor(ans_encoding)
    ans@flippedQuery[is_wrong_strand] <- TRUE
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCompatibleWithSkippedExons().
###

### FIXME: Revisit this and make sure it's doing the right thing for
### paired-end encodings. Maybe look at isCompatibleWithSplicing() for some
### inspiration (used to suffer from similar issue on paired-end encodings but
### was refactored).

setGeneric("isCompatibleWithSkippedExons", signature="x",
    function(x, max.skipped.exons=NA)
        standardGeneric("isCompatibleWithSkippedExons")
)

.build_CompatibleWithSkippedExons_pattern0 <- function(max.njunc1,
                                                       max.Lnjunc, max.Rnjunc,
                                                       max.skipped.exons=NA)
{
    if (!identical(max.skipped.exons, NA))
        stop("only 'max.skipped.exons=NA' is supported for now, sorry")

    ## Subpattern for single-end reads.
    skipped_exons_subpatterns <- c(":(.:)*", ":(..:)*",
                                   ":(...:)*", ":(....:)*")
    subpattern1 <- sapply(0:max.njunc1,
                     function(njunc)
                       paste0(build_compatible_encoding_subpatterns(njunc),
                              collapse=skipped_exons_subpatterns[njunc+1L]))
    subpattern1 <- paste0(":(", paste0(subpattern1, collapse="|"), "):")

    ## Subpattern for paired-end reads.
    Lsubpattern <- sapply(0:max.Lnjunc,
                     function(njunc)
                       paste0(":",
                              build_compatible_encoding_subpatterns(njunc),
                              "-",
                              collapse=".*"))
    Lsubpattern <- paste0("(", paste0(Lsubpattern, collapse="|"), ")")

    Rsubpattern <- sapply(0:max.Rnjunc,
                     function(njunc)
                       paste0("-",
                              build_compatible_encoding_subpatterns(njunc),
                              ":",
                              collapse=".*"))
    Rsubpattern <- paste0("(", paste0(Rsubpattern, collapse="|"), ")")

    LRsubpattern <- paste0(Lsubpattern, ".*", Rsubpattern)

    ## Final pattern.
    paste0("(", subpattern1, "|", LRsubpattern, ")")
}

.build_CompatibleWithSkippedExons_pattern <- function(x, max.skipped.exons=NA)
{
    njunc <- njunc(x)
    Lnjunc <- Lnjunc(x)
    Rnjunc <- Rnjunc(x)
    max.njunc1 <- max(c(0L, njunc[is.na(Lnjunc)]))
    max.Lnjunc <- max(c(0L, Lnjunc), na.rm=TRUE)
    max.Rnjunc <- max(c(0L, Rnjunc), na.rm=TRUE)
    .build_CompatibleWithSkippedExons_pattern0(max.njunc1,
                    max.Lnjunc, max.Rnjunc,
                    max.skipped.exons=max.skipped.exons)
}

.isCompatibleWithSkippedExons <- function(x, max.skipped.exons=NA)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    pattern1 <- .build_CompatibleWithSkippedExons_pattern(x,
                                                          max.skipped.exons)
    grepl(pattern1, x) & !isCompatibleWithSplicing(x)
}

setMethod("isCompatibleWithSkippedExons", "character",
    .isCompatibleWithSkippedExons
)

setMethod("isCompatibleWithSkippedExons", "factor",
    function(x, max.skipped.exons=NA)
    {
        ok <- isCompatibleWithSkippedExons(levels(x),
                                           max.skipped.exons=max.skipped.exons)
        ok[as.integer(x)]
    }
)

setMethod("isCompatibleWithSkippedExons", "OverlapEncodings",
    function(x, max.skipped.exons=NA)
        isCompatibleWithSkippedExons(encoding(x),
                        max.skipped.exons=max.skipped.exons)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractSteppedExonRanks().
###

.extract_njunc_from_encoding <- function(x)
{
    as.integer(unlist(strsplit(sub(":.*", "", x), "--", fixed=TRUE),
                      use.names=FALSE)) - 1L
}

.extractSteppedExonRanksFromEncodingBlocks <- function(encoding_blocks,
                                                       encoding_patterns)
{
    patterns <- paste0("^", encoding_patterns, "$")
    ii <- lapply(patterns, grep, encoding_blocks)
    ii_eltNROWS <- elementNROWS(ii)
    if (any(ii_eltNROWS == 0L))
        return(integer(0))
    if (any(ii_eltNROWS != 1L))
        stop("cannot unambiguously extract stepped exon ranks from ",
             "encoding \"", paste0(encoding_blocks, collapse=":"), "\"")
    ans <- unlist(ii, use.names=FALSE)
    diff_ans <- diff(ans)
    if (any(diff_ans <= 0L))
        return(integer(0))
    ans
}

### 'encoding' must be a single encoding.
### Returns an integer vector. If the encoding is single-end then the vector
### is unnamed and strictly sorted. If it's paired-end then it's made of the
### concatenation of the 2 strictly sorted integer vectors that correspond to
### each end. In that case the vector is named.
.extractSteppedExonRanks <- function(encoding, for.query.right.end=FALSE)
{
    if (!isTRUEorFALSE(for.query.right.end))
        stop("'for.query.right.end' must be TRUE or FALSE")
    encoding_blocks <- strsplit(encoding, ":", fixed=TRUE)[[1L]]
    njunc <- .extract_njunc_from_encoding(encoding_blocks[1L])
    encoding_blocks <- encoding_blocks[-1L]
    if (length(njunc) == 1L) {
        ## Single-end read.
        if (for.query.right.end)
            stop("cannot use 'for.query.right.end=TRUE' ",
                 "on single-end encoding: ", encoding)
        encoding_patterns <- build_compatible_encoding_subpatterns(njunc)
        return(.extractSteppedExonRanksFromEncodingBlocks(encoding_blocks,
                                                          encoding_patterns))
    }
    if (length(njunc) != 2L)  # should never happen
        stop(encoding, ": invalid encoding")
    ## Paired-end read.
    encoding_blocks <- strsplit(encoding_blocks, "--", fixed=TRUE)
    if (!all(elementNROWS(encoding_blocks) == 2L))  # should never happen
        stop(encoding, ": invalid encoding")
    encoding_blocks <- matrix(unlist(encoding_blocks, use.names=FALSE), nrow=2L)
    Lencoding_patterns <- build_compatible_encoding_subpatterns(njunc[1L])
    Lranks <- .extractSteppedExonRanksFromEncodingBlocks(encoding_blocks[1L, ],
                                                         Lencoding_patterns)
    Rencoding_patterns <- build_compatible_encoding_subpatterns(njunc[2L])
    Rranks <- .extractSteppedExonRanksFromEncodingBlocks(encoding_blocks[2L, ],
                                                         Rencoding_patterns)
    if (for.query.right.end)
        return(Rranks)  # unnamed! (like for a single-end read)
    names(Rranks) <- rep.int("R", length(Rranks))
    names(Lranks) <- rep.int("L", length(Lranks))
    c(Lranks, Rranks)
}

setGeneric("extractSteppedExonRanks",
    function(x, for.query.right.end=FALSE)
        standardGeneric("extractSteppedExonRanks")
)

setMethod("extractSteppedExonRanks", "character",
    function(x, for.query.right.end=FALSE)
    {
        lapply(x, .extractSteppedExonRanks, for.query.right.end)
    }
)

setMethod("extractSteppedExonRanks", "factor",
    function(x, for.query.right.end=FALSE)
    {
        if (length(x) == 0L)
            return(list())
        ranks <- extractSteppedExonRanks(levels(x),
                     for.query.right.end=for.query.right.end)
        ranks[as.integer(x)]
    }
)

setMethod("extractSteppedExonRanks", "OverlapEncodings",
    function(x, for.query.right.end=FALSE)
    {
        ranks <- extractSteppedExonRanks(encoding(x),
                     for.query.right.end=for.query.right.end)
        ranks_eltNROWS <- elementNROWS(ranks)
        tmp <- unlist(unname(ranks), use.names=TRUE)  # we want the inner names
        tmp <- tmp + rep.int(Loffset(x), ranks_eltNROWS)
        flevels <- seq_len(length(ranks))
        f <- factor(rep.int(flevels, ranks_eltNROWS), levels=flevels)
        unname(split(tmp, f))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractSpannedExonRanks().
###

setGeneric("extractSpannedExonRanks",
    function(x, for.query.right.end=FALSE)
        standardGeneric("extractSpannedExonRanks")
)

setMethod("extractSpannedExonRanks", "character",
    function(x, for.query.right.end=FALSE)
    {
        .extractRanks <- function(encoding) {
            ranks <- .extractSteppedExonRanks(encoding,
                          for.query.right.end=for.query.right.end)
            if (length(ranks) == 0L)
                return(c(NA_integer_, NA_integer_))
            c(ranks[1L], ranks[length(ranks)])
        }
        ranks <- lapply(x, .extractRanks)
        if (length(ranks) == 0L) {
            firstSpannedExonRank <- lastSpannedExonRank <- integer(0)
        } else {
            ranks <- unlist(ranks, use.names=FALSE)
            firstSpannedExonRank <- ranks[c(TRUE, FALSE)]
            lastSpannedExonRank <- ranks[c(FALSE, TRUE)]
        }
        data.frame(firstSpannedExonRank=firstSpannedExonRank,
                   lastSpannedExonRank=lastSpannedExonRank,
                   check.names=FALSE, stringsAsFactors=FALSE)
    }
)

setMethod("extractSpannedExonRanks", "factor",
    function(x, for.query.right.end=FALSE)
    {
        if (length(x) == 0L)
            return(list())
        ranks <- extractSpannedExonRanks(levels(x),
                     for.query.right.end=for.query.right.end)
        ans <- ranks[as.integer(x), , drop=FALSE]
        rownames(ans) <- NULL
        ans
    }
)

setMethod("extractSpannedExonRanks", "OverlapEncodings",
    function(x, for.query.right.end=FALSE)
    {
        ranks <- extractSpannedExonRanks(encoding(x),
                     for.query.right.end=for.query.right.end)
        ranks$firstSpannedExonRank <- ranks$firstSpannedExonRank + Loffset(x)
        ranks$lastSpannedExonRank <- ranks$lastSpannedExonRank + Loffset(x)
        ranks
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractSkippedExonRanks().
###

setGeneric("extractSkippedExonRanks",
    function(x, for.query.right.end=FALSE)
        standardGeneric("extractSkippedExonRanks")
)

setMethod("extractSkippedExonRanks", "character",
    function(x, for.query.right.end=FALSE)
    {
        .extractRanks <- function(encoding) {
            ranks <- .extractSteppedExonRanks(encoding,
                          for.query.right.end=for.query.right.end)
            if (length(ranks) == 0L)
                return(ranks)
            ranks_names <- names(ranks)
            if (is.null(ranks_names))  # single-end read
                return(setdiff(ranks[1L]:ranks[length(ranks)], ranks))
            ## Paired-end read.
            ranks <- split(unname(ranks), ranks_names)
            Lranks <- ranks$L
            Lranks <- setdiff(Lranks[1L]:Lranks[length(Lranks)], Lranks)
            Rranks <- ranks$R
            Rranks <- setdiff(Rranks[1L]:Rranks[length(Rranks)], Rranks)
            names(Lranks) <- rep.int("L", length(Lranks))
            names(Rranks) <- rep.int("R", length(Rranks))
            c(Lranks, Rranks)
        }
        lapply(x, .extractRanks)
    }
)

setMethod("extractSkippedExonRanks", "factor",
    function(x, for.query.right.end=FALSE)
    {
        if (length(x) == 0L)
            return(list())
        ranks <- extractSkippedExonRanks(levels(x),
                     for.query.right.end=for.query.right.end)
        ranks[as.integer(x)]
    }
)

setMethod("extractSkippedExonRanks", "OverlapEncodings",
    function(x, for.query.right.end=FALSE)
    {
        ranks <- extractSkippedExonRanks(encoding(x),
                     for.query.right.end=for.query.right.end)
        ranks_eltNROWS <- elementNROWS(ranks)
        tmp <- unlist(unname(ranks), use.names=TRUE)  # we want the inner names
        tmp <- tmp + rep.int(Loffset(x), ranks_eltNROWS)
        flevels <- seq_len(length(ranks))
        f <- factor(rep.int(flevels, ranks_eltNROWS), levels=flevels)
        unname(split(tmp, f))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractQueryStartInTranscript().
###

### TODO: Maybe put this in IRanges and rename it setElementNROWS, or even
### better introduce an "elementNROWS<-" generic and make this the method
### for CompressedList objects?
.setElementNROWS <- function(x, eltNROWS)
{
    if (!is(x, "CompressedList"))
        stop("'x' must be a CompressedList object")
    if (!is.numeric(eltNROWS) || length(eltNROWS) != length(x))
        stop("'eltNROWS' must be an integer vector of the same length as 'x'")
    if (!is.integer(eltNROWS))
        eltNROWS <- as.integer(eltNROWS)
    if (S4Vectors:::anyMissingOrOutside(eltNROWS, lower=0L))
        stop("'eltNROWS' cannot contain NAs or negative values")
    x_eltNROWS <- elementNROWS(x)
    if (!all(eltNROWS <= x_eltNROWS))
        stop("'all(eltNROWS <= elementNROWS(x))' must be TRUE")
    offset <- cumsum(c(0L, x_eltNROWS[-length(x_eltNROWS)]))
    ii <- S4Vectors:::fancy_mseq(eltNROWS, offset=offset)
    x@unlistData <- x@unlistData[ii]
    x@partitioning@end <- unname(cumsum(eltNROWS))
    x
}

### Returns a data.frame with 1 row per overlap, and 3 integer columns:
###     1. startInTranscript
###     2. firstSpannedExonRank
###     3. startInFirstSpannedExon
### Rows for overlaps that are not "compatible" or "almost compatible"
### contain NAs.
extractQueryStartInTranscript <- function(query, subject,
                                          hits=NULL, ovenc=NULL,
                                          flip.query.if.wrong.strand=FALSE,
                                          for.query.right.end=FALSE)
{
    if (!is(query, "GRangesList") || !is(subject, "GRangesList"))
        stop("'query' and 'subject' must be GRangesList objects")
    seqinfo <- merge(seqinfo(query), seqinfo(subject))
    seqlevels(query) <- seqlevels(subject) <- seqlevels(seqinfo)
    if (is.null(hits)) {
        if (length(query) != length(subject))
            stop("'query' and 'subject' must have the same length")
    } else {
        if (!is(hits, "Hits"))
            stop("'hits' must be a Hits object or NULL")
        if (queryLength(hits) != length(query) ||
            subjectLength(hits) != length(subject))
            stop("'hits' is not compatible with 'query' and 'subject' ",
                 "('queryLength(hits)' and 'subjectLength(hits)' don't ",
                 "match the lengths of 'query' and 'subject')")
        query <- query[queryHits(hits)]
        subject <- subject[subjectHits(hits)]
    }
    if (is.null(ovenc)) {
        ovenc <- encodeOverlaps(query, subject,
                         flip.query.if.wrong.strand=flip.query.if.wrong.strand)
    } else {
        if (!is(ovenc, "OverlapEncodings"))
            stop("'ovenc' must be an OverlapEncodings object")
        if (length(ovenc) != length(query))
            stop("when not NULL, 'ovenc' must have the same length ",
                 "as 'hits', if specified, otherwise as 'query'")
    }
    if (!isTRUEorFALSE(for.query.right.end))
        stop("'for.query.right.end' must be TRUE or FALSE")

    query <- flipQuery(query, flippedQuery(ovenc))

    ## Extract first range from each list element of 'query'.
    if (for.query.right.end) {
        query.break <- mcols(query)$query.break
        if (is.null(query.break))
            stop("using 'for.query.right.end=TRUE' requires that ",
                 "'mcols(query)' has a \"query.break\" column ",
                 "indicating for each paired-end read the position of the ",
                 "break between the ranges coming from one end and those ",
                 "coming from the other end")
        query <- tails(query, n=-query.break)
    }
    query1 <- unlist(heads(query, n=1L), use.names=FALSE)
    ## A sanity check.
    if (length(query1) != length(query)) 
        stop("some list elements in 'query' are empty")

    query_start1 <- start(query1)
    query_end1 <- end(query1)
    query_strand1 <- as.factor(strand(query1))

    ## Extract start/end/strand of the first spanned exon
    ## in each top-level element of 'subject'.
    exrank <- extractSpannedExonRanks(ovenc,
                  for.query.right.end=for.query.right.end)$firstSpannedExonRank
    sii1 <- start(subject@partitioning) + exrank - 1L
    subject_start1 <- start(subject@unlistData)[sii1]
    subject_end1 <- end(subject@unlistData)[sii1]
    subject_strand1 <- as.factor(strand(subject@unlistData))[sii1]

    ## A sanity check.
    if (any(!is.na(exrank) & (query_strand1 != subject_strand1))) {
        ## TODO: Error message needs to take into account whether 'hits'
        ## and/or 'ovenc' was supplied or not.
        stop("'ovenc' is incompatible with the supplied 'query' ",
             "and/or 'subject' and/or 'hits'")
    }

    ## Compute the "query start in first spanned exon".
    startInFirstSpannedExon <- rep.int(NA_integer_, length(query))
    is_on_plus <- query_strand1 == "+"
    idx <- which(!is.na(exrank) & is_on_plus)
    startInFirstSpannedExon[idx] <- query_start1[idx] - subject_start1[idx] + 1L
    idx <- which(!is.na(exrank) & !is_on_plus)
    startInFirstSpannedExon[idx] <- subject_end1[idx] - query_end1[idx] + 1L

    ## Truncate each transcript in 'subject' right before the first spanned
    ## exon and compute the cumulated width of the truncated object.
    subject2_eltNROWS <- exrank - 1L
    subject2_eltNROWS[is.na(exrank)] <- 0L
    subject2 <- .setElementNROWS(subject, subject2_eltNROWS)
    subject2_cumwidth <- unname(sum(width(subject2)))
    subject2_cumwidth[is.na(exrank)] <- NA_integer_

    ## Compute the "query start in transcript".
    startInTranscript <- subject2_cumwidth + startInFirstSpannedExon

    data.frame(startInTranscript=startInTranscript,
               firstSpannedExonRank=exrank,
               startInFirstSpannedExon=startInFirstSpannedExon,
               check.names=FALSE, stringsAsFactors=FALSE)
}

