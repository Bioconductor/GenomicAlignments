### =========================================================================
### Extract junctions from genomic alignments
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### junctions() generic and methods.
###

setGeneric("junctions", signature="x",
    function(x, use.mcols=FALSE, ...) standardGeneric("junctions")
)

setMethod("junctions", "GAlignments",
    function(x, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        grl <- grglist(x, order.as.in.query=TRUE)
        ans <- psetdiff(granges(x), grl)
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)

setMethod("junctions", "GAlignmentPairs",
    function(x, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        first_junctions <- junctions(x@first)
        last_junctions <- junctions(invertRleStrand(x@last))
        ## Fast way of doing mendoapply(c, first_junctions, last_junctions)
        ## on 2 CompressedList objects.
        ans <- c(first_junctions, last_junctions)
        collate_subscript <-
            IRanges:::make_XYZxyz_to_XxYyZz_subscript(length(x))
        ans <- ans[collate_subscript]
        ans <- shrinkByHalf(ans)
        names(ans) <- names(x)
        if (use.mcols) {
            mcols(ans) <- mcols(x)
        } else {
            mcols(ans) <- NULL
        }
        ans
    }
)

setMethod("junctions", "GAlignmentsList",
    function(x, use.mcols=FALSE, ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand)
            strand(x@unlistData) <- "*"
        grl <- junctions(x@unlistData)
        ans_breakpoints <- end(grl@partitioning)[end(x@partitioning)]
        ans_partitioning <- PartitioningByEnd(ans_breakpoints, names=names(x))
        ans <- relist(grl@unlistData, ans_partitioning)
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### summarizeJunctions()
###

.extract_unoriented_intron_motif <- function(genome, junctions)
{
    mcols(junctions) <- NULL
    junctions_len <- length(junctions)
    Ldinucl_gr <- Rdinucl_gr <- junctions
    end(Ldinucl_gr) <- start(Ldinucl_gr) + 1L
    start(Rdinucl_gr) <- end(Rdinucl_gr) - 1L
    all_dinucl <- getSeq(genome, c(Ldinucl_gr, Rdinucl_gr))
    Ldinucl <- head(all_dinucl, n=junctions_len)
    Rdinucl <- tail(all_dinucl, n=junctions_len)
    xscat(Ldinucl, "-", Rdinucl)
}

### Natural intron motifs taken from:
###   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC84117/
.NATURAL_INTRON_MOTIFS <- c("GT-AG", "GC-AG", "AT-AC", "AT-AA", "AT-AG")

.infer_intron_strand <- function(unoriented_intron_motif)
{
    natural_intron_motifs <- DNAStringSet(.NATURAL_INTRON_MOTIFS)
    intron_strand <- rep.int(NA, length(unoriented_intron_motif))
    idx <- which(unoriented_intron_motif %in%
                 natural_intron_motifs)
    intron_strand[idx] <- FALSE
    idx <- which(unoriented_intron_motif %in%
                 reverseComplement(natural_intron_motifs))
    intron_strand[idx] <- TRUE
    if (any(is.na(intron_strand)))
        warning("strand of some introns could not be determined")
    strand(intron_strand)
}

.orient_intron_motif <- function(unoriented_intron_motif, intron_strand)
{
    ans <- unoriented_intron_motif
    idx <- which(intron_strand == "-")
    ans[idx] <- reverseComplement(ans[idx])
    ans <- factor(as.character(ans), levels=.NATURAL_INTRON_MOTIFS)
}

summarizeJunctions <- function(x, with.revmap=FALSE, genome=NULL)
{
    if (!isTRUEorFALSE(with.revmap))
        stop("'with.revmap' must be TRUE or FALSE")
    if (!is.null(genome))
        genome <- getBSgenome(genome)

    x_junctions <- junctions(x)
    unlisted_junctions0 <- unlist(x_junctions, use.names=FALSE)
    unlisted_junctions <- unstrand(unlisted_junctions0)
    ans <- sort(unique(unlisted_junctions))
    unq2dups <- as(findMatches(ans, unlisted_junctions), "List")
    ans_score <- elementLengths(unq2dups)
    tmp <- extractList(strand(unlisted_junctions0), unq2dups)
    ans_plus_score <- sum(tmp == "+")
    ans_minus_score <- sum(tmp == "-")
    ans_mcols <- DataFrame(score=ans_score,
                           plus_score=ans_plus_score,
                           minus_score=ans_minus_score)
    if (with.revmap) {
        supported_by <- togroup(x_junctions)
        ans_mcols$revmap <- extractList(supported_by, unq2dups)
    }
    if (!is.null(genome)) {
        unoriented_intron_motif <- .extract_unoriented_intron_motif(genome,
                                                                    ans)
        ans_intron_strand <- .infer_intron_strand(unoriented_intron_motif)
        ans_intron_motif <- .orient_intron_motif(unoriented_intron_motif,
                                                 ans_intron_strand)
        ans_mcols <- cbind(ans_mcols,
                           DataFrame(intron_strand=ans_intron_strand,
                                     intron_motif=ans_intron_motif))
    }
    mcols(ans) <- ans_mcols
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readTopHatJunctions()
###
### Read splice junctions file (junctions.bed) generated by the TopHat
### aligner into a GRanges object.
### Usage:
###   readTopHatJunctions("junctions.bed")
###
### Comparing with output of bed_to_juncs script (assuming the
### 'new_list.juncs' file was obtained by passing 'junctions.bed' thru
### bed_to_juncs):
###   junctions1 <- readTopHatJunctions("junctions.bed")
###   junctions2 <- readTopHatJunctions("new_list.juncs",
###                                     file.is.bed_to_juncs.output=TRUE)
###   stopifnot(all(junctions1 == junctions2))
###

.bed_to_Juncs <- function(x)
{
    if (!is(x, "GRanges"))
        stop("'x' must be a GRanges object")
    if (!identical(mcols(x)$thick, ranges(x)))
        stop("this BED file doesn't look like the junctions.bed file ",
             "generated by TopHat")
    blocks <- mcols(x)$blocks
    stopifnot(all(elementLengths(blocks) == 2L))
    unlisted_blocks <- unlist(blocks, use.names=FALSE)
    even_idx <- 2L * seq_along(x)
    odd_idx <- even_idx - 1L
    ans_start <- start(x) + end(unlisted_blocks)[odd_idx]
    ans_end <- start(x) + start(unlisted_blocks)[even_idx] - 2L
    ans <- GRanges(seqnames(x), IRanges(ans_start, ans_end), strand=strand(x))
    mcols(ans) <- DataFrame(name=mcols(x)$name,
                            score=as.integer(mcols(x)$score))
    ans
}

### 'file' must be the path or a connection object to a junctions.bed file as
### generated by TopHat, or to a tab-delimited file obtained by running
### TopHat's bed_to_juncs script on a junctions.bed file.
### Returns the junctions in a GRanges object.
### IMPORTANT NOTE: readTopHatJunctions() does NOT follow the convention used
### by TopHat that describes a junction by the position of the nucleotide
### immediately before and after the intron. In the GRanges object returned
### by readTopHatJunctions(), a junction is considered to start at the
### left-most and to end at the right-most nucleotide of the intron.
readTopHatJunctions <- function(file, file.is.bed_to_juncs.output=FALSE)
{
    if (!isTRUEorFALSE(file.is.bed_to_juncs.output))
        stop("'file.is.bed_to_juncs.output' must be TRUE or FALSE")
    if (is.character(file)) {
        if (!isSingleString(file))
            stop("'file' must be a single string")
        file_ext0 <- ".bed"
        file_ext <- substr(file, start=nchar(file) - nchar(file_ext0) + 1L,
                                 stop=nchar(file))
        if (file.is.bed_to_juncs.output) {
            if (file_ext == file_ext0)
                stop("'file.is.bed_to_juncs.output=TRUE' is not aimed to be ",
                     "used on a file\n  with the .bed extension")
            df <- read.table(file)
            ## The 2nd and 3rd columns in 'new_list.juncs' are the left and
            ## right positions of the junctions, respectively. The convention
            ## used by TopHat is that these are NOT the positions of the
            ## left-most and right-most nucleotides of the intron, but rather
            ## the positions immediately before and after, respectively, that
            ## is, the last and the first positions of the flanking exons.
            ## Also these positions are *both* 0-based.
            ans_ranges <- IRanges(df[[2L]] + 2L, df[[3L]])
            ans <- GRanges(df[[1L]], ans_ranges, strand=df[[4L]])
            return(ans)
        }
        if (file_ext != file_ext0)
            warning("'file' has no .bed extension, suggesting it may not ",
                    "be a junctions.bed\n  file as generated by TopHat. ",
                    "I will assume it is this file anyway (or a BED\n  file ",
                    "with similar content). If 'file' is a tab-delimited ",
                    "file obtained\n  by running TopHat's bed_to_juncs script ",
                    "on a junctions.bed file, you\n  should use ",
                    "'file.is.bed_to_juncs.output=TRUE'")
    }
    junctions_bed <- rtracklayer::import(file)
    .bed_to_Juncs(junctions_bed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readSTARJunctions()
###
### Read splice junctions file (SJ.out.tab) generated by the STAR aligner
### into a GRanges object.
### Usage:
###   readSTARJunctions("SJ.out.tab")
###

readSTARJunctions <- function(file)
{
    df <- read.table(file)
    GRanges(df[[1L]],
            IRanges(df[[2L]], df[[3L]]),
            strand=strand(df[[4L]] == 2L),
            um_reads=df[[7L]],
            mm_reads=df[[8L]])
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Old stuff (deprecated & defunct)
###

introns <- function(...)
{
    .Deprecated("junctions")
    junctions(...)
}

