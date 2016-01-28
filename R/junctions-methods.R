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
        first_junctions <- junctions(first(x, real.strand=TRUE))
        last_junctions <- junctions(last(x, real.strand=TRUE))
        ## pc() is a fast "parallel c()" for list-like objects.
        ## In the case below, it's equivalent to (but faster than) doing
        ## 'mendoapply(c, first_junctions, last_junctions)'.
        ans <- pc(first_junctions, last_junctions)
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
### Natural intron motifs taken from:
###   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC84117/

NATURAL_INTRON_MOTIFS <- c("GT-AG", "GC-AG", "AT-AC", "AT-AA", "AT-AG")


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

.infer_intron_strand <- function(unoriented_intron_motif)
{
    natural_intron_motifs <- DNAStringSet(NATURAL_INTRON_MOTIFS)
    intron_strand <- rep.int(NA, length(unoriented_intron_motif))
    idx <- which(unoriented_intron_motif %in%
                 natural_intron_motifs)
    intron_strand[idx] <- FALSE
    idx <- which(unoriented_intron_motif %in%
                 reverseComplement(natural_intron_motifs))
    intron_strand[idx] <- TRUE
    if (any(is.na(intron_strand)))
        warning("For some junctions, the dinucleotides found at the intron ",
                "boundaries don't\n  match any of the natural intron motifs ",
                "stored in predefined character vector\n  ",
                "'NATURAL_INTRON_MOTIFS'. For these junctions, the ",
                "intron_motif and\n  intron_strand metadata columns ",
                "were set to NA and *, respectively.")
    strand(intron_strand)
}

.orient_intron_motif <- function(unoriented_intron_motif, intron_strand)
{
    ans <- unoriented_intron_motif
    idx <- which(intron_strand == "-")
    ans[idx] <- reverseComplement(ans[idx])
    ans <- factor(as.character(ans), levels=NATURAL_INTRON_MOTIFS)
}

summarizeJunctions <- function(x, with.revmap=FALSE, genome=NULL)
{
    if (!isTRUEorFALSE(with.revmap))
        stop("'with.revmap' must be TRUE or FALSE")
    if (!is.null(genome)) {
        if (!suppressWarnings(require(BSgenome, quietly=TRUE)))
            stop("you need to install the BSgenome package in order ",
                 "to use the 'genome' argument")
        genome <- BSgenome::getBSgenome(genome)
    }

    x_junctions <- junctions(x)
    unlisted_junctions <- unlist(x_junctions, use.names=FALSE)
    unstranded_unlisted_junctions <- unstrand(unlisted_junctions)
    ans <- sort(unique(unstranded_unlisted_junctions))
    unq2dups <- as(findMatches(ans, unstranded_unlisted_junctions), "List")
    ans_score <- elementNROWS(unq2dups)
    tmp <- extractList(strand(unlisted_junctions), unq2dups)
    ans_plus_score <- sum(tmp == "+")
    ans_minus_score <- sum(tmp == "-")
    ans_mcols <- DataFrame(score=ans_score,
                           plus_score=ans_plus_score,
                           minus_score=ans_minus_score)
    if (with.revmap) {
        crossed_by <- togroup(x_junctions)
        ans_revmap <- extractList(crossed_by, unq2dups)
        ## 'ans_revmap' should never contain duplicates when 'x' is a
        ## GAlignments object, because a given junction can show up at most
        ## once per SAM/BAM record (i.e. per element in 'x', or per alignment).
        ## This doesn't hold anymore if the elements in 'x' consist of more
        ## than 1 SAM/BAM record (or alignment) e.g. if 'x' is a
        ## GAlignmentPairs or GAlignmentsList object, because, in that case,
        ## the same junction can show up more than once per element in 'x'.
        if (!is(x, "GAlignments"))
            ans_revmap <- unique(ans_revmap)
        ans_mcols$revmap <- ans_revmap
    }
    if (!is.null(genome)) {
        unoriented_intron_motif <- .extract_unoriented_intron_motif(genome,
                                                                    ans)
        ans_intron_strand <- .infer_intron_strand(unoriented_intron_motif)
        ans_intron_motif <- .orient_intron_motif(unoriented_intron_motif,
                                                 ans_intron_strand)
        ans_mcols <- cbind(ans_mcols,
                           DataFrame(intron_motif=ans_intron_motif,
                                     intron_strand=ans_intron_strand))
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
###                                     file.is.raw.juncs=TRUE)
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
    stopifnot(all(elementNROWS(blocks) == 2L))
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
readTopHatJunctions <- function(file, file.is.raw.juncs=FALSE)
{
    if (!isTRUEorFALSE(file.is.raw.juncs))
        stop("'file.is.raw.juncs' must be TRUE or FALSE")
    if (is.character(file)) {
        if (!isSingleString(file))
            stop("'file' must be a single string")
        file_ext0 <- ".bed"
        file_ext <- substr(file, start=nchar(file) - nchar(file_ext0) + 1L,
                                 stop=nchar(file))
        if (file.is.raw.juncs) {
            if (file_ext == file_ext0)
                stop("'file.is.raw.juncs=TRUE' is not aimed to be ",
                     "used on a file\n  with the .bed extension")
            df <- read.table(file, stringsAsFactors=FALSE)
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
                    "'file.is.raw.juncs=TRUE'")
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

.STAR_INTRON_MOTIFS <- c("GT-AG", "CT-AC",
                         "GC-AG", "CT-GC",
                         "AT-AC", "GT-AT")

.get_STAR_intron_motif_levels <- function()
{
    ans <- .STAR_INTRON_MOTIFS[(1:3)*2L - 1L]
    stopifnot(all(ans == NATURAL_INTRON_MOTIFS[1:3]))
    rev_motifs <- .STAR_INTRON_MOTIFS[(1:3)*2L]
    stopifnot(all(rev_motifs == reverseComplement(DNAStringSet(ans))))
    ans
}

readSTARJunctions <- function(file)
{
    motif123 <- .get_STAR_intron_motif_levels()
    df <- read.table(file, stringsAsFactors=FALSE)
    ans_seqnames <- df[[1L]]
    ans_start <- df[[2L]]
    ans_end <- df[[3L]]
    ans_strand <- strand(df[[4L]] == 2L)
    STAR_intron_motif_code <- df[[5L]]
    if (!is.integer(ans_start) || !is.integer(ans_end)
     || !is.integer(STAR_intron_motif_code)
     || S4Vectors:::anyMissingOrOutside(STAR_intron_motif_code,
                                      lower=0L, upper=6L))
        stop("'file' does not look like a junction file generated ",
             "by the STAR aligner (normally the SJ.out.tab file)")
    STAR_intron_motif_code[STAR_intron_motif_code == 0L] <- NA_integer_
    code1 <- STAR_intron_motif_code + 1L
    ans_intron_motif <- factor(motif123[code1 %/% 2L], levels=motif123)
    ans_intron_strand <- strand(as.logical(code1 %% 2L))
    has_code_zero <- is.na(ans_intron_motif)
    stopifnot(identical(has_code_zero, ans_intron_strand == "*"))
    idx0 <- which(!has_code_zero)
    if (!identical(ans_strand[idx0], ans_intron_strand[idx0]))
        warning("For some junctions, the strand reported in the motif_strand ",
                "metadata column\n  (which was inferred from the STAR intron ",
                "motif code stored in column 5 of\n  'file') is conflicting ",
                "with the strand of the junction reported in column 4\n  ",
                "of 'file'. Bug in STAR? Obscure feature? Or corrupted file? ",
                "Please ask on the\n  STAR general user mailing list ",
                "(https://groups.google.com/d/forum/rna-star)\n  for ",
                "clarifications about this (only if you're confident that ",
                "your SJ.out.tab\n  file is not corrupted though).")
    GRanges(ans_seqnames,
            IRanges(ans_start, ans_end),
            strand=ans_strand,
            intron_motif=ans_intron_motif,
            intron_strand=ans_intron_strand,
            um_reads=df[[7L]],
            mm_reads=df[[8L]])
}

