### =========================================================================
### pileLettersAt()
### -------------------------------------------------------------------------


### .pileLettersOnSingleRefAt() is the workhorse behind pileLettersAt().
### 'x', 'pos', 'cigar': 3 parallel vectors describing N strings aligned
### to the same reference sequence. 'x' must be an XStringSet (typically
### DNAStringSet) object containing the unaligned strings (a.k.a. the
### query sequences) reported with respect to the + strand. 'pos' must
### be an integer vector where 'pos[i]' is the 1-based position on the
### reference sequence of the first aligned letter in 'x[[i]]'. 'cigar'
### must be a character vector containing the extended CIGAR strings.
### 'at': must be an integer vector containing the individual positions of
### interest with respect to the reference sequence.
### Returns an XStringSet (typically DNAStringSet) object parallel to
### 'at' (i.e. with 1 string per position of interest).
.pileLettersOnSingleRefAt <- function(x, pos, cigar, at)
{
    stopifnot(is(x, "XStringSet"))
    N <- length(x)  # nb of alignments
    stopifnot(is.integer(pos) && length(pos) == N)
    stopifnot(is.character(cigar) && length(cigar) == N)
    stopifnot(is.integer(at))

    ops <- c("M", "=", "X")
    ranges_on_ref <- cigarRangesAlongReferenceSpace(cigar, pos=pos, ops=ops)
    ranges_on_query <- cigarRangesAlongQuerySpace(cigar, ops=ops)

    ## 'ranges_on_ref' and 'ranges_on_query' are IRangesList objects parallel
    ## to 'x', 'pos', and 'cigar'. In addition, the 2 IRangesList objects
    ## have the same "shape" (i.e. same elementLengths()), so, after
    ## unlisting, the 2 unlisted objects are parallel IRanges objects.
    unlisted_ranges_on_ref <- unlist(ranges_on_ref, use.names=FALSE)
    unlisted_ranges_on_query <- unlist(ranges_on_query, use.names=FALSE)

    ## 2 integer vectors parallel to IRanges objects 'unlisted_ranges_on_ref'
    ## and 'unlisted_ranges_on_query' above.
    range_group <- togroup(ranges_on_ref)
    query2ref_shift <- start(unlisted_ranges_on_ref) -
                       start(unlisted_ranges_on_query)

    hits <- findOverlaps(at, unlisted_ranges_on_ref)
    hits_at_in_x <- at[queryHits(hits)] - query2ref_shift[subjectHits(hits)]
    hits_group <- range_group[subjectHits(hits)]
    unlisted_piles <- subseq(x[hits_group], start=hits_at_in_x, width=1L)
    piles_skeleton <- PartitioningByEnd(queryHits(hits), NG=length(at),
                                       names=names(at))
    piles <- relist(unlisted_piles, piles_skeleton)
    unstrsplit(piles)
}

### 'x', 'seqnames', 'pos', 'cigar': 4 parallel vectors describing N
### aligned strings. 'x', 'pos', and 'cigar' as above. 'seqnames' must
### be a factor-Rle where 'seqnames[i]' is the name of the reference
### sequence of the i-th alignment.
### 'at': must be a GRanges object containing the individual genomic
### positions of interest. 'seqlevels(at)' must be identical to
### 'levels(seqnames)'.
### Returns an XStringSet (typically DNAStringSet) object parallel to
### 'at' (i.e. with 1 string per genomic position of interest).
pileLettersAt <- function(x, seqnames, pos, cigar, at)
{
    stopifnot(is(x, "XStringSet"))
    N <- length(x)  # nb of alignments
    stopifnot(is(seqnames, "Rle"))
    stopifnot(is.factor(runValue(seqnames)))
    stopifnot(length(seqnames) == N)
    stopifnot(is.integer(pos) && length(pos) == N)
    stopifnot(is.character(cigar) && length(cigar) == N)
    stopifnot(is(at, "GRanges"))
    stopifnot(all(width(at) == 1L))
    stopifnot(identical(seqlevels(at), levels(seqnames)))

    ## We process 1 chromosome at a time. So we start by splitting
    ## 'x', 'pos', 'cigar', and 'start(at)' by chromosome. The 4
    ## resulting list-like objects have 1 list element per chromosome
    ## in 'seqlevels(at)' (or in 'levels(seqnames)', which is identical
    ## to 'seqlevels(at)').
    x_by_chrom <- split(x, seqnames)               # XStringSetList
    pos_by_chrom <- split(pos, seqnames)           # IntegerList
    cigar_by_chrom <- split(cigar, seqnames)       # CharacterList
    at_by_chrom <- split(start(at), seqnames(at))  # IntegerList

    ## Unsplit index.
    split_idx <- unlist(split(seq_along(at), seqnames(at)),
                        use.names=FALSE)
    unsplit_idx <- integer(length(at))
    unsplit_idx[split_idx] <- seq_along(at)

    do.call("c", lapply(seq_along(seqlevels(at)),
        function(i)
            .pileLettersOnSingleRefAt(
                x_by_chrom[[i]],
                pos_by_chrom[[i]],
                cigar_by_chrom[[i]],
                at_by_chrom[[i]])))[unsplit_idx]
}

