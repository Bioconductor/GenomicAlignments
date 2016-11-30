### =========================================================================
### "findSpliceOverlaps" methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helper functions
###

.rangeForSorted <- function(x)
{
    part <- PartitioningByWidth(x)
    xflat <- unlist(x, use.names=FALSE)
    IRanges(start(xflat)[start(part)], end(xflat)[end(part)])
}

.compatibleTranscription <- function(query, subject, splice, clip=0L)
{
    qrng <- ranges(query)
    srng <- ranges(subject)
    sprng <- ranges(splice)
    ## FIXME : should clip be a modifiable parameter?
    ## Or we could consider taking it out entirely... the aligner should clip
    if (clip != 0L) {
        bounds <- .rangeForSorted(qrng) - clip
        qrange <- restrict(qrng, start(bounds), end(bounds),
                           keep.all.ranges = TRUE)
    }
    bnds <- elementNROWS(setdiff(qrng, srng)) == 0L
    splc <- elementNROWS(intersect(srng, sprng)) == 0L
    bnds & splc
}

.novelBounds <- function(query, subject, qhits)
{
    qrange <- range(query)
    qstrand <- as.character(strand(qrange))
    qrange <- unlist(qrange, use.names=FALSE)
    srange <- unlist(range(subject), use.names=FALSE)
    ## bounds violation for each match 
    lviol <- start(qrange) < start(srange)
    rviol <- end(qrange) > end(srange)
    TSSviol <- as.logical(ifelse(qstrand != "-", lviol, rviol))
    TSEviol <- as.logical(ifelse(qstrand != "+", lviol, rviol))
    ## bounds violation across all subjects hit (grouped)
    lviol <- start(qrange) == start(srange)
    rviol <- end(qrange) == end(srange)
    TSSmatch <- as.logical(ifelse(qstrand != "-", lviol, rviol))
    TSEmatch <- as.logical(ifelse(qstrand != "+", lviol, rviol))
    TSSgroup <- (rowsum(as.integer(TSSmatch), qhits)[,1] > 0L)[factor(qhits)]
    TSEgroup <- (rowsum(as.integer(TSEmatch), qhits)[,1] > 0L)[factor(qhits)]

    TSS <- TSSviol & !TSSgroup
    TSE <- TSEviol & !TSEgroup

    DataFrame(TSS=TSS, TSE=TSE)
}

.allMatch <- function(x, idx)
{
    (rowsum(as.integer(x), idx)[,1] == table(idx))[idx]
}

.oneMatch <- function(x, idx)
{
    xcnt <- rowsum(as.integer(x), idx)[,1]
    oneMatch <- rep((xcnt == 1L), table(idx))
    unname(x & oneMatch)
}

.novelExon <- function(splice, intronRegion)
{
    splice_eltNROWS <- elementNROWS(splice)
    if (sum(splice_eltNROWS) == 0L)
        return(logical(length(splice)))

    ## subset on elements with splices
    ans <- logical(length(splice))
    idx <- splice_eltNROWS > 0L
    splice <- splice[idx]
    internal <- unlist(.gaps(splice), use.names=FALSE)
    if (sum(length(internal)) == 0L)
        return(ans)

    ## FIXME : not competely "within"
    hits <- findOverlaps(internal, intronRegion, ignore.strand=TRUE)
    if (length(hits) > 0L) {
        ans0 <- logical(length(splice))
        togroup <- togroup(PartitioningByWidth(splice), j=queryHits(hits))
        ne <- table(togroup) == 1L
        ans0[unique(togroup)] <- ne
        ans[idx] <- ans0
        ans
    } 
    ans
}

.novelSpliceEvent <- function(splice, intron)
{
    splice_eltNROWS <- elementNROWS(splice)
    intron_eltNROWS <- elementNROWS(intron)
    if (sum(splice_eltNROWS) == 0L |
        sum(intron_eltNROWS) == 0L)
        DataFrame(Site=rep.int(FALSE, length(splice)),
                  Junction=rep.int(FALSE, length(splice)))

    ## subset on elements with splices
    site <- junction <- rep.int(FALSE, length(splice))
    idx <- splice_eltNROWS > 0L
    splice <- splice[idx]
    intron <- intron[idx]

    iflat <- unlist(intron, use.names=FALSE)
    sflat <- unlist(splice, use.names=FALSE)
    site[idx] <- .spliceEvent("site", iflat, sflat, splice)
    junction[idx] <- .spliceEvent("junction", iflat, sflat, splice)
    DataFrame(Site=site, Junction=junction)
}

.spliceEvent <- function(type, iflat, sflat, splice)
{
    elt <- togroup(PartitioningByWidth(splice))
    if (type == "site") {
        combiner <- c
        elt <- rep(elt, each=2)
    } else if (type == "junction") {
        combiner <- paste
    }
    ikeys <- paste(seqnames(iflat), combiner(start(iflat),
                   end(iflat)), strand(iflat), sep = ":")
    skeys <- paste(seqnames(sflat), combiner(start(sflat),
                   end(sflat)), strand(sflat), sep = ":")
    novel <- !(skeys %in% ikeys)
    rowsum(as.integer(novel), elt) > 0L
}

.novelRetention <- function(query, intronRegion)
{
    ans <- logical(length(query))
    if (length(query) == 0L)
        return(ans)

    hits <- findOverlaps(unlist(query, use.names=FALSE), intronRegion,
                        ignore.strand=TRUE)
    if (length(hits) > 0L) {
        togroup <- togroup(PartitioningByWidth(query), j=queryHits(hits))
        ans[unique(togroup)] <- TRUE
    }
    ans
}

.intronicRegions <- function(tx, intron) {
    txflt <- unlist(tx, use.names = FALSE)
    intronflt <- unlist(intron, use.names = FALSE)
    regions <- setdiff(intronflt, txflt, ignore.strand = TRUE)
    #map <- findOverlaps(regions, intronflt)
    #mcols(regions)$tx_id <-
    #  splitAsList(names(tx)[togroup(introns)][subjectHits(intronic_to_tx)],
    #              queryHits(intronic_to_tx))
    regions
}

.result <- function(hits, nc=NULL, compatible=NULL, unique=NULL, coding=NULL,
                    strandSpecific=NULL, novelTSS=NULL, novelTSE=NULL,
                    novelSite=NULL, novelJunction=NULL, novelExon=NULL,
                    novelRetention=NULL)
{
    nms <- c("compatible", "unique", "coding", "strandSpecific",
             "novelTSS", "novelTSE", "novelSite", "novelJunction",
             "novelExon", "novelRetention")
    ## full result
    if (!is.null(nc)) {
        mcols(hits) <- DataFrame(compatible, unique, coding, strandSpecific,
                                 novelTSS, novelTSE, novelSite, novelJunction,
                                 novelExon, novelRetention)
        hits
    ## no overlaps 
    } else if (is.null(compatible)) {
        mat <- matrix(logical(0), length(hits), length(nms))
        mcols(hits) <- DataFrame(mat)
        names(mcols(hits)) <- nms
        hits
    ## no compatible overlaps 
    } else {
        mat <- matrix(FALSE, length(hits), length(nms))
        mcols(hits) <- DataFrame(cbind(compatible, unique, coding,
                                  strandSpecific, mat))
        names(mcols(hits)) <- nms
        hits
    }
}

.insertGaps <- function(reads)
{
    query.break <- mcols(reads)$query.break
    if (is.null(query.break))
      stop("missing 'query.break' metadata variable: reads not paired?")
    reads_flat <- unlist(reads, use.names = FALSE)
    reads_part <- PartitioningByWidth(reads)
    left_end <- start(reads_part) + query.break - 1L
    right_start <- left_end + 1L
    start <- end(reads_flat)[left_end]
    end <- pmax(start(reads_flat)[right_start], start - 1L)
    if (any(seqnames(reads_flat)[left_end] !=
                                 seqnames(reads_flat)[right_start]))
      stop("reads are on different chromosomes")
    GRanges(seqnames(reads_flat)[left_end],
            IRanges(start, end), strand(reads_flat)[left_end])
}

## Until we have the formal 'gaps' method for GRangeList
.isNumericOrNAs <- S4Vectors:::isNumericOrNAs
.gaps <- function(x, start=NA, end=NA)
{
    if (!.isNumericOrNAs(start))
        stop("'start' must be an integer vector or NA")
    if (!is.integer(start))
        start <- as.integer(start)
    if (!.isNumericOrNAs(end))
        stop("'end' must be an integer vector or NA")
    if (!is.integer(end))
        end <- as.integer(end)

    ## seqname and strand consistent in list elements
    if (all(elementNROWS(runValue(seqnames(x))) == 1L) &&
        all(elementNROWS(runValue(strand(x))) == 1L)) {
        flat <- unlist(x, use.names=FALSE)
        gaps <- gaps(ranges(x), start, end)
### FIXME: this makes this function more of an 'introns' than a .gaps.
### FIXME: this breaks when the GRangesList is not ordered by position
        if (!is.null(mcols(x)$query.break)) {
          insert_gaps <- as(ranges(.insertGaps(x)), "RangesList")
          gaps <- setdiff(gaps, insert_gaps)
        }

        idx <- elementNROWS(gaps) != 0
        ## FIXME : can't handle lists with empty elements 
        ##         'start' and 'end' not quite right here
        firstseg <- start(PartitioningByWidth(x))
        seqnms <- rep(seqnames(flat)[firstseg], elementNROWS(gaps))
        strand <- rep(strand(flat)[firstseg], elementNROWS(gaps))
        gr <- relist(GRanges(seqnms, unlist(gaps, use.names=FALSE), strand), gaps)
        gr
    } else {
### FIXME: does not handle query.break column yet
        setdiff(range(x), x)
    }

}

.findSpliceOverlaps <- function(query, subject, ignore.strand=FALSE, cds=NULL)
{
    ## adjust strand based on 'XS'
    if (!is.null(xs <- mcols(query)$XS)) {
        strand <- ifelse(!is.na(xs), xs, "*")
        strand(query) <- relist(Rle(strand, elementNROWS(query)),
                                query)
    }
    ## NOTE: this misses reads completely within an intron, but this
    ## is intentional: a read is only assigned to a transcript if it
    ## hits an exon. Otherwise, it could be from another gene inside
    ## an intron (happens frequently).
    olap <- findOverlaps(query, subject, ignore.strand=ignore.strand)
    if (length(olap) == 0L)
        return(.result(olap))
    if (!is.null(cds)) {
        coding <- logical(length(olap))
        hits <- findOverlaps(query, cds, ignore.strand=ignore.strand)
        coding[queryHits(olap) %in% queryHits(hits)] <- TRUE
    } else {
        coding <- rep.int(NA, length(olap))
    }

    query <- query[queryHits(olap)]
    subject <- subject[subjectHits(olap)]
    splice <- .gaps(query)

    compatible <- .compatibleTranscription(query, subject, splice)
    unique <- .oneMatch(compatible, queryHits(olap))
    strandSpecific <- all(strand(query) != "*")
    mcols(olap) <- DataFrame(compatible, unique, coding, strandSpecific)
    olap
}


setGeneric("findSpliceOverlaps", signature=c("query", "subject"),
    function(query, subject, ignore.strand=FALSE, ...)
        standardGeneric("findSpliceOverlaps")
)

setMethod("findSpliceOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    .findSpliceOverlaps(query, subject, ignore.strand, cds=cds)
})

setMethod("findSpliceOverlaps", c("GAlignments", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    findSpliceOverlaps(grglist(query, order.as.in.query=TRUE), subject,
                       ignore.strand, ..., cds=cds)
})

setMethod("findSpliceOverlaps", c("GAlignmentPairs", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
### FIXME:
### instead of relying on query.break column, maybe we should add a
### 'splice = .gaps(query)' argument to .findSpliceOverlaps that we
### set to junctions(query) here. The downside is that a GRangesList
### derived from GAlignmentPairs will no longer work.
    findSpliceOverlaps(grglist(query), subject,
                       ignore.strand, ..., cds=cds)
})

setMethod("findSpliceOverlaps", c("character", "ANY"),
          function(query, subject,
                   ignore.strand=FALSE, ...,
                   param=ScanBamParam(), singleEnd=TRUE)
{
    findSpliceOverlaps(BamFile(query), subject,
                       ignore.strand, ...,
                       param=param, singleEnd=singleEnd)
})

setMethod("findSpliceOverlaps", c("BamFile", "ANY"),
    function(query, subject,
             ignore.strand=FALSE, ...,
             param=ScanBamParam(), singleEnd=TRUE)
{
    findSpliceOverlaps(.readRanges(query, param, singleEnd), subject,
                       ignore.strand, ...)
})

.readRanges <- function(bam, param, singleEnd)
{
    if (!"XS" %in% bamTag(param))
        bamTag(param) <- c(bamTag(param), "XS")
    if (singleEnd)
        reads <- readGAlignments(bam, param=param)
    else {
        reads <- readGAlignmentPairs(path(bam), param=param)
        first_xs <- mcols(first(reads))$XS
        last_xs <- mcols(last(reads))$XS
        if (!is.null(first_xs) && !is.null(last_xs)) {
            xs <- first_xs
            xs[is.na(xs)] <- last_xs[is.na(xs)]
            mcols(reads)$XS <- xs
        }
    }

    metadata(reads)$bamfile <- bam

    reads
}

