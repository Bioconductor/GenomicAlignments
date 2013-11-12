### =========================================================================
### summarizeOverlaps() generic and methods
### -------------------------------------------------------------------------


setGeneric("summarizeOverlaps", signature=c("features", "reads"),
    function(features, reads, mode=Union, ignore.strand=FALSE, ...)
{
    standardGeneric("summarizeOverlaps")
})


### -------------------------------------------------------------------------
### Methods for GAlignments, GAlignmentsList and GAlignmentPairs objects
###

.dispatchOverlaps <-
    function(features, reads, mode, ignore.strand, inter.feature, ...)
{
    if (ignore.strand) {
        if (class(features) == "GRangesList") {
            r <- unlist(features)
            strand(r) <- "*"
            features@unlistData <- r
        } else {
            strand(features) <- "*"
        }
    }
    mode(features, reads, ignore.strand, inter.feature=inter.feature)
}

.summarizeOverlaps <- function(features, reads, mode, ignore.strand,
                               ..., inter.feature=inter.feature)
{
    if (all(strand(reads) == "*"))
        ignore.strand <- TRUE
    mode <- match.fun(mode)
    counts <- .dispatchOverlaps(features, reads, mode, ignore.strand,
                                inter.feature=inter.feature)
    colData <- DataFrame(object=class(reads),
                         records=length(reads),
                         row.names="reads")
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
                         rowData=features, colData=colData)
}

setMethod("summarizeOverlaps", c("GRanges", "GAlignments"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, reads, mode, ignore.strand,
                       ..., inter.feature=inter.feature)
})

setMethod("summarizeOverlaps", c("GRangesList", "GAlignments"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, reads, mode, ignore.strand,
                       ..., inter.feature=inter.feature)
})

setMethod("summarizeOverlaps", c("GRanges", "GAlignmentsList"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, grglist(reads, ignore.strand=TRUE),
                       mode, ignore.strand, ...,
                       inter.feature=inter.feature)
})

setMethod("summarizeOverlaps", c("GRangesList", "GAlignmentsList"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, grglist(reads, ignore.strand=TRUE),
                       mode, ignore.strand, ...,
                       inter.feature=inter.feature)
})

setMethod("summarizeOverlaps", c("GRanges", "GAlignmentPairs"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, grglist(reads), mode, ignore.strand,
                       ..., inter.feature=inter.feature)
})

setMethod("summarizeOverlaps", c("GRangesList", "GAlignmentPairs"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, grglist(reads), mode, ignore.strand,
                       ..., inter.feature=inter.feature)
})


### -------------------------------------------------------------------------
### 'mode' functions 
###

Union <- function(features, reads, ignore.strand=FALSE, inter.feature=TRUE)
{
    ov <- findOverlaps(features, reads, ignore.strand=ignore.strand)
    if (inter.feature) {
        ## Remove ambigous reads.
        reads_to_keep <- which(countSubjectHits(ov) == 1L)
        ov <- ov[subjectHits(ov) %in% reads_to_keep]
    }
    countQueryHits(ov)
}

### Drop from 'reads' circular seqlevels that are in use in *both*: 'reads'
### and 'features'.
.dropCircularSeqlevelsInUse <- function(reads, features)
{
    seqlevels_in_use <- intersect(seqlevelsInUse(reads),
                                  seqlevelsInUse(features))
    seqinfo <- merge(seqinfo(reads), seqinfo(features))
    is_circ <- isCircular(seqinfo)
    circular_seqlevels <- names(is_circ)[is_circ]
    seqlevels_to_drop <- intersect(seqlevels_in_use, circular_seqlevels)
    if (length(seqlevels_to_drop) != 0L) {
        warning("reads on circular sequence(s) '",
                paste(seqlevels_to_drop, sep="', '"),
                "' were ignored")
        seqlevels(reads, force=TRUE) <- setdiff(seqlevels(reads),
                                                seqlevels_to_drop)
    }
    reads
}

IntersectionStrict <- function(features, reads, ignore.strand=FALSE,
                               inter.feature=TRUE)
{
    reads <- .dropCircularSeqlevelsInUse(reads, features)
    ov <- findOverlaps(reads, features, type="within",
                       ignore.strand=ignore.strand)
    if (inter.feature) {
        ## Remove ambigous reads.
        reads_to_keep <- which(countQueryHits(ov) == 1L)
        ov <- ov[queryHits(ov) %in% reads_to_keep]
    }
    countSubjectHits(ov)
}

.removeSharedRegions <- function(features, ignore.strand=FALSE)
{
    if (is(features, "GRanges")) {
        regions <- disjoin(features, ignore.strand=ignore.strand)
    } else if (is(features, "GRangesList")) {
        regions <- disjoin(features@unlistData, ignore.strand=ignore.strand)
    } else {
        stop("internal error")  # should never happen
    }
    ov <- findOverlaps(features, regions, ignore.strand=ignore.strand)
    regions_to_keep <- which(countSubjectHits(ov) == 1L)
    ov <- ov[subjectHits(ov) %in% regions_to_keep]
    ans_flesh <- regions[subjectHits(ov)]
    ## Using 'countQueryHits(ov)' to compute the skeleton of 'ans' relies
    ## on the assumption that the hits returned by findOverlaps() are always
    ## ordered by query.
    ans_eltlens <- countQueryHits(ov)
    ans_skeleton <- PartitioningByEnd(cumsum(ans_eltlens))
    relist(ans_flesh, ans_skeleton)
}

IntersectionNotEmpty <-  function(features, reads, ignore.strand=FALSE,
                                  inter.feature=TRUE)
{
    features <- .removeSharedRegions(features, ignore.strand=ignore.strand)
    Union(features, reads, ignore.strand=ignore.strand,
           inter.feature=inter.feature)
}


### -------------------------------------------------------------------------
### Methods for BamFiles and BamViews objects
###

.checkArgs <- function(bam, singleEnd, fragments)
{
    if (singleEnd) {
        if (all(isTRUE(asMates(bam))))
            stop("cannot specify both 'singleEnd=TRUE' and 'asMates=TRUE'")
        if (fragments)
            stop("when 'fragments=TRUE', 'singleEnd' should be FALSE")
    ## all paired-end reading now goes through new C algo
    } else {
        asMates(bam) <- TRUE
    }
}

.getReadFunction <- function(singleEnd, fragments)
{
    if (singleEnd) {
        FUN <- readGAlignmentsFromBam
    } else {
        if (fragments)
            FUN <- readGAlignmentsListFromBam
        else
            FUN <- readGAlignmentPairsFromBam
    }
    FUN
}

.countWithYieldSize <- function(FUN, features, bf, mode, ignore.strand,
                                inter.feature, param)
{
    if (is.na(yieldSize(bf))) {
        x <- FUN(bf, param=param)
        .dispatchOverlaps(features, x, mode, ignore.strand, inter.feature)
    } else {
        if (!isOpen(bf)) {
            open(bf)
            on.exit(close(bf))
        }
        ct <- integer(length(features))
        while (length(x <- FUN(bf, param=param))) {
            ct <- ct + .dispatchOverlaps(features, x, mode, ignore.strand,
                                         inter.feature)
        }
        ct
    }
}

.dispatchBamFiles <-
    function(features, reads, mode, ignore.strand, ...,
             count.mapped.reads=FALSE,
             inter.feature=TRUE, singleEnd=TRUE, fragments=FALSE,
             param=ScanBamParam())
{
    FUN <- .getReadFunction(singleEnd, fragments)

    if ("package:parallel" %in% search() & .Platform$OS.type != "windows")
        lapply <- parallel::mclapply

    cts <- lapply(setNames(seq_along(reads), names(reads)),
               function(i, FUN, reads, features, mode, ignore.strand,
                        inter.feature, param) {
                   bf <- reads[[i]]
                   .countWithYieldSize(FUN, features, bf, mode, ignore.strand,
                                       inter.feature, param)
               }, FUN, reads, features, mode=match.fun(mode), ignore.strand,
               inter.feature, param
           )

    counts <- as.matrix(do.call(cbind, cts))
    if (count.mapped.reads) {
        countBam <- countBam(reads)
        flag <- scanBamFlag(isUnmappedQuery=FALSE)
        param <- ScanBamParam(flag=flag, what="seq")
        colData <- DataFrame(countBam[c("records", "nucleotides")],
                             mapped=countBam(reads, param=param)$records,
                             row.names=colnames(counts))
    } else {
        colData <- DataFrame(row.names=colnames(counts))
    }
    SummarizedExperiment(assays=SimpleList(counts=counts),
                         rowData=features, colData=colData)
}

setMethod("summarizeOverlaps", c("GRanges", "BamFile"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE, singleEnd=TRUE, fragments=FALSE,
             param=ScanBamParam())
{
    .checkArgs(reads, singleEnd, fragments)
    .dispatchBamFiles(features, BamFileList(reads), mode, ignore.strand, ...,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments, param=param)
})

setMethod("summarizeOverlaps", c("GRangesList", "BamFile"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE, singleEnd=TRUE, fragments=FALSE,
             param=ScanBamParam())
{
    .checkArgs(reads, singleEnd, fragments)
    .dispatchBamFiles(features, BamFileList(reads), mode, ignore.strand, ...,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments, param=param)
})

.summarizeOverlaps_character <-
    function(features, reads, mode, ignore.strand=FALSE, ...,
        yieldSize=1000000L, inter.feature=TRUE, singleEnd=TRUE,
        fragments=FALSE, param=ScanBamParam())
{
    if (is.null(names(reads))) {
        if (any(duplicated(reads)))
            stop("duplicate 'reads' paths not allowed; use distinct names()")
    } else if (any(duplicated(names(reads))))
        stop("duplicate 'names(reads)' file paths not allowed")
    reads <- BamFileList(reads, yieldSize=yieldSize, obeyQname=FALSE,
                         asMates=!singleEnd)
    summarizeOverlaps(features, reads, mode, ignore.strand, ...,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments)
}

setMethod("summarizeOverlaps", c("GRanges", "character"),
    .summarizeOverlaps_character)

setMethod("summarizeOverlaps", c("GRangesList", "character"),
    .summarizeOverlaps_character)

.summarizeOverlaps_BamFileList <-
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE, singleEnd=TRUE, fragments=FALSE,
             param=ScanBamParam())
{
    if (any(duplicated(names(reads))))
        stop("duplicate 'names(reads)' not allowed")
    .checkArgs(reads, singleEnd, fragments)
    .dispatchBamFiles(features, reads, mode, ignore.strand, ...,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments, param=param)
}

setMethod("summarizeOverlaps", c("GRanges", "BamFileList"),
    .summarizeOverlaps_BamFileList)

setMethod("summarizeOverlaps", c("GRangesList", "BamFileList"),
    .summarizeOverlaps_BamFileList)

setMethod("summarizeOverlaps", c("BamViews", "missing"),
function(features, reads, mode, ignore.strand=FALSE,
         ..., inter.feature=TRUE, singleEnd=TRUE, fragments=FALSE,
         param=ScanBamParam())
{
    se <- callGeneric(bamRanges(features), BamFileList(bamPaths(features)),
                      mode, ignore.strand, ..., inter.feature=inter.feature,
                      singleEnd=singleEnd, fragments=fragments, param=param)
    colData(se)$bamSamples <- bamSamples(features)
    colData(se)$bamIndicies <- bamIndicies(features)
    exptData(se)$bamExperiment <- bamExperiment(features)
    se
})

