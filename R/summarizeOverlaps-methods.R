### =========================================================================
### summarizeOverlaps() generic and methods
### -------------------------------------------------------------------------


setGeneric("summarizeOverlaps", signature=c("features", "reads"),
    function(features, reads, mode=Union, ignore.strand=FALSE, ...)
        standardGeneric("summarizeOverlaps")
)


### -------------------------------------------------------------------------
### Methods for GAlignments, GAlignmentsList and GAlignmentPairs objects
###

.dispatchOverlaps <- function(features, reads, mode,
                              ignore.strand,
                              inter.feature, preprocess.reads, ...)
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
    if (!is.null(preprocess.reads))
        reads <- preprocess.reads(reads, ...)
    mode(features, reads,
         ignore.strand=ignore.strand, inter.feature=inter.feature)
}

.summarizeOverlaps <- function(features, reads, mode=Union,
                               ignore.strand=FALSE,
                               inter.feature=TRUE, preprocess.reads=NULL, ...)
{
    if (class(reads) == "GRangesList") {
        if (all(unlist(strand(reads), use.names=FALSE) == "*"))
            ignore.strand <- TRUE
    } else {
        if (all(strand(reads) == "*"))
            ignore.strand <- TRUE
    }
    mode <- match.fun(mode)
    counts <- .dispatchOverlaps(features, reads, mode,
                                ignore.strand,
                                inter.feature, preprocess.reads, ...)
    colData <- DataFrame(object=class(reads),
                         records=length(reads),
                         row.names="reads")
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
                         rowRanges=features, colData=colData)
}

setMethod("summarizeOverlaps", c("GRanges", "GAlignments"),
    .summarizeOverlaps
)

setMethod("summarizeOverlaps", c("GRangesList", "GAlignments"),
    .summarizeOverlaps
)

setMethod("summarizeOverlaps", c("GRanges", "GAlignmentsList"),
    .summarizeOverlaps
)

setMethod("summarizeOverlaps", c("GRangesList", "GAlignmentsList"),
    .summarizeOverlaps
)

setMethod("summarizeOverlaps", c("GRanges", "GAlignmentPairs"),
    .summarizeOverlaps
)

setMethod("summarizeOverlaps", c("GRangesList", "GAlignmentPairs"),
    .summarizeOverlaps
)

setMethod("summarizeOverlaps", c("GRanges", "GRanges"),
    .summarizeOverlaps
)

setMethod("summarizeOverlaps", c("GRangesList", "GRanges"),
    .summarizeOverlaps
)

setMethod("summarizeOverlaps", c("GRanges", "GRangesList"),
    .summarizeOverlaps
)

setMethod("summarizeOverlaps", c("GRangesList", "GRangesList"),
    .summarizeOverlaps
)

### -------------------------------------------------------------------------
### 'mode' functions 
###

Union <- function(features, reads,
                  ignore.strand=FALSE, inter.feature=TRUE)
{
    ov <- findOverlaps(features, reads,
                       ignore.strand=ignore.strand)
    if (inter.feature) {
        ## Remove ambigous reads.
        reads_to_keep <- which(countSubjectHits(ov) == 1L)
        ov <- ov[subjectHits(ov) %in% reads_to_keep]
    }
    countQueryHits(ov)
}

IntersectionStrict <- function(features, reads,
                               ignore.strand=FALSE, inter.feature=TRUE)
{
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
    ov <- findOverlaps(features, regions,
                       ignore.strand=ignore.strand)
    regions_to_keep <- which(countSubjectHits(ov) == 1L)
    ov <- ov[subjectHits(ov) %in% regions_to_keep]
    unlisted_ans <- regions[subjectHits(ov)]
    ans_partitioning <- as(ov, "PartitioningByEnd")
    relist(unlisted_ans, ans_partitioning)
}

IntersectionNotEmpty <-  function(features, reads,
                                  ignore.strand=FALSE, inter.feature=TRUE)
{
    features <- .removeSharedRegions(features, ignore.strand=ignore.strand)
    Union(features, reads,
          ignore.strand=ignore.strand, inter.feature=inter.feature)
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
        FUN <- readGAlignments
    } else {
        if (fragments)
            FUN <- readGAlignmentsList
        else
            FUN <- readGAlignmentPairs
    }
    FUN
}

.countWithYieldSize <- function(FUN, features, bf, mode,
                                ignore.strand,
                                inter.feature, param, preprocess.reads, ...)
{
    if (is.na(yieldSize(bf))) {
        x <- FUN(bf, param=param, ...)
        .dispatchOverlaps(features, x, mode,
                          ignore.strand,
                          inter.feature, preprocess.reads, ...)
    } else {
        if (!isOpen(bf)) {
            open(bf)
            on.exit(close(bf))
        }
        ct <- integer(length(features))
        while (length(x <- FUN(bf, param=param, ...))) {
            ct <- ct + .dispatchOverlaps(features, x, mode,
                                         ignore.strand,
                                         inter.feature, preprocess.reads, ...)
        }
        ct
    }
}

.dispatchBamFiles <-
    function(features, reads, mode, ignore.strand,
             count.mapped.reads=FALSE, inter.feature=TRUE, 
             singleEnd=TRUE, fragments=FALSE,
             param=ScanBamParam(), preprocess.reads=NULL, ...)
{
    exist <- sapply(reads, function(bf) file.exists(path(bf)))
    if (!all(exist))
        stop(paste0("file(s): ", paste(path(reads)[!exist], collapse=","), 
                    " do not exist"))
    FUN <- .getReadFunction(singleEnd, fragments)

    cts <- bplapply(setNames(seq_along(reads), names(reads)),
               function(i, FUN, reads, features, mode,
                        ignore.strand,
                        inter.feature, param, preprocess.reads, ...) {
                   bf <- reads[[i]]
                   .countWithYieldSize(FUN, features, bf, mode,
                                       ignore.strand,
                                       inter.feature, param, 
                                       preprocess.reads, ...)
               }, FUN, reads, features, mode=match.fun(mode),
               ignore.strand=ignore.strand,
               inter.feature=inter.feature, param=param, 
               preprocess.reads=preprocess.reads, ...
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
                         rowRanges=features, colData=colData)
}

setMethod("summarizeOverlaps", c("GRanges", "BamFile"),
    function(features, reads, mode=Union,
             ignore.strand=FALSE,
             inter.feature=TRUE, singleEnd=TRUE,
             fragments=FALSE, param=ScanBamParam(), preprocess.reads=NULL, ...)
{
    .checkArgs(reads, singleEnd, fragments)
    .dispatchBamFiles(features, BamFileList(reads), mode,
                      ignore.strand,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments, param=param, 
                      preprocess.reads=preprocess.reads, ...)
})

setMethod("summarizeOverlaps", c("GRangesList", "BamFile"),
    function(features, reads, mode=Union,
             ignore.strand=FALSE,
             inter.feature=TRUE, singleEnd=TRUE,
             fragments=FALSE, param=ScanBamParam(), preprocess.reads=NULL, ...)
{
    .checkArgs(reads, singleEnd, fragments)
    .dispatchBamFiles(features, BamFileList(reads), mode,
                      ignore.strand,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments, param=param,
                      preprocess.reads=preprocess.reads, ...)
})

.summarizeOverlaps_character <-
    function(features, reads, mode=Union,
             ignore.strand=FALSE,
             yieldSize=1000000L, inter.feature=TRUE, singleEnd=TRUE,
             fragments=FALSE, param=ScanBamParam(), preprocess.reads=NULL, ...)
{
    
    if (!all(file.exists(reads)))
        stop("file(s) do not exist:\n  ",
             paste(reads[!file.exists(reads)], collapse="\n  "))
    if (is.null(names(reads))) {
        if (any(duplicated(reads)))
            stop("duplicate 'reads' paths not allowed; use distinct names()")
    } else if (any(duplicated(names(reads))))
        stop("duplicate 'names(reads)' file paths not allowed")
    reads <- BamFileList(reads, yieldSize=yieldSize, obeyQname=FALSE,
                         asMates=!singleEnd)
    summarizeOverlaps(features, reads, mode,
                      ignore.strand=ignore.strand,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments, param=param,
                      preprocess.reads=preprocess.reads, ...)
}

setMethod("summarizeOverlaps", c("GRanges", "character"),
    .summarizeOverlaps_character
)

setMethod("summarizeOverlaps", c("GRangesList", "character"),
    .summarizeOverlaps_character
)

.summarizeOverlaps_BamFileList <-
    function(features, reads, mode=Union,
             ignore.strand=FALSE,
             inter.feature=TRUE, singleEnd=TRUE,
             fragments=FALSE, param=ScanBamParam(), preprocess.reads=NULL, ...)
{
    if (any(duplicated(names(reads))))
        stop("duplicate 'names(reads)' not allowed")
    .checkArgs(reads, singleEnd, fragments)
    .dispatchBamFiles(features, reads, mode,
                      ignore.strand,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments, param=param, 
                      preprocess.reads=preprocess.reads, ...)
}

setMethod("summarizeOverlaps", c("GRanges", "BamFileList"),
    .summarizeOverlaps_BamFileList
)

setMethod("summarizeOverlaps", c("GRangesList", "BamFileList"),
    .summarizeOverlaps_BamFileList
)

setMethod("summarizeOverlaps", c("BamViews", "missing"),
    function(features, reads, mode=Union,
             ignore.strand=FALSE,
             inter.feature=TRUE, singleEnd=TRUE,
             fragments=FALSE, param=ScanBamParam(), preprocess.reads=NULL, ...)
{
    se <- callGeneric(bamRanges(features), BamFileList(bamPaths(features)),
                      mode=mode,
                      ignore.strand=ignore.strand,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments, param=param,
                      preprocess.reads=preprocess.reads, ...)
    colData(se)$bamSamples <- bamSamples(features)
    colData(se)$bamIndices <- bamIndicies(features)
    metadata(se)$bamExperiment <- bamExperiment(features)
    se
})

