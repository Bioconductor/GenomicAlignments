### =========================================================================
### "summarizeOverlaps" methods
### -------------------------------------------------------------------------


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

.dispatchOverlaps <- GenomicRanges:::.dispatchOverlaps
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

