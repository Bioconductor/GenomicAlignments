test_readGAlignments <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    which <- RangesList(seq1=IRanges(1, width=100))
    param <- ScanBamParam(which=which)
    result <- readGAlignments(fl, param=param)
    checkTrue(validObject(result))
    checkIdentical(c(seq1=1575L, seq2=1584L), seqlengths(result))
}

test_readGAlignments_missing_param <- function()
{
    fl <- system.file("unitTests", "cases", "ex1_noindex.bam",
                      package="Rsamtools")
    result0 <- readGAlignments(fl)
    checkTrue(validObject(result0))

    bf <- open(BamFile(fl, character()))
    result1 <- readGAlignments(bf)
    checkIdentical(result1, result0)
}

test_readGAlignments_length0 <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

    which <- RangesList(seq1=IRanges(100000, width=100))
    param <- ScanBamParam(which=which)
    result <- readGAlignments(fl, param=param)
    checkTrue(validObject(result))

    which <- RangesList(seq1=IRanges(c(1, 100000), width=100))
    param <- ScanBamParam(which=which)
    result <- readGAlignments(fl, param=param)
    checkTrue(validObject(result))
}

test_readGAlignments_tag <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

    ## valid
    param <- ScanBamParam(tag=("NM"))
    gal <- readGAlignments(fl, param=param)
    checkIdentical(924L, sum(mcols(gal)[["NM"]]))

    ## empty
    param <- ScanBamParam(tag=("FO"))
    gal <- readGAlignments(fl, param=param)
    checkIdentical(rep.int(NA, length(gal)), mcols(gal)[["FO"]])
}

test_readGAlignments_BamViews <- function()
{
    checkTrue(validObject(readGAlignments(BamViews())))
    src <- system.file("unitTests", "cases", package="Rsamtools")
    fl <- c(system.file("extdata", "ex1.bam", package="Rsamtools"),
            file.path(src, "ex1_shuf1000.bam"))
    bv <- BamViews(fl, auto.range=TRUE)
    rng <- bamRanges(bv)
    aln <- readGAlignments(bv)
    checkEquals(length(bamPaths(bv)), length(aln))

    fl <- c(fl, tempfile())
    bv <- BamViews(fl, bamRanges=rng)
    current <- suppressWarnings({
        tryCatch({
            aln <- readGAlignments(bv)
        }, error=identity)
    })
    checkTrue(is(current, "bperror"))
}
