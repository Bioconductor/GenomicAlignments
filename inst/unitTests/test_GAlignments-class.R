test_GAlignments_constructor <- function()
{
    checkTrue(validObject(GAlignments()))
    checkTrue(validObject(GAlignments("A", pos=1, cigar="1M")))
    checkTrue(validObject(GAlignments("A", pos=1, cigar="1M", strand="-")))
}

test_GAlignments_seqlevels <- function()
{
    gal0 <- GAlignments(seqnames=c("chr1", "chr2"),
                        pos=c(10, 100),
                        cigar=c("50M", "50M"))

    ## Drop
    gal <- gal0
    seqlevels(gal, pruning.mode="coarse") <- "chr2"
    checkIdentical("chr2", seqlevels(gal))

    ## Rename
    gal <- gal0
    seqlevels(gal)[seqlevels(gal) == "chr2"] <- "2"
    checkIdentical(c("chr1", "2"),  seqlevels(gal))
}

test_GAlignments_combine <- function() 
{
    galn <- GAlignments(seqnames="A", pos=1, cigar="1M", strand="-")
    galn_c <- GAlignments(seqnames=Rle(factor("A"), 2),
                          pos=rep(1, 2), cigar=rep("1M", 2),
                          strand=Rle(strand("-"), 2))
    checkIdentical(galn_c, c(galn, galn))
}

