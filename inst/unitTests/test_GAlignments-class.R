test_GAlignments_constructor <- function()
{
    checkTrue(validObject(GAlignments()))
    checkTrue(validObject(GAlignments(seqnames=factor("A"),
                                      pos=1L, cigar="1M",
                                      strand=strand("-"))))
}

test_GAlignments_seqlevels <- function()
{
    gal0 <- GAlignments(seqnames=Rle(c("chr1", "chr2")),
                        pos=as.integer(c(10, 100)),
                        cigar=c("50M", "50M"),
                        strand=strand(c("*", "*")))

    ## Drop
    gal <- gal0
    seqlevels(gal, force=TRUE) <- "chr2"
    checkIdentical("chr2", seqlevels(gal))

    ## Rename
    gal <- gal0
    seqlevels(gal)[seqlevels(gal) == "chr2"] <- "2"
    checkIdentical(c("chr1", "2"),  seqlevels(gal))
}

test_GAlignments_combine <- function() 
{
    galn <- GAlignments(seqnames=factor("A"),
                        pos=1L, cigar="1M",
                        strand=strand("-"))
    galn_c <- GAlignments(seqnames=rep(factor("A"), 2),
                          pos=rep(1L, 2), cigar=rep("1M", 2),
                          strand=rep(strand("-"), 2))
    checkIdentical(galn_c, c(galn, galn))
}

