.noGaps <- GAlignments(
    Rle(factor(c("chr1", "chr2", "chr1", "chr3")),
        c(1, 3, 2, 4)),
    pos=1:10, cigar=paste0(10:1, "M"),
    strand=Rle(strand(c("-", "+", "*", "+", "-")),
        c(1, 2, 2, 3, 2)),
    names=head(letters, 10), score=1:10)
.Gaps <- GAlignments(
    Rle(factor(c("chr2", "chr4")), c(3, 4)), pos=1:7,
    cigar=c("5M", "3M2N3M2N3M", "5M", "10M", "5M1N4M", "8M2N1M", "5M"),
    strand=Rle(strand(c("-", "+")), c(4, 3)),
    names=tail(letters, 7), score=1:7)

test_GAlignments_qnarrow <- function()
{
    gal <- GAlignments(seqnames=rep(factor("A"), 8),
                       pos=10:17,
                       cigar=c("5M", "5X", "3M2I3M", "3M2D3M",
                               "3M2N3M", "3M2S3M", "3M2H3M", "3M2P3M"),
                       strand=Rle(strand(rep("+", 8))))
    n1 <- narrow(gal, start=3)
    q1 <- qnarrow(gal, start=3)
    checkIdentical(qwidth(n1), qwidth(q1))
    checkIdentical(width(n1), width(q1))

    n2 <- narrow(gal, start=4)
    q2 <- qnarrow(gal, start=4)
    checkIdentical(width(n2), width(q2))
    ## M and X 
    checkIdentical(qwidth(n2[1:2]), qwidth(q2[1:2]))
    ## I 
    checkIdentical(qwidth(q2[3]), width(q2[3]) + 2L)
    ## D, N and P
    checkIdentical(qwidth(q2[c(4,5,8)]), width(q2[c(4,5,8)]))
    ## S and H
    checkIdentical(qwidth(q2[6]), width(q2[6]) + 2L)
    checkIdentical(qwidth(q2[7]), width(q2[7]))
}

test_GAlignmentsList_qnarrow <- function()
{
    galist <- GAlignmentsList(.noGaps[1:6], .Gaps)
    qn <- qnarrow(galist, end=-4)
    checkIdentical(qnarrow(galist[[1]], end=-4), qn[[1]])
    checkIdentical(qnarrow(galist[[2]], end=-4), qn[[2]])
}

