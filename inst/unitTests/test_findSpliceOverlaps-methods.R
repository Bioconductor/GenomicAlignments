## TODO: Add tests for "findSpliceOverlaps" methods defined in the
## GenomicAlignments package i.e.:
##   findSpliceOverlaps,GAlignments,GRangesList
##   findSpliceOverlaps,GAlignmentPairs,GRangesList
##   findSpliceOverlaps,character,ANY
##   findSpliceOverlaps,BamFile,ANY

.extract <- function(x, col) as.logical(mcols(x)[[col]])

#test_findSpliceOverlaps_novelJunction <- function()
#{
#    ## novel junction, no novel sites
#    genes <- GRangesList(
#        GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
#        GRanges("chr1", IRanges(c(5, 22), c(15, 25)), "+"))
#
#    ## query = GAlignments
#    gal <- GAlignments("chr1", 5L, "11M4N6M", strand("+"))
#    GALres <- findSpliceOverlaps(gal, genes)
#    checkIdentical(c(TRUE, TRUE), .extract(GALres, "novelJunction"))
#    checkIdentical(c(FALSE, FALSE), .extract(GALres, "novelSite"))
#
#    ## query = GAlignmentPairs
#    gal1 <- GAlignments("chr1", 5L, "11M4N6M", strand("+"))
#    gal2 <- GAlignments("chr1", 50L, "6M", strand("-"))
#    galp <- GAlignmentPairs(gal1, gal2, TRUE)
#    GALPres <- findSpliceOverlaps(galp, genes)
#    checkIdentical(c(TRUE, TRUE), .extract(GALPres, "novelJunction"))
#    checkIdentical(c(FALSE, FALSE), .extract(GALPres, "novelSite"))
#}

