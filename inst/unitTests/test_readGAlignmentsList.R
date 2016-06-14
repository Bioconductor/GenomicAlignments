library(pasillaBamSubset)
chr4 <- untreated3_chr4()

test_readGAlignmentsList_construction <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    bf <- BamFile(fl, asMates=TRUE)
    galist <- readGAlignmentsList(fl)
    checkTrue(is.null(names(galist)))
    galist <- readGAlignmentsList(fl, use.names=TRUE)
    target <- c("EAS54_61:4:143:69:578", "EAS219_FC30151:7:51:1429:1043")
    checkIdentical(names(galist)[1:2], target)

    ## first segment first
    param <- ScanBamParam(what="flag")
    galist <- readGAlignmentsList(fl, param=param)
    mates <- galist[mcols(galist)$mate_status == "mated"]
    flagBit <- bamFlagAsBitMatrix(mcols(unlist(mates))$flag,
                                  bitnames="isFirstMateRead") 
    m <- matrix(flagBit, nrow=2)
    checkIdentical(c(1572L, 0), rowSums(m))
}

test_readGAlignmentsList_noYieldSize <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    bf <- BamFile(fl, asMates=TRUE)
    galist <- readGAlignmentsList(fl)
    checkTrue(validObject(galist))
}

test_readGAlignmentsList_yieldSize <- function()
{
    bf <- BamFile(chr4, asMates=TRUE, yieldSize=1)
    scn1 <- scanBam(bf)
    galist1 <- readGAlignmentsList(bf)
    checkTrue(length(scn1[[1]]$qname) == 2)
    checkTrue(length(unique(scn1[[1]]$qname)) == 1)
    checkTrue(length(unique(scn1[[1]]$qname)) == length(galist1))

    bf <- BamFile(chr4, asMates=TRUE, yieldSize=2)
    scn2 <- scanBam(bf)
    galist2 <- readGAlignmentsList(bf)
    checkTrue(length(scn2[[1]]$qname) == 4)
    checkTrue(length(unique(scn2[[1]]$qname)) == 2)
    checkTrue(length(unique(scn2[[1]]$qname)) == length(galist2))
}

test_readGAlignmentsList_mcols <- function()
{
    bf <- BamFile(chr4, asMates=TRUE, yieldSize=100)
    param <- ScanBamParam(tag=("NM"))
    galist <- readGAlignmentsList(bf, param=param)
    checkIdentical(colnames(mcols(unlist(galist))), "NM")
    checkTrue(names(mcols(galist)) == "mate_status")

    param <- ScanBamParam(tag=("FO"))
    galist <- readGAlignmentsList(bf, param=param)
    checkIdentical(rep.int(NA, length(unlist(galist))), 
                   mcols(unlist(galist))[["FO"]])
}

test_readGAlignmentsList_compare_pairs <- function()
{
    bamfile <- BamFile(untreated3_chr4(), asMates=TRUE)
    galist <- readGAlignmentsList(bamfile)
    mates <- galist[mcols(galist)$mate_status == "mated"]
    galp <- readGAlignmentPairs(bamfile)
    checkIdentical(length(galp), 75409L)
    tbl <- table(mcols(galist))
    checkIdentical(tbl[["mated"]], 75409L)
    checkIdentical(tbl[["ambiguous"]], 0L)
    checkIdentical(tbl[["unmated"]], 21227L)
}

test_readGAlignmentsList_flags <- function()
{
    bamfile <- BamFile(untreated3_chr4(), asMates=TRUE)
    param <- ScanBamParam(flag=scanBamFlag(isProperPair=TRUE))
    galist <- readGAlignmentsList(bamfile, param=param)
    status <- table(mcols(galist)$mate_status)
    checkIdentical(status[["mated"]], 45828L)
    checkIdentical(status[["ambiguous"]], 0L)
    checkIdentical(status[["unmated"]], 0L)
}

## toy_bamfile read summary: 
## --------------------------

## single-end
## s001: 1 primary alignment
## s002: 1 primary alignment + 3 secondary alignments
## s003: unmapped

## paired-end
## p991: 1 pair with a missing mate (can happen if file was subsetted)
## p992: 1 pair with 1st mate unmapped and 2nd mate mapped
## p993: 1 pair with both mates unmapped

## multi-segments
## m001: 3 segments in the template (index of each segment is known)
## m002: 3 segments in the template (index of each segment was lost) 

## mapped pairs ('pi' tag only exists for these mapped pairs)
## p001: 1 primary proper pair
## p002: 1 primary non proper pair
## p003: 2 proper pairs: 1 primary + 1 secondary
## p004: 2 non proper pairs: 1 primary + 1 secondary
## p005: 2 primary pairs
## p006: 3 pairs: 1 primary proper + 1 secondary proper + 
##                1 secondary non proper
## p007: 2 pairs mapped to the same position: 
##       1 primary proper + 1 secondary proper
## p008: 3 pairs mapped to the same position: 1 primary proper + 
##       1 secondary proper + 1 secondary non proper
## p009: 3 pairs mapped to the same position: 
##       1 primary proper + 2 secondary proper. 
source(system.file("unitTests", "test_readGAlignmentPairs.R", 
                   package="GenomicAlignments"))
bf <- BamFile(toy_bamfile, asMates=TRUE)

test_readGAlignmentsList_toybamfile <- function()
{
    param <- ScanBamParam(tag="pi")
    galp <- readGAlignmentPairs(toy_bamfile, use.names=TRUE, param=param)
    galist <- readGAlignmentsList(bf, use.names=TRUE, param=param)

    ## 'mated' 
    mated_galist <- unlist(galist[mcols(galist)$mate_status == "mated"]) 
    pi_target <- c("p001", "p002", "p003a", "p003b", "p004a", "p004b",
                   "p005a", "p005b", "p006a", "p006b", "p006c",
                   "p007a", "p007b", "p008a", "p008b", "p008c", "p009a")
    checkTrue(all(mcols(mated_galist)$pi %in% pi_target))

    ## 'ambiguous' GAList match 'dumped' GAPairs
    ambig_galist <- unlist(galist[mcols(galist)$mate_status == "ambiguous"]) 
    dumped_galp <- getDumpedAlignments()
    pi_target <- rep(c("p009b", "p009c"), each=2)
    checkIdentical(pi_target, sort(mcols(dumped_galp)$pi))
    checkIdentical(pi_target, sort(mcols(ambig_galist)$pi))

    ## 'unmated':
    unmated_galist <- unlist(galist[mcols(galist)$mate_status == "unmated"]) 
    ## unmated single-end, paired-end or multi-segment (no pi tags) 
    name_target <- c("m001", "m002", "p991", "p992", "s001", "s002")
    unmated <- names(unmated_galist)[is.na(mcols(unmated_galist)$pi)] 
    checkTrue(all(unmated %in% name_target)) 
    ## non-proper mapped-pairs (have pi tags) 
    pi_target <- c("p002", "p004a", "p004b", "p006c", "p008c")
    unmated <- na.omit(unique(mcols(unmated_galist)$pi)) 
    checkTrue(all(unmated %in% pi_target)) 

    ## Reads of this type cannot be filtered out wrt readGAlignmentsList.
    ## They are always returned by readGAlignmentsList but never returned by
    ## readGAlignmentPairs.
    bamFlag(param) <- scanBamFlag(isProperPair=TRUE,
                                  hasUnmappedMate=FALSE,
                                  isUnmappedQuery=FALSE,
                                  isPaired=TRUE)
    galist <- readGAlignmentsList(bf, use.names=TRUE, param=param)
    unmated_galist <- unlist(galist[mcols(galist)$mate_status == "unmated"]) 
    unmated <- names(unmated_galist)[is.na(mcols(unmated_galist)$pi)] 
    name_target <- c("m001", "m002", "p991")
    checkTrue(all(unmated %in% name_target))
}

test_readGAlignmentsList_which <- function()
{
    ## 4 non-overlapping regions of interest: first two regions only overlap
    ## with first p001 mate and last two regions only with last p001 mate.
    my_ROI <- GRanges("chr2", IRanges(c(10, 15, 110, 115), width=1))
    my_ROI_labels <- c("chr2:10-10", "chr2:15-15",
                       "chr2:110-110", "chr2:115-115")
    param <- ScanBamParam(tag="pi", which=my_ROI[c(1, 4)])
    target1 <- readGAlignmentsList(toy_bamfile, use.names=TRUE,
                                   param=param, with.which_label=TRUE)
    ## Duplicate results with distinct 'which_label'
    checkIdentical(2L, length(target1))
    checkIdentical(as.vector(mcols(target1)$mate_status), c("mated", "mated"))
    rng1 <- as.vector(mcols(unlist(target1[1]))$which_label)
    checkTrue(all(rng1 %in% my_ROI_labels[1]))
    rng2 <- as.vector(mcols(unlist(target1[2]))$which_label)
    checkTrue(all(rng2 %in% my_ROI_labels[4]))
}
