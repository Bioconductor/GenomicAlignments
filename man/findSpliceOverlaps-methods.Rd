\name{findSpliceOverlaps-methods}

\alias{findSpliceOverlaps-methods}

\alias{findSpliceOverlaps,GAlignments,GRangesList-method}
\alias{findSpliceOverlaps,GAlignmentPairs,GRangesList-method}
\alias{findSpliceOverlaps,character,ANY-method}
\alias{findSpliceOverlaps,BamFile,ANY-method}


\title{findSpliceOverlaps methods for GAlignments and GAlignmentPairs objects}

\description{
  The \code{findSpliceOverlaps} function identifies ranges (reads) that are 
  compatible with a specific transcript isoform. The non-compatible ranges are
  analyzed for the presence of novel splice events.

  NOTE: The \code{\link[GenomicRanges]{findSpliceOverlaps}} generic function
  and method for \link[GenomicRanges]{GRangesList} objects are defined and
  documented in the \pkg{GenomicRanges} package.
  Methods for \link[GenomicFeatures]{TranscriptDb} objects are defined and
  documented in the \pkg{GenomicFeatures} package.
}

\usage{
\S4method{findSpliceOverlaps}{GAlignments,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)

\S4method{findSpliceOverlaps}{GAlignmentPairs,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)

\S4method{findSpliceOverlaps}{BamFile,ANY}(query, subject, ignore.strand=FALSE, ...,
    param=ScanBamParam(), singleEnd=TRUE)
}

\arguments{
  \item{query}{
    A \link[GenomicRanges]{GAlignments}, \link[GenomicRanges]{GAlignmentPairs},
    or \link[Rsamtools]{BamFile} object containing the reads.
    Can also be a single string containing the path to a BAM file.

    Single or paired-end reads are specified with the \code{singleEnd} 
    argument (default FALSE). Paired-end reads can be supplied in a BAM
    file or \linkS4class{GAlignmentPairs} object. Single-end are
    expected to be in a BAM file, \linkS4class{GAlignments} or 
    \linkS4class{GRanges} object.
  }
  \item{subject}{
    A \link{GRangesList} containing the annotations. This list is
    expected to be exons by transcripts. 
  }
  \item{ignore.strand}{
    When set to \code{TRUE}, strand information is ignored in the
    overlap calculations.
  }
  \item{...}{
    Additional arguments such as \code{param} and \code{singleEnd} used
    in the method for \link[Rsamtools]{BamFile} objects. See below.
  }
  \item{cds}{
    See \code{\link[GenomicRanges]{findSpliceOverlaps}} in the
    \pkg{GenomicRanges} package.
  }
  \item{param}{
    An optional \code{\link[Rsamtools]{ScanBamParam}} instance to
    further influence scanning, counting, or filtering.
  }
  \item{singleEnd}{
    A logical value indicating if reads are single or paired-end.
    See \code{\link{summarizeOverlaps}} for more information.
  }
}

\details{
  See \code{\link[GenomicRanges]{findSpliceOverlaps}} in the
  \pkg{GenomicRanges} package for the details.
}

\value{
  See \code{\link[GenomicRanges]{findSpliceOverlaps}} in the
  \pkg{GenomicRanges} package for the returned value.
}

\author{
  Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org> 
}
 
\seealso{
  \itemize{
    \item \link[GenomicRanges]{findSpliceOverlaps} in the
          \pkg{GenomicRanges} package.
    \item The \link[GenomicRanges]{GRangesList},
          \link[GenomicRanges]{GAlignments}, and
          \link[GenomicRanges]{GAlignmentPairs} classes.
  }
}

\examples{
## ----------------------------------------------------------------------- 
## Isoform expression :
## ----------------------------------------------------------------------- 
## findSpliceOverlaps() can assist in quantifying isoform expression
## by identifying reads that map compatibly and uniquely to a 
## transcript isoform.
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(pasillaBamSubset)
se <- untreated1_chr4() ## single-end reads
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
exbytx <- exonsBy(txdb, "tx")
cdsbytx <- cdsBy(txdb, "tx")
param <- ScanBamParam(which=GRanges("chr4", IRanges(1e5,3e5)))
sehits <- findSpliceOverlaps(se, exbytx, cds=cdsbytx, param=param)

## Tally the reads by category to get an idea of read distribution.
lst <- lapply(mcols(sehits), table)
nms <- names(lst) %in% c("compatible", "unique")
tbl <- do.call(rbind, lst[nms]) 
tbl

## Reads compatible with one or more transcript isoforms.
rnms <- rownames(tbl)
tbl[rnms == "compatible","TRUE"]/sum(tbl[rnms == "compatible",])

## Reads compatible with a single isoform. 
tbl[rnms == "unique","TRUE"]/sum(tbl[rnms == "unique",])

## All reads fall in a coding region as defined by 
## the txdb annotation.
lst[["coding"]] 

## Check : Total number of reads should be the same across categories.
lapply(lst, sum)
 
## ----------------------------------------------------------------------- 
## Paired-end reads :
## ----------------------------------------------------------------------- 
## 'singleEnd' is set to FALSE for a BAM file with paired-end reads.
pe <- untreated3_chr4()
hits2 <- findSpliceOverlaps(pe, exbytx, singleEnd=FALSE, param=param)
 
## In addition to BAM files, paired-end reads can be supplied in a 
## GAlignmentPairs object. 
genes <- GRangesList(
    GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
    GRanges("chr1", IRanges(c(5, 22), c(15, 25)), "+"))
galp <- GAlignmentPairs(
    GAlignments("chr1", 5L, "11M4N6M", strand("+")),
    GAlignments("chr1", 50L, "6M", strand("-")),
    isProperPair=TRUE)
findSpliceOverlaps(galp, genes)
}

\keyword{methods}
\keyword{utilities}
