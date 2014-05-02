\name{findSpliceOverlaps-methods}

\alias{findSpliceOverlaps-methods}

\alias{findSpliceOverlaps}
\alias{findSpliceOverlaps,GRangesList,GRangesList-method}
\alias{findSpliceOverlaps,GAlignments,GRangesList-method}
\alias{findSpliceOverlaps,GAlignmentPairs,GRangesList-method}
\alias{findSpliceOverlaps,character,ANY-method}
\alias{findSpliceOverlaps,BamFile,ANY-method}


\title{Classify ranges (reads) as compatible with existing genomic annotations
       or as having novel splice events}

\description{
  The \code{findSpliceOverlaps} function identifies ranges (reads) that are
  compatible with a specific transcript isoform. The non-compatible ranges are
  analyzed for the presence of novel splice events.
}

\usage{
findSpliceOverlaps(query, subject, ignore.strand=FALSE, ...)

\S4method{findSpliceOverlaps}{GRangesList,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)

\S4method{findSpliceOverlaps}{GAlignments,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)

\S4method{findSpliceOverlaps}{GAlignmentPairs,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)

\S4method{findSpliceOverlaps}{BamFile,ANY}(query, subject, ignore.strand=FALSE, ...,
                   param=ScanBamParam(), singleEnd=TRUE)
}

\arguments{
  \item{query}{
    A \link[GenomicRanges]{GRangesList}, \link{GAlignments},
    \link{GAlignmentPairs}, or \link[Rsamtools]{BamFile} object
    containing the reads.
    Can also be a single string containing the path to a BAM file.

    Single or paired-end reads are specified with the \code{singleEnd} 
    argument (default FALSE). Paired-end reads can be supplied in a BAM
    file or \link{GAlignmentPairs} object. Single-end are
    expected to be in a BAM file, \link{GAlignments} or 
    \link[GenomicRanges]{GRanges} object.
  }
  \item{subject}{
    A \link[GenomicRanges]{GRangesList} containing the annotations.
    This list is expected to contain exons grouped by transcripts. 
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
    Optional \link[GenomicRanges]{GRangesList} of coding regions for
    each transcript in the \code{subject}. If provided, the "coding"
    output column will be a \code{logical} vector indicating if the
    read falls in a coding region. When not provided, the "coding"
    output is \code{NA}.
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
  When a read maps compatibly and uniquely to a transcript isoform we
  can quantify the expression and look for shifts in the balance of
  isoform expression. If a read does not map in compatible way, novel
  splice events such as splice junctions, novel exons or retentions
  can be quantified and compared across samples.

  \code{findSpliceOverlaps} detects which reads (query) match to
  transcripts (subject) in a compatible fashion. Compatibility is based
  on both the transcript bounds and splicing pattern. Assessing the
  splicing pattern involves comparision of the read splices (i.e., the
  N operations in the CIGAR) with the transcript introns. For paired-end
  reads, the inter-read gap is not considered a splice junction. The analysis
  of non-compatible reads for novel splice events is under construction.
}

\value{
  The output is a \link[S4Vectors]{Hits} object with
  the metadata columns defined below. Each column is a \code{logical}
  indicating if the read (query) met the criteria.

  \itemize{
    \item{compatible: }{Every splice (N) in a read alignment matches
      an intron in an annotated transcript. The read does not
      extend into an intron or outside the transcript bounds.
    }
    \item{unique: }{The read is compatible with only one annotated
      transcript.
    }
    \item{strandSpecific: }{The query (read) was stranded.
    }
  }
}

\author{
  Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org> 
}
 
\seealso{
  \itemize{
    \item \link[GenomicRanges]{GRangesList} objects in the
          \pkg{GenomicRanges} package.

    \item \link{GAlignments} and \link{GAlignmentPairs} objects.

    \item \link[Rsamtools]{BamFile} objects in the \pkg{Rsamtools} package.
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
