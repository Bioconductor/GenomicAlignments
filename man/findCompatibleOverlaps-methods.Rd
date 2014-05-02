\name{findCompatibleOverlaps-methods}
\alias{findCompatibleOverlaps-methods}

\alias{findCompatibleOverlaps}
\alias{findCompatibleOverlaps,GAlignments,GRangesList-method}
\alias{findCompatibleOverlaps,GAlignmentPairs,GRangesList-method}
\alias{countCompatibleOverlaps}


\title{Finding hits between reads and transcripts that are \emph{compatible}
       with the splicing of the transcript}

\description{
  In the context of an RNA-seq experiment, \code{findCompatibleOverlaps}
  (or \code{countCompatibleOverlaps}) can be used for finding (or counting)
  hits between reads and transcripts that are \emph{compatible}
  with the splicing of the transcript.
}

\usage{
findCompatibleOverlaps(query, subject)
countCompatibleOverlaps(query, subject)
}

\arguments{
  \item{query}{
    A \link{GAlignments} or \link{GAlignmentPairs} object representing
    the aligned reads.
  }
  \item{subject}{
    A \link{GRangesList} object representing the transcripts.
  }
}

\details{
  \code{findCompatibleOverlaps} is a specialized version of
  \code{\link[IRanges]{findOverlaps}} that uses
  \code{\link{encodeOverlaps}} internally to keep only
  the hits where the junctions in the aligned read are \emph{compatible}
  with the splicing of the annotated transcript.

  The topic of working with overlap encodings is covered in details
  in the "OverlapEncodings" vignette located this package
  (\pkg{GenomicAlignments}) and accessible with
  \code{vignette("OverlapEncodings")}.
}

\value{
  A \link[S4Vectors]{Hits} object for \code{findCompatibleOverlaps}.

  An integer vector \emph{parallel} to (i.e. same length as) \code{query}.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \code{\link[IRanges]{findOverlaps}} generic function defined
          in the \pkg{IRanges} package.

    \item The \code{\link{encodeOverlaps}} generic function and
          \link{OverlapEncodings} class.

    \item The "OverlapEncodings" vignette in this package.

    \item \link{GAlignments} and \link{GAlignmentPairs} objects.

    \item \link[GenomicRanges]{GRangesList} objects in the
          \pkg{GenomicRanges} package.
  }
}

\examples{
## Here we only show a simple example illustrating the use of
## countCompatibleOverlaps() on a very small data set. Please
## refer to the "OverlapEncodings" vignette in the GenomicAlignments
## package for a comprehensive presentation of "overlap
## encodings" and related tools/concepts (e.g. "compatible"
## overlaps, "almost compatible" overlaps etc...), and for more
## examples.

## sm_treated1.bam contains a small subset of treated1.bam, a BAM
## file containing single-end reads from the "Pasilla" experiment
## (RNA-seq, Fly, see the pasilla data package for the details)
## and aligned to reference genome BDGP Release 5 (aka dm3 genome on
## the UCSC Genome Browser):
sm_treated1 <- system.file("extdata", "sm_treated1.bam",
                           package="GenomicAlignments", mustWork=TRUE)

## Load the alignments:
flag0 <- scanBamFlag(isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
param0 <- ScanBamParam(flag=flag0)
gal <- readGAlignments(sm_treated1, use.names=TRUE, param=param0)

## Load the transcripts (IMPORTANT: Like always, the reference genome
## of the transcripts must be *exactly* the same as the reference
## genome used to align the reads):
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
exbytx <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, by="tx", use.names=TRUE)

## Number of "compatible" transcripts per alignment in 'gal':
gal_ncomptx <- countCompatibleOverlaps(gal, exbytx)
mcols(gal)$ncomptx <- gal_ncomptx
table(gal_ncomptx)
mean(gal_ncomptx >= 1)
## --> 33% of the alignments in 'gal' are "compatible" with at least
## 1 transcript in 'exbytx'.

## Keep only alignments compatible with at least 1 transcript in
## 'exbytx':
compgal <- gal[gal_ncomptx >= 1]
head(compgal)
}

\keyword{methods}
\keyword{utilities}
