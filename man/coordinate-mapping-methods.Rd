\name{mapToAlignment}
\alias{mapToAlignment}
\alias{mapToAlignment,Ranges,GAlignments-method}
\alias{mapToAlignment,GenomicRanges,GAlignments-method}
\alias{pmapToAlignment}
\alias{pmapToAlignment,Ranges,GAlignments-method}
\alias{pmapToAlignment,GenomicRanges,GAlignments-method}
\alias{coordinate-mapping-methods}

\title{Map range coordinates between local and genome space using alignments}

\description{
  Map range coordinates between local and genome (reference) space using the
  CIGAR in a \code{GAlignments} object.
 
  See \code{?`\link[GenomicFeatures]{mapToTranscriptome-methods}`} in the 
  \pkg{GenomicRanges} package for mapping ranges between transcriptome
  and genome space.
}

\usage{
\S4method{mapToAlignment}{Ranges,GAlignments}(x, alignment, 
         reverse = FALSE, ...)
\S4method{mapToAlignment}{GenomicRanges,GAlignments}(x, alignment, 
         reverse = FALSE, ...)

\S4method{pmapToAlignment}{Ranges,GAlignments}(x, alignment,
         reverse = FALSE, ...)
\S4method{pmapToAlignment}{GenomicRanges,GAlignments}(x, alignment, 
         reverse = FALSE, ...)
}

\arguments{
  \item{x}{
    \code{\linkS4class{GenomicRanges}} object of positions to be mapped.
  }
  \item{alignment}{
    A \code{\linkS4class{GAlignments}} object that represents the alignment 
    (aka mapping) between positions in \code{x} and the result.
  }
  \item{reverse}{
    TRUE or FALSE. When FALSE (default) the range in \code{x} is mapped to the
    reference space using the CIGAR in \code{alignment}. TRUE indicates the
    reverse mapping; the CIGAR in \code{alignment} is used to map the range in
    \code{x} to the local position. See examples.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  \itemize{
    \item{\code{mapToAlignment}}{
      These methods use a \code{GAlignments} object to represent the alignment 
      between the ranges in \code{x} and the output. The following CIGAR 
      operations in the "Extended CIGAR format" are used in the mapping
      algorithm:
      \itemize{
        \item{M, X, =} Sequence match or mismatch
        \item{I} Insertion to the reference 
        \item{D} Deletion from the reference 
        \item{N} Skipped region from the reference
        \item{S} Soft clip on the read
        \item{H} Hard clip on the read
        \item{P} Silent deletion from the padded reference
      } 
      When \code{reverse=FALSE} (default) \code{x} is considered a `local`
      position (e.g., nucleotide position of interest where multiple reads had
      poor alignments). The return value is the reference position corresponding
      to \code{x}. The mapped position starts at \code{start(x) +
      start(alignment) + 1} and is incremented or decremented as the algorithm
      walks the length of the CIGAR.

      When \code{reverse=TRUE}, \code{x} is considered to be a reference
      position. In this case the CIGAR is used to map \code{x} to the `local`
      position. The mapped position starts at \code{x - start(alignment) + 1}
      and is incremented or decremented as the algorithm walks the length of the
      CIGAR.
    }
    \item{\code{pmapToAlignment}}{
      Element-wise (aka `parallel`) versions of \code{mapToAlignment}. The 
      i-th range in \code{x} is mapped to the i-th range in \code{alignment}; 
      \code{x} and \code{alignment} must have the same length. The
      \code{reverse} argument has the same meaning as for \code{mapToAlignment}.
    }
    \item{strand}{
      By SAM convention, the CIGAR string is reported for mapped reads on the 
      forward genomic strand. There is no need to consider strand in these
      methods. When the output is a \code{GenomicRanges} object the strand
      will always be "*".
    }
  }
}

\value{
  An object the same class as \code{x}.

  Parallel methods return an object the same shape as \code{x}. Ranges that
  are non-hits (out of bounds or strand-mismatch) are returned as zero-width
  ranges starting at 1. The seqname of the non-hits is "unmapped".

  Non-parallel methods return an object that varies in length similar to a
  Hits object. The result only contains mapped records, strand mismatch
  and non-hits are not returned. \code{xHits} and \code{alignmentHits}
  metadata columns record the elements of \code{x} and \code{alignment} used
  in the mapping.
}

\seealso{
  \itemize{
    %\item \code{?`\link[GenomicFeatures]{mapToTranscriptome-methods}`} in the 
    %  in the \pkg{GenomicFeatures} package for methods mapping to and from
    %  the transcriptome.

    \item \url{http://samtools.sourceforge.net/} for a description of
      the Extended CIGAR format.
  }
}

\author{V. Obenchain, M. Lawrence and H. Pages}

\examples{
## ---------------------------------------------------------------------
## A. Basic use 
## ---------------------------------------------------------------------

## 1. Map to genome space (reverse = FALSE)
## --------------------------------------------------------
## The simplest example is an alignment with no junctions, insertions 
## or deletions in the CIGAR.
x <- GRanges("chr1", IRanges(c(5, 20), width = 2))
alignment <- GAlignments("chr1", 10L, "11M", strand("*"))

## Mapping to genome space returns a result for the range starting at 
## 5 (10 + 5 - 1 = 14). The range starting at 20 is not within the
## alignment range (20 + 10 - 1 = 29).
mapToAlignment(x, alignment, reverse = FALSE)

## Mapping the same range starting at 5 to an alignment with an insertion
## to the reference ('I' CIGAR operation) shortens the result by 1.
x <- GRanges("chr1", IRanges(5, width = 2))
alignment <- GAlignments("chr1", 10L, "5M10I4M", strand("+"))
mapToAlignment(x, alignment, reverse = FALSE)

## A skipped region from the reference ('N' CIGAR operation) extends the 
## length of the result by the number of skips.
x <- GRanges("chr1", IRanges(5, width = 2))
alignment <- GAlignments("chr1", 10L, "5M10N4M", strand("+"))
mapToAlignment(x, alignment, reverse = FALSE)


## 2. Map to local space (reverse = TRUE)
## ---------------------------------------------------------------------

## For a successful mapping in the opposite direction, 'x' must be 
## located within the 'alignment' range. In this 'x', the start 
## of both ranges fall within 'alignment' but the width of the 
## second makes it too large.
x <- GRanges("chr1", IRanges(c(12, 12), width = c(6, 20)))
alignment <- GAlignments("chr1", 10L, "11M", strand("*"))
mapToAlignment(x, alignment, reverse = TRUE)

## The element-wise version of the function maps the i-th element of 'x'
## to the i-th element of 'alignment'. Here we can see the affect of
## mapping the same range with different CIGAR ops.
x <- GRanges(rep("chr1", 4), IRanges(rep(12, 4), width = rep(6, 4)))
cigar <- c("11M", "5M2I4M", "5M2N4M", "5M2D4M") 
strand <- strand(rep("*", 4))
alignment <- GAlignments(rep("chr1", 4), rep(10L, 4) , cigar, strand)
map <- pmapToAlignment(x, alignment, reverse = TRUE)
names(map) <- c("reference", "insertion", "junction", "deletion")
map
}

\keyword{methods}
\keyword{utilities}
