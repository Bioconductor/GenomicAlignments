\name{mapToAlignments}
\alias{coordinate-mapping-methods}

\alias{mapToAlignments}
\alias{mapToAlignments,Ranges,GAlignments-method}
\alias{mapToAlignments,GenomicRanges,GAlignments-method}

\alias{pmapToAlignments}
\alias{pmapToAlignments,Ranges,GAlignments-method}
\alias{pmapToAlignments,GenomicRanges,GAlignments-method}

\title{Map range coordinates between reads and genome space using
       CIGAR alignments}

\description{
  Map range coordinates between reads (local) and genome (reference) space 
  using the CIGAR in a \code{GAlignments} object.
 
  See \code{?\link[GenomicFeatures]{mapToTranscripts}} in the 
  \pkg{GenomicRanges} package for mapping coordinates between features
  in the transcriptome and genome space.
}

\usage{
\S4method{mapToAlignments}{Ranges,GAlignments}(x, alignments, 
         reverse = FALSE, ...)
\S4method{mapToAlignments}{GenomicRanges,GAlignments}(x, alignments, 
         reverse = FALSE, ...)

\S4method{pmapToAlignments}{Ranges,GAlignments}(x, alignments,
         reverse = FALSE, ...)
\S4method{pmapToAlignments}{GenomicRanges,GAlignments}(x, alignments, 
         reverse = FALSE, ...)
}

\arguments{
  \item{x}{
    \code{\linkS4class{GenomicRanges}} object of positions to be mapped.
  }
  \item{alignments}{
    A \code{\linkS4class{GAlignments}} object that represents the alignment 
    of \code{x} to the genome.
  }
  \item{reverse}{
    TRUE or FALSE. When FALSE (default) the range in \code{x} is mapped to the
    local space, when TRUE \code{x} is mapped to genomic space. Both mappings
    use the CIGAR in \code{alignments}. See examples.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  \itemize{
    \item{\code{mapToAlignments}}{
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
      \itemize{
        \item{\code{reverse = FALSE} (default)}{
          The CIGAR is used to map the genomic (reference) position \code{x} to 
          local coordinates. The mapped position starts at 
          \preformatted{
          start(x) - start(alignments) + 1
          }
          and is incremented or decremented as the algorithm walks the length of 
          the CIGAR. A successful mapping in this direction requires that
          \code{x} fall within \code{alignments}.
        }
        \item{\code{reverse = TRUE}}{
          The CIGAR is used to map the local position \code{x} to genomic
          (reference) coordinates. The mapped position starts at
          \preformatted{
          start(x) + start(alignments) + 1
          } 
          and is incremented or decremented as the algorithm walks the length of 
          the CIGAR. A successful mapping in this direction requires that the
          width of \code{alignments} is <= the width of \code{x}.
        }
      }
    }
    \item{\code{pmapToAlignments}}{
      Element-wise (aka `parallel`) version of \code{mapToAlignments}. The 
      i-th range in \code{x} is mapped to the i-th range in \code{alignments}; 
      \code{x} and \code{alignments} must have the same length. The
      \code{reverse} argument has the same meaning as described in
      code{mapToAlignments}.

      Ranges in \code{x} that do not map (out of bounds or strand mismatch) 
      are returned as zero-width ranges starting at 1. These ranges are given
      the special seqname of "unmapped". Note the non-parallel methods do not 
      return unmapped ranges so the "unmapped" seqname is unique to 
      \code{pmapToAlignments}.
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
  and non-hits are not returned. \code{xHits} and \code{alignmentsHits}
  metadata columns record the elements of \code{x} and \code{alignments} used
  in the mapping.
}

\seealso{
  \itemize{
    \item \code{?\link[GenomicFeatures]{mapToTranscriptome}} in the 
      in the \pkg{GenomicFeatures} package for methods mapping between
      transcriptome and genome space.

    \item \url{http://samtools.sourceforge.net/} for a description of
      the Extended CIGAR format.
  }
}

\author{V. Obenchain, M. Lawrence and H. Pages}

\examples{
## ---------------------------------------------------------------------
## A. Basic use 
## ---------------------------------------------------------------------

## 1. Map to local space (reverse = FALSE)
## ---------------------------------------------------------------------

## For a successful map to local coordinates 'x' must be within 'alignments'.
alignments <- GAlignments("chr1", 10L, "11M", strand("*"))
x <- GRanges("chr1", IRanges(c(12, 12), width=c(6, 20)))

## The width of the second range in 'x' is too large and is a 'non-hit'.
mapToAlignments(x, alignments)

## Non-hits are returned as "unmapped" in the element-wise version:
pmapToAlignments(x, c(alignments, alignments))

## Here we can see the affect of mapping the same range through different 
## alignments (aka CIGAR operations):
ops <- c("no-op", "junction", "insertion", "deletion")
x <- GRanges(rep("chr1", 4), IRanges(rep(12, 4), width=rep(6, 4), names=ops)) 
alignments <- GAlignments(rep("chr1", 4), rep(10L, 4), 
                         cigar = c("11M", "5M2N4M", "5M2I4M", "5M2D4M"),
                         strand = strand(rep("*", 4)))
pmapToAlignments(x, alignments)

## 2. Map to genome space (reverse = TRUE)
## --------------------------------------------------------

## Mapping to genome space returns a result for the range in 'x'
## that starts at position 5 (10 + 5 - 1 = 14). The range starting 
## at 20 does not fall within the alignments range (20 + 10 - 1 = 29).
x <- GRanges("chr1", IRanges(c(5, 20), width = 2))
alignments <- GAlignments("chr1", 10L, "11M", strand("*"))
mapToAlignments(x, alignments, reverse = TRUE)

## Mapping the same range starting at 5 to an alignment with an insertion
## to the reference ('I' CIGAR operation) shortens the result by 1.
x <- GRanges("chr1", IRanges(5, width = 2))
alignments <- GAlignments("chr1", 10L, "5M10I4M", strand("+"))
mapToAlignments(x, alignments, reverse = TRUE)

## A skipped region from the reference ('N' CIGAR operation) extends the 
## length of the result by the number of skips.
x <- GRanges("chr1", IRanges(5, width = 2))
alignments <- GAlignments("chr1", 10L, "5M10N4M", strand("+"))
mapToAlignments(x, alignments, reverse = TRUE)
}

\keyword{methods}
\keyword{utilities}
