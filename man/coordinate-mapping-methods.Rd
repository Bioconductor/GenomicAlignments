\name{mapToAlignments}
\alias{coordinate-mapping-methods}

\alias{mapToAlignments}
\alias{mapToAlignments,Ranges,GAlignments-method}
\alias{mapToAlignments,GenomicRanges,GAlignments-method}
\alias{pmapToAlignments}
\alias{pmapToAlignments,Ranges,GAlignments-method}
\alias{pmapToAlignments,GenomicRanges,GAlignments-method}

\alias{mapFromAlignments}
\alias{mapFromAlignments,Ranges,GAlignments-method}
\alias{mapFromAlignments,GenomicRanges,GAlignments-method}
\alias{pmapFromAlignments}
\alias{pmapFromAlignments,Ranges,GAlignments-method}
\alias{pmapFromAlignments,GenomicRanges,GAlignments-method}

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
\S4method{mapToAlignments}{GenomicRanges,GAlignments}(x, alignments, ...) 
\S4method{pmapToAlignments}{GenomicRanges,GAlignments}(x, alignments, ...) 

\S4method{mapFromAlignments}{GenomicRanges,GAlignments}(x, alignments, ...) 
\S4method{pmapFromAlignments}{GenomicRanges,GAlignments}(x, alignments, ...) 
}

\arguments{
  \item{x}{
    \code{\linkS4class{GenomicRanges}} object of positions to be mapped.
  }
  \item{alignments}{
    A \code{\linkS4class{GAlignments}} object that represents the alignment 
    of \code{x} to the genome.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
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
    \item{\code{mapToAlignments}, \code{pmapToAlignments}}{
      The CIGAR is used to map the genomic (reference) position \code{x} to 
      local coordinates. The mapped position starts at 
      \preformatted{
      start(x) - start(alignments) + 1
      }
      and is incremented or decremented as the algorithm walks the length of 
      the CIGAR. A successful mapping in this direction requires that
      \code{x} fall within \code{alignments}.

       The seqlevels of the return object are taken from the 
       \code{alignments} object and will be a name descriptive of the read
       or aligned region. In this direction, mapping is attempted between all 
       elements of \code{x} and all elements of \code{alignments}.
    }
    \item{\code{mapFromAlignments}, \code{pmapFromAlignments}}{
      The CIGAR is used to map the local position \code{x} to genomic
      (reference) coordinates. The mapped position starts at
      \preformatted{
      start(x) + start(alignments) - 1
      } 
      and is incremented or decremented as the algorithm walks the length of 
      the CIGAR. A successful mapping in this direction requires that the
      width of \code{alignments} is <= the width of \code{x}.

      When mapping to the genome, name matching is used to determine the
      mapping pairs (vs attempting to match all possible pairs). Ranges in
      \code{x} are only mapped to ranges in \code{alignments} with the
      same name. Name matching is motivated by use cases such as
      differentially expressed regions where the expressed regions in
      \code{x} would only be related to a subset of regions in
      \code{alignments}, which may contains gene or transcript ranges.
    }
    \item{element-wise versions}{
      \code{pmapToAlignments} and \code{pmapFromAlignments} are element-wise
      (aka `parallel`) versions of \code{mapToAlignments} and 
      \code{mapFromAlignments}. The i-th range in \code{x} is mapped to the 
      i-th range in \code{alignments}; \code{x} and \code{alignments} must 
      have the same length.

      Ranges in \code{x} that do not map (out of bounds or strand mismatch) 
      are returned as zero-width ranges starting at 1. These ranges are given
      the special seqname of "unmapped". Note the non-parallel methods do not 
      return unmapped ranges so the "unmapped" seqname is unique to 
      \code{pmapToAlignments} and \code{pmapFromAlignments}.
    }
    \item{strand}{
      By SAM convention, the CIGAR string is reported for mapped reads on the 
      forward genomic strand. There is no need to consider strand in these
      methods. The output of these methods will always be unstranded 
      (i.e., "*").
    }
  }
}

\value{
  An object the same class as \code{x}.

  Parallel methods return an object the same shape as \code{x}. Ranges that
  cannot be mapped (out of bounds) are returned as zero-width ranges starting 
  at 1 with a seqname of "unmapped".

  Non-parallel methods return an object that varies in length similar to a
  Hits object. The result only contains mapped records, out of bound ranges 
  are not returned. \code{xHits} and \code{alignmentsHits} metadata columns 
  indicate the elements of \code{x} and \code{alignments} used in the mapping.

  When present, names from \code{x} are propagated to the output. When
  mapping locally, the seqlevels of the output are the names of 
  \code{alignment}; when mapping globally the seqlevels are the seqlevels 
  of \code{alignment} (chromosome names).
}

\seealso{
  \itemize{
    \item \code{?\link[GenomicFeatures]{mapToTranscripts}} in the 
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

## 1. Map to local space with mapToAlignments()
## ---------------------------------------------------------------------

## Mapping to local coordinates requires 'x' to be within the
## 'alignments' range. In this example only the first 'x' range
## can be mapped because the width of the second is too large.
alignments <- GAlignments("chr1", 10L, "11M", strand("*"), names="read_A")
x <- GRanges("chr1", IRanges(c(12, 12), width=c(6, 20)))
mapToAlignments(x, alignments)

## The element-wise version of the function returns the second range
## with a seqlevel of "unmapped":
pmapToAlignments(x, c(alignments, alignments))

## Mapping the same range through different alignments demonstrates 
## how the CIGAR operations affect the outcome.
ops <- c("no-op", "junction", "insertion", "deletion")
x <- GRanges(rep("chr1", 4), IRanges(rep(12, 4), width=rep(6, 4), names=ops)) 
alignments <- GAlignments(rep("chr1", 4), rep(10L, 4), 
                         cigar = c("11M", "5M2N4M", "5M2I4M", "5M2D4M"),
                         strand = strand(rep("*", 4)),
                         names = paste0("region_", 1:4))
pmapToAlignments(x, alignments)

## 2. Map to genome space with mapFromAlignments()
## ---------------------------------------------------------------------

## Mapping to genomic coordinates requires the shifted 'x' to fall
## within the 'alignments' range. The first 'x' range has a shifted
## start value of 14 (5 + 10 - 1 = 14) with a width of 2 and is
## successfully mapped. The second has a shifted start of 29
## (20 + 10 - 1 = 29) which is outside the range of 'alignments'.
x <- GRanges("chr1", IRanges(c(5, 20), width = 2, 
             names = rep("region_A", 2)))
alignments <- GAlignments("chr1", 10L, "11M", strand("*"), names = "region_A")
mapFromAlignments(x, alignments)

## When mapping to the genome, a map is only attempted between ranges with 
## the same name. Changing the name of the first 'x' range makes it
## incompatible with 'alignments' so we get an empty result.
names(x) <- c("region_B", "region_A")
mapFromAlignments(x, alignments)

## In this direction a CIGAR insertion reduces the width of the result 
## while junctions and deletions increase it.
ops <- c("no-op", "junction", "insertion", "deletion")
x <- GRanges(rep("chr1", 4), IRanges(rep(3, 4), width=rep(5, 4), names=ops)) 
alignments <- GAlignments(rep("chr1", 4), rep(10L, 4), 
                         cigar = c("11M", "5M2N4M", "5M2I4M", "5M2D4M"),
                         strand = strand(rep("*", 4)))
pmapFromAlignments(x, alignments)

## ---------------------------------------------------------------------
## B. Round-trip from read -> genome -> transcript 
## ---------------------------------------------------------------------

## TODO

}

\keyword{methods}
\keyword{utilities}
