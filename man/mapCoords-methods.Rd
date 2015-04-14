\name{mapCoords-methods}

\alias{mapCoords-methods}
\alias{mapCoords,GenomicRanges,GAlignments-method}

\title{Mapping ranges between sequences}

\description{
  DEPRECATED! Use \code{\link{mapToAlignments}} instead.

  A method for translating a set of input ranges through a
  \link{GAlignments} object. Returns a \linkS4class{GenomicRanges}
  object.

  NOTE: The \code{\link[IRanges]{mapCoords}} generic function is 
  defined and documented in the \pkg{IRanges} package. A method for 
  translating a set of input ranges through a
  \link[GenomicRanges]{GRangesList} object is defined and documented 
  in the \pkg{GenomicRanges} package.
}

\usage{
\S4method{mapCoords}{GenomicRanges,GAlignments}(from, to, ...)
}

\arguments{
  \item{from}{The input ranges to map, usually a \code{\linkS4class{GRanges}}.}
 
  \item{to}{The alignment between the sequences in \code{from} and the
    sequences in the result.}

  \item{\dots}{Arguments passed to other methods.}
}

\value{
  A \code{GRanges} object of mapped coordinates with matching data
  as metadata columns (`fromHits` and `toHits`).

  Matching data are the result of calling \code{findOverlaps} with
  type `within` on  ranges in \code{from} (the query) and the 
  ranges in \code{to} (the subject). Matching can be many-to-one 
  or one-to-many; one row is reported for each match. 
}

\details{
  DEPRECATED! Use \code{\link{mapToAlignments}} instead.

  Each element in \code{to} is taken to represent the
  alignment of a (read) sequence. The CIGAR string is used to
  translate the input ranges to be relative to the read start. This
  is useful, for example, when determining the cycle (read position)
  at which a particular genomic mismatch occurs.
}

\seealso{
  The generic \link[IRanges]{mapCoords-methods} in the IRanges package.
  Additional methods in the GenomicRanges package 
    \link[GenomicRanges]{mapCoords-methods}.
}

\examples{
## DEPRECATED! See ?mapToAlignments for a replacement.
}

\author{M. Lawrence}
