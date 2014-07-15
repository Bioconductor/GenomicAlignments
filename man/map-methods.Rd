\name{map-methods}

\alias{map-methods}
\alias{map,GenomicRanges,GAlignments-method}

\title{Mapping ranges between sequences}
\description{
  \code{map} and \code{pmap} are defunct. Use \code{mapCoords}
  and \code{pmapCoords} instead.

  A method for translating a set of input ranges through a
  \link{GAlignments} object. Returns a GRangesMapping object.

  NOTE: The \code{\link[IRanges]{map}} generic function is defined and
  documented in the \pkg{IRanges} package.
  A method for translating a set of input ranges through a
  \link[GenomicRanges]{GRangesList} object is defined and documented
  in the \pkg{GenomicRanges} package.
}
\usage{
\S4method{map}{GenomicRanges,GAlignments}(from, to)
}
\arguments{
  \item{from}{The input ranges to map, usually
    a \code{\linkS4class{GenomicRanges}}}
  \item{to}{The alignment between the sequences in \code{from} and the
    sequences in the result.}
}
\value{
  An object of class GRangesMapping, which consists of a
  \code{Hits} instance and \code{GRanges} instance. 
}
\details{
  Each element in \code{to} is taken to represent the
  alignment of a (read) sequence. The CIGAR string is used to
  translate the input ranges to be relative to the read start. This
  is useful, for example, when determining the cycle (read position)
  at which a particular genomic mismatch occurs.
}
\seealso{
  \code{\link{mapCoords}} and \code{\link{pmapCoords}} methods.
}

\author{M. Lawrence}

