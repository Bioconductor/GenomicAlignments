\name{map-methods}

\alias{map,GenomicRanges,GAlignments-method}

\title{Mapping ranges between sequences}
\description{
  A method for translating a set of input ranges through a
  \link{GAlignments} object. Returns a \linkS4class{RangesMapping}
  object.

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
  An object of class \code{RangesMapping}. The \pkg{GenomicRanges}
  package provides some additional methods on this object.
  See \code{?`\link[GenomicRanges]{map-methods}`} in the \pkg{GenomicRanges}
  package for more information.
}
\details{
  Each element in \code{to} is taken to represent the
  alignment of a (read) sequence. The CIGAR string is used to
  translate the input ranges to be relative to the read start. This
  is useful, for example, when determining the cycle (read position)
  at which a particular genomic mismatch occurs.
}
\seealso{
  The \linkS4class{RangesMapping} class is the typical return value.
}

\author{M. Lawrence}

