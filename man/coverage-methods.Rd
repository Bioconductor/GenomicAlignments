\name{coverage-methods}
\alias{coverage-methods}

\alias{coverage,BamFile-method}

\title{Coverage of a BamFile object}

\description{
  \code{\link[IRanges]{coverage}} method for \link[Rsamtools]{BamFile}
  objects.

  NOTE: The \code{\link[IRanges]{coverage}} generic function and methods
  for \link[IRanges]{Ranges} and \link[IRanges]{RangesList} objects
  are defined and documented in the \pkg{IRanges} package.
  The methods for \link[GenomicRanges]{GRanges},
  \link[GenomicRanges]{GRangesList}, \link[GenomicRanges]{GAlignments},
  and \link[GenomicRanges]{GAlignmentPairs} objects are defined and
  documented in the \pkg{GenomicRanges} package.
}

\usage{
\S4method{coverage}{BamFile}(x, shift=0L, width=NULL, weight=1L,
            ..., param=ScanBamParam())
}

\arguments{
  \item{x}{
    A \link[Rsamtools]{BamFile} object.
  }
  \item{shift, width, weight}{
    See \code{coverage} method for GAlignments objects.
  }
  \item{...}{
    Additional arguments passed to the \code{coverage} method for
    \link[GenomicRanges]{GAlignments} objects.
  }
  \item{param}{
    An optional \link[Rsamtools]{ScanBamParam} object passed to
    \code{\link{readGAlignmentsFromBam}}.
  }
}

\details{
  The method for \link[Rsamtools]{BamFile} objects does:
\preformatted{
  coverage(readGAlignmentsFromBam(x, param=param),
           shift=shift, width=width, weight=weight, ...)
}
}

\value{
  A named \link[IRanges]{RleList} object with one coverage vector per
  seqlevel in \code{x}.
}

\seealso{
  \itemize{
    \item \code{\link[IRanges]{coverage}} in the \pkg{IRanges} package
          and \link[GenomicRanges]{coverage-methods} in the
          \pkg{GenomicRanges} package.

    \item \code{\link{readGAlignmentsFromBam}}.
  }
}

\examples{
##TODO
}

\keyword{methods}
\keyword{utilities}
