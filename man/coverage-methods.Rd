\name{coverage-methods}
\alias{coverage-methods}

\alias{coverage,GAlignments-method}
\alias{coverage,GAlignmentPairs-method}
\alias{coverage,BamFile-method}


\title{Coverage of a GAlignments or GAlignmentPairs object}

\description{
  \code{\link[IRanges]{coverage}} methods for \link{GAlignments},
  \link{GAlignmentPairs}, and \link[Rsamtools]{BamFile} objects.

  NOTE: The \code{\link[IRanges]{coverage}} generic function and methods
  for \link[IRanges]{Ranges} and \link[IRanges]{RangesList} objects
  are defined and documented in the \pkg{IRanges} package.
  Methods for \link[GenomicRanges]{GRanges} and
  \link[GenomicRanges]{GRangesList} objects are defined and
  documented in the \pkg{GenomicRanges} package.
}

\usage{
\S4method{coverage}{GAlignments}(x, shift=0L, width=NULL, weight=1L,
         method=c("auto", "sort", "hash"), drop.D.ranges=FALSE)

\S4method{coverage}{GAlignmentPairs}(x, shift=0L, width=NULL, weight=1L,
         method=c("auto", "sort", "hash"), drop.D.ranges=FALSE)

\S4method{coverage}{BamFile}(x, shift=0L, width=NULL, weight=1L,
         ..., param=ScanBamParam())
}

\arguments{
  \item{x}{
    A \link{GAlignments}, \link{GAlignmentPairs}, or \link[Rsamtools]{BamFile}
    object.
  }
  \item{shift, width, weight}{
    See \code{coverage} method for \link[GenomicRanges]{GRanges} objects
    in the \pkg{GenomicRanges} package.
  }
  \item{method}{
    See \code{?\link[IRanges]{coverage}} in the \pkg{IRanges} package
    for a description of this argument.
  }
  \item{drop.D.ranges}{
    Whether the coverage calculation should ignore ranges corresponding
    to D (deletion) in the CIGAR string.
  }
  \item{...}{
    Additional arguments passed to the \code{coverage} method for
    \link{GAlignments} objects.
  }
  \item{param}{
    An optional \link[Rsamtools]{ScanBamParam} object passed to
    \code{\link{readGAlignmentsFromBam}}.
  }
}

\details{
  The methods for \link{GAlignments} and \link{GAlignmentPairs} objects do:
\preformatted{
  coverage(grglist(x), ...)
}
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
    \item \code{\link[IRanges]{coverage}} in the \pkg{IRanges} package.

    \item \link[GenomicRanges]{coverage-methods} in the
          \pkg{GenomicRanges} package.

    \item \link[IRanges]{RleList} objects in the \pkg{IRanges} package.

    \item \link{GAlignments} and \link{GAlignmentPairs} objects.

    \item \code{\link{readGAlignmentsFromBam}}.
  }
}

\examples{
ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")

gal <- readGAlignments(ex1_file)
stopifnot(identical(coverage(gal), coverage(as(gal, "GRangesList"))))

galp <- readGAlignmentPairs(ex1_file)
stopifnot(identical(coverage(galp), coverage(as(galp, "GRangesList"))))
}

\keyword{methods}
\keyword{utilities}
