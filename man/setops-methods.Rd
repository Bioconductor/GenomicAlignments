\name{setops-methods}
\alias{setops-methods}

\alias{pintersect,GAlignments,GRanges-method}
\alias{pintersect,GRanges,GAlignments-method}


\title{Set operations on GAlignments objects}

\description{
  Performs set operations on \link{GAlignments} objects.
}

\usage{
\S4method{pintersect}{GAlignments,GRanges}(x, y, ...)
\S4method{pintersect}{GRanges,GAlignments}(x, y, ...)
}

\arguments{
  \item{x, y}{
    A \link{GAlignments} object and a \link[GenomicRanges]{GRanges} object.
    They must have the same length.
  }
  \item{...}{
    Further arguments to be passed to or from other methods.
  }
}

\value{
  A \link{GAlignments} object \emph{parallel} to (i.e. same length as)
  \code{x} and \code{y}.
}

\seealso{
  \itemize{
    \item The \link{GAlignments} class.

    \item The \link[GenomicRanges]{setops-methods} man page in the
          \pkg{GenomicRanges} package.
  }
}

\examples{
## Parallel intersection of a GAlignments and a GRanges object:
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
gal <- readGAlignments(bamfile)
pintersect(gal, shift(as(gal, "GRanges"), 6L))
}

\keyword{methods}
\keyword{utilities}
