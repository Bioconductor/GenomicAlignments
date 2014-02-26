\name{junctions-methods}

\alias{junctions-methods}

\alias{junctions}
\alias{junctions,GAlignments-method}
\alias{junctions,GAlignmentPairs-method}
\alias{junctions,GAlignmentsList-method}

% Old stuff:
\alias{introns}


\title{Extract junctions from genomic alignments}

\description{
  Tools for extracting junctions from a \link{GAlignments},
  \link{GAlignmentPairs}, or \link{GAlignmentsList} object.
}

\usage{
junctions(x, ...)

\S4method{junctions}{GAlignments}(x)

\S4method{junctions}{GAlignmentPairs}(x)

\S4method{junctions}{GAlignmentsList}(x, ignore.strand=FALSE)
}

\arguments{
  \item{x}{
    A \link{GAlignments}, \link{GAlignmentPairs}, or \link{GAlignmentsList}
    object.
  }
  \item{...}{
    Additional arguments supported by specific \code{junctions} methods.
  }
  \item{ignore.strand}{
  }
}

\details{
  An N operation in the CIGAR of a genomic alignment is interpreted as a
  junction. \code{junctions(x)} will return the genomic ranges of all
  junctions found in \code{x}.

  More precisely, if \code{x} is a \link{GAlignments} object,
  \code{junctions(x)} is equivalent to:
\preformatted{
  psetdiff(granges(x), grglist(x, order.as.in.query=TRUE))
}
  On a \code{x} is a \link{GAlignmentPairs} object, it's equivalent to (but
  faster than):
\preformatted{
  junctions1 <- junctions(first(x))
  junctions2 <- junctions(last(x, invert.strand=TRUE))
  mendoapply(c, junctions1, junctions2)
}
}

\value{
  \code{junctions(x)} returns the genomic ranges of the junctions in a
  \link[GenomicRanges]{GRangesList} object \emph{parallel} to \code{x}
  (i.e. with 1 list element per element in \code{x}).
}

\author{H. Pages}

\seealso{
  \itemize{
    \item \link{GAlignments}, \link{GAlignmentPairs}, and
          \link{GAlignmentsList} objects.

    \item The \link[GenomicRanges]{GRanges} and
          \link[GenomicRanges]{GRangesList} classes defined and documented
          in the \pkg{GenomicRanges} package.
  }
}

\examples{
library(Rsamtools)  # for the ex1.bam file
ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")

gal <- readGAlignments(ex1_file, param=ScanBamParam(what="flag"))
junctions(gal)
stopifnot(identical(unname(elementLengths(junctions(gal))), njunc(gal)))

galp <- readGAlignmentPairs(ex1_file, use.names=TRUE)
junctions(galp)
stopifnot(identical(unname(elementLengths(junctions(galp))), njunc(galp)))
}

\keyword{methods}
\keyword{manip}
