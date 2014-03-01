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
junctions(x, use.mcols=FALSE, ...)

\S4method{junctions}{GAlignments}(x, use.mcols=FALSE)

\S4method{junctions}{GAlignmentPairs}(x, use.mcols=FALSE)

\S4method{junctions}{GAlignmentsList}(x, use.mcols=FALSE, ignore.strand=FALSE)
}

\arguments{
  \item{x}{
    A \link{GAlignments}, \link{GAlignmentPairs}, or \link{GAlignmentsList}
    object.
  }
  \item{use.mcols}{
    \code{TRUE} or \code{FALSE} (the default).
    Whether the metadata columns on \code{x} (accessible with \code{mcols(x)})
    should be propagated to the returned object or not.
  }
  \item{...}{
    Additional arguments, for use in specific methods.
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
  \link[GenomicRanges]{GRangesList} object \emph{parallel} to \code{x},
  (i.e. with 1 list element per element in \code{x}).
  If \code{x} has names on it, they're propagated to the returned object.
  If \code{use.mcols} is TRUE and \code{x} has metadata columns on it
  (accessible with \code{mcols(x)}), they're propagated to the returned object.
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
library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[1]

## ---------------------------------------------------------------------
## junctions()
## ---------------------------------------------------------------------

gal <- readGAlignments(bamfile, param=ScanBamParam(what="flag"))
table(njunc(gal))  # some alignments have 3 junctions!
junctions(gal)
stopifnot(identical(unname(elementLengths(junctions(gal))), njunc(gal)))

galp <- readGAlignmentPairs(bamfile)
junctions(galp)
stopifnot(identical(unname(elementLengths(junctions(galp))), njunc(galp)))

## ---------------------------------------------------------------------
## summarizeJunctions()
## ---------------------------------------------------------------------

}

\keyword{methods}
\keyword{manip}
