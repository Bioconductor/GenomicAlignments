\name{intra-range-methods}

\alias{intra-range-methods}

\alias{narrow,GAlignments-method}
\alias{narrow,GAlignmentsList-method}
\alias{narrow,GappedReads-method}
\alias{qnarrow}
\alias{qnarrow,GAlignments-method}
\alias{qnarrow,GAlignmentsList-method}
\alias{qnarrow,GappedReads-method}


\title{Intra range transformations of a GAlignments or GAlignmentsList object}

\description{
  This man page documents intra range transformations of a
  \link[GenomicRanges]{GAlignments} or \link[GenomicRanges]{GAlignmentsList}
  object.

  See \code{?`\link[IRanges]{intra-range-methods}`} and
  \code{?`\link[IRanges]{inter-range-methods}`} in the \pkg{IRanges}
  package for a quick introduction to intra range and inter range
  transformations.

  Intra range methods for \link{GRanges} and \link{GRangesList}
  objects are defined and documented in the \pkg{GenomicRanges} package.
}

\usage{
\S4method{narrow}{GAlignments}(x, start=NA, end=NA, width=NA, use.names=TRUE)
\S4method{narrow}{GAlignmentsList}(x, start=NA, end=NA, width=NA, use.names=TRUE)

\S4method{qnarrow}{GAlignments}(x, start=NA, end=NA, width=NA)
\S4method{qnarrow}{GAlignmentsList}(x, start=NA, end=NA, width=NA)
}

\arguments{
  \item{x}{
    A \link[GenomicRanges]{GAlignments} or
    \link[GenomicRanges]{GAlignmentsList} object.
  }
  \item{start, end, width}{
    Vectors of integers.
    NAs and negative values are accepted and "solved" according to the
    rules of the SEW (Start/End/Width) interface (see
    \code{?\link[IRanges]{solveUserSEW}} for more information about the
    SEW interface).

    See \code{?`\link[IRanges]{intra-range-methods}`} for more information
    about the \code{start}, \code{end}, and \code{width} arguments.
  }
  \item{use.names}{
    See \code{?`\link[IRanges]{intra-range-methods}`}.
  }
}

\details{
  \itemize{
    \item(){
      \code{narrow} on a \link[GenomicRanges]{GAlignments} object behaves
      like on a \link[IRanges]{Ranges} object. See
      \code{?`\link[IRanges]{intra-range-methods}`} for the details.

      A major difference though is that it returns a
      \link[GenomicRanges]{GAlignments} object instead of a
      \link[IRanges]{Ranges} object.

      Unlike with \code{qnarrow} (see below), the
      \code{start}/\code{end}/\code{width} arguments here describe
      the narrowing on the reference side, not the query side.
    }
    \item(){
      \code{qnarrow} on a \link[GenomicRanges]{GAlignments} object
      behaves like \code{narrow} except that the
      \code{start}/\code{end}/\code{width} arguments here specify
      the narrowing with respect to the query sequences.

      \code{qnarrow} on a \link[GenomicRanges]{GAlignmentsList} object
      returns a \link[GenomicRanges]{GAlignmentsList} object.
    }
  }
}

\value{
  An object of the same class as, and \emph{parallel} to (i.e. same length
  and names as), the original object \code{x}.
}

\note{
  There is no difference between \code{narrow} and \code{qnarrow} when
  all the alignments have a simple CIGAR (i.e. no indels or gaps).
}

\author{H. Pages and V. Obenchain <vobencha@fhcrc.org>}

\seealso{
  \itemize{
    \item \link[GenomicRanges]{GAlignments} and
          \link[GenomicRanges]{GAlignmentsList} objects.
    \item The \link[IRanges]{intra-range-methods} man page in the
          \pkg{IRanges} package.
    \item The \link[GenomicRanges]{intra-range-methods} man page in the
          \pkg{GenomicRanges} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. ON A GAlignments OBJECT
## ---------------------------------------------------------------------
ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
gal <- readGAlignments(ex1_file, param=ScanBamParam(what="flag"))
gal

## This trims 3 nucleotides on the left and 5 nucleotides on the right
## of each alignment:
qnarrow(gal, start=4, end=-6)
## Note that the 'start' and 'end' arguments specify what part of each
## query sequence should be kept (negative values being relative to the
## right end of the query sequence), not what part should be trimmed.

## Trimming on the left doesn't change the "end" of the queries.
qnarrow(gal, start=21)
stopifnot(identical(end(qnarrow(gal, start=21)), end(gal)))

## ---------------------------------------------------------------------
## B. ON A GAlignmentsList OBJECT
## ---------------------------------------------------------------------
gal1 <- GAlignments(
    seqnames=Rle(factor(c("chr1", "chr2", "chr1", "chr3")),
        c(1, 3, 2, 4)),
    pos=1:10, cigar=paste0(10:1, "M"),
    strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    names=head(letters, 10), score=1:10)

gal2 <- GAlignments(
    seqnames=Rle(factor(c("chr2", "chr4")), c(3, 4)), pos=1:7,
    cigar=c("5M", "3M2N3M2N3M", "5M", "10M", "5M1N4M", "8M2N1M", "5M"),
    strand=Rle(strand(c("-", "+")), c(4, 3)),
    names=tail(letters, 7), score=1:7)

galist <- GAlignmentsList(noGaps=gal1, Gaps=gal2)
galist

qnarrow(galist)
}

\keyword{methods}
\keyword{utilities}
