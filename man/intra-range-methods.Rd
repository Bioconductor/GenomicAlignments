\name{intra-range-methods}

\alias{intra-range-methods}

\alias{narrow}
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
  \link{GAlignments} or \link{GAlignmentsList} object.

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
    A \link{GAlignments} or \link{GAlignmentsList} object.
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
      \code{narrow} on a \link{GAlignments} object behaves
      like on a \link[IRanges]{IntegerRanges} object. See
      \code{?`\link[IRanges]{intra-range-methods}`} for the details.

      A major difference though is that it returns a \link{GAlignments}
      object instead of a \link[IRanges]{IntegerRanges} object.

      Unlike with \code{qnarrow} (see below), the
      \code{start}/\code{end}/\code{width} arguments here describe
      the narrowing on the reference side, not the query side.
    }
    \item(){
      \code{qnarrow} on a \link{GAlignments} object behaves like \code{narrow}
      except that the \code{start}/\code{end}/\code{width} arguments here
      specify the narrowing with respect to the query sequences.

      \code{qnarrow} on a \link{GAlignmentsList} object
      returns a \link{GAlignmentsList} object.
    }
  }
}

\value{
  An object of the same class as, and \emph{parallel} to (i.e. same length
  and names as), the original object \code{x}.
}

\note{
  There is no difference between \code{narrow} and \code{qnarrow} when
  all the alignments have a simple CIGAR (i.e. no indels or junctions).
}

\author{Hervé Pagès and V. Obenchain <Valerie.Obenchain@RoswellPark.org>}

\seealso{
  \itemize{
    \item \link{GAlignments} and \link{GAlignmentsList} objects.

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
param <- ScanBamParam(what=c("seq", "qual"))
gal <- readGAlignments(ex1_file, param=param)
gal

## This trims 3 nucleotides on the left and 5 nucleotides on the right
## of each alignment:
gal2 <- qnarrow(gal, start=4, end=-6)
gal2

## Note that the 'start' and 'end' values are relative to the query
## sequences and specify the query substring that must be kept for each
## alignment. Negative values are relative to the right end of the query
## sequence.

## Also note that the metadata columns on 'gal' are propagated as-is so
## the "seq" and "qual" matadata columns must be adjusted "by hand" with
## narrow();
mcols(gal2)$seq <- narrow(mcols(gal)$seq, start=4, end=-6)
mcols(gal2)$qual <- narrow(mcols(gal)$qual, start=4, end=-6)
gal2

## Sanity checks:
stopifnot(identical(qwidth(gal2), width(mcols(gal2)$seq)))
stopifnot(identical(qwidth(gal2), width(mcols(gal2)$qual)))

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
