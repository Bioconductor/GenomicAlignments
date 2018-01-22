\name{findOverlaps-methods}
\alias{findOverlaps-methods}

\alias{findOverlaps}
\alias{findOverlaps,GAlignments,Vector-method}
\alias{findOverlaps,Vector,GAlignments-method}
\alias{findOverlaps,GAlignments,GAlignments-method}
\alias{findOverlaps,GAlignmentPairs,Vector-method}
\alias{findOverlaps,Vector,GAlignmentPairs-method}
\alias{findOverlaps,GAlignmentPairs,GAlignmentPairs-method}
\alias{findOverlaps,GAlignmentsList,Vector-method}
\alias{findOverlaps,Vector,GAlignmentsList-method}
\alias{findOverlaps,GAlignmentsList,GAlignmentsList-method}


\title{Finding overlapping genomic alignments}

\description{
  Finds range overlaps between a \link{GAlignments}, \link{GAlignmentPairs}, 
  or \link{GAlignmentsList} object, and another range-based object.

  NOTE: The \code{\link[IRanges]{findOverlaps}} generic function and methods
  for \link[IRanges]{IntegerRanges} and \link[IRanges]{IntegerRangesList}
  objects are defined and documented in the \pkg{IRanges} package.
  The methods for \link[GenomicRanges]{GRanges} and
  \link[GenomicRanges]{GRangesList} objects are defined and
  documented in the \pkg{GenomicRanges} package.

  \link{GAlignments}, \link{GAlignmentPairs}, and \link{GAlignmentsList}
  objects also support \code{countOverlaps}, \code{overlapsAny}, and
  \code{subsetByOverlaps} thanks to the default methods defined in the
  \pkg{IRanges} package and to the \code{findOverlaps} method defined in
  this package and documented below.
}

\usage{
\S4method{findOverlaps}{GAlignments,GAlignments}(query, subject,
    maxgap=-1L, minoverlap=0L,
    type=c("any", "start", "end", "within"),
    select=c("all", "first", "last", "arbitrary"),
    ignore.strand=FALSE)
}

\arguments{
  \item{query, subject}{
    A \link{GAlignments}, \link{GAlignmentPairs}, or \link{GAlignmentsList}
    object for either \code{query} or \code{subject}.
    A vector-like object containing ranges for the other one.
  }
  \item{maxgap, minoverlap, type, select}{
    See \code{?\link[IRanges]{findOverlaps}} in the \pkg{IRanges} package
    for a description of these arguments.
  }
  \item{ignore.strand}{
    When set to \code{TRUE}, the strand information is ignored in the
    overlap calculations.
  }
}

\details{
  When the query or the subject (or both) is a \link{GAlignments}
  object, it is first turned into a \link{GRangesList} object (with
  \code{as( , "GRangesList")}) and then the rules described previously
  apply. \link{GAlignmentsList} objects are coerced to \link{GAlignments}
  then to a \link{GRangesList}. Feature indices are mapped back to the 
  original \link{GAlignmentsList} list elements.

  When the query is a \link{GAlignmentPairs} object, it is first
  turned into a \link{GRangesList} object (with \code{as( , "GRangesList")})
  and then the rules described previously apply.
}

\value{
  A \link[S4Vectors]{Hits} object when \code{select="all"} or an integer
  vector otherwise. 
}

\seealso{
  \itemize{
    \item \code{\link[IRanges]{findOverlaps}}.
    \item \link[S4Vectors]{Hits-class}.
    \item \link{GRanges-class}.
    \item \link{GRangesList-class}.
    \item \link{GAlignments-class}.
    \item \link{GAlignmentPairs-class}.
    \item \link{GAlignmentsList-class}.
  }
}

\examples{
ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn <- readGAlignments(ex1_file)

subject <- granges(galn)[1]

## Note the absence of query no. 9 (i.e. 'galn[9]') in this result:
as.matrix(findOverlaps(galn, subject))

## This is because, by default, findOverlaps()/countOverlaps() are
## strand specific:
galn[8:10]
countOverlaps(galn[8:10], subject)
countOverlaps(galn[8:10], subject, ignore.strand=TRUE)

## Count alignments in 'galn' that DO overlap with 'subject' vs those
## that do NOT:
table(overlapsAny(galn, subject))
## Extract those that DO:
subsetByOverlaps(galn, subject)

## GAlignmentsList
galist <- GAlignmentsList(galn[8:10], galn[3000:3002])
gr <- GRanges(c("seq1", "seq1", "seq2"), 
              IRanges(c(15, 18, 1233), width=1),
              strand=c("-", "+", "+"))

countOverlaps(galist, gr)
countOverlaps(galist, gr, ignore.strand=TRUE)
findOverlaps(galist, gr)
findOverlaps(galist, gr, ignore.strand=TRUE)
}

\keyword{methods}
\keyword{utilities}
