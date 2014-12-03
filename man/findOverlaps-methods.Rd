\name{findOverlaps-methods}
\alias{findOverlaps-methods}

\alias{findOverlaps}
\alias{findOverlaps,GAlignments,Vector-method}
\alias{findOverlaps,Vector,GAlignments-method}
\alias{findOverlaps,GAlignments,GAlignments-method}
\alias{findOverlaps,GAlignments,GRangesList-method}
\alias{findOverlaps,GRangesList,GAlignments-method}
\alias{findOverlaps,GAlignmentPairs,Vector-method}
\alias{findOverlaps,Vector,GAlignmentPairs-method}
\alias{findOverlaps,GAlignmentPairs,GAlignmentPairs-method}
\alias{findOverlaps,GAlignmentsList,Vector-method}
\alias{findOverlaps,Vector,GAlignmentsList-method}
\alias{findOverlaps,GAlignmentsList,GAlignmentsList-method}

\alias{countOverlaps}
\alias{countOverlaps,GAlignments,Vector-method}
\alias{countOverlaps,Vector,GAlignments-method}
\alias{countOverlaps,GAlignments,GAlignments-method}
\alias{countOverlaps,GAlignments,GenomicRanges-method}
\alias{countOverlaps,GenomicRanges,GAlignments-method}
\alias{countOverlaps,GAlignments,GRangesList-method}
\alias{countOverlaps,GRangesList,GAlignments-method}
\alias{countOverlaps,GAlignmentPairs,Vector-method}
\alias{countOverlaps,Vector,GAlignmentPairs-method}
\alias{countOverlaps,GAlignmentPairs,GAlignmentPairs-method}
\alias{countOverlaps,GAlignmentsList,Vector-method}
\alias{countOverlaps,Vector,GAlignmentsList-method}
\alias{countOverlaps,GAlignmentsList,GAlignmentsList-method}

\alias{overlapsAny}
\alias{overlapsAny,GAlignments,Vector-method}
\alias{overlapsAny,Vector,GAlignments-method}
\alias{overlapsAny,GAlignments,GAlignments-method}
\alias{overlapsAny,GAlignmentPairs,Vector-method}
\alias{overlapsAny,Vector,GAlignmentPairs-method}
\alias{overlapsAny,GAlignmentPairs,GAlignmentPairs-method}
\alias{overlapsAny,GAlignmentsList,Vector-method}
\alias{overlapsAny,Vector,GAlignmentsList-method}
\alias{overlapsAny,GAlignmentsList,GAlignmentsList-method}

\alias{subsetByOverlaps}
\alias{subsetByOverlaps,GAlignments,Vector-method}
\alias{subsetByOverlaps,Vector,GAlignments-method}
\alias{subsetByOverlaps,GAlignments,GAlignments-method}
\alias{subsetByOverlaps,GAlignmentPairs,Vector-method}
\alias{subsetByOverlaps,Vector,GAlignmentPairs-method}
\alias{subsetByOverlaps,GAlignmentPairs,GAlignmentPairs-method}
\alias{subsetByOverlaps,GAlignmentsList,Vector-method}
\alias{subsetByOverlaps,Vector,GAlignmentsList-method}
\alias{subsetByOverlaps,GAlignmentsList,GAlignmentsList-method}


\title{Finding overlapping genomic alignments}

\description{
  Finds range overlaps between a \link{GAlignments}, \link{GAlignmentPairs}, 
  or \link{GAlignmentsList} object, and another range-based object.

  NOTE: The \code{\link[IRanges]{findOverlaps}} generic function and methods
  for \link[IRanges]{Ranges} and \link[IRanges]{RangesList} objects
  are defined and documented in the \pkg{IRanges} package.
  The methods for \link[GenomicRanges]{GRanges} and
  \link[GenomicRanges]{GRangesList} objects are defined and
  documented in the \pkg{GenomicRanges} package.
}

\usage{
\S4method{findOverlaps}{GAlignments,GAlignments}(query, subject,
    maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within"),
    select=c("all", "first", "last", "arbitrary"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=FALSE)

\S4method{countOverlaps}{GAlignments,GAlignments}(query, subject,
    maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=FALSE)

\S4method{overlapsAny}{GAlignments,GAlignments}(query, subject,
    maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=FALSE)

\S4method{subsetByOverlaps}{GAlignments,GAlignments}(query, subject,
    maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=FALSE)
}

\arguments{
  \item{query, subject}{
    A \link{GAlignments}, \link{GAlignmentPairs}, or \link{GAlignmentsList}
    object for either \code{query} or \code{subject}.
    A vector-like object containing ranges for the other one.
  }
  \item{maxgap, minoverlap, type, select}{
    See \code{\link[IRanges]{findOverlaps}} in the \pkg{IRanges} package
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
  For \code{findOverlaps} either a \link[S4Vectors]{Hits} object when
  \code{select = "all"} or an integer vector otherwise. 

  For \code{countOverlaps} an integer vector containing the tabulated
  query overlap hits.

  For \code{overlapsAny} a logical vector of length equal to the number of
  ranges in \code{query} indicating those that overlap any of the ranges
  in \code{subject}.

  For \code{subsetByOverlaps} an object of the same class as \code{query}
  containing the subset that overlapped at least one entity in \code{subject}.
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
