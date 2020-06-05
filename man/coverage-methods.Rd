\name{coverage-methods}
\alias{coverage-methods}

\alias{coverage}
\alias{coverage,GAlignments-method}
\alias{coverage,GAlignmentPairs-method}
\alias{coverage,GAlignmentsList-method}
\alias{coverage,BamFile-method}
\alias{coverage,character-method}


\title{Coverage of a GAlignments, GAlignmentPairs, or GAlignmentsList object}

\description{
  \code{\link[IRanges]{coverage}} methods for \link{GAlignments},
  \link{GAlignmentPairs}, \link{GAlignmentsList}, and
  \link[Rsamtools]{BamFile} objects.

  NOTE: The \code{\link[IRanges]{coverage}} generic function and methods
  for \link[IRanges]{IntegerRanges} and \link[IRanges]{IntegerRangesList}
  objects are defined and documented in the \pkg{IRanges} package.
  Methods for \link[GenomicRanges]{GRanges} and
  \link[GenomicRanges]{GRangesList} objects are defined and
  documented in the \pkg{GenomicRanges} package.
}

\usage{
\S4method{coverage}{GAlignments}(x, shift=0L, width=NULL, weight=1L,
         method=c("auto", "sort", "hash", "naive"), drop.D.ranges=FALSE)

\S4method{coverage}{GAlignmentPairs}(x, shift=0L, width=NULL, weight=1L,
         method=c("auto", "sort", "hash", "naive"), drop.D.ranges=FALSE)

\S4method{coverage}{GAlignmentsList}(x, shift=0L, width=NULL, weight=1L, ...)

\S4method{coverage}{BamFile}(x, shift=0L, width=NULL, weight=1L, ...,
         param=ScanBamParam())

\S4method{coverage}{character}(x, shift=0L, width=NULL, weight=1L, ...,
         yieldSize=2500000L)
}

\arguments{
  \item{x}{
    A \link{GAlignments}, \link{GAlignmentPairs}, \link{GAlignmentsList},
    or \link[Rsamtools]{BamFile} object, or the path to a BAM file.
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
    \code{\link{readGAlignments}}.
  }
  \item{yieldSize}{
    An optional argument controlling how many records are input
    when iterating through a \link[Rsamtools]{BamFile}.
  }
}

\details{
  The methods for \link{GAlignments} and \link{GAlignmentPairs} objects
  do:
\preformatted{  coverage(grglist(x, drop.D.ranges=drop.D.ranges), ...)
}

  The method for \link{GAlignmentsList} objects does:
\preformatted{  coverage(unlist(x), ...)
}

  The method for \link[Rsamtools]{BamFile} objects iterates through a
  BAM file, reading \code{yieldSize(x)} records (or all records, if
  \code{is.na(yieldSize(x))}) and calculating:
\preformatted{  gal <- readGAlignments(x, param=param)
  coverage(gal, shift=shift, width=width, weight=weight, ...)
}

  The method for \code{character} vectors of length 1 creates a
  \link[Rsamtools]{BamFile} object from \code{x} and performs the
  calculation for \code{coverage,BamFile-method}.
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

    \item \code{\link{readGAlignments}}.

    \item \link[Rsamtools]{BamFile} objects in the \pkg{Rsamtools} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. EXAMPLE WITH TOY DATA
## ---------------------------------------------------------------------

ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")

## Coverage of a GAlignments object:
gal <- readGAlignments(ex1_file)
cvg1 <- coverage(gal)
cvg1

## Coverage of a GAlignmentPairs object:
galp <- readGAlignmentPairs(ex1_file)
cvg2 <- coverage(galp)
cvg2

## Coverage of a GAlignmentsList object:
galist <- readGAlignmentsList(ex1_file)
cvg3 <- coverage(galist)
cvg3

table(mcols(galist)$mate_status)
mated_idx <- which(mcols(galist)$mate_status == "mated")
mated_galist <- galist[mated_idx]
mated_cvg3 <- coverage(mated_galist)
mated_cvg3

## Sanity checks:
stopifnot(identical(cvg1, cvg3))
stopifnot(identical( cvg2, mated_cvg3))

## ---------------------------------------------------------------------
## B. EXAMPLE WITH REAL DATA
## ---------------------------------------------------------------------

library(pasillaBamSubset)
## See '?pasillaBamSubset' for more information about the 2 BAM files
## included in this package.
reads <- readGAlignments(untreated3_chr4())
table(njunc(reads))  # data contains junction reads

## Junctions do NOT contribute to the coverage:
read1 <- reads[which(njunc(reads) != 0L)[1]]  # 1st read with a junction
read1  # cigar shows a "skipped region" of length 15306
grglist(read1)[[1]]  # the junction is between pos 4500 and 19807
coverage(read1)$chr4  # junction is not covered

## Sanity checks:
cvg <- coverage(reads)
read_chunks <- unlist(grglist(reads), use.names=FALSE)
read_chunks_per_chrom <- split(read_chunks, seqnames(read_chunks))
stopifnot(identical(sum(cvg), sum(width(read_chunks_per_chrom))))

galist <- readGAlignmentsList(untreated3_chr4())
stopifnot(identical(cvg, coverage(galist)))
}
\keyword{methods}
\keyword{utilities}
