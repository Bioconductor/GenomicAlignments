\name{GappedReads-class}
\docType{class}

% Class:
\alias{GappedReads-class}
\alias{GappedReads}

% Constructors:
\alias{readGappedReads}

% Accessors:
\alias{qseq}
\alias{qseq,GappedReads-method}
\alias{qwidth,GappedReads-method}

% Combining:
\alias{c,GappedReads-method}

% Other methods:
\alias{qnarrow,GappedReads-method}
\alias{narrow,GappedReads-method}


\title{(Legacy) GappedReads objects}

\description{
  The GappedReads class extends the \link[GenomicRanges]{GAlignments}
  class defined in the GenomicRanges package.
  A GappedReads object contains all the information contained in a
  \link[GenomicRanges]{GAlignments} object plus the sequences
  of the queries. Those sequences can be accessed via the \code{qseq}
  accessor.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{readGappedReads(file, format="BAM", use.names=FALSE, ...)}:
      Read a file as a GappedReads object.
      Like with the \code{\link[GenomicRanges]{readGAlignments}}
      constructor, by default (i.e. \code{use.names=FALSE}) the resulting
      object has no names. If \code{use.names} is \code{TRUE}, then the
      names are constructed from the query template names (QNAME field in
      a SAM/BAM file).

      Note that this function is just a front-end that delegates to the
      format-specific back-end function specified via the \code{format}
      argument. The \code{use.names} argument and any extra argument are
      passed to the back-end function.
      Only the BAM format is supported for now. Its back-end is the
      \code{\link[Rsamtools]{readGappedReadsFromBam}} function defined
      in the Rsamtools package.
      See \code{?\link[Rsamtools]{readGappedReadsFromBam}} for more
      information (you might need to install and load the Rsamtools package
      first).
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GappedReads object.

  \describe{
    \item{}{
      \code{qseq(x)}:
      Extracts the sequences of the queries as a
      \code{\link[Biostrings]{DNAStringSet}} object.
    }
  }
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{
  H. Pages and P. Aboyoun
}

\seealso{
  \link[GenomicRanges]{GAlignments-class},
  \code{\link[Rsamtools]{readGappedReadsFromBam}}
}

\examples{
greads_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
greads <- readGappedReads(greads_file)
greads
qseq(greads)
}

\keyword{methods}
\keyword{classes}
