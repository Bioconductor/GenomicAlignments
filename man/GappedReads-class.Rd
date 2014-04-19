\name{GappedReads-class}
\docType{class}

% Class:
\alias{class:GappedReads}
\alias{GappedReads-class}

% Constructor:
\alias{GappedReads}

% Accessors:
\alias{qseq}
\alias{qseq,GappedReads-method}
\alias{qwidth,GappedReads-method}

% Combining:
\alias{c,GappedReads-method}


\title{(Legacy) GappedReads objects}

\description{
  The GappedReads class extends the \link{GAlignments} class.

  A GappedReads object contains all the information contained in a
  \link{GAlignments} object plus the sequences of the queries.
  Those sequences can be accessed via the \code{qseq} accessor.
}

\section{Constructor}{
  GappedReads objects are typically created when reading a file
  containing aligned reads with the \code{\link{readGappedReads}}
  function.
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
  \itemize{
    \item \link{GAlignments} objects.

    \item \code{\link{readGappedReads}}.
  }
}

\examples{
greads_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
greads <- readGappedReads(greads_file)
greads
qseq(greads)
}

\keyword{methods}
\keyword{classes}
