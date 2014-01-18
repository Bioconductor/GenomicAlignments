\name{encodeOverlaps-methods}
\alias{encodeOverlaps-methods}

\alias{flipQuery}

\alias{encodeOverlaps,GRangesList,GRangesList-method}

\alias{selectEncodingWithCompatibleStrand}

\alias{isCompatibleWithSplicing}
\alias{isCompatibleWithSplicing,character-method}
\alias{isCompatibleWithSplicing,factor-method}
\alias{isCompatibleWithSplicing,OverlapEncodings-method}

\alias{isCompatibleWithSkippedExons}
\alias{isCompatibleWithSkippedExons,character-method}
\alias{isCompatibleWithSkippedExons,factor-method}
\alias{isCompatibleWithSkippedExons,OverlapEncodings-method}

\alias{extractSteppedExonRanks}
\alias{extractSteppedExonRanks,character-method}
\alias{extractSteppedExonRanks,factor-method}
\alias{extractSteppedExonRanks,OverlapEncodings-method}

\alias{extractSpannedExonRanks}
\alias{extractSpannedExonRanks,character-method}
\alias{extractSpannedExonRanks,factor-method}
\alias{extractSpannedExonRanks,OverlapEncodings-method}

\alias{extractSkippedExonRanks}
\alias{extractSkippedExonRanks,character-method}
\alias{extractSkippedExonRanks,factor-method}
\alias{extractSkippedExonRanks,OverlapEncodings-method}

\alias{extractQueryStartInTranscript}


\title{Encode the overlaps between RNA-seq reads and the transcripts of
       a gene model}

\description{
  In the context of an RNA-seq experiment, encoding the overlaps between
  the aligned reads and the transcripts of a given gene model can be used
  for detecting those overlaps that are \emph{compatible} with the splicing
  of the transcript.

  The central tool for this is the \code{encodeOverlaps} method for
  \link[GenomicRanges]{GRangesList} objects.

  Othe related utilities are also documented in this man page.
}

\usage{
\S4method{encodeOverlaps}{GRangesList,GRangesList}(query, subject, hits=NULL,
               flip.query.if.wrong.strand=FALSE)

## Related utilities:

flipQuery(x, i)

selectEncodingWithCompatibleStrand(ovencA, ovencB,
                                   query.strand, subject.strand, hits=NULL)

isCompatibleWithSplicing(x)
isCompatibleWithSkippedExons(x, max.skipped.exons=NA)

extractSteppedExonRanks(x, for.query.right.end=FALSE)
extractSpannedExonRanks(x, for.query.right.end=FALSE)
extractSkippedExonRanks(x, for.query.right.end=FALSE)

extractQueryStartInTranscript(query, subject, hits=NULL, ovenc=NULL,
                              flip.query.if.wrong.strand=FALSE,
                              for.query.right.end=FALSE)
}

\arguments{
  \item{query}{
    A \link[GenomicRanges]{GRangesList} object representing the aligned
    reads.
  }
  \item{subject}{
    A \link[GenomicRanges]{GRangesList} object representing the transcripts
    of a given gene model.
  }
  \item{hits}{
    A \link[IRanges]{Hits} object.
    See \code{?`\link[IRanges]{encodeOverlaps}`} for a description of
    how a supplied \link[IRanges]{Hits} object is handled.
  }
  \item{flip.query.if.wrong.strand}{
    See the "OverlapEncodings" vignette in this package
    (\pkg{GenomicAlignments}).
  }
  \item{x}{
    For \code{flipQuery}: a \link[GenomicRanges]{GRangesList} object.

    For \code{isCompatibleWithSplicing}, \code{isCompatibleWithSkippedExons},
    \code{extractSteppedExonRanks}, \code{extractSpannedExonRanks}, and
    \code{extractSkippedExonRanks}:
    an \link[IRanges]{OverlapEncodings} object, a factor or a character
    vector.
  }
  \item{i}{
    Subscript specifying the elements in \code{x} to flip. If missing, all
    the elements are flipped.
  }
  \item{ovencA, ovencB, ovenc}{
    \link[IRanges]{OverlapEncodings} objects.
  }
  \item{query.strand, subject.strand}{
    Vector-like objects containing the strand of the query and subject,
    respectively.
  }
  \item{max.skipped.exons}{
    Not supported yet. If \code{NA} (the default), the number of skipped
    exons must be 1 or more (there is no max).
  }
  \item{for.query.right.end}{
    If \code{TRUE}, then the information reported in the output is for
    the right ends of the paired-end reads.
    Using \code{for.query.right.end=TRUE} with single-end reads is an error.
  }
}

\details{
  The topic of working with overlap encodings is covered in details
  in the "OverlapEncodings" vignette located this package
  (\pkg{GenomicAlignments}) and accessible with
  \code{vignette("OverlapEncodings")}.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \code{\link[GenomicAlignments]{findCompatibleOverlaps}} for a
          specialized version of \code{\link[IRanges]{findOverlaps}} that
          uses \code{\link[IRanges]{encodeOverlaps}} internally to keep only
          the hits where the gaps in the aligned read are \emph{compatible}
          with the splicing of the transcript.

    \item The "OverlapEncodings" vignette in this package.

    \item The \code{\link[IRanges]{encodeOverlaps}} generic function and
          \link[IRanges]{OverlapEncodings} class defined in the \pkg{IRanges}
          package.

    \item The \code{\link[IRanges]{findOverlaps}} generic function defined
          in the \pkg{IRanges} package.

    \item \link{GRangesList} objects.
  }
}

\examples{
## Please refer to the "OverlapEncodings" vignette in the
## GenomicAlignments package for a comprehensive presentation
## of "overlap encodings" and related tools/concepts (e.g.
## "compatible" overlaps, "almost compatible" overlaps etc...),
## and for examples.
}

\keyword{methods}
\keyword{utilities}
