\name{encodeOverlaps-methods}
\alias{encodeOverlaps-methods}

\alias{encodeOverlaps}
\alias{encodeOverlaps,RangesList,RangesList-method}
\alias{encodeOverlaps,RangesList,Ranges-method}
\alias{encodeOverlaps,Ranges,RangesList-method}

\alias{encodeOverlaps1}

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
  \link[GenomicRanges]{GRangesList} objects, which computes the "overlap
  encodings" between a \code{query} and a \code{subject}, both list-like
  objects with list elements containing multiple ranges.

  Other related utilities are also documented in this man page.
}

\usage{
encodeOverlaps(query, subject, hits=NULL, ...)

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
  \item{query, subject}{
    Typically \link[GenomicRanges]{GRangesList} objects representing the
    the aligned reads and the transcripts of a given gene model, respectively.
    If the 2 objects don't have the same length, and if the \code{hits}
    argument is not supplied, then the shortest is recycled to the length
    of the longest (the standard recycling rules apply).

    More generally speaking, \code{query} and \code{subject} must be
    list-like objects with list elements containing multiple ranges e.g.
    \link[IRanges]{RangesList} or \link[GenomicRanges]{GRangesList} objects.
  }
  \item{hits}{
    An optional \link[IRanges]{Hits} object typically obtained from a
    previous call to \code{\link[IRanges]{findOverlaps}(query, subject)}.

    Strictly speaking, \code{hits} only needs to be compatible with
    \code{query} and \code{subject}, that is,
    \code{\link[IRanges]{queryLength}(hits)} and
    \code{\link[IRanges]{subjectLength}(hits)} must be equal to
    \code{length(query)} and \code{length(subject)}, respectively.

    Supplying \code{hits} is a convenient way to do
    \code{encodeOverlaps(query[queryHits(hits)], subject[subjectHits(hits)])},
    that is, calling \code{encodeOverlaps(query, subject, hits)} is equivalent
    to the above, but is much more efficient, especially when \code{query}
    and/or \code{subject} are big. Of course, when \code{hits} is supplied,
    \code{query} and \code{subject} are not expected to have the same length
    anymore.
  }
  \item{...}{
    Additional arguments for methods.
  }
  \item{flip.query.if.wrong.strand}{
    See the "OverlapEncodings" vignette located in this package
    (\pkg{GenomicAlignments}).
  }
  \item{x}{
    For \code{flipQuery}: a \link[GenomicRanges]{GRangesList} object.

    For \code{isCompatibleWithSplicing}, \code{isCompatibleWithSkippedExons},
    \code{extractSteppedExonRanks}, \code{extractSpannedExonRanks}, and
    \code{extractSkippedExonRanks}: an \link{OverlapEncodings} object,
    a factor, or a character vector.
  }
  \item{i}{
    Subscript specifying the elements in \code{x} to flip. If missing, all
    the elements are flipped.
  }
  \item{ovencA, ovencB, ovenc}{
    \link{OverlapEncodings} objects.
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
  See \code{?OverlapEncodings} for a short introduction to "overlap encodings".

  The topic of working with overlap encodings is covered in details
  in the "OverlapEncodings" vignette located this package
  (\pkg{GenomicAlignments}) and accessible with
  \code{vignette("OverlapEncodings")}.
}

\value{
  For \code{encodeOverlaps}: An \link{OverlapEncodings} object.
  If \code{hits} is not supplied, this object is \emph{parallel} to the
  longest of \code{query} and \code{subject}, that is, it has the length
  of the longest and the i-th encoding in it corresponds to the i-th element
  in the longest.
  If \code{hits} is supplied, then the returned object is \emph{parallel}
  to it, that is, it has one encoding per hit.

  For \code{flipQuery}: TODO

  For \code{selectEncodingWithCompatibleStrand}: TODO

  For \code{isCompatibleWithSplicing} and \code{isCompatibleWithSkippedExons}:
  A logical vector \emph{parallel} to \code{x}.

  For \code{extractSteppedExonRanks}, \code{extractSpannedExonRanks}, and
  \code{extractSkippedExonRanks}: TODO

  For \code{extractQueryStartInTranscript}: TODO
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \link{OverlapEncodings} class for a brief introduction to
          "overlap encodings".

    \item The \link[IRanges]{Hits} class defined and documented in the
          \pkg{IRanges} package.

    \item The "OverlapEncodings" vignette in this package.

    \item \code{\link{findCompatibleOverlaps}} for a specialized version
          of \code{\link[IRanges]{findOverlaps}} that uses
          \code{encodeOverlaps} internally to keep only the hits where
          the junctions in the aligned read are \emph{compatible} with
          the splicing of the annotated transcript.

    \item The \link[GenomicRanges]{GRangesList} class defined and documented
          in the \pkg{GenomicRanges} package.

    \item The \code{\link[IRanges]{findOverlaps}} generic function defined
          in the \pkg{IRanges} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. BETWEEN 2 RangesList OBJECTS
## ---------------------------------------------------------------------
## In the context of an RNA-seq experiment, encoding the overlaps
## between 2 GRangesList objects, one containing the reads (the query),
## and one containing the transcripts (the subject), can be used for
## detecting hits between reads and transcripts that are "compatible"
## with the splicing of the transcript. Here we illustrate this with 2
## RangesList objects, in order to keep things simple:

## 4 aligned reads in the query:
read1 <- IRanges(c(7, 15, 22), c(9, 19, 23))  # 2 junctions
read2 <- IRanges(c(5, 15), c(9, 17))  # 1 junction
read3 <- IRanges(c(16, 22), c(19, 24))  # 1 junction
read4 <- IRanges(c(16, 23), c(19, 24))  # 1 junction
query <- IRangesList(read1, read2, read3, read4)

## 1 transcript in the subject:
tx <- IRanges(c(1, 4, 15, 22, 38), c(2, 9, 19, 25, 47))  # 5 exons
subject <- IRangesList(tx)

## Encode the overlaps:
ovenc <- encodeOverlaps(query, subject)
ovenc
encoding(ovenc)

## Reads that are "compatible" with the transcript can be detected with
## a regular expression (the regular expression below assumes that
## reads have at most 2 junctions):
regex0 <- "(:[fgij]:|:[jg].:.[gf]:|:[jg]..:.g.:..[gf]:)"
grepl(regex0, encoding(ovenc))  # read4 is NOT "compatible"

## This was for illustration purpose only. In practise you don't need
## (and should not) use this regular expression, but use instead the
## isCompatibleWithSplicing() utility function:
isCompatibleWithSplicing(ovenc)

## ---------------------------------------------------------------------
## B. BETWEEN 2 GRangesList OBJECTS
## ---------------------------------------------------------------------
## With real RNA-seq data, the reads and transcripts will typically be
## stored in GRangesList objects. Please refer to the "OverlapEncodings"
## vignette in this package for realistic examples.
}

\keyword{methods}
\keyword{utilities}
