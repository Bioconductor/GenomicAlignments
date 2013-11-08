\name{readGAlignments}

\alias{readGAlignments}
\alias{readGAlignmentPairs}
\alias{readGAlignmentsList}
\alias{readGappedReads}

\alias{readGAlignmentsFromBam}
\alias{readGAlignmentsFromBam,BamFile-method}
\alias{readGAlignmentsFromBam,character-method}
\alias{readGAlignmentsFromBam,BamViews-method}

\alias{readGAlignmentPairsFromBam}
\alias{readGAlignmentPairsFromBam,BamFile-method}
\alias{readGAlignmentPairsFromBam,character-method}

\alias{readGAlignmentsListFromBam}
\alias{readGAlignmentsListFromBam,BamFile-method}
\alias{readGAlignmentsListFromBam,character-method}

\alias{readGappedReadsFromBam}
\alias{readGappedReadsFromBam,BamFile-method}
\alias{readGappedReadsFromBam,character-method}

% Old stuff:
\alias{readGappedAlignments}
\alias{readGappedAlignmentPairs}

\alias{readBamGappedAlignments}
\alias{readBamGappedReads}
\alias{readBamGappedAlignmentPairs}
\alias{readBamGAlignmentsList}


\title{Reading alignments from a file}

\description{
  Read alignments from a file (typically a BAM file) into a
  \link[GenomicRanges]{GAlignments}, \link[GenomicRanges]{GAlignmentPairs},
  \link[GenomicRanges]{GAlignmentsList}, or \link{GappedReads}
  object.
}

\usage{
## Front-ends
readGAlignments(file, format="BAM", use.names=FALSE, ...)
readGAlignmentPairs(file, format="BAM", use.names=FALSE, ...)
readGAlignmentsList(file, format="BAM", use.names=FALSE, ...)
readGappedReads(file, format="BAM", use.names=FALSE, ...)

## BAM specific back-ends
readGAlignmentsFromBam(file, index=file, ..., use.names=FALSE,
                       param=NULL, with.which_label=FALSE)

readGAlignmentPairsFromBam(file, index=file, use.names=FALSE,
                           param=NULL, with.which_label=FALSE)

readGAlignmentsListFromBam(file, index=file, ..., use.names=FALSE,
                           param=ScanBamParam(), with.which_label=FALSE)

readGappedReadsFromBam(file, index=file, use.names=FALSE,
                       param=NULL, with.which_label=FALSE)
}

\arguments{
  \item{file, index}{The path to the BAM file to read, and to the index
    file of the BAM file to read, respectively. The latter is given
    \emph{without} the '.bai' extension. See \code{\link{scanBam}} for more
    information.}

  \item{\dots}{Arguments passed to other methods.}

  \item{use.names}{Use the query template names (QNAME field) as the
    names of the returned object? If not (the default), then the returned
    object has no names.}

  \item{param}{\code{NULL} or a \link[Rsamtools]{ScanBamParam} object.
    Like for \code{\link[Rsamtools]{scanBam}}, this influences what fields
    and which records are imported. However, note that the fields specified
    thru this \link[Rsamtools]{ScanBamParam} object will be loaded
    \emph{in addition} to any field required for generating the returned
    object (\link[GenomicRanges]{GAlignments},
    \link{GappedReads}, or \link[GenomicRanges]{GAlignmentPairs} object),
    but only the fields requested by the user will actually be kept as
    metadata columns of the object.

    By default (i.e. \code{param=NULL} or \code{param=ScanBamParam()}), no 
    additional field is loaded. The flag used is 
    \code{scanBamFlag(isUnmappedQuery=FALSE)} for
    \code{readGAlignmentsFromBam}, \code{readGappedReadsFromBam} and
    \code{readGAlignmentsListFromBam}
    (i.e. only records corresponding to mapped reads are loaded),
    and \code{scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE,
    hasUnmappedMate=FALSE)} for \code{readGAlignmentPairsFromBam}
    (i.e. only records corresponding to paired-end reads with both ends
    mapped are loaded).}

  \item{with.which_label}{\code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE} and if \code{param} has a \code{which} component,
    a \code{"which_label"} metadata column is added to the returned
    \link[GenomicRanges]{GAlignments} or \link{GappedReads} object,
    or to the \code{\link{first}} and \code{\link{last}} components
    of the returned \link[GenomicRanges]{GAlignmentPairs} object.
    In the case of \code{readGAlignmentsListFromBam}, it's added as an
    \emph{inner} metadata column, that is, the metadata column is placed
    on the \link[GenomicRanges]{GAlignments} object obtained by unlisting
    the returned \link[GenomicRanges]{GAlignmentsList} object.

    The purpose of this metadata column is to unambiguously identify
    the range in \code{which} where each element in the returned object
    originates from. The labels used to identify the ranges are normally
    of the form \code{"seq1:12250-246500"}, that is, they're the same as
    the names found on the outer list that \code{\link{scanBam}} would
    return if called with the same \code{param} argument. If some ranges
    are duplicated, then the labels are made unique by appending a unique
    suffix to all of them. The \code{"which_label"} metadata column is
    represented as a factor-\link[IRanges]{Rle}.}
}

\details{
  See \code{?\link[GenomicRanges]{GAlignments}} for a
  description of \link[GenomicRanges]{GAlignments} objects.

  See \code{?\link{GappedReads}} for a
  description of \link{GappedReads} objects.

  \subsection{Front-ends}{
      \code{readGAlignments} reads a file containing aligned reads as a
      \link[GenomicRanges]{GAlignments} object.

      \code{readGAlignmentPairs} reads a file containing paired-end reads as a
      \link[GenomicRanges]{GAlignmentPairs} object.

      \code{readGAlignmentsList} reads a file containing aligned reads as a
      \link[GenomicRanges]{GAlignmentsList} object.

      By default (i.e. \code{use.names=FALSE}), the resulting object has no
      names. If \code{use.names} is \code{TRUE}, then the names are
      constructed from the query template names (QNAME field in a SAM/BAM
      file). Note that the 2 records in a pair (when using
      \code{readGAlignmentPairs} or the records in a group (when using
      \code{readGAlignmentsList}) have the same QNAME.

      These functions are just front-ends that delegate to a
      format-specific back-end function depending on the supplied \code{format}
      argument. The \code{use.names} argument and any extra argument are
      passed to the back-end function.
      Only the BAM format is supported for now via the \code{read*FromBam}
      back-end functions.
  }

  \subsection{BAM specific back-ends}{

    When \code{file} is \link[Rsamtools]{BamViews} object
    \code{readGAlignmentsFromBam} visits each path in \code{bamPaths(file)},
    returning the result of \code{readGAlignmentsFromBam} applied to the
    specified path. When \code{index} is missing, it is set equal to
    \code{bamIndicies(file)}. Only reads in \code{bamRanges(file)} are
    returned (if \code{param} is supplied, \code{bamRanges(file)} takes
    precedence over \code{bamWhich(param)}).
    The return value is a \link[IRanges]{SimpleList} object, with elements
    of the list corresponding to each path. \code{bamSamples(file)} is
    available as metadata columns (accessed with \code{mcols}) of the
    returned \link[IRanges]{SimpleList} object.

    \code{readGAlignmentPairsFromBam} proceeds in 2 steps:
    \enumerate{
      \item Load the BAM file into a \link[GenomicRanges]{GAlignments}
            object with \code{readGAlignmentsFromBam};
      \item Turn this \link[GenomicRanges]{GAlignments} object into a
            \link[GenomicRanges]{GAlignmentPairs} object by pairing
            its elements.
    }
    See \code{?\link[GenomicRanges]{GAlignmentPairs-class}} for a
    description of \link[GenomicRanges]{GAlignmentPairs} objects,
    and \code{?\link{findMateAlignment}} for a description of the
    pairing algorithm (including timing and memory requirement).

    \code{readGAlignmentsListFromBam} pairs records into \sQuote{mates}
    acording to the criteria below. A \code{GAlignmentsList} is 
    returned with a \sQuote{mates} metadata column which indicates mate 
    status. The mates are returned first followed by non-mates. When
    the \sQuote{file} argument is a BamFile, \sQuote{asMates=TRUE}
    must be set, otherwise the data are treated as single-end reads. 
    See the \sQuote{asMates} section of \code{?\link{BamFile}} for details. 

      Mate criteria:
      \itemize{
        \item Bit 0x1 (multiple segments) is 1.
        \item Bit 0x4 (segment unmapped) is 0.
        \item Bit 0x8 (next segment unmapped) is 0.
        \item Bit 0x40 and 0x80 (first/last segment): 
               Segments are a pair of first/last OR
               neither segment is marked first/last.
        \item Bit 0x100 (secondary alignment): 
              Both segments are secondary OR both not secondary
        \item Bit 0x2 (properly aligned):
              Both segments are properly aligned 
        \item \sQuote{qname} match.
        \item \sQuote{tid} match.
        \item segment1 \sQuote{mpos} matches segment2 \sQuote{pos} AND
               segment2 \sQuote{mpos} matches segment1 \sQuote{pos}
      } 

      Records not passing these criteria are returned with mate status
      FALSE. Flags, tags and ranges may be specified in the
      \code{ScanBamParam} for fine tuning of results.

    See \code{?\link[GenomicRanges]{GAlignmentsList-class}} for a 
    description of \link[GenomicRanges]{GAlignmentsList} objects.
  }
}

\value{
  A \link[GenomicRanges]{GAlignments} object for
  \code{readGAlignmentsFromBam}.

  A \link{GappedReads} object for \code{readGappedReadsFromBam}.

  A \link[GenomicRanges]{GAlignmentPairs} object for
  \code{readGAlignmentPairsFromBam}.
  Note that a BAM (or SAM) file can in theory contain a mix of single-end
  and paired-end reads, but in practise it seems that single-end and
  paired-end are not mixed. In other words, the value of flag bit 0x1
  (\code{isPaired}) is the same for all the records in a file.
  So if \code{readGAlignmentPairsFromBam} returns a
  \link[GenomicRanges]{GAlignmentPairs} object of length zero,
  this almost certainly means that the BAM (or SAM) file contains
  alignments for single-end reads (although it could also mean that the
  user-supplied \code{\linkS4class{ScanBamParam}} is filtering out everything,
  or that the file is empty, or that all the records in the file correspond
  to unmapped reads).

  A \link[GenomicRanges]{GAlignmentsList} object for
  \code{readGAlignmentsListFromBam}. 
}

\note{
  BAM records corresponding to unmapped reads are always ignored.

  Starting with Rsamtools 1.7.1 (BioC 2.10), PCR or optical duplicates
  are loaded by default (use \code{scanBamFlag(isDuplicate=FALSE)} to
  drop them).
}

\author{H. Pages <hpages@fhcrc.org> and
        Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \link[GenomicRanges]{GAlignments-class},
  \link[GenomicRanges]{GAlignmentsList-class},
  \link{GappedReads-class},
  \link[GenomicRanges]{GAlignmentPairs-class},
  \code{\link{findMateAlignment}},
  \code{\link{scanBam}},
  \code{\link{ScanBamParam}}
}

\examples{
## ---------------------------------------------------------------------
## A. readGAlignmentsFromBam()
## ---------------------------------------------------------------------

## Simple use:
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                       mustWork=TRUE)
gal1 <- readGAlignmentsFromBam(bamfile)
gal1
names(gal1)

## Using the 'use.names' arg:
gal2 <- readGAlignmentsFromBam(bamfile, use.names=TRUE)
gal2
head(names(gal2))

## Using the 'param' arg to drop PCR or optical duplicates as well as
## secondary alignments, and to load additional BAM fields:
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
                                       isNotPrimaryRead=FALSE),
                      what=c("qual", "flag"))
gal3 <- readGAlignmentsFromBam(bamfile, param=param)
gal3
mcols(gal3)

## Using the 'param' arg to load reads from particular regions.
## Note that if we weren't providing a 'what' argument here, all the
## BAM fields would be loaded:
which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
param <- ScanBamParam(which=which)
gal4 <- readGAlignmentsFromBam(bamfile, param=param)
gal4

## Note that a given record is loaded one time for each region it
## belongs to (this is a scanBam() feature, readGAlignmentsFromBam()
## is based on scanBam()):
which <- IRangesList(seq2=IRanges(c(1563, 1567), width=1))
param <- ScanBamParam(which=which)
gal5 <- readGAlignmentsFromBam(bamfile, param=param)
gal5

## Use 'with.which_label=TRUE' to identify the range in 'which'
## where each element in 'gal5' originates from.
gal5 <- readGAlignmentsFromBam(bamfile, param=param,
                               with.which_label=TRUE)
gal5

## Using the 'param' arg to load tags. Except for MF and Aq, the tags
## specified below are predefined tags (see the SAM Spec for the list
## of predefined tags and their meaning).
param <- ScanBamParam(tag=c("MF", "Aq", "NM", "UQ", "H0", "H1"),
                      what="isize")
gal6 <- readGAlignmentsFromBam(bamfile, param=param)
mcols(gal6)  # "tag" cols always after "what" cols

## With a BamViews object:
fls <- system.file("extdata", "ex1.bam", package="Rsamtools",
                   mustWork=TRUE)
bv <- BamViews(fls,
               bamSamples=DataFrame(info="test", row.names="ex1"),
               auto.range=TRUE)
aln <- readGAlignmentsFromBam(bv)
aln
aln[[1]]
aln[colnames(bv)]
mcols(aln)

## ---------------------------------------------------------------------
## B. readGappedReadsFromBam()
## ---------------------------------------------------------------------
greads1 <- readGappedReadsFromBam(bamfile)
greads1
names(greads1)
qseq(greads1)
greads2 <- readGappedReadsFromBam(bamfile, use.names=TRUE)
head(greads2)
head(names(greads2))

## ---------------------------------------------------------------------
## C. readGAlignmentPairsFromBam()
## ---------------------------------------------------------------------
galp1 <- readGAlignmentPairsFromBam(bamfile)
head(galp1)
names(galp1)
## Using the 'param' arg to drop PCR or optical duplicates as well as
## secondary alignments (dropping secondary alignments can help make the
## pairing algorithm run significantly faster, see ?findMateAlignment):
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
                                       isNotPrimaryRead=FALSE))
galp2 <- readGAlignmentPairsFromBam(bamfile, use.names=TRUE, param=param)
galp2
head(galp2)
head(names(galp2))

## ---------------------------------------------------------------------
## D. readGAlignmentsListFromBam()
## ---------------------------------------------------------------------

library(pasillaBamSubset)

## 'file' as character.
fl <- untreated3_chr4() 
galist1 <- readGAlignmentsListFromBam(fl)
galist1[1:3]
length(galist1)
table(elementLengths(galist1))

## When 'file' is a BamFile, 'asMates' must be TRUE. If FALSE,
## the data are treated as single-end and each list element of the
## GAlignmentsList will be of length 1. For single-end data 
## use readGAlignments().
bf <- BamFile(fl, yieldSize=3, asMates=TRUE)
readGAlignmentsList(bf)

## Use a 'param' to fine tune the results.
param <- ScanBamParam(flag=scanBamFlag(isProperPair=TRUE))
galist2 <- readGAlignmentsListFromBam(fl, param=param)
length(galist2)
}

\keyword{manip}
