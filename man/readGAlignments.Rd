\name{readGAlignments}

\alias{readGAlignments}
\alias{readGAlignments,BamFile-method}
\alias{readGAlignments,character-method}
\alias{readGAlignments,BamViews-method}

\alias{readGAlignmentPairs}
\alias{readGAlignmentPairs,BamFile-method}
\alias{readGAlignmentPairs,character-method}

\alias{readGAlignmentsList}
\alias{readGAlignmentsList,BamFile-method}
\alias{readGAlignmentsList,character-method}

\alias{readGappedReads}
\alias{readGappedReads,BamFile-method}
\alias{readGappedReads,character-method}

% old stuff (Deprecated or Defunct)
\alias{readGAlignmentsFromBam}
\alias{readGAlignmentPairsFromBam}
\alias{readGAlignmentsListFromBam}
\alias{readGappedReadsFromBam}


\title{Reading genomic alignments from a file}

\description{
  Read genomic alignments from a file (typically a BAM file) into a
  \link{GAlignments}, \link{GAlignmentPairs}, \link{GAlignmentsList},
  or \link{GappedReads} object.
}

\usage{
readGAlignments(file, index=file, use.names=FALSE, param=NULL,
                      with.which_label=FALSE)

readGAlignmentPairs(file, index=file, use.names=FALSE, param=NULL,
                          with.which_label=FALSE)

readGAlignmentsList(file, index=file, use.names=FALSE, param=ScanBamParam(),
                          with.which_label=FALSE)

readGappedReads(file, index=file, use.names=FALSE, param=NULL,
                      with.which_label=FALSE)
}

\arguments{
  \item{file}{
    The path to the file to read or a \link[Rsamtools]{BamFile} object.
    Can also be a \link[Rsamtools]{BamViews} object for
    \code{readGAlignments}.
  }
  \item{index}{
    The path to the index file of the BAM file to read.
    Must be given \emph{without} the '.bai' extension.
    See \code{\link[Rsamtools]{scanBam}} in the \pkg{Rsamtools} packages
    for more information.
  }
  \item{use.names}{
    \code{TRUE} or \code{FALSE}. By default (i.e. \code{use.names=FALSE}),
    the resulting object has no names. If \code{use.names} is \code{TRUE},
    then the names are constructed from the query template names (QNAME
    field in a SAM/BAM file). Note that the 2 records in a pair (when using
    \code{readGAlignmentPairs} or the records in a group (when using
    \code{readGAlignmentsList}) have the same QNAME.
  }
  \item{param}{\code{NULL} or a \link[Rsamtools]{ScanBamParam} object.
    Like for \code{\link[Rsamtools]{scanBam}}, this influences what fields
    and which records are imported. However, note that the fields specified
    thru this \link[Rsamtools]{ScanBamParam} object will be loaded
    \emph{in addition} to any field required for generating the returned
    object (\link{GAlignments}, \link{GAlignmentPairs},
    or \link{GappedReads} object),
    but only the fields requested by the user will actually be kept as
    metadata columns of the object.

    By default (i.e. \code{param=NULL} or \code{param=ScanBamParam()}), no 
    additional field is loaded. The flag used is 
    \code{scanBamFlag(isUnmappedQuery=FALSE)} for
    \code{readGAlignments}, \code{readGAlignmentsList}, and
    \code{readGappedReads}.
    (i.e. only records corresponding to mapped reads are loaded),
    and \code{scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE,
    hasUnmappedMate=FALSE)} for \code{readGAlignmentPairs}
    (i.e. only records corresponding to paired-end reads with both ends
    mapped are loaded).}

  \item{with.which_label}{\code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE} and if \code{param} has a \code{which} component,
    a \code{"which_label"} metadata column is added to the returned
    \link{GAlignments} or \link{GappedReads} object,
    or to the \code{\link{first}} and \code{\link{last}} components
    of the returned \link{GAlignmentPairs} object.
    In the case of \code{readGAlignmentsList}, it's added as an
    \emph{inner} metadata column, that is, the metadata column is placed
    on the \link{GAlignments} object obtained by unlisting
    the returned \link{GAlignmentsList} object.

    The purpose of this metadata column is to unambiguously identify
    the range in \code{which} where each element in the returned object
    originates from. The labels used to identify the ranges are normally
    of the form \code{"seq1:12250-246500"}, that is, they're the same as
    the names found on the outer list that \code{\link{scanBam}} would
    return if called with the same \code{param} argument. If some ranges
    are duplicated, then the labels are made unique by appending a unique
    suffix to all of them. The \code{"which_label"} metadata column is
    represented as a factor-\link[S4Vectors]{Rle}.}
}

\details{
  \itemize{
    \item \code{readGAlignments} reads a file containing aligned reads as a
          \link{GAlignments} object. See \code{?\link{GAlignments}} for a
          description of \link{GAlignments} objects.

          When \code{file} is a \link[Rsamtools]{BamViews} object,
          \code{readGAlignments} visits each path in \code{bamPaths(file)},
          returning the result of \code{readGAlignments} applied to the
          specified path. When \code{index} is missing, it is set equal to
          \code{bamIndicies(file)}. Only reads in \code{bamRanges(file)} are
          returned (if \code{param} is supplied, \code{bamRanges(file)} takes
          precedence over \code{bamWhich(param)}). The return value is a
          \link[S4Vectors]{SimpleList} object, with elements of the list
          corresponding to each path. \code{bamSamples(file)} is available
          as metadata columns (accessed with \code{mcols}) of the returned
          \link[S4Vectors]{SimpleList} object.

    \item \code{readGAlignmentPairs} reads a file containing aligned paired-end
          reads as a \link{GAlignmentPairs} object.
          See \code{?\link{GAlignmentPairs}} for a description of
          \link{GAlignmentPairs} objects.

    \item \code{readGAlignmentsList} reads a file containing aligned reads as
          a \link{GAlignmentsList} object. See \code{?\link{GAlignmentsList}} 
          for a description of \link{GAlignmentsList} objects.
          \code{readGAlignmentsList} pairs records into mates according to the
          pairing criteria described below. The 1st mate will always be 1st in
          the \link{GAlignmentsList} list elements that have mate_status set to
          \code{"mated"}, and the 2nd mate will always be 2nd.

          A \code{GAlignmentsList} is returned with a \sQuote{mate_status}
          metadata column on the outer list elements. \code{mate_status} is a
          factor with 3 levels indicating mate status, \sQuote{mated},
          \sQuote{ambiguous} or \sQuote{unmated}:
          \itemize{
            \item{mated:} primary or non-primary pairs 
            \item{ambiguous:} multiple segments matching to the
                              same location (indistinguishable) 
            \item{unmated:} mate does not exist or is unmapped
          }

          When the \sQuote{file} argument is a BamFile, \sQuote{asMates=TRUE}
          must be set, otherwise the data are treated as single-end reads. 
          See the \sQuote{asMates} section of \code{?\link[Rsamtools]{BamFile}}
          in the \pkg{Rsamtools} package for details. 

    \item \code{readGappedReads} reads a file containing aligned reads as a
          \link{GappedReads} object. See \code{?\link{GappedReads}} for a
          description of \link{GappedReads} objects.
  }

  For all these functions, flags, tags and ranges may be specified in the
  supplied \link[Rsamtools]{ScanBamParam} object for fine tuning of results.
}

\section{Pairing criteria}{
  This section describes the pairing criteria used by
  \code{readGAlignmentsList} and \code{readGAlignmentPairs}.
  \itemize{
    \item First, only records with flag bit 0x1 (multiple segments) set to 1,
          flag bit 0x4 (segment unmapped) set to 0, and flag bit 0x8 (next
          segment in the template unmapped) set to 0, are candidates for
          pairing (see the SAM Spec for a description of flag bits and fields).
          Records that correspond to single-end reads, or records that
          correspond to paired-end reads where one or both ends are unmapped,
          will remain unmated.

    \item Then the following fields and flag bits are considered:
          \itemize{
            \item (A) QNAME
            \item (B) RNAME, RNEXT
            \item (C) POS, PNEXT
            \item (D) Flag bits Ox10 (segment aligned to minus strand)
                      and 0x20 (next segment aligned to minus strand)
            \item (E) Flag bits 0x40 (first segment in template) and 0x80 (last
                      segment in template)
            \item (F) Flag bit 0x2 (proper pair)
            \item (G) Flag bit 0x100 (secondary alignment)
          }
          2 records rec1 and rec2 are considered mates iff all the following
          conditions are satisfied:
          \itemize{
            \item (A) QNAME(rec1) == QNAME(rec2)
            \item (B) RNEXT(rec1) == RNAME(rec2) and RNEXT(rec2) == RNAME(rec1)
            \item (C) PNEXT(rec1) == POS(rec2) and PNEXT(rec2) == POS(rec1)
            \item (D) Flag bit 0x20 of rec1 == Flag bit 0x10 of rec2 and
                      Flag bit 0x20 of rec2 == Flag bit 0x10 of rec1
            \item (E) rec1 corresponds to the first segment in the template and
                      rec2 corresponds to the last segment in the template, OR,
                      rec2 corresponds to the first segment in the template and
                      rec1 corresponds to the last segment in the template
            \item (F) rec1 and rec2 have same flag bit 0x2
            \item (G) rec1 and rec2 have same flag bit 0x100
          }
  }
  Note that this is actually the pairing criteria used by
  \code{\link[Rsamtools]{scanBam}} (when the \link[Rsamtools]{BamFile}
  passed to it has the \code{asMates} toggle set to \code{TRUE}), which
  \code{readGAlignmentsList} and \code{readGAlignmentPairs} call behind
  the scene. It is also the pairing criteria used by
  \code{\link{findMateAlignment}}.
}

\value{
  A \link{GAlignments} object for \code{readGAlignments}.

  A \link{GAlignmentPairs} object for \code{readGAlignmentPairs}.
  Note that a BAM (or SAM) file can in theory contain a mix of single-end
  and paired-end reads, but in practise it seems that single-end and
  paired-end are not mixed. In other words, the value of flag bit 0x1
  (\code{isPaired}) is the same for all the records in a file.
  So if \code{readGAlignmentPairs} returns a \link{GAlignmentPairs} object
  of length zero, this almost always means that the BAM (or SAM) file
  contains alignments for single-end reads (although it could also mean that
  the user-supplied \code{\linkS4class{ScanBamParam}} is filtering out
  everything, or that the file is empty, or that all the records in the file
  correspond to unmapped reads).

  A \link{GAlignmentsList} object for \code{readGAlignmentsList}. 
  When the list contains paired-end reads a metadata data column of
  \code{mate_status} is added to the object. See details in the 
  `Bam specific back-ends' section on this man page. 

  A \link{GappedReads} object for \code{readGappedReads}.
}

\note{
  BAM records corresponding to unmapped reads are always ignored.

  Starting with Rsamtools 1.7.1 (BioC 2.10), PCR or optical duplicates
  are loaded by default (use \code{scanBamFlag(isDuplicate=FALSE)} to
  drop them).
}

\author{H. Pages <hpages@fredhutch.org> and
        Valerie Obenchain <vobencha@fredhutch.org>}

\seealso{
  \itemize{
    \item \code{\link[Rsamtools]{scanBam}} and
          \code{\link[Rsamtools]{ScanBamParam}} in the \pkg{Rsamtools}
          package.

    \item \link{GAlignments}, \link{GAlignmentPairs}, \link{GAlignmentsList},
          and \link{GappedReads} objects.

    \item \link[IRanges]{RangesList} objects (used in the examples below
          to specify the \code{which} regions) in the \pkg{IRanges} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. readGAlignments()
## ---------------------------------------------------------------------

## Simple use:
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                       mustWork=TRUE)
gal1 <- readGAlignments(bamfile)
gal1
names(gal1)

## Using the 'use.names' arg:
gal2 <- readGAlignments(bamfile, use.names=TRUE)
gal2
head(names(gal2))

## Using the 'param' arg to drop PCR or optical duplicates as well as
## secondary alignments, and to load additional BAM fields:
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
                                       isSecondaryAlignment=FALSE),
                      what=c("qual", "flag"))
gal3 <- readGAlignments(bamfile, param=param)
gal3
mcols(gal3)

## Using the 'param' arg to load alignments from particular regions.
## Note that if we weren't providing a 'what' argument here, all the
## BAM fields would be loaded:
which <- RangesList(seq1=IRanges(1000, 1100),
                    seq2=IRanges(c(1546, 1555, 1567), width=10))
param <- ScanBamParam(which=which)
gal4 <- readGAlignments(bamfile, use.names=TRUE, param=param)
gal4

## IMPORTANT NOTE: A given record is loaded one time for each region
## it overlaps with. We call this "duplicated record selection" (this
## is a scanBam() feature, readGAlignments() is based on scanBam()):
which <- IRangesList(seq2=IRanges(c(1555, 1567), width=10))
param <- ScanBamParam(which=which)
gal5 <- readGAlignments(bamfile, use.names=TRUE, param=param)
gal5  # record EAS114_26:7:37:79:581 was loaded twice

## This becomes clearer if we use 'with.which_label=TRUE' to identify
## the region in 'which' where each element in 'gal5' originates from.
gal5 <- readGAlignments(bamfile, use.names=TRUE, param=param,
                                 with.which_label=TRUE)
gal5

## Not surprisingly, we also get "duplicated record selection" when
## 'which' contains repeated or overlapping regions. Using the same
## regions as we did for 'gal4' above, except that now we're
## repeating the region on seq1:
which <- RangesList(seq1=rep(IRanges(1000, 1100), 2),
                    seq2=IRanges(c(1546, 1555, 1567), width=10))
param <- ScanBamParam(which=which)
gal4b <- readGAlignments(bamfile, use.names=TRUE, param=param)
length(gal4b)  # > length(gal4), because all the records overlapping
               # with bases 1000 to 1100 on seq1 are now duplicated

## The "duplicated record selection" will artificially increase the
## coverage or affect other downstream results. It can be mitigated
## (but not completely eliminated) by first "reducing" the set of
## regions:
which <- reduce(which)
which
param <- ScanBamParam(which=which)
gal4c <- readGAlignments(bamfile, use.names=TRUE, param=param)
length(gal4c)  # < length(gal4), because the 2 first original regions
               # on seq2 were merged into a single one

## Note that reducing the set of regions didn't completely eliminate
## "duplicated record selection". Records that overlap the 2 reduced
## regions on seq2 (which$seq2) are loaded twice (like for 'gal5'
## above).

## Using the 'param' arg to load tags. Except for MF and Aq, the tags
## specified below are predefined tags (see the SAM Spec for the list
## of predefined tags and their meaning).
param <- ScanBamParam(tag=c("MF", "Aq", "NM", "UQ", "H0", "H1"),
                      what="isize")
gal6 <- readGAlignments(bamfile, param=param)
mcols(gal6)  # "tag" cols always after "what" cols

## With a BamViews object:
fls <- system.file("extdata", "ex1.bam", package="Rsamtools",
                   mustWork=TRUE)
bv <- BamViews(fls,
               bamSamples=DataFrame(info="test", row.names="ex1"),
               auto.range=TRUE)
## Note that the "readGAlignments" method for BamViews objects
## requires the ShortRead package to be installed.
aln <- readGAlignments(bv)
aln
aln[[1]]
aln[colnames(bv)]
mcols(aln)

## ---------------------------------------------------------------------
## B. readGAlignmentPairs()
## ---------------------------------------------------------------------
galp1 <- readGAlignmentPairs(bamfile)
head(galp1)
names(galp1)
## Using the 'param' arg to drop PCR or optical duplicates as well as
## secondary alignments (dropping secondary alignments can help make the
## pairing algorithm run significantly faster):
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
                                       isSecondaryAlignment=FALSE))
galp2 <- readGAlignmentPairs(bamfile, use.names=TRUE, param=param)
galp2
head(galp2)
head(names(galp2))

## ---------------------------------------------------------------------
## C. readGAlignmentsList()
## ---------------------------------------------------------------------
library(pasillaBamSubset)

## 'file' as character.
fl <- untreated3_chr4() 
galist1 <- readGAlignmentsList(fl)
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
galist2 <- readGAlignmentsList(fl, param=param)
length(galist2)

## ---------------------------------------------------------------------
## D. readGappedReads()
## ---------------------------------------------------------------------
greads1 <- readGappedReads(bamfile)
greads1
names(greads1)
qseq(greads1)
greads2 <- readGappedReads(bamfile, use.names=TRUE)
head(greads2)
head(names(greads2))
}

\keyword{manip}
