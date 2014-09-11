\name{stackStringsFromBam}

\alias{stackStringsFromBam}
\alias{alphabetFrequencyFromBam}


\title{Stack the read sequences stored in a BAM file on a region of interest}

\description{
  \code{stackStringsFromBam} stacks the read sequences (or their quality
  strings) stored in a BAM file over a user-specified region.

  \code{alphabetFrequencyFromBam} computes the alphabet frequency of the
  reads over a user-specified region.

  Both functions take into account the CIGAR of each read to "lay" the read
  sequence (or its quality string) alongside the reference space.
  This step ensures that each nucleotide in a read is associated with the
  correct position on the reference sequence.
}

\usage{
stackStringsFromBam(file, index=file, param,
                    what="seq", use.names=FALSE,
                    D.letter="-", N.letter=".",
                    Lpadding.letter="+", Rpadding.letter="+")

alphabetFrequencyFromBam(file, index=file, param, what="seq", ...)
}

\arguments{
  \item{file, index}{
    The path to the BAM file to read, and to the index file of the BAM file
    to read, respectively. The latter is given \emph{without} the '.bai'
    extension. See \code{\link{scanBam}} for more information.
  }
  \item{param}{
    A \link{ScanBamParam} object containing exactly 1 genomic region
    (i.e. \code{unlist(bamWhich(param))} must have length 1).
    Alternatively, \code{param} can be a \link[GenomicRanges]{GRanges} or
    \link[IRanges]{RangesList} object containing exactly 1 genomic region
    (the strand will be ignored in case of a \link[GenomicRanges]{GRanges}
    object), or a character string specifying a single genomic region (in
    the \code{"chr14:5201-5300"} format).
  }
  \item{what}{
    A single string. Either \code{"seq"} or \code{"qual"}. If \code{"seq"}
    (the default), the read sequences will be stacked. If \code{"qual"},
    the read quality strings will be stacked.
  }
  \item{use.names}{
    Use the query template names (QNAME field) as the names of the returned
    object? If not (the default), then the returned object has no names.
  }
  \item{D.letter, N.letter}{
    A single letter used as a filler for injections. The 2 arguments are
    passed down to the \code{\link{sequenceLayer}} function.
    See \code{?\link{sequenceLayer}} for more details.
  }
  \item{Lpadding.letter, Rpadding.letter}{
    A single letter to use for padding the sequences on the left, and another
    one to use for padding on the right. The 2 arguments are passed down to
    the \code{\link[Biostrings]{stackStrings}} function defined in the
    \pkg{Biostrings} package.
    See \code{?\link[Biostrings]{stackStrings}} in the \pkg{Biostrings}
    package for more details.
  }
  \item{...}{
    Further arguments to be passed to \link[Biostrings]{alpahabetFrequency}.
  }
}

\details{
  \code{stackStringsFromBam} performs the 3 following steps:
  \enumerate{
    \item Load the read sequences (or their quality strings) from the BAM
          file. Only the read sequences that overlap with the specified region
          are loaded. This is done with the
          \code{\link{readGAlignmentsFromBam}} function. Note that if the
          file contains paired-end reads, the pairing is ignored.

    \item Lay the sequences alongside the reference space, using their CIGARs.
          This is done with the \code{\link{sequenceLayer}} function.

    \item Stack them on the specified region. This is done with the
          \code{\link[Biostrings]{stackStrings}} function defined in the
          \pkg{Biostrings} package.
  }

  \code{alphabetFrequencyFromBam} also performs steps 1. and 2. but,
  instead of stacking the sequences at step 3., it computes the nucleotide
  frequencies for each genomic position in the specified region.
}

\value{
  For \code{stackStringsFromBam}: A rectangular (i.e. constant-width)
  \link[Biostrings]{DNAStringSet} object (if \code{what} is \code{"seq"})
  or \link[Biostrings]{BStringSet} object (if \code{what} is \code{"qual"}).

  For \code{alphabetFrequencyFromBam}: By default a matrix like one returned
  by \code{\link[Biostrings]{alphabetFrequency}}. The matrix has 1 row per
  nucleotide position in the specified region.
}

\note{
  TWO IMPORTANT CAVEATS ABOUT \code{stackStringsFromBam}:

  Specifying a big genomic region, say >= 100000 bp, can require a lot of
  memory (especially with high coverage reads) and is not recommended.
  See the \code{\link{pileLettersAt}} function for piling the read letters
  on top of a set of individual genomic positions, which is more
  flexible and more memory efficient.

  Paired-end reads are treated as single-end reads (i.e. they're not paired).
}

\author{H. Pages}

\seealso{
  \itemize{
    \item The \code{\link{pileLettersAt}} function for piling the letters
          of a set of aligned reads on top of a set of individual genomic
          positions.

    \item The \code{\link{readGAlignmentsFromBam}} function for loading read
          sequences (or their quality strings) from a BAM file (via a
          \link{GAlignments} object).

    \item The \code{\link{sequenceLayer}} function for laying read sequences
          alongside the reference space, using their CIGARs.

    \item The \code{\link[Biostrings]{stackStrings}} function in the
          \pkg{Biostrings} package for stacking an arbitrary
          \link[Biostrings]{XStringSet} object.

    \item The \code{\link[Biostrings]{alphabetFrequency}} function in the
          \pkg{Biostrings} package.

    \item The SAMtools mpileup command available at
          \url{http://samtools.sourceforge.net/} as part of the
          SAMtools project.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. EXAMPLE WITH TOY DATA
## ---------------------------------------------------------------------

bamfile1 <- BamFile(system.file("extdata", "ex1.bam", package="Rsamtools"))

region1 <- GRanges("seq1", IRanges(1, 60))  # region of interest

## Stack the read sequences:
stackStringsFromBam(bamfile1, param=region1)
af <- alphabetFrequencyFromBam(bamfile1, param=region1, baseOnly=TRUE)
cm1a <- t(af[ , DNA_BASES])  # consensus matrix
cm1a

## Stack their quality strings:
stackStringsFromBam(bamfile1, param=region1, what="qual")

## Control the number of reads to display:
options(showHeadLines=18)
options(showTailLines=6)
stackStringsFromBam(bamfile1, param=GRanges("seq1", IRanges(61, 120)))

stacked_qseq <- stackStringsFromBam(bamfile1, param="seq2:1509-1519")
stacked_qseq  # deletion in read 13
af <- alphabetFrequencyFromBam(bamfile1, param="seq2:1509-1519",
                                baseOnly=TRUE)
cm1b <- t(af[ , DNA_BASES])
cm1b

## Sanity check:
stopifnot(identical(consensusMatrix(stacked_qseq)[DNA_BASES, ], cm1b))

stackStringsFromBam(bamfile1, param="seq2:1509-1519", what="qual")

## ---------------------------------------------------------------------
## B. EXAMPLE WITH REAL DATA
## ---------------------------------------------------------------------

library(RNAseqData.HNRNPC.bam.chr14)
bamfile2 <- BamFile(RNAseqData.HNRNPC.bam.chr14_BAMFILES[1])

## Region of interest:
region2 <- GRanges("chr14", IRanges(19650095, 19650159))

readGAlignments(bamfile2, param=ScanBamParam(which=region2))

stackStringsFromBam(bamfile2, param=region2)

af <- alphabetFrequencyFromBam(bamfile2, param=region2, baseOnly=TRUE)
cm2 <- t(af[ , DNA_BASES])  # consensus matrix
cm2
}

\keyword{methods}
\keyword{manip}
