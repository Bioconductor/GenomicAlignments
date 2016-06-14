\name{pileLettersAt}

\alias{pileLettersAt}


\title{Pile the letters of a set of aligned reads on top of a set
       of genomic positions}

\description{
  \code{pileLettersAt} extracts the letters/nucleotides of a set of
  reads that align to a set of genomic positions of interest.
  The extracted letters are returned as "piles of letters" (one per
  genomic position of interest) stored in an \link[Biostrings]{XStringSet}
  (typically \link[Biostrings]{DNAStringSet}) object.
}

\usage{
pileLettersAt(x, seqnames, pos, cigar, at)
}

\arguments{
  \item{x}{
    An \link[Biostrings]{XStringSet} (typically
    \link[Biostrings]{DNAStringSet}) object containing N \emph{unaligned}
    read sequences (a.k.a. the query sequences) reported with respect to
    the + strand. 
  }
  \item{seqnames}{
    A factor-\link[S4Vectors]{Rle} \emph{parallel} to \code{x}.
    For each \code{i}, \code{seqnames[i]} must be the name of the reference
    sequence of the i-th alignment.
  }
  \item{pos}{
    An integer vector \emph{parallel} to \code{x}.
    For each \code{i}, \code{pos[i]} must be the 1-based position
    on the reference sequence of the first aligned letter in \code{x[[i]]}.
  }
  \item{cigar}{
    A character vector \emph{parallel} to \code{x}. Contains the extended
    CIGAR strings of the alignments.
  }
  \item{at}{
    A \link[GenomicRanges]{GPos} object containing the genomic positions
    of interest. \code{seqlevels(at)} must be identical to
    \code{levels(seqnames)}.

    If \code{at} is not a \link[GenomicRanges]{GPos} object,
    \code{pileLettersAt} will first try to turn it into one by calling
    the \code{\link[GenomicRanges]{GPos}()} constructor function on it.
    So for example \code{at} can be a \link[GenomicRanges]{GRanges} object
    (or any other \link[GenomicRanges]{GenomicRanges} derivative), and, in
    that case, each range in it will be interpreted as a run of adjacent
    genomic positions. See \code{?\link[GenomicRanges]{GPos}} in the
    \pkg{GenomicRanges} package for more information.
  }
}

\details{
  \code{x}, \code{seqnames}, \code{pos}, \code{cigar} must be 4 \emph{parallel}
  vectors describing N aligned reads.
}

\value{
  An \link[Biostrings]{XStringSet} (typically \link[Biostrings]{DNAStringSet})
  object \emph{parallel} to \code{at} (i.e. with 1 string per genomic
  position).
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item The \code{pileup} and \code{applyPileups} functions defined in the
          \pkg{Rsamtools} package, as well as the SAMtools mpileup command
          (available at \url{http://samtools.sourceforge.net/} as part of the
          SAMtools project), for more powerful flexible alternatives.

    \item The \code{\link{stackStringsFromBam}} function
          for stacking the read sequences (or their quality strings)
          stored in a BAM file on a region of interest.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item \link[GenomicRanges]{GPos} objects in the \pkg{GenomicRanges}
          package.

    \item \link{GAlignments} objects.

    \item \link{cigar-utils} for the CIGAR utility functions used internally
          by \code{pileLettersAt}.
  }
}

\examples{
## Input

##   - A BAM file:
bamfile <- BamFile(system.file("extdata", "ex1.bam", package="Rsamtools"))
seqinfo(bamfile)  # to see the seqlevels and seqlengths
stackStringsFromBam(bamfile, param="seq1:1-21")  # a quick look at
                                                 # the reads

##   - A GPos object containing Genomic Positions Of Interest:
my_GPOI <- GPos(c("seq1:1-5", "seq1:21-21", "seq1:1575-1575",
                  "seq2:1513-1514"))

## Some preliminary massage on 'my_GPOI'

seqinfo(my_GPOI) <- merge(seqinfo(my_GPOI), seqinfo(bamfile))
seqlevels(my_GPOI) <- seqlevelsInUse(my_GPOI)

## Load the BAM file in a GAlignments object. Note that we load only
## the reads aligned to the sequences in 'seqlevels(my_GPOI)'. Also,
## in order to be consistent with applyPileups() and SAMtools (m)pileup,
## we filter out the following BAM records:
##   - secondary alignments (flag bit 0x100);
##   - reads not passing quality controls (flag bit 0x200);
##   - PCR or optical duplicates (flag bit 0x400).
## See ?ScanBamParam and the SAM Spec for more information. 

which <- as(seqinfo(my_GPOI), "GRanges")
flag <- scanBamFlag(isSecondaryAlignment=FALSE,
                    isNotPassingQualityControls=FALSE,
                    isDuplicate=FALSE)
what <- c("seq", "qual")
param <- ScanBamParam(flag=flag, what=c("seq", "qual"), which=which)
gal <- readGAlignments(bamfile, param=param)
seqlevels(gal) <- seqlevels(my_GPOI) 

## Extract the read sequences (a.k.a. query sequences) and quality
## strings. Both are reported with respect to the + strand.

qseq <- mcols(gal)$seq
qual <- mcols(gal)$qual

nucl_piles <- pileLettersAt(qseq, seqnames(gal), start(gal), cigar(gal),
                            my_GPOI)
qual_piles <- pileLettersAt(qual, seqnames(gal), start(gal), cigar(gal),
                            my_GPOI)
mcols(my_GPOI)$nucl_piles <- nucl_piles
mcols(my_GPOI)$qual_piles <- qual_piles
my_GPOI 

## Finally, to summarize A/C/G/T frequencies at each position:
alphabetFrequency(nucl_piles, baseOnly=TRUE)

## Note that the pileup() function defined in the Rsamtools package
## can be used to obtain a similar result:
scanbam_param <- ScanBamParam(flag=flag, which=my_GPOI)
pileup_param <- PileupParam(max_depth=5000,
                            min_base_quality=0,
                            distinguish_strands=FALSE)
pileup(bamfile, scanBamParam=scanbam_param, pileupParam=pileup_param)
}

\keyword{methods}
\keyword{manip}
