\name{pileLettersAt}

\alias{pileLettersAt}


\title{Pile the letters of a set of aligned reads on top of a set
       of individual genomic positions}

\description{
  \code{pileLettersAt} extracts the letters/nucleotides of a set of
  reads that align to a set of individual genomic positions of interest.
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
    A factor-\link[IRanges]{Rle} \emph{parallel} to \code{x}.
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
    A \link[GenomicRanges]{GRanges} object containing the individual genomic
    positions of interest. \code{seqlevels(at)} must be identical to
    \code{levels(seqnames)}.
  }
}

\details{
  \code{x}, \code{seqnames}, \code{pos}, \code{cigar} must be 4 \emph{parallel}
  vectors describing N aligned reads.
}

\value{
  An \link[Biostrings]{XStringSet} (typically \link[Biostrings]{DNAStringSet})
  object \emph{parallel} to \code{at} (i.e. with 1 string per individual
  genomic position).
}

\author{H. Pages}

\seealso{
  \itemize{
    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item \link[GenomicRanges]{GRanges} objects in the \pkg{GenomicRanges}
          package.

    \item The \code{\link{stackStringsFromBam}} function
          for stacking the read sequences (or their quality strings)
          stored in a BAM file on a region of interest.

    \item \link{GAlignments} objects.

    \item \link{cigar-utils} for the CIGAR utility functions used internally
          by \code{pileLettersAt}.

    \item The SAMtools mpileup command available at
          \url{http://samtools.sourceforge.net/} as part of the
          SAMtools project.
  }
}

\examples{
## Input

##   - A BAM file:
bamfile <- BamFile(system.file("extdata", "ex1.bam", package="Rsamtools"))
seqinfo(bamfile)  # to see the seqlevels and seqlengths
stackStringsFromBam(bamfile, param="seq1:1-21")  # a quick look at
                                                 # the reads

##   - A GRanges object containing Individual Genomic Positions Of
##     Interest:
my_IGPOI <- GRanges(Rle(c("seq1", "seq2"), c(7, 2)),
                    IRanges(c(1:5, 21, 1575, 1513:1514), width=1))

## Some preliminary massage on 'my_IGPOI'

seqinfo(my_IGPOI) <- merge(seqinfo(my_IGPOI), seqinfo(bamfile))
seqlevels(my_IGPOI) <- seqlevelsInUse(my_IGPOI)

## Load the BAM file in a GAlignments object. We load only the reads
## aligned to the sequences in 'seqlevels(my_IGPOI)' and we filter out
## reads not passing quality controls (flag bit 0x200) and PCR or
## optical duplicates (flag bit 0x400). See ?ScanBamParam and the SAM
## Spec for more information. 

which <- as(seqinfo(my_IGPOI), "GRanges")
flag <- scanBamFlag(isNotPassingQualityControls=FALSE,
                    isDuplicate=FALSE)
what <- c("seq", "qual")
param <- ScanBamParam(flag=flag, what=c("seq", "qual"), which=which)
gal <- readGAlignmentsFromBam(bamfile, param=param)
seqlevels(gal) <- seqlevels(my_IGPOI) 

## Extract the read sequences (a.k.a. query sequences) and quality
## strings. Both are reported with respect to the + strand.

qseq <- mcols(gal)$seq
qual <- mcols(gal)$qual

nucl_piles <- pileLettersAt(qseq, seqnames(gal), start(gal), cigar(gal),
                            my_IGPOI)
qual_piles <- pileLettersAt(qual, seqnames(gal), start(gal), cigar(gal),
                            my_IGPOI)
mcols(my_IGPOI)$nucl_piles <- nucl_piles
mcols(my_IGPOI)$qual_piles <- qual_piles
my_IGPOI 

## Finally, to summarize A/C/G/T frequencies at each position:
alphabetFrequency(nucl_piles, baseOnly=TRUE)
}

\keyword{methods}
\keyword{manip}
