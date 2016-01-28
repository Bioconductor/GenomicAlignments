\name{mapToAlignments}
\alias{coordinate-mapping-methods}

\alias{mapToAlignments}
\alias{mapToAlignments,Ranges,GAlignments-method}
\alias{mapToAlignments,GenomicRanges,GAlignments-method}
\alias{pmapToAlignments}
\alias{pmapToAlignments,Ranges,GAlignments-method}
\alias{pmapToAlignments,GenomicRanges,GAlignments-method}

\alias{mapFromAlignments}
\alias{mapFromAlignments,Ranges,GAlignments-method}
\alias{mapFromAlignments,GenomicRanges,GAlignments-method}
\alias{pmapFromAlignments}
\alias{pmapFromAlignments,Ranges,GAlignments-method}
\alias{pmapFromAlignments,GenomicRanges,GAlignments-method}

\title{Map range coordinates between reads and genome space using
       CIGAR alignments}

\description{
  Map range coordinates between reads (local) and genome (reference) space 
  using the CIGAR in a \code{GAlignments} object.
 
  See \code{?\link[GenomicFeatures]{mapToTranscripts}} in the 
  \pkg{GenomicRanges} package for mapping coordinates between features
  in the transcriptome and genome space.
}

\usage{
\S4method{mapToAlignments}{GenomicRanges,GAlignments}(x, alignments, ...) 
\S4method{pmapToAlignments}{GenomicRanges,GAlignments}(x, alignments, ...) 

\S4method{mapFromAlignments}{GenomicRanges,GAlignments}(x, alignments, ...) 
\S4method{pmapFromAlignments}{GenomicRanges,GAlignments}(x, alignments, ...) 
}

\arguments{
  \item{x}{
    \code{\linkS4class{GenomicRanges}} object of positions to be mapped.
    \code{x} must have names when mapping to the genome.
  }
  \item{alignments}{
    A \code{\linkS4class{GAlignments}} object that represents the alignment of
    \code{x} to the genome. The \code{aligments} object must have names. When
    mapping to the genome names are used to determine mapping pairs and in the
    reverse direction they are used as the seqlevels of the output object.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  These methods use a \code{GAlignments} object to represent the alignment 
  between the ranges in \code{x} and the output. The following CIGAR 
  operations in the "Extended CIGAR format" are used in the mapping
  algorithm:
  \itemize{
    \item{M, X, =} Sequence match or mismatch
    \item{I} Insertion to the reference 
    \item{D} Deletion from the reference 
    \item{N} Skipped region from the reference
    \item{S} Soft clip on the read
    \item{H} Hard clip on the read
    \item{P} Silent deletion from the padded reference
  }
  \itemize{
    \item{\code{mapToAlignments}, \code{pmapToAlignments}}{
      The CIGAR is used to map the genomic (reference) position \code{x} to 
      local coordinates. The mapped position starts at 
      \preformatted{
      start(x) - start(alignments) + 1
      }
      and is incremented or decremented as the algorithm walks the length of 
      the CIGAR. A successful mapping in this direction requires that
      \code{x} fall within \code{alignments}.

       The seqlevels of the return object are taken from the 
       \code{alignments} object and will be a name descriptive of the read
       or aligned region. In this direction, mapping is attempted between all 
       elements of \code{x} and all elements of \code{alignments}.
    }
    \item{\code{mapFromAlignments}, \code{pmapFromAlignments}}{
      The CIGAR is used to map the local position \code{x} to genomic
      (reference) coordinates. The mapped position starts at
      \preformatted{
      start(x) + start(alignments) - 1
      } 
      and is incremented or decremented as the algorithm walks the length of 
      the CIGAR. A successful mapping in this direction requires that the
      width of \code{alignments} is <= the width of \code{x}.

      When mapping to the genome, name matching is used to determine the
      mapping pairs (vs attempting to match all possible pairs). Ranges in
      \code{x} are only mapped to ranges in \code{alignments} with the
      same name. Name matching is motivated by use cases such as
      differentially expressed regions where the expressed regions in
      \code{x} would only be related to a subset of regions in
      \code{alignments}, which may contains gene or transcript ranges.
    }
    \item{element-wise versions}{
      \code{pmapToAlignments} and \code{pmapFromAlignments} are element-wise
      (aka `parallel`) versions of \code{mapToAlignments} and 
      \code{mapFromAlignments}. The i-th range in \code{x} is mapped to the 
      i-th range in \code{alignments}; \code{x} and \code{alignments} must 
      have the same length.

      Ranges in \code{x} that do not map (out of bounds) are returned as
      zero-width ranges starting at 0. These ranges are given the special
      seqname of "UNMAPPED". Note the non-parallel methods do not return
      unmapped ranges so the "UNMAPPED" seqname is unique to
      \code{pmapToAlignments} and \code{pmapFromAlignments}.
    }
    \item{strand}{
      By SAM convention, the CIGAR string is reported for mapped reads on the 
      forward genomic strand. There is no need to consider strand in these
      methods. The output of these methods will always be unstranded 
      (i.e., "*").
    }
  }
}

\value{
  An object the same class as \code{x}.

  Parallel methods return an object the same shape as \code{x}. Ranges that
  cannot be mapped (out of bounds) are returned as zero-width ranges starting 
  at 0 with a seqname of "UNMAPPED".

  Non-parallel methods return an object that varies in length similar to a
  Hits object. The result only contains mapped records, out of bound ranges 
  are not returned. \code{xHits} and \code{alignmentsHits} metadata columns 
  indicate the elements of \code{x} and \code{alignments} used in the mapping.

  When present, names from \code{x} are propagated to the output. When
  mapping locally, the seqlevels of the output are the names on the 
  \code{alignment} object. When mapping globally, the output seqlevels are 
  the seqlevels of \code{alignment} which are usually chromosome names.
}

\seealso{
  \itemize{
    \item \code{?\link[GenomicFeatures]{mapToTranscripts}} in the 
      in the \pkg{GenomicFeatures} package for methods mapping between
      transcriptome and genome space.

    \item \url{http://samtools.sourceforge.net/} for a description of
      the Extended CIGAR format.
  }
}

\author{V. Obenchain, M. Lawrence and H. Pagès}

\examples{
## ---------------------------------------------------------------------
## A. Basic use 
## ---------------------------------------------------------------------

## 1. Map to local space with mapToAlignments()
## ---------------------------------------------------------------------

## Mapping to local coordinates requires 'x' to be within 'alignments'.
## In this 'x', the second range is too long and can't be mapped.
alignments <- GAlignments("chr1", 10L, "11M", strand("*"), names="read_A")
x <- GRanges("chr1", IRanges(c(12, 12), width=c(6, 20)))
mapToAlignments(x, alignments)

## The element-wise version of the function returns unmapped ranges
## as zero-width ranges with a seqlevel of "UNMAPPED":
pmapToAlignments(x, c(alignments, alignments))

## Mapping the same range through different alignments demonstrates 
## how the CIGAR operations affect the outcome.
ops <- c("no-op", "junction", "insertion", "deletion")
x <- GRanges(rep("chr1", 4), IRanges(rep(12, 4), width=rep(6, 4), names=ops)) 
alignments <- GAlignments(rep("chr1", 4), rep(10L, 4), 
                          cigar = c("11M", "5M2N4M", "5M2I4M", "5M2D4M"),
                          strand = strand(rep("*", 4)),
                          names = paste0("region_", 1:4))
pmapToAlignments(x, alignments)

## 2. Map to genome space with mapFromAlignments()
## ---------------------------------------------------------------------

## One of the criteria when mapping to genomic coordinates is that the
## shifted 'x' range falls within 'alignments'. Here the first 'x' 
## range has a shifted start value of 14 (5 + 10 - 1 = 14) with a width of 
## 2 and so is successfully mapped. The second has a shifted start of 29
## (20 + 10 - 1 = 29) which is outside the range of 'alignments'.
x <- GRanges("chr1", IRanges(c(5, 20), width=2, names=rep("region_A", 2)))
alignments <- GAlignments("chr1", 10L, "11M", strand("*"), names="region_A")
mapFromAlignments(x, alignments)

## Another characteristic of mapping this direction is the name matching
## used to determine pairs. Mapping is only attempted between ranges in 'x' 
## and 'alignments' with the same name. If we change the name of the first 'x' 
## range, only the second will be mapped to 'alignment'. We know the second
## range fails to map so we get an empty result.
names(x) <- c("region_B", "region_A")
mapFromAlignments(x, alignments)

## CIGAR operations: insertions reduce the width of the output while
## junctions and deletions increase it.
ops <- c("no-op", "junction", "insertion", "deletion")
x <- GRanges(rep("chr1", 4), IRanges(rep(3, 4), width=rep(5, 4), names=ops)) 
alignments <- GAlignments(rep("chr1", 4), rep(10L, 4), 
                         cigar = c("11M", "5M2N4M", "5M2I4M", "5M2D4M"),
                         strand = strand(rep("*", 4)))
pmapFromAlignments(x, alignments)

## ---------------------------------------------------------------------
## B. TATA box motif: mapping from read -> genome -> transcript
## ---------------------------------------------------------------------

## The TATA box motif is a conserved DNA sequence in the core promoter
## region. Many eukaryotic genes have a TATA box located approximately
## 25-35 base pairs upstream of the transcription start site. The motif is 
## the binding site of general transcription factors or histones and
## plays a key role in transcription.

## In this example, the position of the TATA box motif (if present) is 
## located in the DNA sequence corresponding to read ranges. The local 
## motif positions are mapped to genome coordinates and then mapped
## to gene features such as promoters regions.

## Load reads from chromosome 4 of D. melanogaster (dm3):
library(pasillaBamSubset)
fl <- untreated1_chr4()
gal <- readGAlignments(fl)

## Extract DNA sequences corresponding to the read ranges:
library(GenomicFeatures)
library(BSgenome.Dmelanogaster.UCSC.dm3)
dna <- extractTranscriptSeqs(BSgenome.Dmelanogaster.UCSC.dm3, grglist(gal))

## Search for the consensus motif TATAAA in the sequences:
box <- vmatchPattern("TATAAA", dna)

## Some sequences had more than one match:
table(elementNROWS(box))

## The element-wise function we'll use for mapping to genome coordinates
## requires the two input argument to have the same length. We need to
## replicate the read ranges to match the number of motifs found.

## Expand the read ranges to match motifs found:
motif <- elementNROWS(box) != 0
alignments <- rep(gal[motif], elementNROWS(box)[motif])

## We make the IRanges into a GRanges object so the seqlevels can
## propagate to the output. Seqlevels are needed in the last mapping step.
readCoords <- GRanges(seqnames(alignments), unlist(box, use.names=FALSE))

## Map the local position of the motif to genome coordinates:
genomeCoords <- pmapFromAlignments(readCoords, alignments) 
genomeCoords

## We are interested in the location of the TATA box motifs in the
## promoter regions. To perform the mapping we need the promoter ranges 
## as a GRanges or GRangesList.

## Extract promoter regions 50 bp upstream from the transcription start site:
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
promoters <- promoters(txdb, upstream=50, downstream=0)

## Map the genome coordinates to the promoters:
names(promoters) <- mcols(promoters)$tx_name  ## must be named 
mapToTranscripts(genomeCoords, promoters) 









}

\keyword{methods}
\keyword{utilities}
