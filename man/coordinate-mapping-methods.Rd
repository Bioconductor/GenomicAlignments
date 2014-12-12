\name{coordinateMapping-methods}
\alias{coordinateMapping-methods}
\alias{coordinateMapping}

\alias{mapToGenome}
\alias{pmapToGenome}
\alias{pmapToGenome,Ranges,GAlignments-method}
\alias{pmapToGenome,GRanges,GAlignments-method}

\alias{mapToTranscript}
\alias{pmapToTranscript}
\alias{pmapToTranscript,Ranges,GAlignments-method}
\alias{pmapToTranscript,GRanges,GAlignments-method}


\title{Map range coordinates between genome and transcript space}

\description{
  Coordinate mapping methods for \link{GAlignments} objects.
 
  See \code{?`\link[GenomicRanges]{coordinate-mapping-methods}`} in the 
  \pkg{GenomicRanges} package for methods on \link[GenomicRanges]{GRanges} and
  \link[GenomicRanges]{GRangesList} objects.
}

\usage{
\S4method{pmapToGenome}{Ranges,GAlignments}(from, to, ...) 

\S4method{pmapToGenome}{GRanges,GAlignments}(from, to, 
         ignore.strand = TRUE, ...)

\S4method{pmapToTranscript}{Ranges,GAlignments}(from, to, ...) 

\S4method{pmapToTranscript}{GRanges,GAlignments}(from, to, 
         ignore.strand = TRUE, ...)
}

\arguments{
  \item{from}{
    A \code{\linkS4class{Ranges}} or \code{\linkS4class{GRanges}}
    object of positions to be mapped.

    \itemize{
      \item{} For \code{mapToTranscript} \code{from} is in genome-based 
        coordinates 
      \item{} For \code{mapToGenome} \code{from} is in transcript-based 
        coordinates 
    }
  }
  \item{to}{
    A \code{\linkS4class{GAlignments}} object that represents the alignment (aka
    mapping) between the sequences in \code{from} and the sequences in the
    result. Ranges in \code{to} are the 1-based positions of the clipped query
    relative to the reference sequence; this is the standard representation of
    ranges in a GAlignments object.

    The current implementation assumes that positions in \code{from} 
    occur within the \code{to} region, outside of any deletions or insertions.
  }
  \item{ignore.strand}{
    When TRUE, strand is ignored in overlap operations.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  NGS data analysis often requires conversion between genome (reference) and
  transcript-based coordinates. Transcript-based coordinates are used in RNA-Seq
  analysis steps such as mapping spliced reads against the genome,
  quantification, and detection of novel exons. In DNA-Seq, they are used to
  predict the effect of variants detected in the sample. Once the anlaysis is
  done it is often of interest to map these locations back to the linear genome.

  In the reverse direction, reference coordinates for a chromosomal region of
  interest or dbSNP variants need to be mapped to the transcriptome for the
  exploration of coding regions or other gene features of interest.

  On this man page the terms `genomic` and `reference` are used interchangeably;
  `transcript` refers to any transcript-based coordinates such as gene ranges,
  exons, introns, etc.  These functions use a \code{GAlignments} object as the
  `mapper` to represent alignment between the input ranges \code{from} and the
  output. The following CIGAR operations in the "Extended CIGAR format" are 
  used in the coordinate mapping:
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
    \item{\code{mapToGenome}}{
      Range coordinates in transcript space are mapped to genome space.
    }
    \item{\code{mapToTranscript}}{
      Range coordinates in genome space are mapped to transcript space.
    }
    \item{\code{pmapToGenome} and \code{pmapToTranscript}}{
      Parallel versions of \code{mapToGenome} and \code{mapToTranscript}. The 
      i-th range in \code{from} is mapped to the i-th range in \code{to}; 
      \code{from} and \code{to} must have the same length.
    }
  }
}

\value{
  An object the same class as \code{from}.
}

\seealso{
  \itemize{
    \item \code{?`\link[GenomicRanges]{mapping-methods}`} in the 
      in the \pkg{GenomicRanges} package for methods on \link{GenomicRanges} 
      objects.

    \item \url{http://samtools.sourceforge.net/} for a description of
      the Extended CIGAR format.
  }
}

\examples{
## ---------------------------------------------------------------------
## Mapping dbSNP variants from genome space to CDS coordinates
## ---------------------------------------------------------------------

## ---------------------------------------------------------------------
## Mapping cDNA positions to the genome
## ---------------------------------------------------------------------

}

\keyword{methods}
\keyword{utilities}
