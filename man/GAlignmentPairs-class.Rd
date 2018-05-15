\name{GAlignmentPairs-class}
\docType{class}

% Class:
\alias{class:GAlignmentPairs}
\alias{GAlignmentPairs-class}

\alias{parallelSlotNames,GAlignmentPairs-method}

% Constructors:
\alias{GAlignmentPairs}

% Accessors:
\alias{strandMode}
\alias{strandMode,GAlignmentPairs-method}
\alias{strandMode<-}
\alias{strandMode<-,GAlignmentPairs-method}
\alias{invertStrand,GAlignmentPairs-method}
\alias{names,GAlignmentPairs-method}
\alias{names<-,GAlignmentPairs-method}
\alias{first}
\alias{first,GAlignmentPairs-method}
\alias{second}
\alias{second,GAlignmentPairs-method}
\alias{last}
\alias{last,GAlignmentPairs-method}
\alias{seqnames,GAlignmentPairs-method}
\alias{strand,GAlignmentPairs-method}
\alias{njunc,GAlignmentPairs-method}
\alias{isProperPair}
\alias{isProperPair,GAlignmentPairs-method}
\alias{seqinfo,GAlignmentPairs-method}
\alias{seqlevelsInUse,GAlignmentPairs-method}
\alias{seqinfo<-,GAlignmentPairs-method}

% List methods:
\alias{[[,GAlignmentPairs,ANY,ANY-method}
\alias{unlist,GAlignmentPairs-method}

% Coercion:
\alias{ranges,GAlignmentPairs-method}
\alias{granges,GAlignmentPairs-method}
\alias{grglist,GAlignmentPairs-method}
\alias{coerce,GAlignmentPairs,IntegerRanges-method}
\alias{coerce,GAlignmentPairs,GRanges-method}
\alias{coerce,GAlignmentPairs,GRangesList-method}

\alias{coerce,GAlignmentPairs,GAlignments-method}

\alias{coerce,GAlignmentPairs,DataFrame-method}
\alias{as.data.frame,GAlignmentPairs-method}

% "show" method:
\alias{show,GAlignmentPairs-method}

% Concatenation:
\alias{bindROWS,GAlignmentPairs-method}

\title{GAlignmentPairs objects}

\description{
  The GAlignmentPairs class is a container for storing \emph{pairs of genomic
  alignments}. These pairs are typically obtained by aligning paired-end reads
  to a reference genome or transcriptome.
}

\details{
  A GAlignmentPairs object is a list-like object where each list element
  represents a pair of genomic alignment.

  An \emph{alignment pair} is made of a "first" and a "last"/"second"
  alignment, and is formally represented by a \link{GAlignments} object of
  length 2. In most applications, an \emph{alignment pair} will represent
  an aligned paired-end read. In that case, the "first" member of the pair
  represents the alignment of the first end of the read (aka "first segment
  in the template", using SAM Spec terminology), and the "last" member of
  the pair represents the alignment of the second end of the read (aka "last
  segment in the template", using SAM Spec terminology).

  In general, a GAlignmentPairs object will be created by loading
  records from a BAM (or SAM) file containing aligned paired-end reads,
  using the \code{readGAlignmentPairs} function (see below).
  Each element in the returned object will be obtained by pairing 2
  records.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GAlignmentPairs(first, last, strandMode=1,
                            isProperPair=TRUE, names=NULL)}:
      Low-level GAlignmentPairs constructor. Generally not used directly.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GAlignmentPairs object.

  \describe{
    \item{}{
      \code{strandMode(x)}, \code{strandMode(x) <- value}:
      The \emph{strand mode} is a per-object switch on GAlignmentPairs objects
      that controls the behavior of the \code{strand} getter.
      More precisely, it indicates how the strand of a pair should be
      inferred from the strand of the first and last alignments in the pair:
      \itemize{
        \item 0: strand of the pair is always *.
        \item 1: strand of the pair is strand of its first alignment.
              This mode should be used when the paired-end data was generated
              using one of the following stranded protocols:
              Directional Illumina (Ligation), Standard SOLiD.
        \item 2: strand of the pair is strand of its last alignment.
              This mode should be used when the paired-end data was generated
              using one of the following stranded protocols:
              dUTP, NSR, NNSR, Illumina stranded TruSeq PE protocol.
      }
      These modes are equivalent to \code{strandSpecific} equal 0, 1, and 2,
      respectively, for the \code{featureCounts} function defined in the
      \pkg{Rsubread} package.

      Note that, by default, the \code{\link{readGAlignmentPairs}} function
      sets the \emph{strand mode} to 1 on the returned GAlignmentPairs object.
      The function has a \code{strandMode} argument to let the user set a
      different \emph{strand mode}. The \emph{strand mode} can also be changed
      any time with the \code{strandMode} setter or with \code{invertStrand}.

      Also note that 3rd party programs TopHat2 and Cufflinks have a
      \code{--library-type} option to let the user specify which protocol
      was used. Please refer to the documentation of these programs for more
      information.
    }
    \item{}{
      \code{length(x)}:
      Return the number of alignment pairs in \code{x}.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or set the names on \code{x}.
      See \code{\link{readGAlignmentPairs}} for how to automatically extract
      and set the names when reading the alignments from a file.
    }
    \item{}{
      \code{first(x, real.strand=FALSE)},
      \code{last(x, real.strand=FALSE)},
      \code{second(x, real.strand=FALSE)}:
      Get the "first" or "last"/"second" alignment for each alignment pair in
      \code{x}. The result is a \link{GAlignments} object of the same length
      as \code{x}.

      If \code{real.strand=TRUE}, then the strand is inverted on-the-fly
      according to the \emph{strand mode} currently set on the object (see
      \code{strandMode(x)} above). More precisely, if \code{strandMode(x)}
      is 0, then the strand is set to * for the \link{GAlignments} object
      returned by both, \code{first()} and \code{last()}.
      If \code{strandMode(x)} is 1, then the strand of the object returned
      by \code{last()} is inverted. If \code{strandMode(x)} is 2, then the
      strand of the object returned by \code{first()} is inverted.
    }
    \item{}{
      \code{seqnames(x)}:
      Get the sequence names of the pairs in \code{x} i.e. the name of the
      reference sequence for each alignment pair in \code{x}. The sequence
      name of a pair is the sequence name of the 2 alignments in the pair
      if they are the same (\emph{concordant seqnames}), or NA if they differ
      (\emph{discordant seqnames}).

      The sequence names are returned in a factor-\link[S4Vectors]{Rle}
      object that is \emph{parallel} to \code{x}, i.e. the i-th element in
      the returned object is the sequence name of the i-th pair in \code{x}.
    }
    \item{}{
      \code{strand(x)}:
      Get the strand for each alignment pair in \code{x}.
      Obey \code{strandMode(x)} above to infer the strand of a pair.
      Return \code{*} for pairs with \emph{discordant strand}, or for
      all pairs if \code{strandMode(x)} is 0.
    }
    \item{}{
      \code{njunc(x)}:
      Equivalent to \code{njunc(first(x)) + njunc(last(x))}.
    }
    \item{}{
      \code{isProperPair(x)}:
      Get the "isProperPair" flag bit (bit 0x2 in SAM Spec) set by
      the aligner for each alignment pair in \code{x}.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Get or set the information about the underlying sequences.
      \code{value} must be a \link[GenomeInfoDb]{Seqinfo} object.
    }
    \item{}{
      \code{seqlevels(x)}, \code{seqlevels(x) <- value}:
      Get or set the sequence levels.
      \code{seqlevels(x)} is equivalent to \code{seqlevels(seqinfo(x))}
      or to \code{levels(seqnames(x))}, those 2 expressions being
      guaranteed to return identical character vectors on a
      GAlignmentPairs object. \code{value} must be a character vector
      with no NAs.
      See \code{?\link{seqlevels}} for more information.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Get or set the sequence lengths.
      \code{seqlengths(x)} is equivalent to \code{seqlengths(seqinfo(x))}.
      \code{value} can be a named non-negative integer or numeric vector
      eventually with NAs.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Get or set the circularity flags.
      \code{isCircular(x)} is equivalent to \code{isCircular(seqinfo(x))}.
      \code{value} must be a named logical vector eventually with NAs.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}:
      Get or set the genome identifier or assembly name for each sequence.
      \code{genome(x)} is equivalent to \code{genome(seqinfo(x))}.
      \code{value} must be a named character vector eventually with NAs.
    }
    \item{}{
      \code{seqnameStyle(x)}:
      Get or set the seqname style for \code{x}.
      Note that this information is not stored in \code{x} but inferred
      by looking up \code{seqnames(x)} against a seqname style database
      stored in the seqnames.db metadata package (required).
      \code{seqnameStyle(x)} is equivalent to \code{seqnameStyle(seqinfo(x))}
      and can return more than 1 seqname style (with a warning)
      in case the style cannot be determined unambiguously.
    }
  }
}

\section{Vector methods}{
  In the code snippets below, \code{x} is a GAlignmentPairs object.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new GAlignmentPairs object made of the selected
      alignment pairs. 
    }
  }
}

\section{List methods}{
  In the code snippets below, \code{x} is a GAlignmentPairs object.

  \describe{
    \item{}{
      \code{x[[i]]}:
      Extract the i-th alignment pair as a \link{GAlignments} object
      of length 2. As expected \code{x[[i]][1]} and \code{x[[i]][2]} are
      respectively the "first" and "last" alignments in the pair.
    }
    \item{}{
      \code{unlist(x, use.names=TRUE)}:
      Return the \link{GAlignments} object conceptually defined
      by \code{c(x[[1]], x[[2]], ..., x[[length(x)]])}.
      \code{use.names} determines whether \code{x} names should be
      propagated to the result or not.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GAlignmentPairs object.

  \describe{
    \item{}{
      \code{granges(x, use.names=TRUE, use.mcols=FALSE,
                       on.discordant.seqnames=c("error", "drop", "split"))},
      \code{ranges(x, use.names=TRUE, use.mcols=FALSE,
                      on.discordant.seqnames=c("error", "drop", "split"))}:

      Return a \link{GRanges} object (for \code{granges()}) or
      \link[IRanges]{IRanges}) object (for \code{ranges()}).

      If \code{x} contains no pairs with \emph{discordant seqnames}, the
      operation is guaranteed to be successful and to return an object
      \emph{parallel} to \code{x}, that is, an object where the i-th element
      is the range of the genomic region spanned by the i-th alignment in
      \code{x} (all gaps in the region are ignored).

      If \code{x} contains pairs with discordant seqnames, then an error
      is raised, unless the \code{on.discordant.seqnames} argument is set
      to \code{"drop"} or \code{"split"}, in which case the pairs with
      discordant seqnames are either dropped or represented with 2 genomic
      ranges (or 2 ranges for \code{ranges()}) in the returned object.
      In that case, the returned object is NOT \emph{parallel} to \code{x}.

      If \code{use.names} is TRUE, then the names on \code{x} 
      (if any) are propagated to the returned object.
      If \code{use.mcols} is TRUE, then the metadata columns on \code{x}
      (if any) are propagated to the returned object.
    }
    \item{}{
      \code{grglist(x, use.mcols=FALSE, drop.D.ranges=FALSE)}:

      Return a \link{GRangesList} object of length \code{length(x)}
      where the i-th element represents the ranges (with respect to the
      reference) of the i-th alignment pair in \code{x}. The strand of
      the returned ranges obeys the \emph{strand mode} currently set on the
      object (see \code{strandMode(x)} above).

      More precisely, if \code{grl1} and \code{grl2} are
      \code{grglist(first(x, real.strand=TRUE), order.as.in.query=TRUE)} and
      \code{grglist(last(x, real.strand=TRUE), order.as.in.query=TRUE)},
      respectively, then the i-th element in the returned \link{GRangesList}
      object is \code{c(grl1[[i]], grl2[[i]])} if \code{strandMode(x)} is 0
      or 1, or \code{c(grl2[[i]], grl1[[i]])} if \code{strandMode(x)} is 2.

      Note that, if \code{strandMode(x)} is 1 or 2, this results in the ranges
      being in consistent order with the original "query template", that is,
      being in the order defined by walking the "query template" from the
      beginning to the end.

      If \code{use.names} is TRUE, then the names on \code{x} 
      (if any) are propagated to the returned object.
      If \code{use.mcols} is TRUE, then the metadata columns on \code{x}
      (if any) are propagated to the returned object.

      If \code{drop.D.ranges} is TRUE, then deletions (Ds in the
      CIGAR) are treated like junctions (Ns in the CIGAR), that is, the
      ranges corresponding to deletions are dropped.
    }
    \item{}{
      \code{as(x, "GRanges")}, \code{as(x, "IntegerRanges")},
      \code{as(x, "GRangesList")}:
      Alternate ways of doing
      \code{granges(x, use.names=TRUE, use.mcols=TRUE)},
      \code{ranges(x, use.names=TRUE, use.mcols=TRUE)}, and
      \code{grglist(x, use.names=TRUE, use.mcols=TRUE)}, respectively.
    }
    \item{}{
      \code{as(x, "GAlignments")}:
      Equivalent of \code{unlist(x, use.names=TRUE)}.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is a GAlignmentPairs object.

  \describe{
    \item{}{
      \code{show(x)}:
      By default the \code{show} method displays 5 head and 5 tail
      elements. This can be changed by setting the global options
      \code{showHeadLines} and \code{showTailLines}. If the object
      length is less than (or equal to) the sum of these 2 options
      plus 1, then the full object is displayed.
      Note that these options also affect the display of \link{GRanges}
      and \link{GAlignments} objects, as well as other \link[S4Vectors]{Vector}
      derivatives defined in the \pkg{IRanges} and \pkg{Biostrings} packages
      (e.g. \link[IRanges]{IRanges} and \link[Biostrings]{XStringSet} objects).
    }
  }
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \code{\link{readGAlignmentPairs}} for reading aligned paired-end
          reads from a file (typically a BAM file) into a GAlignmentPairs
          object.

    \item \link{GAlignments} objects for handling aligned single-end reads.

    \item \code{\link{makeGAlignmentPairs}} for pairing the elements of a
          \link{GAlignments} object into a GAlignmentPairs object.

    \item \link{junctions-methods} for extracting and summarizing junctions
          from a GAlignmentPairs object.

    \item \link[GenomicAlignments]{coverage-methods} for computing the
          coverage of a GAlignmentPairs object.

    \item \link[GenomicAlignments]{findOverlaps-methods} for finding range
          overlaps between a GAlignmentPairs object and another range-based
          object.

    \item \code{\link[GenomeInfoDb]{seqinfo}} in the \pkg{GenomeInfoDb}
          package for getting/setting/modifying the sequence information
          stored in an object.

    \item The \link[GenomicRanges]{GRanges} and
          \link[GenomicRanges]{GRangesList} classes defined and documented
          in the \pkg{GenomicRanges} package.
  }
}

\examples{
library(Rsamtools)  # for the ex1.bam file
ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galp <- readGAlignmentPairs(ex1_file, use.names=TRUE, strandMode=1)
galp

length(galp)
head(galp)
head(names(galp))

first(galp)
last(galp)
# or
second(galp)

strandMode(galp)
first(galp, real.strand=TRUE)
last(galp, real.strand=TRUE)
strand(galp)

strandMode(galp) <- 2
first(galp, real.strand=TRUE)
last(galp, real.strand=TRUE)
strand(galp)

seqnames(galp)

head(njunc(galp))
table(isProperPair(galp))
seqlevels(galp)

## Rename the reference sequences:
seqlevels(galp) <- sub("seq", "chr", seqlevels(galp))
seqlevels(galp)

galp[[1]]
unlist(galp)

grglist(galp)  # a GRangesList object

strandMode(galp) <- 1
grglist(galp)

## Alternatively the strand mode can be changed with invertStrand():
invertStrand(galp)

stopifnot(identical(unname(elementNROWS(grglist(galp))), njunc(galp) + 2L))

granges(galp)  # a GRanges object
}

\keyword{methods}
\keyword{classes}
