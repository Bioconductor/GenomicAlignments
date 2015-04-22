\name{GAlignmentPairs-class}
\docType{class}

% Class:
\alias{class:GAlignmentPairs}
\alias{GAlignmentPairs-class}

% Constructors:
\alias{GAlignmentPairs}

% Accessors:
\alias{length,GAlignmentPairs-method}
\alias{names,GAlignmentPairs-method}
\alias{names<-,GAlignmentPairs-method}
\alias{first}
\alias{first,GAlignmentPairs-method}
\alias{last}
\alias{last,GAlignmentPairs-method}
\alias{left}
\alias{left,GAlignmentPairs-method}
\alias{right}
\alias{right,GAlignmentPairs-method}
\alias{seqnames,GAlignmentPairs-method}
\alias{strandMode}
\alias{strandMode,GAlignmentPairs-method}
\alias{strandMode<-}
\alias{strandMode<-,GAlignmentPairs-method}
\alias{strand,GAlignmentPairs-method}
\alias{strand<-,GAlignmentPairs-method}
\alias{njunc,GAlignmentPairs-method}
\alias{isProperPair}
\alias{isProperPair,GAlignmentPairs-method}
\alias{elementMetadata<-,GAlignmentPairs-method}
\alias{seqinfo,GAlignmentPairs-method}
\alias{seqlevelsInUse,GAlignmentPairs-method}
\alias{seqinfo<-,GAlignmentPairs-method}

% List methods:
\alias{[[,GAlignmentPairs,ANY,ANY-method}
\alias{unlist,GAlignmentPairs-method}

% Coercion:
\alias{grglist,GAlignmentPairs-method}
\alias{granges,GAlignmentPairs-method}
\alias{coerce,GAlignmentPairs,GRangesList-method}
\alias{coerce,GAlignmentPairs,GRanges-method}
\alias{coerce,GAlignmentPairs,GAlignments-method}

% "show" method:
\alias{show,GAlignmentPairs-method}

% Combining:
\alias{c,GAlignmentPairs-method}


\title{GAlignmentPairs objects}

\description{
  The GAlignmentPairs class is a container for "genomic alignment pairs".

  IMPORTANT NOTE: The GAlignmentPairs container can only hold \emph{concordant}
  pairs at the moment, that is, pairs where the 2 alignments are on the
  opposite strands of the same reference sequence.
}

\details{
  A GAlignmentPairs object is a list-like object where each element
  describes a pair of genomic alignment.

  An "alignment pair" is made of a "first" and a "last" alignment,
  and is formally represented by a \link{GAlignments} object of
  length 2. It is typically representing a hit of a paired-end read to
  the reference genome that was used by the aligner. More precisely,
  in a given pair, the "first" alignment represents the hit of the first
  end of the read (aka "first segment in the template", using SAM Spec
  terminology), and the "last" alignment represents the hit of the second
  end of the read (aka "last segment in the template", using SAM Spec
  terminology).

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
      \code{first(x, invert.strand=FALSE)},
      \code{last(x, invert.strand=FALSE)}:
      Get the "first" or "last" alignment for each alignment pair in
      \code{x}.
      The result is a \link{GAlignments} object of the same length
      as \code{x}.
      If \code{invert.strand=TRUE}, then the strand is inverted on-the-fly,
      i.e. "+" becomes "-", "-" becomes "+", and "*" remains unchanged.
    }
    \item{}{
      \code{left(x)}:
      Get the "left" alignment for each alignment pair in \code{x}.
      By definition, the "left" alignment in a pair is the alignment that
      is on the + strand. If this is the "first" alignment, then it's returned
      as-is by \code{left(x)}, but if this is the "last" alignment, then it's
      returned by \code{left(x)} with the strand inverted.
    }
    \item{}{
      \code{right(x)}:
      Get the "right" alignment for each alignment pair in \code{x}.
      By definition, the "right" alignment in a pair is the alignment that
      is on the - strand. If this is the "first" alignment, then it's returned
      as-is by \code{right(x)}, but if this is the "last" alignment, then it's
      returned by \code{right(x)} with the strand inverted.
    }
    \item{}{
      \code{seqnames(x)}:
      Get the name of the reference sequence for each alignment pair
      in \code{x}. When reading the alignments from a BAM file, this comes
      from the RNAME field which has the same value for the 2 records in a
      pair (\code{\link{readGAlignmentPairs}}, the function used for reading
      paired-end reads from a BAM file as a GAlignmentPairs object, rejects
      pairs with incompatible RNAME values).
    }
    \item{}{
      \code{strandMode(x)}, \code{strandMode(x) <- value}:
      The strand mode indicates how to infer the strand of a pair from the
      strand of the first and last alignments in the pair:
      \itemize{
        \item 0: strand of the pair is always *.
        \item 1: strand of the pair is strand of its first alignment.
              This mode should be used when the paired-end data was obtained
              using one of the following stranded protocols:
              Directional Illumina (Ligation), Standard SOLiD.
        \item 2: strand of the pair is strand of its last alignment.
              This mode should be used when the paired-end data was obtained
              using one of the following stranded protocols:
              dUTP, NSR, NNSR, Illumina stranded TruSeq PE protocol.
      }
      These modes are equivalent to \code{strandSpecific} equal 0, 1, and 2,
      respectively, for the \code{featureCounts} function defined in the
      \pkg{Rsubread} package.

      Also note that 3rd party software TopHat2 and Cufflinks have a
      \code{--library-type} option to let the user specify which protocol
      was used. Please refer to their documentation for more information.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Get or set the strand for each alignment pair in \code{x}.
      Obeys \code{strandMode(x)} above to infer the strand of a pair.
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
      \code{grglist(x, use.mcols=FALSE,
                       order.as.in.query=FALSE,
                       drop.D.ranges=FALSE)}:

      Return a \link{GRangesList} object of length \code{length(x)}
      where the i-th element represents the ranges (with respect to the
      reference) of the i-th alignment pair in \code{x}.

      If \code{use.mcols} is TRUE and \code{x} has metadata columns on it
      (accessible with \code{mcols(x)}), they're propagated to the returned
      object.

      IMPORTANT: If \code{strandMode(x)} is 1, the strand of the ranges
      coming from the "first" alignment in a pair is preserved and the
      strand of those coming from the "last" alignment in a pair is
      inverted. If \code{strandMode(x)} is 2, it's the other way around.

      The \code{order.as.in.query} toggle affects the order of the ranges
      \emph{within} each top-level element of the returned object.

      If \code{FALSE} (the default), then the "left" ranges are placed before
      the "right" ranges, and, within each left or right group, are ordered
      from 5' to 3' in elements associated with the plus strand and from 3'
      to 5' in elements associated with the minus strand.
      More formally, the i-th element in the returned \link{GRangesList}
      object can be defined as \code{c(grl1[[i]], grl2[[i]])}, where
      \code{grl1} is \code{grglist(left(x))} and \code{grl2} is
      \code{grglist(right(x))}.

      If \code{TRUE}, then the "first" ranges are placed before the "last"
      ranges, and, within each first or last group, are \emph{always}
      ordered from 5' to 3', whatever the strand is.
      More formally, the i-th element in the returned \link{GRangesList}
      object can be defined as \code{c(grl1[[i]], grl2[[i]])}, where
      \code{grl1} is \code{grglist(first(x),
                                   order.as.in.query=TRUE)}
      and
      \code{grl2} is \code{grglist(last(x, invert.strand=TRUE),
                                   order.as.in.query=TRUE)}.

      Note that the relationship between the 2 \link{GRangesList} objects
      obtained with \code{order.as.in.query} being respectively
      \code{FALSE} or \code{TRUE} is simpler than it sounds: the only
      difference is that the order of the ranges in elements associated
      with the \emph{minus} strand is reversed.

      Finally note that, in the latter, the ranges are \emph{always} ordered
      consistently with the original "query template", that is, in the order
      defined by walking the "query template" from the beginning to the end.

      If \code{drop.D.ranges} is \code{TRUE}, then deletions (Ds in the
      CIGAR) are treated like junctions (Ns in the CIGAR), that is, the
      ranges corresponding to deletions are dropped.
    }
    \item{}{
      \code{granges(x, use.mcols=FALSE)}: Return a \link{GRanges} object
      of length \code{length(x)} where each range is obtained by merging
      all the ranges within the corresponding top-level element in
      \code{grglist(x)}.

      If \code{use.mcols} is TRUE and \code{x} has metadata columns on it
      (accessible with \code{mcols(x)}), they're propagated to the returned
      object.
    }
    \item{}{
      \code{as(x, "GRangesList")}, \code{as(x, "GRanges")}:
      Alternate ways of doing \code{grglist(x, use.mcols=TRUE)} and
      \code{granges(x, use.mcols=TRUE)}, respectively.
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
      and \link{GAlignments} objects, as well as other objects defined
      in the IRanges and Biostrings packages (e.g. \link[IRanges]{Ranges}
      and \link[Biostrings]{XStringSet} objects).
    }
  }
}

\author{
  H. Pages
}

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
last(galp, invert.strand=TRUE)
left(galp)
right(galp)
seqnames(galp)

strandMode(galp)
strand(galp)

strandMode(galp) <- 2
strand(galp)

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

grglist(galp, order.as.in.query=TRUE)
stopifnot(identical(unname(elementLengths(grglist(galp))), njunc(galp) + 2L))

granges(galp)  # a GRanges object
}

\keyword{methods}
\keyword{classes}
