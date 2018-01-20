\name{GAlignmentsList-class}
\docType{class}

% Class
\alias{class:GAlignmentsList}
\alias{GAlignmentsList-class}
\alias{GAlignmentsList}

% Constructors:
\alias{GAlignmentsList}

% Accessors:
\alias{names,GAlignmentsList-method}
\alias{names<-,GAlignmentsList-method}
\alias{seqnames,GAlignmentsList-method}
\alias{seqnames<-,GAlignmentsList-method}
\alias{rname,GAlignmentsList-method}
\alias{rname<-,GAlignmentsList-method}
\alias{strand,GAlignmentsList-method}
\alias{strand<-,GAlignmentsList-method}
\alias{strand<-,GAlignmentsList,character-method}
\alias{cigar,GAlignmentsList-method}
\alias{qwidth,GAlignmentsList-method}
\alias{njunc,GAlignmentsList-method}
\alias{elementMetadata,GAlignmentsList-method}
\alias{elementMetadata<-,GAlignmentsList-method}
\alias{seqinfo,GAlignmentsList-method}
\alias{seqinfo<-,GAlignmentsList-method}
\alias{start,GAlignmentsList-method}
\alias{end,GAlignmentsList-method}
\alias{width,GAlignmentsList-method}

% Coercion:
\alias{ranges,GAlignmentsList-method}
\alias{granges,GAlignmentsList-method}
\alias{grglist,GAlignmentsList-method}
\alias{rglist,GAlignmentsList-method}
\alias{coerce,GAlignmentsList,IntegerRanges-method}
\alias{coerce,GAlignmentsList,GRanges-method}
\alias{coerce,GAlignmentsList,GRangesList-method}
\alias{coerce,GAlignmentsList,RangesList-method}

\alias{coerce,GAlignmentPairs,GAlignmentsList-method}
\alias{coerce,GAlignmentsList,GAlignmentPairs-method}

\alias{coerce,list,GAlignmentsList-method}

% extractList() and family:
\alias{relistToClass,GAlignments-method}

% show:
\alias{show,GAlignmentsList-method}


\title{GAlignmentsList objects}

\description{
  The GAlignmentsList class is a container for storing a collection of
  \link{GAlignments} objects.
}

\details{
  A GAlignmentsList object contains a list of \link{GAlignments} objects.
  The majority of operations on this page are described in more detail
  on the GAlignments man page, see ?\code{GAlignments}.
} 

\section{Constructor}{
  \describe{
    \item{}{
      \code{GAlignmentsList(...)}:
      Creates a GAlignmentsList from a list of \link{GAlignments} objects.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GAlignmentsList object.

  \describe{
    \item{}{
      \code{length(x)}:
      Return the number of elements in \code{x}.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or set the names of the elements of \code{x}.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Get or set the name of the reference sequences of the
      alignments in each element of \code{x}.
    }
    \item{}{
      \code{rname(x)}, \code{rname(x) <- value}:
      Same as \code{seqnames(x)} and \code{seqnames(x) <- value}.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Get or set the strand of the alignments in each element 
      of \code{x}.
    }
    \item{}{
      \code{cigar(x)}:
      Returns a character list of length \code{length(x)}
      containing the CIGAR string for the alignments in
      each element of \code{x}.
    }
    \item{}{
      \code{qwidth(x)}:
      Returns an integer list of length \code{length(x)}
      containing the length of the alignments in each element of
      \code{x} *after* hard clipping (i.e. the length of the 
      query sequence that is stored in the corresponding SAM/BAM record).
    }
    \item{}{
      \code{start(x)}, \code{end(x)}:
      Returns an integer list of length \code{length(x)}
      containing the "start" and "end" (respectively) of the 
      alignments in each element of \code{x}. 
    }
    \item{}{
      \code{width(x)}:
      Returns an integer list of length \code{length(x)} containing
      the "width" of the alignments in each element of \code{x}.
    }
    \item{}{
      \code{njunc(x)}:
      Returns an integer list of length \code{x} containing the number
      of junctions (i.e. N operations in the CIGAR) for the alignments
      in each element of \code{x}.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Get or set the information about the underlying sequences.
      \code{value} must be a \link[GenomeInfoDb]{Seqinfo} object.
    }
    \item{}{
      \code{seqlevels(x)}, \code{seqlevels(x) <- value}:
      Get or set the sequence levels of the alignments in each element
      of \code{x}.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Get or set the sequence lengths for each element of \code{x}.
      \code{seqlengths(x)} is equivalent to \code{seqlengths(seqinfo(x))}.
      \code{value} can be a named non-negative integer or numeric vector
      eventually with NAs.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Get or set the circularity flags for the alignments in each
      element in \code{x}. \code{value} must be a named logical list 
      eventually with NAs.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}:
      Get or set the genome identifier or assembly name for the alignments 
      in each element of \code{x}. \code{value} must be a named character 
      list eventually with NAs.
    }
    \item{}{
      \code{seqnameStyle(x)}:
      Get or set the seqname style for alignments in each element of \code{x}.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GAlignmentsList object.

  \describe{
    \item{}{
      \code{granges(x, use.names=TRUE, use.mcols=FALSE, ignore.strand=FALSE)}, 
      \code{ranges(x, use.names=TRUE, use.mcols=FALSE)}:

      Return either a \link{GRanges} or a \link[IRanges]{IRanges}
      object of length \code{length(x)}. Note this coercion IGNORES 
      the cigar information. The resulting ranges span the entire
      range, including any junctions or spaces between paired-end reads.

      If \code{use.names} is TRUE, then the names on \code{x}
      (if any) are propagated to the returned object.
      If \code{use.mcols} is TRUE, then the metadata columns on \code{x}
      (if any) are propagated to the returned object.

      \code{granges} coercion supports \code{ignore.strand} to allow 
      ranges of opposite strand to be combined (see examples). All
      ranges in the resulting GRanges will have strand \sQuote{*}.
    }
    \item{}{
      \code{grglist(x, use.names=TRUE, use.mcols=FALSE, ignore.strand=FALSE)}, 
      \code{rglist(x, use.names=TRUE, use.mcols=FALSE)}:

      Return either a \link{GRangesList} or an \link[IRanges]{IRangesList}
      object of length \code{length(x)}. This coercion RESPECTS the cigar 
      information. The resulting ranges are fragments of the original ranges 
      that do not include junctions or spaces between paired-end reads.
 
      If \code{use.names} is TRUE, then the names on \code{x}
      (if any) are propagated to the returned object.
      If \code{use.mcols} is TRUE, then the metadata columns on \code{x}
      (if any) are propagated to the returned object.

      \code{grglist} coercion supports \code{ignore.strand} to allow 
      ranges of opposite strand to be combined (see examples). When
      \code{ignore.strand} is TRUE all ranges in the resulting 
      GRangesList have strand \sQuote{*}.
    }
    \item{}{
      \code{as(x, "GRanges")}, \code{as(x, "IntegerRanges")},
      \code{as(x, "GRangesList")}, \code{as(x, "RangesList")}:
      Alternate ways of doing
      \code{granges(x, use.names=TRUE, use.mcols=TRUE)},
      \code{ranges(x, use.names=TRUE, use.mcols=TRUE)},
      \code{grglist(x, use.names=TRUE, use.mcols=TRUE)}, and
      \code{rglist(x, use.names=TRUE, use.mcols=TRUE)}, respectively.
    }
    \item{}{
      \code{as.data.frame(x, row.names = NULL, optional = FALSE,
         ..., value.name = "value", use.outer.mcols = FALSE, 
         group_name.as.factor = FALSE)}: 
      Coerces \code{x} to a \code{data.frame}. See as.data.frame on the 
      \code{List} man page for details (?\code{List}). 
    }
    \item{}{
      \code{as(x, "GALignmentsList")}: Here \code{x} is a
      \link{GAlignmentPairs} object. Return a GAlignmentsList object of length
      \code{length(x)} where the i-th list element represents the ranges of
      the i-th alignment pair in \code{x}.
    }
  }
}

\section{Subsetting and related operations}{
  In the code snippets below, \code{x} is a GAlignmentsList object.

  \describe{
    \item{}{
      \code{x[i]}, \code{x[i] <- value}:
      Get or set list elements \code{i}. \code{i} can be a numeric 
      or logical vector. \code{value} must be a GAlignments.
    }
    \item{}{
      \code{x[[i]]}, \code{x[[i]] <- value}:
      Same as \code{x[i]}, \code{x[i] <- value}. 
    }
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Get or set list elements \code{i} with optional metadata columns
      \code{j}. \code{i} can be a numeric, logical or missing. 
      \code{value} must be a GAlignments.
    }
  }
}

\section{Concatenation}{
  \describe{
    \item{}{
      \code{c(x, ..., ignore.mcols=FALSE)}:
      Concatenate GAlignmentsList object \code{x} and the GAlignmentsList
      objects in \code{...} together.
      See \code{?\link[S4Vectors]{c}} in the \pkg{S4Vectors} package for
      more information about concatenating Vector derivatives.
    }
  }
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{Valerie Obenchain <Valerie.Obenchain@RoswellPark.org>}

\seealso{
  \itemize{
    \item \code{\link{readGAlignmentsList}} for reading genomic alignments
          from a file (typically a BAM file) into a GAlignmentsList object.

    \item \link{GAlignments} and \link{GAlignmentPairs} objects for handling
          aligned single- and paired-end reads, respectively.

    \item \link{junctions-methods} for extracting and summarizing junctions
          from a GAlignmentsList object.

    \item \link[GenomicAlignments]{findOverlaps-methods} for finding range
          overlaps between a GAlignmentsList object and another range-based
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
gal1 <- GAlignments(
    seqnames=Rle(factor(c("chr1", "chr2", "chr1", "chr3")),
        c(1, 3, 2, 4)),
    pos=1:10, cigar=paste0(10:1, "M"),
    strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    names=head(letters, 10), score=1:10)

gal2 <- GAlignments(
    seqnames=Rle(factor(c("chr2", "chr4")), c(3, 4)), pos=1:7,
    cigar=c("5M", "3M2N3M2N3M", "5M", "10M", "5M1N4M", "8M2N1M", "5M"),
    strand=Rle(strand(c("-", "+")), c(4, 3)),
    names=tail(letters, 7), score=1:7)

galist <- GAlignmentsList(noGaps=gal1, Gaps=gal2)


## ---------------------------------------------------------------------
## A. BASIC MANIPULATION
## ---------------------------------------------------------------------

length(galist)
names(galist)
seqnames(galist)
strand(galist)
head(cigar(galist))
head(qwidth(galist))
head(start(galist))
head(end(galist))
head(width(galist))
head(njunc(galist))
seqlevels(galist)

## Rename the reference sequences:
seqlevels(galist) <- sub("chr", "seq", seqlevels(galist))
seqlevels(galist)

grglist(galist)  # a GRangesList object
rglist(galist)   # an IRangesList object

## ---------------------------------------------------------------------
## B. SUBSETTING
## ---------------------------------------------------------------------

galist[strand(galist) == "-"]
has_junctions <- sapply(galist,
                        function(x) any(grepl("N", cigar(x), fixed=TRUE)))
galist[has_junctions]

## Different ways to subset:
galist[2]             # a GAlignments object of length 1
galist[[2]]           # a GAlignments object of length 1
grglist(galist[2])  # a GRangesList object of length 1
rglist(galist[2])   # a NormalIRangesList object of length 1

## ---------------------------------------------------------------------
## C. mcols()/elementMetadata()
## ---------------------------------------------------------------------

## Metadata can be defined on the individual GAlignment elements
## and the overall GAlignmentsList object. By default, 'level=between' 
## extracts the GALignmentsList metadata. Using 'level=within' 
## will extract the metadata on the individual GAlignments objects.

mcols(galist) ## no metadata on the GAlignmentsList object
mcols(galist, level="within")


## ---------------------------------------------------------------------
## D. readGAlignmentsList()
## ---------------------------------------------------------------------

library(pasillaBamSubset)

## 'file' as character.
fl <- untreated3_chr4()
galist1 <- readGAlignmentsList(fl)

galist1[1:3]
length(galist1)
table(elementNROWS(galist1))

## When 'file' is a BamFile, 'asMates' must be TRUE. If FALSE,
## the data are treated as single-end and each list element of the
## GAlignmentsList will be of length 1. For single-end data 
## use readGAlignments() instead of readGAlignmentsList().
bf <- BamFile(fl, yieldSize=3, asMates=TRUE)
readGAlignmentsList(bf)

## Use a 'param' to fine tune the results.
param <- ScanBamParam(flag=scanBamFlag(isProperPair=TRUE))
galist2 <- readGAlignmentsList(fl, param=param)
length(galist2)


## ---------------------------------------------------------------------
## E. COERCION 
## ---------------------------------------------------------------------

## The granges() and grlist() coercions support 'ignore.strand' to 
## allow ranges from different strands to be combined. In this example 
## paired-end reads aligned to opposite strands were read into a 
## GAlignmentsList. If the desired operation is to combine these ranges, 
## regardless of junctions or the space between pairs, 'ignore.strand'
## must be TRUE.
granges(galist[1])
granges(galist[1], ignore.strand=TRUE)

## grglist()
galist <- GAlignmentsList(noGaps=gal1, Gaps=gal2)
grglist(galist)
grglist(galist, ignore.strand=TRUE)
}

\keyword{methods}
\keyword{classes}

