\name{summarizeOverlaps-methods}

\alias{summarizeOverlaps-methods}

\alias{summarizeOverlaps}
\alias{summarizeOverlaps,GRanges,GAlignments-method}
\alias{summarizeOverlaps,GRangesList,GAlignments-method}
\alias{summarizeOverlaps,GRanges,GAlignmentsList-method}
\alias{summarizeOverlaps,GRangesList,GAlignmentsList-method}
\alias{summarizeOverlaps,GRanges,GAlignmentPairs-method}
\alias{summarizeOverlaps,GRangesList,GAlignmentPairs-method}
\alias{summarizeOverlaps,GRanges,GRanges-method}
\alias{summarizeOverlaps,GRangesList,GRanges-method}
\alias{summarizeOverlaps,GRanges,GRangesList-method}
\alias{summarizeOverlaps,GRangesList,GRangesList-method}

\alias{Union}
\alias{IntersectionStrict}
\alias{IntersectionNotEmpty}

\alias{summarizeOverlaps,GRanges,BamFile-method}
\alias{summarizeOverlaps,GRangesList,BamFile-method}
\alias{summarizeOverlaps,GRanges,character-method}
\alias{summarizeOverlaps,GRangesList,character-method}
\alias{summarizeOverlaps,GRanges,BamFileList-method}
\alias{summarizeOverlaps,GRangesList,BamFileList-method}
\alias{summarizeOverlaps,BamViews,missing-method}


\title{Perform overlap queries between reads and genomic features} 

\description{
  \code{summarizeOverlaps} extends \code{findOverlaps} by providing 
  options to resolve reads that overlap multiple features. 
}

\usage{
\S4method{summarizeOverlaps}{GRanges,GAlignments}(
  features, reads, mode=Union,
  ignore.strand=FALSE, inter.feature=TRUE, preprocess.reads=NULL, ...)

\S4method{summarizeOverlaps}{GRangesList,GAlignments}(
  features, reads, mode=Union,
  ignore.strand=FALSE, inter.feature=TRUE, preprocess.reads=NULL, ...)

\S4method{summarizeOverlaps}{GRanges,GRanges}(
  features, reads, mode=Union,
  ignore.strand=FALSE, inter.feature=TRUE, preprocess.reads=NULL, ...)

\S4method{summarizeOverlaps}{GRangesList,GRanges}(
  features, reads, mode=Union,
  ignore.strand=FALSE, inter.feature=TRUE, preprocess.reads=NULL, ...)

\S4method{summarizeOverlaps}{GRanges,GAlignmentPairs}(
  features, reads, mode=Union,
  ignore.strand=FALSE, inter.feature=TRUE, preprocess.reads=NULL, ...)

\S4method{summarizeOverlaps}{GRangesList,GAlignmentPairs}(
  features, reads, mode=Union,
  ignore.strand=FALSE, inter.feature=TRUE, preprocess.reads=NULL, ...)

## mode funtions
Union(features, reads, ignore.strand=FALSE,
                       inter.feature=TRUE)
IntersectionStrict(features, reads, ignore.strand=FALSE,
                                    inter.feature=TRUE)
IntersectionNotEmpty(features, reads, ignore.strand=FALSE,
                                      inter.feature=TRUE)

\S4method{summarizeOverlaps}{GRanges,BamFile}(
  features, reads, mode=Union,
  ignore.strand=FALSE, inter.feature=TRUE, singleEnd=TRUE,
  fragments=FALSE, param=ScanBamParam(), preprocess.reads=NULL, ...)

\S4method{summarizeOverlaps}{BamViews,missing}(
  features, reads, mode=Union,
  ignore.strand=FALSE, inter.feature=TRUE, singleEnd=TRUE,
  fragments=FALSE, param=ScanBamParam(), preprocess.reads=NULL, ...)
}

\arguments{
  \item{features}{
    A \link[GenomicRanges]{GRanges} or a \link[GenomicRanges]{GRangesList}
    object of genomic regions of interest. When a \link[GenomicRanges]{GRanges}
    is supplied, each row is considered a feature. When a
    \link[GenomicRanges]{GRangesList} is supplied, each higher list-level is
    considered a feature. This distinction is important when defining
    overlaps.

    When \code{features} is a \link[Rsamtools]{BamViews} the \code{reads}
    argument is missing. Features are extracted from the \code{bamRanges} 
    and the \code{reads} from \code{bamPaths}. Metadata from \code{bamPaths} 
    and \code{bamSamples} are stored in the \code{colData} of the 
    resulting \link[SummarizedExperiment]{RangedSummarizedExperiment} object.
    \code{bamExperiment} metadata are stored in the \code{metadata} slot. 
  } 
  \item{reads}{
    A \link[GenomicRanges]{GRanges}, \link[GenomicRanges]{GRangesList}
    \link{GAlignments}, \link{GAlignmentsList}, 
    \link{GAlignmentPairs}, \link[Rsamtools]{BamViews} or 
    \link[Rsamtools]{BamFileList} object that represents the data to be 
    counted by \code{summarizeOverlaps}.

    \code{reads} is missing when a \link[Rsamtools]{BamViews} object is 
    the only argument supplied to \code{summarizeOverlaps}.
    \code{reads} are the files specified in \code{bamPaths} of the
    \link[Rsamtools]{BamViews} object.
  }
  \item{mode}{
    \code{mode} can be one of the pre-defined count methods such as
    "Union", "IntersectionStrict", or "IntersectionNotEmpty" or
    it a user supplied count function. For a custom count
    function, the input arguments must match those of the pre-defined 
    options and the function must return a vector of counts the same
    length as the annotation ('features' argument). See examples for
    details.

    The pre-defined options are designed after the counting modes 
    available in the HTSeq package by Simon Anders (see references).

    \itemize{
      \item "Union" : (Default) Reads that overlap any portion of exactly one
            feature are counted. Reads that overlap multiple features are
            discarded. This is the most conservative of the 3 modes.
      \item "IntersectionStrict" : A read must fall completely "within" the
            feature to be counted. If a read overlaps multiple features but
            falls "within" only one, the read is counted for that feature.
            If the read is "within" multiple features, the read is discarded.
      \item "IntersectionNotEmpty" : A read must fall in a unique disjoint
            region of a feature to be counted. When a read overlaps multiple
            features, the features are partitioned into disjoint intervals.
            Regions that are shared between the features are discarded leaving
            only the unique disjoint regions. If the read overlaps one of
            these remaining regions, it is assigned to the feature the
            unique disjoint region came from.
      \item user supplied function : A function can be supplied as the
            \code{mode} argument. It must (1) have arguments that correspond
            to \code{features}, \code{reads}, \code{ignore.strand} and
            \code{inter.feature} arguments (as in the defined mode functions)
            and (2) return a vector of counts the same length as
            \code{features}. 
    }
  }
  \item{ignore.strand}{
    A logical indicating if strand should be considered when matching.
  }
  \item{inter.feature}{
    (Default TRUE) A logical indicating if the counting \code{mode} should 
    be aware of overlapping features. When TRUE (default), reads mapping to 
    multiple features are dropped (i.e., not counted). When FALSE, these 
    reads are retained and a count is assigned to each feature they map to.

    There are 6 possible combinations of the \code{mode} and
    \code{inter.feature} arguments. When \code{inter.feature=FALSE} the
    behavior of modes \sQuote{Union} and \sQuote{IntersectionStrict} are 
    essentially \sQuote{countOverlaps} with \sQuote{type=any} and
    \code{type=within}, respectively. \sQuote{IntersectionNotEmpty} does
    not reduce to a simple countOverlaps because common (shared) regions 
    of the annotation are removed before counting.
  }
  \item{preprocess.reads}{
    A function applied to the reads before counting. The first argument
    should be \code{reads} and the return value should be an object 
    compatible with the \code{reads} argument to the counting modes,
    Union, IntersectionStrict and IntersectionNotEmpty.

    The distinction between a user-defined 'mode' and user-defined
    'preprocess.reads' function is that in the first case the user
    defines how to count; in the second case the reads are preprocessed
    before counting with a pre-defined mode. See examples.
  }
  \item{...}{
    Additional arguments passed to functions or methods called
    from within \code{summarizeOverlaps}. For BAM file methods 
    arguments may include \code{singleEnd}, \code{fragments} or 
    \code{param} which apply to reading records from a file 
    (see below). Providing \code{count.mapped.reads=TRUE} include 
    additional passes through the BAM file to collect statistics 
    similar to those from \code{countBam}.

    A \code{BPPARAM} argument can be passed down to the \code{bplapply} 
    called by \code{summarizeOverlaps}. The argument can be MulticoreParam(),
    SnowParam(), BatchJobsParam() or DoparParam(). See the
    BiocParallel package for details in specifying the params.
  }
  \item{singleEnd}{
    (Default TRUE) A logical indicating if reads are single or 
    paired-end. In Bioconductor > 2.12 it is not necessary to sort
    paired-end BAM files by \code{qname}. When counting with 
    \code{summarizeOverlaps}, setting \code{singleEnd=FALSE} will trigger 
    paired-end reading and counting. It is fine to also set 
    \code{asMates=TRUE} in the \code{BamFile} but is not necessary when
    \code{singleEnd=FALSE}.
  }
  \item{fragments}{
    (Default FALSE) A logical; applied to paired-end data only.
    \code{fragments} controls which function is used to read the data which
    subsequently affects which records are included in counting.

    When \code{fragments=FALSE}, data are read with
    \code{\link{readGAlignmentPairs}} and returned in a \code{GAlignmentPairs}
    class. This class only holds \sQuote{mated pairs} from opposite strands;
    same-strand pairs singletons, reads with unmapped pairs and other fragments
    are dropped.

    When \code{fragments=TRUE}, data are read with
    \code{\link{readGAlignmentsList}} and returned in a \code{GAlignmentsList}
    class. This class holds \sQuote{mated pairs} as well as same-strand pairs,
    singletons, reads with unmapped pairs and other fragments. Because more
    records are kept, generally counts will be higher when
    \code{fragments=TRUE}.

    The term \sQuote{mated pairs} refers to records paired with the algorithm
    described on the \code{?\link{readGAlignmentsList}} man page.
  }
  \item{param}{An optional \code{\link[Rsamtools]{ScanBamParam}} instance to
    further influence scanning, counting, or filtering.

    See \code{?\link{BamFile}} for details of how records are returned
    when both \code{yieldSize} is specified in a \code{\link{BamFile}} and
    \code{which} is defined in a \code{\link{ScanBamParam}}.
  }
}

\details{
  \describe{
    \item{}{\code{summarizeOverlaps} offers counting modes to resolve reads that
      overlap multiple features. The \code{mode} argument defines a set of rules
      to resolve the read to a single feature such that each read is counted a
      maximum of once. New to GenomicRanges >= 1.13.9 is the
      \code{inter.feature} argument which allows reads to be counted for each
      feature they overlap.  When \code{inter.feature=TRUE} the counting modes
      are aware of feature overlap; reads that overlap multiple features are
      dropped and not counted.  When \code{inter.feature=FALSE} multiple feature
      overlap is ignored and reads are counted once for each feature they map
      to. This essentially reduces modes \sQuote{Union} and
      \sQuote{IntersectionStrict} to \code{countOverlaps} with
      \code{type="any"}, and \code{type="within"}, respectively.
      \sQuote{IntersectionNotEmpty} is not reduced to a derivative of
      \code{countOverlaps} because the shared regions are removed before
      counting.

      The \code{BamViews}, \code{BamFile} and \code{BamFileList} methods
      summarize overlaps across one or several files. The latter uses
      \code{bplapply}; control parallel evaluation using the
      \code{\link{register}} interface in the \pkg{BiocParallel} package.
    }
    \item{features :}{
      A \sQuote{feature} can be any portion of a genomic region such as a gene, 
      transcript, exon etc. When the \code{features} argument is a 
      \link[GenomicRanges]{GRanges} the rows define the features. The result
      will be the same length as the \link[GenomicRanges]{GRanges}. When 
      \code{features} is a \link[GenomicRanges]{GRangesList} the highest
      list-level defines the features and the result will be the same length
      as the \link[GenomicRanges]{GRangesList}. 

      When \code{inter.feature=TRUE}, each count \code{mode} attempts to 
      assign a read that overlaps multiple features to a single feature. If 
      there are ranges that should be considered together (e.g., exons by 
      transcript or cds regions by gene) the \link[GenomicRanges]{GRangesList}
      would be appropriate. If there is no grouping in the data then a 
      \link[GenomicRanges]{GRanges} would be appropriate. 
    }
    \item{paired-end reads :}{ 
      Paired-end reads are counted as a single hit if one or both parts
      of the pair are overlapped. Paired-end records can be counted in 
      a \link{GAlignmentPairs} container or BAM file. 

      Counting pairs in BAM files:
      \itemize{
        \item{The \code{singleEnd} argument should be FALSE.}
        \item{When \code{reads} are supplied as a BamFile or BamFileList, 
              the \code{asMates} argument to the BamFile should be TRUE.}
        \item{When \code{fragments} is FALSE, a \code{GAlignmentPairs}
              object is used in counting (pairs only).}
        \item{When \code{fragments} is TRUE, a \code{GAlignmentsList}
              object is used in counting (pairs, singletons, unmapped
              mates, etc.)}
      }
    }
  }
}

\value{
  A \link[SummarizedExperiment]{RangedSummarizedExperiment} object. The
  \code{assays} slot holds the counts, \code{rowRanges} holds the annotation
  from \code{features}.

  When \code{reads} is a \code{BamFile} or \code{BamFileList} \code{colData} 
  is an empty DataFrame with a single row named \sQuote{counts}. If 
  \code{count.mapped.reads=TRUE}, \code{colData} holds the output of 
  \code{countBam} in 3 columns named \sQuote{records} (total records), 
  \sQuote{nucleotides} and \sQuote{mapped} (mapped records).

  When \code{features} is a \code{BamViews} \code{colData} includes
  2 columns named \code{bamSamples} and \code{bamIndices}. 

  In all other cases, \code{colData} has columns of \sQuote{object} 
  (class of reads) and \sQuote{records} (length of \code{reads}). 
}

\references{
  HTSeq :
  \url{http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html}

  htseq-count :
  \url{http://www-huber.embl.de/users/anders/HTSeq/doc/count.html}
}

\author{Valerie Obenchain <vobencha@fredhutch.org>}

\seealso{
  \itemize{
    \item The \pkg{DESeq2}, \pkg{DEXSeq} and \pkg{edgeR} packages.

    \item The \link[SummarizedExperiment]{RangedSummarizedExperiment} class
          defined in the \pkg{SummarizedExperiment} package.

    \item The \link{GAlignments} and \link{GAlignmentPairs} classes.

    \item The \link[Rsamtools]{BamFileList} and \link[Rsamtools]{BamViews}
          classes in the \pkg{Rsamtools} package.

    \item The \link{readGAlignments} and \link{readGAlignmentPairs} functions.
  }
}

\examples{
reads <- GAlignments(
    names = c("a","b","c","d","e","f","g"),
    seqnames = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
    pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
    cigar = c("500M", "100M", "300M", "500M", "300M", 
              "50M200N50M", "50M150N50M"),
    strand = strand(rep("+", 7)))

gr <- GRanges(
    seqnames = c(rep("chr1", 7), rep("chr2", 4)), strand = "+", 
    ranges = IRanges(c(1000, 3000, 3600, 4000, 4000, 5000, 5400, 
                       2000, 3000, 7000, 7500), 
                     width = c(500, 500, 300, 500, 900, 500, 500, 
                               900, 500, 600, 300),
                     names=c("A", "B", "C1", "C2", "D1", "D2", "E", "F",
                             "G", "H1", "H2"))) 
groups <- factor(c(1,2,3,3,4,4,5,6,7,8,8))
grl <- splitAsList(gr, groups)
names(grl) <- LETTERS[seq_along(grl)]

## ---------------------------------------------------------------------
## Counting modes. 
## ---------------------------------------------------------------------

## First count with a GRanges as the 'features'. 'Union' is the
## most conservative counting mode followed by 'IntersectionStrict' 
## then 'IntersectionNotEmpty'.
counts1 <- 
    data.frame(union=assays(summarizeOverlaps(gr, reads))$counts, 
               intStrict=assays(summarizeOverlaps(gr, reads, 
                                mode="IntersectionStrict"))$counts,
               intNotEmpty=assays(summarizeOverlaps(gr, reads,
                                  mode="IntersectionNotEmpty"))$counts)

colSums(counts1)

## Split the 'features' into a GRangesList and count again.
counts2 <- 
    data.frame(union=assays(summarizeOverlaps(grl, reads))$counts, 
               intStrict=assays(summarizeOverlaps(grl, reads, 
                                mode="IntersectionStrict"))$counts,
               intNotEmpty=assays(summarizeOverlaps(grl, reads,
                                  mode="IntersectionNotEmpty"))$counts)
colSums(counts2)

## The GRangesList ('grl' object) has 8 features whereas the GRanges 
## ('gr' object) has 11. The affect on counting can be seen by looking
## at feature 'H' with mode 'Union'. In the GRanges this feature is 
## represented by ranges 'H1' and 'H2',
gr[c("H1", "H2")]

## and by list element 'H' in the GRangesList, 
grl["H"]
 
## Read "d" hits both 'H1' and 'H2'. This is considered a multi-hit when
## using a GRanges (each range is a separate feature) so the read was 
## dropped and not counted.
counts1[c("H1", "H2"), ]

## When using a GRangesList, each list element is considered a feature.
## The read hits multiple ranges within list element 'H' but only one 
## list element. This is not considered a multi-hit so the read is counted.
counts2["H", ]

## ---------------------------------------------------------------------
## Counting multi-hit reads.
## ---------------------------------------------------------------------

## The goal of the counting modes is to provide a set of rules that
## resolve reads hitting multiple features so each read is counted
## a maximum of once. However, sometimes it may be desirable to count 
## a read for each feature it overlaps. This can be accomplished by 
## setting 'inter.feature' to FALSE.

## When 'inter.feature=FALSE', modes 'Union' and 'IntersectionStrict'
## essentially reduce to countOverlaps() with type="any" and 
## type="within", respectively.

## When 'inter.feature=TRUE' only features "A", "F" and "G" have counts.
se1 <- summarizeOverlaps(gr, reads, mode="Union", inter.feature=TRUE)
assays(se1)$counts

## When 'inter.feature=FALSE' all 11 features have a count. There are 
## 7 total reads so one or more reads were counted more than once.
se2 <- summarizeOverlaps(gr, reads, mode="Union", inter.feature=FALSE)
assays(se2)$counts

## ---------------------------------------------------------------------
## Counting BAM files.
## ---------------------------------------------------------------------

library(pasillaBamSubset)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
exbygene <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "gene")

## (i) Single-end :

## Large files can be iterated over in chunks by setting a
## 'yieldSize' on the BamFile.
bf_s <- BamFile(untreated1_chr4(), yieldSize=50000)
se_s <- summarizeOverlaps(exbygene, bf_s, singleEnd=TRUE)
table(assays(se_s)$counts > 0)

## When a character (file name) is provided as 'reads' instead 
## of a BamFile object summarizeOverlaps() will create a BamFile
## and set a reasonable default 'yieldSize'.

## (ii) Paired-end :

## A paired-end file may contain singletons, reads with unmapped
## pairs or reads with more than two fragments. When 'fragments=FALSE'
## only reads paired by the algorithm are included in the counting. 
nofrag <- summarizeOverlaps(exbygene, untreated3_chr4(), 
                            singleEnd=FALSE, fragments=FALSE)
table(assays(nofrag)$counts > 0)

## When 'fragments=TRUE' all singletons, reads with unmapped pairs 
## and other fragments will be included in the counting.
bf <- BamFile(untreated3_chr4(), asMates=TRUE)
frag <- summarizeOverlaps(exbygene, bf, singleEnd=FALSE, fragments=TRUE)
table(assays(frag)$counts > 0)

## As expected, using 'fragments=TRUE' results in a larger number 
## of total counts because singletons, unmapped pairs etc. are 
## included in the counting.

## Total reads in the file:
countBam(untreated3_chr4())

## Reads counted with 'fragments=FALSE':
sum(assays(nofrag)$counts)

## Reads counted with 'fragments=TRUE':
sum(assays(frag)$counts)

## ---------------------------------------------------------------------
## Use ouput of summarizeOverlaps() for differential expression analysis
## with DESeq2 or edgeR.
## ---------------------------------------------------------------------

fls <- list.files(system.file("extdata", package="GenomicAlignments"),
                  recursive=TRUE, pattern="*bam$", full=TRUE)
names(fls) <- basename(fls)
bf <- BamFileList(fls, index=character(), yieldSize=1000)
genes <- GRanges(
    seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
    ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 
                       4000, 7500, 5000, 5400), 
                     width=c(rep(500, 3), 600, 900, 500, 300, 900, 
                             300, 500, 500))) 
se <- summarizeOverlaps(genes, bf)

## When the reads are BAM files, the 'colData' contains summary 
## information from a call to countBam().
colData(se)

## Start differential expression analysis with the DESeq2 or edgeR
## package:
library(DESeq2)
deseq <- DESeqDataSet(se, design= ~ 1)
library(edgeR)
edger <- DGEList(assays(se)$counts, group=rownames(colData(se)))

## ---------------------------------------------------------------------
## Filter records by map quality before counting. 
## (user-supplied 'mode' function) 
## ---------------------------------------------------------------------

## The 'mode' argument can take a custom count function whose
## arguments are the same as those in the current counting modes
## (i.e., Union, IntersectionNotEmpty, IntersectionStrict). 
## In this example records are filtered by map quality before counting.

mapq_filter <- function(features, reads, ignore.strand, inter.feature)
{ 
    require(GenomicAlignments) # needed for parallel evaluation
    Union(features, reads[mcols(reads)$mapq >= 20],
          ignore.strand, inter.feature) 
}

genes <- GRanges("seq1", IRanges(seq(1, 1500, by=200), width=100))
param <- ScanBamParam(what="mapq")
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
se <- summarizeOverlaps(genes, fl, mode=mapq_filter, param=param) 
assays(se)$counts

## The count function can be completely custom (i.e., not use the
## pre-defined count functions at all). Requirements are that
## the input arguments match the pre-defined modes and the output
## is a vector of counts the same length as 'features'. 
 
my_count <- function(features, reads,  ignore.strand, inter.feature) { 
    ## perform filtering, or subsetting etc. 
    require(GenomicAlignments) # needed for parallel evaluation
    countOverlaps(features, reads)
}

## ---------------------------------------------------------------------
## Preprocessing reads before counting with a standard count mode.
## (user-supplied 'preprocess.reads' function) 
## ---------------------------------------------------------------------

## The 'preprocess.reads' argument takes a function that is
## applied to the reads before counting with a pre-defined mode.

ResizeReads <- function(reads, width=1, fix="start", ...) {
    reads <- as(reads, "GRanges")
    stopifnot(all(strand(reads) != "*"))
    resize(reads, width=width, fix=fix, ...)
}

## By default ResizeReads() counts reads that overlap on the 5' end:
summarizeOverlaps(grl, reads, mode=Union, preprocess.reads=ResizeReads)

## Count reads that overlap on the 3' end by passing new values
## for 'width' and 'fix':
summarizeOverlaps(grl, reads, mode=Union, preprocess.reads=ResizeReads,
                  width=1, fix="end")

## ---------------------------------------------------------------------
## summarizeOverlaps() with BamViews.
## ---------------------------------------------------------------------

## bamSamples and bamPaths metadata are included in the colData.
## bamExperiment metadata is put into the metadata slot.
fl <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
rngs <- GRanges(c("seq1", "seq2"), IRanges(1, c(1575, 1584)))
samp <- DataFrame(info="test", row.names="ex1")
view <- BamViews(fl, bamSamples=samp, bamRanges=rngs)
se <- summarizeOverlaps(view, mode=Union, ignore.strand=TRUE)
colData(se)
metadata(se)
}

\keyword{methods}
\keyword{utilities}
