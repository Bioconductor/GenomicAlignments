\name{junctions-methods}

\alias{junctions-methods}

\alias{junctions}
\alias{junctions,GAlignments-method}
\alias{junctions,GAlignmentPairs-method}
\alias{junctions,GAlignmentsList-method}

\alias{summarizeJunctions}

% Old stuff:
\alias{introns}


\title{Extract junctions from genomic alignments}

\description{
  Tools for extracting junctions from a \link{GAlignments},
  \link{GAlignmentPairs}, or \link{GAlignmentsList} object.
}

\usage{
## junctions() and summarizeJunctions():

junctions(x, use.mcols=FALSE, ...)

\S4method{junctions}{GAlignments}(x, use.mcols=FALSE)

\S4method{junctions}{GAlignmentPairs}(x, use.mcols=FALSE)

\S4method{junctions}{GAlignmentsList}(x, use.mcols=FALSE, ignore.strand=FALSE)

summarizeJunctions(x, with.revmap=FALSE, genome=NULL)

## Utilities for importing the junction file generated
## by some aligners:

readTopHatJunctions(file, file.is.bed_to_juncs.output=FALSE)

readSTARJunctions(file)
}

\arguments{
  \item{x}{
    A \link{GAlignments}, \link{GAlignmentPairs}, or \link{GAlignmentsList}
    object.
  }
  \item{use.mcols}{
    \code{TRUE} or \code{FALSE} (the default).
    Whether the metadata columns on \code{x} (accessible with \code{mcols(x)})
    should be propagated to the returned object or not.
  }
  \item{...}{
    Additional arguments, for use in specific methods.
  }
  \item{ignore.strand}{
  }
  \item{with.revmap}{
  }
  \item{genome}{
  }
  \item{file}{
  }
  \item{file.is.bed_to_juncs.output}{
  }
}

\details{
  An N operation in the CIGAR of a genomic alignment is interpreted as a
  junction. \code{junctions(x)} will return the genomic ranges of all
  junctions found in \code{x}.

  More precisely, if \code{x} is a \link{GAlignments} object,
  \code{junctions(x)} is equivalent to:
\preformatted{
  psetdiff(granges(x), grglist(x, order.as.in.query=TRUE))
}
  On a \code{x} is a \link{GAlignmentPairs} object, it's equivalent to (but
  faster than):
\preformatted{
  junctions1 <- junctions(first(x))
  junctions2 <- junctions(last(x, invert.strand=TRUE))
  mendoapply(c, junctions1, junctions2)
}
}

\value{
  \code{junctions(x)} returns the genomic ranges of the junctions in a
  \link[GenomicRanges]{GRangesList} object \emph{parallel} to \code{x},
  (i.e. with 1 list element per element in \code{x}).
  If \code{x} has names on it, they're propagated to the returned object.
  If \code{use.mcols} is TRUE and \code{x} has metadata columns on it
  (accessible with \code{mcols(x)}), they're propagated to the returned object.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item \link{GAlignments}, \link{GAlignmentPairs}, and
          \link{GAlignmentsList} objects.

    \item The \link[GenomicRanges]{GRanges} and
          \link[GenomicRanges]{GRangesList} classes defined and documented
          in the \pkg{GenomicRanges} package.
  }
}

\examples{
library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[1]

## ---------------------------------------------------------------------
## A. junctions()
## ---------------------------------------------------------------------

gal <- readGAlignments(bamfile)
table(njunc(gal))  # some alignments have 3 junctions!
juncs <- junctions(gal)
juncs

stopifnot(identical(unname(elementLengths(juncs)), njunc(gal)))

galp <- readGAlignmentPairs(bamfile)
juncs <- junctions(galp)
juncs

stopifnot(identical(unname(elementLengths(juncs)), njunc(galp)))

## ---------------------------------------------------------------------
## B. summarizeJunctions()
## ---------------------------------------------------------------------

## By default, only the "score", "plus_score", and "minus_score"
## metadata columns are returned:
junc_summary <- summarizeJunctions(gal)
junc_summary

## The "score" metadata column reports the total number of alignments
## that "support" each junction, i.e., that have the junction encoded
## in their CIGAR:
median(mcols(junc_summary)$score)

## The "plus_score" and "minus_score" metadata columns report the number
## of "supporting" alignments on each strand:
stopifnot(identical(mcols(junc_summary)$score,
                    mcols(junc_summary)$plus_score +
                    mcols(junc_summary)$minus_score))

## If 'with.revmap' is TRUE, the "revmap" metadata column is added to
## the result. It contains the mapping from elements in the ouput (the
## junctions) to elements in the input ('x'). We use this to compute a
## 'score2' for each junction. We obtain this score by summing the
## mapping qualities of the alignments supporting the junction:
gal <- readGAlignments(bamfile, param=ScanBamParam(what="mapq"))
junc_summary <- summarizeJunctions(gal, with.revmap=TRUE)
junc_score2 <- sum(extractList(mcols(gal)$mapq,
                               mcols(junc_summary)$revmap))
mcols(junc_summary)$score2 <- junc_score2

## If a genome is specified thru the 'genome' argument (in which case
## the corresponding BSgenome data package needs to be installed), then
## summarizeJunctions() returns the intron strand and motif for each
## junction. Since the reads in RNAseqData.HNRNPC.bam.chr14 were aligned
## to the hg19 genome, the following requires that you have
## BSgenome.Hsapiens.UCSC.hg19 installed:
junc_summary <- summarizeJunctions(gal, with.revmap=TRUE, genome="hg19")
mcols(junc_summary)$score2 <- junc_score2  # putting 'score2' back

## The "intron_motif" metadata column is a factor whose levels are the
## 5 natural intron motifs described at
## http://www.ncbi.nlm.nih.gov/pmc/articles/PMC84117/
table(mcols(junc_summary)$intron_motif)

## ---------------------------------------------------------------------
## C. STRANDED RNA-seq PROTOCOL
## ---------------------------------------------------------------------

## Here is a simple test for checking whether the RNA-seq protocol was
## stranded or not:
strandedTest <- function(plus_score, minus_score)
    (sum(plus_score ^ 2) + sum(minus_score ^ 2)) /
        sum((plus_score + minus_score) ^ 2)

## The result of this test is guaranteed to be >= 0.5 and <= 1.
## If, for each junction, the strand of the supporting alignments looks
## random (i.e. "plus_score" and "minus_score" are close), then
## strandedTest() will return a result close to 0.5. If it doesn't look
## random (i.e. for each junction, one of "plus_score" and "minus_score"
## is much bigger than the other), then strandedTest() will return a
## result close to 1.
strandedTest(mcols(junc_summary)$plus_score,
             mcols(junc_summary)$minus_score)

## ---------------------------------------------------------------------
## UTILITIES FOR IMPORTING THE JUNCTION FILE GENERATED BY SOME ALIGNERS
## ---------------------------------------------------------------------

## The TopHat aligner generates a junctions.bed file where it reports
## all the junctions satisfying some "quality" criteria (see the TopHat
## manual at http://tophat.cbcb.umd.edu/manual.shtml for more
## information). We can load this file with readTopHatJunctions():
runname <- names(RNAseqData.HNRNPC.bam.chr14_BAMFILES)[1]
junctions_file <- system.file("extdata", "tophat2_out", runname,
                              "junctions.bed",
                              package="RNAseqData.HNRNPC.bam.chr14")
th_junctions <- readTopHatJunctions(junctions_file)

## Comparing the "TopHat junctions" with the result of
## summarizeJunctions():
th_junctions14 <- th_junctions
seqlevels(th_junctions14, force=TRUE) <- "chr14"
mcols(th_junctions14)$intron_strand <- strand(th_junctions14)
strand(th_junctions14) <- "*"

## All the "TopHat junctions" are in 'junc_summary':
stopifnot(all(th_junctions14 \%in\% junc_summary))

## But not all the junctions in 'junc_summary' are reported by TopHat
## (that's because TopHat reports only junctions that satisfy some
## "quality" criteria):
is_in_th_junctions14 <- junc_summary \%in\% th_junctions14
table(is_in_th_junctions14)  # 32 junctions are not in TopHat's
                             # junctions.bed file
junc_summary2 <- junc_summary[is_in_th_junctions14]

## 'junc_summary2' and 'th_junctions14' contain the same junctions in
## the same order:
stopifnot(all(junc_summary2 == th_junctions14))

## Let's merge their metadata columns. We use our own version of
## merge() for this, which is stricter (it checks that the common
## columns are the same in the 2 data frames to merge) and also
## simpler:
merge2 <- function(df1, df2)
{
    common_colnames <- intersect(colnames(df1), colnames(df2))
    lapply(common_colnames,
           function(colname)
             stopifnot(all(df1[ , colname] == df2[ , colname])))
    extra_mcolnames <- setdiff(colnames(df2), colnames(df1))
    cbind(df1, df2[ , extra_mcolnames, drop=FALSE])
}

mcols(th_junctions14) <- merge2(mcols(th_junctions14),
                                mcols(junc_summary2))

## Here is a peculiar junction reported by TopHat:
idx0 <- which(mcols(th_junctions14)$score2 == 0L)
th_junctions14[idx0]
gal[mcols(th_junctions14)$revmap[[idx0]]]
## The junction is supported by 5 alignments (score is 5), all of
## which have a mapping quality of 0!
}

\keyword{methods}
\keyword{manip}
