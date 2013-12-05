### =========================================================================
### GappedReads objects
### -------------------------------------------------------------------------
###

setClass("GappedReads",
    contains="GAlignments",
    representation(
        qseq="DNAStringSet"
        ## TODO: Maybe add the read quality? mismatch information?
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setGeneric("qseq", function(x) standardGeneric("qseq"))

setMethod("qseq", "GappedReads", function(x) x@qseq)

### Overriding "qwidth" method for GAlignments objects with a faster
### method.
setMethod("qwidth", "GappedReads", function(x) width(qseq(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GappedReads.qseq <- function(x)
{
    x_qseq <- qseq(x)
    if (class(x_qseq) != "DNAStringSet" || !is.null(names(x_qseq)))
        return("'qseq(x)' must be an unnamed DNAStringSet instance")
    if (length(x_qseq) != length(cigar(x)))
        return("'qseq(x)' and 'cigar(x)' must have the same length")
    if (!identical(width(x_qseq), cigarWidthAlongQuerySpace(cigar(x))))
        return(paste("'width(qseq(x))' and",
                     "'cigarWidthAlongQuerySpace(cigar(x))'",
                     "must be identical"))
    NULL
}

.valid.GappedReads <- function(x)
{
    .valid.GappedReads.qseq(x)
}

setValidity2("GappedReads", .valid.GappedReads,
             where=asNamespace("GenomicAlignments"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GappedReads <- function(seqnames=Rle(factor()), pos=integer(0),
                        cigar=character(0), strand=NULL,
                        qseq=DNAStringSet(),
                        names=NULL, seqlengths=NULL)
{
    galn <- GAlignments(seqnames=seqnames, pos=pos,
                        cigar=cigar, strand=strand,
                        names=names, seqlengths=seqlengths)
    new("GappedReads", galn, qseq=qseq)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod(IRanges:::extractROWS, "GappedReads",
    function(x, i)
    {
        if (missing(i) || !is(i, "Ranges"))
            i <- IRanges:::normalizeSingleBracketSubscript(i, x)
        x@qseq <- IRanges:::extractROWS(x@qseq, i)
        callNextMethod()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining.
###

setMethod("c", "GappedReads",
    function (x, ..., recursive = FALSE)
    {
        stop("coming soon")
    }
)

