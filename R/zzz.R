.onUnload <- function(libpath)
{
    library.dynam.unload("GenomicAlignments", libpath)
}

.test <- function() BiocGenerics:::testPackage("GenomicAlignments")

