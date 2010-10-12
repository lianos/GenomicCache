setClassUnion("MaybeTranscriptDb", c('TranscriptDb', 'NULL'))
setClass("GenomicFeaturesX", 
         representation(.genome='character',
                        .cache='environment'),
         contains="VIRTUAL")

setClass("GenomicCache",
         representation(.path='character',
                        .txdb='MaybeTranscriptDb'),
         contains="GenomicFeaturesX")

setClass("GFGene",
         representation(.id='character',
                        .entrez.id='character',
                        .symbol='character',
                        .chromosome='factor',
                        .strand='factor',
                        .exons="GRangesList",
                        .cds="GRangesList",
                        .utr5="GRangesList",
                        .utr3="GRangesList"),
         contains="GenomicFeaturesX")

setClass("RefSeqGene", contains="GFGene")
setClass("EnsemblGene", contains="GFGene")
setClass("AceviewGene", contains="GFGene")
setClass("UcscGene", contains="GFGene")

################################################################################
## Methods : Generic (work on all GenomicFeaturesX-type objects)

setGeneric("duplicate", function(x, ...) standardGeneric("duplicate"))

setGeneric("annotationSource", function(object) {
  standardGeneric("annotationSource")
})

setGeneric("annotationPackage", function(x, ...) {
  standardGeneric("annotationPackage")
})

setMethod("annotationPackage", c(x="GenomicFeaturesX"),
function(x, package=NULL) {
  ## You can pass the annotation package to NULL case you have a corner case
  ## that I haven't thought of.
  annotationPackages(genome(x), package=package)
})

setMethod("annotationPackage", c(x="character"),
function(x, package=NULL) {
  if (is.null(package)) {
    package <- switch(x,
      hg18='org.Hs.eg.db', hg19='org.Hs.eg.db',     ## Human
      mm8='org.Mm.eg.db', mm9='org.Mm.eg.db',       ## Mouse
      rn4='org.Rn.eg.db', rn3='org.Rn.eg.db',       ## Rat
      dm3='org.Dm.eg.db',                           ## Fly
      ce4='org.Ce.eg.db', ce6='org.Ce.eg.db',       ## Worm
      stop("Unknown genome -- provide the annotation package in `package`")
    )
  }
  package
})

setGeneric("genome", function(x, ...) standardGeneric("genome"))

setGeneric("cacheFetch", function(x, what, expr, force.eval=FALSE) {
  standardGeneric("cacheFetch")
})
setGeneric("clearCache", function(x, ...) standardGeneric("clearCache"))
setGeneric("dispose", function(x, ...) standardGeneric("dispose"))
setGeneric("getBsGenome", function(x, ...) standardGeneric("getBsGenome"))

################################################################################
## Methods: Gene
setGeneric("entrezId", function(x, ...) standardGeneric("entrezId"))
setGeneric("cdsBounds", function(x, ...) standardGeneric("cdsBounds"))
setGeneric("chromosome", function(x, ...) standardGeneric("chromosome"))
setGeneric("id", function(x, ...) standardGeneric("id"))
setGeneric("symbol", function(x, ...) standardGeneric("symbol"))
setGeneric("txBounds", function(x, ...) standardGeneric("txBounds"))
setGeneric("txNames", function(x, ...) {
  standardGeneric("txNames")
})
setGeneric("utr5", function(x, ...) standardGeneric("utr5"))
setGeneric("utr3", function(x, ...) standardGeneric("utr3"))
setGeneric("isProteinCoding", function(x, ...) {
  standardGeneric("isProteinCoding")
})

################################################################################
## Methods: GenomicCache + TranscriptDb
setGeneric("chromosomes", function(x, ...) {
  standardGeneric("chromosomes")
})


################################################################################
## Methods: GenomicCache
setGeneric("cacheDir", function(x, ...) standardGeneric("cacheDir"))
setGeneric("txdb", function(x, ...) standardGeneric("txdb"))

## setGeneric("txdbConn", function(txdb) standardGeneric('txdbConn'))
