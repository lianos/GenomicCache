setClassUnion("MaybeTranscriptDb", c('TranscriptDb', 'NULL'))
setClass("GenomicFeaturesX", contains="VIRTUAL")

setClass("GenomicCache",
         representation(.txdb='MaybeTranscriptDb',
                        .cache="environment"))

setClass("GFGene",
         representation(.id='character',
                        .entrez.id='character',
                        .symbol='character',
                        .chromosome='factor',
                        .strand='factor',
                        .exons="GRangesList",
                        .cds="GRangesList",
                        .utr5="GRangesList",
                        .utr3="GRangesList",
                        .transcript.names="character",
                        .genome='character',
                        .cache='environment'))

setClass("RefSeqGene", contains="GFGene")
setClass("EnsemblGene", contains="GFGene")
setClass("AceviewGene", contains="GFGene")

################################################################################
## Methods : Generic (work on all GenomicFeaturesX-type objects)

setGeneric("duplicate", function(x, ...) standardGeneric("duplicate"))
setGeneric("getEgAnnotationMap", function(x, what, ...) {
  standardGeneric("getEgAnnotationMap")
})
setGeneric("annotationSource", function(x, ...) {
  standardGeneric("annotationSource")
})
setGeneric("genome", function(x, ...) standardGeneric("genome"))

##setGeneric("cacheFetch", function(x, what, expr) standardGeneric("cachFetch"))
setGeneric("clearCache", function(x, ...) standardGeneric("clearCache"))
setGeneric("dispose", function(x, ...) standardGeneric("dispose"))

################################################################################
## Methods: TranscriptDb + GenomicCache + Gene
setGeneric("exons", function(x, vals=NULL) standardGeneric("exons"))
setGeneric("cds", function(x, ...) standardGeneric("cds"))
setGeneric("transcripts", function(x, ...) standardGeneric("transcripts"))

################################################################################
## Methods: Gene
setGeneric("chromosome", function(object, ...) standardGeneric("chromosome"))

################################################################################
## Methods: GenomicCache + TranscriptDb
setGeneric("chromosomes", function(x, ...) {
  standardGeneric("chromosomes")
})

setGeneric("dataSource", function(x, ...) standardGeneric("dataSource"))
setGeneric("getBsGenome", function(x, ...) standardGeneric("getBsGenome"))
setGeneric("getMetadata", function(x, ...) standardGeneric("getMetadata"))
setGeneric("setMetadata", function(x, ...) standardGeneric("setMetadata"))

################################################################################
## Methods: GenomicCache
