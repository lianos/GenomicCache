setClassUnion("MaybeTranscriptDb", c('TranscriptDb', 'NULL'))

##' A virtual class for annotation objects/tracks to extend
##'
##' @exportClass GenomicFeaturesX
##'
##' @slot .genome The genome abbreviation (hg18, etc.) this
##' \code{GenomicFeatureX} is annotted against.
##' @slot .cache An \code{environment} to hold calculated objects for quick
##' retrieval.
setClass("GenomicFeaturesX",
         contains=c("VIRTUAL"),
         representation=representation(
           .genome='character',
           .cache='environment'),
         prototype=prototype(
           .genome=character(),
           .cache=new.env()))

##' The universal/super object that contains all annotations for a given genome.
##' This will minimally have a \code{\link{TranscriptDb}} object.
##'
##' @exportClass GenomicCache
##'
##' @slot .path The absolute path the parent directory of the container
##' @slot .txdb The \code{link{TranscriptDb}} object for the genome.
setClass("GenomicCache",
         contains=c("GenomicFeaturesX"),
         representation=representation(
           .path='character',
           .txdb='MaybeTranscriptDb'),
         prototype=prototype(
           .path=character(),
           .txdb=NULL))

##' The base class for a gene which stores transcript, cds, utr3/5, transcript
##' bounds, etc.
##'
##' @exportClass GFGene
##'
##' @slot .id The gene id (primary key) in the \code{\link{TranscriptDb}}
##' @slot .entrez.id The entrez id for the gene
##' @slot .symbol The hgnc gene symbol
##' @slot .chromosome The chromosome(s) the transcripts of the gene are
##' annotated to
##' @slot .strand The strand(s) the transcripts are annotated to.
##' @slot .exons A \code{\link{GRangesList}} containing as many items as there
##' are annotated transcripts. Each item defines the exon structure of the
##' transcript.
##' @slot .cds A \code{\link{GRangesList}} containing the cds exons for each
##' transcript.
##' @slot .utr5 A \code{\link{GRangesList}} containing the 5' UTRs for each
##' transcript
##' @slot .utr3 A \code{\link{GRangesList}} containgin the 3' UTRs for each
##' transcript.
setClass("GFGene",
         contains=c("GenomicFeaturesX"),
         representation=representation(
           .id='character',
           .entrez.id='character',
           .symbol='character',
           .chromosome='factor',
           .strand='factor',
           .exons="GRangesList",
           .cds="GRangesList",
           .utr5="GRangesList",
           .utr3="GRangesList"),
         prototype=prototype(
           .id=character(),
           .entrez.id=character(),
           .symbol=character(),
           .chromosome=factor(),
           .strand=strand(),
           .cds=GRangesList(),
           .utr5=GRangesList(),
           .utr3=GRangesList()))

##' @exportClass RefSeqGene
setClass("RefSeqGene", contains=c("GFGene"))

##' @exportClass EnsemblGene
setClass("EnsemblGene", contains=c("GFGene"))

##' @exportClass AceviewGene
setClass("AceviewGene", contains=c("GFGene"))

##' @exportClass UcscGene
setClass("UcscGene", contains=c("GFGene"))

##' A subclass of \code{\linkS4class{GRanges}} which stores the complete
##' intervals spanning the length of a chromosome. Each range has a
##' corresponding entry in the \code{values()} of the ranges with
##' \code{exon.anno}, \code{symbol}, \code{entrez.id} columns providing the
##' annotation for the interval.
##'
##' These classes are generated from the \code{\link{annotateChromosome}}
##' function, and its driver function
##' \code{\link{generateAnnotatedChromosomesByGenes}}. Different annotations
##' are created given the \code{\linkS4class{GenomicCache}}, and the desired
##' length of the flanking regions around the annotated transcription start/stop
##' boundaries given.
##'
##' @exportClass AnnotatedChromosome
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
setClass('AnnotatedChromosome', contains=c("GRanges"))

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

## imported from GenomicRanges
## setGeneric("genome", function(x) standardGeneric("genome"))

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

## == Imported from ShortRead ==
## setGeneric("chromosome", function(x, ...) standardGeneric("chromosome"))
## setGeneric("id", function(x, ...) standardGeneric("id"))

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

setGeneric("getGenesOnChromosome",
function(x, chromosome, start=NULL, end=NULL, maxgap=0L, minoverlap=1L,
         overlap.type=c('any', 'start', 'end', 'within', 'equal'),
         use.cache=TRUE, ...) {
  standardGeneric("getGenesOnChromosome")
})


################################################################################
## Methods: GenomicCache
setGeneric("cacheDir", function(x, ...) standardGeneric("cacheDir"))
setGeneric("txdb", function(x, ...) standardGeneric("txdb"))
setGeneric("fiveUTRsBy",
function(x, by=c('tx', 'gene'), use.names=FALSE, ...) {
  standardGeneric("fiveUTRsBy")
})
setGeneric("threeUTRsBy",
function(x, by=c('tx', 'gene'), use.names=FALSE, ...) {
  standardGeneric("threeUTRsBy")
})

################################################################################
## Methods: (G)GappedRanges
setGeneric("gwidth", function(x, ...) standardGeneric("gwidth"))

################################################################################
## Uniqueness

setGeneric("getUniquenessMap",
function(x, chromosome, k, ...) {
  standardGeneric("getUniquenessMap")
})

##' ##' Returns a logical vector that is the same length as ranges indicating
##' whether that range is uniquelly mappable in the given genome
##'
##' @param ranges Either a \code{GRanges} object from which the chromosome(s)
##' are inferred, or an \code{IRanges} object
##' or an IRanges object,hg18, mm9, etc.
##' @param genome hg18, mm9, etc.
##' @param chromosome the name of the chromosome these ranges belong
setGeneric("flagUniqueRanges",
function(ranges, genome, chromosome=NULL, ...) {
  standardGeneric("flagUniqueRanges")
})

