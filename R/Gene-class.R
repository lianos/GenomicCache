setMethod("initialize", "GFGene",
function(.Object, ...,
         .id=character(),
         .entrez.id=character(),
         .symbol=character(),
         .chromosome=factor(),
         .strand=strand(),
         .transcripts=GRangesList(),
         .transcript.names=character(),
         .txBounds=GRanges(),
         .cdsBounds=GRanges(),
         .genome=character(),
         .cache=new.env()) {
  callNextMethod(.Object,
                 .id=.id,
                 .entrez.id=.entrez.id,
                 .symbol=.symbol,
                 .chromosome=.chromosome,
                 .strand=.strand,
                 .transcripts=.transcripts,
                 .transcript.names=.transcript.names,
                 .txBounds=.txBounds,
                 .cdsBounds=.cdsBounds,
                 .genome=.genome,
                 .cache=.cache,
                 ...)
})

##' Create a Gene object.
##' This requires the use of a \code{\link{GenomicCache}} object \code{.gc},
##' which is, for convenience, automatically "pulled out" from the argument list
##' (\code{...}) if not explicitly passed.
GFGene <- function(..., .gc=NULL) {
  args <- list(...)
  arg.names <- names(args)

  if (is.null(.gc)) {
    .gc <- takeFromListByType(args, 'GenomicCache')
    if (is.null(.gc)) {
      stop("GenomicCache (.gc) object expected.")
    }
  }
  
  a.source <- annotationSource(.gc)
  class.name <- switch(a.source, refGene="RefSeqGene",
                       ensGene="EnsemblGene", aceGene="AceviewGene",
                       stop("Unknown source: ", a.source))
  
  id <- takeFromListByType(args, 'character', index=TRUE)
  if (is.null(id)) {
    stop("No gene name/id passed")
  }
  
  id.type <- names(args)[id]
  id <- args[[id]]

  if (is.null(id.type)) {
    id.type <- 'symbol'
  }
  
  if (id.type == 'symbol') {
    symbol <- id
    entrez.id <- getEntrezIdFromSymbol(.gc, symbol)
  } else if (id.type == 'id') {
    ## Will throw an error if its not an ensembl source
    entrez.id <- getEntrezIdFromGeneId(.gc, id)
    symbol <- getSymbolFromEntrezId(.gc, entrez.id)
  } else if (id.type == 'tx.id') {
    entrez.id <- getTranscriptIdFromEntrezId(.gc, id)
    symbol <- getSymbolFromEntrezId(.gc, entrez.id)
  }
  
  if (a.source == 'refGene') {
    id <- character()
  } else if (a.source == 'ensGene') {
    id <- getGeneIdFromEntrezId(.gc, entrez.id)
  } else {
    stop("Unknown source: ", a.source)
  }

  .exons <- exonsBy(.gc, 'tx')
  .transcripts <- transcripts(.gc)
    
  tx.name <- getTranscriptIdFromEntrezId(.gc, entrez.id)
  xcripts <- subset(.transcripts, values(.transcripts)$tx_name %in% tx.name)
  xm <- values(xcripts)
  
  new(class.name,
      .id=id,
      .entrez.id=entrez.id,
      .symbol=symbol,
      .strand=as.vector(strand(xcripts)[1]),
      .chromosome=as.vector(seqnames(xcripts)[1]),
      .transcripts=.exons[as.character(xm$tx_id)],
      .transcript.names=xm$tx_name,
      .genome=genome(.gc)
      )
}

setMethod("show", c(object="GFGene"),
function(object) {
  asource <- switch(substring(annotationSource(object), 1, 3),
                    ens="Ensembl", ref="RefSeq", ace="AceView",
                    "unknown source")
  xcripts <- transcripts(object)
  cat(sprintf("%s (%s), %d transcripts [%s]\n",
              symbol(object), strand(object), length(xcripts), asource))
  for (idx in 1:length(xcripts)) {
    bounds <- range(ranges(xcripts[[idx]]))
    cat(" ", object@.transcript.names[idx], ": ")
    cat(chromosome(object), ":", sep="")
    cat(format(start(bounds), big.mark=","), "-", sep="")
    cat(format(end(bounds), big.mark=","), "\n")
  }
  
})

setGeneric("txBounds", function(x, ...) standardGeneric("txBounds"))
setMethod("txBounds", c(x="GFGene"),
function(x, ...) {
  x@.txBounds
})

setGeneric("cdsBounds", function(x, ...) standardGeneric("cdsBounds"))
setMethod("cdsBounds", c(x="GFGene"),
function(x, ...) {
  x@.cdsBounds
})

#' Returns RangesList
setMethod("cds", c(x="GFGene"),
function(x, ...) {
  stop("Not implemented yet")
  cd.s <- cdsStart(x, ...)
  cd.e <- cdsEnd(x, ...)
  ranges <- RangesList(lapply(seq_along(cd.s), function(idx) {
    IRanges(start=cd.s[idx], end=cd.e[idx])
  }))
  names(ranges) <- names(cd.s)
  ranges
})

setMethod("ranges", c(x="GFGene"),
function(x, type=c('tx', 'transcript', 'cds', 'coding'),
         summary=TRUE, ...) {
  stop("Not implemented yet")
  ## Return the bounds of the gene's transcripts, or coding regions
})

setMethod("getBsGenome", c(x="GFGene"),
function(x, ...) {
  getBsGenomeFromVersion(genome(x))
})

############################################################## Simple Accessors
setGeneric("id", function(x, ...) standardGeneric("id"))
setMethod("id", c(x="GFGene"),
function(x, ...) {
  xt@.id
})


setMethod("genome", c(x="GFGene"),
function(x, ...) {
  x@.genome
})

setGeneric("symbol", function(x, ...) standardGeneric("symbol"))
setMethod("symbol", c(x="GFGene"),
function(x, ...) {
  x@.symbol
})

setGeneric("chromosome", function(x, ...) standardGeneric("chromosome"))
setMethod("chromosome", c(x="GFGene"), 
function(x, as.DNAString=FALSE, unmasked=TRUE, ...) {
  chr <- as.character(x@.chromosome)[1]
  if (as.DNAString) {
    genome <- getBsGenome(x)
    chr <- genome[[chr]]
    if (unmasked) {
      chr <- unmasked(chr)
    }
  }
  chr
})

setMethod("strand", c(x="GFGene"),
function(x) {
  x@.strand
})

setMethod('start', c(x="GFGene"),
function(x, ...) {
  start(ranges(x, ...))
})

setMethod('end', c(x="GFGene"),
function(x, ...) {
  end(ranges(x, ...))
})

setMethod("annotationSource", c(x="GFGene"),
function(x) {
  switch(class(x)[1],
    RefSeqGene='refGene',
    EnsemblGene='ensGene',
    AceviewGene='aceGene')
})

setMethod("transcripts", c(x="GFGene"),
function(x, ...) {
  x@.transcripts
})

###############################################################################
# Summary
###############################################################################
setMethod("summary", c(object="GFGene"), function(object, ...) {
  stop("Not implemented")
})
