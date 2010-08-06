## This GFGene class is really handy when interactively working with your data.
## Imagine you want to see what's happening over a certain gene, or some such.

setMethod("initialize", "GFGene",
function(.Object, ...,
         .id=character(),
         .entrez.id=character(),
         .symbol=character(),
         .chromosome=factor(),
         .strand=strand(),
         .exons=GRangesList(),
         .cds=GRangesList(),
         .utr5=GRangesList(),
         .utr3=GRangesList(),
         .transcript.names=character(),
         .genome=character(),
         .cache=new.env()) {
  callNextMethod(.Object,
                 .id=.id,
                 .entrez.id=.entrez.id,
                 .symbol=.symbol,
                 .chromosome=.chromosome,
                 .strand=.strand,
                 .exons=.exons,
                 .cds=.cds,
                 .utr5=.utr5,
                 .utr3=.utr3,
                 .transcript.names=.transcript.names,
                 .genome=.genome,
                 .cache=.cache,
                 ...)
})

.guessGeneIdType <- function(id, anno.source) {
  if (anno.source == 'refGene') {
    id.type <- switch(substring(id, 1, 3), NM_='tx.id', NR_='tx.id', "symbol")
  } else if (anno.source == 'ensGene') {
    id.type <- switch(substring(id, 1, 4), ENST='tx.id', ENSG='id', 'symbol')
  } else {
    stop("Don't know how to parse gene names of this type: ", anno.source)
  }

  id.type
}

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
    id.type <- .guessGeneIdType(id, annotationSource(.gc))
  }
  
  if (is.null(id.type)) {
    id.type <- .guessGeneIdType(id, .gc)
  }

  ## Necessary monkey business to "do the right thing" in order to get
  ## the correct entrez.id and symbol for this gene
  if (id.type == 'symbol') {
    entrez.id <- getEntrezIdFromSymbol(.gc, id)
  } else {
    if (class.name == 'EnsemblGene') {
      if (id.type == 'id') {
        entrez.id <- getEntrezIdFromGeneId(.gc, id)
      } else {
        entrez.id <- getEntrezIdFromTranscriptId(.gc, id)
      }
    } else {
      entrez.id <- getEntrezIdFromTranscriptId(.gc, id)      
    }
  }

  ## DEBUG: Why does more than one entrez id come back -- and even for the
  ##        wrong transcript id? Try ENST00000270722, it returns results for
  ##        ENST000002707221 and ENST000002707222
  if (!is.null(entrez.id) && !(id %in% names(entrez.id))) {
    entrez.id <- NULL
  }
  if (is.null(entrez.id)) {
    stop("Unknown gene identifier: ", id, " (is it a ", id.type, "?)")
  }
  
  symbol <- getSymbolFromEntrezId(.gc, entrez.id)

  
  .exons <- exonsBy(.gc, 'tx')
  .transcripts <- transcripts(.gc)

  ## Get all transcript IDs associated with this gene
  tx.name <- getTranscriptIdFromEntrezId(.gc, entrez.id)
  xcripts <- subset(.transcripts, values(.transcripts)$tx_name %in% tx.name)
  xm <- values(xcripts)
  tx.ids <- as.character(xm$tx_id)

  ## Errors fly when you try to index a GRangesList with a multiple keys, and
  ## one of which isn't present. This is why I lapply of the tx.ids instead
  ## of cdsBy(.gd, 'tx)[tx.ids]
  take <- function(grl, idxs) {
    x <- GRangesList(lapply(idxs, function(idx) {
      xx <- grl[[idx]]
      if (is.null(xx)) xx <- GRanges()
      xx
    }))
    names(x) <- idxs
    x
  }
  
  g.exons <- .exons[tx.ids]
  g.cds <- take(cdsBy(.gc, 'tx'), tx.ids)
  g.utr5 <- take(fiveUTRsByTranscript(.gc), tx.ids)
  g.utr3 <- take(threeUTRsByTranscript(.gc), tx.ids)
  
  new(class.name,
      .id=id,
      .entrez.id=entrez.id,
      .symbol=symbol,
      .strand=as.vector(strand(xcripts)[1]),
      .chromosome=as.vector(seqnames(xcripts)[1]),
      .exons=g.exons,
      .cds=g.cds,
      .utr5=g.utr5,
      .utr3=g.utr3,
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
  cat(sprintf("%s[%s], %s(%s) %d transcripts\n",
              symbol(object), asource, chromosome(object), strand(object),
              length(xcripts)))
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
  bounds <- unlist(endoapply(transcripts(x), function(xx) range(xx)),
                   use.names=FALSE)
  values(bounds) <- DataFrame(tx_name=x@.transcript.names,
                              tx_id=names(transcripts(x)))
  bounds
})

##' Returns a GRanges object
setGeneric("cdsBounds", function(x, ...) standardGeneric("cdsBounds"))
setMethod("cdsBounds", c(x="GFGene"),
function(x, ...) {
  bounds <- unlist(endoapply(cds(x), function(xx) range(xx)),
                   use.names=FALSE)
  values(bounds) <- DataFrame(tx_name=x@.transcript.names,
                              tx_id=names(transcripts(x)))
  bounds
})

#' Returns RangesList
setMethod("cds", c(x="GFGene"),
function(x, ...) {
  x@.cds
})

setGeneric("utr5", function(x, ...) standardGeneric("utr5"))
setMethod("utr5", c(x="GFGene"),
function(x, ...) {
  x@.utr5
})

setGeneric("utr3", function(x, ...) standardGeneric("utr3"))
setMethod("utr3", c(x="GFGene"),
function(x, ...) {
  x@.utr3
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
setGeneric("isProteinCoding", function(x, ...) {
  standardGeneric("isProteinCoding")
})
setMethod("isProteinCoding", c(x="GFGene"),
function(x, ...) {
  sapply(cds(x), function(exons) length(exons) > 0)
})

setGeneric("id", function(x, ...) standardGeneric("id"))
setMethod("id", c(x="GFGene"),
function(x, ...) {
  x@.id
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
  x@.exons
})

setMethod("exons", c(x="GFGene"),
function(x) {
  x@.exons
})
