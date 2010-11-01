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
                 .genome=.genome,
                 .cache=.cache,
                 ...)
})

.guessGeneIdType <- function(id, anno.source) {
  if (!is.na(suppressWarnings(as.integer(id)))) {
    return("entrez")
  }
  
  if (anno.source == 'refGene') {
    ## refseq
    id.type <- switch(substring(id, 1, 3), NM_='tx.id', NR_='tx.id', "symbol")
  } else if (anno.source == 'ensGene') {
    ## ensemble
    id.type <- switch(substring(id, 1, 4), ENST='tx.id', ENSG='id', 'symbol')
  } else if (anno.source == 'knownGene') {
    ## UCSC genes
    id.type <- switch(substring(id, 1, 2), uc='tx.id', 'symbol')
  } else if (anno.source == 'acembly') {
    ## Aceview
    if (length(grep('.', id, fixed=TRUE)) == 0) {
      id.type <- 'symbol'
    } else {
      id.type <- 'tx.id'
    }
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
  
  anno.source <- annotationSource(.gc)
  class.name <- switch(anno.source, refGene="RefSeqGene",
                       ensGene="EnsemblGene", acembly="AceviewGene",
                       knownGene="UcscGene",
                       stop("Unknown source: ", anno.source))
  genome <- genome(.gc)
  
  id <- takeFromListByType(args, 'character', index=TRUE)
  if (is.null(id)) {
    stop("No gene name/id passed")
  }
  
  id.type <- names(args)[id]
  id <- args[[id]]
  if (is.null(id.type)) {
    id.type <- .guessGeneIdType(id, anno.source)
  }
  
  ## Necessary monkey business to "do the right thing" in order to get
  ## the correct entrez.id and symbol for this gene
  symbol <- NULL
  gene.id <- NULL
  if (id.type == 'symbol') {
    symbol <- id
    entrez.id <- getEntrezIdFromSymbol(.gc, id, anno.source)
  } else if (id.type == 'entrez') {
    entrez.id <- id
  } else {
    if (class.name == 'EnsemblGene') {
      if (id.type == 'id') {
        entrez.id <- getEntrezIdFromGeneId(.gc, id, anno.source)
      } else {
        entrez.id <- getEntrezIdFromTranscriptId(.gc, id, anno.source)
      }
    } else {
      entrez.id <- getEntrezIdFromTranscriptId(.gc, id, anno.source)
    }
  }

  ## DEBUG: Why does more than one entrez id come back -- and even for the
  ##        wrong transcript id? Try ENST00000270722, it returns results for
  ##        ENST000002707221 and ENST000002707222
  ## if (!is.null(entrez.id) && !(id %in% names(entrez.id))) {
  ##   entrez.id <- NULL
  ## }
  if (is.null(entrez.id)) {
    stop("Unknown gene identifier: ", id, " (is it a ", id.type, "?)")
  }
  
  if (is.null(symbol)) {
    symbol <- getSymbolFromEntrezId(.gc, entrez.id, anno.source)
  }
  
  if (is.null(gene.id)) {
    gene.id <- switch(anno.source, refGene=symbol,
                      ensGene=getGeneIdFromEntrezId(.gc, entrez.id),
                      knownGene=symbol, acembly=symbol,
                      stop("Unknown anno.source"))
  }
  
  .exons <- exonsBy(.gc, 'tx')
  .transcripts <- transcripts(.gc)

  ## Get all transcript IDs associated with this gene
  if (class.name != 'AceviewGene') {
    tx.name <- getTranscriptIdFromEntrezId(.gc, entrez.id, anno.source)
    xcripts <- subset(.transcripts, values(.transcripts)$tx_name %in% tx.name)
  } else {
    tx.names <- values(.transcripts)$tx_name
    take <- grep(sprintf('%s\\.', symbol), tx.names, ignore.case=TRUE)
    if (length(take) == 0) {
      stop("Can not find transcripts in Aceview")
    }
    xcripts <- .transcripts[take]      
  }
  
  xm <- values(xcripts)
  tx.ids <- as.character(xm$tx_id)
  meta <- DataFrame(tx_name=xm$tx_name)
  
  ## Errors fly when you try to index a GRangesList with a multiple keys, and
  ## one of which isn't present. This is why I lapply of the tx.ids instead
  ## of cdsBy(.gd, 'tx)[tx.ids]. This results in an empty GRanges object in the
  ## the slot the expected cds, utr, etc. should have been if the transcript
  ## had one.
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
  values(g.exons) <- meta
  
  g.cds <- take(cdsBy(.gc, 'tx'), tx.ids)
  values(g.cds) <- meta
  
  g.utr5 <- take(fiveUTRsByTranscript(.gc), tx.ids)
  values(g.utr5) <- meta
  
  g.utr3 <- take(threeUTRsByTranscript(.gc), tx.ids)
  values(g.utr3) <- meta

  new(class.name,
      .entrez.id=entrez.id,
      .id=gene.id,
      .symbol=symbol,
      .strand=as.factor(strand(xcripts)[1]),
      .chromosome=as.factor(unique(seqnames(xcripts))),
      .exons=g.exons,
      .cds=g.cds,
      .utr5=g.utr5,
      .utr3=g.utr3,
      .genome=genome
      )
}

setMethod("show", c(object="GFGene"),
function(object) {
  asource <- switch(substring(annotationSource(object), 1, 3),
                    ens="Ensembl", ref="RefSeq", ace="AceView",
                    kno="UCSC", "unknown source")
  xcripts <- transcripts(object)
  meta <- values(xcripts)
  
  cat(sprintf("%s|%s [%s]: %d transcripts\n",
              symbol(object), entrezId(object), asource,
              length(xcripts)))
  for (idx in 1:length(xcripts)) {
    xcript <- xcripts[[idx]]
    bounds <- range(ranges(xcript))
    cat(" ", meta$tx_name[idx], ": ")
    cat(as.character(seqnames(xcript)[1]))
    cat(" (", as.character(strand(xcript)[1]), ") : ", sep="")
    cat(format(start(bounds), big.mark=","), "-", sep="")
    cat(format(end(bounds), big.mark=","), "\n")
  }
  
})

setMethod("length", c(x="GFGene"),
function(x) {
  length(x@.exons)
})

setMethod("txBounds", c(x="GFGene"),
function(x, ...) {
  bounds <- unlist(endoapply(transcripts(x, ...), function(xx) range(xx)),
                   use.names=FALSE)
  values(bounds) <- DataFrame(tx_name=txNames(x, ...),
                              tx_id=names(transcripts(x, ...)))
  bounds
})

setMethod("txNames", c(x="GFGene"),
function(x, ...) {
  values(transcripts(x, ...))$tx_name
})

##' Returns a GRanges object
setMethod("cdsBounds", c(x="GFGene"),
function(x, ...) {
  bounds <- unlist(endoapply(cds(x, ...), function(xx) range(xx)),
                   use.names=FALSE)
  values(bounds) <- DataFrame(tx_name=txNames(x, ...),
                              tx_id=names(transcripts(x, ...)))
  bounds
})

#' Returns RangesList
setMethod("cds", c(x="GFGene"),
function(x, which.chr=NULL, ...) {
  filterByChr(x@.cds, which.chr=which.chr)
})

setMethod("utr5", c(x="GFGene"),
function(x, which.chr=NULL, ...) {
  filterByChr(x@.utr5, which.chr=which.chr)
})

setMethod("utr3", c(x="GFGene"),
function(x, which.chr=NULL, ...) {
  filterByChr(x@.utr3, which.chr=which.chr)
})

setMethod("range", c(x="GFGene"),
function(x, by=c('gene', 'tx', 'cds'), ..., na.rm=TRUE) {
  by <- match.arg(by)
  switch(by,
         gene=range(unlist(transcripts(x, ...))),
         tx=endoapply(transcripts(x, ...), range),
         cds=endoapply(cds(x, ...), range))
})

setMethod("ranges", c(x="GFGene"),
function(x, type=c('tx', 'transcript', 'cds', 'coding'),
         summary=TRUE, ...) {
  stop("Not implemented yet")
  ## Return the bounds of the gene's transcripts, or coding regions
})

setMethod("getBsGenome", c(x="GFGene"),
function(x, ...) {
  getBsGenome(genome(x))
})

############################################################## Simple Accessors
setMethod("entrezId", c(x="GFGene"),
function(x, ...) {
  x@.entrez.id
})

setMethod("id", c(x="GFGene"),
function(x, ...) {
  x@.id
})

setMethod("isProteinCoding", c(x="GFGene"),
function(x, which.chr=NULL, ...) {
  sapply(cds(x, which.chr=which.chr, ...), function(exons) length(exons) > 0)
})

setMethod("genome", c(x="GFGene"),
function(x, ...) {
  x@.genome
})

setMethod("symbol", c(x="GFGene"),
function(x, ...) {
  x@.symbol
})

setMethod("chromosome", c(x="GFGene"), 
function(x, as.DNAString=FALSE, unmasked=TRUE, ...) {
  chr <- as.character(x@.chromosome)
  if (as.DNAString) {
    stop("Manually retreive the chromosome from a BSgenome object")
    ## genome <- getBsGenome(x)
    ## chr <- genome[[chr]]
    ## if (unmasked) {
    ##   chr <- unmasked(chr)
    ## }
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

setMethod("annotationSource", c(object="GFGene"),
function(object) {
  switch(class(object)[1],
    RefSeqGene='refGene',
    EnsemblGene='ensGene',
    AceviewGene='acembly',
    UcscGene='knownGene')
})

setMethod("transcripts", c(x="GFGene"),
function(x, which.chr=NULL, ...) {
  filterByChr(x@.exons, which.chr=which.chr)
})

setMethod("exons", c(x="GFGene"),
function(x, ...) {
  transcripts(x, ...)
})
