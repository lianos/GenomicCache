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

  ## Errors fly when you try to index a GRangesList with a key that it doesn't
  ## have. This happens when slicing cds()[...] with an RNA, for example.
  none <- function(e) GRangesList()
  
  g.exons <- .exons[as.character(xm$tx_id)]
  g.cds <- tryCatch(cdsBy(.gc, 'tx')[tx.ids], error=none)
  g.utr5 <- tryCatch(fiveUTRsByTranscript(.gc)[tx.ids], error=none)
  g.utr3 <- tryCatch(threeUTRsByTranscript(.gc)[tx.ids], error=none)
  
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


matchGFGeneSummaryType <- function(what) {
  opts <- c('constitutiveExons', 'cover', 'first')
   if (!is.numeric(what)) {
    what <- match.arg(what, opts)
  } else {
    what <- as.integer(what)
  }
  if (what == 'first') {
    what <- 1L
  }
  what
}

##' Returns an "annotated" GRanges object by "compressing" the exons
##' across isoforms
setGeneric("summarized", function(x, ...) standardGeneric("summarized"))
setMethod("summarized", c(x="GFGene"),
function(x, exon.type='constitutiveExons', cds.cover=c('min', 'max'),
         ...) {
  cds.cover <- match.arg(cds.cover)
  exon.type <- matchGFGeneSummaryType(exon.type)
  
  if (is.integer(exon.type)) {
    if (exon.type > length(transcripts(x)) || exon.type < 1L) {
      stop("Transcript chosen for `exon.type` is out of bounds")
    }
  }

  reducef <- switch(exon.type, cover='union', 'intersect')
  
  if (all(isProteinCoding(x))) {
    if (is.integer(exon.type)) {
      .cds <- cds(x)[[exon.type]]
      .utr3 <- utr3(x)[[exon.type]]
      .utr5 <- utr5(x)[[exon.type]]
    } else {
      .cds <- Reduce(reducef, lapply(cds(x), ranges))
      .utr3 <- Reduce(reducef, lapply(utr3(x), ranges))
      .utr5 <- Reduce(reducef, lapply(utr5(x), ranges))

      if (cds.cover == 'min') {
        .cds <- setdiff(.cds, union(.utr3, .utr5))
      } else {
        .utr3 <- setdiff(.utr3, .cds)
        .utr5 <- setdiff(.utr5, .cds)
      }
    }
    
    ranges <- c(.utr5, .cds, .utr3)
    exon.anno <- c(rep('utr5', length(.utr5)),
                   rep('cds', length(.cds)),
                   rep('utr3', length(.utr3)))
  } else {
    ranges <- Reduce(reducef, lapply(transcripts(x), ranges))
    exon.anno <- c(rep('utr', length(ranges)))
  }
  
  gr <- GRanges(seqnames=rep(chromosome(x), length(ranges)),
                ranges=ranges,
                strand=rep(strand(x), length(ranges)))
  values(gr) <- DataFrame(exon.anno=exon.anno)
  gr <- gr[order(start(gr))]
  gr
})

############################################################## Simple Accessors
setGeneric("isProteinCoding", function(x, ...) standardGeneric("isProteinCoding"))
setMethod("isProteinCoding", c(x="GFGene"),
function(x, ...) {
  xcript.names <- names(transcripts(x))
  cds.names <- names(cds(x))
  if (is.null(cds.names)) {
    cds.names <- character()
  }
  ispc <- xcript.names %in% cds.names
  names(ispc) <- xcript.names
  ispc
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
