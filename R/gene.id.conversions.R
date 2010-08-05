setGeneric("getEntrezIdFromSymbol", function(x, symbol) {
  standardGeneric("getEntrezIdFromSymbol")
})
setMethod("getEntrezIdFromSymbol", c(x="TranscriptDb"),
function(x, symbol) {
  symbol2eg <- getEgAnnotationMapFromVersion("SYMBOL2EG", genome(x))
  entrez.id <- symbol2eg[[symbol]]
  
  if (is.null(entrez.id)) {
    alias2eg <- getEgAnnotationMapFromVersion("ALIAS2EG", genome(x))
    entrez.id <- alias2eg[[symbol]]    
  }
  
  if (is.null(entrez.id)) {
    stop("Unknown gene symbol: ", symbol)
  }
  ## if (length(entrez.id) != length(symbol)) {
  ##   warning("More than one entrez id for ", symbol)
  ## }
  entrez.id
})
setMethod("getEntrezIdFromSymbol", c(x="GenomicCache"),
function(x, symbol) {
  getEntrezIdFromSymbol(x@.txdb, symbol)
})

setGeneric("getEntrezIdFromTranscriptId", function(x, id) {
  standardGeneric("getEntrezIdFromTranscriptId")
})
setMethod("getEntrezIdFromTranscriptId", c(x="TranscriptDb"),
function(x, id) {
  asource <- annotationSource(x)
  if (asource == 'ensGene') {
    map <- revmap(getEgAnnotationMapFromVersion('ENSEMBLTRANS', genome(x)))
  } else if (asource == 'refGene') {
    map <- getEgAnnotationMapFromVersion('REFSEQ2EG', genome(x))
  } else {
    stop("Unknown annotation source (UCSC Table): ", asource)
  }
  ids <- mget(id, map, ifnotfound=NA)
  unk <- which(sapply(ids, is.na))
  if (length(unk) > 0) {
    warning("Unknown transcript ids: ", paste(id[unk], collapse=","))
    ids[unk] <- NULL
  }
  if (length(id) == 1) {
    ids <- unlist(ids)
  }
  ids
})
setMethod("getEntrezIdFromTranscriptId", c(x="GenomicCache"),
function(x, id) {
  getEntrezIdFromTranscriptId(x@.txdb, id)
})


setGeneric("getEntrezIdFromGeneId", function(x, id) {
  standardGeneric("getEntrezIdFromGeneId")
})
setMethod("getEntrezIdFromGeneId", c(x="TranscriptDb"),
function(x, id) {
  asource <- annotationSource(x)
  if (asource == 'refGene') {
    warning("Assuming gne id is its symbol for RefSeq")
    ids <- getEntrezIdFromSymbol(x, id)
  } else if (asource == 'ensGene') {
    map <- getEgAnnotationMapFromVersion('ENSEMBL2EG', genome(x))
    ids <- mget(id, map, ifnotfound=NA)
  }
  unk <- which(sapply(ids, is.na))
  if (length(unk) > 0) {
    warning("Unknown transcript ids: ", paste(id[unk], collapse=","))
    ids[unk] <- NULL
  }
  if (length(ids) == 1) {
    ids <- unlist(ids)
  }
  ids
})
setMethod("getEntrezIdFromGeneId", c(x="GenomicCache"),
function(x, id) {
  getEntrezIdFromGeneId(x@.txdb, id)
})


setGeneric("getTranscriptIdFromEntrezId", function(x, id) {
  standardGeneric("getTranscriptIdFromEntrezId")
})
setMethod("getTranscriptIdFromEntrezId", c(x="TranscriptDb"),
function(x, id) {
  asource <- annotationSource(x)
  if (asource == 'ensGene') {
    map <- getEgAnnotationMapFromVersion('ENSEMBLTRANS', genome(x))
  } else if (asource == 'refGene') {
    map <- revmap(getEgAnnotationMapFromVersion('REFSEQ2EG', genome(x)))
  } else {
    stop("Unknown annotation source (UCSC Table): ", asource)
  }
  ids <- mget(id, map, ifnotfound=NA)
  unk <- which(sapply(ids, is.na))
  if (length(unk) > 0) {
    warning("Unknown entrez ids: ", paste(id[unk], collapse=","))
    ids[unk] <- NULL
  }
  
  if (asource == 'refGene') {
    ## Remove protein refseq IDs
    ids <- lapply(ids, function(vals) {
      vals[substring(vals, 1, 3) != "NP_"]
    })
  }
  
  if (length(ids) == 1) {
    ids <- unlist(ids)
  }
  
  ids
})
setMethod("getTranscriptIdFromEntrezId", c(x="GenomicCache"),
function(x, id) {
  getTranscriptIdFromEntrezId(x@.txdb, id)
})


setGeneric("getSymbolFromEntrezId", function(x, id) {
  standardGeneric("getSymbolFromEntrezId")
})
setMethod("getSymbolFromEntrezId", c(x="TranscriptDb"),
function(x, id) {
  map <- getEgAnnotationMapFromVersion('SYMBOL', genome(x))
  ids <- mget(id, map, ifnotfound=NA)
  unk <- which(sapply(ids, is.na))
  if (length(unk) > 0) {
    warning("Unknown entrez ids: ", paste(id[unk], collapse=","))
    ids[unk] <- NULL
  }
  if (length(ids) == 1) {
    ids <- unlist(ids)
  }
  ids
})
setMethod("getSymbolFromEntrezId", c(x="GenomicCache"),
function(x, id) {
  getSymbolFromEntrezId(x@.txdb, id)
})

setGeneric("getGeneIdFromEntrezId", function(x, id) {
 standardGeneric("getGeneIdFromEntrezId") 
})
setMethod("getGeneIdFromEntrezId", c(x="TranscriptDb"),
function(x, id) {
  asource <- annotationSource(x)
  if (asource == 'refGene') {
    warning("RefSeq doesn't have 'gene ids' since their IDs are ",
            "returning the gene symbol")
   ids <- getSymbolFromEntrezId(x, id)
  } else if (asource == 'ensGene') {
    map <- revmap(getEgAnnotationMapFromVersion('ENSEMBL2EG', genome(x)))
    ids <- mget(id, map, ifnotfound=NA)
    unk <- which(sapply(ids, is.na))
    if (length(unk) > 0) {
      warning("Unknown entrez ids: ", paste(id[unk], collapse=","))
      ids[unk] <- NULL
    }
  }
  
  if (length(ids) == 1) {
    ids <- unlist(ids)
  }
  ids
})
setMethod("getGeneIdFromEntrezId", c(x="GenomicCache"),
function(x, id) {
  getGeneIdFromEntrezId(x@.txdb, id)
})


###############################################################################
## Old

##' Uses the ID from the annotation type (ensembl, refseq, etc.), such as
##' ENSGXXXXXXXX to get the "Symbol" people use
## getSymbolFromGeneId <- function(txdb, id) {
##   asource <- annotationSource(txdb)
##   name <- switch(asource,
##                  ensGene="ENSEMBL2EG",
##                  stop("Unknown annotation source (UCSC Table): ", asource))
##   id2eg <- getEgAnnotationMapFromVersion(name, genome(txdb))
##   eg2symbol <- revmap(getEgAnnotationMapFromVersion("SYMBOL2EG", genome(txdb)))
##   entrez.id <- id2eg[[id]]
##   if (is.null(entrez.id)) {
##     stop("Unknown Gene ID: ", id, " for ", asource, ".")
##   }
##   if (length(entrez.id) != length(id)) {
##     warning("More than one entrez id for ", id)
##   }

##   eg2symbol[[entrez.id[1]]]
## }


## ##' Returns the appropriate Gene ID for the given symbol. For example, if the
## ##' txdb was created from ensembl annotations, the gene ids are ENSGxxxxxx
## getGeneIdFromSymbol <- function(txdb, symbol) {
##   getEntrezIdFromSymbol(txdb, symbol)
##   getGeneIdFromEntrezId(txdb, entrez.id)
## }

## getGeneIdFromTranscriptId <- function(txdb, id) {
##   asource <- annotationSource(txdb)
##   name <- switch(asource,
##                  ensGene="ENSEMBLTRANS",
##                  refGene="REFSEQ2EG")
##   xcript2eg <- revmap(getGeneIdFromTranscriptId(name, genome(txdb)))
##   entrez.id <- xcript2eg[[id]]
##   if (is.null(entrez.id)) {
##     stop("Unknown transcript id ", id)
##   }
##   if (length(entrez.id) != length(id)) {
##     warning("More than one entrez id maps to transcript ", id)
##   }
##   getGeneIdFromEntrez(txdb, entrez.id[1])
## }
