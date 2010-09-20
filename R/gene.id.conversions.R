.getAnnoPackageName <- function(from, package=NULL) {
  is.anno.package <- length(grep('^org\\..*\\.db$', from) == 1L)
  if (is.anno.package) {
    ## this is probably the package name itself
    if (!require(from, character.only=TRUE)) {
      stop("Unknown package: ", from)
    }
    from
  } else {
    ## probably the genome
    annotationPackage(from, package=package)
  }
}


###############################################################################
# getEntrezIdFromSymbol
###############################################################################
setGeneric("getEntrezIdFromSymbol",
function(x, id, anno.source=NULL, rm.unknown=TRUE) {
  standardGeneric("getEntrezIdFromSymbol")
})

setMethod("getEntrezIdFromSymbol", c(x="character"),
function(x, id, anno.source=NULL, rm.unknown) {
  x <- .getAnnoPackageName(x)
  symbol2eg <- getAnnMap("SYMBOL2EG", x)
  ids <- mget(id, symbol2eg, ifnotfound=NA)
  
  if (is.null(ids)) {
    unk <- seq_along(id)
  } else {
    unk <- which(sapply(ids, is.na))    
  }
  
  if (length(unk) > 0) {
    ## Look in alias
    alias2eg <- getAnnMap("ALIAS2EG", x)
    ids2 <- mget(id[unk], alias2eg, ifnotfound=NA)

    if (rm.unknown) {
      if (is.null(ids2)) {
        unk2 <- seq_along(ids2)
      } else {
        unk2 <- which(sapply(ids2, is.na))
      }

      if (length(unk2) > 0) {
        warning("Unknown transcript ids: ", paste(id[unk][unk2], collapse=","))
        ids2 <- ids2[-unk2]
      }
      
      ids2[unk2] <- NULL
    }
    
    ids <- c(ids[-unk], ids2)
  }
  
  if (length(id) == 1) {
    ids <- unlist(ids)
  }
  
  ids
})


setMethod("getEntrezIdFromSymbol", c(x="TranscriptDb"),
function(x, id, rm.unknown) {
  genome <- subset(metadata(txdb), name == "Genome")$value
  getEntrezIdFromSymbol(genome, id, rm.unknown)
})

setMethod("getEntrezIdFromSymbol", c(x="GenomicCache"),
function(x, id, rm.unknown) {
  getEntrezIdFromSymbol(genome(x), id, rm.unknown, genome=genome)
})


###############################################################################
# getEntrezIdFromTranscriptId
###############################################################################
setGeneric("getEntrezIdFromTranscriptId",
function(x, id, anno.source, rm.unknown=TRUE) {
  standardGeneric("getEntrezIdFromTranscriptId")
})

setMethod("getEntrezIdFromTranscriptId", c(x="character"),
function(x, id, anno.source, rm.unknown) {
  x <- .getAnnoPackageName(x)
  
  if (anno.source == 'ensGene') {
    map <- revmap(getAnnMap('ENSEMBLTRANS', x))
  } else if (anno.source == 'refGene') {
    map <- getAnnMap('REFSEQ2EG', x)
  } else {
    stop("Unknown annotation source (UCSC Table): ", anno.source)
  }
  
  ids <- mget(id, map, ifnotfound=NA)
  
  if (rm.unknown) {
    if (is.null(ids)) {
      unk <- seq_along(id)
    } else {
      unk <- which(sapply(ids, is.na))    
    }
    
    if (length(unk) > 0) {
      warning("Unknown transcript ids: ", paste(id[unk], collapse=","))
      ids[unk] <- NULL
    }
  }

  if (length(id) == 1) {
    ids <- unlist(ids)
  }

  ids
})

setMethod("getEntrezIdFromTranscriptId", c(x="GenomicCache"),
function(x, id, anno.source=annotationSource(x), rm.unknown) {
  getEntrezIdFromTranscriptId(genome(x), id, anno.source, rm.unknown)
})

setMethod("getEntrezIdFromTranscriptId", c(x="TranscriptDb"),
function(x, id, anno.source, rm.unknown) {
  if (missing(anno.source)) {
    anno.source <- subset(metadata(object@.txdb), name == "UCSC Table")$value
  }
  genome <- subset(metadata(txdb), name == "Genome")$value
  getEntrezIdFromTranscriptId(genome, id, anno.source, rm.unknown)
})


###############################################################################
# getEntrezIdFromGeneId
###############################################################################
setGeneric("getEntrezIdFromGeneId",
function(x, id, anno.source, rm.unknown=TRUE) {
  standardGeneric("getEntrezIdFromGeneId")
})

setMethod("getEntrezIdFromGeneId", c(x="character"),
function(x, id, anno.source, rm.unknown) {
  x <- .getAnnoPackageName(x)
  
  if (anno.source == 'refGene') {
    warning("Assuming gne id is its symbol for RefSeq")
    ids <- getEntrezIdFromSymbol(x, id, anno.source, rm.unknown)
  } else if (anno.source == 'ensGene') {
    map <- getAnnMap('ENSEMBL2EG', x)
    ids <- mget(id, map, ifnotfound=NA)
  }

  if (rm.unknown) {
    if (is.null(ids)) {
      unk <- seq_along(id)
    } else {
      unk <- which(sapply(ids, is.na))    
    }
    
    if (length(unk) > 0) {
      warning("Unknown transcript ids: ", paste(id[unk], collapse=","))
      ids[unk] <- NULL
    }
  }
  
  if (length(ids) == 1) {
    ids <- unlist(ids)
  }
  
  ids
})

setMethod("getEntrezIdFromGeneId", c(x="GenomicCache"),
function(x, id, anno.source=annotationSource(x), rm.unknown) {
  getEntrezIdFromGeneId(genome(x), id, anno.source, rm.unknown)
})

setMethod("getEntrezIdFromGeneId", c(x="TranscriptDb"),
function(x, id, anno.source, rm.unknown) {
  if (missing(anno.source)) {
    anno.source <- subset(metadata(object@.txdb), name == "UCSC Table")$value
  }
  genome <- subset(metadata(txdb), name == "Genome")$value
  getEntrezIdFromGeneId(genome, id, anno.source, rm.unknown)
})


###############################################################################
# getTranscriptIdFromEntrezId
###############################################################################
setGeneric("getTranscriptIdFromEntrezId",
function(x, id, anno.source, rm.unknown=TRUE) {
  standardGeneric("getTranscriptIdFromEntrezId")
})

setMethod("getTranscriptIdFromEntrezId", c(x="character"),
function(x, id, anno.source, rm.unknown) {
  x <- .getAnnoPackageName(x)
  
  if (anno.source == 'ensGene') {
    map <- getEgAnnotationMapFromVersion('ENSEMBLTRANS', genome(x))
  } else if (anno.source == 'refGene') {
    map <- revmap(getEgAnnotationMapFromVersion('REFSEQ2EG', genome(x)))
  } else {
    stop("Unknown annotation source (UCSC Table): ", asource)
  }

  ids <- mget(id, map, ifnotfound=NA)
  if (rm.unknown) {
    if (is.null(ids)) {
      unk <- seq_along(id)
    } else {
      unk <- which(sapply(ids, is.na))    
    }
    
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
  }
  
  if (length(ids) == 1) {
    ids <- unlist(ids)
  }
  
  ids
})

setMethod("getTranscriptIdFromEntrezId", c(x="GenomicCache"),
function(x, id, rm.unknown) {
  getTranscriptIdFromEntrezId(x@.txdb, id)
})


###############################################################################
# getSymbolFromEntrezId
###############################################################################
setGeneric("getSymbolFromEntrezId", function(x, id, rm.unknown=TRUE) {
  standardGeneric("getSymbolFromEntrezId")
})

setMethod("getSymbolFromEntrezId", c(x="TranscriptDb"),
function(x, id, rm.unknown) {
  map <- getEgAnnotationMapFromVersion('SYMBOL', genome(x))
  ids <- mget(id, map, ifnotfound=NA)
  
  if (rm.unknown) {
    if (is.null(ids)) {
      unk <- seq_along(id)
    } else {
      unk <- which(sapply(ids, is.na))    
    }
    
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

setMethod("getSymbolFromEntrezId", c(x="GenomicCache"),
function(x, id, rm.unknown) {
  getSymbolFromEntrezId(x@.txdb, id, rm.unknown)
})


###############################################################################
# getGeneIdFromEntrezId
###############################################################################
setGeneric("getGeneIdFromEntrezId", function(x, id, rm.unknown=TRUE) {
 standardGeneric("getGeneIdFromEntrezId") 
})

setMethod("getGeneIdFromEntrezId", c(x="TranscriptDb"),
function(x, id, rm.unknown) {
  asource <- annotationSource(x)
  if (asource == 'refGene') {
    warning("RefSeq doesn't have 'gene ids' since their IDs are ",
            "returning the gene symbol")
   ids <- getSymbolFromEntrezId(x, id)
  } else if (asource == 'ensGene') {
    map <- revmap(getEgAnnotationMapFromVersion('ENSEMBL2EG', genome(x)))
    ids <- mget(id, map, ifnotfound=NA)
  }

  if (rm.unknown) {
    if (is.null(ids)) {
      unk <- seq_along(id)
    } else {
      unk <- which(sapply(ids, is.na))    
    }
    
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
function(x, id, rm.unknown) {
  getGeneIdFromEntrezId(x@.txdb, id, rm.unknown)
})

