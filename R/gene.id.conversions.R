setGeneric("getEntrezIdFromSymbol", function(x, id, rm.unknown=TRUE) {
  standardGeneric("getEntrezIdFromSymbol")
})

setMethod("getEntrezIdFromSymbol", c(x="TranscriptDb"),
function(x, id, rm.unknown) {
  symbol2eg <- getEgAnnotationMapFromVersion("SYMBOL2EG", genome(x))
  ids <- mget(id, symbol2eg, ifnotfound=NA)
  
  if (is.null(ids)) {
    unk <- seq_along(id)
  } else {
    unk <- which(sapply(ids, is.na))    
  }
  
  if (length(unk) > 0) {
    ## Look in alias
    alias2eg <- getEgAnnotationMapFromVersion("ALIAS2EG", genome(x))
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

setMethod("getEntrezIdFromSymbol", c(x="GenomicCache"),
function(x, id, rm.unknown) {
  getEntrezIdFromSymbol(x@.txdb, id, rm.unknown)
})

setGeneric("getEntrezIdFromTranscriptId", function(x, id, rm.unknown=TRUE) {
  standardGeneric("getEntrezIdFromTranscriptId")
})
setMethod("getEntrezIdFromTranscriptId", c(x="TranscriptDb"),
function(x, id, rm.unknown) {
  asource <- annotationSource(x)
  if (asource == 'ensGene') {
    map <- revmap(getEgAnnotationMapFromVersion('ENSEMBLTRANS', genome(x)))
  } else if (asource == 'refGene') {
    map <- getEgAnnotationMapFromVersion('REFSEQ2EG', genome(x))
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
function(x, id, rm.unknown) {
  getEntrezIdFromTranscriptId(x@.txdb, id, rm.unknown)
})


setGeneric("getEntrezIdFromGeneId", function(x, id, rm.unknown=TRUE) {
  standardGeneric("getEntrezIdFromGeneId")
})
setMethod("getEntrezIdFromGeneId", c(x="TranscriptDb"),
function(x, id, rm.unknown) {
  asource <- annotationSource(x)
  if (asource == 'refGene') {
    warning("Assuming gne id is its symbol for RefSeq")
    ids <- getEntrezIdFromSymbol(x, id)
  } else if (asource == 'ensGene') {
    map <- getEgAnnotationMapFromVersion('ENSEMBL2EG', genome(x))
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
function(x, id, rm.unknown) {
  getEntrezIdFromGeneId(x@.txdb, id, rm.unknown)
})

setGeneric("getTranscriptIdFromEntrezId", function(x, id, rm.unknown=TRUE) {
  standardGeneric("getTranscriptIdFromEntrezId")
})

setMethod("getTranscriptIdFromEntrezId", c(x="TranscriptDb"),
function(x, id, rm.unknown) {
  asource <- annotationSource(x)
  if (asource == 'ensGene') {
    map <- getEgAnnotationMapFromVersion('ENSEMBLTRANS', genome(x))
  } else if (asource == 'refGene') {
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

