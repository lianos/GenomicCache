## These funcions do "cross reference" look ups for transcript ids, symbol
## names gene ids, etc. The "base" function for each receives a character
## value for `x` which should be the genome name (ie. hg18, mm9, etc.)
## 
## All functions have the same signature:
## @param x The name of the genome the query is against
## @param id The identifiers that need to be translated.
## @param anno.source What type of annotation are we dealing with? Valid values
## include refGene, ensGene, aceGene (?) ...
## @param rm.unknown A logical specifying whether cross references which could
## not be identified should be returned as NA (\code{FALSE}) or removed
## (\code{TRUE})

###############################################################################
## getEntrezIdFromSymbol
###############################################################################

##' Fetch the entrez id given a gene symbol from a given genome.
##'
##' @importFrom annotate getAnnMap
##' @exportMethod getEntrezIdFromSymbol
##' @rdname getEntrezIdFromSymbol-methods
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param x A character vector specifying the genome (ie. "hg18"), a
##' \code{linkS4class{TranscriptDb}} object, or a
##' \code{\linkS4class{GenomicCache}} that the genome information can be
##' extracted from.
##' @param id A character vector of symbols to get the id from
##' @param anno.source Explicitly specify \code{refGene}, etc.
##' @param rm.unknown If \code{TRUE} symbols with unknown entrez id's are
##' removed, otherwise \code{NA} kept for them
##'
##' @return If only one \code{id} is given, then a character vector with the
##' entrez id(s) is returned, otherwise a named list is returned with the
##' entrez id(s) as elements, and the symbol(s) as names
setGeneric("getEntrezIdFromSymbol",
function(x, id, anno.source=NULL, rm.unknown=TRUE) {
  standardGeneric("getEntrezIdFromSymbol")
})

setMethod("getEntrezIdFromSymbol", c(x="character"),
function(x, id, anno.source, rm.unknown) {
  x <- getAnnoPackageName(x)
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
  
  if (length(id) == 1L) {
    ids <- unlist(ids)
  }
  
  ids
})


setMethod("getEntrezIdFromSymbol", c(x="TranscriptDb"),
function(x, id, anno.source, rm.unknown) {
  genome <- subset(metadata(txdb), name == "Genome")$value
  getEntrezIdFromSymbol(genome, id, anno.source, rm.unknown)
})

setMethod("getEntrezIdFromSymbol", c(x="GenomicCache"),
function(x, id, anno.source, rm.unknown) {
  getEntrezIdFromSymbol(genome(x), id, annotationSource(x), rm.unknown)
})


###############################################################################
## getEntrezIdFromTranscriptId
###############################################################################
setGeneric("getEntrezIdFromTranscriptId",
function(x, id, anno.source, rm.unknown=TRUE) {
  standardGeneric("getEntrezIdFromTranscriptId")
})

setMethod("getEntrezIdFromTranscriptId", c(x="character"),
function(x, id, anno.source, rm.unknown) {
  x <- getAnnoPackageName(x)
  
  if (anno.source == 'ensGene') {
    map <- revmap(getAnnMap('ENSEMBLTRANS', x))
  } else if (anno.source == 'refGene') {
    map <- getAnnMap('REFSEQ2EG', x)
  } else if (anno.source == 'knownGene') {
    map <- revmap(getAnnMap('UCSCKG', x))
  } else if (anno.source == 'acembly') {
    warning("No transcript ID map in org.*.db Annotation maps, parsing symbol")
    symbol <- strsplit(id, '.', fixed=TRUE)[[1]][1]
    return(getEntrezIdFromSymbol(x, symbol, anno.source, rm.unknown))
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

  if (length(id) == 1L) {
    ids <- unlist(ids)
  }

  ids
})

setMethod("getEntrezIdFromTranscriptId", c(x="GenomicCache"),
function(x, id, anno.source, rm.unknown) {
  getEntrezIdFromTranscriptId(txdb(x), id, annotationSource(x), rm.unknown)
})

setMethod("getEntrezIdFromTranscriptId", c(x="TranscriptDb"),
function(x, id, anno.source=annotationSource(x), rm.unknown) {
  if (anno.source == 'ensGene') {
    return(getEntrezIdFromTranscriptId(genome(x), id, anno.source, rm.unknown))
  }
  query <- sprintf("SELECT _tx_id FROM transcript WHERE tx_name='%s'", id)
  tx.id <- dbGetQuery(txdbc(x), query)[,1]
  if (length(tx.id) == 0L) {
    stop("Unknown transcripts: ", paste(id, collapse=","))
  }
  df <- data.frame(tx.id=tx.id)
  query <- "SELECT gene_id,_tx_id FROM gene WHERE _tx_id=?"
  df <- dbGetPreparedQuery(txdbc(x), query, df)
  res <- df[["gene_id"]]
  names(res) <- df[['_tx_id']]
  res
})


###############################################################################
## getEntrezIdFromGeneId
###############################################################################
setGeneric("getEntrezIdFromGeneId",
function(x, id, anno.source, rm.unknown=TRUE) {
  standardGeneric("getEntrezIdFromGeneId")
})

setMethod("getEntrezIdFromGeneId", c(x="character"),
function(x, id, anno.source, rm.unknown) {
  x <- getAnnoPackageName(x)
  
  if (anno.source %in% c('refGene', 'knownGene', 'acembly')) {
    message("Assuming gene id is its symbol")
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
  
  if (length(ids) == 1L) {
    ids <- unlist(ids)
  }
  
  ids
})

setMethod("getEntrezIdFromGeneId", c(x="GenomicCache"),
function(x, id, anno.source, rm.unknown) {
  getEntrezIdFromGeneId(genome(x), id, annotationSource(x), rm.unknown)
})

setMethod("getEntrezIdFromGeneId", c(x="TranscriptDb"),
function(x, id, anno.source, rm.unknown) {
  meta <- metadata(x@.txdb)
  if (missing(anno.source)) {
    anno.source <- subset(meta, name == "UCSC Table")$value
  }
  genome <- subset(meta, name == "Genome")$value
  getEntrezIdFromGeneId(genome, id, anno.source, rm.unknown)
})


###############################################################################
## getTranscriptIdFromEntrezId
###############################################################################
setGeneric("getTranscriptIdFromEntrezId",
function(x, id, anno.source, rm.unknown=TRUE) {
  standardGeneric("getTranscriptIdFromEntrezId")
})

setMethod("getTranscriptIdFromEntrezId", c(x="character"),
function(x, id, anno.source, rm.unknown) {
  x <- getAnnoPackageName(x)
  
  if (anno.source == 'ensGene') {
    map <- getAnnMap('ENSEMBLTRANS', x)
  } else if (anno.source == 'refGene') {
    map <- revmap(getAnnMap('REFSEQ2EG', x))
  } else if (anno.source == 'knownGene') {
    map <- getAnnMap('UCSCKG', x)
  } else if (anno.source == 'acembly') {
    stop("No transcript map for Aceview genes in org.*.db pacakges")
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
      warning("Unknown entrez ids: ", paste(id[unk], collapse=","))
      ids[unk] <- NULL
    }
    
    if (anno.source == 'refGene') {
      ## Remove protein refseq IDs
      ids <- lapply(ids, function(vals) {
        vals[substring(vals, 1, 3) != "NP_"]
      })
    }
  }
  
  if (length(ids) == 1L) {
    ids <- unlist(ids)
  }
  
  ids
})

setMethod("getTranscriptIdFromEntrezId", c(x="TranscriptDb"),
function(x, id, anno.source=annotationSource(x), rm.unknown) {
  if (anno.source == 'ensGene') {
    return(getTranscriptIdFromEntrezId(genome(x), id, anno.source, rm.unknown))
  }
  query <- sprintf("SELECT _tx_id FROM gene WHERE gene_id='%s'", id)
  tx.id <- dbGetQuery(txdbc(x), query)[,1]
  if (length(tx.id) == 0L) {
    stop("Unknown entrez id ", id)
  }
  query <- sprintf("SELECT tx_name FROM transcript WHERE _tx_id=?")
  df <- data.frame("_tx_id"=tx.id)
  tx.name <- dbGetPreparedQuery(txdbc(x), query, df)[,1]
  tx.name
})

setMethod("getTranscriptIdFromEntrezId", c(x="GenomicCache"),
function(x, id, anno.source, rm.unknown) {
  ## getTranscriptIdFromEntrezId(genome(x), id, annotationSource(x), rm.unknown)
  getTranscriptIdFromEntrezId(txdb(x), id, annotationSource(x), rm.unknown)  
})


###############################################################################
## getSymbolFromEntrezId
###############################################################################
setGeneric("getSymbolFromEntrezId", function(x, id, anno.source, rm.unknown=TRUE) {
  standardGeneric("getSymbolFromEntrezId")
})

setMethod("getSymbolFromEntrezId", c(x="character"),
function(x, id, anno.source, rm.unknown) {
  x <- getAnnoPackageName(x)
  map <- getAnnMap('SYMBOL', x)
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
function(x, id, anno.source, rm.unknown) {
  getSymbolFromEntrezId(genome(x), id, annotationSource(x), rm.unknown)
})

setMethod("getSymbolFromEntrezId", c(x="TranscriptDb"),
function(x, id, anno.source, rm.unknown) {
  meta <- metadata(x@.txdb)
  if (missing(anno.source)) {
    anno.source <- subset(meta, name == "UCSC Table")$value
  }
  genome <- subset(meta, name == "Genome")$value
  getSymbolFromEntrezId(genome, id, anno.source, rm.unknown)
})

###############################################################################
## getGeneIdFromEntrezId
###############################################################################
setGeneric("getGeneIdFromEntrezId",
function(x, id, anno.source, rm.unknown=TRUE) {
 standardGeneric("getGeneIdFromEntrezId") 
})

setMethod("getGeneIdFromEntrezId", c(x="character"),
function(x, id, anno.source, rm.unknown) {
  x <- getAnnoPackageName(x)
  
  if (anno.source %in% c('refGene', 'knownGene', 'acembly')) {
    warning("RefSeq doesn't have 'gene ids' since their IDs are ",
            "returning the gene symbol")
   ids <- getSymbolFromEntrezId(x, id, anno.source, rm.unknown)
  } else if (anno.source == 'ensGene') {
    map <- revmap(getAnnMap('ENSEMBL2EG', x))
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
function(x, id, anno.source, rm.unknown) {
  getGeneIdFromEntrezId(genome(x), id, annotationSource(x), rm.unknown)
})

setMethod("getGeneIdFromEntrezId", c(x="TranscriptDb"),
function(x, id, anno.source, rm.unknown) {
  meta <- metadata(x@.txdb)
  if (missing(anno.source)) {
    anno.source <- subset(meta, name == "UCSC Table")$value
  }
  genome <- subset(meta, name == "Genome")$value
  getGeneIdFromEntrezId(genome, id, anno.source, rm.unknown)
})

