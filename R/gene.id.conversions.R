## These functions only support passing in one symbol/id/etc at a time
map.get <- function(what, map) {
  if (length(what) == 1) {
    return(map[[what]])
  } else {
    return(mget(what, map, ifnotfound=NA))
  }
}

getEntrezIdFromSymbol <- function(txdb, symbol) {
  symbol2eg <- getEgAnnotationMapFromVersion("SYMBOL2EG", genome(txdb))
  entrez.id <- symbol2eg[[symbol]]
  
  if (is.null(entrez.id)) {
    alias2eg <- getEgAnnotationMapFromVersion("ALIAS2EG", genome(txdb))
    entrez.id <- alias2eg[[symbol]]    
  }
  
  if (is.null(entrez.id)) {
    stop("Unknown gene symbol: ", symbol)
  }
  
  ## if (length(entrez.id) != length(symbol)) {
  ##   warning("More than one entrez id for ", symbol)
  ## }
  
  entrez.id
}

getEntrezIdFromTranscriptId <- function(txdb, id) {
  asource <- annotationSource(txdb)
  if (asource == 'ensGene') {
    map <- revmap(getEgAnnotationMapFromVersion('ENSEMBLTRANS', genome(txdb)))
  } else if (asource == 'refGene') {
    map <- getEgAnnotationMapFromVersion('REFSEQ2EG', genome(txdb))
  } else {
    stop("Unknown annotation source (UCSC Table): ", asource)
  }
  x <- mget(id, map, ifnotfound=NA)
  unk <- which(sapply(x, is.na))
  if (length(unk) > 0) {
    warning("Unknown transcript ids: ", paste(id[unk], collapse=","))
    x[unk] <- NULL
  }
  if (length(id) == 1) {
    x <- unlist(x)
  }
  x
}

getEntrezIdFromGeneId <- function(txdb, id) {
  asource <- annotationSource(txdb)
  if (asource != 'ensGene') {
    stop("Can only reteive entrez id from a gene id using ensembl db")
  }
  map <- getEgAnnotationMapFromVersion('ENSEMBL2EG', genome(txdb))
  x <- mget(id, map, ifnotfound=NA)
  unk <- which(sapply(x, is.na))
  if (length(unk) > 0) {
    warning("Unknown transcript ids: ", paste(id[unk], collapse=","))
    x[unk] <- NULL
  }
  if (length(id) == 1) {
    x <- unlist(x)
  }
  x
}

getTranscriptIdFromEntrezId <- function(txdb, id) {
  asource <- annotationSource(txdb)
  if (asource == 'ensGene') {
    map <- getEgAnnotationMapFromVersion('ENSEMBLTRANS', genome(txdb))
  } else if (asource == 'refGene') {
    map <- revmap(getEgAnnotationMapFromVersion('REFSEQ2EG', genome(txdb)))
  } else {
    stop("Unknown annotation source (UCSC Table): ", asource)
  }
  x <- mget(id, map, ifnotfound=NA)
  unk <- which(sapply(x, is.na))
  if (length(unk) > 0) {
    warning("Unknown entrez ids: ", paste(id[unk], collapse=","))
    x[unk] <- NULL
  }
  
  if (asource == 'refGene') {
    ## Remove protein refseq IDs
    x <- lapply(x, function(vals) {
      vals[substring(vals, 1, 3) != "NP_"]
    })
  }
  
  if (length(id) == 1) {
    x <- unlist(x)
  }
  
  x
}

getSymbolFromEntrezId <- function(txdb, id) {
  map <- getEgAnnotationMapFromVersion('SYMBOL', genome(txdb))
  x <- mget(id, map, ifnotfound=NA)
  unk <- which(sapply(x, is.na))
  if (length(unk) > 0) {
    warning("Unknown entrez ids: ", paste(id[unk], collapse=","))
    x[unk] <- NULL
  }
  if (length(id) == 1) {
    x <- unlist(x)
  }
  x
}
# getGeneIdFromEntrezId <- function(txdb, entrez) {
#   asource <- annotationSource(txdb)
#   if (asource == 'refGene') {
#     return(get)
#   }
#   name <- switch(asource,
#                  ensGene="ENSEMBL2EG",
#                  refGene="REFSEQ2EG",
#                  stop("Unknown annotation source (UCSC Table): ", asource))
#   e2id <- revmap(getEgAnnotationMapFromVersion(name, genome(txdb)))
#   gene.id <- e2id[[entrez]]
#   if (length(gene.id) != length(entrez)) {
#     warning("Multiple matches to entrez id ", entrez)
#   }
#   gene.id
# }

###############################################################################
## Old

##' Uses the ID from the annotation type (ensembl, refseq, etc.), such as
##' ENSGXXXXXXXX to get the "Symbol" people use
getSymbolFromGeneId <- function(txdb, id) {
  asource <- annotationSource(txdb)
  name <- switch(asource,
                 ensGene="ENSEMBL2EG",
                 stop("Unknown annotation source (UCSC Table): ", asource))
  id2eg <- getEgAnnotationMapFromVersion(name, genome(txdb))
  eg2symbol <- revmap(getEgAnnotationMapFromVersion("SYMBOL2EG", genome(txdb)))
  entrez.id <- id2eg[[id]]
  if (is.null(entrez.id)) {
    stop("Unknown Gene ID: ", id, " for ", asource, ".")
  }
  if (length(entrez.id) != length(id)) {
    warning("More than one entrez id for ", id)
  }

  eg2symbol[[entrez.id[1]]]
}


##' Returns the appropriate Gene ID for the given symbol. For example, if the
##' txdb was created from ensembl annotations, the gene ids are ENSGxxxxxx
getGeneIdFromSymbol <- function(txdb, symbol) {
  getEntrezIdFromSymbol(txdb, symbol)
  getGeneIdFromEntrezId(txdb, entrez.id)
}

getGeneIdFromTranscriptId <- function(txdb, id) {
  asource <- annotationSource(txdb)
  name <- switch(asource,
                 ensGene="ENSEMBLTRANS",
                 refGene="REFSEQ2EG")
  xcript2eg <- revmap(getGeneIdFromTranscriptId(name, genome(txdb)))
  entrez.id <- xcript2eg[[id]]
  if (is.null(entrez.id)) {
    stop("Unknown transcript id ", id)
  }
  if (length(entrez.id) != length(id)) {
    warning("More than one entrez id maps to transcript ", id)
  }
  getGeneIdFromEntrez(txdb, entrez.id[1])
}
