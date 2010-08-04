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

getGeneIdFromEntrez <- function(txdb, entrez) {
  asource <- annotationSource(txdb)
  name <- switch(asource,
                 ensGene="ENSEMBL2EG",
                 stop("Unknown annotation source (UCSC Table): ", asource))
  e2id <- revmap(getEgAnnotationMapFromVersion(name, genome(txdb)))
  gene.id <- e2id[[entrez]]
  if (length(gene.id) != length(entrez)) {
    warning("Multiple matches to entrez id ", entrez)
  }
  gene.id
}

##' Returns the appropriate Gene ID for the given symbol. For example, if the
##' txdb was created from ensembl annotations, the gene ids are ENSGxxxxxx
getGeneIdFromSymbol <- function(txdb, symbol) {
  asource <- annotationSource(txdb)
  symbol2eg <- getEgAnnotationMapFromVersion("SYMBOL2EG", genome(txdb))
  entrez.id <- symbol2eg[[symbol]]
  
  if (is.null(entrez.id)) {
    alias2eg <- getEgAnnotationMapFromVersion("ALIAS2EG", genome(txdb))
    entrez.id <- alias2eg[[symbol]]
    if (is.null(entrez.id)) {
      stop("Unknown gene symbol: ", symbol)
    }
  }
  if (length(entrez.id) != length(symbol)) {
    warning("More than one entrez id for ", symbol)
  }

  getGeneIdFromEntrez(txdb, entrez.id)
}

getGeneIdFromTranscriptId <- function(txdb, id) {
  asource <- annotationSource(txdb)
  name <- switch(asource,
                 ensGene="org.Hs.egENSEMBLTRANS")
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
