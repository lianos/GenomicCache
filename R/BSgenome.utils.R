##' Returns the genome object for the appropriate BSgenome library
getBsGenomeFromVersion <- function(genome='hg18') {
  ## BSgenome.Celegans.UCSC.ce2
  ## BSgenome.Hsapiens.UCSC.hg18
  ## BSgenome.Mmusculus.UCSC.mm9
  lib.name <- 'BSgenome.:bioc.genome:.UCSC.:genome:'
  
  bioc.genome <- switch(substring(genome, 1, 2),
    hg="Hsapiens",
    mm="Mmusculus",
    sa="Scerevisiae",
    stop("Unknown Genome Version: ", genome)
  )
  
  lib.name <- gsub(':genome:', genome, lib.name)
  lib.name <- gsub(':bioc.genome:', bioc.genome, lib.name)

  suppressPackageStartupMessages({
    found <- require(lib.name, character.only=TRUE)
  })
  
  if (!found) {
    stop(lib.name, " package required.")
  }
  
  get(bioc.genome)
}

##' Gets the annotation library name for the given genome version
getAnnotationLibraryName <- function(genome='hg18') {  
  lib.name <- 'org.:abbrev:.eg.db'
  abbrev <- switch(substring(genome, 1, 2),
                   hg="Hs",
                   mm="Mm",
                   sa="Sc",
                   stop("Unknown genome version: ", version))
  lib.name <- gsub(':abbrev:', abbrev, lib.name)
  lib.name
}

getEgAnnotationMapFromVersion <- function(what, version='hg18') {
  what <- toupper(what)
  valid.maps <- c(
    "CHRLOC",     "ENSEMBLTRANS2EG",  "MAP",        "PATH2EG",    "SYMBOL",     "UCSCKG",
    "CHRLOCEND",  "ENZYME",           "MAP2EG",     "PFAM",       "SYMBOL2EG",  "UNIGENE",
    "ACCNUM",     "ENSEMBL",          "ENZYME2EG",  "MAPCOUNTS",  "PMID",       "UNIGENE2EG",
    "ACCNUM2EG",  "ENSEMBL2EG",       "GENENAME",   "OMIM",       "PMID2EG",    "UNIPROT",
    "ALIAS2EG",   "ENSEMBLPROT",      "GO",         "OMIM2EG",    "PROSITE",
    "CHR",        "ENSEMBLPROT2EG",   "GO2ALLEGS",  "ORGANISM",   "REFSEQ",
    "CHRLENGTHS", "ENSEMBLTRANS",     "GO2EG",      "PATH",       "REFSEQ2EG"
  )
  if (!what %in% valid.maps) {
    stop("Illegal map name for org.XX.eg.db's : ", what)
  }

  lib.name <- getAnnotationLibraryName(version)
  suppressPackageStartupMessages({
    found <- require(lib.name, character.only=TRUE)
  })
  if (!found) {
    stop("This functionality requires the ", lib.name, " package")
  }
  
  var.name <- paste(gsub('\\.db', '', lib.name), what, sep="")
  get(var.name)
}
