## S4 Definitions I feel should be added to base GenomicFeatures package
setMethod("metadata", c(x="TranscriptDb"),
function(x, what=NULL, ...) {
  if (is.null(what)) {
    dbReadTable(GenomicFeatures:::txdbConn(x), "metadata")
  } else {
    query <- "SELECT value FROM metadata WHERE name='%s';"
    dbGetQuery(GenomicFeatures:::txdbConn(x), sprintf(query, what))$value
  }
})

setMethod("getBsGenome", c(x="TranscriptDb"),
function(x, ...) {
  getBsGenomeFromVersion(genome(x))
})

##' Returns the version of the genome (hg18, mm9, etc.)
setMethod("genome", c(x="TranscriptDb"),
function(x, ...) {
  metadata(x, 'Genome')
})

##' The Data Source -- "UCSC"
setMethod("dataSource", c(x="TranscriptDb"),
function(x, ...) {
  metadata(x, 'Data source')
})

##' The data table that was pulled from the data source, ensGene, refGene
setMethod("annotationSource", c(x="TranscriptDb"),
function(x, ...) {
  metadata(x, 'UCSC Table')
})

## if (!isGeneric("transcripts")) {
##   if (is.function("transcripts")) {
##     fun <- transcripts
##   } else {
##     fun <- function(x, ...) standardGeneric("transcripts")
##   }
##   setGeneric("transcripts", fun)
## }
