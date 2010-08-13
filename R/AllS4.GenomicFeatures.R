## S4 Definitions I feel should be added to base GenomicFeatures package
setMethod("getMetadata", c(x="TranscriptDb"),
function(x, what, ...) {
  query <- "SELECT value FROM metadata WHERE name='%s';"
  dbGetQuery(txdbConn(x), sprintf(query, what))$value
})

setMethod("setMetadata", c(x="TranscriptDb"),
function(x, what, value, ...) {
  query <- "INSERT OR REPLACE INTO metadata (name,value) VALUES ('%s', '%s');"
  invisible(dbGetQuery(txdbConn(x), sprintf(query, what, value)))
})

setMethod("getBsGenome", c(x="TranscriptDb"),
function(x, ...) {
  getBsGenomeFromVersion(genome(x))
})

##' Returns the version of the genome (hg18, mm9, etc.)
setMethod("genome", c(x="TranscriptDb"),
function(x, ...) {
  getMetadata(x, 'Genome')
})

##' The Data Source -- "UCSC"
setMethod("dataSource", c(x="TranscriptDb"),
function(x, ...) {
  getMetadata(x, 'Data source')
})

##' The data table that was pulled from the data source, ensGene, refGene
setMethod("annotationSource", c(x="TranscriptDb"),
function(x, ...) {
  getMetadata(x, 'UCSC Table')
})

setMethod("transcripts", c(x="TranscriptDb"),
function(x, vals=NULL, columns=c("tx_id", "tx_name")) {
  GenomicFeatures::transcripts(x, vals, columns)
})

setGeneric("exonsBy", function(x, by=c('tx', 'gene')) {
  standardGeneric("exonsBy")
})
setMethod("exonsBy", c(x="TranscriptDb"),
function(x, by) {
  GenomicFeatures::exonsBy(x, by)
})

setMethod("exons", c(x="TranscriptDb"),
function(x, vals) {
  GenomicFeatures::exons(x, vals)
})

setGeneric("cdsBy", function(x, by=c('tx', 'gene')) {
  standardGeneric("cdsBy")
})
setMethod("cdsBy", c(x="TranscriptDb"),
function(x, by) {
  GenomicFeatures::cdsBy(x, by)
})

setMethod("cds", c(x="TranscriptDb"),
function(x, vals=NULL) {
  GenomicFeatures::cds(x, vals)
})

setGeneric("threeUTRsByTranscript", function(x) {
  standardGeneric("threeUTRsByTranscript")
})
setMethod("threeUTRsByTranscript", c(x="TranscriptDb"),
function(x) {
  GenomicFeatures::threeUTRsByTranscript(x)
})

setGeneric("fiveUTRsByTranscript", function(x) {
  standardGeneric("fiveUTRsByTranscript")
})
setMethod("fiveUTRsByTranscript", c(x="TranscriptDb"),
function(x) {
  GenomicFeatures::fiveUTRsByTranscript(x)
})

## if (!isGeneric("transcripts")) {
##   if (is.function("transcripts")) {
##     fun <- transcripts
##   } else {
##     fun <- function(x, ...) standardGeneric("transcripts")
##   }
##   setGeneric("transcripts", fun)
## }
