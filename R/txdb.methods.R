txdbConn <- GenomicFeatures:::txdbConn  ## maped for interactive dev purposes

setGeneric("getMetadata", function(x, ...) standardGeneric("getMetadata"))
setMethod("getMetadata", c(x="TranscriptDb"),
function(x, what, ...) {
  query <- "SELECT value FROM metadata WHERE name='%s';"
  dbGetQuery(txdbConn(x), sprintf(query, what))$value
})


setGeneric("setMetadata", function(x, ...) standardGeneric("setMetadata"))
setMethod("setMetadata", c(x="TranscriptDb"),
function(x, what, value, ...) {
  query <- "INSERT OR REPLACE INTO metadata (name,value) VALUES ('%s', '%s');"
  invisible(dbGetQuery(txdbConn(x), sprintf(query, what, value)))
})

setGeneric("getBsGenome", function(x, ...) standardGeneric("getBsGenome"))
setMethod("getBsGenome", c(x="TranscriptDb"),
function(x, ...) {
  getBsGenomeFromVersion(genome(x))
})

setGeneric("genome", function(x, ...) standardGeneric("genome"))
setMethod("genome", c(x="TranscriptDb"),
function(x, ...) {
  getMetadata(x, 'Genome')
})

setGeneric("dataSource", function(x, ...) standardGeneric("dataSource"))
setMethod("dataSource", c(x="TranscriptDb"),
function(x, ...) {
  getMetadata(x, 'Data source')
})

setGeneric("annotationSource", function(x, ...) standardGeneric("annotationSource"))
setMethod("annotationSource", c(x="TranscriptDb"),
function(x, ...) {
  getMetadata(x, 'UCSC Table')
})
