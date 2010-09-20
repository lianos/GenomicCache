## S4 Definitions I feel should be added to base GenomicFeatures package
# setMethod("metadata", c(x="TranscriptDb"),
# function(x, what=NULL, ...) {
#   if (is.null(what)) {
#     dbReadTable(GenomicFeatures:::txdbConn(x), "metadata")
#   } else {
#     query <- "SELECT value FROM metadata WHERE name='%s';"
#     dbGetQuery(GenomicFeatures:::txdbConn(x), sprintf(query, what))$value
#   }
# })

setMethod("getBsGenome", c(x="TranscriptDb"),
function(x, ...) {
  getBsGenomeFromVersion(genome(x))
})

## if (!isGeneric("transcripts")) {
##   if (is.function("transcripts")) {
##     fun <- transcripts
##   } else {
##     fun <- function(x, ...) standardGeneric("transcripts")
##   }
##   setGeneric("transcripts", fun)
## }
