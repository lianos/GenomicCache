## txdbConn <- GenomicFeatures:::txdbConn  ## maped for interactive dev purposes
txdbConn <- function(x) {
  if (inherits(x, 'GenomicCache')) {
    x <- x@.txdb
  }
  stopifnot(is(x, 'TranscriptDb'))
  x@.xData$conn
}

setMethod("genome", "TranscriptDb", function(x) {
  subset(metadata(x), name == "Genome")$value
})

