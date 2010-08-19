setMethod("initialize", "GenomicCache",
function(.Object, ...,
         .genome=character(),
         .path=character(),
         .txdb=NULL,
         .cache=new.env()) {
  callNextMethod(.Object,
                 .genome=.genome,
                 .path=.path,
                 .txdb=.txdb,
                 .cache=.cache,
                 ...)
})

GenomicCache <- function(path, pre.load=c('transcripts', 'exons')) {
  if (!dir.exists(path)) {
    stop("Cannot read directory: ", path)
  }
  features.path <- file.path(path, 'features') 
  if (!dir.exists(features.path)) {
    stop("Invalid GenomicCache directory -- no 'features' subdirectory found.")
  }
  
  ## Get the TranscriptDb object
  txdb.path <- list.files(features.path, 'TranscriptDb', ignore.case=TRUE,
                          full.names=TRUE)
  if (length(txdb.path) != 1) {
    stop("Need 1 and only 1 TranscriptDb object in the GenomicCache/features ",
         "directory, found: ", length(txdb.path))
  }
  txdb <- loadFeatures(txdb.path)
  attr(txdb, 'path') <- txdb.path
  
  gc <- new('GenomicCache',
            .genome=genome(txdb),
            .path=path,
            .txdb=txdb)
  
  ## Preload objects
  can.load <- c('transcripts', 'exons', 'utr3', 'utr5')
  for (what in pre.load) {
    if (!what %in% can.load) warning("Don't know how to load: ", what, "\n")
    cat("Preloading", what, "...\n")
    getFunction(what)(gc)
  }
  
  gc
}

setMethod("seqnames", c(x="GenomicCache"),
function(x) {
  seqnames(x@.txdb)
})

setMethod("chromosomes", c(x="GenomicCache"),
function(x, ...) {
  seqnames(x@.txdb)
})

setMethod("transcripts", c(x="GenomicCache"),
function(x, vals=NULL, columns=c("tx_id", "tx_name")) {
  var <- generateCacheName('transcripts', vals=vals, columns=columns)
  cacheFetch(x, var, {
    transcripts(x@.txdb, vals, columns)
  })
})


setMethod("exons", c(x="GenomicCache"),
function(x, vals=NULL, columns="exon_id") {
  var <- generateCacheName('exons', vals=NULL, columns=columns)
  cacheFetch(x, var, {
    exons(x@.txdb, vals, columns)
  })
})


setMethod("exonsBy", c(x="GenomicCache"),
function(x, by=c('tx', 'gene'), use.names=FALSE, ...) {
  by <- match.arg(by)
  var <- generateCacheName('exonsBy', by=by, use.names=use.names)
  cacheFetch(x, var, {
    exonsBy(x@.txdb, by, use.names=use.names)
  })
})


setMethod("cds", c(x="GenomicCache"),
function(x, vals=NULL, columns="cds_id") {
  var <- generateCacheName('cds', vals=vals, columns=columns)
  cacheFetch(x, var, {
    cds(x@.txdb, vals, columns)
  })
})

setMethod("cdsBy", c(x="GenomicCache"),
function(x, by=c('tx', 'gene'), use.names=FALSE, ...) {
  by <- match.arg(by)
  var <- generateCacheName('cdsBy', by=by, use.names=use.names)
  cacheFetch(x, var, {
    cdsBy(x@.txdb, by, use.names=use.names)
  })
})

setMethod("fiveUTRsByTranscript", c(x="GenomicCache"),
function(x, use.names=FALSE, ...) {
  var <- generateCacheName('fiveUTRsByTranscript', use.names=FALSE)
  cacheFetch(x, 'utr5', {
    fiveUTRsByTranscript(x@.txdb, use.names=use.names)    
  })
})

setMethod("threeUTRsByTranscript", c(x="GenomicCache"),
function(x, use.names=FALSE, ...) {
  var <- generateCacheName('threeUTRsByTranscript', use.names=FALSE)
  cacheFetch(x, 'utr3', {
    threeUTRsByTranscript(x@.txdb, use.names=use.names)
  })
})

setMethod("genome", c(x="GenomicCache"),
function(x, ...) {
  x@.genome
})

setMethod("dataSource", c(x="GenomicCache"),
function(x, ...) {
  dataSource(x@.txdb)
})

setMethod("annotationSource", c(x="GenomicCache"),
function(x, ...) {
  annotationSource(x@.txdb)
})

