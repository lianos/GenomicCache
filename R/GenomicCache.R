setMethod("initialize", "GenomicCache",
function(.Object, ...,
         .txdb=NULL,
         .cache=new.env()) {
  callNextMethod(.Object,
                 .txdb=.txdb,
                 .cache=.cache,
                 ...)
})

GenomicCache <- function(txdb, pre.load=c('transcripts', 'exons', 'utr3')) {
  can.load <- c('transcripts', 'exons', 'utr3', 'utr5')
  bad.load <- which(!pre.load %in% can.load)
  if (length(bad.load) > 0) {
    warning("Illegal things to load: ", paste(pre.load[bad.load], collapse=","))
    pre.load <- pre.load[-bad.load]
  }
  
  if (is.character(txdb)) {
    if (!file.exists(txdb)) {
      stop("Can't read TranscriptDb file: ", txdb)
    }
    txdb <- loadFeatures(txdb)
  }
  if (!inherits(txdb, 'TranscriptDb')) {
    stop("Need a valid TranscriptDb object to continue")
  }

  for (what in pre.load) {
    
  }
  
  new('GenomicCache',
      .txdb=txdb)
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
  cacheFetch(x, 'transcripts', {
    transcripts(x@.txdb)
  })
})


setMethod("exons", c(x="GenomicCache"),
function(x, vals) {
  cacheFetch(x, 'exons', {
    exons(x@.txdb, vals)
  })
})


setMethod("exonsBy", c(x="GenomicCache"),
function(x, by) {
  by <- match.arg(by)
  var.name <- sprintf("exonsBy.%s", by)
  cacheFetch(x, var.name, {
    exonsBy(x@.txdb, by)
  })
})


setMethod("cds", c(x="GenomicCache"),
function(x, ...) {
  cacheFetch(x, 'cds', {
    cds(x@.txdb)
  })
})

setMethod("cdsBy", c(x="GenomicCache"),
function(x, by) {
  by <- match.arg(by)
  var.name <- sprintf('cdsBy.%s', by)
  cacheFetch(x, var.name, {
    cdsBy(x@.txdb, by)
  })
})

setMethod("fiveUTRsByTranscript", c(x="GenomicCache"),
function(x) {
  cacheFetch(x, 'utr5', {
    fiveUTRsByTranscript(x@.txdb)    
  })
})

setMethod("threeUTRsByTranscript", c(x="GenomicCache"),
function(x) {
  cacheFetch(x, 'utr3', {
    threeUTRsByTranscript(x@.txdb)
  })
})

setMethod("genome", c(x="GenomicCache"),
function(x, ...) {
  genome(x@.txdb)
})

setMethod("dataSource", c(x="GenomicCache"),
function(x, ...) {
  dataSource(x@.txdb)
})

setMethod("annotationSource", c(x="GenomicCache"),
function(x, ...) {
  annotationSource(x@.txdb)
})

