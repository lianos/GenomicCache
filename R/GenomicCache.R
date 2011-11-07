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

##' Create a GenomicCache object and save it to filesystem.
##'
##' This function could be driven from the command line via the
##' \code{gfx-create-cache} script.
##'
##' @export
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param genome The assembly of the genome to use (hg18, mm9, etc.)
##' @param annotation Which annotation source to use (refseq, ensembl, aceview)
##' @param path The directory to create the \code{GenomicCache} directory in.
##' @param gc.name The name of the GenomicCache directory to create. Leave this
##' \code{NULL} for a "reasonable" one to be generated for you.
##' @param is.table.name Is \code{annotation} a UCSC table name (refGene, etc.)?
##' @param table2name A named character vector the user can pass to inform the
##' UCSC table name to "normal" name mapping. The \code{names()} are the
##' UCSC table names, the values are the "normal" names. This can be left
##' \code{NULL} and the default (refseq, ucsc, aceview) one will be used.
##' @return A \code{GenomicCache} object (invisbly)
createGenomicCache <- function(genome, annotation, path='.', gc.name=NULL,
                               is.table.name=FALSE, table2name=NULL) {
  table2name <- c(refGene='refseq', knownGene='ucsc', ensGene="ensembl",
                  acembly='aceview')
  if (is.table.name) {
    ucsc.table <- match.arg(annotation, names(table2name))
    annotation <- table2name[ucsc.table]
  } else {
    annotation <- match.arg(annotation, table2name)
    ucsc.table <- names(table2name)[table2name == annotation]
  }

  if (is.null(gc.name)) {
    gc.name <- paste("GenomicCache", genome, annotation, sep=".")
  }
  gc.path <- file.path(path, gc.name)

  ## Setup directory structure
  fi <- file.info(gc.path)
  if (!is.na(fi$isdir)) {
    stop("GenomicCache directory already exists, ", gc.path)
  }

  dirs <- c(gc.path, paste(gc.path, c('cache', 'features'), sep="/"))
  for (dir in dirs) {
    if (!dir.create(dir)) {
      stop("Could not create directory:", dir, "\nCheck permissions?")
    }
  }

  txdb <- makeTranscriptDbFromUCSC(genome=genome, tablename=ucsc.table)
  fn <- paste('TranscriptDb', genome, annotation, 'sqlite', sep=".")
  saveFeatures(txdb, file.path(gc.path, 'features', fn))

  invisible(GenomicCache(gc.path, pre.load=NULL))
}

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

  genome <- subset(metadata(txdb), name == "Genome")$value

  gc <- new('GenomicCache',
            .genome=genome,
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

## Get the gene symbol for these genes
id2symbol <- function(x, ids=NULL) {
  stopifnot(inherits(x, 'GenomicCache') || inherits(x, 'TranscriptDb'))
  .genome <- genome(x)

  if (is.null(ids)) {
    SQL <- "SELECT * FROM gene"
    ## ids <- GenomicFeatures:::dbEasyQuery(GenomicFeatures:::txdbConn(x), SQL)
    ## Changed in bioc 2.9:
    ids <- GenomicFeatures:::dbEasyQuery(txdbConn(x), SQL)
    ids <- ids$gene_id
  }

  symbols <- getSymbolFromEntrezId(.genome, ids, rm.unknown=FALSE)
  data.frame(entrez=ids, symbol=sapply(symbols, '[', 1))
}

setMethod("show", c(object="GenomicCache"),
function(object) {
  cat("GenomicCache object:\n")
  cat("  Genome:", genome(object), "\n")
  cat("  Annotation source:", annotationSource(object), "\n")
})

## Makes a copy of the GenomicCache ensuring that it has a separate connection
## to the transcript database. \code{pre.load} is set to \code{NULL} because we
## often don't want to waste time loading things since this is likely called
## when running within
setMethod("duplicate", c(x="GenomicCache"),
function(x, pre.load=NULL, ...) {
  GenomicCache(x@.path, pre.load=pre.load)
})

setMethod("dispose", c(x="GenomicCache"),
function(x, ...) {
  clearCache(x)
  # sqliteCloseConnection(GenomicFeatures:::txdbConn(txdb(x)))
  sqliteCloseConnection(GenomicFeatures:::txdbConn(txdb(x)))
})

setMethod("seqnames", c(x="GenomicCache"),
function(x) {
  seqlevels(txdb(x))
})

setMethod("seqlengths", c(x="GenomicCache"), function(x) {
  seqlengths(txdb(x))
})

setMethod('seqlevels', c(x="GenomicCache"),
function(x) {
  seqlevels(txdb(x))
})

setMethod("seqinfo", c(x="GenomicCache"),
function(x) {
  seqinfo(txdb(x))
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

setMethod("transcriptsBy", c(x="GenomicCache"),
function(x, by=c("gene", "exon", "cds"), use.names=FALSE) {
  by <- match.arg(by)
  var <- generateCacheName('transcriptsBy', by=by, use.names=use.names)

  cacheFetch(x, var, {
    xc <- transcriptsBy(x@.txdb, by, use.names)
    if (by == 'gene') {
      ## Associate gene symbols to the id's
      symbols <- id2symbol(x, names(xc))
      values(xc) <- DataFrame(symbol=as.character(symbols$symbol))
    }
    xc
  })
})


setMethod("exons", c(x="GenomicCache"),
function(x, vals=NULL, columns="exon_id") {
  var <- generateCacheName('exons', vals=NULL, columns=columns)
  cacheFetch(x, var, exons(x@.txdb, vals, columns))
})


setMethod("exonsBy", c(x="GenomicCache"),
function(x, by=c('tx', 'gene'), use.names=FALSE, ...) {
  by <- match.arg(by)
  var <- generateCacheName('exonsBy', by=by, use.names=use.names)
  cacheFetch(x, var, exonsBy(x@.txdb, by, use.names=use.names))
})


setMethod("cds", c(x="GenomicCache"),
function(x, vals=NULL, columns="cds_id") {
  var <- generateCacheName('cds', vals=vals, columns=columns)
  cacheFetch(x, var, cds(x@.txdb, vals, columns))
})

setMethod("cdsBy", c(x="GenomicCache"),
function(x, by=c('tx', 'gene'), use.names=FALSE, ...) {
  by <- match.arg(by)
  var <- generateCacheName('cdsBy', by=by, use.names=use.names)
  cacheFetch(x, var, cdsBy(x@.txdb, by, use.names=use.names))
})

setMethod("fiveUTRsBy", c(x="GenomicCache"),
function(x, by, use.names, flank.up=1000, flank.down=1000, ...) {
  by <- match.arg(by)
  if (by == 'tx') {
    utrs <- fiveUTRsByTranscript(x, use.names, ...)
  } else {

  }
  if (flank.up != 0) {
    utrs <- endoapply(utrs, function(x) {
      resize(x, width=width(x) + flank.up, fix='start')
    })
  }
  if (flank.down != 0) {
    utrs <- endoapply(utrs, function(x) {
      resize(x, width=width(x) + flank.down, fix='end')
    })
  }

  utrs
})

setMethod("fiveUTRsByTranscript", c(x="GenomicCache"),
function(x, use.names=FALSE, ...) {
  var <- generateCacheName('fiveUTRsByTranscript', use.names=FALSE)
  cacheFetch(x, var, fiveUTRsByTranscript(x@.txdb, use.names=use.names))
})


setMethod("threeUTRsBy", c(x="GenomicCache"),
function(x, by, use.names, flank.up=100, flank.down=1000, ...) {
  by <- match.arg(by)
  if (by == 'tx') {
    utrs <- threeUTRsByTranscript(x, use.names, ...)
  } else {

  }

  ## These endoapply's take way too long
  ## if (flank.up != 0) {
  ##   utrs <- endoapply(utrs, function(x) {
  ##     resize(x, width=width(x) + flank.up, fix='start')
  ##   })
  ## }
  ## if (flank.down != 0) {
  ##   utrs <- endoapply(utrs, function(x) {
  ##     resize(x, width=width(x) + flank.down, fix='end')
  ##   })
  ## }

  utrs
})

setMethod("threeUTRsByTranscript", c(x="GenomicCache"),
function(x, use.names=FALSE, ...) {
  var <- generateCacheName('threeUTRsByTranscript', use.names=FALSE)
  cacheFetch(x, var, threeUTRsByTranscript(x@.txdb, use.names=use.names))
})

setMethod("genome", c(x="GenomicCache"),
function(x, ...) {
  x@.genome
})

setMethod("genome", c(x="TranscriptDb"),
function(x, ...) {
  subset(metadata(x), name == "Genome")$value
})

##' @importFrom annotate dataSource
setMethod("dataSource", c(object="GenomicCache"),
function(object) {
  subset(metadata(object@.txdb), name == "Data source")$value
})

setMethod("annotationSource", c(object="GenomicCache"),
function(object) {
  subset(metadata(object@.txdb), name == "UCSC Table")$value
})

setMethod("annotationSource", c(object="TranscriptDb"),
function(object) {
  subset(metadata(object), name == "UCSC Table")$value
})

setMethod("cacheDir", c(x="GenomicCache"),
function(x, ..., global=FALSE) {
  if (global) {
    cdir <- Sys.getenv(.GFX$cache$environment.key)
  } else {
    cdir <- file.path(x@.path, 'cache')
  }
  more <- list(...)
  if (length(more) > 0L) {
    cdir <- do.call(file.path, c(list(cdir), more))
  }
  cdir
})

setMethod("txdb", c(x="GenomicCache"),
function(x, ...) {
  x@.txdb
})
