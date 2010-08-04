setClass("GFGene", representation(id='character',
                                  symbol='character',
                                  transcripts="GRanges",
                                  exons="GRangesList",
                                  chromosome='character',
                                  strand='character',
                                  start='integer',
                                  end='integer',
                                  cache='environment'))
setClass("RefSeqGene", contains="GFGene")
setClass("EnsemblGene", contains="GFGene")
setClass("AceviewGene", contains="GFGene")

setMethod("initialize", "GFGene",
function(.Object, ...,
         id=character(),
         symbol=character(),
         transcripts=GRanges(),
         exons=GRangesList(),
         chromosome=character(),
         strand=character(),
         start=integer(),
         end=integer(),
         cache=new.env()) {
  callNextMethod(.Object,
                 id=id,
                 symbol=symbol,
                 transcripts=transcripts,
                 exons=exons,
                 chromosome=chromosome,
                 strand=strand,
                 start=start,
                 end=end,
                 cache=cache,
                 ...)
})

##' Returns an object of type \code{type} from a list
##' If this object isn't found, or other error, returns \code{NULL}
.getGFGeneConstructorObject <- function(args, type, multi=FALSE, index=FALSE) {
  take <- which(sapply(args, function(arg) inherits(arg, type)))
  if (length(take) == 0L) {
    return(NULL)
  }
  
  if (length(take) > 1) {
    if (is.logical(multi[1]) && !multi) {
      warning("Multiple objects of type ", type, " found.")
      take <- '..NOTHING..'
    } else if (is.integer(multi)) {
      if (any(multi > length(take)) || any(multi < 0L)) {
        warning("multi take subscript(s) out of bounds")
        take <- '..NOTHING..'
      }
    } else {
      warning("Illegal type of multi argument: ", is(multi)[1])
      take <- '..NOTHING..'
    }
  }

  if (index) {
    ret.val <- take
  } else {
    ret.val <- if (length(take) > 1) args[take] else args[[take]]
  }

  ret.val
}

GFGene <- function(..., .transcripts=NULL, .exons=NULL, .txdb=NULL) {
  args <- list(...)
  arg.names <- names(args)

  if (is.null(.txdb)) {
    .txdb <- .getGFGeneConstructorObject(args, 'TranscriptDb')
    if (is.null(.txdb)) {
      stop("TranscriptDb (.txdb) expected")
    }
  }
  
  id <- .getGFGeneConstructorObject(args, 'character', index=TRUE)
  if (is.null(id)) {
    stop("No gene name/id passed")
  }
  
  id.type <- names(args)[id]
  id <- args[[id]]

  if (is.null(id.type)) {
    id.type <- 'symbol'
  }

  id <- switch(id.type,
               symbol=getGeneIdFromSymbol(.txdb, id),
               id=id,
               tx.id=getGeneIdFromTranscriptId(.txdb, id),
               stop("Unknown id type: ", id.type))
  
  ## Try to fetch cached objects from args
  if (is.null(.transcripts)) {
    .transcripts <- .getGFGeneConstructorObject(args, 'GRangesList')
  }
  if (is.null(.exons)) {
    .exons <- .getGFGeneConstructorObject(args, 'GRanges')
  }
  
  if (is.null(.transcripts)) {
    .transcripts <- transcriptsBy(.txdb, 'gene')
  }
  
  if (is.null(.exons)) {
    .exons <- transcripts(.txdb)
  }
  
  ## Get down to business
  browser()
  xcript.ranges <- .transcripts[[id]]
  if (length(xcript.ranges) == 0) {
    stop("Gene not found in txdb: ", id)
  }
  
}
