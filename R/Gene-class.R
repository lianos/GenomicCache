setClass("GFGene", representation(id='character',
                                  symbol='character',
                                  chromosome='character',
                                  strand='factor',
                                  transcripts="GRangesList",
                                  start='integer',
                                  end='integer',
                                  cache='environment'))

# setClass("RefSeqGene", contains="GFGene")
# setClass("EnsemblGene", contains="GFGene")
# setClass("AceviewGene", contains="GFGene")

setMethod("initialize", "GFGene",
function(.Object, ...,
         id=character(),
         symbol=character(),
         chromosome=character(),
         strand=strand(),
         transcripts=GRanges(),
         start=integer(),
         end=integer(),
         cache=new.env()) {
  callNextMethod(.Object,
                 id=id,
                 symbol=symbol,
                 chromosome=chromosome,
                 strand=strand,
                 transcripts=transcripts,
                 start=start,
                 end=end,
                 cache=cache,
                 ...)
})

GFGene <- function(..., .txdb=NULL, .exons=NULL) {
  args <- list(...)
  arg.names <- names(args)

  if (is.null(.txdb)) {
    .txdb <- takeFromListByType(args, 'TranscriptDb')
    if (is.null(.txdb)) {
      stop("TranscriptDb (.txdb) expected")
    }
  }
  
  a.source <- annotationSource(.txdb)
  class.name <- switch(a.source, refGene="RefSeqGene",
                       ensGene="EnsemblGene", aceGene="AceviewGene",
                       stop("Unknown source: ", a.source))
  
  id <- takeFromListByType(args, 'character', index=TRUE)
  if (is.null(id)) {
    stop("No gene name/id passed")
  }
  
  id.type <- names(args)[id]
  id <- args[[id]]

  if (is.null(id.type)) {
    id.type <- 'symbol'
  }
  
  if (id.type == 'symbol') {
    symbol <- id
    entrez.id <- getEntrezIdFromSymbol(.txdb, symbol)
  } else if (id.type == 'id') {
    ## Will throw an error if its not an ensembl source
    entrez.id <- getEntrezIdFromGeneId(.txdb, id)
    symbol <- getSymbolFromEntrezId(.txdb, entrez.id)
  } else if (id.type == 'tx.id') {
    entrez.id <- getTranscriptIdFromEntrezId(.txdb, id)
    symbol <- getSymbolFromEntrezId(.txdb, entrez.id)
  }
  
  if (a.source == 'refGene') {
    id <- symbol
  } else if (a.source == 'ensGene') {
    id <- getGeneIdFromEntrezId(.txdb, entrez.id)
  } else {
    stop("Unknown source: ", a.source)
  }
  
  tx.id <- getTranscriptIdFromEntrezId(.txdb, entrez.id)
  
  if (is.null(.exons)) {
    .exons <- takeFromListByType(args, 'GRangesList')
  }
  if (is.null(.exons)) {
    .exons <- exonsBy(.txdb, by='tx')
  }
  
  xcripts <- .exons[tx.id]
  #new(class.name,
  new("GFGene",
      id=id,
      symbol=symbol,
      strand=as.vector(strand(xcripts[[1]])[1]),
      chromosome=as.character(seqnames(xcripts[[1]])[1]),
      transcripts=xcripts
      )
}
