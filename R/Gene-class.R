
setMethod("initialize", "GFGene",
function(.Object, ...,
         .id=character(),
         .symbol=character(),
         .chromosome=factor(),
         .transcripts=GRangesList(),
         .transcript.names=character(),
         .start=integer(),
         .end=integer(),
         .strand=strand(),
         .cache=new.env()) {
  callNextMethod(.Object,
                 .id=.id,
                 .symbol=.symbol,
                 .chromosome=.chromosome,
                 .transcripts=.transcripts,
                 .transcript.names=.transcript.names,
                 .start=.start,
                 .end=.end,
                 .cache=.cache,
                 .strand=.strand,
                 ...)
})

# GFGene <- function(..., .transcripts=NULL, .exons=NULL, .txdb=NULL) {
#   args <- list(...)
#   arg.names <- names(args)
# 
#   if (is.null(.txdb)) {
#     .txdb <- takeFromListByType(args, 'TranscriptDb')
#     if (is.null(.txdb)) {
#       stop("TranscriptDb (.txdb) expected")
#     }
#   }
#   
#   a.source <- annotationSource(.txdb)
#   class.name <- switch(a.source, refGene="RefSeqGene",
#                        ensGene="EnsemblGene", aceGene="AceviewGene",
#                        stop("Unknown source: ", a.source))
#   
#   id <- takeFromListByType(args, 'character', index=TRUE)
#   if (is.null(id)) {
#     stop("No gene name/id passed")
#   }
#   
#   id.type <- names(args)[id]
#   id <- args[[id]]
# 
#   if (is.null(id.type)) {
#     id.type <- 'symbol'
#   }
#   
#   if (id.type == 'symbol') {
#     symbol <- id
#     entrez.id <- getEntrezIdFromSymbol(.txdb, symbol)
#   } else if (id.type == 'id') {
#     ## Will throw an error if its not an ensembl source
#     entrez.id <- getEntrezIdFromGeneId(.txdb, id)
#     symbol <- getSymbolFromEntrezId(.txdb, entrez.id)
#   } else if (id.type == 'tx.id') {
#     entrez.id <- getTranscriptIdFromEntrezId(.txdb, id)
#     symbol <- getSymbolFromEntrezId(.txdb, entrez.id)
#   }
#   
#   if (a.source == 'refGene') {
#     id <- character()
#   } else if (a.source == 'ensGene') {
#     id <- getGeneIdFromEntrezId(.txdb, entrez.id)
#   } else {
#     stop("Unknown source: ", a.source)
#   }
#   
#   if (is.null(.exons)) {
#     .exons <- takeFromListByType(args, 'GRangesList')
#   }
#   if (is.null(.exons)) {
#     .exons <- exonsBy(.txdb, by='tx')
#   }
#   
#   if (is.null(.transcripts)) {
#     .transcripts <- takeFromListByType(args, 'GRanges')
#   }
#   if (is.null(.transcripts)) {
#     .transcripts <- transcripts(.txdb)
#   }
#   
#   tx.name <- getTranscriptIdFromEntrezId(.txdb, entrez.id)
#   xcripts <- subset(.transcripts, values(.transcripts)$tx_name %in% tx.name)
#   xm <- values(xcripts)
#   
#   #new(class.name,
#   new("GFGene",
#       id=id,
#       symbol=symbol,
#       strand=as.vector(strand(.transcripts)[1]),
#       chromosome=as.vector(seqnames(.transcripts)[1]),
#       transcripts=.exons[as.character(xm$tx_id)],
#       transcript.names=xm$tx_name
#       )
# }
