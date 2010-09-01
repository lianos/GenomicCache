## Functions that mimic my GenomeAnnotation calls
## (for package migration purposes)
.geneCacheFileName <- function(gcache, chromosome) {
  paste(genome(gcache), annotationSource(gcache), "genes", chromosome, "rda",
        sep=".")
}

loadGFXGeneModels <- function(gcache, chromosome, cache.dir=NULL) {
  cache.dir <- GFXCacheDir(cache.dir, 'gene.models')
  file.name <- .geneCacheFileName(gcache, chromosome)
  fpath <- file.path(cache.dir, file.name)
  if (is.na(file.info(fpath)$isdir)) {
    return(NULL)
  }
  obj <- load(fpath)
  get(obj)
}

## ~ 6 Hours for RefSeq hg18
## ~ 8.3 hours for Ensembl hg18
generateGFXGeneModels <- function(gcache, chromosomes=NULL, cache.dir=NULL) {
  if (!require(plyr)) stop("Plyr is used here")
  cache.dir <- GFXCacheDir(cache.dir, 'gene.models')
  
  if (is.null(chromosomes)) {
    chromosomes <- seqnames(gcache)
  }

  xcripts <- transcripts(gcache)
  
  for (chr in chromosomes) {
    cat(chr, "...\n")
    chr.xcripts <- xcripts[which(seqnames(xcripts) == chr)]
    .so.far <- character(length(chr.xcripts))
    .idx <- 1L
    genes <- llply(values(chr.xcripts)$tx_name, function(tx) {
      if (tx %in% .so.far) {
        return(NULL)
      }
      g <- tryCatch(GFGene(tx.id=tx, gcache), error=function(e) NULL)
      if (!is.null(g)) {
        to <- .idx + length(transcripts(g)) - 1
        .so.far[.idx:to] <<- g@.transcript.names
        .idx <<- to + 1
      }
      g
    }, .progress='text')

    genes <- genes[!sapply(genes, is.null)]
    names(genes) <- uniquefy(sapply(genes, symbol))
    save(genes, file=file.path(cache.dir, .geneCacheFileName(gcache, chr)))
    chr
  }
  
}

setGeneric("genesOnChromosome",
function(x, chromosome, start=NULL, end=NULL, maxgap=0L, minoverlap=1L,
         overlap.type=c('any', 'start', 'end', 'within', 'equal'),
         use.cache=FALSE, cache.dir=NULL, ...) {
  standardGeneric("getGenesOnChromosome")
})


## NOTE: genesOnChromosome is not done
setMethod("genesOnChromosome", c(x="GenomicCache"),
function(x, chromosome, start, end, maxgap, minoverlap, overlap.type,
         use.cache, cache.dir, ...) {
  xcripts <- transcripts(x)
  if (is.null(end)) {
    end <- seqlengths(xcripts)[[chromosome]]
  }
  
  bounds <- switch(is(start)[1],
    "NULL"=GRanges(seqname=chromosome, ranges=IRanges(start=1, end=end)),
    numeric=GRanges(seqname=chromosome, ranges=IRanges(start=start, end=end)),
    integer=GRanges(seqname=chromosome, ranges=IRanges(start=start, end=end)),
    IRanges=GRanges(seqname=chromosome, ranges=start),
    stop("Illegal value for `start`"))
  
  if (use.cache) {
    genes <- loadGFXGeneModels(x, chromosome, cache.dir=cache.dir)
    if (is.null(genes)) {
      stop("Build the gene models first!")
    }
  } else {
    possibly <- subsetByOverlaps(xcripts, bounds, maxgap, minoverlap,
                                 overlap.type)
    if (length(possibly) > 0) {
      entrez <- getEntrezIdFromTranscriptId(x, values(possibly)$tx_name)
      entrez <- unique(unlist(entrez))
      symbols <- unlist(getSymbolFromEntrezId(x, entrez))
      genes <- lapply(symbols, GFGene, .gc=x)
      names(genes) <- sapply(genes, symbol)
    }
  }
  

})
