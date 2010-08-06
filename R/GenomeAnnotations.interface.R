## Functions that mimic my GenomeAnnotation calls
## (for package migration purposes)
.geneCacheFileName <- function(gcache, chromosome) {
  paste(genome(gcache), annotationSource(gcache), "genes", chromosome, "rda", sep=".")
}

loadGFXGeneModels <- function(gcache, chromosome, cache.dir=NULL) {
  cache.dir <- .getCacheDir(cache.dir)
  file.name <- .geneCacheFileName(gcache, chromosome)
  fpath <- file.path(cache.dir, file.name)
  if (is.na(file.info(fpath)$isdir)) {
    return(NULL)
  }
  obj <- load(fpath)
  get(obj)
}

generateGFXGeneModels <- function(gcache, chromosomes=NULL, cache.dir=NULL) {
  require(plyr)
  cache.dir <- .getCacheDir(cache.dir)
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

setGeneric("getGenesOnChromosome",
function(x, chromosome, start=NULL, end=NULL, cache.dir=NULL) {
  standardGeneric("getGenesOnChromosome")
})
setMethod("getGenesOnChromosome", c(x="GenomicCache"),
function(x, chromosome, start, end, cache.dir=NULL) {
  genes <- loadGFXGeneModels(x, chromosome, cache.dir=cache.dir)
  if (is.null(genes)) {
    stop("Build the gene models first!")
  }
  genes
})
