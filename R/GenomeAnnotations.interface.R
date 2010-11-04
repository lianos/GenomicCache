## Functions that mimic my GenomeAnnotation calls
## (for package migration purposes)
.geneCacheFileName <- function(gcache, chromosome) {
  paste(genome(gcache), annotationSource(gcache), "genes", chromosome, "rda",
        sep=".")
}

loadGFXGeneModels <- function(gcache, chromosome, cache.dir=NULL) {
  cache.dir <- cacheDir(gcache, 'gene.models')
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
generateGFXGeneModels <- function(gcache, chromosomes=NULL,
                                  flank.up=c(0, 500, 1000),
                                  flank.down=c(0, 500, 1000)) {
  if (!require(plyr)) stop("Plyr is used here")
  cache.dir <- cacheDir(gcache, 'gene.models')
  if (!dir.exists(cache.dir)) {
    cat("Making cache directory:", cache.dir, "\n")
    if (!dir.create(cache.dir)) {
      stop("Error creating directory")
    }
  }
  
  if (is.null(chromosomes)) {
    chromosomes <- seqnames(gcache)
  }

  xcripts <- transcripts(gcache)

  foreach(chr=chromosomes, .packages="GenomicFeaturesX", .inorder=FALSE) %dopar% {
    cat("===", chr, "===...\n")
    .gc <- duplicate(gcache)
    on.exit(dispose(.gc))

    chr.xcripts <- xcripts[which(seqnames(xcripts) == chr)]
    .so.far <- character(length(chr.xcripts))
    .idx <- 1L
    
    genes <- llply(values(chr.xcripts)$tx_name, function(tx) {
      if (tx %in% .so.far) {
        return(NULL)
      }
      g <- tryCatch(GFGene(tx.id=tx, .gc), error=function(e) NULL)
      if (!is.null(g)) {
        ## cat(chr, symbol(g), "\n")
        xcripts <- transcripts(g, which.chr=chr)
        to <- .idx + length(xcripts) - 1
        .so.far[.idx:to] <<- values(xcripts)$tx_name
        .idx <<- to + 1
        ## Generate idealized models
        for (i in 1:length(flank.up)) {
          idealized(g, which.chr=chr, flank.up=flank.up[i],
                    flank.down=flank.down[i])
        }
      }
      g
    }, .progress='none')

    genes <- genes[!sapply(genes, is.null)]
    names(genes) <- make.unique(sapply(genes, symbol))
    save(genes, file=file.path(cache.dir, .geneCacheFileName(.gc, chr)))
    chr
  }
  
}

## NOTE: genesOnChromosome is not done
setMethod("getGenesOnChromosome", c(x="GenomicCache"),
function(x, chromosome, start, end, maxgap, minoverlap, overlap.type,
         use.cache, cache.dir=cacheDir(x), ...) {
  xcripts <- transcripts(x)
  if (use.cache) {
    genes <- loadGFXGeneModels(x, chromosome, cache.dir=cache.dir)
    if (is.null(genes)) {
      stop("Build the gene models first!")
    }
    
    if (!is.null(start) || !is.null(end)) {
      ## NOTE: This is poorly implemented
      if (is.null(start)) {
        start <- 1L
      }
      if (is.null(end)) {
        end <- xcripts[[chromosome]]
      }
      tx.bounds <- txBounds(x, chromosome, .genes=genes)
      subset.bounds <- IRangest(start, end)
      o <- findOverlaps(tx.bounds, subset.bounds, maxgap=maxgap,
                        minoverlap=minoverlap, type=overlap.type)
      genes <- genes[queryHits(o)]
    }
    
    return(genes)
  } else {
    stop("getGenesOnChromosome w/o cache not implemented yet.")

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
  }
})

setMethod("txBounds", c(x="GenomicCache"),
function(x, chromosome, .genes=NULL, ...) {
  var.name <- paste('txBounds', chromosome, sep=".")
  cacheFetch(x, var.name, {
    if (is.null(.genes)) {
      .genes <- getGenesOnChromosome(x, chromosome)
    }
    txBoundsFromGeneList(.genes)
  })
})

txBoundsFromGeneList <- function(genes) {
  tx.bounds <- RangesList(lapply(genes, function(g) {
    range(unlist(ranges(transcripts(g))))
  }))
  unlist(tx.bounds)
}
