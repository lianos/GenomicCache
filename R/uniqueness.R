## TODO: Debug this -- it was written cold, ie. not tested/run at all.
## 
## I'm planning to get unique reads from uniquely-aligned BAMs, but
## this could still be useful.
withinUniquenessBounds <- function(k) {
  .GFX$uniqueness$min.length <= k && k <= .GFX$uniqueness$max.length
}

## you can pass in a uniqueness.map to use one that's already been loaded.
## if this is done and use.universal is TRUE, we assume it is the universall/all
## map, otherwise we assume that the map is for the same width that the ranges
## and all of the ranges must be the same width
setMethod("flagUniqueRanges", c(ranges="IRanges"),
function(ranges, genome, chromosome, use.universal=FALSE,
         uniqueness.map=NULL, ...) {
  if (missing(genome)) {
    stop("Need genome")
  }
  if (is.null(chromosome)) {
    stop("Need chromosome")
  }
  if (is.null(use.universal)) {
    use.universal <- FALSE
  }
  
  are.unique <- logical(length(ranges))  
  widths <- unique(width(ranges))

  if (!is.null(uniqueness.map)) {
    if (!use.universal && length(width) > 1) {
      warning("We can't use the map passed in, will load it on the fly")
      uniqueness.map <- NULL
    }
  }

  if (use.universal && is.null(uniqueness.map)) {
    uniqueness.map <- getUniquenessMap(genome, chromosome, 'all')
  }
  
  for (w in widths) {
    these <- width(ranges) == w
    starts <- start(ranges[these])
    if (use.universal) {
      if (withinUniquenessBounds(w)) {
        is.unique <- as.logical(uniqueness.map[starts] <= w)
      } else {
        warning("Width of range (", w, ") is not within bounds, ",
                "all ranges marked as FALSE")
        is.unique <- logical(length(starts))
      }
    } else {
      uniqueness.map <- tryCatch(getUniquenessMap(genome, chromosome, w),
                                 error=function(e) NULL)
      if (is.null(uniqueness.map)) {
        warning("No unique map for ",
                sprintf("%s.%s %d-mer", genome, chromosome, w),
                " all ranges of this length are FALSE")
      }
      is.unique <- as.logical(uniqueness.map[starts] == 1L)
    }
    are.unique[these] <- is.unique
  }
  
  are.unique
})

setMethod("flagUniqueRanges", c(ranges="GRanges"),
function(ranges, genome, chromosome, use.universal=FALSE, ...) {
  chromosomes <- unique(as.character(seqnames(ranges)))
  result <- logical(length(ranges))
  
  for (chr in chromosomes) {
    these <- which(seqnames(ranges) == chr)
    u <- flagUniqueRanges(ranges(ranges[these]), genome, chr, use.universal, ...)
    result[these] <- u
  }

  result
})

################################################################################
## Utility functions

##' @param genome hg18, mm9, etc
##' @param chromosome chr1, chrX, etc.
##' @param k the length of the kmer to check, or 'all' to use the "universal"
##' uniqueness map
.uniquenessMapFilePath <- function(genome, chromosome, k, cache.dir=NULL) {
  if (suppressWarnings(is.na(as.integer(k)))) {
    if (!is.character(k) && (k != "all")) {
      stop("k can only be an integer, or 'all'")
    }
  } else {
    if (!withinUniquenessBounds(k)) {
      stop("Uniqueness maps only exists for kmers >= ",
           .GFX$uniqueness$min.length, " and <= ",
           .GFX$uniqueness$max.length)
    }
    k <- sprintf("k%d", k)
  }
  
  cache.dir <- GFXCacheDir(cache.dir, 'uniqueness', genome, k)
  name <- paste(chromosome, "Rle.rda", sep=".")
  file.path(cache.dir, name)
}


## x : genome
setMethod("getUniquenessMap", c(x="character"),
function(x, chromosome, k, ...) {
  file.name <- .uniquenessMapFilePath(x, chromosome, k)
  if (!file.exists(file.name)) {
    stop("No uniqueness map available for genome ", x,
         " and read length ", k)
  }
  load(file.name)
  uniqueness
})

setMethod("getUniquenessMap", c(x="GenomicCache"),
function(x, chromosome, k, ...) {
  getUniquenessMap(genome(x), chromosome, k, allow.universal, ...)
})

###############################################################################
# Functions to convert Anshul's uniqueness maps to Rle objects
# These functions aren't meant to be used by the general public
###############################################################################
convertUniqueness2Rle <- function(u.dir, genome='hg18') {
  files <- dir(u.dir, pattern="uint8.unique\\.*")
  gdb <- GenomeDB('hg18')
  chromosomes <- sapply(strsplit(files, ".", fixed=TRUE), '[', 1)
  
  uniqueness <- lapply(seq_along(files), function(idx) {
    file.name <- files[idx]
    chr <- chromosomes[idx]
    chr.length <- getChromosomeLength(gdb, chr)
    cat("Chromosome", chr, "-- length:", chr.length, "\n")
    # cat("  ", file.path(u.dir, file.name), "\n")
    u <- readBin(file.path(u.dir, file.name), what='integer', size=1,
                 signed=FALSE, n=chr.length+100)
    if (length(u) != chr.length) {
      cat("  Lengths don't match")
    }
    uniqueness <- Rle(u)
    save(uniqueness, file=file.path(u.dir, paste(chr, "Rle.rda", sep=".")))
    Rle(u)
  })
  names(uniqueness) <- chromosomes
  
  invisible(uniqueness)
}
