## Used to represent contiguous "transcriptome space" while keeping track of
## potential gaps in genomic space (ie. as would occur from introns)
setClass("GappedRanges",
         representation(seqnames="Rle",
                        ranges="IRanges",
                        gaps="IRangesList",
                        strand="Rle",
                        metadata="list"),
         prototype(seqnames=Rle(factor()),
                   ranges=IRanges(),
                   gaps=IRangesList(),
                   strand=Rle(strand()),
                   elementMetadata=NULL),
         contains="Sequence")

GappedRanges <- function(seqnames='*', ranges=IRanges(), gaps=NULL,
                         strand=Rle("*", length(ranges)), values=NULL, ...) {
  if (!is(seqnames, 'Rle')) {
    seqnames <- Rle(seqnames)
  }
  if (!is.factor(runValue(seqnames))) {
    runValue(seqnames) <- factor(runValue(seqnames))
  }
  
  if (!is(ranges, 'IRanges')) {
    ranges <- as(ranges, 'IRanges')
  }
  
  if (!is(strand, "Rle")) {
    strand <- Rle(strand)
  }
  
  if (!is.factor(runValue(strand)) ||
      !identical(levels(runValue(strand)), levels(strand()))) {
    runValue(strand) <- strand(runValue(strand))
  }
  
  if (IRanges:::anyMissing(runValue(strand))) {
    warning("missing values in strand converted to \"*\"")
    runValue(strand)[is.na(runValue(strand))] <- "*"
  }

  lx <- max(length(seqnames), length(ranges), length(strand))
  if (lx > 1) {
    if (length(seqnames) == 1)
      seqnames <- rep(seqnames, lx)
    if (length(ranges) == 1)
      ranges <- rep(ranges, lx)
    if (length(strand) == 1)
      strand <- rep(strand, lx)
  }

  ## Set and check the gaps information for the ranges
  if (is.null(gaps)) {
    gaps <- do.call("IRangesList", replicate(length(ranges), IRanges()))
  }
  if (length(gaps) < lx) {
    xg <- do.call('IRangesList', replicate(lx - length(gaps), IRanges()))
    gaps <- c(gaps, xg)
  }
  
  if (!validGappedRangesGaps(ranges, gaps)) {
    stop("Invalid gaps object")
  }
  
  ## Create a default DataFrame to drop into the elementMetadata slot
  if (is.null(values)) {
    values <- DataFrame(...)
    if (ncol(values) == 0) {
      values <- new("DataFrame", nrows=length(ranges))
    }
    if (!is.null(rownames(values))) {
      if (!is.null(names(ranges))) {
        names(ranges) <- rownames(values)
      }
      rownames(values) <- NULL
    }
  }

  new('GappedRanges', seqnames=seqnames, ranges=ranges, gaps=gaps,
      strand=strand, elementMetadata=values)
}

validGappedRangesGaps <- function(ranges, gaps=NULL) {
  if (is(ranges, 'GappedRanges')) {
    gaps <- gaps(ranges)
    ranges <- ranges(ranges)
  }
  if (!is(ranges, 'IRanges')) {
    stop("Invalid ranges object")
  }
  if (!is(gaps, 'IRangesList')) {
    stop("Invalid gaps object")
  }
  if (length(ranges) != length(gaps)) {
    message("Differing lengths for ranges and gaps")
    return(FALSE)
  }

  ## Check that all gaps lie within their conjugate ranges
  for (i in 1:length(ranges)) {
    grange <- range(gaps[[i]])
    if (length(grange) > 0) {
      rrange <- ranges[i]
      if (start(rrange) > start(grange) || end(rrange) < end(grange)) {
        message("Out of bounds gaps for range ", i)
        return(FALSE)
      }
    }
  }

  TRUE
}

setMethod("show", c(object="GappedRanges"),
function(object) {
  cat("GappedRanges of length", length(object), "\n")
  df <- as.data.frame(object)
  df$gaps <- sapply(gaps(object), length)
  print(df, quote=FALSE, right=TRUE)
})

setGeneric("gappedWidth", function(x, ...) standardGeneric("gappedWidth"))
setMethod("gappedWidth", c(x="GappedRanges"),
function(x, ...) {
  .ranges <- ranges(x)
  .gaps <- gaps(x)
  width(.ranges) - sapply(width(.gaps), sum)
})


setMethod("length", c(x="GappedRanges"),
function(x) {
  length(ranges(x))
})

setMethod("names", c(x="GappedRanges"),
function(x) {
  names(ranges(x))
})

setMethod("ranges", c(x="GappedRanges"),
function(x, ...) {
  x@ranges
})

setMethod("strand", c(x="GappedRanges"),
function(x) {
  x@strand
})

setMethod("start", c(x="GappedRanges"),
function(x, ...) {
  start(ranges(x))
})

setMethod("end", c(x="GappedRanges"),
function(x, ...) {
  end(ranges(x))
})

setMethod("seqnames", c(x="GappedRanges"),
function(x) {
  x@seqnames
})

setMethod("gaps", c(x="GappedRanges"),
function(x, start=NA, end=NA) {
  x@gaps
})

setMethod("width", c(x="GappedRanges"),
function(x) {
  width(ranges(x))
})

setMethod("as.data.frame", c(x="GappedRanges"),
function(x, row.names=NULL, optional=FALSE) {
  ranges <- ranges(x)
  if (missing(row.names))
    row.names <- names(x)
  if (!is.null(names(x)))
    names(x) <- NULL
  data.frame(seqnames=as.factor(seqnames(x)),
             start=start(x),
             end=end(x),
             width=width(x),
             strand=as.factor(strand(x)),
             as.data.frame(elementMetadata(x)),
             row.names=row.names,
             stringsAsFactors=FALSE)
})
