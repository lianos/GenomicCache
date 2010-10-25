GappedRanges <- function(irl=IRangesList(), ...) {
  if (is(irl, 'IRanges')) {
    irl <- IRangesList(irl)
  }
  args <- list(...)
  if (length(args) > 0) {
    more <- lapply(args, function(arg) {
      if (is(arg, 'IRanges')) arg else NULL
    })
    more <- more[!sapply(more, is.null)]
    if (length(more) > 0) {
      irl <- c(irl, IRangesList(more))
    }
  }
  
  as(irl, 'GappedRanges')
}

setReplaceMethod("[", "GappedRanges",
function(x, i, j, ..., value) {
  ## if (!missing(i)) {
  ##   iInfo <- IRanges:::.bracket.Index(i, length(x), names(x))
  ##   if (!is.null(iInfo[['msg']])) {
  ##     stop(iInfo[['msg']])
  ##   }
  ## }
  ## if (missing(i) || !iInfo[['useIdx']]) {
  ## 
  ## }
  if (is(value, 'GappedRanges')) {
    value <- ranges(value, mind.the.gap=TRUE)
  }

  if (!inherits(value, 'IRangesList')) {
    stop("Illegal replacement value")
  }

  if (missing(i)) {
    i <- 1:length(x)
  }
  
  if (!missing(j)) {
    warning("Not sure how to handle the j argument gracefully")
    ## if (is.null(elementMetadata(x))) {
    ##   j <- NULL
    ## } else {
    ##   j <- 1:ncol(elementMetadata(x))
    ## }
  }

  ## Replacing elements in a compressed list is all-sorts-of whacky, so I'm
  ## "uncompressing" and manipulating "Normal" lists, then compress it back
  ## into the expected CompressedNormalIRangesList
  nirl <- as.list(x@cnirl)
  value <- lapply(value, as, 'NormalIRanges')
  nirl[i] <- value
  cnirl <- as(IRangesList(nirl), 'CompressedNormalIRangesList')
  x@cnirl <- cnirl
  x
})

setMethod("gwidth", c(x="GappedRanges"),
function(x, mind.the.gap=TRUE, ...) {
  if (mind.the.gap) {
    sum(width(x@cnirl))
  } else {
    width(x)
  }
})

setMethod("ranges", c(x="GappedRanges"),
function(x, mind.the.gap=FALSE, ...) {
  if (mind.the.gap) {
    x@cnirl
  } else {
    ## range(x@cnirl) doesn't work because of some "not NormalIRanges" result
    inner <- x@cnirl
    class(inner) <- 'CompressedIRangesList'
    unlist(range(inner))
  }
})

setMethod("findOverlaps", c("Ranges", "GappedRanges"),
function(query, subject, maxgap=0L, minoverlap=1L,
         type=c("any", "start", "end", "within", "equal"),
         select=c("all", "first", "last", "arbitrary"),
         mind.the.gap=TRUE, usage.warning=TRUE, ...) {
  findOverlaps(query, ranges(subject, mind.the.gap), maxgap=maxgap,
               minoverlap=minoverlap, type=type, select=select)
})

## An empty intersection returns an IRanges of start=0, end=0. This isn't
## correct, but since this is for only genomic coordinates, I'm doing this
## for consistency's sake (for some definition of consistency!).
setMethod("intersect", c(x="IRanges", y="GappedRanges"),
function(x, y) {
  intersected <- seqapply(ranges(y, mind.the.gap=TRUE), function(.ranges) {
    i <- intersect(x, .ranges)
    if (length(i) == 0L) {
      i <- IRanges(0, 0)
    }
    i
  })
  y@cnirl <- as(intersected, 'CompressedNormalIRangesList')
  y
})

setMethod("as.data.frame", c(x="GappedRanges"),
function(x, row.names=NULL, optional=FALSE, ...) {
  DF <- as.data.frame(ranges(x))
  DF$width.nogaps <- gwidth(x)
  DF$ngap <- ngap(x)
  DF
})

setMethod("show", c(object="GappedRanges"),
function(object) {
  show(as.data.frame(object))
})

