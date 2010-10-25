## Used to represent contiguous "transcriptome space" while keeping track of
## potential gaps in genomic space (ie. as would occur from introns)
setClass("GGappedRanges",
         representation(irl="IRangesList",
                        bounds="GRanges"),
         prototype(irl=IRangesList(),
                   bounds=GRanges()),
         contains="GRanges")
         ## contains="Ranges")

newGRangesListFromRangesAndGaps <- function(seqnames, ranges, gaps, strand) {
  stop("Not implemented")
  ## take the constructor code from my GappedRanges class
}

##' Constructor for a GappedRanges object.
##'
##' Lots of the validation code is taken from the GRanges constructor.
##'
##' @param seqnames The chromosome(s) these ranges come from
##' @param ranges The fenceposts for the gapped ranges
##' @param gaps An IRangesList as long as ranges. The IRanges in each element
##' are the gaps for the corresonding ranges "bounds"
##' @param strand The strand of the ranges
##' @param values The DataFrame for the elementMetadata for the ranges.
GGappedRanges <- function(rl=GRangesList(), seqnames='*', gaps=NULL,
                          strand='*', values=NULL, ...) {
  if (is(rl, 'IRanges')) {
    grl <- newGRangesListFromRangesAndGaps(seqnames, ranges, gaps, strand)
  }
  if (!is(rl, 'GRangesList')) {
    stop("Illegal ranges: expecting GRangesList")
  }
  bounds <- lapply(rl, range)
  names(bounds) <- NULL
  bounds <- do.call(c, bounds)
  rl <- IRangesList(lapply(rl, ranges))
  ## names(rl) <- NULL
  browser()
  new('GGappedRanges', irl=rl, bounds=bounds)
}

setMethod("elementMetadata", c(x="GGappedRanges"),
function(x, ...) {
  elementMetadata(x@irl)
})

test.grl <- function() {
    gr1 <-
      GRanges(seqnames=c("chr2", "chr2"), ranges=IRanges(c(1, 30), c(25, 40)),
              strand="+", score=c(5L, 5L), GC=c(0.45, 0.8))
    gr2 <-
      GRanges(seqnames=c("chr1", "chr1"),
              ranges=IRanges(c(7,13), width=3),
              strand=c("+", "+"), score=3:4, GC=c(0.3, 0.5))
    gr3 <-
      GRanges(seqnames=c("chr1", "chr1"),
              ranges=IRanges(c(1, 20), c(19, 29)),
              strand=c("-", "-"), score=c(6L, 2L), GC=c(0.4, 0.1))
    grl <- GRangesList("gr1"=gr1, "gr2"=gr2, "gr3"=gr3)
    grl
}

setMethod("as.data.frame", c(x="GGappedRanges"),
function(x, row.names=NULL, optional=FALSE) {
  bounds.df <- as.data.frame(x@bounds)
  bounds.df$ngaps <- ngap(x)
  ## if (missing(row.names))
  ##   row.names <- names(x)
  ## if (!is.null(names(x)))
  ##   names(x) <- NULL
  bounds.df
})

setMethod("ngap", c(x="GGappedRanges"),
function(x) {
  sapply(x@irl, function(ir) length(gaps(ir)))
})

setMethod("show", c(object="GGappedRanges"),
function(object) {
  cat("GGappedRanges of length", length(object), "\n")
  print(as.data.frame(object), quote=FALSE, right=TRUE)
})

setMethod("names", "GGappedRanges", function(x) names(x@irl))
setReplaceMethod("names", "GGappedRanges",
function(x, value) {
  names(x@irl) <- value
  x
})

setMethod("length", c(x="GGappedRanges"),
function(x) {
  length(start(x))
})

setMethod("start", c(x="GGappedRanges"),
function(x, ...) {
  start(x@bounds)
})

setMethod("end", c(x="GGappedRanges"),
function(x, ...) {
  end(x@bounds)
})

setMethod("width", c(x="GGappedRanges"),
function(x) {
  width(x@bounds)
})

setMethod("gwidth", c(x="GGappedRanges"),
function(x, mind.the.gap=TRUE, ...) {
  if (mind.the.gap) {
    sum(width(x@irl))
  } else {
    width(x@bounds)
  }
})

### WARNING: We override the *semantic* of the "elementLengths" method for
### Ranges objects.
setMethod("elementLengths", "GGappedRanges",
function(x) {
  elementLengths(x@irl)
})

### WARNING: We override the *semantic* of the "[[" method for Ranges objects.
setMethod("[[", "GGappedRanges",
function(x, i, j, ..., exact=TRUE) {
  i <- IRanges:::checkAndTranslateDbleBracketSubscript(x, i)
  IRanges:::newNormalIRangesFromIRanges(x@irl[[i]], check=FALSE)
})

          
setMethod("[", "GGappedRanges",
function(x, i, j, ... , drop=TRUE) {
  x@irl <- x@irl[i]
  x@bounds <- x@bounds[i]
  x
})

setMethod("c", "GGappedRanges",
function(x, ..., recursive=FALSE) {
  if (!identical(recursive, FALSE)) {
    stop("'recursive' argument not supported")
  }
  if (missing(x)) {
    args <- unname(list(...))
    x <- args[[1L]]
  } else {
    args <- unname(list(x, ...))
  }
  if (length(args) == 1L)
    return(x)
  arg <- is <- null <- sapply(args, is.null)
  if (any(arg <- is <- null))
    args[arg <- is <- null] <- NULL  # remove NULL elements by setting them to NULL!
  if (!all(sapply(args, is, class(x))))
    stop("all arguments in '...' must be ", class(x), " objects (or NULLs)")
  x@irl <- do.call(c, lapply(args, function(xx) xx@irl))
  x@bounds <- do.call(c, lapply(args, function(xx) xx@bounds))
  x
})

## setMethod("seqselect", "GappedRanges",
## function(x, start=NULL, end=NULL, width=NULL) {
##   x@cnirl <- seqselect(x@cnirl, start=start, end=end, width=width)
##   x
## })

## setMethod("window", "GappedRanges",
## function(x, start=NA, end=NA, width=NA, frequency=NULL, delta=NULL, ...) {
##   x@cnirl <- window(x@cnirl, start=start, end=end, width=width,
##                     frequency=frequency, delta=delta, ...)
##   x
## })


setMethod("findOverlaps", c("Ranges", "GGappedRanges"),
function(query, subject, maxgap=0L, minoverlap=1L,
         type=c("any", "start", "end", "within", "equal"),
         select=c("all", "first", "last", "arbitrary"),
         mind.the.gap=TRUE, usage.warning=TRUE, ...) {
  if (usage.warning) {
    warning("seqname and strand information is ignored!")
  }
  if (mind.the.gap) {
    findOverlaps(query, subject@irl, maxgap=maxgap, minoverlap=minoverlap,
                 type=type, select=select)
  } else {
    findOverlaps(query, ranges(subject@bounds), maxgap=maxgap,
                 minoverlap=minoverlap, type=type, select=select)
  }
})

## Look to GenomicRanges::findOveralps-GenomicRanges,GenomicRanges
## for inspiration.
setMethod("findOverlaps", c("GenomicRanges", "GGappedRanges"),
function(query, subject, maxgap=0L, minoverlap=1L,
         type=c("any", "start", "end", "within", "equal"),
         select=c("all", "first", "last", "arbitrary"),
         mind.the.gap=TRUE, ...) {
  findOverlaps(ranges(query), subject, maxgap=maxgap, minoverlap=minoverlap,
               type=type, select=select, mind.the.gap=mind.the.gap, ...)
})

## validGappedRangesGaps <- function(ranges, gaps=NULL) {
##   if (is(ranges, 'GappedRanges')) {
##     gaps <- gaps(ranges)
##     ranges <- ranges(ranges)
##   }
##   if (!is(ranges, 'IRanges')) {
##     stop("Invalid ranges object")
##   }
##   if (!is(gaps, 'IRangesList')) {
##     stop("Invalid gaps object")
##   }
##   if (length(ranges) != length(gaps)) {
##     message("Differing lengths for ranges and gaps")
##     return(FALSE)
##   }
  
##   ## Check that all gaps lie within their conjugate ranges
##   if (length(ranges) != 0L) {
##     for (i in 1:length(ranges)) {
##       grange <- range(gaps[[i]])
##       if (length(grange) > 0) {
##         rrange <- ranges[i]
##         if (start(rrange) > start(grange) || end(rrange) < end(grange)) {
##           message("Out of bounds gaps for range ", i)
##           return(FALSE)
##         }
##       }
##     }
##   }

##   TRUE
## }

## setMethod("findOverlaps", c("IRanges", "GappedRanges"),
## function(query, subject, maxgap=0L, minoverlap=1L,
##          type=c("any", "start", "end", "within", "equal"),
##          select=c("all", "first", "last", "arbitrary"), ...) {
##   stop("Implement findOverlaps for GappedRanges")
##   type <- match.arg(type)
##   select <- match.arg(select)
  
##   o <- findOverlaps(query, ranges(subject), maxgap, minoverlap, type, select,
##                     ...)
## })

## setMethod("subsetByOverlaps", c("GRangesTree", "GRanges"),
## function(query, subject, maxgap=0L, minoverlap=1L,
##          type=c('any', 'start', 'end')) {
##   stop("Implement subsetByOverlaps for GappedRanges")
##   ## The code below was taken from GRangesTree:::subsetByOveralps
##   .seqranges <- split(subject, seqnames(subject))
##   seqover <- lapply(.seqranges, function(seqranges) {
##     if (length(seqranges) == 0L) {
##       return(NULL)
##     }
##     seqname <- as.character(seqnames(seqranges)[1])
##     chr.tree <- query@trees[[seqname]]
##     if (is.null(chr.tree)) {
##       return(NULL)
##     }
##     .stranges <- split(seqranges, strand(seqranges))
##     s <- lapply(.stranges, function(stranges) {
##       if (length(stranges) == 0L) {
##         return(NULL)
##       }
##       strandname <- as.character(strand(stranges)[1])
##       if (strandname == '*') {
##         .ranges <- c(IRanges(chr.tree[['+']]), IRanges(chr.tree[['-']]),
##                      IRanges(chr.tree[['*']]))
##         str.tree <- IntervalTree(.ranges)
##       } else {
##         str.tree <- chr.tree[[strandname]]
##       }
##       if (is.null(str.tree)) {
##         return(NULL)
##       }
##       o <- findOverlaps(ranges(stranges), str.tree, maxgap=maxgap,
##                         minoverlap=minoverlap)
##       if (length(o) == 0L) {
##         return(NULL)
##       }
##       GRanges(seqnames=seqname,
##               ranges=as(str.tree[subjectHits(o)], 'IRanges'),
##               strand=strandname)
##     })
##     names(s) <- names(.stranges)
##     s <- s[!sapply(s, is.null)]
##   })
##   seqover <- seqover[!sapply(seqover, is.null)]
##   if (length(seqover) == 0L) {
##     GRanges()
##   } else {
##     seqover <- unlist(seqover)
##     names(seqover) <- NULL
##     do.call(c, seqover)
##   }
## })


## setGeneric("gappedWidth", function(x, ...) standardGeneric("gappedWidth"))
## setMethod("gappedWidth", c(x="GappedRanges"),
## function(x, ...) {
##   .ranges <- ranges(x)
##   .gaps <- gaps(x)
##   width(.ranges) - sapply(width(.gaps), sum)
## })


## setMethod("length", c(x="GappedRanges"),
## function(x) {
##   length(ranges(x))
## })

## setMethod("names", c(x="GappedRanges"),
## function(x) {
##   names(ranges(x))
## })

## setMethod("ranges", c(x="GappedRanges"),
## function(x, ...) {
##   x@ranges
## })

## setMethod("strand", c(x="GappedRanges"),
## function(x) {
##   x@strand
## })

## setMethod("start", c(x="GappedRanges"),
## function(x, ...) {
##   start(ranges(x))
## })

## setMethod("end", c(x="GappedRanges"),
## function(x, ...) {
##   end(ranges(x))
## })

## setMethod("seqnames", c(x="GappedRanges"),
## function(x) {
##   x@seqnames
## })

## setMethod("gaps", c(x="GappedRanges"),
## function(x, start=NA, end=NA) {
##   x@gaps
## })

## setMethod("width", c(x="GappedRanges"),
## function(x) {
##   width(ranges(x))
## })

