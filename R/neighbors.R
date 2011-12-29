##' Provides information about up/downstream neighbors between ranges objects.
##'
##' For each event in `from`, the index into `to` for the closest up and
##' downstream events is returned, as well as their respective distances.
##'
##' @param from A *Ranges object
##' @param to A *Ranges object. Defaults to \code{from}
setGeneric("neighbors", signature=c('from', 'to'),
function(from, to, ...) {
  standardGeneric("neighbors")
})

setMethod("neighbors", c(from="GRanges", to="missing"),
function(from, to, both.strands=TRUE, ...) {
  neighbors(from, from, both.strands=both.strands, to.missing=TRUE, ...)
})

setMethod("neighbors", c(from="GRanges", to="GRanges"),
function(from, to, both.strands=TRUE, ...) {
  ## downstream
  ## downstream.from has the index of `to` that are immediately downstream from
  ## `from`.
  args <- list(...)
  downstream.from <- precede(from, to)
  calc <- !is.na(downstream.from)

  down.dists <- rep(NA_integer_, length(from))
  use.to <- to[downstream.from[calc]]
  use.from <- from[calc]
  down.dists[calc] <- ifelse(as.logical(strand(use.from) == '-'),
                             start(use.from) - end(use.to),
                             start(use.to) - end(use.from))

  ## upstream
  upstream.from <- follow(from, to)
  calc <- !is.na(upstream.from)

  up.dists <- rep(NA_integer_, length(from))
  use.to <- to[upstream.from[calc]]
  use.from <- from[calc]

  up.dists[calc] <- ifelse(as.logical(strand(use.from) == '-'),
                           start(use.to) - end(use.from),
                           start(use.from) - end(use.to))

  dists <- data.frame(upstream=up.dists, upstream.idx=upstream.from,
                      downstream=down.dists, downstream.idx=downstream.from)

  if (both.strands) {
    new.to <- to
    strand(new.to) <- ifelse(as.logical(strand(new.to) != '-'), '-', '+')
    opp.dists <- neighbors(from, new.to, both.strands=FALSE, ...)
    colnames(opp.dists) <- paste('opp', colnames(opp.dists), sep=".")
    dists <- cbind(dists, opp.dists)
  }

  if (is.null(args$skip.check) || !args$skip.check) {
    if (!is.null(args$to.missing) && args$to.missing) {
      o <- findOverlaps(from, type='any', ignoreSelf=TRUE,
                        ignoreRedundant=TRUE)
    } else {
      o <- findOverlaps(from, to, type='any')
    }
    if (length(o) > 0) {
      warning("Overlaps were found between from/to, but they are ignored",
              immediate.=TRUE)
    }
  }

  dists
})

setMethod("neighbors", c(from="GRanges", to="IRanges"),
function(from, to, ...) {
  neighbors(ranges(from), to, ...)
})

setMethod("neighbors", c(from="IRanges", to="GRanges"),
function(from, to, ...) {
  seqname <- unique(as.character(seqnames(to)))
  if (length(seqname) != 1) {
    stop("A `to` object across multiple seqnames isn't allowed")
  }
  from <- GRanges(seqname, from, '+')
  d.fwd <- neighbors(from, to, both.strands=FALSE, ...)
  colnames(d.fwd) <- paste('fwd', colnames(d.fwd), sep='.')

  browser()
  strand(from) <- '-'
  d.rev <- neighbors(from, to, both.strands=FALSE, ...)
  colnames(d.rev) <- paste('rev', colnames(d.rev), sep=".")

  cbind(d.fwd, d.rev)
})

setMethod("neighbors", c(from="IRanges", to="missing"),
function(from, to, ...) {
  neighbors(from, from, to.missing=TRUE, ...)
})

setMethod("neighbors", c(from="IRanges", to="IRanges"),
function(from, to, ...) {
  ## downstream
  ## downstream.from has the index of `to` that are immediately downstream from
  ## `from`.
  downstream.from <- precede(from, to)
  calc <- !is.na(downstream.from)

  down.dists <- rep(NA_integer_, length(from))
  use.to <- to[downstream.from[calc]]
  use.from <- from[calc]
  down.dists[calc] <- start(use.to) - end(use.from)

  ## upstream
  upstream.from <- follow(from, to)
  calc <- !is.na(upstream.from)

  up.dists <- rep(NA_integer_, length(from))
  use.to <- to[upstream.from[calc]]
  use.from <- from[calc]
  up.dists[calc] <- start(use.from) - end(use.to)

  dists <- data.frame(upstream=up.dists, upstream.idx=upstream.from,
                      downstream=down.dists, downstream.idx=downstream.from)

  args <- list(...)
  if (is.null(args$skip.check) || !args$skip.check) {
    if (!is.null(args$to.missing) && args$to.missing) {
      o <- findOverlaps(from, type='any', ignoreSelf=TRUE,
                        ignoreRedundant=TRUE)
    } else {
      o <- findOverlaps(from, to, type='any')
    }
    if (length(o) > 0) {
      warning("Overlaps were found between from/to, but they are ignored",
              immediate.=TRUE)
    }
  }


  dists
})
