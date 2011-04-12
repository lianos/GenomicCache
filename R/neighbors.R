##' Provides information about up/downstream neighbors between ranges objects.
##'
##' For each event in `from`, the index into `to` for the closest up and
##' downstream events is returned, as well as their respective distances.
##'
##' @param from A *Ranges object
##' @param to A *Ranges object. Defaults to \code{from}
setGeneric("neighbors",
function(from, to, include.overlap=!is.null(to), ...) {
  standardGeneric("neighbors")
})

setMethod("neighbors", c(from="GRanges", to="missing"),
function(from, to, include.overlap=FALSE, both.strands=TRUE, ...) {
  neighbors(from, from, include.overlap=include.overlap,
            both.strands=both.strands, ...)
})

setMethod("neighbors", c(from="GRanges", to="GRanges"),
function(from, to, include.overlaps, both.strands=TRUE, ...) {
  ## downstream
  ## downstream.from has the index of `to` that are immediately downstream from
  ## `from`.
  downstream.from <- precede(from, to)
  calc <- !is.na(downstream.from)

  down.dists <- rep(NA_integer_, length(from))
  use.to <- to[downstream.from[calc]]
  use.from <- from[calc]
  down.dists[calc] <- ifelse(as.logical(strand(use.from) == '-'),
                             start(use.from) - end(use.to),
                             start(use.to) - end(use.from))
  ## } else {
  ##   down.dists[calc] <- start(use.to) - end(use.from)
  ## }


  ## upstream
  upstream.from <- follow(from, to)
  calc <- !is.na(upstream.from)

  up.dists <- rep(NA_integer_, length(from))
  use.to <- to[upstream.from[calc]]
  use.from <- from[calc]

  up.dists[calc] <- ifelse(as.logical(strand(use.from) == '-'),
                           start(use.to) - end(use.from),
                           start(use.from) - end(use.to))
  ## } else {
  ##   up.dists[calc] <- start(use.from) - end(use.to)
  ## }

  dists <- data.frame(upstream=up.dists, upstream.idx=upstream.from,
                      downstream=down.dists, downstream.idx=downstream.from)

  if (both.strands) {
    new.to <- to
    strand(new.to) <- ifelse(as.logical(strand(new.to) != '-'), '-', '+')
    opp.dists <- neighbors(from, new.to, include.overlap=include.overlap,
                           both.strands=FALSE)
    colnames(opp.dists) <- paste('opp', colnames(opp.dists), sep=".")
    dists <- cbind(dists, opp.dists)
  }

  dists
})

setMethod("neighbors", c(from="GRanges", to="IRanges"),
function(from, to, include.overlaps,...) {
  neighbors(ranges(from), to, include.overlaps, ...)
})

setMethod("neighbors", c(from="IRanges", to="GRanges"),
function(from, to, include.overlaps, ...) {
  seqname <- as.character(unique(seqnames(to)))
  if (length(seqname) != 1) {
    stop("A `to` object across multiple seqnames isn't allowed")
  }
  from <- GRanges(seqname, from, '+')
  d.fwd <- neighbors(from, to, include.overlaps=include.overlaps,
                     both.strands=FALSE)
  colnames(d.fwd) <- paste('fwd', colnames(d.fwd), sep='.')

  browser()
  strand(from) <- '-'
  d.rev <- neighbors(from, to, include.overlaps=include.overlaps,
                     both.strands=FALSE)
  colnames(d.rev) <- paste('rev', colnames(d.rev), sep=".")

  cbind(d.fwd, d.rev)
})

setMethod("neighbors", c(from="IRanges", to="missing"),
function(from, to, include.overlaps, ...) {
  neighbors(from, from, include.overlaps, ...)
})

setMethod("neighbors", c(from="IRanges", to="IRanges"),
function(from, to, include.overlaps, ...) {
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
  dists
})
