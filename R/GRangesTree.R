## An IntervalTree for GRanges objects
setClass("GRangesTree",
  representation(trees='list'),
  prototype=list())

setValidity("GRangesTree", function(object) {
  if (is.null(names(object@trees))) {
    return(paste("`trees` list must be names"))
  }
  TRUE
})

setAs("GRanges", "GRangesTree", function(from) {
  .seqranges <- split(from, seqnames(from))
  trees <- lapply(.seqranges, function(seqranges) {
    if (length(seqranges) == 0L) {
      return(NULL)
    }
    .stranges <- split(seqranges, strand(seqranges))
    s <- lapply(.stranges, function(stranges) {
      if (length(stranges) == 0L) {
        return(NULL)
      }
      IntervalTree(ranges(stranges))
    })
    names(s) <- names(.stranges)
    s <- s[!sapply(s, is.null)]
  })
  names(trees) <- names(.seqranges)
  trees <- trees[!sapply(trees, is.null)]
  new('GRangesTree', trees=trees)
})

setAs("CompressedReads", "GRangesTree", function(from) {
  from <- as(from, "GRanges")
  as(from, "GRangesTree")
})

setMethod("findOverlaps", c("GRanges", "GRangesTree"),
function(query, subject, maxgap=0L, minoverlap=1L,
         type=c("any", "start", "end"),
         select=c("all", "first")) {
  ##
  type <- match.arg(type)
  select <- match.arg(select)
  stop("Not implemented")
})

setMethod("subsetByOverlaps", c("GRangesTree", "GRanges"),
function(query, subject, maxgap=0L, minoverlap=1L,
         type=c('any', 'start', 'end')) {
  .seqranges <- split(subject, seqnames(subject))
  seqover <- lapply(.seqranges, function(seqranges) {
    if (length(seqranges) == 0L) {
      return(NULL)
    }
    seqname <- as.character(seqnames(seqranges)[1])
    chr.tree <- query@trees[[seqname]]
    if (is.null(chr.tree)) {
      return(NULL)
    }
    .stranges <- split(seqranges, strand(seqranges))
    s <- lapply(.stranges, function(stranges) {
      if (length(stranges) == 0L) {
        return(NULL)
      }
      strandname <- as.character(strand(stranges)[1])
      if (strandname == '*') {
        .ranges <- c(IRanges(chr.tree[['+']]), IRanges(chr.tree[['-']]),
                     IRanges(chr.tree[['*']]))
        str.tree <- IntervalTree(.ranges)
      } else {
        str.tree <- chr.tree[[strandname]]
      }
      if (is.null(str.tree)) {
        return(NULL)
      }
      o <- findOverlaps(ranges(stranges), str.tree, maxgap=maxgap,
                        minoverlap=minoverlap)
      if (length(o) == 0L) {
        return(NULL)
      }
      GRanges(seqnames=seqname,
              ranges=as(str.tree[subjectHits(o)], 'IRanges'),
              strand=strandname)
    })
    names(s) <- names(.stranges)
    s <- s[!sapply(s, is.null)]
  })
  seqover <- seqover[!sapply(seqover, is.null)]
  if (length(seqover) == 0L) {
    GRanges()
  } else {
    seqover <- unlist(seqover)
    names(seqover) <- NULL
    do.call(c, seqover)
  }
})

