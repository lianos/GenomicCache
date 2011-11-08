#' @nord
matchGFGeneCollapse <- function(collapse, gene=NULL) {
  choices <- c('cover', 'constitutive', 'first', 'longest', 'shortest')
  if (is.numeric(collapse)) {
    collapse <- as.integer(collapse)
    if (collapse != 1L && !is.null(gene) && length(transcripts(x) < collapse)) {
      stop("Transcript chosen for `exon.type` is out of bounds")
    }
  } else {
    collapse <- match.arg(collapse, choices)
  }
  if (collapse == 'first') {
    collapse <- 1L
  }
  collapse
}

##' Returns an annotated GRanges object by "compressing" the exons
##' across isoforms into one prototypical / idealized version of the gene..
##'
##' @param x The \code{GFGene} object
##' @param by How the caller wants the "idealized" gene to be calculated.
##' \code{all} exons? Just the coding (\code{cds}) exons? etc.
##' @param collapse How are the exons summarized? Is it the union of all exons
##' accross transcripts (\code{cover}), the constitutively spliced portions?
##' \code{constitutive}, just the \code{first} transcript? The caller can
##' also set this to an integer indicating which transcript to pick of them all.
##' @param cds.cover Some genes have isoforms with coding regions over utr's
##' over other isoforms. If set the \code{min} then the \code{cds} boundaries
##' are the narrowest possible (never overlapping a \code{utr}). \code{max} 
##' will allos the \code{cds} to overlap \code{utr}'s.
##'
##' @return An annotated \code{GRanges} object, whose accomponied \code{values}
##' has an \code{exon.anno} collumn identifying which corresponding rantges are
##' \code{cds, utr3, utr5, utr}
##'
##' If \code{by} is either \code{utr3, utr5}, then it returns a
##' \code{GRangesList} containing ranges that are split by intermediadary
##' coding regions. The utr's are split this way to provide a "best guess"
##' of which utrs can be considered as a whole (and not spanning separate
##' transcripts).
##'
##' If this last point seems strange to you, you should just use \code{by='all'}
##' and \code{subset} the returned object by, for example,
##' \code{values(...)$exon.ann == 'utr3}.
setGeneric("idealized", 
function(x, by=c('all', 'cds', 'utr5', 'utr3'),
         collapse=c('cover', 'constitutive', 'first', 'longest', 'shortest'),
         cds.cover=c('min', 'max'), which.chr=NULL, flank.up=0L, flank.down=0L,
         ...) {
  standardGeneric("idealized")
})

setMethod("idealized", c(x="GFGene"),
function(x, by, collapse, cds.cover, which.chr, flank.up, flank.down, ...) {
  by <- match.arg(by)
  cds.cover <- match.arg(cds.cover)
  collapse <- matchGFGeneCollapse(collapse, x)
  reducef <- switch(collapse, constitutive='interset', 'union')
  metaD <- list()
  
  ## g.strand <- strand(x)
  g.strand <- as.vector(strand(transcripts(x, which.chr=which.chr)[[1]][1]))
  
  args <- list(...)
  force.eval <- args$force.eval
  
  var.name <- paste('idealized', collapse,
                    sprintf('by-%s', by),
                    sprintf('cdscover-%s', cds.cover),
                    sprintf('utr5*-%d', flank.up),
                    sprintf('utr3*-%d', flank.down),
                    sep=".")
  if (!is.null(which.chr)) {
    var.name <- paste(var.name, which.chr, sep=".")
  }
  
  cacheFetch(x, var.name, {
    .idealizedBy <- switch(by, all=.idealizedByAll, cds=.idealizedByCds,
                           utr5=.idealizedByUtr5, utr3=.idealizedByUtr3,
                           longest=.idealizedByLongest)
    is.pc <- isProteinCoding(x, which.chr=which.chr)
    is.pc <- names(is.pc)[is.pc]
    probably.rna <- length(is.pc) == 0

    xcripts <- transcripts(x, which.chr=which.chr)
    
    if (probably.rna) {
      ## A gene with non-protein coding transcripts is tagged as `utr` across
      ## across all of its exons in this scenario.
      ## TODO: Cleanup duplicated code in probably.rna if/else block
      if (collapse %in% c('longest', 'shortest')) {
        lens <- sapply(xcripts, function(xc) {
          width(range(ranges(xc)))
        })
        collapse <- switch(collapse, longest=which.max(lens), which.min(lens))
        metaD$tx_name <- values(xcripts)$tx_name[collapse]
      }
      if (is.numeric(collapse)) {
        ranges <- ranges(xcripts[[collapse]])
      } else {
        ranges <- Reduce(reducef, lapply(xcripts, ranges))        
      }
      exon.anno <- c(rep('utr', length(ranges)))
    } else {
      take.pc <- function(.list, namez) {
        .list[names(.list) %in% namez]
      }
      takeOrEmpty <- function(x, idx) {
        if (length(x) < idx) IRanges() else x[[idx]]
      }
      
      .cds <- lapply(take.pc(cds(x, which.chr=which.chr), is.pc), ranges)
      .utr3 <- lapply(take.pc(utr3(x, which.chr=which.chr), is.pc), ranges)
      .utr5 <- lapply(take.pc(utr5(x, which.chr=which.chr), is.pc), ranges)

      if (collapse %in% c('longest', 'shortest')) {
        lens <- sapply(take.pc(xcripts, is.pc), function(xc) {
          width(range(ranges(xc)))
        })
        collapse <- switch(collapse, longest=which.max(lens), which.min(lens))
        metaD$tx_name <- values(xcripts)$tx_name[collapse]
      }
      
      if (is.numeric(collapse)) {
        ## Perhaps the Nth transcript doesn't have all of the expected segments
        ## in which case, we want to return an empty IRanges
        .exons <- list(utr5=takeOrEmpty(.utr5, collapse),
                       cds=takeOrEmpty(.cds, collapse),
                       utr3=takeOrEmpty(.utr3, collapse))
      } else {
        ## `Reduce`-ing over an empty list (which can result from a gene not
        ## having a .utr3) will produce a NULL in insteand of an empty IRanges
        ## object, which will throw errors in the setdiff operations below
        null2empty <- function(x) if (is.null(x)) IRanges() else x
        .cds <- null2empty(Reduce(reducef, .cds))
        .utr3 <- null2empty(Reduce(reducef, .utr3))
        .utr5 <- null2empty(Reduce(reducef, .utr5))
        
        if (cds.cover == 'min') {
          .cds <- setdiff(.cds, union(.utr3, .utr5))
        } else {
          .utr3 <- setdiff(.utr3, .cds)
          .utr5 <- setdiff(.utr5, .cds)
        }

        .exons <- .idealizedBy(.utr5, .cds, .utr3, cds.cover)
      }
      ranges <- c(.exons$utr5, .exons$cds, .exons$utr3)
      exon.anno <- c(rep('utr5', length(.exons$utr5)),
                     rep('cds', length(.exons$cds)),
                     rep('utr3', length(.exons$utr3)))
    }

    chr <- if (!is.null(which.chr)) which.chr else chromosome(x)[1]
    
    if (probably.rna && by == 'cds') {
      gr <- GRanges()
      values(gr) <- DataFrame(exon.anno=character())
    } else if (probably.rna || by == 'all') {
      gr <- GRanges(seqnames=rep(chr, length(ranges)),
                    ranges=ranges,
                    strand=rep(g.strand, length(ranges)))
      values(gr) <- DataFrame(exon.anno=exon.anno)
      gr <- gr[order(ranges(gr))]
    } else {
      gr <- lapply(.exons, function(e) {
        g <- GRanges(seqnames=rep(chr, length(e)),
                     ranges=e,
                     strand=rep(g.strand, length(e)))
        values(g) <- DataFrame(exon.anno=rep(by, length(g)))
        if (length(g) > 0) {
          g[order(ranges(g))]
        }
        g
      })
      gr <- GRangesList(gr)
    }

    if (flank.up > 0) {
      if (probably.rna) {
        
      } else {
        take <- which(values(gr)$exon.anno == 'utr5')
        if (length(take) > 0) {
          if (g.strand == '+') {
            ext.end <- start(gr[take[1L]]) - 1L
            ext.start <- ext.end - flank.up + 1
            utr.ext <- GRanges(seqnames=chr, strand=g.strand,
                               ranges=IRanges(start=ext.start, end=ext.end))
            values(utr.ext) <- DataFrame(exon.anno='utr5*')
            gr <- c(utr.ext, gr)
          } else {
            ext.start <- end(gr[take[length(take)]]) + 1L
            ext.end <- ext.start + flank.up - 1L
            utr.ext <- GRanges(seqnames=chr, strand=g.strand,
                               ranges=IRanges(start=ext.start, end=ext.end))
            values(utr.ext) <- DataFrame(exon.anno='utr5*')
            gr <- c(gr, utr.ext)
          }
        }
      }
    }
    
    if (flank.down > 0) {
      if (probably.rna) {
        
      } else {
        take <- which(values(gr)$exon.anno == 'utr3')
        if (length(take) > 0) {
          if (g.strand == '+') {
            ext.start <- end(gr[take[length(take)]]) + 1L
            ext.end <- ext.start + flank.down - 1L
            utr.ext <- GRanges(seqnames=chr, strand=g.strand,
                               ranges=IRanges(start=ext.start, end=ext.end))
            values(utr.ext) <- DataFrame(exon.anno='utr3*')
            gr <- c(gr, utr.ext)
          } else {
            ext.end <- start(gr[take[1L]]) - 1L
            ext.start <- ext.end - flank.down + 1
            utr.ext <- GRanges(seqnames=chr, strand=g.strand,
                               ranges=IRanges(start=ext.start, end=ext.end))
            values(utr.ext) <- DataFrame(exon.anno='utr3*')
            gr <- c(utr.ext, gr)
          }
        }
      }
    }
    metadata(gr) <- metaD
    gr
  }, force.eval=force.eval)
})


.idealizedByUtr <- function(utr, .cds) {
  if (length(utr) == 0) {
    return(list(IRanges()))
  }

  cgaps <- gaps(.cds)
  mm <- matchMatrix(findOverlaps(cgaps, utr))
  
  if (nrow(mm) == 0) {
    utr <- list(utr)
  } else {
    qh <- unique(mm[, 'query'])
    sh <- unique(mm[, 'subject'])
    utr.split <- lapply(qh, function(i) {
      take <- mm[, 'query'] == i
      utr[mm[take, 'subject']]
    })
    utr.rest <- setdiff(1:length(utr), sh)
    utr <- list(c(utr.split, list(utr[utr.rest])))
  }
  
  utr
}

## Will split consecutive 3' UTRs into 2 GRanges if there is a cds region
## between the "idealized" 3'utr stretch
.idealizedByUtr3 <- function(.utr5, .cds, .utr3, cds.cover) {
  .idealizedByUtr(.utr3, .cds)
}

.idealizedByUtr5 <- function(.utr5, .cds, .utr3, cds.cover) {
  .idealizedByUtr(.utr5, .cds)
}

.idealizedByCds <- function(.utr5, .cds, .utr3, cds.cover) {
  warning("idealized(gene, by='cds') is not really implemented yet.")
  list(utr5=IRanges(), cds=.cds, utr3=IRanges())
}

.idealizedByAll <- function(.utr5, .cds, .utr3, cds.cover) {
  list(utr5=.utr5, cds=.cds, utr3=.utr3)
}


