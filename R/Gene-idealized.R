matchGFGeneCollapse <- function(x, collapse) {
  if (is.numeric(collapse)) {
    collapse <- as.integer(collapse)
    if (collapse > length(transcripts(x)) || collapse < 1L) {
      stop("Transcript chosen for `exon.type` is out of bounds")
    }
  } else {
    collapse <- match.arg(collapse, c('cover', 'constitutive', 'first'))
  }
  if (collapse == 'first') {
    collapse <- 1L
  }
  collapse
}

##' Returns an annotated GRanges object by "compressing" the exons
##' across isoforms into one prototypical / idealized version of the gene..
##'
##' @param The \code{GFGene} object
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
         collapse=c('cover', 'constitutive', 'first'),
         cds.cover=c('min', 'max'), flank.up=0L, flank.down=0L, ...) {
  standardGeneric("idealized")
})

setMethod("idealized", c(x="GFGene"),
function(x, by, collapse, cds.cover, flank.up, flank.down, ...) {
  by <- match.arg(by)
  cds.cover <- match.arg(cds.cover)
  collapse <- matchGFGeneCollapse(x, collapse)
  reducef <- switch(collapse, cover='union', 'intersect')
  g.strand <- strand(x)
  
  var.name <- paste('idealized', collapse,
                    sprintf('by-%s', by),
                    sprintf('cdscover-%s', cds.cover),
                    sprintf('utr5*-%d', flank.up),
                    sprintf('utr3*-%d', flank.down),
                    sep=".")
  
  cacheFetch(x, var.name, {
    .idealizedBy <- switch(by, all=.idealizedByAll, cds=.idealizedByCds,
                           utr5=.idealizedByUtr5, utr3=.idealizedByUtr3)

    probably.rna <- (all(!isProteinCoding(x)) ||
                     (any(!isProteinCoding(x)) && cds.cover == 'min'))
    if (probably.rna) {
      ## A gene with non-protein coding transcripts is tagged as `utr` across
      ## across all of its exons in this scenario
      ranges <- Reduce(reducef, lapply(transcripts(x), ranges))
      exon.anno <- c(rep('utr', length(ranges)))
    } else {
      if (is.integer(collapse)) {
        ## Perhaps the Nth transcript doesn't have all of the expected segments
        ## in which case, we want to return an empty IRanges
        takeOrEmpty <- function(x, idx) {
          if (length(x) < idx) IRanges() else x[[idx]]
        }
        .exons <- list(utr5=takeOrEmpty(utr5(x), collapse),
                       cds=takeOrEmpty(cds(x), collapse),
                       utr3=takeOrEmpty(utr3(x), collapse))
      } else {
        ## `Reduce`-ing over an empty list (which can result from a gene not
        ## having a .utr3) will produce a NULL in insteand of an empty IRanges
        ## object, which will throw errors in the setdiff operations below
        null2empty <- function(x) if (is.null(x)) IRanges() else x
        .cds <- null2empty(Reduce(reducef, lapply(cds(x), ranges)))
        .utr3 <- null2empty(Reduce(reducef, lapply(utr3(x), ranges)))
        .utr5 <- null2empty(Reduce(reducef, lapply(utr5(x), ranges)))

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

    if (probably.rna && by == 'cds') {
      gr <- GRanges()
      values(gr) <- DataFrame(exon.anno=character())
    } else if (probably.rna || by == 'all') {
      gr <- GRanges(seqnames=rep(chromosome(x), length(ranges)),
                    ranges=ranges,
                    strand=rep(g.strand, length(ranges)))
      values(gr) <- DataFrame(exon.anno=exon.anno)
      gr <- gr[order(start(gr))]
    } else {
      gr <- lapply(.exons, function(e) {
        g <- GRanges(seqnames=rep(chromosome(x), length(e)),
                     ranges=e,
                     strand=rep(g.strand, length(e)))
        values(g) <- DataFrame(exon.anno=rep(by, length(g)))
        if (length(g) > 0) {
          g[order(start(g))]
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
            utr.ext <- GRanges(seqnames=chromosome(x), strand=g.strand,
                               ranges=IRanges(start=ext.start, end=ext.end))
            values(utr.ext) <- DataFrame(exon.anno='utr5*')
            gr <- c(utr.ext, gr)
          } else {
            ext.start <- end(gr[take[length(take)]]) + 1L
            ext.end <- ext.start + flank.up - 1L
            utr.ext <- GRanges(seqnames=chromosome(x), strand=g.strand,
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
            utr.ext <- GRanges(seqnames=chromosome(x), strand=g.strand,
                               ranges=IRanges(start=ext.start, end=ext.end))
            values(utr.ext) <- DataFrame(exon.anno='utr3*')
            gr <- c(gr, utr.ext)
          } else {
            ext.end <- start(gr[take[1L]]) - 1L
            ext.start <- ext.end - flank.down + 1
            utr.ext <- GRanges(seqnames=chromosome(x), strand=g.strand,
                               ranges=IRanges(start=ext.start, end=ext.end))
            values(utr.ext) <- DataFrame(exon.anno='utr3*')
            gr <- c(utr.ext, gr)
          }
        }
      }
    }
    
    gr
  })
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


