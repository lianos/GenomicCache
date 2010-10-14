setClass('AnnotatedChromosome', contains="GRanges")

.annotatedChromosomeFileName <- function(gcache, seqname, flank.up, flank.down,
                                         stranded) {
  stranded <- if (stranded) 'stranded' else 'not-stranded'
  fn <- sprintf("%s.annotated.up-%d.down-%d.%s.rda", seqname, flank.up,
                flank.down, stranded)
  file.path(cacheDir(gcache), 'annotated.chromosomes', fn)
}

##' Returns the annotated chromosome object from the given parameters.
##'
##' If the file has not been built, an error will be thrown prompting the
##' caller to genereate this file first.
getAnnotatedChromosome <- function(gcache, seqname, flank.up=1000L,
                                   flank.down=1000L, stranded=TRUE) {
  fn <- .annotatedChromosomeFileName(gcache, seqname, flank.up, flank.down,
                                     stranded)
  if (!file.exists(fn)) {
    v <- sprintf('gcache, flank.up=%d, flank.down=%d, stranded=%s, chrs=%s)',
                 flank.up, flank.down, stranded, seqname)
    stop(basename(fn), " file not found. Generate it first via:\n",
         sprintf('  annotateChromosomeByGenes(%s, ...)', v))
  }
  var.name <- load(fn)
  anno <- get(var.name, inherits=FALSE)
  class(anno) <- 'AnnotatedChromosome'
  anno
}

##' Calculates a GRanges object for a chromosome, with internal ranges
##' corresponding to annotated exon boundaries for genes.
##'
##' NOTE: It's not clear how a genomic locus that belongs to an intron of
##' two genes is assigned the "gene owner". Investigate further!
##' 
##' @param gene.list list of GRanges objects, each object in the list indicates
##' the set of exons (across all isoforms) for a gene -- such as you might get
##' from \code{\link{idealized}(GFGene)}. This can also be a \code{GRangesList},
##' but a normal list is preffered for now since GRangesList object are slow
##' to iterate over.
##' @param flank.up The number of basepairs to extend the 5'utr annotation
##' @param flank.down The number of basepairs to extend the 3'utr annotation
##' @param seqname The name of the chromosome we are building annotations for
##' @param seqlength The length of the chromosome
annotateChromosome <- function(gene.list, flank.up=0L, flank.down=flank.up,
                               seqname='NA', seqlength=NA, stranded=TRUE) {
  if (!stranded) {
    stop("Unstranded annotation not fully functional ... look to fix code ",
         "from 'buildFlank...' onwards ...")
  }
  ## Parameter Bureaucracy
  sl <- integer()
  sl[[seqname]] <- seqlength
  seqlength <- sl
  
  if (is(gene.list, 'GRangesList')) {
    exons <- unlist(gene.list)
  } else {
    gl <- gene.list
    names(gl) <- NULL
    exons <- do.call(c, gl)
  }
  
  ## I am handling strand issues outside of the GRanges framework
  if (!stranded) {
    strand(exons) <- '*'
  }
  exons <- split(exons, strand(exons))
  
  ## Locate intervals that are "excluively" annotated, and others that
  ## have two+ annotations on the same region.
  interval.annos <- lapply(exons, function(.exons) {
    if (length(.exons) == 0L) {
      exclusive.itree <- IntervalTree(IRanges())
      overlap.itree <- IntervalTree(IRanges())
    } else {
      d.exons <- disjoin(.exons)
      d.matches <- subjectHits(findOverlaps(.exons, d.exons))
      d.tab <- tabulate(d.matches)
      exclusive <- d.tab == 1L
      exclusive.itree <- IntervalTree(ranges(d.exons[exclusive]))
      overlap.itree <- IntervalTree(reduce(ranges(d.exons[!exclusive])))
    }
    list(exclusive=exclusive.itree, overlap=overlap.itree)
  })

  ## Use the "exclusive intervals" to pull out the strictly exclusive exons
  ## per gene.
  clean.list <- lapply(1:length(gene.list), function(idx) {
    g.name <- names(gene.list)[idx]
    g.exons <- gene.list[[idx]]
    ref.strand <- if (!stranded) '*' else as.character(strand(g.exons)[1])
    itree <- interval.annos[[ref.strand]]$exclusive
    mm <- matchMatrix(findOverlaps(ranges(g.exons), itree))
    if (nrow(mm) > 0) {
      clean <- GRanges(seqnames=seqname,
                       ranges=IRanges(itree[mm[, 2]]),
                       strand=ref.strand,
                       seqlengths=seqlength)
      values(clean) <- values(g.exons)[mm[, 1],,drop=FALSE]
      values(clean)$symbol <- g.name
    } else {
      clean <- GRanges()
    }
    clean
  })

  annotated <- do.call(c, clean.list)
  
  ## Convert overlap regions to GRanges and combine
  overlaps <- lapply(names(interval.annos), function(istrand) {
    itree <- interval.annos[[istrand]]$overlap
    if (length(itree) > 0L) {
      GRanges(seqnames=seqname, ranges=IRanges(itree), strand=istrand)
    } else {
      GRanges()
    }
  })
  overlaps <- do.call(c, overlaps)
  values(overlaps) <- DataFrame(exon.anno='overlap', symbol=NA)
  
  annotated <- c(annotated, overlaps)
  annotated <- annotated[order(start(annotated))]

  ## Annotated extended/flanking utrs. If the extended flank runs into
  ## a region that is already annotated, we only take the region that
  ## starts the flank up until the first annotation.
  if (flank.up > 0) {
    up.fwd <- buildFlankAnnotation(annotated, flank.up, 'up', '+')
    up.rev <- buildFlankAnnotation(annotated, flank.up, 'up', '-')
    annotated <- c(annotated, up.fwd, up.rev)
    resort <- TRUE
  }
  
  if (flank.down > 0) {
    down.fwd <- buildFlankAnnotation(annotated, flank.down, 'down', '-')
    down.rev <- buildFlankAnnotation(annotated, flank.down, 'down', '+')
    annotated <- c(annotated, down.fwd, down.rev)
    resort <- TRUE
  }

  if (resort) {
    annotated <- annotated[order(start(annotated))]
    resort <- FALSE
  }
  
  ## Annotate introns
  introns <- buildIntronAnnotation(annotated, stranded=stranded)
  annotated <- c(annotated, introns)
  
  ## Whatever isn't marked by now must be intergenic
  intergenic <- buildIntergenicRegions(annotated, stranded=stranded)
  annotated <- c(annotated, intergenic)
  annotated <- annotated[order(start(annotated))]
  class(annotated) <- 'AnnotatedChromosome'
  annotated
}

## if (FALSE) {
##   dt <- data.table(symbol=values(annotated)$symbol,
##                    strand=as.character(strand(annotated)),
##                    start=start(annotated),
##                    end=end(annotated),
##                    exon.anno=values(annotated)$exon.anno)
##   key(dt) <- 'symbol'
## }

buildIntergenicRegions <- function(annotated, stranded=TRUE) {
  intergenic <- gaps(annotated)
  take <- as.logical(seqnames(intergenic) == seqnames(annotated)[1])
  intergenic <- intergenic[take]
  if (stranded) {
    intergenic <- intergenic[strand(intergenic) != '*']
  }
  values(intergenic) <- DataFrame(exon.anno='intergenic', symbol=NA)
  intergenic
}

buildIntronAnnotation <- function(annotated, stranded=TRUE) {
  bounds <- annotatedTxBounds(annotated)
  unannotated <- gaps(annotated)
  if (stranded) {
    unannotated <- unannotated[strand(unannotated) != '*']
  }
  o <- findOverlaps(unannotated, bounds)
  mm <- matchMatrix(o)

  if (nrow(mm) > 0) {
    ## redundant matches can happen -- we ignore them for now (danger!)
    ## and just pick the first
    m2 <- mm[!duplicated(mm[, 1]),,drop=FALSE]
    introns <- unannotated[m2[,1]]
    values(introns) <- DataFrame(exon.anno='intron',
                                 symbol=values(bounds)$symbol[m2[, 2]])
  } else {
    introns <- GRanges()
  }
  introns
}

.fstart.val <- function(direction, strand) {
  ((direction == 'up') && (strand %in% c('+', '*'))) ||
  ((direction == 'down' && (strand == '-')))
}

buildFlankAnnotation <- function(annotated, fdist, fdir, fstrand) {
  fdir <- match.arg(fdir, c('up', 'down'))
  fstrand <- match.arg(fstrand, levels(strand()))
  if (fstrand == '*') {
    warning("Building flanks not supported for strand: *")
    return(GRanges())
  }
  
  is.rev <- strand(annotated) == '-'
  atake <- if (fstrand == '-') is.rev else !is.rev
  annotated <- annotated[atake]
  seqname <- seqnames(annotated)[1]
  
  bounds <- annotatedTxBounds(annotated)
  
  .flank <- flank(bounds, width=fdist, start=.fstart.val(fdir, fstrand))

  d <- setdiff(ranges(.flank), ranges(annotated))
  o <- findOverlaps(d, ranges(bounds), maxgap=1L)
  mm <- matchMatrix(o)
  if (nrow(mm) > 0) {
    ## There will be duplicate matches here when there shouldn't be
    ## I think this is due to small (1) width annotations ... ignore
    ## these
    
    ## axe <- unique(mm[,2][duplicated(mm[,2])])
    ## axe <- mm[,2] %in% axe
    ## mm <- mm[!axe,,drop=FALSE]

    ## There should only be max of 1 assignment per new flank,
    ## uniqueness over txBound ids doesn't work
    axe <- unique(mm[, 1][duplicated(mm[, 1])])
    keep <- !(mm[, 1] %in% axe)
    mm <- mm[keep, , drop=FALSE]
    
    new.flanks <- GRanges(seqnames=seqname, ranges=d[mm[,1]], strand=fstrand)
    values(new.flanks) <- values(bounds)[mm[,2],]
    exon.anno <- if (fdir == 'up') 'utr5*' else 'utr3*'
    values(new.flanks)$exon.anno <- exon.anno
  } else {
    new.flanks <- GRanges()
  }

  new.flanks
}

annotatedTxBounds <- function(annotated, flank.up=0L, flank.down=0L) {
  ## Calculate inferredmax-bounds by symbol
  dt <- data.table(start=start(annotated), end=end(annotated),
                   symbol=values(annotated)$symbol,
                   strand=as.character(strand(annotated)))
  key(dt) <- 'symbol'
  bounds <- dt[, {
    list(start=min(start), end=max(end), strand=strand[1])
  }, by='symbol']
  bounds <- bounds[!is.na(bounds$symbol),]
  ## calculate flanks
  flanks <- GRanges(seqnames=seqnames(annotated[1]),
                    ranges=IRanges(start=bounds$start, end=bounds$end),
                    strand=bounds$strand)
  values(flanks) <- DataFrame(exon.anno='flank',
                              symbol=as.character(bounds$symbol))
  is.rev <- strand(flanks) == '-'
  if (flank.up > 0) {
    ranges(flanks) <- resize(ranges(flanks), width(flanks) + flank.up,
                             fix=ifelse(is.rev, 'start', 'end'))
  }
  if (flank.down > 0) {
    ranges(flanks) <- resize(ranges(flanks), width(flanks) + flank.down,
                             fix=ifelse(is.rev, 'end', 'start'))
  }

  flanks <- flanks[order(start(flanks))]
  flanks
}

## Flagging genes as "bad" if:
## (i)   they map to more than one chromosome
##       (my cached idealized version is hosed when this happens)
## (ii)  the intersection of all transcript boundaries is empty
##       (This hoses the utr{3|5}* annotation logic)
.goodGene <- function(gene, which.chr) {
  xcripts <- transcripts(gene, which.chr=which.chr)

  ## Does it map to more than one chromosome?
  ## chrs <- unique(as.character(unique(seqnames(xcripts))))
  ## if (length(chrs) > 1) {
  ##   return(FALSE)
  ## }

  ## Are the txBounds disjoint?
  tx.bounds <- lapply(xcripts, function(x) range(ranges(x)))
  if (length(Reduce(intersect, tx.bounds)) == 0) {
    return(FALSE)
  }
  
  TRUE
}


##' The crank that turns the annotateChromosome function over the chromosomes
##' of a \code{GenomicCache}
##'
##' I'm depending on the getGenesOnChromosom function to load a cached list
##' of gene objects, which have an "idealized" version of the gene we are
##' using, otherwise the speed will be terrible.
##' 
##' @param gcache A \code{\link{GenomicCache}} object
##' @param flank.up Number of base pairs upstream to extend the 5' UTR
##' @param flank.down Number of base pairs downstream to extend 3' UTR
##' @param gene.by The \code{by} parameter for the \code{GFGene::idealized}
##' function
##' @param gene.collapse The \code{collapse} parameter for the
##' \code{GFGene::idealized} function
##' @param gene.collapse The \code{collapse} parameter for the
##' \code{GFGene::idealized} function
annotateChromosomeByGenes <- function(gcache, flank.up=1000L, flank.down=flank.up,
                                      stranded=TRUE, gene.by='all',
                                      gene.collapse='cover',
                                      gene.cds.cover='min', chrs=NULL,
                                      do.save=TRUE, ...) {
  bsg <- getBsGenome(gcache)
  if (is.null(chrs)) {
    chrs <- chromosomes(gcache)
  }
  
  illegal.chr <- !chrs %in% names(seqlengths(bsg))
  if (any(illegal.chr)) {
    stop("Bad chromosome names: ", paste(chrs[illegal.chr], collapse=","))
  }

  annos <- lapply(chrs, function(chr) {
    cat(chr, "...\n")
    chr.length <- seqlengths(bsg)[chr]
    genes <- getGenesOnChromosome(gcache, chr)

    cat("... cleaning gene models ...\n")
    models <- lapply(genes, function(gene) {
      ## Do not include genes whose transcripts do not overlap at all
      ## or exist on a different chromosome
      if (.goodGene(gene, chr)) {
        gm <- idealized(gene, by=gene.by, collapse=gene.collapse,
                        cds.cover=gene.cds.cover, flank.up=0, flank.down=0,
                        which.chr=chr)
        ## remove the utr{3|5}*
        ## axe <- grep("*", values(gm)$exon.anno, fixed=TRUE)
        ## if (length(axe) > 0L) {
        ##   gm <- gm[-axe]
        ## }
        gm
      } else {
        NULL
      }
    })
    models <- models[!sapply(models, is.null)]
    
    cat("... annotating chromosome ...")
    st <- proc.time()['elapsed']
    chr.anno <- annotateChromosome(models, flank.up, flank.down,
                                   seqname=chr, seqlength=chr.length,
                                   stranded=stranded)
    cat(proc.time()['elapsed'] - st, "secons\n")

    if (do.save) {
      fn <- .annotatedChromosomeFileName(gcache, chr, flank.up, flank.down,
                                         stranded)
      save(chr.anno, file=fn)
    }

    chr.anno
  })

  annos
}
