##' Calculates a GRanges object for a chromosome, with internal ranges
##' corresponding to annotated exon boundaries for genes.
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
                               seqname='NA', seqlength=NA,stranded=TRUE) {
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
  
  ## is.rev <- strand(annotated) == '-'
  
  ## Annotated extended/flanking utrs. If the extended flank runs into
  ## a region that is already annotated, we only take the region that
  ## starts the flank up until the first annotation.
  ## bounds <- .txBounds(annotated)
  ## bounds.is.rev <- strand(bounds) == '-'
  
  ## if (flank.up > 0) {
  ##   flank.up <- flank(bounds, width=flank.up, start=TRUE)
    
  ##   o <- findOverlaps(annotated, flank.up)
  ##   problem.idx <- unique(subjectHits(o))
  ##   .flank.up <- flank.up[-problem.idx]
  ##   problem <- flank.up[problem.idx]
  ##   values(.flank.up)$exon.anno <- 'utr5*'

  ##   pieces <- setdiff(problem, annotated)
  ##   pieces <- pieces[order(start(pieces))]
    
  ##   o.fwd <- findOverlaps(ranges(problem[strand(problem) != '-']),
  ##                         ranges(pieces),
  ##                         select='last')
  ##   take.flank <- which(!is.na(o.fwd))
  ##   take.piece <- o.fwd[take.flank]

  ##   ## fwd
  ##   fu <- flank.up[!bounds.is.rev]
  ##   fu <- fu[order(start(fu))]
  ##   fanno <- annotated[!is.rev]
  ##   bf <- bounds[!bounds.is.rev]
  ##   d <- setdiff(ranges(fu), ranges(fanno))
  ##   o <- findOverlaps(d, ranges(bf), maxgap=1L)
  ##   mm <- matchMatrix(o)
  ## }
  ## if (flank.down > 0) {
  ##   flank.down <- flank(bounds, width=flank.down, start=FALSE)
  ## }
  
  ## ## annotate remaining introns intergenic regions
  ## unannotated <- gaps(cleaned)
  ## if (stranded) {
  ##   unannotated <- unannotated[!strand(unannotated) == '*']
  ## } else {
  ##   unannotated <- unannotated[strand(unannotated) == '*']
  ## }
  
  annotated
}

.txBounds <- function(annotated, flank.up=0L, flank.down=0L) {
  ## Calculate inferredmax-bounds by symbol
  dt <- data.table(start=start(annotated), end=end(annotated),
                   symbol=values(annotated)$symbol,
                   strand=as.character(strand(annotated)))
  key(dt) <- 'symbol'
  bounds <- dt[, {
    list(start=min(start), end=max(end), strand=strand[1])
  }, by='symbol']
  bounds <- bounds[-is.na(bounds$symbol),]
  ## calculate flanks
  flanks <- GRanges(seqnames=seqnames(annotated[1]),
                    ranges=IRanges(start=bounds$start, end=bounds$end),
                    strand=bounds$strand)
  values(flanks) <- DataFrame(exon.anno='flank', symbol=bounds$symbol)
  is.rev <- strand(flanks) == '-'
  if (flank.up > 0) {
    ranges(flanks) <- resize(ranges(flanks), width(flanks) + flank.up,
                             fix=ifelse(is.rev, 'start', 'end'))
  }
  if (flank.down > 0) {
    ranges(flanks) <- resize(ranges(flanks), width(flanks) + flank.down,
                             fix=ifelse(is.rev, 'end', 'start'))
  }
  flanks
}

## Flagging genes as "bad" if:
## (i)   they map to more than one chromosome
##       (my cached idealized version is hosed when this happens)
## (ii)  the intersection of all transcript boundaries is empty
##       (This hoses the utr{3|5}* annotation logic)
.goodGene <- function(gene) {
  xcripts <- transcripts(gene)

  ## Does it map to more than one chromosome?
  chrs <- unique(as.character(unique(seqnames(xcripts))))
  if (length(chrs) > 1) {
    return(FALSE)
  }

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
annotateGenomeByGenes <- function(gcache, flank.up=1000L, flank.down=flank.up,
                                  gene.by='all', gene.collapse='cover',
                                  gene.cds.cover='min', chrs=NULL,
                                  do.save=TRUE, path=cacheDir(gcache), ...) {
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
    genes <- getGenesOnChromosome(gcr, chr)

    cat("... cleaning gene models ...\n")
    models <- lapply(genes, function(gene) {
      ## Do not include genes whose transcripts do not overlap at all
      ## or exist on a different chromosome
      if (.goodGene(gene)) {
        gm <- idealized(gene, by=gene.by, collapse=gene.collapse,
                        cds.cover=gene.cds.cover, flank.up=flank.up,
                        flank.down=flank.down)
        ## remove the utr{3|5}*
        axe <- grep("*", values(gm)$exon.anno, fixed=TRUE)
        if (length(axe) > 0L) {
          gm <- gm[-axe]
        }
        gm
      } else {
        NULL
      }
    })
    models <- models[!sapply(models, is.null)]
    
    cat("... annotating chromosome ...")
    st <- proc.time()['elapsed']
    chr.anno <- annotateChromosome(models, chr, flank.up, flank.down,
                                   seqname=chr, seqlength=chr.length,
                                   stranded=TRUE)
    cat(proc.time()['elapsed'] - st, "secons\n")

    if (do.save) {
      fname <- sprintf("%s.annotated.stranded.rda", chr)
      save(chr.anno, file=file.path(path, fname))
    }

    chr.anno
  })

  annos
}
