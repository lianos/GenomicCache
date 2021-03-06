## TODO: Building flank annotation is still screwed up! Look at RefSeq EIF4A1.
## there are many INTERNAL utr3* locations!
## setAs("GRanges", "AnnotatedChromosome", function(from) {
##   class(from) <- "AnnotatedChromosome"
##   from
## })

#' Checks a GenomicRanges-like object if it looks like an annotated genome
isValidAnnotatedGenome <- function(x, check.extensions=FALSE,
                                   key.by=c('seqnames','strand','entrez.id')) {
  if (!inherits(x, "GenomicRanges")) {
    return("This does not inherit from GenomicRanges")
  }
  errs <- character()
  ##############################################################################
  ## Ensure object is densely annotated and has no overlaps
  x.gaps <- gaps(x)
  x.gaps <- x.gaps[strand(x.gaps) != '*']
  has.gaps <- length(x.gaps) > 0
  has.overlaps <- length(findOverlaps(x, ignoreSelf=TRUE)) > 0
  if (has.gaps) {
    errs <- "Object has gaps in the annotation"
  }
  if (has.overlaps) {
    errs <- c(errs, "Object has overlapping annotations")
  }

  ###############################################################################
  ## Check that there is a max of 1 utr3* and utr5* for entrez.id
  if (check.extensions) {
    ag.dt <- as(anno, 'data.table')
    setkeyv(ag.dt, key.by)
    n.ext <- ag.dt[, {
      list(utr5e=sum(exon.anno == 'utr5*'), utr3e=sum(exon.anno == 'utr3*'))
    }, by=key.by]
    if (any(n.ext$utr5e > 1)) {
      errs <- c(errs, "utr5* annotations are not consistent")
    }
    if (any(n.ext$utr3e > 1)) {
      errs <- c(errs, "utr3* annotations are not consistent")
    }
  }

  if (length(errs) > 0) {
    return(errs)
  } else {
    return(TRUE)
  }
}


annotateIntronUtr3 <- function(ag, si.object=ag, nuke.factors=TRUE) {
  if (inherits(try(seqinfo(si.object)), "try-error")) {
    stop("An object with `seqinfo` is required for `si.object`")
  }

  agdt <- as.data.table(as.data.frame(ag))
  agdt[, exon.anno := as.character(exon.anno)]

  ## Hammer out any factors to characters -- you will thank me later
  if (nuke.factors) {
    for (wut in names(agdt)) {
      if (is.factor(agdt[[wut]])) {
        agdt[, (wut) := as.character(agdt[[wut]])]
      }
    }
  }

  is.genic <- !is.na(agdt$entrez.id)
  other <- agdt[!is.genic]
  genic <- agdt[is.genic]
  setkeyv(genic, c('seqnames', 'strand', 'entrez.id', 'start'))

  re <- genic[, {
    is.utr3 <- which(exon.anno == 'utr3')
    breaks.utr3 <- which(exon.anno %in% c('cds', 'utr', 'utr5'))
    intron.utr3 <- integer()
    if (strand == '+') {
      if (length(is.utr3) > 0 & length(breaks.utr3) > 0) {
        intron.utr3 <- is.utr3[is.utr3 < max(breaks.utr3)]
      }
    } else {
      if (length(is.utr3) > 0 & length(breaks.utr3) > 0) {
        intron.utr3 <- is.utr3[is.utr3 > min(breaks.utr3)]
      }
    }

    if (length(intron.utr3) > 0) {
      sd <- copy(.SD)
      sd$exon.anno[intron.utr3] <- 'intron.utr3'
    } else {
      sd <- .SD
    }

    sd
  }, by=c('seqnames', 'strand', 'entrez.id')]
  setcolorder(re, names(other))
  out <- rbind(other, re)

  gr <- as(out, "GRanges")
  gr <- rematchSeqinfo(gr, si.object)
  gr <- gr[order(gr)]

  gr
}

setValidity("AnnotatedChromosome", function(object) {
  errs <- character()
  expected.meta <- c('exon.anno', 'symbol', 'entrez.id')
  meta <- values(object)
  is.missing <- !expected.meta %in% colnames(meta)
  if (any(is.missing)) {
    errs <- sprintf("Annotation columns missing: %s",
                    paste(expected.meta[is.missing], collapse=","))
  }

  if (length(errs) == 0L) {
    TRUE
  } else {
    errs
  }

})

##' Plot the number of basepairs assigned to different annotations across genome
genomeAnnotationDistribution <- function(agenome) {
  stopifnot(inherits(agenome, "GenomicRanges"))

  agdt <- as.data.frame(agenome)[, c('width', 'exon.anno')]
  agdt <- data.table(agdt, key="exon.anno")
  counts <- agdt[, list(count=sum(width)), by="exon.anno"]

  g <- ggplot(as.data.frame(counts), aes(exon.anno, count)) + theme_bw() +
    geom_bar(aes(fill=exon.anno), stat='identity') +
      ylab("Count") + xlab("Exon Annotation") +
        opts(axis.text.x=theme_text(angle=-45, hjust=0, vjust=1),
             title=title)

}

##' Returns the annotated chromosome object from the given parameters.
##'
##' If the file has not been built, an error will be thrown prompting the
##' caller to genereate this file first.
##'
##' NOTE: 2010-10-23 -- Integrating entrez.id into the values() of annotated
##' chromosomes. We should switch to them as the primary key instead of using
##' the gene symbol. All entrez.id in here stuff hasn't been tested yet.
##'
##' @export
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param gcache \code{\linkS4class{GenomicCache}} object
##' @param seqnames A character vector of seqnames/chromosomes to specifying
##' which annotated chromosomes to get, or a \code{\linkS4class{GRanges}} object
##' which the unique seqnames are pulled out of
##' @param flank.up The flanking paremeters to specify which annotated chromosome
##' to retrieve
##' @param flank.down Same as \code{flank.up}
##' @param stranded Logical indicating if the stranded annotated chromosome is
##' desired
##'
##' @return An \code{\linkS4class{AnnotatedChromosome}} object
getAnnotatedChromosome <- function(gcache, seqnames, gene.collapse='longest',
                                   flank.up=1000L, flank.down=flank.up,
                                   stranded=TRUE) {
  if (inherits(seqnames, 'GRanges')) {
    seqnames <- as.character(seqnames(seqnames))
  }
  gene.collapse <- matchGFGeneCollapse(gene.collapse)
  seqnames <- unique(as.character(seqnames))
  annotated <- lapply(seqnames, function(seqname) {
    fn <- annotatedChromosomeFN(gcache, seqname, gene.collapse, flank.up=flank.up,
                                flank.down=flank.down, stranded=stranded)
    if (!file.exists(fn)) {
      do.try <- paste('gcache, collapse=%s, flank.up=%d, flank.down=%d,',
                      'stranded=%s, chrs=%s')
      do.try <- sprintf(do.try, gene.collapse, flank.up, flank.down, stranded,
                        seqname)
      stop(basename(fn), " file not found. Generate it first via:\n",
           sprintf('  annotateChromosomeByGenes(%s, ...)', do.try))
    }
    var.name <- load(fn)
    anno <- get(var.name, inherits=FALSE)
    ## as(anno, 'AnnotatedChromosome')
    class(anno) <- "AnnotatedChromosome"
    anno
  })
  suppressWarnings(do.call(c, unname(annotated)))
}

getAnnotatedGenome <- function(gcache, gene.collapse='longest', flank.up=1000L,
                               flank.down=flank.up, stranded=TRUE) {
  gene.collapse <- matchGFGeneCollapse(gene.collapse)
  genome.file <- annotatedGenomeFN(gcache, gene.collapse, flank.up=flank.up,
                                   flank.down=flank.down, stranded=stranded)
  if (file.exists(genome.file)) {
    return(load.it(genome.file))
  }

  ## The genome file isn't there, try loading it one chromosome at a time
  annotated <- lapply(seqnames(gcache), function(seqname) {
    fn <- annotatedChromosomeFN(gcache, seqname, gene.collapse, flank.up=flank.up,
                                flank.down=flank.down, stranded=stranded)
    if (file.exists(fn)) load.it(fn) else NULL
  })
  annotated <- annotated[!sapply(annotated, is.null)]
  if (length(annotated) == 0L) {
    stop("No annotated chromoomes found")
  }
  suppressWarnings(annotated <- do.call(c, unname(annotated)))
  as(annotated, 'AnnotatedChromosome')
}


##' The crank that turns the annotateChromosome function over the chromosomes
##' of a \code{GenomicCache}
##'
##' This function depends on the call to \code{\link{getGenesOnChromosome}}
##' function to load a cached list of gene objects, which have an "idealized"
##' version of the gene we are using, otherwise the speed will be terrible.
##'
##' @export
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##' @seealso \code{\link{getGenesOnChromosome}}
##' @seealso \code{\link{generateGFXGeneModels}}
##'
##' @param gcache A \code{\link{GenomicCache}} object to use as a reference
##' for the chromosome annotations
##' @param flank.up Number of base pairs upstream to extend the 5' UTR
##' @param flank.down Number of base pairs downstream to extend 3' UTR
##' @param gene.by The \code{by} parameter for the \code{GFGene::idealized}
##' function
##' @param gene.collapse The \code{collapse} parameter for the
##' \code{\link{idealized}} function
##' @param gene.collapse The \code{collapse} parameter for the
##' \code{\link{idealized}} function
##' @param fusion.filter Annotated fusion transcripts trip us up because the
##' fusion product will mask the two individual product and most all
##' exonic regions are identified as 'overlap'. This parameter uses a regex
##' to identify genes as fusions. Genes that match will be exlcuded from
##' annotation building. Set to \code{NULL} to keep all.
##' @return Invisibly returns an \code{\linkS4class{AnnnotatedChromosome}}
##' object
generateAnnotatedChromosomesByGenes <-
  function(gcache, flank.up=1000L, flank.down=flank.up, stranded=TRUE,
           gene.by='all', gene.collapse='cover', gene.cds.cover='min',
           add.ensembl.id=TRUE, chrs=NULL, do.save=TRUE,
           fusion.filter="[[:alnum:]]-[^0-9]", bsg.seqlengths=NULL,
           return.anno=TRUE, ...) {
  verbose <- checkVerbose(...)
  if (is.null(bsg.seqlengths)) {
    bsg.seqlengths <- seqlengths(gcache)
  }
  checkOrCreateDirectory(cacheDir(gcache, 'annotated.chromosomes'), TRUE)
  if (is.null(return.anno)) {
    return.anno <- !do.save
  }
  if (is.null(chrs)) {
    ## chrs <- seqlevels(gcache)
    chr.files <- gsub('.rda', '', dir(cacheDir(gcache, 'gene.models')))
    chrs <- sapply(strsplit(chr.files, '.', fixed=TRUE), tail, 1)
  }

  if (length(chrs) == 0) {
    stop("No chromosomes found to build gene models from")
  }

  illegal.chr <- !chrs %in% names(bsg.seqlengths)
  if (any(illegal.chr)) {
    stop("Bad chromosome names: ", paste(chrs[illegal.chr], collapse=","))
  }

  ## annos <- foreach(chr=chrs,
  ##                  .inorder=FALSE, .options.multicore=list(preschedule=FALSE),
  ##                  .verbose=verbose) %dopar% {
  annos <- mclapply(chrs, function(chr) {
    cat(chr, "...\n")
    seqlength <- bsg.seqlengths[chr]
    ## .gc <- duplicate(gcache, pre.load=NULL)

    ## genes <- tryCatch(getGenesOnChromosome(.gc, chr), error=function(e) NULL)
    genes <- tryCatch(getGenesOnChromosome(gcache, chr), error=function(e) NULL)
    if (is.null(genes)) {
      return(NULL)
    }
    if (is.character(fusion.filter)) {
      symbol <- sapply(genes, symbol)
      axe <- grep(fusion.filter, symbol)
      if (length(axe) > 0) {
        cat("... removing", length(axe), "genes as detected fusions")
        genes <- genes[-axe]
      }
    }

    entrez.id <- sapply(genes, entrezId)

    cat("... (", chr, ") cleaning gene models\n", sep="")
    models <- lapply(genes, function(gene) {
      ## Do not include genes whose transcripts do not overlap at all
      ## or exist on a different chromosome
      if (.goodGene(gene, chr)) {
        gm <- idealized(gene, by=gene.by, collapse=gene.collapse,
                        cds.cover=gene.cds.cover, flank.up=0L,
                        flank.down=0L, which.chr=chr)
        metadata(gm) <- c(metadata(gm), list(entrez.id=entrezId(gene)))
        gm
      } else {
        NULL
      }
    })
    keep <- !sapply(models, is.null)
    models <- models[keep]
    entrez.id <- entrez.id[keep]

    if (length(models) > 0) {
      cat("... (", chr, ") annotating chromosome\n", sep="")
      st <- proc.time()['elapsed']
      chr.anno <- annotateChromosome(models, entrez.id, flank.up, flank.down,
                                     seqname=chr, seqlength=seqlength,
                                     stranded=stranded)

      if (gene.collapse %in% c('longest', 'shortest')) {
        xref.txname <- unname(sapply(models, function(x) metadata(x)$tx_name))
        xref.entrez <- unname(sapply(models, function(x) metadata(x)$entrez.id))
        xref <- data.frame(entrez.id=xref.entrez, txname=xref.txname)
        xrf <- match(values(chr.anno)$entrez.id, xref$entrez.id)
        values(chr.anno)$tx_name <- xref$txname[xrf]
      }

      if (add.ensembl.id) {
        meta <- insertEnsemblId(as.data.frame(values(chr.anno)),
                                genome(gcache))
        values(chr.anno) <- DataFrame(meta)
      }

      cat(proc.time()['elapsed'] - st, "seconds\n")

      if (do.save) {
        ## fn <- annotatedChromosomeFN(.gc, chr, gene.collapse, flank.up=flank.up,
        ##                             flank.down=flank.down, stranded=stranded)
        fn <- annotatedChromosomeFN(gcache, chr, gene.collapse, flank.up=flank.up,
                                    flank.down=flank.down, stranded=stranded)

        cat("... (", chr, ") Saving to", fn, "\n")
        save(chr.anno, file=fn)
      }
    } else {
      cat("...(", chr, ") No genes found ... skipping\n")
      chr.anno <- NULL
    }

    ## dispose(.gc)
    if (return.anno) chr.anno else chr
  }, mc.preschedule=FALSE)

  annos <- annos[!sapply(annos, is.null)]
  annos <- suppressWarnings(do.call(c, unname(annos)))

  reanno <- tryCatch({
    rematchSeqinfo(annos, gcache)
  }, error=function(e) NULL)
  if (is.null(reanno)) {
    warning("Can't load BSgenome to reorder seqinfo")
  } else {
    annos <- reanno
    annos <- annos[order(annos)]
  }
  anno.fn <- annotatedGenomeFN(gcache, gene.collapse, gene.cds.cover, flank.up,
                               flank.down, stranded)
  save(annos, file=anno.fn)
  invisible(annos)
}

insertEnsemblId <- function(xdf, genome.id) {
  stopifnot(inherits(xdf, 'data.frame'))
  entrez.col <- which(colnames(xdf) == 'entrez.id')
  if (length(entrez.col) != 1L) {
    stop("entrez.id required")
  }
  ens.id <- matchEnsId(xdf[[entrez.col]], genome.id)
  ans <- transform(xdf[, 1:entrez.col, drop=FALSE], ens.id=ens.id)
  if (ncol(ans) < ncol(xdf) + 1L) {
    ans <- cbind(ans, xdf[,-(1:entrez.col),drop=FALSE])
  }
  for (cname in colnames(ans)) {
    if (is.factor(ans[[cname]])) {
      ans[[cname]] <- as.character(ans[[cname]])
    }
  }
  ans
}

matchEnsId <- function(x, genome.id) {
  x <- as.character(x)
  xentrez <- unique(x)
  xentrez <- xentrez[!is.na(xentrez)]

  xens <- getGeneIdFromEntrezId(genome.id, xentrez, 'ensGene', rm.unknown=FALSE)
  xref <- sapply(xens, '[[', 1L)
  names(xref) <- names(xens)

  xref[x]
}


##' Calculates a GRanges object for a chromosome, with internal ranges
##' corresponding to annotated exon boundaries for genes.
##'
##' NOTE: There is ambiguity to which gene an "intronic" locus belong to if
##' it is intronic to two genes that overlap the same genomic locus.
##' Investigate further!
##'
##' @export
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param gene.list list of \code{\linkS4class{GRanges}} objects, each object
##' in the list indicates the set of exons (across all isoforms) for a gene --
##' such as you might get from \code{\link{idealized}}. This can also be a
##' \code{GRangesList}, but a normal list is preffered for now since GRangesList
##' object are slow to iterate over.
##' @param entrez.id A vector of entrez id's that correspond to the genes in
##' gene.list
##' @param flank.up The number of basepairs to extend the 5'utr annotation
##' @param flank.down The number of basepairs to extend the 3'utr annotation
##' @param seqname The name of the chromosome we are building annotations for
##' @param seqlength The length of the chromosome
##'
##' @return An \code{\linkS4class{AnnotatedChromosome}} object.
annotateChromosome <- function(gene.list, entrez.id, flank.up=0L,
                               flank.down=flank.up, seqname=NULL,
                               seqlength=NA_integer_, stranded=TRUE) {
  if (length(entrez.id) != length(gene.list)) {
    stop("Must have entrezIds for all genes in gene.list")
  }
  if (!stranded) {
    stop("Unstranded annotation not fully functional ... look to fix code ",
         "from 'buildFlank...' onwards ...")
  }

  ## Parameter Bureaucracy
  if (is.null(seqname)) {
    seqname <- as.character(seqnames(gene.list[[1]][1]))
  }
  if (is.null(seqlength) || is.na(seqlength)) {
    names(seqlength) <- seqname
  }

  if (is(gene.list, 'GRangesList')) {
    exons <- unlist(unname(gene.list))
  } else {
    if (!all(sapply(gene.list, inherits, 'GRanges'))) {
      stop("Illegal object passed into gene.list")
    }
    exons <- do.call(c, unname(gene.list))
  }

  ## I am handling strand issues outside of the GRanges framework
  if (!stranded) {
    strand(exons) <- '*'
  }
  exons <- split(exons, strand(exons))

  ## Locate intervals that are "excluively" annotated, and others that
  ## have two+ annotations on the same region.
  interval.annos <- lapply(exons, function(.exons) {
    lapply(splitRangesByOverlap(ranges(.exons)), IntervalTree)
  })

  ## Use the "exclusive intervals" to redefine the exclusive portions of exons
  ## in each gene model by using the exclusive intervals that overlap with
  ## each genes exon.
  clean.list <- lapply(1:length(gene.list), function(idx) {
    g.exons <- gene.list[[idx]]
    if (length(g.exons) == 0L) {
      cat("idealized gene hosed for entrez", entrez.id[idx], "\n")
      return(GRanges())
    }
    ref.strand <- if (!stranded) '*' else as.character(strand(g.exons)[1])
    ## ref.strand <- tryCatch({
    ##  if (!stranded) '*' else as.character(strand(g.exons)[1])
    ## }, error=function(e) browser())
    itree <- interval.annos[[ref.strand]]$exclusive
    mm <- as.matrix(findOverlaps(ranges(g.exons), itree))
    if (nrow(mm) > 0) {
      ## use exclusive ranges only for exon boundaries
      clean <- GRanges(seqnames=seqname, ranges=IRanges(itree[mm[, 2]]),
                       strand=ref.strand, seqlengths=seqlength)
      ## Take annotations from original exons
      values(clean) <- values(g.exons)[mm[, 1], , drop=FALSE]
      values(clean)$symbol <- names(gene.list)[idx]
      values(clean)$entrez.id <- entrez.id[idx]
    } else {
      clean <- GRanges()
    }
    clean
  })

  annotated <- do.call(c, unname(clean.list))

  ## Convert overlap regions to GRanges and combine
  overlaps <- lapply(names(interval.annos), function(istrand) {
    itree <- interval.annos[[istrand]]$overlap
    if (length(itree) > 0L) {
      GRanges(seqnames=seqname, ranges=IRanges(itree), strand=istrand)
    } else {
      GRanges()
    }
  })
  overlaps <- do.call(c, unname(overlaps))
  if (length(overlaps) > 0) {
    values(overlaps) <- DataFrame(exon.anno='overlap', symbol=NA, entrez.id=NA)
    annotated <- c(annotated, overlaps)
  }
  ## annotated <- annotated[order(ranges(annotated))]
  annotated <- sort(annotated)

  ## Annotated extended/flanking utrs. If the extended flank runs into
  ## a region that is already annotated, we only take the region that
  ## starts the flank up until the first annotation.
  if (flank.down > 0) {
    down <- buildFlankAnnotation(annotated, flank.down, 'down', seqlength)
    annotated <- c(annotated, down)
    resort <- TRUE
  }

  if (flank.up > 0) {
    up <- buildFlankAnnotation(annotated, flank.up, 'up', seqlength)
    annotated <- c(annotated, up)
    resort <- TRUE
  }

  if (resort) {
    ## annotated <- annotated[order(ranges(annotated))]
    annotated <- sort(annotated)
    resort <- FALSE
  }

  ## Annotate introns
  introns <- buildIntronAnnotation(annotated, stranded=stranded)
  annotated <- c(annotated, introns)

  ## Whatever isn't marked by now must be intergenic
  intergenic <- buildIntergenicRegions(annotated, stranded=stranded)
  annotated <- c(annotated, intergenic)
  ## annotated <- annotated[order(ranges(annotated))]
  annotated <- sort(annotated)

  if (length(findOverlaps(annotated, ignoreSelf=TRUE, type='any')) > 0) {
    warning("Annotation for chromosome", seqname, "is not clean", sep=" ")
  }

  as(annotated, 'AnnotatedChromosome')
}

##' Returns the exclusive portion of the ranges in .ranges as $exclusive.
##' The regions that ovlerap in .ranges are removed and put into $overlap
##'
##' @return An \code{\link{IRangesList}} object with $exsluive and $overlap
splitRangesByOverlap <- function(.ranges) {
  if (!inherits(.ranges, 'IRanges')) {
    stop("IRanges object required")
  }
  if (length(.ranges) == 0L) {
    exclusive <- IRanges()
    overlap <- IRanges()
  } else {
    disjoint.ranges <- disjoin(.ranges)
    disjoint.matches <- subjectHits(findOverlaps(.ranges, disjoint.ranges))
    ## A range is exclusive if it doesn't overlap more than one disjoint.range
    is.exclusive <- tabulate(disjoint.matches) == 1L
    exclusive <- disjoint.ranges[is.exclusive]
    overlap <- reduce(disjoint.ranges[!is.exclusive])
  }
  IRangesList(exclusive=exclusive, overlap=overlap)
}

##' @nord
buildIntergenicRegions <- function(annotated, stranded=TRUE) {
  intergenic <- gaps(annotated)
  take <- as.logical(seqnames(intergenic) == seqnames(annotated)[1])
  intergenic <- intergenic[take]
  if (stranded) {
    intergenic <- intergenic[strand(intergenic) != '*']
  }
  values(intergenic) <- DataFrame(exon.anno='intergenic', symbol=NA,
                                  entrez.id=NA)
  intergenic
}

##' @nord
buildIntronAnnotation <- function(annotated, stranded=TRUE) {
  bounds <- annotatedTxBounds(annotated)
  unannotated <- gaps(annotated)
  if (stranded) {
    unannotated <- unannotated[strand(unannotated) != '*']
  }
  o <- findOverlaps(unannotated, bounds)
  mm <- as.matrix(o)

  if (nrow(mm) > 0) {
    ## redundant matches can happen -- we ignore them for now (danger!)
    ## and just pick the first
    m2 <- mm[!duplicated(mm[, 1]),,drop=FALSE]
    introns <- unannotated[m2[,1]]
    values(introns) <- DataFrame(exon.anno='intron',
                                 symbol=values(bounds)$symbol[m2[, 2]],
                                 entrez.id=values(bounds)$entrez.id[m2[, 2]])
  } else {
    introns <- GRanges()
  }
  introns
}

##' @nord
buildFlankAnnotation <- function(annotated, distance, direction, seqlength=NA) {
  direction <- match.arg(direction, c('up', 'down'))
  resize.fix <- if (direction == 'up') 'start' else 'end'
  flank.start <- direction == 'up'
  exon.anno <- if (direction == 'up') 'utr5*' else 'utr3*'
  ## DEBUG: 2012-01-08, utr5* region for APITD1-CORT (and rest of annotations!),
  ##         - chr1 [10490626, 10490803] +) seems "wrong" and is due to upstream
  ##           annotations being marked as overlap.
  ##         -  why is chr1 10507873 10508775 (+) intergenic!?
  ## DEBUG: 2012-01-08: Why is NPPA not showing up at all? Look at the CLCN6,
  ##        NPPA-AS1 and NPPA locus.
  bounds <- annotatedTxBounds(annotated)
  chr.mask <- reduce(union(annotated, bounds))
  flanks <- flank(bounds, width=distance, start=flank.start)
  ## unique.flanks <- setdiff(flanks, annotated)
  ## unique.flanks <- setdiff(flanks, bounds)
  unique.flanks <- setdiff(flanks, chr.mask)
  unique.flanks <- resize(unique.flanks, width=width(unique.flanks) + 1,
                          fix=resize.fix)
  o <- findOverlaps(unique.flanks, bounds)
  mm <- as.matrix(o)

  if (nrow(mm) > 0) {
    ## There will be duplicate matches here due to txbounds that overlap
    ## with eachother (think genes inside of other genes). Keep only flanks
    ## that overlap with one annotated tx bound, the others are tossed.
    ## keep <- tabulate(mm[, 1]) == 1L
    keep <- !(duplicated(mm[, 1]) | duplicated(mm[, 1], fromLast=TRUE))
    mm <- mm[keep, , drop=FALSE]
    new.flanks <- unique.flanks[mm[, 1]]
    new.flanks <- resize(new.flanks, width=width(new.flanks)-1, fix=resize.fix)
    values(new.flanks) <- values(bounds)[mm[, 2], , drop=FALSE]
    values(new.flanks)$exon.anno <- exon.anno
  } else {
    new.flanks <- GRanges()
  }

  new.flanks <- trimRangesToSeqlength(new.flanks, seqlength)
  new.flanks
}

##' @nord
trimRangesToSeqlength <- function(granges, seqlength=NA) {
  ## TODO: Handle isCircular=TRUE chromosomes
  too.low <- start(granges) < 1
  if (any(too.low)) {
    granges[too.low] <- 1L
  }

  if (!is.na(seqlength)) {
    too.high <- end(granges) > seqlength
    if (any(too.high)) {
      end(granges[too.high]) <- seqlength
    }
  }

  granges
}

##' Determine the transcription bounds of annotated regions (by entrez.id).
##'
##' Calculates the extremes of regions annotated for a given entrez.id.
##'
##' This function assumes that all ranges in \code{annotated} are on the
##' same chromosome.
##'
##' @param annotated A data.frame/data.table/AnnotatedChromsome
##' @param flank.up Number of bases to extend the upstream tx bound
##' @param flank.down Number of bases to extedn the downstream tx bound
##' @return A \code{GRanges} with the transcription bounds
annotatedTxBounds <- function(annotated, flank.up=0L, flank.down=0L,
                              seqlength=NA) {
  if ((inherits(annotated, 'GenomicRanges') && length(annotated) == 0L) ||
      (inherits(annotated, 'data.frame') && nrow(annotated) == 0L)) {
    stop("No annotations in `annotated`")
  }

  ## Calculate inferredmax-bounds by symbol
  if (is.na(seqlength) && is(annotated, 'GRanges')) {
    seqlength <- seqlengths(annotated)
  }

  dt <- subset(as(annotated, 'data.table'), !is.na(entrez.id))
  setkeyv(dt, c('seqnames', 'strand', 'entrez.id'))

  ## by entrez.id, chr, strand
  bounds <- dt[, {
    transform(.SD[1L], start=min(.SD$start), end=max(.SD$end), exon.anno='txbound')
  }, by=key(dt)]

  bounds <- as(bounds, 'GRanges')
  if (flank.up > 0) {
    bounds <- resize(bounds, width=width(bounds) + flank.up, fix='end')
  }
  if (flank.down > 0) {
    bounds <- resize(bounds, width=width(bounds) + flank.down, fix='start')
  }
  if (length(unique(as.character(seqnames(bounds)))) == 1) {
    bounds <- bounds[order(ranges(bounds))]
  }

  if (inherits(annotated, 'GRanges')) {
    old.meta <- values(annotated)
    new.meta <- resortColumns(values(bounds), values(annotated))

    for (col in colnames(old.meta)) {
      if (is.character(old.meta[[col]])) {
        new.meta[[col]] <- as.character(new.meta[[col]])
      }
    }
    values(bounds) <- new.meta
  }

  if (flank.up + flank.down > 0 && !is.na(seqlength)) {
    bounds <- trimRangesToSeqlength(bounds, seqlength)
  }

  if (!all(is.na(seqlength))) {
    seqlengths(bounds) <- seqlength
  }

  bounds
}

##' @nord
resortColumns <- function(from, to) {
  common.names <- intersect(colnames(from), colnames(to))
  if (length(common.names) != length(colnames(to))) {
    stop("Need matching column names")
  }

  xref <- match(colnames(to), colnames(from))
  from[, xref]
}

##' Flagging genes as "bad" if:
##' (i)   they map to more than one chromosome
##'       (my cached idealized version is hosed when this happens)
##' (ii)  the intersection of all transcript boundaries is empty
##'       (This hoses the utr{3|5}* annotation logic)
##'
##' TODO: Fix (ii) so that genes like DNMT3A don't get hosed. Its
##'       intersection turns up empty.
##'
##'     :::||||||||||----------------//---------------||||||||:::::
##'
##'         :::||||---|||:::                 :::||||----||||::::
##' @nord
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
