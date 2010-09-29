##' Calculates a GRanges object for a chromosome, with internal ranges
##' corresponding to annotated exon boundaries for genes.
##' 
##' @param gene.list list of GRanges objects, each object in the list indicates
##' the set of exons (across all isoforms) for a gene -- such as you might get
##' from \code{\link{idealized}(GFGene)}. This can also be a \code{GRangesList},
##' but a normal list is preffered for now since GRangesList object are slow
##' to iterate over.
##' @param seqname The name of the chromosome we are building annotations for
##' @param seqlength The length of the chromosome
annotateChromosome <- function(gene.list, seqname, seqlength=NA,
                               stranded=TRUE) {
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

  cleaned <- do.call(c, clean.list)
  
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

  cleaned <- c(cleaned, overlaps)
  cleaned <- cleaned[order(start(cleaned))]
  cleaned
}



##' The crank that turns the annotateChromosome function over the chromosomes
##' of a \code{GenomicCache}
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
                                  gene.cds.cover='min', chromosomes=NULL,...) {
  stop("Not implemented ...")
  if (is.null(chromosomes)) {
    chromosomes <- seqnames(gcache)
  }
  
  annos <- foreach(chr=chromosomes) %dopar% {
    genes <- getGenesOnChromosome(gcache, chr)
    ideal <- llply(genes, function(gene) {
      ig <- idealized(gene, by=gene.by, collapse=gene.collapse,
                      cds.cover=gene.cds.cover)
                      
    }, .progress='text')
  }
  
  annos <- do.call(GRangesList(annos))
}

