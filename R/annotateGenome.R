##' Creates a GRangesList covering all of the genome that has annotations
##' for each interval along the chromosome
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