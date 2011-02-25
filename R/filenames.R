## Centralized place to keep filename logic. When writing to or reading from
## a file on disk, call appropriate methods here instead of hard wiring the
## names of files in the code
##
## Functions should be vectorized to return results for several file queries
## at a time.

##' @param gc.path The path to the cache directory for the genomic cache
annotatedChromosomeFN <- function(gc.path, seqname, gene.collapse='cover',
                                  cds.cover='min', flank.up=1000L,
                                  flank.down=flank.up, stranded=TRUE, ...) {
  if (inherits(gc.path, 'GenomicCache')) {
    gc.path <- cacheDir(gc.path)
  }
  message("cds.cover parameter not used yet.")
  stranded <- if (stranded) 'stranded' else 'not-stranded'
  fn <- sprintf("%s.annotated.collapse-%s.up-%d.down-%d.%s.rda", seqname,
                gene.collapse, flank.up, flank.down, stranded)
  file.path(gc.path, 'annotated.chromosomes', fn)
}

internalPrimingCacheFN <- function(gc.path, seqname, window, strand,
                                   type=c('wig', 'rda', 'ff')) {
  type <- match.arg(type)
  if (inherits(gc.path, 'GenomicCache')) {
    gc.path <- cacheDir(gc.path)
  }
  strand <- if (strand == '-' || strand == 'rev') 'rev' else 'fwd'
  fp <- file.path(gc.path, 'internal.priming', window)

  suffix <- switch(type, wig='bigWig', rda='rda', ff='ffdata')
  file.path(fp, type, paste(seqname, strand, suffix, sep="."))
}
