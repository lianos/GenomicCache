getInternalPrimingScores <- function(gcache, seqname, strand, window.size,
                                     type=c('rda', 'wig')) {
  type <- match.arg(type)
  file <- internalPrimingCacheFN(gcache, seqname, window.size, strand,
                                 type=type)
  if (type == 'wig'){
    BigWigFile(file)
  } else {
    load.it(file)
  }
}
