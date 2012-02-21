##' Adds `intron.utr3` category for 3'UTRs "inside" the longest isofrom.
##'
##' This is a combination of a region being defined as `utr3` with a
##' utr3.index != 1
annotateWithIntronUtr3 <- function(ag, ...) {
  verbose <- checkVerbose(...)
  dt <- data.table(as.data.frame(ag),
                   key=c('seqnames', 'strand', 'entrez.id', 'start'))
  dt$exon.anno <- as.character(dt$exon.anno)

  if (verbose) {
    cat("... annotating intron.utr3's\n")
  }

  ag.dt <- dt[, {
    SD <- copy(.SD)
    if (!is.na(entrez.id) && any(utr3.index > 1)) {
      change.utr3 <- utr3.index != max(utr3.index) & exon.anno == 'utr3'
      SD$exon.anno[change.utr3] <- 'intron.utr3'
    }
    SD
  }, by=head(key(dt), -1)]

  ag <- as(ag.dt, 'GRanges')
  values(ag) <- DataFrame(lapply(values(ag), function(x) {
    if (is.factor(x)) x <- as.character(x)
    x
  }))

  class(ag) <- 'AnnotatedChromosome'
  ag
}

