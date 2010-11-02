## setOldClass(c("data.table", "data.frame"))
setOldClass(c('data.table', 'data.frame'))

setAs("GRanges", "data.table", function(from) {
  if (length(from) == 0L) {
    return(data.table())
  }
  data.table(as.data.frame(from))
})

setAs("IRanges", "data.table", function(from) {
  if (length(from) == 0L) {
    return(data.table())
  }
  data.table(cbind(as.data.frame(from), as.data.frame(values(from))))
})

setAs("data.table", "GRanges", function(from) {
  if (nrow(from) == 0L || all(is.na(from))) {
    return(GRanges())
  }
  if (!all(c('seqnames', 'start', 'end') %in% colnames(from))) {
    stop("seqnames, start, end required")
  }
  .strand <- if ('strand' %in% colnames(from)) from$strand else '*'
  gr <- GRanges(seqnames=from$seqnames, ranges=IRanges(from$start, from$end),
                strand=.strand)
  gr.colnames <- c('seqnames', 'start', 'end', 'strand', 'width')
  meta.cols <- colnames(from)[!colnames(from) %in% gr.colnames]
  if (length(meta.cols) > 0) {
    values(gr) <- DataFrame(from[, meta.cols, with=FALSE])
  }
  gr
})

setAs("data.table", "IRanges", function(from) {
  if (nrow(from) == 0L || all(is.na(from))) {  
    return(IRanges())
  }
  if (!all(c('start', 'end') %in% colnames(from))) {
    stop("seqnames, start, end required")
  }
  ir <- IRanges(from$start, from$end)
  ir.colnames <- c('start', 'end', 'width')
  meta.cols <- colnames(from)[!colnames(from) %in% ir.colnames]
  if (length(meta.cols) > 0) {
    values(ir) <- DataFrame(from[, meta.cols, with=FALSE])
  }
  ir
})
