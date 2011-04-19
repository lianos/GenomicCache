context("Neighbor Info")

test_that("GRagnes <-> GRanges neighbors (same strands)", {
  events <- GRanges(c(rep("chr1", 5), rep('chr2', 3)),
                    c(IRanges(c(1, 5, 10, 30, 20), width=2),
                      IRanges(c(10, 5, 50), width=2)),
                    strand='+')
  is.chr1 <- as.logical(seqnames(events) == 'chr1')
  annos <- GRanges(c('chr1', 'chr1', 'chr2', 'chr2'),
                   c(IRanges(c(15, 80), c(18, 82)),
                     IRanges(c(18, 100), c(20, 105))),
                   strand='+')
  n <- neighbors(events, annos, include.overlaps=TRUE)

  up.idx <- as.integer(c(NA, NA, NA, 1, 1,   NA, NA, 3))
  dn.idx <- as.integer(c(1, 1, 1, 2, 2,      3, 3, 4))

  expect_identical(n$upstream.idx, up.idx,
                   info="Upstream neighbors miscalculated")
  expect_identical(n$downstream.idx, dn.idx,
                   info="Downstream neighbors miscalculated")
})
