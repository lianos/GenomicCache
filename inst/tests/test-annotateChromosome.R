################################################################################
## Test data
clengths <- c(chr1=10000, chr2=5000)
chr1 <- factor('chr1', levels=c('chr1', 'chr2'))
chr2 <- factor('chr2', levels=c('chr1', 'chr2'))
g1 <- GRanges(seqname=chr1,
              ranges=IRanges(c(10, 40, 80, 120), width=21),
              strand='+',
              seqlengths=clengths)
values(g1) <- DataFrame(exon.anno=c('utr5', 'cds', 'cds', 'utr3'))

g2 <- GRanges(seqname=chr1,
              ranges=IRanges(c(200, 250, 300, 400), width=31),
              strand='+',
              seqlengths=clengths)
values(g2) <- DataFrame(exon.anno=c('utr5', 'cds', 'cds', 'utr3'))

## First exon overlaps with parts of exons 2,3 in g1
## Secon exon is in intron 3
o1 <- GRanges(seqname=chr1,
              ranges=IRanges(start=c(50, 105), end=c(85, 115)),
              strand='+',
              seqlengths=clengths)
values(o1) <- DataFrame(exon.anno=c('cds', 'cds'))

g3r <- GRanges(seqnames=chr1,
               ranges=IRanges(start=c(15, 50), end=c(45, 100)),
               strand='-',
               seqlengths=clengths)
values(g3r) <- DataFrame(exon.anno=c('cds', 'utr5'))

gchr1 <- GRangesList(GeneA=g1, GeneB=g2, GeneC=o1, GeneD=g3r)
lchr1 <- list(GeneA=g1, GeneB=g2, GeneC=o1, GeneD=g3r)

################################################################################
## Run tests
context("Annotate Chromosome")
test_that("Overlapping regions are disjoint(ed)", {
  
})

test_that("Chromosome is annotated", {
  
})
