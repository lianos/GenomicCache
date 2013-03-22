##' Create an ExonCountSet for DEXSeq
##'
##' @param x An annotated genome object
##' @param input List of BAM files
##' @param A vector listing the conditions of each group, or
##' a data.frame for multivariate designs
createExonCountSet <- function(x, input, groups, ...) {
  stopifnot(inherits(x, "GenomicRanges"))
  if (is(input, 'BamFile')) {
    input <- list(input)
  }
  stopifnot(is.list(input) && all(sapply(input, is, "BamFile")))
  stopifnot(all(file.exists(sapply(input, path))))

  if (is.vector(groups)) {
    groups <- data.frame(condition=groups)
  }
  if (all(rownames(groups) == as.character(1:nrow(groups)))) {
    if (!is.null(names(input))) {
      rownames(groups) <- names(input)
    } else {
      rownames(groups) <- make.unique(groups$condition)
    }
  }
  stopifnot(is(groups, "data.frame"))
  for (wut in colnames(groups)) {
    groups[[wut]] <- factor(as.character(groups[[wut]]))
  }
  se <- tabulateIntoSummarizedExperiment(x, input)
}