setMethod("subset", "GRanges", function(x, subset, select, drop=FALSE, ...) {
  if (missing(subset)) {
    i <- TRUE
  } else {
    i <- eval(substitute(subset), values(x), parent.frame())
    i <- try(as.logical(i), silent = TRUE)
    if (inherits(i, "try-error"))
      stop("'subset' must be coercible to logical")
    i <- i & !is.na(i)
  }
  if (missing(select))
    j <- TRUE
  else {
    nl <- as.list(seq_len(ncol(values(x))))
    names(nl) <- colnames(values(x))
    j <- eval(substitute(select), nl, parent.frame())
  }
  xx <- x[i]
  values(xx) <- values(xx)[, j, drop=drop]
  xx
})
