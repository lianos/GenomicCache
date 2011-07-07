## Most (all) objects defined in this package have `.cache` slots, which are
## environments that are used to (surprise )store/cache things

##' Checks for the existense of a variable `what` in the cache, if it doesn't
##' exist, it computes the appropriate value by evaluating \code{expr}, and
##' sets its value into the \code{cache} to variable named \code{what}. The
##' value is then returned.
##'
##' @usage
##' cacheFetch(x, 'transcripts', {
##'   transcripts(x@.transcripts, by='tx)
##' })

setMethod("cacheFetch", c(x="ANY"),
function(x, what, expr, force.eval) {
  if (!exists(what, x@.cache, inherits=FALSE) ||
      (is.logical(force.eval) && force.eval)) {
    value <- eval.parent(expr)
    assign(what, value, x@.cache, inherits=FALSE)
  }
  get(what, x@.cache, inherits=FALSE)
})

setMethod("clearCache", c(x="ANY"),
function(x, ...) {
  verbose <- checkVerbose(...)
  clear <- ls(x@.cache)
  more <- list(...)
  
  if (length(more) > 0) {
    more <- unlist(more)
    more <- more[sapply(more, is.character)]
    clear <- more[more %in% clear]
  }

  if (length(clear) == 0L) {
    return(NULL)
  }
  
  if (verbose) {
    cat("  Clearing cache:", paste(clear, collapse=", "), "\n")
  }
  rm(list=clear, envir=x@.cache)
  invisible(gc(reset=TRUE))  
})

generateCacheName <- function(base, ...) {
  params <- list(...)
  if (length(params) > 0) {
    values <- sapply(params, function(p) {
      name <- if (is.null(p)) "NULL" else as.character(p)
      if (length(name) > 1) {
        name <- paste(name, collapse="|")
      }
      name
    })
  } else {
    values <- 'missing'
  }
  
  params <- paste(names(params), values, sep=":", collapse=",")
  paste(base, params, sep="/")
}
