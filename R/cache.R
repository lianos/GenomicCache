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
cacheFetch <- function(x, what, expr) {
  if (!exists(what, x@.cache, inherits=FALSE)) {
    value <- eval.parent(expr)
    assign(what, value, x@.cache, inherits=FALSE)
  }
  get(what, x@.cache, inherits=FALSE)
}

setMethod("clearCache", c(x="GenomicCache"),
function(x, ...) {
  verbose <- checkVerbose(...)
  clear <- ls(x@.cache)
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
################################################################################
## These are used to store cached objects for a given GenomicCache, such as
## stored GFGene objects / chromosome, see getGenesOnChromosome


## Pass in extra sub directories through \code{...}
GFXCacheDir <- function(path=NULL, ...) {
  if (is.null(path)) {
    path <- Sys.getenv(.GFX$cache$environment.key)
    if (is.null(path)) {
      stop("No cache directory provided, and one is not set in the environment")
    }
    subdirs <- list(...)
    if (length(subdirs) > 0) {
      path <- do.call(file.path, c(list(path), subdirs))
    }
  }
  assert.dir.exists(path)
  path
}
