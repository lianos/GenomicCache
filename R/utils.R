na.logical <- function(the.logical) {
  the.logical <- as.logical(the.logical)
  the.logical[is.na(the.logical)] <- FALSE
  the.logical
}

filterByChr <- function(grl, which.chr=NULL) {
  if (!is.null(which.chr)) {
    keep <- sapply(grl, function(g) {
      length(g) > 0 && all(seqnames(g) == which.chr)
    })
    grl <- grl[keep]
  }
  grl
}

checkVerbose <- function(...) {
  verbose <- list(...)$verbose
  if (is.null(verbose)) verbose <- options()$verbose
  verbose
}

## Loads one item from the rda file. if what is null, it will
## load the first item 
load.it <- function(rda.file, what=NULL) {
  if (!file.exists(rda.file)) {
    stop("Can't find data file ", rda.file)
  }
  e <- new.env()
  vars <- load(rda.file, e)
  if (length(vars) == 0L) {
    stop("No objects found in ", rda.file)
  }
  if (is.null(what)) {
    what <- vars[1]
  }
  if (!what %in% vars) {
    stop("Object `", what, "` not found in ", rda.file)
  }
  get(what, e, inherits=FALSE)
}

##' Returns the bioconductor annoation package name for the given genome.
##' 
##' @param from A character string naming the genome, ie. hg18, mm9, etc.
##' The function also checks to see if it is the name of the package itself.
##' @param package Passed through to the \code{\link{annotationPackage}}
##' function. 
getAnnoPackageName <- function(from, package=NULL) {
  is.anno.package <- length(grep('^org\\..*\\.db$', from) == 1L)
  if (is.anno.package) {
    ## this is probably the package name itself
    if (!require(from, character.only=TRUE)) {
      stop("Unknown package: ", from)
    }
    from
  } else {
    ## probably the genome
    annotationPackage(from, package=package)
  }
}

##' Takes a character vector and tacks on a .1, .2, etc. to any duplicate
##' names until they are all made unique
##'
##' uniquefy(c('A', 'B', 'C', 'D', 'A', 'A', 'B', 'C', 'C')) will return:
##'     "A"   "B"   "C"   "D"   "A.1" "A.2" "B.1" "C.1" "C.2"
uniquefy <- function(values) {
  return(make.unique(values))
  ## dups <- duplicated(values)
  ## if (any(dups)) {
  ##   values[dups] <- paste(values[dups], '1', sep='.')
  ## }
  ## dups <- duplicated(values)
  ## while (any(dups)) {
  ##   count <- regexpr("\\.(\\d+)$", values[dups], perl=TRUE)
  ##   repl <- substring(values[dups], count + 1,
  ##                     count + attr(count, 'match.length'))
  ##   repl <- as.integer(repl) + 1
  ##   values[dups] <- paste(substring(values[dups], 1, count -1),
  ##                         repl, sep=".")
  ##   dups <- duplicated(values)
  ## }
  ## values
}

dir.exists <- function(path) {
  path <- as.character(path)
  if (is.na(file.info(path)$isdir) || !file.info(path)$isdir) {
    FALSE
  } else {
    TRUE
  }
}

assert.dir.exists <- function(path) {
  if (!dir.exists(path)) stop("Can't access directory: ", path)
}


## setAs("GRanges", "IntervalTree", function(from) {
##   as(ranges(from), "IntervalTree")
## })

##' Returns an object of type \code{type} from a list, this is most useful
##' when \code{the.list} has one object of \code{type} in it.
##' 
##' Primarily used to get arguments out of function calls with \code{(...)}
##' assumes tha
##'
##' If this object isn't found, or other error, returns \code{NULL}
takeFromListByType <- function(the.list, type, multi=FALSE, index=FALSE) {
  take <- which(sapply(the.list, function(arg) inherits(arg, type)))
  if (length(take) == 0L) {
    return(NULL)
  }
  
  if (length(take) > 1) {
    if (is.logical(multi[1])) {
      if (!multi) {
        warning("Multiple objects of type ", type, " found.")
        take <- '..NOTHING..'
      }
    } else if (is.numeric(multi)) {
      if (any(multi > length(take)) || any(multi < 0L)) {
        warning("multi take subscript(s) out of bounds")
        take <- '..NOTHING..'
      } else {
        take <- take[multi]
      }
    } else {
      warning("Illegal type of multi argument: ", is(multi)[1])
      take <- '..NOTHING..'
    }
  }

  if (index) {
    ret.val <- take
  } else {
    ret.val <- if (length(take) > 1) the.list[take] else the.list[[take]]
  }

  ret.val
}
