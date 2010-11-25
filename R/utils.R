##' Checks that a directory exists and will create it if not.
##'
##' If the directory does not exist, and the caller does not want to create it
##' an error will be thrown
##'
##' @export
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##' 
##' @param path The path to the directory to check.
##' @param create A logical indicating whether or not the directory should be
##' created if it doesn't exist
##' @param verbose Let us know what's going on
##'
##' @return \code{TRUE} if everything is kosher, otherwise an error is thrown.
checkOrCreateDirectory <- function(path, create=FALSE, verbose=TRUE) {
  if (!dir.exists(path)) {
    if (!create) {
      stop("Directory", path, "does not exist", sep=" ")
    } else {
      if (verbose) cat("Creating directory", path, "...\n")
      if (!dir.create(save.path)) {
        stop("Error! Check permissions? Parent directory exists?")
      }
    }
  }

  TRUE
}

##' Convenience method to sets \code{NA}'s in a logical vector to \code{FALSE}.
##' 
##' @export
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##' 
##' @param the.logical A logical vector/Rle
##' @return A \code{logical} with \code{NA} values set to \code{FALSE}
na.logical <- function(the.logical) {
  the.logical <- as.logical(the.logical)
  the.logical[is.na(the.logical)] <- FALSE
  the.logical
}

##' Convert NA values in vectors and data.frames to a default value
##'
##' @param wut The thing to convert
##' @param to The value convert NA to
##' @return The same type as \code{wut}
convert.na <- function(wut, to=".defaults.") {
  if (is.character(to) && to[1] == ".defaults.") {
    to <- list(logical=FALSE, numeric=0, integer=0L, character="",
               factor="")
  }
  if (is.vector(wut) || is.factor(wut)) {
    wut.type <- is(wut)[1]
    if (is.list(to)) {
      if (!wut.type %in% names(to)) {
        stop("Unknown default conversion value for", wut.type, sep=" ")
      }
      to <- to[[wut.type]]
    }
    if (wut.type == 'factor') {
      levels(wut) <- c(levels(wut), to)
    }
    wut[is.na(wut)] <- to
  } else if (inherits(wut, 'data.frame') || inherits(wut, 'DataFrame')) {
    cols <- 1:ncol(wut)
    if (is(wut, 'data.table')) {
      ## Don't change key columns
      cols <- setdiff(cols, which(colnames(wut) %in% key(wut)))
    }
    for (idx in cols) {
      wut[[idx]] <- convert.na(wut[[idx]], to=to)
    }
  }
  
  wut
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

dir.exists <- function(path) {
  path <- as.character(path)
  info <- file.info(path)
  !is.na(info$isdir) && info$isdir
}

assert.dir.exists <- function(path) {
  if (!dir.exists(path)) stop("Can't access directory: ", path)
}

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
