##' Takes a character vector and tacks on a .1, .2, etc. to any duplicate
##' names until they are all made unique
##'
##' uniquefy(c('A', 'B', 'C', 'D', 'A', 'A', 'B', 'C', 'C')) will return:
##'     "A"   "B"   "C"   "D"   "A.1" "A.2" "B.1" "C.1" "C.2"
uniquefy <- function(values) {
  dups <- duplicated(values)
  if (any(dups)) {
    values[dups] <- paste(values[dups], '1', sep='.')
  }
  dups <- duplicated(values)
  while (any(dups)) {
    count <- regexpr("\\.(\\d+)$", values[dups], perl=TRUE)
    repl <- substring(values[dups], count + 1,
                      count + attr(count, 'match.length'))
    repl <- as.integer(repl) + 1
    values[dups] <- paste(substring(values[dups], 1, count -1),
                          repl, sep=".")
    dups <- duplicated(values)
  }
  values
}

dir.exists <- function(path) {
  if (!is.character(path)) {
    stop("Illegal variable type for path: ", is(path)[1])
  }
  
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
    if (is.logical(multi[1]) && !multi) {
      warning("Multiple objects of type ", type, " found.")
      take <- '..NOTHING..'
    } else if (is.integer(multi)) {
      if (any(multi > length(take)) || any(multi < 0L)) {
        warning("multi take subscript(s) out of bounds")
        take <- '..NOTHING..'
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
