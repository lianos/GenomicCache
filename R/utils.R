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
