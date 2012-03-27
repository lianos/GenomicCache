# setClass("BigWigFile",
#          representation=representation(
#            path='character'),
#          prototype=prototype(
#            path=character()))
# 
# setMethod("initialize", "BigWigFile",
# function(.Object, ..., path=character()) {
#   callNextMethod(.Object, path=path, ...)
# })
# 
# setMethod("show", c(object="BigWigFile"),
# function(object) {
#   sep <- .Platform$file.sep
#   paths <- strsplit(dirname(object@path), sep, fixed=TRUE)[[1]]
#   if (length(paths) > 3) {
#     path <- paste(tail(paths, 3), collapse=sep)
#     path <- paste("...", path, sep=sep)
#   }
#   cat("BigWigFile [", file.path(path, basename(object@path)), "]\n", sep="")
# })
# 
# BigWigFile <- function(path) {
#   if (!file.exists(path)) {
#     stop("Can't find bigWig file: ", path)
#   }
# 
#   new("BigWigFile", path=path)
# }
# 
# ##' Query ranges for their bigWig score
# ##'
# ##' @importFrom rtracklayer import.bw
# setMethod("[", c(x="BigWigFile"),
# function(x, i, j, ..., drop) {
#   if (missing(i) || !missing(j) || length(list(...)) > 0L) {
#     stop("invalid subsetting")
#   }
#   if (!(inherits(i, 'GRanges') || inherits(i, 'RangesList'))) {
#     stop("can only subset with a GRanges-like object")
#   }
#   selection <- BigWigSelection(ranges=i, colnames='score')
#   scores <- import.bw(x@path, selection)
# 
#   if (inherits(i, 'GRanges')) {
#     strand(i) <- '*'
#     suppressWarnings(scores <- as(scores, 'GRanges'))
#     o <- findOverlaps(i, scores)
#     values(scores)$queryIdx <- queryHits(o)
#   }
# 
#   scores
# })
