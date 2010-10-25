## Taken from findOverlaps-GenomicRanges,GRangesList
## much better than I could ever do ...
setMethod("findOverlaps", c("Ranges", "RangesList"),
function(query, subject, maxgap=0L, minoverlap=1L,
         type=c("any", "start", "end", "within", "equal"),
         select=c("all", "first", "last", "arbitrary"), ...) {
  if (!IRanges:::isSingleNumber(maxgap) || maxgap < 0) {
    stop("'maxgap' must be a non-negative integer")
  }
  type <- match.arg(type)
  select <- match.arg(select)
  
  unlistSubject <- unlist(subject, use.names=FALSE)
  subjectGroups <- togroup(subject@partitioning)
  if (type == "start") {
    keep <- which(IRanges:::diffWithInitialZero(subjectGroups) != 0L)
    unlistSubject <-  unlistSubject[keep]
    subjectGroups <- subjectGroups[keep]
  } else if (type == "end") {
    keep <- end(subject@partitioning)[elementLengths(subject) > 0L]
    unlistSubject <-  unlistSubject[keep]
    subjectGroups <- subjectGroups[keep]
  }
  ans <- callGeneric(query, unlistSubject, maxgap=maxgap, type=type,
                     select="all")
  matchMatrix <- ans@matchMatrix
  if (minoverlap > 1L && nrow(ans@matchMatrix) > 0) {
    matchMatrix <- ans@matchMatrix[FALSE,]
    intrsct <- pintersect(ranges(query)[queryHits(ans)],
                          ranges(unlistSubject[subjectHits(ans)]))
    df <- data.frame(query = queryHits(ans),
                     s = subjectHits(ans),
                     subject =subjectGroups[subjectHits(ans)],
                     w=width(intrsct))
    mat <-  with(df, aggregate( w, list(query=query, subject=subject), sum))
    indx <- mat$x >= minoverlap
    if(any(indx))
      matchMatrix <- as.matrix(mat[indx, c("query", "subject"), drop=FALSE],
                               rownames.force = FALSE)
  } else {
    matchMatrix[, 2L] <- subjectGroups[matchMatrix[, 2L, drop=TRUE]]
    matchMatrix <- GenomicRanges:::.cleanMatchMatrix(matchMatrix)
  }
  if (select == "all") {
    DIM <- c(length(query), length(subject))
    initialize(ans, matchMatrix = matchMatrix, DIM = DIM)
  } else {
    IRanges:::.matchMatrixToVector(matchMatrix, length(query))
  }
})
