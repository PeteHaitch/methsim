### =========================================================================
### findMostOverlapping
### -------------------------------------------------------------------------
###

#' Find "most overlapping" range in subject for each element of query.
#' @param query Any object for which a `findOverlaps` method is defined, e.g.,
#' a \code{\link[GenomicRange]{GRanges}}.
#' @param subject Any object for which a `findOverlaps` method is defined, e.g.,
#' a \code{\link[GenomicRange]{GRanges}}.
#' @details Ties are broken at random.
#' @return An \code{integer} vector containing for each \code{query} the index
#' of the \code{subject} that it most overlaps.
#' @export
findMostOverlapping <- function(query, subject) {
  stopifnot(isDisjoint(subject))
  ol <- findOverlaps(query, subject)
  wol <- ranges(ol, ranges(query), ranges(subject))
  x <- data.table(qh = queryHits(ol), sh = subjectHits(ol), w = width(wol),
                  key = "qh")
  x[, sh[nnet::which.is.max(w)], by = qh][, V1]
}
