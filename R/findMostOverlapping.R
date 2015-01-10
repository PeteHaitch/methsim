### =========================================================================
### findMostOverlapping
### -------------------------------------------------------------------------
###

# TODO: Works when copy-pasted but not as part of the package.
# TODO: Read and understand
# http://stackoverflow.com/questions/16573995/subset-by-group-with-data-table
# TODO: This code might not handle ties.

#' Find "most overlapping" range in subject for each element of query
#' @param query Any object for which a `findOverlaps` method is defined, e.g.,
#' a \code{\link[GenomicRange]{GRanges}}.
#' @param subject Any object for which a `findOverlaps` method is defined, e.g.,
#' a \code{\link[GenomicRange]{GRanges}}.
#' @return An \code{integer} vector.
#' @export
findMostOverlapping <- function(query, subject) {
  stopifnot(isDisjoint(subject))
  ol <- findOverlaps(query, subject)
  wol <- ranges(ol, ranges(query), ranges(subject))
  x <- data.table("sh" = subjectHits(ol), "qh" = queryHits(ol),
                  "w" = width(wol))
  # The following line works interactively but not in a function.
  # x[x[, .I[w == max(w)], by = qh]$V1][, sh]

  # Here are a few attempts to make "quote()-ed" versions of the above, but
  # without success.
  #   x[x[, eval(quote(.I[w == max(w)])),
  #       by = eval(quote(qh))]$V1][, eval(quote(sh))]  q1 <- substitute(.I[w == max(w)])
#   q1 <- substitute(.I[w == max(w)])
#   q2 <- substitute(qh)
#   q3 <- substitute(sh)
#   print(head(x))
#   x[x[, eval(q1), by = eval(q2)]$V1][, eval(q3)]

#   setkey(x, "sh")
#   x[x[, .I[w == max(w)], by = qh][, V1]][, sh]

  # This is slightly faster but still doesn't work in a function
  x[, sh[w == max(w)], by = qh][, V1]


}
