# Copyright (C) 2015 Peter Hickey
#
# This file is part of methsim.
#
# methsim is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# methsim is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with methsim  If not, see <http://www.gnu.org/licenses/>.

### =========================================================================
### findMostOverlapping
### -------------------------------------------------------------------------
###

# TODO (long term): Make this a generic method and suggest the addition of
# findMostOverlapping,IRanges-method to IRanges,
# findMostOverlapping,GRanges-method to GenomicRanges, etc. Perhaps do this via
# the 'select' parameter in findOverlaps, e.g.,
# findOverlaps(x, y, select = 'most').

#' Find most overlapping element in subject for each element of query.
#'
#' Runs the appropriate \code{findOverlaps} method and then for elements of
#' \code{query} that overlap multiple elements of \code{subject} chooses the
#' element of \code{subject} that has the greater overlap. Ties are broken at
#' random.
#'
#' @param query Any object for which a `findOverlaps` method is defined, e.g.,
#' a \code{\link[GenomicRanges]{GRanges}}.
#' @param subject Any object for which a `findOverlaps` method is defined, e.g.,
#' a \code{\link[GenomicRanges]{GRanges}}.
#' @param ignore.strand When set to \code{TRUE}, the strand information is
#' ignored in the overlap calculations.
#'
#' @note This was written for working with \code{\link[GenomicRanges]{GRanges}}
#' objects. While it may work with other objects for which a
#' \code{findOverlaps} method exists, these are untested and correctness is not
#' guaranteed. In particular, this will return an error if either the
#' \code{query} or \code{subject} are \code{\link[GenomicTuples]{GTuples}} or
#' \code{\link[GenomicTuples]{GTuplesList}} objects.
#'
#' @return An \code{integer} vector containing for each \code{query} the index
#' of the \code{subject} that it most overlaps (\code{NA} if there is no
#' overlap).
#'
#' @export
findMostOverlapping <- function(query, subject, ignore.strand = FALSE) {
  if (inherits(query, "GTuples") || inherits(query, "GTuplesList") ||
      inherits(subject, "GTuples") || inherits(subject, "GTuplesList")) {
    stop(paste0("'findMostOverlapping' not yet implemented for 'GTuples'- and ",
                "'GTuplesList'-based objects."))
  }
  ol <- findOverlaps(query, subject, ignore.strand = ignore.strand)
  mo <- rep(NA_integer_, )
  wol <- ranges(ol, ranges(query), ranges(subject))

  # Only run do this if at least one element of query has an overlap.
  if (any(countQueryHits(ol) != 0L)) {
    mo[countQueryHits(ol) > 0L] <- mapply(function(x, i) x[i],
                                          x = split(subjectHits(ol),
                                                    queryHits(ol)),
                                          i = tapply(width(wol),
                                                     queryHits(ol),
                                                     nnet::which.is.max),
                                          SIMPLIFY = TRUE, USE.NAMES = FALSE)
  }
  mo
}
