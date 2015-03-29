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

#' Find "most overlapping" range in subject for each element of query.
#' @param query Any object for which a `findOverlaps` method is defined, e.g.,
#' a \code{\link[GenomicRanges]{GRanges}}.
#' @param subject Any object for which a `findOverlaps` method is defined, e.g.,
#' a \code{\link[GenomicRanges]{GRanges}}.
#' @details Ties are broken at random.
#' @return An \code{integer} vector containing for each \code{query} the index
#' of the \code{subject} that it most overlaps.
#' @export
findMostOverlapping <- function(query, subject) {
  ol <- findOverlaps(query, subject)
  wol <- ranges(ol, ranges(query), ranges(subject))
  mapply(function(x, i) x[i], x = split(subjectHits(ol), queryHits(ol)),
         i = tapply(width(wol), queryHits(ol), nnet::which.is.max),
         SIMPLIFY = TRUE, USE.NAMES = FALSE)
}
