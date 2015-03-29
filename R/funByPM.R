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
### funByPM: Compute a function on a MethPat object for each sample stratified
### by its PartitionedMethylome.
### -------------------------------------------------------------------------
###

# TODO: Support other FUN, e.g., methLevel and methLevelCor
# TODO: Reconcile funByPM() with function-level stratification, e.g.
# cometh(feature).
#
#' Apply a function to a
#' \code{MethylationTuples::\link[MethylationTuples]{MethPat}} object,
#' with results stratified by the region type of the accompanying
#' \code{\link{PartitionedMethylome}}.
#'
#' @param FUN The function to be applied to each element of \code{pm}: see
#' 'Details'. Only a limited number of functions are currently supported,
#' specifically \code{MethylationTuples::\link[MethylationTuples]{methLevel}()},
#' \code{MethylationTuples::\link[MethylationTuples]{cometh}()}, and
#' \code{MethylationTuples::\link[MethylationTuples]{patternFreqs}()}.
#' @param pm A \code{\link{PartitionedMethylome}} object.
#' @param methpat A \code{MethylationTuples::\link[MethylationTuples]{MethPat}}
#' object.
#' @param min_cov An \code{integer} specifying the minimum coverage required
#' in order to use an m-tuple in the analysis.
#' @param ... Optional arguments to \code{FUN}.
#'
#' @return A \code{data.table::\link[data.table]{data.table}}.
#'
#' @export
funByPM <- function(FUN, pm, methpat, min_cov, ...) {
  # Argument checks
  stopifnot(is(pm, "PartitionedMethylome"))
  stopifnot(is(methpat, "MethPat"))
  if (ncol(methpat) > 1L) {
    stop("'methpat' must only contain data on a single sample.")
  }
  stopifnot(min_cov >= 0L)
  FUN <- match.fun(FUN)
  # TODO: Make more robust method for checking whether it is possible to use
  # FUN with methsim:::funByPM()
  AVAILABLE_FUNCTIONS <- c(MethylationTuples::methLevel,
                           MethylationTuples::cometh,
                           MethylationTuples::patternFreqs)
  stopifnot(any(vapply(X = AVAILABLE_FUNCTIONS, FUN = identical,
                       FUN.VALUE = logical(1), FUN)))

  # Drop m-tuples with insufficient coverage
  cov <- getCoverage(methpat)
  methpat <- methpat[as.vector(!is.na(cov) & (cov >= min_cov)), ]
  # Find the region type of each m-tuple
  ol <- findOverlaps(methpat, pm, type = "within")
  # Drop those m-tuples spanning a region boundary
  methpat <- methpat[queryHits(ol)]
  type <- regionType(pm[subjectHits(ol)])
  val <- FUN(methpat, min_cov = min_cov, ...)
  if (!is(val, "data.table")) {
    val <- data.table::as.data.table(val)
  }
  val[, type := type]
  val
}
