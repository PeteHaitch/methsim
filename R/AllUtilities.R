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
### Utility functions. These functions are not exported.
### NB: This filename is a misnomer since there are other utility functions
### defined in R/sampling-utilities.R
### -------------------------------------------------------------------------


### -------------------------------------------------------------------------
### .makePosAndCounts(): A helper function called by asMethPat
###

#' A helper function called by asMethPat
#' @keywords internal
#' @export
.makePosAndCounts <- function(zz, size) {

  if (size > 1L) {

    # Special handling for the case where there are no reads mapped to this
    # particular seqlevel.
    if (nrow(zz) == 0L) {
      pos <- matrix(integer(0), ncol = size)
      # TODO: How to return counts
      counts <- matrix(integer(0), ncol = 2 ^ size,
                       dimnames =
                         list(NULL,
                              MethylationTuples:::.makeMethPatNames(size)))
    } else {
      # Remove reads with less than 'size' methylation loci.
      setkey(zz, readID)
      zz_reduced <- zz[, n := .N, by = key(zz)][n >= size, ][, n := NULL]

      # Special handling of the case when there are no reads with sufficient
      # methylation loci.
      # NB: Slightly different to the case of there being no reads mapped to
      # this particular seqlevel.
      if (nrow(zz_reduced) == 0L) {
        pos <- matrix(integer(0), ncol = size)
        # TODO: How to return counts
        counts <- matrix(integer(0), ncol = 2 ^ size,
                         dimnames =
                           list(NULL,
                                MethylationTuples:::.makeMethPatNames(size)))
      } else {
        # Create m-tuples from each read (where m = size).
        setkey(zz_reduced, readID, pos)
        # Tabulate methylation patterns at each m-tuple.
        # UP TO HERE: .tabulatez() is causing coredump.
        counts <- .tabulatez(zz_reduced[, readID],
                             zz_reduced[, z],
                             zz_reduced[, pos],
                             size)
        pos <- strsplit(names(counts), ",")
        # Re-structure counts into a list of counts for each methylation
        # pattern.
        counts <- lapply(seq_len(2 ^ size), function(i, counts) {
          unlist(lapply(counts, "[", i), use.names = FALSE)
        }, counts = counts)
        names(counts) <- MethylationTuples:::.makeMethPatNames(size)
        # Drop zero 'rows', i.e., those without counts
        nonzero_rows <- Reduce(`+`, counts) > 0
        counts <- lapply(counts, "[", nonzero_rows)
        pos <- matrix(as.integer(unlist(pos[nonzero_rows],
                                        use.names = FALSE)),
                      ncol = size,
                      byrow = TRUE)
        # Order data by positions (pos is currently sorted lexicographically)
        order_idx <- do.call(order, lapply(seq_len(ncol(pos)),
                                           function(i) pos[, i]))
        pos <- pos[order_idx, ]
        counts <- lapply(counts, "[", order_idx)
      }
    }
  } else {
    setkey(zz, pos)
    zz_reduced <- zz[, list(M = sum(z), U = sum(!z)), by = key(zz)]

    pos <- as.matrix(zz_reduced[, pos])
    # TODO: How to return counts?
    counts <- list(M = as.matrix(zz_reduced[, M]),
                   U = as.matrix(zz_reduced[, U]))
  }

  list(pos = pos, counts = counts)
}
