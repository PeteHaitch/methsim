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
### SimulatedMethylome2: An S4 class to store the transition probabilities
###                      to simulate a true methylome. This is a re-thinking
###                      of the SimulateMethylome class.
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### (MTuples | P)
###     MTuples: A MethylationTuples::MTuples object
###     P: A n \times 2 matrix of the transition probabilities.
###        Let P_{i, .} be the i-th row of (i = 1, ..., n),
###        then P_{i, .} =
###         [Pr(Z_{i} = 1 | Z_{i - 1} = 0), Pr(Z_{i} = 1 | Z_{i - 1} = 1)].
###       More generally, we might use a scheme where the j-th column of P
###       (j = 1, ..., 2^{m - 1}, with m = 2 under the default) stores
###       Pr(Z_{i} = 1 | Z_{i - 1}, Z_{i - 2}, ...).

#' SimulatedMethylome class.
#' @aliases SimulatedMethylome
#' @details An S4 class to store a single simulated methylome.
#' @export
setClass("SimulatedMethylome2",
         contains = "SummarizedExperiment")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulatedMethylome2.rowData <- function(object) {
  msg <- NULL

  if (!is(object@rowData, "MTuples")) {
    msg <- Biobase::validMsg(msg, paste0("'rowData' slot of a ",
                                         "'SimulatedMethylome2' object must ",
                                         "be a 'MTuples' object."))
  } else if (size(object@rowData) != 1L) {
    # Only run this check if the rowData is indeed an MTuples object.
    msg <- Biobase::validMsg(msg, paste0("'size' of 'MTuples' in 'rowData' ",
                                         "slot must be 1."))
  }
  msg
}

.valid.SimulatedMethylome2.assays <- function(object) {
  msg <- NULL

  if (!identical(names(object@assays$data), "P")) {
    msg <- Biobase::validMsg(msg, paste0("'assayNames' of a ",
                                         "'SimulatedMethylome2' object must ",
                                         "be 'P'."))
  }
  msg
}

.valid.SimulatedMethylome2 <- function(object) {

  # First need to check that rowData is an MTuples object with 1-tuples.
  # Otherwise some of the .valid.SimulatedMethylome2.* functions won't work
  msg <- .valid.SimulatedMethylome2.rowData(object)
  if (is.null(msg)) {

    # Include all other .valid.SimulatedMethylome.* functions in this vector
    msg <- c(.valid.SimulatedMethylome2.assays(object))
  }

  if (is.null(msg)) {
    return(TRUE)
  } else{
    msg
  }
}

S4Vectors::setValidity2("SimulatedMethylome2", .valid.SimulatedMethylome2)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# None because I don't want the user constructing these manually, rather they
# should be constructed by simulate2,SimulateMethylomeParam-method.
#
#

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methLevel()
###

# TODO

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### cometh()
###

# TODO

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methinfo()
###

# TODO: Document
#' Extract MethInfo
#'
#' @aliases methinfo
#'
#' @export
setMethod("methinfo",
          "SimulatedMethylome2",
          function(object) {
            methinfo(object@rowData)
          }
)
