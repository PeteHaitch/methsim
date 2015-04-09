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
### SimulatedMethylome2: This is a complete re-thinking of the
###                      SimulatedMethylome class and requires changes to all
###                      downstream methods and functions.
###                      An S4 class to store the marginal probabilities and
###                      log odds-ratios that are required to compute the
###                      transition probabilities. These transition
###                      probabilities are ultimately used to sample from the
###                      'true' methylome in order to simulate an assay.
###                      A SimulatedMethylome2 object can contain data on
###                      multiple methylomes, provided all methylomes have
###                      the same set of methylation loci (or NA marginalProb
###                      at locus for sample without that locus) and the
###                      same seqinfo.
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### (MTuples | marginal, LOR)
###          MTuples: A MethylationTuples::MTuples object
###     marginalProb: An vector of the marginal probabilities that the i-th
###                   methylation locus is methylated, i.e., Pr(Z_{i} = 1).
###              LOR: A vector of the log odds ratio (LOR) for the (i-1)-th and
###                   i-th methylation locus. The LOR is 0 for the first
###                   locus on each seqlevel, corresponding to independence.

#' SimulatedMethylome2 class.
#' @aliases SimulatedMethylome2
#' @details An S4 class to store a simulated methylomes.
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

  if (!identical(assayNames(object), c("marginalProb", "LOR"))) {
    msg <- Biobase::validMsg(msg, paste0("'assayNames' of a ",
                                         "'SimulatedMethylome2' object must ",
                                         "be 'marginalProb' and 'LOR'."))
  }
  if (!is.numeric(assay(object, "marginalProb", withDimnames = FALSE)) ||
      any(is.na(assay(object, "marginalProb", withDimnames = FALSE))) ||
      min(assay(object, "marginalProb", withDimnames = FALSE) < 0) ||
      max(assay(object, "marginalProb", withDimnames = FALSE) > 1)) {
    msg <- Biobase::validMsg(msg, paste0("'marginalProb' assay must not ",
                                         "contain 'NA' and all values must be ",
                                         "probabilities (i.e., between zero ",
                                         " and one)."))
  }
  if (!is.numeric(assay(object, "LOR", withDimnames = FALSE)) ||
      any(is.na(assay(object, "LOR", withDimnames = FALSE)))) {
    msg <- Biobase::validMsg(paste0("'LOR' assay must be numeric and not ", "
                                    contain 'NA'."))
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
