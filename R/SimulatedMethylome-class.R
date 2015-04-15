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
### SimulatedMethylome: An S4 class to store the marginal probabilities,
###                     log odds-ratios and weights that represent a simulated
###                     'true' methylome. A SimulatedMethylome may only
###                     contain data for a single methylome.
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### (MTuples | MarginalProb, LOR, MixtureWeights)
###          MTuples: A MethylationTuples::MTuples object
###     MarginalProb: An vector of the marginal probabilities that the i-th
###                   methylation locus is methylated, i.e., Pr(Z_{i} = 1).
###              LOR: A vector of the log odds ratio (LOR) for the (i-1)-th and
###                   i-th methylation locus. The LOR is 0 for the first
###                   locus on each seqlevel, corresponding to independence.
###   MixtureWeights: A DataFrame specifying the weights of the mixture
###                   distribution. Should be identical be for all loci.

#' SimulatedMethylome class.
#' @aliases SimulatedMethylome
#' @details An S4 class to store a simulated methylomes.
#' @export
setClass("SimulatedMethylome",
         contains = "SummarizedExperiment")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulatedMethylome.rowData <- function(object) {
  msg <- NULL

  if (!is(object@rowData, "MTuples")) {
    msg <- Biobase::validMsg(msg, paste0("'rowData' slot of a ",
                                         "'SimulatedMethylome' object must ",
                                         "be a 'MTuples' object."))
  } else if (size(object@rowData) != 1L) {
    # Only run this check if the rowData is indeed an MTuples object.
    msg <- Biobase::validMsg(msg, paste0("'size' of 'MTuples' in 'rowData' ",
                                         "slot must be 1."))
  }
  msg
}

.valid.SimulatedMethylome.assays <- function(object) {
  msg <- NULL

  if (!identical(assayNames(object),
                 c("MarginalProb", "LOR", "MixtureWeights"))) {
    msg <- Biobase::validMsg(msg, paste0("'assayNames' of a ",
                                         "'SimulatedMethylome' object must ",
                                         "be 'marginalProb', 'LOR' and ",
                                         "'MixtureWeights'."))
  }
  if (!is.numeric(assay(object, "MarginalProb", withDimnames = FALSE)) ||
      any(is.na(assay(object, "MarginalProb", withDimnames = FALSE))) ||
      min(assay(object, "MarginalProb", withDimnames = FALSE) < 0) ||
      max(assay(object, "MarginalProb", withDimnames = FALSE) > 1)) {
    msg <- Biobase::validMsg(msg, paste0("'MarginalProb' assay must not ",
                                         "contain 'NA' and all values must be ",
                                         "probabilities (i.e., between zero ",
                                         " and one)."))
  }
  if (!is.numeric(assay(object, "LOR", withDimnames = FALSE)) ||
      any(is.na(assay(object, "LOR", withDimnames = FALSE)))) {
    msg <- Biobase::validMsg(msg, paste0("'LOR' assay must be numeric and ",
                                         "not contain 'NA'."))
  }
  if (!is(assay(object, "MixtureWeights", withDimnames = FALSE), "DataFrame") ||
      any(sapply(assay(object, "MixtureWeights", withDimnames = FALSE),
                 function(x) length(S4Vectors::runValue(x)) != 1L)) ||
      !all.equal(sum(sapply(assay(object, "MixtureWeights",
                                  withDimnames = FALSE),
                 function(x) as.vector(x[1]))), 1)) {
    msg <- Biobase::validMsg(msg, "'MixtureWeights' must sum to 1.")
  }
  msg
}

.valid.SimulatedMethylome <- function(object) {

  # First need to check that rowData is an MTuples object with 1-tuples.
  # Otherwise some of the .valid.SimulatedMethylome2.* functions won't work
  msg <- .valid.SimulatedMethylome.rowData(object)
  if (is.null(msg)) {

    # Include all other .valid.SimulatedMethylome.* functions in this vector
    msg <- c(.valid.SimulatedMethylome.assays(object))
  }

  if (is.null(msg)) {
    return(TRUE)
  } else{
    msg
  }
}

S4Vectors::setValidity2("SimulatedMethylome", .valid.SimulatedMethylome)

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

# TODO: Is it possible to compute methLevel from SimulatedMethylome without
# explicity simulating Z?

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### cometh()
###

# TODO: Is it possible to compute cometh from SimulatedMethylome without
# explicity simulating Z?

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methinfo()
###

#' @export
setMethod("methinfo",
          "SimulatedMethylome",
          function(object) {
            methinfo(object@rowData)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getters/setters
###
# TODO (long term)
