### =========================================================================
### SimulatedMethylome: An S4 class to store a single simulated methylome.
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### (MTuples | Z, W)
###     MTuples: A MethylationTuples::MTuples object
###     Z: A matrix (TBD) of the methylation haplotypes.
###     W: A DataFrame of the weights (relative frequency) of each haplotype.

#' SimulatedMethylome class.
#' @aliases SimulatedMethylome
#' @details An S4 class to store a single simulated methylome.
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
                                         "'SimulatedMethylome' object must be ",
                                         "a 'MTuples' object."))
  } else if (size(object@rowData) != 1L) {
    # Only run this check if the rowData is indeed an MTuples object.
    msg <- Biobase::validMsg(msg, paste0("'size' of 'MTuples' in 'rowData' ",
                                         "slot must be 1."))
  }
  msg
}

.valid.SimulatedMethylome.assays <- function(object) {
  msg <- NULL

  if (!identical(names(object@assays$data), c("Z", "W"))) {
    msg <- Biobase::validMsg(msg, paste0("'assayNames' of a ",
                                         "'SimulatedMethylome' object must be ",
                                         "'Z' and 'W'"))
  } else if (ncol(object@assays$data$Z) != ncol(object@assays$data$W)) {
    # Only run this check if the assay names are valid
    msg <- Biobase::validMsg(msg, paste0("'ncol(Z)' must equal 'ncol(W)'."))
  }
  msg
}

.valid.SimulatedMethylome <- function(object) {

  # First need to check that rowData is an MTuples object with 1-tuples.
  # Otherwise some of the .valid.SimulatedMethylome.* functions won't work
  msg <- .valid.SimulatedMethylome.rowData(object)
  if (is.null(msg)){

    # Include all other .valid.SimulatedMethylome.* functions in this vector
    msg <- c(.valid.SimulatedMethylome.assays(object))
  }

  if (is.null(msg)) {
    return(TRUE)
  } else{
    msg
  }
}

setValidity2("SimulatedMethylome", .valid.SimulatedMethylome)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# None because I don't want the user constructing these manually, rather they
# should be constructed by simulate,SimulateMethylomeParam-method.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methLevel()
###

## TODO: Document.
## TODO: Unit tests.
#' Compute methylation levels.
#' @param object A \code{\link{SimulatedMethylome}} object.
#' @param statistic A \code{character} string indicating which methylation
#' level statistic is to be computed. One of "\code{beta-values}" or
#' "\code{M-values}" (see below).
#' @param offset A \code{numeric} vector with length 1 used when computing
#' M-values (default: 1).
#'
#' @details
#' TODO: Define beta-values and M-values. Note any differences with how others
#' define beta-values or M-values, e.g., minfi.
#'
#' @return A \code{\link[base]{matrix}}, with the same dimensions and dimension
#' names as \code{x}, of methylation levels at each methylation loci in each
#' sample.
#'
#' @aliases methLevel
#'
#' @export
setMethod("methLevel",
          "SimulatedMethylome",
          function(object, statistic = c("beta-values", "M-values"),
                   offset = 1L) {
            statistic <- match.arg(statistic)
            if (statistic == "beta-values") {
              meth_level <- rowSums(assay(object, "W", withDimnames = FALSE) *
                                      assay(object, "Z", withDimnames = FALSE))
            } else if (statistic == "M-values") {
              # TODO: Figure out how to compute M-values from beta-values
              # using the appropriate 'offset'.
              stop("Sorry, M-values not yet implemented.")
            }

            # TODO: Try to do this in a way that avoids a copy.
            matrix(meth_level, ncol = 1L)
          }
)

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
          "SimulatedMethylome",
          function(object) {
            methinfo(object@rowData)
          }
)
