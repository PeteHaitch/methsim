### =========================================================================
### SimulatedMethylome: An S4 class to store a single simulated methylome.
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### (MTuples | Z, H)
###     MTuples: A MethylationTuples::MTuples object
###     Z: A DataFrame/sparseMatrix/matrix (TBD) of the methylation haplotypes.
###     H: A DataFrame of the frequency of each haplotype.

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
  } else if (size(object@rowData != 1L)) {
    # Only run this check if the rowData is indeed an MTuples object.
    msg <- Biobase::validMsg(msg, paste0("'size' of 'MTuples' in 'rowData' ",
                                         "slot must be 1."))
  }
  msg
}

.valid.SimulatedMethylome.assays <- function(object) {
  msg <- NULL

  if (!identical(GenomicRanges::assayNames(object), c("Z", "H"))) {
    msg <- Biobase::validMsg(msg, paste0("'assayNames' of a ",
                                         "'SimulatedMethylome' object must be ",
                                         "'Z' and 'H'"))
  } else if (ncol(Z) != ncol(H)) {
    # Only run this check if the assay names are valid
    msg <- Biobase::validMsg(msg, paste0("'ncol(Z)' must equal 'ncol(H)'."))
  }
  msg
}

.valid.SimulatedMethylome <- function(object) {

  # First need to check that rowData is an MTuples object with 1-tuples.
  # Otherwise some of the .valid.SimulatedMethylome.* functions won't work
  msg <- .valid.SimulatedMethylome.rowData(object)
  if (is.null(msg)){

    # Include all other .valid.MethPat.* functions in this vector
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
