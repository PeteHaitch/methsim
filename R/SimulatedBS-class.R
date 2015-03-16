### =========================================================================
### SimulatedBS: An S4 class to store data on a single simulated
### bisulfite-sequencing assay; BS = [WGBS | RRBS | eRRBS]
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### data.table("chr", readID", "pos", "z")
###   "seqname": The name of the seqlevel, i.e., chromosome.
###   "queryID": A unique (within a chromosome) read ID.
###       "pos": The position along the chromosome of the methylation locus.
###         "z": The methylation state at that locus.
###
### An object of this class is returned by the
### simulate,SimulateBSParam-method.
### This class will only retain those reads overlapping at least one
### methylation locus.

#' SimulatedBS class.
#' @aliases SimulatedBS
#' @details An S4 class to store a single simulated bisulfite-sequencing
#' experiment.
#' @export
setClass("SimulatedBS",
         contains = "data.table")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulatedBS.data.table <- function(object) {
  msg <- NULL

  if (!inherits(object, "data.table")) {
    msg <- Biobase::validMsg(msg, "'SimulatedBS' must inherit from data.table")
  } else {
    if (!identical(colnames(object), c("seqnames", "pos", "readID", "z"))) {
      msg <- Biobase::validMsg(msg, paste0("'SimulatedBS' colnames must be ",
                                           "'seqnames', 'pos', 'readID' and ",
                                           "'z'."))
    } else {
      if (class(object[["seqnames"]]) != "factor") {
        msg <- Biobase::validMsg(msg, "'seqnames' must be a'factor'.")
      }
      if (class(object[["pos"]]) != "integer") {
        msg <- Biobase::validMsg(msg, "'pos' must be an 'integer'.")
      }
      if (class(object[["readID"]]) != "integer") {
        msg <- Biobase::validMsg(msg, "'readID' must be an 'integer'.")
      }
      if (class(object[["z"]]) != "integer") {
        msg <- Biobase::validMsg(msg, "'z' must be an 'integer'.")
      }
    }
  }
  msg
}

.valid.SimulatedBS <- function(object) {

  # Include all other .valid.SimulatedBS.* functions in this vector
  msg <- c(.valid.SimulatedBS.data.table(object))

  if (is.null(msg)) {
    return(TRUE)
  } else{
    msg
  }
}

setValidity2("SimulatedBS", .valid.SimulatedBS)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# None because I don't want the user constructing these manually, rather they
# should be constructed by simulate,SimulateBSParam-method.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# TODO: Coerction from SimulatedBS object to MethPat, i.e.,
# as(SimulateBS, "MethPat", size = 1L).

