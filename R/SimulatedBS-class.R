### =========================================================================
### SimulatedBS: An S4 class to store data on a single simulated
### bisulfite-sequencing assay; BS = [WGBS | RRBS | eRRBS]
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### list(z, seqinfo)
###        z: A list of data.tables, each element of the list is data for a
###           single seqlevel (and the names of the list are the seqlevels).
###        Each data.table has the following column names
###          "queryID": A unique (within a chromosome) read ID.
###              "pos": The position along the chromosome of the methylation
###                     locus.
###                "z": The methylation state at that locus.
###  seqinfo: A Seqinfo object for the SimulatedMethylome from which this
###           SimulatedBS object was simulated.
### methinfo: A MethInfo object for the SimulatedMethyome from which this
###           SimulatedBS object was simulated.
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
         slots = list(
           z = "list",
           seqinfo = "Seqinfo",
           methinfo = "MethInfo"
         )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulatedBS.z <- function(object) {
  msg <- NULL

  if (!is(object@z, "list")) {
    msg <- Biobase::validMsg(msg, "'SimulatedBS' must be a list.")
  } else {
    if (any(sapply(object@z, function(x) !is(x, "data.table")))) {
      msg <- Biobase::validMsg(msg, paste0("All elements of 'z' slot ",
                                           "must be 'data.table' objects."))
    }
    if (any(sapply(object@z, function(x) {
      !identical(colnames(x), c("pos", "readID", "z"))
    }))) {
      msg <- Biobase::validMsg(msg, paste0("colnames of all elements of ",
                                           "'SimulatedBS' must be 'pos', ", "
                                           'readID' and 'z'."))
    } else {
      if (any(sapply(object@z, function(x) class(x[["pos"]]) != "integer"))) {
        msg <- Biobase::validMsg(msg, "'pos' columns must be 'integer' type.")
      }
      if (any(sapply(object@z, function(x) class(x[["readID"]]) != "integer")))
        {
        msg <- Biobase::validMsg(msg,
                                 "'readID' columns must be 'integer' type.")
      }
      if (any(sapply(object@z, function(x) class(x[["z"]]) != "integer"))) {
        msg <- Biobase::validMsg(msg, "'z' columns must be 'integer' type.")
      }
    }
  }
  msg
}

.valid.SimulatedBS.seqinfo <- function(object) {
  msg <- NULL

  if (!is(object@seqinfo, "Seqinfo")) {
    msg <- Biobase::validMsg(msg, "'seqinfo' slot must be a Seqinfo object.")
  }
  if (!identical(seqlevels(object@seqinfo), names(object@z))) {
    msg <- Biobase::validMsg(msg, paste0("'seqinfo' slot missing names ",
                                         "(seqlevels) of 'z' slot"))
  }
  msg
}

.valid.SimulatedBS.methinfo <- function(object) {
  msg <- NULL

  if (!is(object@methinfo, "MethInfo")) {
    msg <- Biobase::validMsg(msg, "'methinfo' slot must be a MethInfo object.")
  }
  msg
}

.valid.SimulatedBS <- function(object) {

  # Include all other .valid.SimulatedBS.* functions in this vector
  msg <- c(.valid.SimulatedBS.z(object),
           .valid.SimulatedBS.seqinfo(object),
           .valid.SimulatedBS.methinfo(object))

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
### seqinfo()
###

setMethod("seqinfo",
          "SimulatedBS",
          function(object) {
            object@seqinfo
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methinfo()
###

setMethod("methinfo",
          "SimulatedBS",
          function(object) {
            object@methinfo
          }
)

### Coercion
###

# TODO: Coerction from SimulatedBS object to MethPat, i.e.,
# as(SimulateBS, "MethPat", size = 1L).

