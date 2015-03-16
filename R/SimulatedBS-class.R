### =========================================================================
### SimulatedBS: An S4 class to store data on a single simulated
### bisulfite-sequencing assay; BS = [WGBS | RRBS | eRRBS]
### -------------------------------------------------------------------------
###

# UP TO HERE: Ah crap, `[` breaks when the parent of an S4 class is
# data.table (see http://r.789695.n4.nabble.com/Weird-behavior-with-S4-subclasses-of-data-table-after-loading-RCurl-td4644611.html).

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### list(dt)
### dt: A data.table object with the following column names.
###   "seqnames": The name of the seqlevel, i.e., the chromosome.
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
         slots = list(
           dt = "data.table")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulatedBS.dt <- function(object) {
  msg <- NULL

  if (!inherits(object@dt, "data.table")) {
    msg <- Biobase::validMsg(msg, "'dt' slot must be a data.table")
  } else {
    if (!identical(colnames(object@dt), c("seqnames", "pos", "readID", "z"))) {
      msg <- Biobase::validMsg(msg, paste0("colnames of 'dt' slot must be ",
                                           "'seqnames', 'pos', 'readID' and ",
                                           "'z'."))
    } else {
      if (class(object@dt[["seqnames"]]) != "factor") {
        msg <- Biobase::validMsg(msg, "'seqnames' must be a 'factor'.")
      }
      if (class(object@dt[["pos"]]) != "integer") {
        msg <- Biobase::validMsg(msg, "'pos' must be an 'integer'.")
      }
      if (class(object@dt[["readID"]]) != "integer") {
        msg <- Biobase::validMsg(msg, "'readID' must be an 'integer'.")
      }
      if (class(object@dt[["z"]]) != "integer") {
        msg <- Biobase::validMsg(msg, "'z' must be an 'integer'.")
      }
    }
  }
  msg
}

.valid.SimulatedBS <- function(object) {

  # Include all other .valid.SimulatedBS.* functions in this vector
  msg <- c(.valid.SimulatedBS.dt(object))

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

