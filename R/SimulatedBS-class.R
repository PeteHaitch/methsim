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
### SimulatedBS: An S4 class to store data on a single simulated
### bisulfite-sequencing assay; BS = [WGBS | RRBS | eRRBS]
### -------------------------------------------------------------------------
###

# TODO: Add sampleName slot.

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
           methinfo = "MethInfo",
           SampleName = "character"
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
  if (!all(names(object@z) %in% seqlevels(object@seqinfo))) {
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

.valid.SimulatedBS.SampleName <- function(object) {
  msg <- NULL

  if (!is.character(object@SampleName) | length(object@SampleName) != 1L) {
    msg <- Biobase::validMsg(msg, paste0("'SampleName' slot must be a ",
                                         "'character' vector with length 1."))
  }
}

.valid.SimulatedBS <- function(object) {

  # Include all other .valid.SimulatedBS.* functions in this vector
  msg <- c(.valid.SimulatedBS.z(object),
           .valid.SimulatedBS.seqinfo(object),
           .valid.SimulatedBS.methinfo(object),
           .valid.SimulatedBS.SampleName(object))

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
# should be constructed by simulate,BSParam-method.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo()
###

setMethod("seqinfo",
          "SimulatedBS",
          function(x) {
            x@seqinfo
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# TODO: Deprecate sampleName argument once added as slot to SimulatedBS.
# Coercion from SimulatedBS object to MethPat
# TODO: Can't use setAs() because it doesn't allow extra arguments (i.e,
# 'size'). Could perhaps make MethPat constructor a generic (ala
# SummarizedExperiment) and allow MethPat(Simulated, size).

#' Coerce a \code{\link{SimulatedBS}} object to a
#' \code{\link[MethylationTuples]{MethPat}} object.
#'
#' @param from a \code{\link{SimulatedBS}} object.
#' @param size an integer specifying the \code{\link[MethylationTuples]{size}}
#' of the returned \code{\link[MethylationTuples]{MethPat}} object.
#' @return A \code{\link[MethylationTuples]{MethPat}} object.
#'
#' @note When \code{size} > 1, only adjacent m-tuples are created. This
#' function will preserve the \code{seed} attribute of the \code{SimulatedBS}
#' object.
#'
#' @rdname SimulatedBS-class
#' @name asMethPat
#' @export
asMethPat <- function(SimulatedBS, size = 1L, BPPARAM = bpparam()) {

  # Argument checks
  stopifnot(is(SimulatedBS, "SimulatedBS"))
  size <- as.integer(size)
  stopifnot(size > 0L)
  if (size > 1L) {
    warning(paste0("Only adjacent ", size, "-tuples are created."))
  }

  z <- bplapply(SimulatedBS@z, .makePosAndCounts, size)
  names(z) <- names(SimulatedBS@z)

  # Create MethPat object from pos and counts
  seqnames <- Rle(names(z), sapply(lapply(z, "[[", "pos"), nrow))
  pos <- do.call(rbind, lapply(z, "[[", "pos"))
  counts <- lapply(seq_len(2 ^ size), function(i, z) {
    matrix(unlist(lapply(lapply(z, "[[", "counts"), "[[", i),
                  use.names = FALSE), ncol = 1L)
  }, z = z)
  names(counts) <- MethylationTuples:::.makeMethPatNames(size)
  # TODO: Use sampleName from SimulatedBS object
  counts <- lapply(counts, `colnames<-`, SimulatedBS@SampleName)
  methpat <- MethPat(assays = counts,
                     rowRanges = MTuples(
                       GTuples(seqnames, pos, "*",
                               seqinfo = seqinfo(SimulatedBS)),
                       methinfo = methinfo(SimulatedBS)))
  methpat
}
