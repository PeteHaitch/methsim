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
### PartitionedMethylome: An S4 class to formalise an informal GRanges-based
### class used by MethylSeekR.
### -------------------------------------------------------------------------
###

#' PartitionedMethylome class.
#' @aliases PartitionedMethylome
#'
#' @export
setClass("PartitionedMethylome",
         contains = "GRanges",
         slots = list(
           regionType = "factor"
         ),
         prototype = prototype(
           regionType = factor(levels = c("UMR", "LMR", "PMR", "other"))
         )
)

setMethod(GenomicRanges:::extraColumnSlotNames,
          "PartitionedMethylome",
          function(x) {
            c("regionType")
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.PartitionedMethylome.regionType <- function(object) {
  msg <- NULL
  if (!is.factor(object@regionType) ||
      !identical(levels(object@regionType),
                 c("UMR", "LMR", "PMR", "other"))) {
    msg <- validMsg(msg, paste0("'regionType' must be a factor with levels ",
                                "'UMR', 'LMR', 'PMR' and 'other'."))
  }
  msg
}

.valid.PartitionedMethylome.GRanges <- function(object) {
  msg <- NULL
  if (!isDisjoint(object)) {
    msg <- validMsg(msg, "'PartitionedMethylome' must be disjoint.")
  }
  gaps <- gaps(granges(object))
  gaps <- gaps[strand(gaps) == "*"]
  if (length(gaps)) {
    msg <- validMsg(msg, "'PartitionedMethylome' must not contain gaps.")
  }
  msg
}

.valid.PartitionedMethylome <- function(object) {
  # Include all .valid.PartitionedMethylome.* functions in this vector
  msg <- c(.valid.PartitionedMethylome.regionType(object),
           .valid.PartitionedMethylome.GRanges(object))

  if (is.null(msg)){
    return(TRUE)
  } else{
    msg
  }
}

S4Vectors::setValidity2("PartitionedMethylome", .valid.PartitionedMethylome)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# TODO: Document
#' @export
PartitionedMethylome <- function(seqnames = Rle(),
                                 ranges = IRanges(),
                                 strand = Rle("*", length(seqnames)),
                                 regionType = factor(levels = c("UMR", "LMR",
                                                                "PMR",
                                                                "other")),
                                 ...,
                                 seqlengths = NULL,
                                 seqinfo = NULL) {
  if (!identical(levels(regionType), c("UMR", "LMR", "PMR", "other"))) {
    stop(paste0("'regionType' must be a factor with levels 'UMR', 'LMR', ",
                "'PMR' and 'other'."))
  }

  # Create GRanges
  gr <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand,
                seqlengths = seqlengths, seqinfo = seqinfo, ...)

  new("PartitionedMethylome", gr, regionType = regionType)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @include AllGenerics.R
#' @export
setMethod("regionType",
          "PartitionedMethylome",
          function(object) {
            object@regionType
          }
)
