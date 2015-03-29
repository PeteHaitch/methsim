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

### -------------------------------------------------------------------------
### simulate
###
### stats::simulate is an S3 generic.
###

### ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
###   Generic
###     simulate(object, nsim, seed, ...)
###
###   Methods
###     (1) simulate,SimulateMethylomeParam-method
###     (2) simulate,SimulateBSParam-method

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Process
###
### (1) Simulate a methylome with simulate,SimulatedMethylomeParam-method.
### (2) Simulate an assay with simulate,SimulatedBSParam-method.
### (3) Analyse simulated data with as(SimulatedBS, MethPat, size) and
###     methods for working with MethylationTuples::MethPat objects.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Long term (low priority)
###
### - Add methods for simulating methylation microarray data, e.g., Illumina
###   450k.
### - Add methods for creating a FASTQ or BAM file from a SimulatedBS object.

#' Simulate responses
#'
#' @description Simulate data according to given parameters.
#'
#' NOTE: This man page is for the \code{simulate} \emph{S4 generic function}
#' defined in the \pkg{methsim} package. See
#' \code{?statistics::\link[stats]{simulate}} for the default method
#' (defined in the \pkg{statistics} package). Other packages can define
#' specific S4 methods for objects not supported by the default method.
#'
#' @param object an object representing the simulation parameters.
#' @param nsim number of responses to simulate. Defaults to 1.
#' @param seed an object specifying if and how the random number generator
#' should be initialized ('seeded'). Either \code{NULL} or an integer that will
#' be used in a call to \code{base::\link[base]{set.seed}} before simulating
#' the responses. If set, the value is saved as the "\code{seed}" attribute of
#' the returned value. The default, \code{NULL} will not change the random
#' generator state, and return \code{base::\link[base]{.Random.seed}} as the
#' "\code{seed}" attribute, see 'Value'.
#' @param ... additional optional arguments.
#'
#' @return See \code{?stats::\link[stats]{simulate}} for the value returned by
#' the default method.
#'
#' Specific methods defined in other packages should behave as consistently as
#' possible with the default method.
#'
#' @seealso \itemize{
#'  \item \code{stats::\link[stats]{simulate}} for the default \code{simulate}
#'    method.
#'  \item \code{\link[methods]{showMethods}} for displaying a summary of the
#'    methods defined for a given generic function.
#'  \item \code{\link[methods]{selectMethod}} for getting the definition of
#'    a specific method.
#'  \item \code{\link{simulate,SimulateMethylomeParam-method}} for an example
#'  of a specific \code{simulate} method (defined for
#'  \code{\link{SimulateMethylomeParam}} objects).
#'  \item \link[BiocGenerics]{BiocGenerics} for a summary of all the generics
#'  defined in the \pkg{BiocGenerics} package.
#' }
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#'
#' @keywords methods
#'
#' @export
setGeneric("simulate", function(object, nsim = 1, seed = NULL, ...) {
  standardGeneric("simulate")
})

### -------------------------------------------------------------------------
### regionType
###

#' @export
setGeneric("regionType", function(object) {
  standardGeneric("regionType")
})
