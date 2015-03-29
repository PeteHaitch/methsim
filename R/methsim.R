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

#' Simulate DNA methylation data.
#'
#' \pkg{methsim} is software to simulate DNA methylation sequencing data.
#' \code{methsim} can currently simulate data from bisulfite-sequencing assays,
#' such as methylC-seq or BS-seq. It uses a non-stationary, inhomogeneous Markov
#' model. This model can incorporate the strong spatial dependence of DNA
#' methylation. Simulation parameters can be estimated from data or specified
#' by the user.
#'
#' Please refer to the vignettes to see how to use the \pkg{methsim} package.
#'
#' @docType package
#' @name methsim
#' @useDynLib methsim, .registration = TRUE
#' @import methods
#' @import GenomicRanges
#' @import Rcpp
#' @import MethylationTuples
#' @import MethylSeekR
#' @import data.table
#' @import BiocParallel
#' @importFrom Biobase validMsg
#' @importFrom mhsmm hmmspec hmmfit
NULL
