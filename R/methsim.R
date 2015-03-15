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
