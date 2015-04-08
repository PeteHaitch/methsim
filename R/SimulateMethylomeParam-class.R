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
### SimulateMethylomeParam: An S4 class to store the parameters used to
### simulate a single SimulatedMethylome.
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### list(BSgenomeName, PartitionedMethylome, MethLevelDT, ComethDT,
###      PatternFreqsDT, SampleName)
###           BSgenomeName: A character vector with the name of the relevant
###                         BSgenome object.
###   PartitionedMethylome: A PartitionedMethylome object.
###            MethLevelDT: A data.table containing the
###                         MethylationTuples::methLevel data.
###               ComethDT: A data.table containing the
###                         MethylationTuples::cometh data.
###         PatternFreqsDT: A data.table containing the
###                         MethylationTuples::PatternFreqs data.
###             SampleName: The sample name.

#' SimulateMethylomeParam class
#'
#' An S4 class for the parameters used by
#' \code{\link{simulate,SimulateMethylome-method}}.
#'
#' @include PartitionedMethylome-class.R
#'
#' @aliases SimulateMethylomeParam
#'
#' @export
setClass("SimulateMethylomeParam",
         slots = list(
           BSgenomeName = "character",
           PartitionedMethylome = "PartitionedMethylome",
           MethLevelDT = "data.table",
           ComethDT = "data.table",
           PatternFreqsDT = "data.table",
           SampleName = "character")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulateMethylomeParam.BSgenomeName <- function(object) {
  msg <- NULL
  if (length(object@BSgenomeName) == 1L && is.character(object@BSgenomeName)) {
    if (!(object@BSgenomeName %in% BSgenome::available.genomes())) {
      msg <- Biobase::validMsg(msg, paste0("'BSgenomeName' must be an element ",
                                           "of ",
                                           "'BSgenome::available.genomes()'."))
    }
  } else {
    msg <- Biobase::validMsg(msg, paste0("'BSgenomeName' must be a length one ",
                                         "character vector."))
  }

}

.valid.SimulateMethylomeParam.PartitionedMethylome <- function(object) {
  msg <- NULL
  if (!is(object@PartitionedMethylome, "PartitionedMethylome")) {
    msg <- Biobase::validMsg(msg, paste0("'PartitionedMethylome' slot must ",
                                         "be a 'PartitionedMethylome'."))
  }
  msg
}

.valid.SimulateMethylomeParam.MethLevelDT <- function(object) {
  msg <- NULL
  if (!is.data.table(object@MethLevelDT)) {
    msg <- Biobase::validMsg(msg, paste0("'MethLevelDT' slot must be a ",
                                         "'data.table'."))
  }
  if (!identical(colnames(object@MethLevelDT), c("type", "beta", "N"))) {
    msg <- Biobase::validMsg(msg, paste0("Column names of the 'MethLevelDT' ",
                                         "slot must be 'type', 'beta' and ",
                                         "'N'."))
  }
  msg
}

.valid.SimulateMethylomeParam.ComethDT <- function(object) {
  msg <- NULL
  if (!is.data.table(object@ComethDT)) {
    msg <- Biobase::validMsg(msg, paste0("'ComethDT' slot must be a ",
                                         "'data.table'."))
  }
  if (!identical(colnames(object@ComethDT), c("IPD", "type", "statistic",
                                              "N"))) {
    msg <- Biobase::validMsg(msg, paste0("Column names of the 'ComethDT' ",
                                         "slot must be 'IPD', 'type', ",
                                         "'statistic' and 'N'." ))
  }
  msg
}

.valid.SimulateMethylomeParam.PatternFreqsDT <- function(object) {
  msg <- NULL
  if (!is.data.table(object@PatternFreqsDT)) {
    msg <- Biobase::validMsg(msg, paste0("'PatternFreqsDT' slot must be a ",
                                         "'data.table'."))
  }
  w <- grep("^w[0-9]+", colnames(object@PatternFreqsDT), value = TRUE)
  if (!identical(colnames(object@PatternFreqsDT), c("type", w, "N"))) {
    msg <- Biobase::validMsg(msg, paste0("Column names of the ",
                                         "'PatternFreqsDT' slot must be, ",
                                         "'type', '",
                                         paste0(w, collapse = "', '"), "', ",
                                         "'statistic' and 'N'." ))
  }
  msg
}

.valid.SimulateMethylomeParam.SampleName <- function(object) {
  msg <- NULL

  if (!is.character(object@SampleName) | length(object@SampleName) != 1L) {
    msg <- Biobase::validMsg(msg, paste0("'SampleName' slot must be a ",
                                         "'character' vector with length 1."))
  }
}

.valid.SimulateMethylomeParam <- function(object) {
  # Include all .valid.SimulateMethylomeParam.* functions in this vector
  msg <- c(.valid.SimulateMethylomeParam.BSgenomeName(object),
           .valid.SimulateMethylomeParam.PartitionedMethylome(object),
           .valid.SimulateMethylomeParam.MethLevelDT(object),
           .valid.SimulateMethylomeParam.ComethDT(object),
           .valid.SimulateMethylomeParam.PatternFreqsDT(object),
           .valid.SimulateMethylomeParam.SampleName(object))

  if (is.null(msg)) {
    return(TRUE)
  } else{
    msg
  }
}

S4Vectors::setValidity2("SimulateMethylomeParam", .valid.SimulateMethylomeParam)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# TODO (long term): This is a barebones constructor. Might want to make
# MethLevelDT, ComethDT and PatternFreqDT formal S4 classes.
#' @export
SimulateMethylomeParam <- function(BSgenomeName,
                                   PartitionedMethylome,
                                   MethLevelDT,
                                   ComethDT,
                                   PatternFreqsDT,
                                   SampleName) {

  # TODO: Argument checking, e.g., probably a good idea to check length(h) > 2,
  # i.e., that there are a sufficient number of possible patterns in
  # PatternFreqsDT.
  new("SimulateMethylomeParam",
      BSgenomeName = BSgenomeName,
      PartitionedMethylome = PartitionedMethylome,
      MethLevelDT = MethLevelDT,
      ComethDT = ComethDT,
      PatternFreqsDT = PatternFreqsDT,
      SampleName = SampleName)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simulate()
###

# TODO: Should comethylation_function and/or epsilon be part of the
# SimulateMethylomeParam object?
# TODO: Should user-messages be suppressible via suppressMessages() or a
# 'verbose' option.
# TODO: Need to document the methsim:::.sampleComethDT parameters (at least
# mean_fun and sd_fun).
#' Simulate a methylome.
#'
#' @note Currently only simulates CpG methylation.
#'
#' @param object A \code{\link{SimulateMethylomeParam}} object.
#' @param nsim The number of samples to simulate using the parameters given in
#' \code{object}.
#' @param seed An object specifying if and how the random number generator
#' should be initialized ('seeded'). For the "MethSimParam" method, either
#' \code{NULL} Or an integer that will be used in a call to
#' \code{base::\link[base]{set.seed}} before simulating the samples. If set,
#' the value is saved as the "\code{seed}" attribute of the returned value. The
#' default, \code{NULL}, will not change the random generator state, and return
#' \code{\link{.Random.seed}} as the "\code{seed}" attribute, see 'Value'.
#' @param BPPARAM An optional
#' \code{BiocParallel::\link[BiocParallel]{BiocParallelParam}} instance
#' determining the parallel back-end to be used during evaluation.
#' @param epsilon An offset added/subtracted to a region's sampled methylation
#' level should it be one/zero values.
#' @param seqlevels A vector of seqlevels at which to simulate a methylome.
#' Will use all available seqlevels if missing.
#' @param comethylation_function A function used to sample from the
#' co-methylation distribution.
#' @param ... Additional arguments passed to the \code{comethylation_function}.
#'
#' @note Currently only support simulation of CpG methylation and unstranded
#' methylomes.
#'
#' @return A SimulatedMethylome, which is the underlying "true" methylome for
#' a single sample.
#'
#' @export
setMethod("simulate",
          "SimulateMethylomeParam",
          function(object,
                   nsim = 1,
                   seed = NULL,
                   BPPARAM = bpparam(),
                   epsilon = 0.01,
                   seqlevels,
                   comethylation_function,
                   ...) {

            # Argument checks
            if (!object@BSgenomeName %in% BSgenome::available.genomes()) {
              stop(paste0("'", object@BSgenomeName, "' package is not ",
                          "available from Bioconductor."))
            }
            if (!requireNamespace(object@BSgenomeName, quietly = TRUE)) {
              stop(paste0("'", object@BSgenomeName, "' package is required.\n",
                          "To install this package, start R and enter:\n",
                          "source('http://bioconductor.org/biocLite.R')\n",
                          "biocLite('", object@BSgenomeName, "')"))

            } else {
              bsgenome <- eval(parse(text = paste0(object@BSgenomeName,
                                                   "::", object@BSgenomeName)))
            }
            # TODO: Is this the best way to set default seqlevels? Can't use
            # seqlevels = seqlevels(object@PartitionedMethylome) in function
            # signature because of 'recursive default argument reference' error.
            if (missing(seqlevels)) {
              seqlevels <- GenomeInfoDb::seqlevels(object@PartitionedMethylome)
            }
            if (!all(seqlevels %in% GenomeInfoDb::seqlevels(bsgenome))) {
              stop(paste0("Unexpected seqlevels.\n",
                          paste0(seqlevels[!seqlevels %in%
                                             GenomeInfoDb::seqlevels(bsgenome)],
                                 collapse = ", "), " are not seqlevels of ",
                          GenomeInfoDb::bsgenomeName(bsgenome)))
            }
            stopifnot(is.numeric(epsilon) & epsilon > 0 & epsilon < 1)
            # TODO: This is a rather clunky way to set a default value of
            # 'comethylation_function'
            if (missing(comethylation_function)) {
              comethylation_function <- .sampleComethDT
            } else {
              comethylation_function <- match.fun(comethylation_function)
            }

            # Non-CpG methylation is unlikely to be implemented.
            message("Currently only simulates CpG methylation.")
            # TODO (long term): Support stranded methylomes.
            message("Currently only simulates unstranded methylomes.")

            # TODO: Will need to revisit how seed is set and (pseudo) random
            # numbers are generated due to the use of BiocParallel and Rcpp*.
            # This chunk for handling RNG generation is based on
            # stats:::simulate.lm.
            if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
              runif(1)
            }
            if (is.null(seed)) {
              rng_state <- get(".Random.seed", envir = .GlobalEnv)
            } else {
              r_seed <- get(".Random.seed", envir = .GlobalEnv)
              set.seed(seed)
              rng_state <- structure(seed, kind = as.list(RNGkind()))
              on.exit(assign(".Random.seed", r_seed, envir = .GlobalEnv))
            }

            if (nsim != 1) {
              # Currently only simulate one methylome for a given
              # SimulateMethylomeParam.
              stop("'nsim' must be equal to 1.")
            }

            # Sample parameters from the SimulateMethylomeParams object.
            message("Simulating ", nsim, " methylome...")

            # Methylation loci at which to simulate a methylation state.
            exclude <- setdiff(GenomeInfoDb::seqlevels(bsgenome), seqlevels)
            one_tuples <- findMTuples(bsgenome,
                                      MethInfo("CG"),
                                      size = 1,
                                      exclude = exclude)
            # Only want unstranded methylomes
            one_tuples <- unstrand(one_tuples[strand(one_tuples) == "+"])
            ol <- findOverlaps(one_tuples, object@PartitionedMethylome)

            # Sample pseudo-haplotype weights
            message("Sampling w...")
            W_by_region <- .samplePatternFreqsDT(
              object@PatternFreqsDT,
              regionType(object@PartitionedMethylome))
            # "Undo" implicit Rles of W_by_region
            W <- matrix(rep(W_by_region, times = rep(countSubjectHits(ol),
                                                     ncol(W_by_region))),
                        ncol = ncol(W_by_region),
                        dimnames = list(NULL,
                                        paste0('W',
                                               seq_len(ncol(W_by_region)))))

            # Sample average methylation levels in each region
            message("Sampling region methylation levels...")
            beta_by_region <- .sampleMethLevelDT(
              object@MethLevelDT,
              regionType(object@PartitionedMethylome))
            # Add (subtract) epsilon to zero (one) elements of beta_by_region.
            # Otherwise the entire region will be zero (one).
            beta_by_region[beta_by_region == 1] <- 1 - epsilon
            beta_by_region[beta_by_region == 0] <- epsilon
            beta_by_region <- Rle(beta_by_region, countSubjectHits(ol))

            # Generate (pseudo) random numbers used by the .simulateZ().
            # Don't generate random numbers in parallel, e.g., via bplapply().
            # It needlessly complicates things (reproducibility of random
            # numbers when generated in parallel is hard) and any speed ups are
            # swamped by the running times of other steps in this function.
            u <- replicate(ncol(W_by_region), runif(length(one_tuples)),
                           simplify = FALSE)

            # Sample within-fragment co-methylation for each IPD-region_type
            # combination.
            message("Sampling LOR...")
            two_tuples <- endoapply(
              split(one_tuples, seqnames(one_tuples)), function(x) {
                n <- length(x)
                MTuples(GTuples(seqnames(x)[seq_len(n - 1)],
                                matrix(c(start(x)[seq.int(1, n - 1)],
                                         start(x)[seq.int(2, n)]), ncol = 2),
                                strand(x)[seq_len(n - 1)],
                                seqinfo = seqinfo(x)),
                        methinfo(x))
              }
            )
            two_tuples <- unlist(two_tuples, use.names = FALSE)
            lor_by_pair <- comethylation_function(two_tuples,
                                                  object@ComethDT,
                                                  object@PartitionedMethylome,
                                                  ...)

            # Simulate Z as a matrix with ncol = ncol(H_by_region). The
            # resulting object is approximately 900 MB in size for a human
            # methylome with 8 pseudo-haplotypes.
            # While storing Z as a Matrix::sparseMatrix might seem appealing,
            # the vast majority of entries (~80%) are non-zero (1), therefore
            # the object is not in fact sparse and is even larger in size
            # (2.1 GB). While I could "bit-flip" to make 1 = unmethylated and
            # 0 = methylated, this would only be a source of confusion.
            # In contrast, a DataFrame solution with Rle columns is ~ 240 MB in
            # size. However, row-column access is unacceptably slow for my
            # subsequent application that involves sampling from Z.
            message("Simulating Z")
            Z <- bplapply(u, function(u,
                                      beta_by_region,
                                      lor_by_pair,
                                      seqnames_one_tuples) {
              .simulateZ(as.vector(beta_by_region),
                         lor_by_pair,
                         as.vector(seqnames_one_tuples),
                         u)
            }, beta_by_region = beta_by_region, lor_by_pair = lor_by_pair,
            seqnames_one_tuples = seqnames(one_tuples), BPPARAM = BPPARAM)
            Z <- simplify2array(Z)
            # TODO: This is inefficient because it forces a copy of Z.
            colnames(Z) <- colnames(W)

            # Create SimulatedMethylome object
            message("Creating SimulatedMethylome object...")
            sm <- new("SimulatedMethylome",
                      SummarizedExperiment(assays = SimpleList(Z = Z, W = W),
                                           rowData = one_tuples))

            # Ensure "seed" is set as an attribute of the returned value.
            attr(sm, "seed") <- rng_state
            sm
          }
)

# TODO: Should comethylation_function and/or epsilon be part of the
# SimulateMethylomeParam object?
# TODO: Should user-messages be suppressible via suppressMessages() or a
# 'verbose' option.
# TODO: Need to document the methsim:::.sampleComethDT parameters (at least
# mean_fun and sd_fun).
#' Simulate a methylome the second method.
#'
#' Rather than using weights (PatternFreqsDT) and simulating the
#' 'true' methylome directly (as a SimulatedMethylome object), simulate2()
#' samples the distribution of average methylation levels (MethLevelDT)
#' and the distributions of within-fragment co-methylation (ComethDT)
#' and computes the resulting transition probabilities under a first-order
#' Markov process.
#'
#' @note Currently only simulates CpG methylation.
#'
#' @param object A \code{\link{SimulateMethylomeParam}} object.
#' @param nsim The number of samples to simulate using the parameters given in
#' \code{object}.
#' @param seed An object specifying if and how the random number generator
#' should be initialized ('seeded'). For the "MethSimParam" method, either
#' \code{NULL} Or an integer that will be used in a call to
#' \code{base::\link[base]{set.seed}} before simulating the samples. If set,
#' the value is saved as the "\code{seed}" attribute of the returned value. The
#' default, \code{NULL}, will not change the random generator state, and return
#' \code{\link{.Random.seed}} as the "\code{seed}" attribute, see 'Value'.
#' @param epsilon An offset added/subtracted to a region's sampled methylation
#' level should it be one/zero values.
#' @param comethylation_function A function used to sample from the
#' co-methylation distribution.
#' @param seqlevels A vector of seqlevels at which to simulate a methylome.
#' Will use all available seqlevels if missing.
#' @param mc_order The order of the Markov chain. Currently only a value of
#' 1 is supported.
#' @param ... Additional arguments passed to the \code{comethylation_function}.
#'
#' @note Currently only support simulation of CpG methylation and unstranded
#' methylomes.
#' @note TODO: Random number generation is not yet properly implemented and so
#' results will likely not be reproducible even if the same seed is used.
#'
#' @return A SimulatedMethylome, which is the underlying "true" methylome for
#' a single sample.
#'
#' @export
setMethod("simulate2",
          "SimulateMethylomeParam",
          function(object,
                   nsim = 1,
                   seed = NULL,
                   epsilon = 0.01,
                   comethylation_function,
                   seqlevels,
                   mc_order = 1L,
                   ...) {


            # Argument checks
            if (!object@BSgenomeName %in% BSgenome::available.genomes()) {
              stop(paste0("'", object@BSgenomeName, "' package is not ",
                          "available from Bioconductor."))
            }
            if (!requireNamespace(object@BSgenomeName, quietly = TRUE)) {
              stop(paste0("'", object@BSgenomeName, "' package is required.\n",
                          "To install this package, start R and enter:\n",
                          "source('http://bioconductor.org/biocLite.R')\n",
                          "biocLite('", object@BSgenomeName, "')"))

            } else {
              bsgenome <- eval(parse(text = paste0(object@BSgenomeName,
                                                   "::", object@BSgenomeName)))
            }
            # TODO: Is this the best way to set default seqlevels? Can't use
            # seqlevels = seqlevels(object@PartitionedMethylome) in function
            # signature because of 'recursive default argument reference' error.
            if (missing(seqlevels)) {
              seqlevels <- GenomeInfoDb::seqlevels(object@PartitionedMethylome)
            }
            if (!all(seqlevels %in% GenomeInfoDb::seqlevels(bsgenome))) {
              stop(paste0("Unexpected seqlevels.\n",
                          paste0(seqlevels[!seqlevels %in%
                                             GenomeInfoDb::seqlevels(bsgenome)],
                                 collapse = ", "), " are not seqlevels of ",
                          GenomeInfoDb::bsgenomeName(bsgenome)))
            }
            stopifnot(is.numeric(epsilon) & epsilon > 0 & epsilon < 1)
            # TODO: This is a rather clunky way to set a default value of
            # 'comethylation_function'
            if (missing(comethylation_function)) {
              comethylation_function <- .sampleComethDT
            } else {
              comethylation_function <- match.fun(comethylation_function)
            }

            # Non-CpG methylation is unlikely to be implemented.
            message("Currently only simulates CpG methylation.")
            # TODO (long term): Support stranded methylomes.
            message("Currently only simulates unstranded methylomes.")

            # TODO: Will need to revisit how seed is set and (pseudo) random
            # numbers are generated due to the use of BiocParallel and Rcpp*.
            # This chunk for handling RNG generation is based on
            # stats:::simulate.lm.
            if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
              runif(1)
            }
            if (is.null(seed)) {
              rng_state <- get(".Random.seed", envir = .GlobalEnv)
            } else {
              r_seed <- get(".Random.seed", envir = .GlobalEnv)
              set.seed(seed)
              rng_state <- structure(seed, kind = as.list(RNGkind()))
              on.exit(assign(".Random.seed", r_seed, envir = .GlobalEnv))
            }

            if (nsim != 1) {
              # Currently only simulate one methylome for a given
              # SimulateMethylomeParam.
              stop("'nsim' must be equal to 1.")
            }

            # Sample parameters from the SimulateMethylomeParams object.
            message("Simulating ", nsim, " methylome...")

            # Methylation loci at which to simulate a methylation state.
            exclude <- setdiff(GenomeInfoDb::seqlevels(bsgenome), seqlevels)
            one_tuples <- findMTuples(bsgenome,
                                      MethInfo("CG"),
                                      size = 1,
                                      exclude = exclude)
            # Only want unstranded methylomes
            one_tuples <- unstrand(one_tuples[strand(one_tuples) == "+"])
            ol <- findOverlaps(one_tuples, object@PartitionedMethylome)

            # Sample average methylation levels in each region
            message("Sampling region methylation levels...")
            beta_by_region <- .sampleMethLevelDT(
              object@MethLevelDT,
              regionType(object@PartitionedMethylome))
            # Add (subtract) epsilon to zero (one) elements of beta_by_region.
            # Otherwise the entire region will be zero (one).
            beta_by_region[beta_by_region == 1] <- 1 - epsilon
            beta_by_region[beta_by_region == 0] <- epsilon
            beta_by_region <- Rle(beta_by_region, countSubjectHits(ol))

            # Sample within-fragment co-methylation for each IPD-region_type
            # combination.
            message("Sampling LOR...")
            two_tuples <- endoapply(
              split(one_tuples, seqnames(one_tuples)), function(x) {
                n <- length(x)
                MTuples(GTuples(seqnames(x)[seq_len(n - 1)],
                                matrix(c(start(x)[seq.int(1, n - 1)],
                                         start(x)[seq.int(2, n)]), ncol = 2),
                                strand(x)[seq_len(n - 1)],
                                seqinfo = seqinfo(x)),
                        methinfo(x))
              }
            )
            two_tuples <- unlist(two_tuples, use.names = FALSE)
            lor_by_pair <- comethylation_function(two_tuples,
                                                  object@ComethDT,
                                                  object@PartitionedMethylome,
                                                  ...)

            # Simulate P.
            message("Simulating P")
            P <- .computeP(as.vector(beta_by_region),
                           lor_by_pair,
                           as.vector(seqnames(one_tuples)),
                           mc_order)

            # Create SimulatedMethylome2 object
            message("Creating SimulatedMethylome2 object...")
            sm2 <- new("SimulatedMethylome2",
                       SummarizedExperiment(assays = SimpleList(P = P),
                                            rowData = one_tuples))

            # Ensure "seed" is set as an attribute of the returned value.
            attr(sm2, "seed") <- rng_state
            sm2
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###
# TODO (long term)
