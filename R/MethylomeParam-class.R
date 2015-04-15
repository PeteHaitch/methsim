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
### MethylomeParam: An S4 class to store the parameters used to simulate
### SimulatedMethylome.
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
###         MixtureWeights: A vector of weights representing the mixture of
###                         methylomes in the sample.
###             SampleName: The sample name.

#' MethylomeParam class
#'
#' An S4 class for the parameters used by
#' \code{\link{simulate,MethylomeParam-method}}.
#'
#' @include PartitionedMethylome-class.R
#'
#' @aliases MethylomeParam
#'
#' @export
setClass("MethylomeParam",
         slots = list(
           BSgenomeName = "character",
           PartitionedMethylome = "PartitionedMethylome",
           MethLevelDT = "data.table",
           ComethDT = "data.table",
           MixtureWeights = "numeric",
           SampleName = "character")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.MethylomeParam.BSgenomeName <- function(object) {
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

.valid.MethylomeParam.PartitionedMethylome <- function(object) {
  msg <- NULL
  if (!is(object@PartitionedMethylome, "PartitionedMethylome")) {
    msg <- Biobase::validMsg(msg, paste0("'PartitionedMethylome' slot must ",
                                         "be a 'PartitionedMethylome'."))
  }
  msg
}

.valid.MethylomeParam.MethLevelDT <- function(object) {
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

.valid.MethylomeParam.ComethDT <- function(object) {
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

.valid.MethylomeParam.MixtureWeights <- function(object) {
  msg <- NULL
  if (!is.numeric(object@MixtureWeights)) {
    msg <- Biobase::validMsg(msg, paste0("'MixtureWeights' slot must be a ",
                                         "numeric vector."))
  } else {
    if (any(object@MixtureWeights > 1 || object@MixtureWeights < 0)) {
      msg <- Biobase::validMsg(msg, paste0("'MixtureWeights' slot must all be ",
                                           "between zero and one."))

    }
    if (!all.equal(sum(object@MixtureWeights), 1)) {
      msg <- Biobase::validMsg(msg, paste0("'MixtureWeights' slot must sum to ",
                                           "one."))
    }
  }
  msg
}

.valid.MethylomeParam.SampleName <- function(object) {
  msg <- NULL

  if (!is.character(object@SampleName) | length(object@SampleName) != 1L) {
    msg <- Biobase::validMsg(msg, paste0("'SampleName' slot must be a ",
                                         "'character' vector with length 1."))
  }
}

.valid.MethylomeParam <- function(object) {
  # Include all .valid.MethylomeParam.* functions in this vector
  msg <- c(.valid.MethylomeParam.BSgenomeName(object),
           .valid.MethylomeParam.PartitionedMethylome(object),
           .valid.MethylomeParam.MethLevelDT(object),
           .valid.MethylomeParam.ComethDT(object),
           .valid.MethylomeParam.MixtureWeights(object),
           .valid.MethylomeParam.SampleName(object))

  if (is.null(msg)) {
    return(TRUE)
  } else{
    msg
  }
}

S4Vectors::setValidity2("MethylomeParam", .valid.MethylomeParam)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# TODO (long term): This is a barebones constructor. Might want to make
# MethLevelDT and ComethDT formal S4 classes.
#' @export
MethylomeParam <- function(BSgenomeName,
                           PartitionedMethylome,
                           MethLevelDT,
                           ComethDT,
                           MixtureWeights,
                           SampleName) {

  # TODO: Argument checking
  new("MethylomeParam",
      BSgenomeName = BSgenomeName,
      PartitionedMethylome = PartitionedMethylome,
      MethLevelDT = MethLevelDT,
      ComethDT = ComethDT,
      MixtureWeights = MixtureWeights,
      SampleName = SampleName)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simulate()
###

# A helper function called by simulate,MethylomeParam-method
.simulateMethylomeParam <- function(i, object, epsilon, ol,
                                    comethylation_function, one_tuples,
                                    two_tuples, ...) {

  # Need to simulate a MarginalProb and LOR for all positive
  # object@MixtureWeights
  n_components <- sum(object@MixtureWeights > 0)

  # Sample average methylation levels in each region
  mldt <- object@MethLevelDT
  rt <- regionType(object@PartitionedMethylome)
  marginal_prob_by_region <- replicate(n_components,
                                       .sampleMethLevelDT(mldt, rt),
                                       simplify = "array")
  # Add (resp. subtract) epsilon to zero (resp. one) elements of
  # marginal_prob_by_region, otherwise the entire region will be
  # zero (resp. one).
  marginal_prob_by_region[marginal_prob_by_region == 1] <- 1 - epsilon
  marginal_prob_by_region[marginal_prob_by_region == 0] <- epsilon
  n_loci_per_region <- countSubjectHits(ol)

  # Sample within-fragment co-methylation for each IPD-region_type
  # combination.
  # dots (...) are passed to comethylation_function().
  # For a seqlevel with n 1-tuples, there are (n - 1) elements of
  # lor_by_pair; there is no value for the first 1-tuple because it
  # has no predecessor (hence the need below to add a value of zero
  # for the first 1-tuple on each seqlevel)
  cdt <- object@ComethDT
  pm <- object@PartitionedMethylome
  # TODO: This is the slowest part, and it would be good to run the
  # n_component replications in parallel. Need to pass down BPPARAM in a clever
  # way to respect inheritance so that it doesn't blow up.
  lor_by_pair <- replicate(n_components,
                           comethylation_function(two_tuples,
                                                  cdt,
                                                  pm,
                                                  ...),
                           simplify = "array")
  # Add a value of LOR = 0 for the first methylation locus of each
  # seqlevel.
  cn <- paste0("component_", seq_len(n_components))
  lor <- matrix(0,
                nrow = sum(as.numeric(n_loci_per_region)),
                ncol = n_components,
                dimnames = list(NULL, cn))
  first_loci <- start(seqnames(one_tuples))
  lor[-c(first_loci), ] <- lor_by_pair

  # Create SimulatedMethylome object
  # NOTE: The column names are 'component1', ..., 'componentW',
  # where W = n_components
  marginal_prob <- matrix(rep(marginal_prob_by_region,
                              rep(n_loci_per_region, n_components)),
                          ncol = n_components,
                          dimnames = list(NULL, cn))
  mixture_weights <- S4Vectors::DataFrame(lapply(object@MixtureWeights, Rle,
                                                 lengths = nrow(marginal_prob)))
  colnames(mixture_weights) <- cn
  assays <- S4Vectors::SimpleList(MarginalProb = marginal_prob,
                                  LOR = lor,
                                  MixtureWeights = mixture_weights)
  # UP TO HERE: Need to re-write SimulatedMethylome class before I can test
  # the rest of simulate,MethylomeParam-method.
  new("SimulatedMethylome",
      SummarizedExperiment(assays = assays, rowRanges = one_tuples))
}

# TODO: NULL vs. missing arguments; what's the best choice? (I remember
# reading some advice from Hadley on this issue).
# TODO: Investigate locus-specific average methylation levels rather than
# region-specifc methylation levels.
# TODO: Should comethylation_function and/or epsilon be part of the
# MethylomeParam object? No. They specify how to construct the
# SimulatedMethylome; the same MethylomeParam can be used to create multiple
# SimulatedMethylome objects by different sampling schemes.
# TODO: Should user-messages be suppressible via suppressMessages() or a
# 'verbose' option.
# TODO: Need to document the methsim:::.sampleComethDT parameters (at least
# mean_fun and sd_fun).
#' Simulate a methylome.
#'
#' \code{simulate()} samples the distribution of average methylation level
#' (\code{MethLevelDT} slot) and the distributions of within-fragment
#' co-methylation (\code{ComethDT} slot) to compute the transition
#' probabilities under a first-order Markov process.
#'
#' @param object A \code{\link{MethylomeParam}} object.
#' @param nsim The number of methylomes to simulate using the
#' parameters given in \code{object}.
#' @param seed An object specifying if and how the random number generator
#' should be initialized ('seeded'). For the \code{MethylomeParam}
#' method, either \code{NULL} or an integer that will be used in a call to
#' \code{base::\link[base]{set.seed}} before simulating the samples
#' (methylomes). If set, the value is saved as the "\code{seed}" attribute of
#' the returned value. The default, \code{NULL}, will not change the random
#' generator state, and return \code{\link{.Random.seed}} as the "\code{seed}"
#' attribute, see 'Value'.
#' @param epsilon An offset added/subtracted to a region's sampled methylation
#' level to avoid zero/one values.
#' @param comethylation_function A function used to sample from the
#' co-methylation distribution. If not specified, uses the default, see
#' 'Co-methylation'.
#' @param seqlevels A character vector of
#' \code{GenomeInfoDb::\link[GenomeInfoDb]{seqlevels}} at which to simulate a
#' methylome. If missing, the default is to use all available seqlevels.
#' @param BPPARAM An optional
#' \code{BiocParallel::\link[BiocParallel]{BiocParallelParam}} instance
#' determining the parallel back-end to be used during evaluation.
#' @param ... Additional arguments passed to the \code{comethylation_function}.
#'
#' @section Co-methylation:
#' \strong{TODO}: Describe \code{comethylation_function}, \code{mean_fun},
#' \code{sd_fun}, etc.
#'
#' @note Currently only support simulation of CpG methylation and unstranded
#' methylomes.
#'
#' @return A list of length \code{nsim} of SimulatedMethylome objects.
#'
#' @export
setMethod("simulate",
          "MethylomeParam",
          function(object,
                   nsim = 1,
                   seed = NULL,
                   epsilon = 0.01,
                   comethylation_function,
                   seqlevels,
                   BPPARAM = bpparam(),
                   ...) {

            # Argument checks
            # Load the namespace for the BSgenome
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
            # Check epsilon
            stopifnot(is.numeric(epsilon) & epsilon > 0 & epsilon < 1)
            # TODO: This is a rather clunky way to set a default value of
            # 'comethylation_function'
            if (missing(comethylation_function)) {
              comethylation_function <- .sampleComethDT
            } else {
              comethylation_function <- match.fun(comethylation_function)
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

            # Non-CpG methylation is unlikely to be implemented.
            message("Currently only simulates CpG methylation.")
            # TODO (long term): Support stranded methylomes.
            message("Currently only simulates unstranded methylomes.")

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

            # Simulate nsim methylomes.
            if (nsim >= 2) {
              message("Simulating ", nsim, " methylomes...")
            } else {
              message("Simulating ", nsim, " methylome...")
            }

            # Methylation loci at which to simulate a methylation state.
            exclude <- setdiff(GenomeInfoDb::seqlevels(bsgenome), seqlevels)
            one_tuples <- findMTuples(bsgenome,
                                      MethInfo("CG"),
                                      size = 1,
                                      exclude = exclude)
            # Only want unstranded methylomes
            one_tuples <- unstrand(one_tuples[strand(one_tuples) == "+"])
            ol <- findOverlaps(one_tuples, object@PartitionedMethylome)

            # sufficient_one_tuples are those seqlevels with at least two
            # 1-tuples.
            sufficient_one_tuples <- elementLengths(
              split(one_tuples, seqnames(one_tuples))) > 1L
            two_tuples <- endoapply(
              split(one_tuples, seqnames(one_tuples))[sufficient_one_tuples],
              function(x) {
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

            # Simulate nsim SimulatedMethylome objects in parallel.
            list_of_sm <- bplapply(
              seq_len(nsim),
              .simulateMethylomeParam(i,
                                      object,
                                      epsilon,
                                      ol,
                                      comethylation_function,
                                      one_tuples,
                                      two_tuples
              ), object = object, epsilon = epsilon,
              ol = ol, one_tuples = one_tuples,
              comethylation_function = comethylation_function,
              two_tuples = two_tuples, BPPARAM = BPPARAM)

            # Ensure "seed" is set as an attribute of the returned value.
            attr(list_of_simulated_methylomes, "seed") <- rng_state
            list_of_simulated_methylomes
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###
# TODO (long term)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getters/setters
###
# TODO (long term)
