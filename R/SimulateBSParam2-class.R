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
### SimulateBSParam2: An S4 class to store the parameters used to simulate
### bisulfite-sequencing data from a SimulatedMethylome2 object;
### BS = [WGBS | RRBS | eRRBS]
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### list(SimulatedMethylome2, aveCov, errorRate, target)
### SimulatedMethylome2: A SimulatedMethylome2 object.
###             aveCov: The average sequencing coverage to simulate.
###          errorRate: The error rate (combined sequencing error and
###                     bisulfite-conversion error) to use in the simulation.
###             target: The co-ordinates of baits/enriched regions used in
###                     simulating RRBS/eRRBS data.

# TODO: Uncomment if not defined elsewhere in the package.
# setClassUnion("GRangesORNULL",
#               members = c("GRanges", "NULL")
# )

#' SimulateBSParam2 class
#'
#' An S4 class for the parameters used by
#' \code{\link{simulate2,SimulateBS2Param-method}}.
#'
#' @include SimulatedMethylome2-class.R
#'
#' @aliases SimulateBSParam2
#'
#' @export
setClass("SimulateBSParam2",
         slots = list(
           SimulatedMethylome2 = "SimulatedMethylome2",
           aveCov = "numeric",
           errorRate = "numeric",
           target = "GRangesORNULL")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulateBSParam2.SimulatedMethylome2 <- function(object) {
  msg <- NULL
  if (!is(object@SimulatedMethylome2, "SimulatedMethylome2")) {
    msg <- Biobase::validMsg(msg, paste0("'SimulatedMethylome2' slot must be ",
                                         "a 'SimulatedMethylome2' object."))
  }
  msg
}

.valid.SimulateBSParam2.aveCov <- function(object) {
  msg <- NULL
  if (!is.numeric(object@aveCov) ||
      object@aveCov <= 0) {
    msg <- Biobase::validMsg(msg, paste0("'aveCov' slot must be a postive ",
                                         "number."))
  }
  msg
}

.valid.SimulateBSParam2.errorRate <- function(object) {
  msg <- NULL
  if (!is.numeric(object@errorRate) ||
      object@errorRate < 0 ||
      object@errorRate > 1) {
    msg <- Biobase::validMsg(msg, paste0("'errorRate' slot must be a number ",
                                         "between 0 and 1."))
  }
  msg
}

.valid.SimulateBSParam2.target <- function(object) {
  msg <- NULL
  if (!is(object@target, "GRangesORNULL")) {
    msg <- Biobase::validMsg(msg, paste0("'target' slot must be a GRanges ",
                                         "object or NULL."))
  } else {
    if (is(object@target, "GRanges")) {
      if (!identical(seqinfo(object@target),
                     seqinfo(object@SimulatedMethylome2))) {
        msg <- Biobase::validMsg(msg, paste0("'target' slot and ",
                                             "'SimulatedMethylome2' slot must ",
                                             " have identical 'seqinfo'."))
      }
    }
  }
  msg
}

.valid.SimulateBSParam2 <- function(object) {
  # Include all .valid.SimulateBSParam2.* functions in this vector
  msg <- c(.valid.SimulateBSParam2.SimulatedMethylome2(object),
           .valid.SimulateBSParam2.aveCov(object),
           .valid.SimulateBSParam2.errorRate(object),
           .valid.SimulateBSParam2.target(object))

  if (is.null(msg)) {
    return(TRUE)
  } else{
    msg
  }
}

S4Vectors::setValidity2("SimulateBSParam2", .valid.SimulateBSParam2)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
SimulateBSParam2 <- function(SimulatedMethylome2,
                             aveCov = 30L,
                             errorRate = 0.01,
                             target = NULL) {

  # TODO: Argument checks
  new("SimulateBSParam2",
      SimulatedMethylome2 = SimulatedMethylome2,
      aveCov = aveCov,
      errorRate = errorRate,
      target = target)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlevels()
###

# TODO

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simulate2()
###

# # TODO: Should information about the type of sequencing (e.g., single-end or
# # paired-end, read length distribution, fragment length distribution, etc.) be
# # a part of the SimulateBSParam2 object rather than passed via ...?
# # TODO: Add message() output with timing information if verbose = TRUE.
# #' Simulate a bisulfite-sequencing experiment using the second method.
# #'
# #' @param object a \code{\link{SimulateBSParam2}} object.
# #' @param nsim the number of samples to simulate using the parameters given in
# #' \code{object}.
# #' @param seed an object specifying if and how the random number generator
# #' should be initialized ('seeded'). For the "MethSimParam" method, either
# #' \code{NULL} or an integer that will be used in a call to
# #' \code{base::\link[base]{set.seed}} before simulating the samples. If set,
# #' the value is saved as the "\code{seed}" attribute of the returned value. The
# #' default, \code{NULL}, will not change the random generator state, and return
# #' \code{\link{.Random.seed}} as the "\code{seed}" attribute, see 'Value'.
# #' @param BPPARAM an optional
# #' \code{BiocParallel::\link[BiocParallel]{BiocParallelParam}} instance
# #' determining the parallel back-end to be used during evaluation.
# #' @param sequencingType a character specifying either single-end ('SE') or
# #' paired-end ('PE'). \strong{NOTE}: Only single-end data is currently
# #' supported.
# #' @param seqlevels A vector of seqlevels at which to bisulfite-sequencing
# #' reads. Will use all available seqlevels if missing, where availability is
# #' defined by the seqlevels "in use" in the SimulatedMethylome slot.
# #' @param readLength an integer specifying the read length.
# #'
# #' @note Currently only simulates whole-genome bisulfite-sequencing data.
# #' \strong{WARNING}: Currently reads are simulated for circular seqlevels
# #' such as 'chrM' (mitochondrial DNA).
# #'
# #' @return a SimulatedBS2 object.
# #'
# #' @export
# setMethod("simulate",
#           "SimulateBSParam2",
#           function(object,
#                    nsim = 1,
#                    seed = NULL,
#                    BPPARAM = bpparam(),
#                    sequencingType = "SE",
#                    readLength = 100L,
#                    seqlevels,
#                    ...) {
#
#             # Argument checks
#             # TODO: Is this the best way to set default seqlevels? Can't use
#             # seqlevels = seqlevels(object@PartitionedMethylome) in function
#             # signature because of 'recursive default argument reference' error.
#             # TODO: Propose a seqlevelsInUse,SummarizedExperiment-method.
#             valid_seqlevels <- GenomeInfoDb::seqlevelsInUse(
#               rowRanges(object@SimulatedMethylome2))
#             if (missing(seqlevels)) {
#               # Only use seqlevels that are "active" in the SimulatedMethylome2
#               # object
#               seqlevels <- valid_seqlevels
#             } else {
#               # Check that supplied seqlevels are valid
#               if (!all(seqlevels %in% valid_seqlevels)) {
#               stop(paste0("Unexpected seqlevels.\n",
#                           paste0(seqlevels[!seqlevels %in% valid_seqlevels],
#                                  collapse = ", "), " are not seqlevels of ",
#                           "'SimulateBSParam2'."))
#               }
#             }
#             # Only single-end sequencing currently supported
#             stopifnot(sequencingType == "SE")
#
#             # TODO: Will need to revisit how seed is set and (pseudo) random
#             # numbers are generated due to the use of BiocParallel and Rcpp*.
#             # This chunk for handling RNG generation is based on
#             # stats:::simulate.lm.
#             if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
#               runif(1)
#             }
#             if (is.null(seed)) {
#               rng_state <- get(".Random.seed", envir = .GlobalEnv)
#             } else {
#               r_seed <- get(".Random.seed", envir = .GlobalEnv)
#               set.seed(seed)
#               rng_state <- structure(seed, kind = as.list(RNGkind()))
#               on.exit(assign(".Random.seed", r_seed, envir = .GlobalEnv))
#             }
#
#             # TODO: Remove this restriction.
#             if (nsim != 1) {
#               # Currently only simulate one bisulfite-sequencing experiment for
#               # a given SimulateBSParam.
#               stop("'nsim' must be equal to 1.")
#             }
#
#             # TODO (long term): Support stranded data.
#             warning("Currently only simulates unstranded data.")
#             # TOOD: Add support for paired-end reads.
#             warning("Currently only simulates single-end reads.")
#
#             message("Simulating ", nsim, " bisulfite-sequencing sample...")
#
#             # Circular chromosomes are hard. While the read automatically
#             # gets wrapped around by GRanges(), it makes subsequent functions,
#             # e.g., asMethPat(), more complicated. So, for now, I do not allow
#             # simulation of reads for circular chromosomes.
#             is_circ <- isCircular(
#               seqinfo(object@SimulatedMethylome2))[seqlevels]
#             if (any(is_circ)) {
#               warning(paste0("No reads simulated for ",
#                              paste0(seqlevels[is_circ], collapse = ", "),
#                              " (circular chromosomes not yet supported)."))
#             }
#             # Update seqlevels
#             seqlevels <- seqlevels[!is_circ]
#             if (length(seqlevels) == 0L) {
#               # TODO: Return the "empty" object instead of a warning message.
#               stop("No reads simulated")
#             }
#             # Sample read start sites based on uniform sampling with given
#             # average sequencing coverage (aveCov).
#             # TODO: This assumes constant readLength; this code will need
#             # modification if this assumption is changed.
#             n_reads <- as.list(trunc(
#               object@aveCov / readLength *
#                 seqlengths(object@SimulatedMethylome2)[seqlevels]))
#
#             # TODO: Perhaps the number of reads per-chromosome should be
#             # sampled from a multinomial(sum(n_reads), n_reads)?
#             # Don't simulate read_start in parallel, e.g., via bpmapply().
#             # It needlessly complicates things (reproducibility of random
#             # numbers when generated in parallel is hard) and any speed ups are
#             # swamped by the running times of other steps in this function.
#             read_start <- mapply(function(n, seqlength) {
#               .sampleReadStart(n, seqlength)
#             }, n = n_reads,
#             seqlength = seqlengths(object@SimulatedMethylome2)[seqlevels],
#             SIMPLIFY = FALSE)
#             # Sorting makes things easier to keep track of, and sorting an
#             # integer vector is easier than sorting more complicated objects,
#             # e.g., GRanges objects.
#             read_start <- bplapply(read_start, sort, BPPARAM = BPPARAM)
#
#             # Compute the transition probabilities.
#             # UP TO HERE: Re-write .computeP
#             P <- .computeP(assay(object@SimulatedMethylome2, "marginalProb",
#                                  withDimnames = FALSE),
#                            assay(object@SimulatedMethylome2, "LOR",
#                                  withDimnames = FALSE),
#                            mc_order)
#
#             # Find reads that overlap methylation loci and then sample a
#             # methylation pattern for each such read.
#             # TODO: Take care if simulate() itself is being run in parallel
#             # (or at least document that it could spawn heaps of processes).
#             # TODO: Switch to bpmapply to avoid passing the entire read_start
#             # list to each parallel process.
#             z <- bplapply(names(read_start), function(seqname, read_start,
#                                                       readLength, sm) {
#
#               # TODO: This may cause warnings (at least when this isn't run
#               # in parallel, which causes warning()s to be suppressed). These
#               # warnings will occur if a read runs "off the end" of the
#               # seqlevel.
#               gr <- GRanges(seqname,
#                             IRanges(read_start[[seqname]], width = readLength),
#                             seqinfo = seqinfo(sm))
#
#               # Only retain reads that overlap at least one methylation locus.
#               gr <- subsetByOverlaps(gr, sm)
#               ol <- findOverlaps(gr, sm)
#
#               # Sample from beta_by_region to determine the probability that
#               # the first methylation locus in the read is methylated.
#               first_hit <- selectHits(ol, "first")
#
#               #  to determine from which each haplotype each read
#               # should be sampled.
#               # NOTE: Make sure this is only done once per read; to ensure this
#               # I only sample a haplotype for the first methylation locus in
#               # each read.
#               first_hit <- selectHits(ol, "first")
#
#               # UP TO HERE
#
#             }, read_start = read_start, readLength = readLength,
#             sm = object@SimulatedMethylome, BPPARAM = BPPARAM)
#
#
#
#
#
#
#
#
#             # Remove reads that don't overlap any methylation loci (and sort
#             # those that remain).
#             read_start <- bpmapply(function(rs, seqname, row_ranges_sm) {
#               # TODO: This may cause warnings (at least when this isn't run
#               # in parallel, which causes warning()s to be suppressed). These
#               # warnings will occur if a read runs "off the end" of the
#               # seqlevel.
#               gr <- GRanges(seqname,
#                             IRanges(rs, width = readLength),
#                             seqinfo = seqinfo(row_ranges_sm))
#               gr <- subsetByOverlaps(gr, row_ranges_sm)
#               sort(gr)
#             }, rs = read_start,
#             seqname = names(read_start),
#             MoreArgs =
#               list(row_ranges_sm = rowRanges(object@SimulatedMethylome2)),
#             BPPARAM = BPPARAM)
