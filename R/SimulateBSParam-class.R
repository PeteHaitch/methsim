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
### SimulateBSParam: An S4 class to store the parameters used to simulate
### bisulfite-sequencing data from a SimulatedMethylome object;
### BS = [WGBS | RRBS | eRRBS]
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### list(SimulatedMethylome, aveCov, errorRate, target)
### SimulatedMethylome: A SimulatedMethylome object.
###             aveCov: The average sequencing coverage to simulate.
###          errorRate: The error rate (combined sequencing error and
###                     bisulfite-conversion error) to use in the simulation.
###             target: The co-ordinates of baits/enriched regions used in
###                     simulating RRBS/eRRBS data.

setClassUnion("GRangesORNULL",
              members = c("GRanges", "NULL")
)

#' SimulateBSParam class
#'
#' An S4 class for the parameters used by
#' \code{\link{simulate,SimulateBS-method}}.
#'
#' @include SimulatedMethylome-class.R
#'
#' @aliases SimulateBSParam
#'
#' @export
setClass("SimulateBSParam",
         slots = list(
           SimulatedMethylome = "SimulatedMethylome",
           aveCov = "numeric",
           errorRate = "numeric",
           target = "GRangesORNULL")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulateBSParam.SimulatedMethylome <- function(object) {
  msg <- NULL
  if (!is(object@SimulatedMethylome, "SimulatedMethylome")) {
    msg <- Biobase::validMsg(msg, paste0("'SimulatedMethylome' slot must be a ",
                                         "'SimulatedMethylome' object."))
  }
  msg
}

.valid.SimulateBSParam.aveCov <- function(object) {
  msg <- NULL
  if (!is.numeric(object@aveCov) ||
      object@aveCov <= 0) {
    msg <- Biobase::validMsg(msg, paste0("'aveCov' slot must be a postive ",
                                         "number."))
  }
  msg
}

.valid.SimulateBSParam.errorRate <- function(object) {
  msg <- NULL
  if (!is.numeric(object@errorRate) ||
      object@errorRate < 0 ||
      object@errorRate > 1) {
    msg <- Biobase::validMsg(msg, paste0("'errorRate' slot must be a number ",
                                         "between 0 and 1."))
  }
  msg
}

.valid.SimulateBSParam.target <- function(object) {
  msg <- NULL
  if (!is(object@target, "GRangesORNULL")) {
    msg <- Biobase::validMsg(msg, paste0("'target' slot must be a GRanges ",
                                         "object or NULL."))
  } else {
    if (is(object@target, "GRanges")) {
      if (!identical(seqinfo(object@target),
                     seqinfo(object@SimulatedMethylome))) {
        msg <- Biobase::validMsg(msg, paste0("'target' slot and ",
                                             "'SimulatedMethylome' slot must ",
                                             " have identical 'seqinfo'."))
      }
    }
  }
  msg
}

.valid.SimulateBSParam <- function(object) {
  # Include all .valid.SimulateBSParam.* functions in this vector
  msg <- c(.valid.SimulateBSParam.SimulatedMethylome(object),
           .valid.SimulateBSParam.aveCov(object),
           .valid.SimulateBSParam.errorRate(object),
           .valid.SimulateBSParam.target(object))

  if (is.null(msg)){
    return(TRUE)
  } else{
    msg
  }
}

S4Vectors::setValidity2("SimulateBSParam", .valid.SimulateBSParam)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
SimulateBSParam <- function(SimulatedMethylome,
                            aveCov = 30L,
                            errorRate = 0.01,
                            target = NULL) {

  # TODO: Argument checks
  new("SimulateBSParam",
      SimulatedMethylome = SimulatedMethylome,
      aveCov = aveCov,
      errorRate = errorRate,
      target = target)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simulate()
###

# TODO: Should information about the type of sequencing (e.g., single-end or
# paired-end, read length distribution, fragment length distribution, etc.) be
# a part of the SimulateBSParam object rather than passed via ...?
# TODO: Add message() output with timing information if verbose = TRUE.
#' Simulate a bisulfite-sequencing experiment.
#'
#' @param object a \code{\link{SimulateBSParam}} object.
#' @param nsim the number of samples to simulate using the parameters given in
#' \code{object}.
#' @param seed an object specifying if and how the random number generator
#' should be initialized ('seeded'). For the "MethSimParam" method, either
#' \code{NULL} or an integer that will be used in a call to
#' \code{base::\link[base]{set.seed}} before simulating the samples. If set,
#' the value is saved as the "\code{seed}" attribute of the returned value. The
#' default, \code{NULL}, will not change the random generator state, and return
#' \code{\link{.Random.seed}} as the "\code{seed}" attribute, see 'Value'.
#' @param BPPARAM an optional
#' \code{BiocParallel::\link[BiocParallel]{BiocParallelParam}} instance
#' determining the parallel back-end to be used during evaluation.
#' @param sequencingType a character specifying either single-end ('SE') or
#' paired-end ('PE'). \strong{NOTE}: Only single-end data is currently
#' supported.
#' @param readLength an integer specifying the read length.
#'
#' @note Currently only simulates whole-genome bisulfite-sequencing data.
#' \strong{WARNING}: Currently reads are simulated for circular seqlevels
#' such as 'chrM' (mitochondrial DNA).
#'
#' @return a SimulatedBS object.
#'
#' @export
setMethod("simulate",
          "SimulateBSParam",
          function(object,
                   nsim = 1,
                   seed = NULL,
                   BPPARAM = bpparam(),
                   sequencingType = "SE",
                   readLength = 100L,
                   ...) {

            # TODO (long term): Support stranded data.
            warning("Currently only simulates unstranded data.")
            # TOOD: Add support for paired-end reads.
            warning("Currently only simulates single-end reads.")

            # Argument checks
            stopifnot(sequencingType == "SE")

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
              on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
            }

            # TODO: Remove this restriction.
            if (nsim != 1) {
              # Currently only simulate one bisulfite-sequencing experiment for
              # a given SimulateBSParam.
              stop("'nsim' must be equal to 1.")
            }

            message("Simulating ", nsim, " bisulfite-sequencing sample...")

            # Sample read start sites based on uniform sampling with given
            # average sequencing coverage (aveCov).
            # TODO: This assumes constant readLength; this code will need
            # modification if this assumption is changed.
            n_reads <- trunc(object@aveCov / readLength *
                               seqlengths(object@SimulatedMethylome))
            # TODO: Perhaps the number of reads per-chromosome should be
            # sampled from a multinomial(sum(n_reads), n_reads)?
            # Don't simulate read_start in parallel, e.g., via bpmapply().
            # It needlessly complicates things (reproducibility of random
            # numbers when generated in parallel is hard) and any speed ups are
            # swamped by the running times of other steps in this function.
            read_start <- mapply(function(n, seqlength) {
              .sampleReadStart(n, seqlength)
            }, n = n_reads, seqlength = seqlengths(object@SimulatedMethylome))
            read_start <- bplapply(read_start, sort, BPPARAM = BPPARAM)

            # Find reads that overlap methylation loci (and then sample a
            # methylation pattern for each such read). It's 2-3x faster to
            # subset by chromosome and use findOverlaps on the resulting
            # IRanges than it is to run findOverlaps on the GRanges (but this
            # the subsetting destroys the common index).
            # TODO: Take care if simulate() itself is being run in parallel
            # (or at least document that it could spawn heaps of processes).
            z <- bplapply(names(read_start), function(seqname, read_start,
                                                         readLength, sm) {

              # Circular chromosomes are hard. While the read automatically
              # gets wrapped around, it makes subsequent functions, e.g.,
              # asMethPat(), more complicated. So, for now, I do not allow
              # simulation of reads for circular chromosomes.
              if (isCircular(seqinfo(sm))[seqname]) {
                # TODO: This warning isn't displayed when run via bplapply()
                # (but is displayed when using lapply()).
                warning(paste0("No reads simulated for ", seqname,
                               " (circular chromosomes not yet supported)."))
                return(data.table("pos" = integer(0),
                                  "readID" = integer(0),
                                  "z" = integer(0)))
              }

              # TODO: This may cause warnings (at least when this isn't run
              # in parallel, which causes warning()s to be suppressed). These
              # warnings will occur if a read runs "off the end" of the
              # seqlevel.
              gr <- GRanges(seqname,
                            IRanges(read_start[[seqname]], width = readLength),
                            seqinfo = seqinfo(sm))

              # Only retain reads that overlap at least one methylation locus.
              gr <- subsetByOverlaps(gr, sm)
              ol <- findOverlaps(gr, sm)

              # Sample from W to determine from which each haplotype each read
              # should be sampled.
              # NOTE: Make sure this is only done once per read; to ensure this
              # I only sample a haplotype for the first methylation locus in
              # each read.
              first_hit <- selectHits(ol, "first")
              sampled_W <- .sampleW(assay(sm, "W",
                                          withDimnames = FALSE)[first_hit, ])

              # For each read, copy the corresponding haplotype.
              z <- .sampleZ(Z = assay(sm, "Z", withDimnames = FALSE),
                            sampled_W = sampled_W,
                            fh = first_hit,
                            cqh = countQueryHits(ol))

              # Return as a data.table
              data.table("pos" = start(sm)[subjectHits(ol)],
                         "readID" = queryHits(ol),
                         "z" = unlist(z, use.names = FALSE))
            }, read_start = read_start, readLength = readLength,
            sm = object@SimulatedMethylome, BPPARAM = BPPARAM)

            # Don't rbindlist(z). Instead, keeping as list will
            # actually save memory (no need to retain seqnames for every row)
            # and allow easier parallelisation by seqlevel.
            # Ensure seqlevels are set as names(z).
            names(z) <- names(read_start)

            # Introduce sequencing error + bisulfite-conversion error i.e.,
            # flip elements of z[[i]]$z s.t. Prob(flip) = object@errorRate.
            # Don't simulate errors in parallel, e.g., via bplapply().
            # It needlessly complicates things (reproducibility of random
            # numbers when generated in parallel is hard) and any speed ups are
            # swamped by the running times of other steps in this function.
            lapply(z, function(z_, errorRate) {
              .simErrorInPlace(z_[["z"]], runif(nrow(z_)), errorRate)
            }, errorRate = object@errorRate)

            # Construct SimulatedBS object.
            new("SimulatedBS",
                z = z,
                seqinfo = seqinfo(object@SimulatedMethylome),
                methinfo = methinfo(object@SimulatedMethylome))
          }
)
