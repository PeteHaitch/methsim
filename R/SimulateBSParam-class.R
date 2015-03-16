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
#' @include SimulateMethylomeParam-class.R
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
#' @param sequencing_type a character specifying either single-end ('SE') or
#' paired-end ('PE'). \strong{NOTE}: Only single-end data is currently
#' supported.
#' @param read_length an integer specifying the read length.
#'
#' @note Currently only simulates whole-genome bisulfite-sequencing data.
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
                   sequencing_type = "SE",
                   read_length = 100L,
                   ...) {

            # TODO (long term): Support stranded data.
            warning("Currently only simulates unstranded data.")
            # TOOD: Add support for paired-end reads.
            warning("Currently only simulates single-end reads.")

            # Argument checks
            stopifnot(sequencing_type == "SE")

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

            message("Simulating ", nsim, " bisulfite-sequencing experiment...")

            # Sample read start sites based on uniform sampling with given
            # average sequencing coverage (aveCov).
            # TODO: This assumes constant read_length; this code will need
            # modification if this assumption is changed.
            n_reads <- trunc(object@aveCov / read_length *
                               seqlengths(object@SimulatedMethylome))
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
            sbs <- bplapply(names(read_start), function(seqname, read_start,
                                                         sm) {

              # TODO: Is it necessary/useful to add seqinfo? All it is likely
              # to do is give a warning if a read runs off the end of the
              # chromosome.
              gr <- GRanges(seqname,
                            IRanges(read_start[[seqname]], width = 100),
                            seqinfo = seqinfo(sm))

              # Only retain reads that overlap at least one methylation locus.
              gr <- subsetByOverlaps(gr, sm)
              ol <- findOverlaps(gr, sm)

              # Sample from H to determine from which each haplotype each read
              # should be sampled.
              # NOTE: Make sure this is only done once per read; to ensure this
              # I only sample a haplotype for the first methylation locus in
              # each read.
              first_hit <- selectHits(ol, "first")
              h <- .sampleH(assay(sm, "H", withDimnames = FALSE)[first_hit, ])

              # For each read, copy the corresponding haplotype.
              z <- .sampleZ(Z = assay(sm, "Z", withDimnames = FALSE),
                            h = h,
                            fh = first_hit,
                            cqh = countQueryHits(ol))

              # Return as a data.table
              data.table("seqnames" = factor(rep(seqname, length(ol)),
                                             levels = seqlevels(sm)),
                         "pos" = start(sm)[subjectHits(ol)],
                         "readID" = queryHits(ol),
                         "z" = unlist(z, use.names = FALSE))
            }, read_start = read_start, sm = object@SimulatedMethylome,
            BPPARAM = BPPARAM)

            sbs <- rbindlist(sbs)
            # Introduce sequencing error + bisulfite-conversion error
            # i.e., flip elements of sbs$z s.t. Prob(flip) = object@errorRate
            .simErrorInPlace(sbs[["z"]], runif(nrow(sbs)), object@errorRate)
            # TODO: Is this the most useful key?
            setkey(sbs, seqnames, readID)

            # Set class of z_dt as SimulatedBS.
            # TODO: By changing the class this forces a copy, I think :(.
            new("SimulatedBS", dt = sbs)
          }
)
