### =========================================================================
### SimulateMethylomeParam: An S4 class to store the parameters used to
### simulate a single SimulatedMethylome.
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### list(BSgenome, PartitionedMethylome, MethLevelDT, ComethDT, PatternFreqsDT,
###      SampleName)
###               BSgenome: A BSgenome object for the genome from which to
###                         simulate a methylome.
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
           BSgenome = "BSgenome",
           PartitionedMethylome = "PartitionedMethylome",
           MethLevelDT = "data.table",
           ComethDT = "data.table",
           PatternFreqsDT = "data.table",
           SampleName = "character")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulateMethylomeParam.BSgenome <- function(object) {
  msg <- NULL
  if (!is(object@BSgenome, "BSgenome")) {
    msg <- Biobase::validMsg(msg, paste0("'BSgenome' slot must be a ",
                                         "'BSgenome' object."))
  }
  msg
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
  h <- grep("^h[0-9]+", colnames(object@PatternFreqsDT), value = TRUE)
  if (!identical(colnames(object@PatternFreqsDT), c("type", h, "N"))) {
    msg <- Biobase::validMsg(msg, paste0("Column names of the ",
                                         "'PatternFreqsDT' slot must be, ",
                                         "'type', '",
                                         paste0(h, collapse = "', '"), "', ",
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
  msg <- c(.valid.SimulateMethylomeParam.BSgenome(object),
           .valid.SimulateMethylomeParam.PartitionedMethylome(object),
           .valid.SimulateMethylomeParam.MethLevelDT(object),
           .valid.SimulateMethylomeParam.ComethDT(object),
           .valid.SimulateMethylomeParam.PatternFreqsDT(object),
           .valid.SimulateMethylomeParam.SampleName(object))

  if (is.null(msg)){
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
SimulateMethylomeParam <- function(BSgenome,
                                   PartitionedMethylome,
                                   MethLevelDT,
                                   ComethDT,
                                   PatternFreqsDT,
                                   SampleName) {

  # TODO: Argument checking, e.g., probably a good idea to check length(h) > 2,
  # i.e., that there are a sufficient number of haplotypes in PatternFreqsDT.
  new("SimulateMethylomeParam",
      BSgenome = BSgenome,
      PartitionedMethylome = PartitionedMethylome,
      MethLevelDT = MethLevelDT,
      ComethDT = ComethDT,
      PatternFreqsDT = PatternFreqsDT,
      SampleName = SampleName)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simulate()
###

# TODO: Probably want to export the various comethylation_function options.
#' Simulate a methylome.
#'
#' @note Currently only simulates CpG methylation.
#'
#' @param object A \code{\link{SimulateMethylomeParam}} object.
#' @param nsim Number of samples to simulate using the parameters given in
#' \code{object}.
#' @param seed An object specifying if and how the random number generator
#' should be initialized ('seeded'). For the "MethSimParam" method, either
#' \code{NULL} or an integer that will be used in a call to
#' \code{base::\link[base]{set.seed}} before simulating the samples. If set,
#' the value is saved as the "\code{seed}" attribute of the returned value. The
#' default, \code{NULL}, will not change the random generator state, and return
#' \code{\link{.Random.seed}} as the "\code{seed}" attribute, see 'Value'.
#' @param comethylation_model the type
#' @param BPPARAM An optional
#' \code{BiocParallel::\link[BiocParallel]{BiocParallelParam}} instance
#' determining the parallel back-end to be used during evaluation.
#' @param ... additional arguments passed to the \code{comethylation_function}.
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
                   comethylation_function = methsim:::.sampleComethDT,
                   BPPARAM = bpparam(),
                   epsilon = 0.01, ...) {

            warning("Currently only supports CpG methylation.")
            warning("Currently only supports unstranded methylomes.")

            # Argument checks
            comethylation_function <- match.fun(comethylation_function)
            stopifnot(is.numeric(epsilon) & epsilon > 0 & epsilon < 1)

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

            if (nsim != 1) {
              # Currently only simulate one methylome for a given
              # SimulateMethylomeParam.
              stop("'nsim' must be equal to 1.")
            }

            message("Simulating ", nsim, " methylome...")

            # Sample parameters

            # Sample average methylation levels in each region
            beta_by_region <- .sampleMethLevelDT(
              object@MethLevelDT,
              regionType(object@PartitionedMethylome))
            # Add (subtract) epsilon to zero (one) elements of beta_by_region.
            # Otherwise the entire region will be zero (one).
            beta_by_region[beta_by_region == 1] <- 1 - epsilon
            beta_by_region[beta_by_region == 0] <- epsilon

            # Sample within-fragment co-methylation for each IPD-region_type
            # combination.
            two_tuples <- findMTuples(object@BSgenome, MethInfo("CG"), size = 2)
            # Only want unstranded methylomes
            two_tuples <- unstrand(two_tuples[strand(two_tuples) == "+"])
            # Drop unusable seqlevels
            two_tuples <- keepSeqlevels(two_tuples,
                                        seqlevels(object@PartitionedMethylome))
            lor_by_pair <- .sampleComethDT(two_tuples,
                                           object@ComethDT,
                                           object@PartitionedMethylome,
                                           ...)

            # Methylation loci at which to simulate a methylation state.
            one_tuples <- findMTuples(object@BSgenome, MethInfo("CG"), size = 1)
            # Only want unstranded methylomes.
            one_tuples <- unstrand(one_tuples[strand(one_tuples) == "+"])
            # Drop unusable seqlevels.
            one_tuples <- keepSeqlevels(one_tuples,
                                        seqlevels(object@PartitionedMethylome))
            ol <- findOverlaps(one_tuples, object@PartitionedMethylome)
            beta_by_region <- Rle(beta_by_region, countSubjectHits(ol))

            # Sample haplotype frequencies
            H_by_region <- .samplePatternFreqsDT(
              object@PatternFreqsDT,
              regionType(object@PartitionedMethylome))
            # TODO: "Undo" implicit Rles of H_by_region
            H <- matrix(rep(H_by_region, times = rep(countSubjectHits(ol),
                                                     ncol(H_by_region))),
                        ncol = ncol(H_by_region),
                        dimnames = list(NULL,
                                        paste0('h',
                                               seq_len(ncol(H_by_region)))))

            # Don't generate random numbers in parallel, e.g., via mclapply().
            # It needlessly complicates things and any speed ups are swamped
            # by the running times of other steps in this function.
            u <- lapply(seq_len(ncol(H_by_region)), function (i) {
              runif(length(one_tuples))
            })

            # Simulate Z as a matrix with ncol = ncol(H_by_region). The
            # resulting object is approximately 900 MB in size for a human
            # methylome with 8 "haplotypes".
            # While storing Z as a Matrix::sparseMatrix might seem appealing,
            # the vast majority of entries (~80%) are non-zero (1), therefore
            # the object is not in fact sparse and is even larger in size
            # (2.1 GB). While I could "bit-flip" to make 1 = unmethylated and
            # 0 = methylated, this would only be a source of confusion.
            # In contrast, a DataFrame solution with Rle columns is ~ 240 MB in
            # size. However, row-column access is unacceptably slow for my
            # subsequent application that involves sampling from Z.
            # UP TO HERE: Switch to BiocParallel
            # bpmapply(, SIMPLIFY = TRUE) takes a ridiculously longer time to
            # run than a straightforward mcmapply(, SIMPLIFY = TRUE); why?
            Z <- bplapply(seq_along(u), function(u, beta_by_region,
                                                 lor_by_pair) {
              .simulateZ(as.vector(beta_by_region),
                         lor_by_region,
                         as.vector(seqnames_one_tuples),
                         u)
            }, BPPARAM = BPPARAM)
            Z <- simplify2array(Z)
            # TODO: This is inefficient because it forces a copy of Z.
            colnames(Z) <- colnames(H)

            # UP TO HERE: Create SimulatedMethylome object
            # This should work, I think, but doesn't
            sm <- new("SimulatedMethylome",
                      SummarizedExperiment(assays = SimpleList(Z = Z, H = H),
                                           rowData = one_tuples))


            # Ensure "seed" is set as an attribute of the returned value.
            attr(sm, "seed") <- rng_state
            sm
          }
)

