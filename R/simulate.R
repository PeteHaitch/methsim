# TODO: Document.
#' Simulate a methylome.
#'
#' @param object A \code{\link{MethSimParam}} object.
#' @param nsim Number of samples to simulate using the parameters given in
#' \code{object}.
#' @param seed An object specifying if and how the random number generator
#' should be initialized ('seeded'). For the "MethSimParam" method, either
#' \code{NULL} or an integer that will be used in a call to
#' \code{base::\link[base]{set.seed}} before simulating the samples. If set,
#' the value is saved as the "\code{seed}" attribute of the returned value. The
#' default, \code{NULL}, will not change the random generator state, and return
#' \code{\link{.Random.seed}} as the "\code{seed}" attribute, see 'Value'.
#' @param BPPARAM An optional
#' \code{BiocParallel::\link[BiocParallel]{BiocParallelParam}} instance
#' determining the parallel back-end to be used during evaluation.
#'
#' @return A SimulatedMethylome, which is the underlying "true" methylome for
#' a single sample.
#'
#' @export
setMethod("simulate",
          "MethSimParam",
          function(object, nsim = 1, seed = NULL, BPPARAM = bpparam()) {

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
              # Can only simulate one "true" methylome for a given MethSimParam
              stop("'nsim' must be equal to 1.")
            }

            message("Simulating ", nsim, " samples...")

            # Set the






            methpat <- .simulate(msp = object, BPPARAM)


            # Ensure "seed" is set as an attribute of the returned value.
            attr(methpat, "seed") <- rng_state
            methpat
          }
)
