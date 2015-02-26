# TODO: Create generic and set method for MethSimParam objects.

simulate <- function(param, nsim, seed = NULL, ...) {

  # This chunk for handling RNG generation is based on stats:::simulate.lm.
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
  # TODO: Call the function(s) that do the actual simulation
  # PSEUDO CODE
  val <- for (i in seq_len(nsim)) {
    .simulate(param, ...)
  }
  val
}
