### =========================================================================
### Utility functions. These functions are not exported.
### -------------------------------------------------------------------------
###

# A copy of MethylSeekR::plotAlphaDistributionOneChr that can be safely run in
# parallel because it does not create the plot but returns it as an object.
# TODO: Document and robust-ify (long-term)
plotAlphaDistributionOneChr <- function(m, chr.sel, nCGbin = 101,
                                        plot = FALSE) {
  message("determining alpha distribution for chromosome: ",
          chr.sel)
  indx <- as.character(seqnames(m)) == chr.sel
  if (sum(indx) < nCGbin)
    stop(sprintf("Error: less than %d covered CpGs on chromosome %s",
                 nCGbin, chr.sel))
  T <- as.numeric(values(m[indx])[, 1])
  M <- as.numeric(values(m[indx])[, 2])
  score <- MethylSeekR:::calculateAlphaDistr(M, T, nCGbin, num.cores = 1L)
  hist(score, probability = TRUE, breaks = 30,
       xlab = sprintf("posterior mean of alpha (%s)", chr.sel), plot)
}

# A copy of MethylSeekR::plotAlphaDistributionOneChr but does some
# parallelisation by sample.
# TODO: Deprecate
plotAlphaDistributionOneChr_2 <- function(l_msrgr, sn, chr.sel, ncol, nrow,
                                          nCGbin = 101) {
  message("determining alpha distribution for chromosome: ",
          chr.sel)
  l_indx <- bplapply(l_msrgr, function(msrgr, chr.sel) {
    indx <- seqnames(msrgr) == chr.sel
    if (sum(indx) < nCGbin) {
      stop(sprintf("Error: less than %d covered CpGs on chromosome %s",
                   nCGbin, chr.sel))
    }
    indx
  }, chr.sel = chr.sel)
  l_score <- bpmapply(function(m, indx) {
    T <- as.numeric(values(m[indx])[, 1])
    M <- as.numeric(values(m[indx])[, 2])
    MethylSeekR:::calculateAlphaDistr(M, T, nCGbin, num.cores = 1L)
  }, m = l_msrgr, indx = l_indx)
  par(mfrow = c(nrow, ncol))
  mapply(function(score, sn, chr.sel) {
    hist(score, probability = TRUE, breaks = 30,
         xlab = sprintf("posterior mean of alpha (%s)", chr.sel), main = sn)
  }, score = l_score, sn = sn, MoreArgs = list(chr.sel = chr.sel))
  return(NULL)
}

# A copy of MethylSeekR::segmentPMDs that can be safely run in
# parallel because it does not create the plot but returns it as an object.
# TODO: Document and robust-ify (long-term)
segmentPMDs <- function(m, chr.sel, seqLengths, nCGbin = 101, plot = FALSE) {
  hmm.model <- methsim:::trainPMDHMM(m, chr.sel, nCGbin, plot = FALSE)
  y.list <- MethylSeekR:::PMDviterbiSegmentation(m, hmm.model$startval,
                                                 nCGbin, num.cores = 1L)
  segments <- MethylSeekR:::createGRangesObjectPMDSegmentation(m, y.list,
                                                               num.cores = 1L,
                                                               seqLengths)
  list(segments = segments, hist = hmm.model$hist, x = hmm.model$x,
       lines1 = hmm.model$lines1, lines2 = hmm.model$lines2)
}

# A copy of MethylSeekR::trainPMDHMM that can be safely run in
# parallel because it does not create the plot but returns it as an object.
# TODO: Document and robust-ify (long-term)
trainPMDHMM <- function(m, chr.sel, nCGbin, plot = FALSE) {
  message("training PMD-HMM on chromosome ", chr.sel)
  indx <- as.character(seqnames(m)) == chr.sel
  if (sum(indx) < nCGbin) {
    stop(sprintf("Error: less than %d covered CpGs on chromosome %s",
                 nCGbin, chr.sel))
  }
  T <- as.numeric(values(m[indx])[, 1])
  M <- as.numeric(values(m[indx])[, 2])
  score <- MethylSeekR:::calculateAlphaDistr(M, T, nCGbin, num.cores = 1L)
  J = 2
  init0 <- c(0, 1)
  P0 <- t(matrix(c(0.998297563, 0.001702437, 0.002393931, 0.997606069),
                 nrow = J, ncol = J))
  b0 <- list(mu = c(0.3867895, 1.1690474), sigma = c(0.01649962,
                                                     0.1437864))
  startval <- mhsmm::hmmspec(init = init0, trans = P0, parms.emission = b0,
                             dens.emission = dnorm.hsmm)
  train <- list(x = score, N = length(score))
  startval <- mhsmm::hmmfit(train, startval, mstep = mstep.norm)$model
  x <- seq(0, 3, by = 0.01)
  hist <- hist(score, probability = TRUE, breaks = 30,
               xlab = sprintf("posterior mean of alpha (%s)", chr.sel),
               plot = plot)
  lines1 <- dnorm(x, mean = startval$parms.emission$mu[1],
                  sd = sqrt(startval$parms.emission$sigma[1]))
  lines2 <- dnorm(x, mean = startval$parms.emission$mu[2],
                  sd = sqrt(startval$parms.emission$sigma[2]))
  if (plot) {
    hist
    lines(x, lines1, type = "l", col = "red")
    lines(x, lines2, type = "l", col = "green")
    return(list(startval = startval, hist = hist, x = x, lines1 = lines1,
           lines2 = lines2))
  } else {
    list(startval = startval, hist = hist, x = x, lines1 = lines1,
         lines2 = lines2)
  }
}
