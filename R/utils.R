### =========================================================================
### Utility functions. These functions are not exported.
### -------------------------------------------------------------------------
###

# A copy of MethylSeekR::plotAlphaDistributionOneChr that can be safely run in
# parallel because it does not plot the histogram but simply returns it.
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
  score <- MethylSeekR:::calculateAlphaDistr(M, T, nCGbin, num.cores = 1)
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
    MethylSeekR:::calculateAlphaDistr(M, T, nCGbin, num.cores = 1)
  }, m = l_msrgr, indx = l_indx)
  par(mfrow = c(nrow, ncol))
  mapply(function(score, sn, chr.sel) {
    hist(score, probability = TRUE, breaks = 30,
         xlab = sprintf("posterior mean of alpha (%s)", chr.sel), main = sn)
  }, score = l_score, sn = sn, MoreArgs = list(chr.sel = chr.sel))
  return(NULL)
}
