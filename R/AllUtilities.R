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

# A copy of MethylSeekR::segmentPMDs that can be safely run in
# parallel because it does not create the plot but returns it as an object.
# TODO: Document and robust-ify (long-term)
segmentPMDs <- function(m, chr.sel, seqLengths, nCGbin = 101, plot = FALSE) {
  hmm.model <- trainPMDHMM(m, chr.sel, nCGbin, plot = FALSE)
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
    plot(hist, freq = FALSE)
    lines(x, lines1, type = "l", col = "red")
    lines(x, lines2, type = "l", col = "green")
    return(list(startval = startval, hist = hist, x = x, lines1 = lines1,
                lines2 = lines2))
  } else {
    list(startval = startval, hist = hist, x = x, lines1 = lines1,
         lines2 = lines2)
  }
}

# A copy of MethylSeekR::calculateFDRs that adds the sample name to the title
# of each plot.
# TODO: Document and robust-ify (long-term)
calculateFDRs <- function(m, sn, CGIs, PMDs = NA, pdfFilename = NULL, num.cores = 1,
                          nCpG.smoothing = 3,
                          meth.cutoffs = seq(0.3, 0.7, by = 0.1),
                          nCpG.cutoffs = seq(1, 6, by = 1), minCover = 5) {
  m = m[values(m)[, 1] >= minCover]
  message("calculating false discovery rate")
  n <- length(meth.cutoffs) * length(nCpG.cutoffs)
  parameters <- vector("list", n)
  count = 1
  for (meth.cutoff in meth.cutoffs) {
    for (k in nCpG.cutoffs) {
      parameters[[count]] <- c(meth.cutoff, k)
      count <- count + 1
    }
  }
  meth.notrand <- m
  if (class(PMDs) == "GRanges") {
    message("removing PMDs for randomization")
    meth.notrand <- subsetByOverlaps(m, PMDs[values(PMDs)$type ==
                                               "notPMD"])
  }
  meth.rand <- meth.notrand
  ov <- findOverlaps(meth.rand, CGIs)
  meth.rand <- meth.rand[-unique(queryHits(ov))]
  values(meth.rand) <- values(meth.rand)[sample(length(meth.rand)),
                                         ]
  res <- mclapply(parameters, function(params) {
    meth.cutoff <- params[1]
    k <- params[2]
    mean.meth <- runmean(Rle(values(meth.notrand)[, 2]/values(meth.notrand)[,
                                                                            1]), k = nCpG.smoothing, endrule = "constant")
    indx <- mean.meth < meth.cutoff
    nSeg <- sum(runValue(indx) == TRUE & runLength(indx) >=
                  k)
    mean.meth = runmean(Rle(values(meth.rand)[, 2]/values(meth.rand)[,
                                                                     1]), k = nCpG.smoothing, endrule = "constant")
    indx <- mean.meth < meth.cutoff
    nSeg.rand <- sum(runValue(indx) == TRUE & runLength(indx) >=
                       k)
    c(nSeg, nSeg.rand)
  }, mc.cores = num.cores)
  FDRs = 100 * sapply(res, function(x) {
    x[2]/x[1]
  })
  tmp = matrix(NA, nrow = length(meth.cutoffs), ncol = length(nCpG.cutoffs))
  rownames(tmp) = meth.cutoffs
  colnames(tmp) = nCpG.cutoffs
  count = 1
  for (meth.cutoff in meth.cutoffs) {
    for (k in nCpG.cutoffs) {
      tmp[as.character(meth.cutoff), as.character(k)] = FDRs[count]
      count = count + 1
    }
  }
  FDRs = tmp
  numSegments = sapply(res, function(x) {
    x[1]
  })
  tmp = matrix(NA, nrow = length(meth.cutoffs), ncol = length(nCpG.cutoffs))
  rownames(tmp) = meth.cutoffs
  colnames(tmp) = nCpG.cutoffs
  count = 1
  for (meth.cutoff in meth.cutoffs) {
    for (k in nCpG.cutoffs) {
      tmp[as.character(meth.cutoff), as.character(k)] = numSegments[count]
      count = count + 1
    }
  }
  numSegments = tmp
  rownames(FDRs) = as.character(100 * as.numeric(rownames(FDRs)))
  rownames(numSegments) = as.character(100 * as.numeric(rownames(numSegments)))
  if (!is.null(pdfFilename)) {
    pdf(pdfFilename, width = 9, height = 4.5)
  }
  par(mfrow = c(1, 2))
  barplot(pmin((t(FDRs)), 20), beside = TRUE, ylab = "FDR (%)",
          ylim = c(0, 20), xlab = "methylation cut-off (%)", main = sn)
  barplot(t(numSegments), beside = TRUE, ylab = "number of segments",
          xlab = "methylation cut-off (%)", main = sn)
  legend("topleft",
         legend = paste(colnames(FDRs),
                        c("CpG", rep("CpGs", ncol(FDRs) - 1)), sep = " "),
         fill = grey.colors(ncol(FDRs)), bty = "n")
  if (!is.null(pdfFilename))
    dev.off()
  rownames(FDRs) = as.character(as.numeric(rownames(FDRs))/100)
  rownames(numSegments) = as.character(as.numeric(rownames(numSegments))/100)
  list(FDRs = FDRs, numSegments = numSegments)
}

# A copy of MethylSeekR::segmentUMRsLMRs that adds the sample name to the title
# of each plot. Actually, it doesn't create the plot but simply returns the
# data necessary for creating it.
# TODO: Document and robust-ify (long-term)
segmentUMRsLMRs <- function(m, meth.cutoff = 0.5, nCpG.cutoff = 3, PMDs = NA,
                            myGenomeSeq, seqLengths, nCpG.smoothing = 3,
                            minCover = 5) {
  nCG.classification <- 30
  message("identifying UMRs and LMRs")
  m = m[values(m)[, 1] >= minCover]
  # TODO: Not sure why I have to explicitly coerce to a character vector but
  # the Rle,table-method isn't working as expected.
  nCGsPerChr = table(as.character(seqnames(m)))
  chrs = names(nCGsPerChr)[nCGsPerChr >= nCpG.smoothing]
  res <- mclapply(chrs, function(chr) {
    sel <- which(as.character(seqnames(m)) == chr)
    mean.meth <- runmean(Rle(values(m)[sel, 2] /
                               values(m)[sel, 1]), k = nCpG.smoothing,
                         endrule = "constant")
    indx <- mean.meth < meth.cutoff
    runValue(indx)[runLength(indx) < nCpG.cutoff & runValue(indx) ==
                     TRUE] = FALSE
    runValue(indx)[runLength(indx) < nCpG.cutoff & runValue(indx) ==
                     FALSE] = TRUE
    tmp.ids <- rep(1:length(runLength(indx)), runLength(indx))
    tmp <- split(1:length(sel), tmp.ids)
    tmp <- tmp[runValue(indx) == TRUE]
    if (length(tmp) > 0) {
      coords <- cbind(sapply(tmp, min), sapply(tmp, max))
      starts <- round((start(m)[sel[pmax(1, coords[, 1] - 1)]] +
                         start(m)[sel[coords[, 1]]]) / 2)
      ends <- round((start(m)[sel[coords[, 2]]] +
                       start(m)[sel[pmin(length(sel), coords[, 2] + 1)]]) / 2)
      hmr.gr = GRanges(seqnames = unique(seqnames(m[sel])),
                       strand = "*", ranges = IRanges(starts, ends),
                       seqlengths = seqLengths)
    }
    else {
      hmr.gr = GRanges(, seqlengths = seqLengths)
    }
    hmr.gr
  }, mc.cores = 1L)
  segments.gr = do.call(c, unname(res))
  if (class(PMDs) == "GRanges") {
    segments.gr = subsetByOverlaps(segments.gr,
                                   PMDs[values(PMDs)$type == "notPMD"])
  }
  nCG = vcountPattern("CG", getSeq(myGenomeSeq,
                                   resize(segments.gr, width(segments.gr),
                                          fix = "start"), as.character = FALSE))
  ov <- findOverlaps(m, segments.gr)
  T = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), sum)
  M = tapply(values(m[queryHits(ov)])[, 2], subjectHits(ov), sum)
  nCG.segmentation = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov),
                            length)
  median.meth = tapply(as.vector(runmean(Rle(values(m[queryHits(ov)])[, 2] /
                                               values(m[queryHits(ov)])[, 1]),
                                         nCpG.smoothing, endrule = "constant")),
                       subjectHits(ov), median)
  median.meth = pmax(0, median.meth)
  if (!all.equal(as.numeric(names(T)), 1:length(segments.gr))) {
    message("error in calculating methylation levels for PMDs")
  }
  type = c("UMR", "LMR")[as.integer(nCG < nCG.classification) + 1]
  values(segments.gr) = DataFrame(nCG.segmentation, nCG, T, M, pmeth = M / T,
                                  median.meth = median.meth, type)
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  list(x = log2(values(segments.gr)$nCG),
       y = 100 * values(segments.gr)$median.meth,
       colramp = jet.colors,
       xlab = "log2 number of CpGs in segment",
       ylab = "median methylation (%)",
       v = log2(nCG.classification),
       segments.gr = segments.gr)
}

# A copy of MethylSeekR::plotFinalSegmentation that can be safely run in
# parallel because it does not create the plot but returns it as an object.
# TODO: Document and robust-ify (long-term)
plotFinalSegmentation <- function(m, sn, segs, PMDs = NA, meth.cutoff,
                                  numRegions = 1, minCover = 5,
                                  nCpG.smoothing = 3, plot = FALSE) {
  m = m[values(m)[, 1] >= minCover]
  cols.seg = c("#377EB8", "#E41A1C")
  names(cols.seg) = c("UMR", "LMR")
  col.PMD = "#4DAF4A"
  pch.seg <- c(LMR = 24, UMR = 22)
  chrs = unique(as.character(seqnames(m)))
  height = 1
  len = 5000
  nCGperChr = table(as.character(seqnames(m)))
  indx = which(nCGperChr[chrs] < 3 * len)
  if (length(indx) > 0) {
    chrs = chrs[-indx]
  }
  exist.PMDs <- FALSE
  if (class(PMDs) == "GRanges") {
    exist.PMDs <- TRUE
  }
  for (ii in 1:numRegions) {
    chr = sample(chrs, 1)
    i = seqnames(m) == chr
    s = sample(1:(length(m[i]) - 3 * len), 1)
    par(mfrow = c(6, 1))
    for (gg in 1:3) {
      indx = (s + (gg - 1) * len):(s + (gg) * len)
      seg.sel = subsetByOverlaps(segs, m[i][indx])
      if (length(seg.sel) > 0) {
        mysegcol <- cols.seg[as.character(values(seg.sel)$type)]
        mysegpch <- pch.seg[as.character(values(seg.sel)$type)]
      }
      PMD.sel = NA
      if (exist.PMDs) {
        PMD.sel = subsetByOverlaps(PMDs[values(PMDs)$type == "PMD"], m[i][indx])
      }
      add.mar <- 2
      par(mar = c(3 - add.mar, 4, 1 + add.mar, 2))
      plot(start(m[i])[indx], 100 * values(m[i])[indx, 2] /
             values(m[i])[indx, 1], pch = "*", ylim = c(-18, 100), cex = 0.9,
           ylab = "methylation", axes = FALSE, xlab = "", main = sn)
      axis(2, at = c(0, 50, 100))
      box()
      text(par("usr")[1] + par("cxy")[1] * 0.2, par("usr")[3] +
             par("cxy")[2] * 0.3, adj = c(0, 0), label = "raw")
      if (exist.PMDs & length(PMD.sel) > 0) {
        rect(start(PMD.sel), -18 - height / 4, end(PMD.sel),
             -18 + height / 4, lwd = 2, col = col.PMD, border = col.PMD)
      }
      if (length(seg.sel) > 0) {
        points(x = mid(ranges(seg.sel)),
               y = rep(-18, length(seg.sel)), pch = mysegpch, cex = 1.3,
               col = "black", lwd = 0.5, bg = mysegcol)
      }
      par(mar = c(2 + add.mar, 4, 2 - add.mar, 2))
      plot(start(m[i])[indx], 100 *
             as.vector(runmean(Rle(values(m[i])[indx, 2] /
                                     values(m[i])[indx, 1]),
                               nCpG.smoothing, endrule = "constant")),
           pch = "*", ylim = c(-18, 100), cex = 0.9, ylab = "methylation",
           xlab = sprintf("position on %s", chr), axes = FALSE)
      axis(1)
      axis(2, at = c(0, 50, 100))
      box()
      text(par("usr")[1] + par("cxy")[1] * 0.2, par("usr")[3] +
             par("cxy")[2] * 0.3, adj = c(0, 0), label = "smoothed")
      if (exist.PMDs & length(PMD.sel) > 0) {
        rect(start(PMD.sel), -18 - height / 4, end(PMD.sel),
             -18 + height / 4, lwd = 2, col = col.PMD, border = col.PMD)
      }
      if (length(seg.sel) > 0) {
        points(x = mid(ranges(seg.sel)),
               y = rep(-18, length(seg.sel)), pch = mysegpch, cex = 1.3,
               col = "black", lwd = 0.5, bg = mysegcol)
      }
      abline(h = 100 * meth.cutoff, lty = 5, col = "darkgrey")
    }
  }
}

# Compute a function on a MethPat object for each sample stratified by its
# PartitionedMethylome

# TODO: Support other FUN, e.g., methLevel and methLevelCor
#
#' Apply a function to a
#' \code{MethylationTuples::\link[MethylationTuples]{MethPat}} object,
#' with results stratified by the region type of the accompanying
#' \code{\link{PartitionedMethylome}}.
#'
#' @param FUN The function to be applied to each element of \code{pm}: see
#' 'Details'. Only a limited number of functions are currently supported,
#' specifically \code{MethylationTuples::\link[MethylationTuples]{methLevel}()},
#' \code{MethylationTuples::\link[MethylationTuples]{cometh}()}, and
#' \code{MethylationTuples::\link[MethylationTuples]{patternFreqs}()}.
#' @param pm A \code{\link{PartitionedMethylome}} object.
#' @param methpat A \code{MethylationTuples::\link[MethylationTuples]{MethPat}}
#' object.
#' @param min_cov An \code{integer} specifying the minimum coverage required
#' in order to use an m-tuple in the analysis.
#' @param ... Optional arguments to \code{FUN}.
#'
#' @return A \code{data.table::\link[data.table]{data.table}}.
#'
#' @keywords internal
funByPM <- function(FUN, pm, methpat, min_cov, ...) {
  # Argument checks
  stopifnot(is(pm, "PartitionedMethylome"))
  if (ncol(methpat) > 1L) {
    stop("'methpat' must only contain data on a single sample.")
  }
  stopifnot(min_cov >= 0L)
  FUN <- match.fun(FUN)
  # TODO: Make more robust method for checking whether it is possible to use
  # FUN with methsim:::funByPM()
  AVAILABLE_FUNCTIONS <- c(MethylationTuples::methLevel,
                           MethylationTuples::cometh,
                           MethylationTuples::patternFreqs)
  stopifnot(any(vapply(X = AVAILABLE_FUNCTIONS, FUN = identical,
                       FUN.VALUE = logical(1), FUN)))

  # Drop m-tuples with insufficient coverage
  cov <- getCoverage(methpat)
  methpat <- methpat[as.vector(!is.na(cov) & (cov >= min_cov)), ]
  # Find the region type of each m-tuple
  ol <- findOverlaps(methpat, pm, type = "within")
  # Drop those m-tuples spanning a region boundary
  methpat <- methpat[queryHits(ol)]
  type <- mcols(pm)$type[subjectHits(ol)]
  val <- FUN(methpat, min_cov = min_cov, ...)
  if (!is(val, "data.table")) {
    val <- data.table::as.data.table(val)
  }
  val[, type := type]
  val
}
