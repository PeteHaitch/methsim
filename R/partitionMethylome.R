### =========================================================================
### partitionMethylome
### -------------------------------------------------------------------------
###

#' Partition a methylome into UMRs, LMRs, PMRs and others.
#'
#' Paritition a methylome into unmethylated regions (UMRs), lowly-methylated
#' regions (LMRs), partially-methylated regions (PMRs) and mostly-methylated
#' regions (others) based on the beta-values.
#'
#' @param umrs_lmrs The output of
#' \code{MethylSeekR::\link[MethylSeekR]{segmentUMRsLMRs}}.
#' @param pmrs The output of \code{MethylSeekR\link[MethylSeekR]{segmentPMDs}}.
#' @return A \code{\link{PartitionedMethylome}}.
#' @export
partitionMethylome <- function(umrs_lmrs, pmrs) {

  if (!isDisjoint(umrs_lmrs)) {
    stop("Expected 'umrs_lmrs' to be disjoint. Please report this error.")
  }
  if (!isDisjoint(pmrs)) {
    stop("Expected 'pmrs' to be disjoint. Please report this error.")
  }

  gr <- c(granges(umrs_lmrs, use.mcols = FALSE),
         granges(pmrs, use.mcols = FALSE))
  mcols(gr)$type <- c(mcols(umrs_lmrs)$type, mcols(pmrs)$type)
  gr <- sort(gr)
  gr <- disjoin(gr)
  # Get the typo of each range in the partitioned methylome
  regionType <- rep(NA_character_, length(gr))
  # Get the type of those with ranges with exact matches in m
  equal_match_ol <- findOverlaps(gr, m, type = "equal")
  regionType[queryHits(equal_match_ol)] <-
    mcols(m)$type[subjectHits(equal_match_ol)]
  # Get the type of those ranges that overlap both UMR/LMR and PMR
  mo_idx <- findMostOverlapping(gr[-queryHits(equal_match_ol)], m)
  regionType[-queryHits(equal_match_ol)] <- mcols(m)$type[mo_idx]
  regionType <- gsub("notPMD", "other", regionType)
  regionType <- gsub("PMD", "PMR", regionType)
  regionType <- factor(regionType, levels = c("UMR", "LMR", "PMR", "other"))
  seqlevels(gr) <- seqlevelsInUse(gr)
  new("PartitionedMethylome", gr, regionType = regionType)
}
