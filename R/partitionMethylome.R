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

  # Combine umrs_lmrs and pmrs while ensuring that only the 'type' mcol is
  # retained.
  gr <- c(granges(umrs_lmrs, use.mcols = FALSE),
         granges(pmrs, use.mcols = FALSE))
  mcols(gr)$type <- c(mcols(umrs_lmrs)$type, mcols(pmrs)$type)
  # TODO: Remove sort() once a unit test is written to ensure that
  # disjoin(sort(gr)) == disjoint(gr).
  disjoint_gr <- disjoin(sort(gr))
  regionType <- rep(NA_character_, length(disjoint_gr))
  # Get the regionType of the ranges in disjoint_gr with exact matches in gr
  equal_match_ol <- findOverlaps(disjoint_gr, gr, type = "equal")
  regionType[queryHits(equal_match_ol)] <-
    mcols(gr)$type[subjectHits(equal_match_ol)]
  # Get the regionType of the ranges in disjoint_gr that overlap both UMR/LMR
  # and PMR.
  mo_idx <- findMostOverlapping(disjoint_gr[-queryHits(equal_match_ol)], gr)
  regionType[-queryHits(equal_match_ol)] <- mcols(gr)$type[mo_idx]
  regionType <- gsub("notPMD", "other", regionType)
  regionType <- gsub("PMD", "PMR", regionType)
  regionType <- factor(regionType, levels = c("UMR", "LMR", "PMR", "other"))
  seqlevels(disjoint_gr) <- seqlevelsInUse(disjoint_gr)
  new("PartitionedMethylome", disjoint_gr, regionType = regionType)
}
