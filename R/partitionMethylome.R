### =========================================================================
### partitionMethylome
### -------------------------------------------------------------------------
###

# TODO: Create the PartionedMethylome class and return an object of this class.
#' Partition a methylome into UMRs, LMRs, PMRs and others.
#'
#' Paritition a methylome into unmethylated regions (UMRs), lowly-methylated
#' regions (LMRs), partially-methylated regions (PMRs) and mostly-methylated
#' regions (others) based on the beta-values.
#'
#' @param umrs_lmrs The output of
#' \code{MethylSeekR::\link[MethylSeekR]{segmentUMRsLMRs}}.
#' @param pmrs The output of \code{MethylSeekR\link[MethylSeekR]{segmentPMDs}}.
#' @param bsgenome A \code{BSgenome::\link[BSgenome]{BSgenome}} object of the
#' relevant organism.
#' @return A \code{\link{PartitionedMethylome}}.
#' @export
partitionMethylome <- function(umrs_lmrs, pmrs) {

  if (!isDisjoint(umrs_lmrs)) {
    stop("Expected 'umrs_lmrs' to be disjoint. Please report this error.")
  }
  if (!isDisjoint(pmrs)) {
    stop("Expected 'pmrs' to be disjoint. Please report this error.")
  }

  m <- c(granges(umrs_lmrs, use.mcols = FALSE),
         granges(pmrs, use.mcols = FALSE))
  m$type <- c(umrs_lmrs$type, pmrs$type)
  pm <- sort(m)
  pm <- disjoin(m)
  # Get the typo of each range in the partitioned methylome
  type <- rep(NA_character_, length(pm))
  # Get the type of those with ranges with exact matches in m
  equal_match_ol <- findOverlaps(pm, m, type = "equal")
  type[queryHits(equal_match_ol)] <- mcols(m)$type[subjectHits(equal_match_ol)]
  # Get the type of those ranges that overlap both UMR/LMR and PMR
  mo_idx <- findMostOverlapping(pm[-queryHits(equal_match_ol)], m)
  type[-queryHits(equal_match_ol)] <- mcols(m)$type[mo_idx]
  type <- gsub("notPMD", "other", type)
  type <- gsub("PMD", "PMR", type)
  type <- factor(type, levels = c("UMR", "LMR", "PMR", "other"))
  mcols(pm) <- DataFrame(type = type)
  seqlevels(pm) <- seqlevelsInUse(pm)
  new("PartitionedMethylome", pm)
}