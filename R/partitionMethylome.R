### =========================================================================
### partitionMethylome
### -------------------------------------------------------------------------
###

#' Partition a methylome into UMRs, LMRs, PMRs and MMRs.
#'
#' Paritition a methylome into unmethylated regions (UMRs), lowly-methylated
#' regions (LMRs), partially-methylated regions (PMRs) and mostly-methylated
#' regions (MMRs) based on the beta-values.
#'
#' This function is a convenience function wrapping around several functions
#' provided by the \code{MethylSeekR} Bioconductor package. It is not a
#' replacement.
#'
#' @param umrs_lmrs The output of
#' \code{MethylSeekR\link[MethylSeekR]{segmentUMRsLMRs}}.
#' @param pmrs The output of \code{MethylSeekR\link[MethylSeekR]{segmentPMDs}}.
#' @return A \code{\link{PartitionedMethylome}}.
#' @export
partitionMethylome <- function(umrs_lmrs, pmrs) {

  if (!isDisjoint(umrs_lmrs)) {
    stop("Expect 'umrs_lmrs' to be disjoint. Please report this error.")
  }
  if (!isDisjoint(pmrs)) {
    stop("Expect 'pmrs' to be disjoint. Please report this error.")
  }

  m <- c(granges(umrs_lmrs, use.mcols = FALSE),
         granges(pmrs, use.mcols = FALSE))
  m$type <- c(umrs_lmrs$type, pmrs$type)
  m <- sort(m)
  m_disjoint <- disjoin(m)
  mo_idx <- findMostOverlapping(m, m_disjoint)
  # UP TO HERE
}