### =========================================================================
### MethylSeekRGR: An S4 class to formalise the informal GRanges-based class
### used by MethylSeekR.
### -------------------------------------------------------------------------
###

#' @export
setClass("MethylSeekRGR",
         contains = "GRanges"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

# TODO: Really, the 'T' and 'M' columns should probably be extraColumnSlots.
# But that's something that needs to be fixed in MethylSeekR, not in methsim.
.valid.MethylSeekRGR.mcols <- function(object) {
  msg <- NULL
  if (!all(c("T", "M") %in% colnames(mcols(object)))) {
    msg <- validMsg(msg, paste0("Must contain metadata columns 'T' ",
                                         "and 'M'"))
  }
  if (!is(object$T, "integer") || !is(object$M, "integer")) {
    msg <- validMsg(msg, "'T' and 'M' must be 'integer' valued.")
  }
  if (any(object$M > object$T)) {
    msg <- validMsg(msg, "'M' > 'T' should not occur.")
  }
}

.valid.MethylSeekRGR.seqlengths <- function(object) {
  msg <- NULL
  if (length(object)) {
    if (any(is.na(seqlengths(object)))) {
      msg <- validMsg(msg, "Require valid seqlengths.")
    }
  }
  msg
}

.valid.MethylSeekRGR <- function(object) {
  # Include all .valid.MethylSeekRGR.* functions in this vector
  msg <- c(.valid.MethylSeekRGR.mcols(object),
           .valid.MethylSeekRGR.seqlengths(object))

  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

S4Vectors::setValidity2("MethylSeekRGR", .valid.MethylSeekRGR)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
MethylSeekRGR <- function(seqnames = Rle(), ranges = IRanges(),
                          strand = Rle("*", length(seqnames)), T = integer(0),
                          M = integer(0), seqinfo = Seqinfo()) {
  gr <- GRanges(seqnames, ranges, strand, T = T, M = M, seqinfo = seqinfo)
  new("MethylSeekRGR", gr)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

#' Coerce a MethPat object of CG 1-tuples to MethylSeekR-compatible GRanges
#' object(s).
#'
#' @param methpat A \code{\link[MethylationTuples]{MethPat}} object containing CpG
#' 1-tuples. The \code{\link[MethylationTuples]{MethPat}} object should already
#' been processed with
#' \code{MethylationTuples::\link[MethylationTuples]{filterOutVariants}} and
#' \code{MethylationTuples::\link[MethylationTuples]{collapseStrand}}.
#' @return A list of \code{\link{MethylSeekRGR}} objects, which are
#' compatible with the \code{MethylSeekR} Bioconductor package.
#' @name as
#' @aliases MethylSeekRGR
setAs("MethPat",
      "MethylSeekRGR",
      function(from) {
        # Adapted from MethylSeekR::readMethylome
        if (size(methpat) != 1L) {
          stop("'methpat' must contain data for 1-tuples.")
        }
        if (!identical(methtype(methpat), "CG")) {
          stop("'methpat' must have CG 'methtype'.")
        }
        if (!all(strand(methpat) == "*")) {
          stop(paste0("'methpat' must have processed by ",
                      "'MethylationTuples::collapseStrand'."))
        }
        list_of_msrgr <- mapply(function(T, M, seqnames, ranges, strand,
                                         seqinfo) {
          msrgr <- MethylSeekRGR(seqnames, ranges, strand, T, M, seqinfo)
          msrgr <- msrgr[!is.na(T)]
          sort(msrgr)
        }, T = split(getCoverage(methpat), rep(1:ncol(methpat),
                                               each = nrow(methpat))),
        M = split(assay(methpat, "M"), rep(1:ncol(methpat),
                                           each = nrow(methpat))),
        MoreArgs = list(seqnames = seqnames(methpat), ranges = ranges(methpat),
                        strand = strand(methpat), seqinfo = seqinfo(methpat)))
        names(list_of_msrgr) <- colnames(methpat)
        mapply(function(msrgr, nm) {
          mean_cov <- mean(msrgr$T)
          if (mean_cov < 10) {
            warning(paste0("For CpGs with at least one read, sample '", nm,
                           "' ", "has mean coverage = ", mean_cov, "\nThe ",
                           "MethylSeekR developers do not the use of ",
                           "MethylSeekR for methylomes with mean coverage < ",
                           "10X."))
          }
        }, msrgr = list_of_msrgr, nm = names(msrgr))
      }
)
