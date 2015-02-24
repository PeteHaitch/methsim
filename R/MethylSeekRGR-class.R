### =========================================================================
### MethylSeekRGR: An S4 class to formalise the informal GRanges-based class
### used by MethylSeekR.
### -------------------------------------------------------------------------
###

#' MethylSeekRGR class
#'
#' An S4 class to formalise the informal
#' \code{GenomicRanges::\link[GenomicRanges]{GRanges}}-based class used by
#' \code{\link[MethylSeekR]{MethylSeekR-package}}.
#' @aliases MethylSeekRGR
#'
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
  msg
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

#' Coerce a MethPat object of CG 1-tuples to a list of MethylSeekR-compatible
#' GRanges object(s), one list element per sample.
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
#' @export
setAs("MethPat",
      "MethylSeekRGR",
      function(from) {
        # Adapted from MethylSeekR::readMethylome
        if (size(from) != 1L) {
          stop("'MethPat' object must contain data for 1-tuples.")
        }
        if (!identical(methtype(from), "CG")) {
          stop("'MethPat' object must have CG 'methtype'.")
        }
        if (!all(strand(from) == "*")) {
          stop(paste0("'MethPat' object must have processed by ",
                      "'MethylationTuples::collapseStrand'."))
        }
        if (ncol(from) != 1L) {
          stop(paste0("'MethPat' object must only contain data on one sample. ",
                      "For multiple samples please use the ",
                      "as(MethPat, 'MethylSeekRGRList') method."))
        }
        idx <- !is.na(getCoverage(from))
        if (!any(idx)) {
          return(MethylSeekRGR())
        } else {
          msrgr <- MethylSeekRGR(seqnames(from)[idx], ranges(from)[idx],
                                 strand(from)[idx], getCoverage(from)[idx],
                                 assay(from, "M")[idx], seqinfo(from))
          msrgr <- sort(msrgr)
        }
        sn <- colnames(from)
        mean_cov <- mean(mcols(msrgr)$T)
        if (mean_cov < 10L) {
          warning(paste0("For CpGs with at least one read, sample '", sn,
                         "' ", "has mean coverage = ", mean_cov, "\nThe ",
                         "MethylSeekR developers do not recommend the use ",
                         "of MethylSeekR for methylomes with mean coverage ",
                         "< 10X."))
        }
        msrgr
      }
)
