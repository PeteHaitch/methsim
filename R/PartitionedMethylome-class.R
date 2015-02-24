### =========================================================================
### PartitionedMethylome: An S4 class to formalise an informal GRanges-based
### class used by MethylSeekR.
### -------------------------------------------------------------------------
###

#' PartitionedMethylome class.
#' @aliases PartitionedMethylome
#'
#' @export
setClass("PartitionedMethylome",
         contains = "GRanges"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###


.valid.PartitionedMethylome.mcols <- function(object) {
  msg <- NULL
  if (!is.factor(mcols(object)$type) ||
        !identical(levels(mcols(object)$type), c("UMR", "LMR", "PMR", "other"))) {
    msg <- validMsg(msg, paste0("'type' must be a factor with levels ",
                                "'UMR', 'LMR', 'PMR' and 'other'."))
  }
  msg
}

.valid.PartitionedMethylome.GRanges <- function(object) {
  msg <- NULL
  if (!isDisjoint(object)) {
    msg <- validMsg(msg, "'PartitionedMethylome' must be disjoint.")
  }
  gaps <- gaps(granges(object))
  gaps <- gaps[strand(gaps) == "*"]
  if (length(gaps)) {
    msg <- validMsg(msg, "'PartitionedMethylome' must not contain gaps.")
  }
  msg
}

.valid.PartitionedMethylome <- function(object) {
  # Include all .valid.PartitionedMethylome.* functions in this vector
  msg <- c(.valid.PartitionedMethylome.mcols(object),
           .valid.PartitionedMethylome.GRanges(object))

  if (is.null(msg)){
    return(TRUE)
  } else{
    msg
  }
}

S4Vectors::setValidity2("PartitionedMethylome", .valid.PartitionedMethylome)
