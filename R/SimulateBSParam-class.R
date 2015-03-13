### =========================================================================
### SimulateBSParam: An S4 class to store the parameters used to simulate
### bisulfite-sequencing data from a SimulatedMethylome object;
### BS = [WGBS | RRBS | eRRBS]
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### list(SimulatedMethylome, aveCov, errorRate, target)
### SimulatedMethylome: A SimulatedMethylome object.
###             aveCov: The average sequencing coverage to simulate.
###          errorRate: The error rate (combined sequencing error and
###                     bisulfite-conversion error) to use in the simulation.
###             target: The co-ordinates of baits/enriched regions used in
###                     simulating RRBS/eRRBS data.

setClassUnion("GRangesORNULL",
              members = c("GRanges", "NULL")
)

#' SimulateBSParam class
#'
#' An S4 class for the parameters used by
#' \code{\link{simulate,SimulateBS-method}}.
#'
#' @include SimulateMethylomeParam-class.R
#'
#' @aliases SimulateBSParam
#'
#' @export
setClass("SimulateBSParam",
         slots = list(
           SimulatedMethylome = "SimulatedMethylome",
           aveCov = "numeric",
           errorRate = "numeric",
           target = "GRangesORNULL")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.SimulateBSParam.SimulatedMethylome <- function(object) {
  msg <- NULL
  if (!is(object@SimulatedMethylome, "SimulatedMethylome")) {
    msg <- Biobase::validMsg(msg, paste0("'SimulatedMethylome' slot must be a ",
                                         "'SimulatedMethylome' object."))
  }
  msg
}

.valid.SimulateBSParam.aveCov <- function(object) {
  msg <- NULL
  if (!is.numeric(object@aveCov) ||
      object@aveCov <= 0) {
    msg <- Biobase::validMsg(msg, paste0("'aveCov' slot must be a postive ",
                                         "number."))
  }
  msg
}

.valid.SimulateBSParam.errorRate <- function(object) {
  msg <- NULL
  if (!is.numeric(object@errorRate) ||
      object@errorRate < 0 ||
      object@errorRate > 1) {
    msg <- Biobase::validMsg(msg, paste0("'errorRate' slot must be a number ",
                                         "between 0 and 1."))
  }
  msg
}

.valid.SimulateBSParam.target <- function(object) {
  msg <- NULL
  if (!is(object@target, "GRangesORNULL")) {
    msg <- Biobase::validMsg(msg, paste0("'target' slot must be a GRanges ",
                                         "object or NULL."))
  } else {
    if (is(object@target, "GRanges")) {
      if (!identical(seqinfo(object@target),
                     seqinfo(object@SimulatedMethylome))) {
        msg <- Biobase::validMsg(msg, paste0("'target' slot and ",
                                             "'SimulatedMethylome' slot must ",
                                             " have identical 'seqinfo'."))
      }
    }
  }
  msg
}

.valid.SimulateBSParam <- function(object) {
  # Include all .valid.SimulateBSParam.* functions in this vector
  msg <- c(.valid.SimulateBSParam.SimulatedMethylome(object),
           .valid.SimulateBSParam.aveCov(object),
           .valid.SimulateBSParam.errorRate(object),
           .valid.SimulateBSParam.target(object))

  if (is.null(msg)){
    return(TRUE)
  } else{
    msg
  }
}

S4Vectors::setValidity2("SimulateBSParam", .valid.SimulateBSParam)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
SimulateBSParam <- function(SimulatedMethylome,
                            aveCov = 30L,
                            errorRate = 0.01,
                            target = NULL) {

  # TODO: Argument checks
  new("SimulateBSParam",
      SimulatedMethylome = SimulatedMethylome,
      aveCov = aveCov,
      errorRate = errorRate,
      target = target)
}
