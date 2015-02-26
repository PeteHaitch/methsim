### =========================================================================
### MethSimParam: An S4 class for the parameters used by methsim::simulate()
### -------------------------------------------------------------------------
###

#' MethSimParam class
#'
#' An S4 class for the parameters used by \code{\link{simulate}}.
#' @aliases MethSimParam
#'
#' @export
setClass("MethSimParam",
         slots = list(
           partitionedMethylome = "PartitionedMethylome",
           methLevel = "data.table",
           cometh = "data.table",
           patternFreqs = "data.table",
           sampleName = "character"),
         prototype = list(
           PartitionedMethylome = new("PartitionedMethylome"),
           methLevel = data.table(type = character(0),
                                  beta = numeric(0),
                                  N = integer(0)),
           cometh = data.table(IPD = integer(0),
                               type = character(0),
                               statistic = numeric(0),
                               N = integer(0)),
           patternFreqs = data.table(type = character(0),
                                     h1 = numeric(0),
                                     h2 = numeric(0),
                                     h3 = numeric(0),
                                     h4 = numeric(0),
                                     h5 = numeric(0),
                                     h6 = numeric(0),
                                     h7 = numeric(0),
                                     h8 = numeric(0),
                                     N = integer(0)),
           sampleName = character(0)

         )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.MethSimParam.PartitionedMethylome <- function(object) {
  msg <- NULL
  if (!is(object@PartitionedMethylome, "PartitionedMethylome")) {
    msg <- BBiobase::validMsg(msg, paste0("'PartitionedMethylome' slot must ",
                                          "be a 'PartitionedMethylome'."))
  }
  msg
}

.valid.MethSimParam.methLevel <- function(object) {
  msg <- NULL
  if (!is.data.table(object@methLevel)) {
    msg <- Biobase::validMsg(msg, paste0("'methLevel' slot must be a ",
                                         "'data.table'"))
  }
  if (!identical(colnames(object@methLevel), c("type", "beta", "N"))) {
    msg <- Biobase::validMsg(msg, paste0("Column names of the 'methLevel' ",
                                         "slot must be 'type', 'beta' and ",
                                         "'N'." ))
  }
  msg
}

.valid.MethSimParam.cometh <- function(object) {
  msg <- NULL
  if (!is.data.table(object@cometh)) {
    msg <- Biobase::validMsg(msg, paste0("'cometh' slot must be a ",
                                         "'data.table'"))
  }
  if (!identical(colnames(object@cometh), c("IPD", "type", "statistic",
                                               "N"))) {
    msg <- Biobase::validMsg(msg, paste0("Column names of the 'cometh' ",
                                         "slot must be 'IPD', 'type', ",
                                         "'statistic' and 'N'." ))
  }
  msg
}

.valid.MethSimParam.patternFreqs <- function(object) {
  msg <- NULL
  if (!is.data.table(object@patternFreqs)) {
    msg <- Biobase::validMsg(msg, paste0("'patternFreqs' slot must be a ",
                                         "'data.table'"))
  }
  h <- grep("^h[0-9]+", colnames(object@patternFreqs), value = TRUE)
  if (!identical(colnames(object@patternFreqs), c("type", h, "N"))) {
    msg <- Biobase::validMsg(msg, paste0("Column names of the 'patternFreqs' ",
                                         "slot must be, 'type', '",
                                         paste0(h, collapse = "', '"), "', ",
                                         "'statistic' and 'N'." ))
  }
  msg
}

.valid.MethSimParam.sampleName <- function(object) {
  msg <- NULL
  if (!is.character(object@sampleName) | length(object@sampleName) != 1L) {
    msg <- Biobase::validMsg(msg, paste0("'sampleName' slot must be a ",
                                         "'character' vector with length 1."))
  }
}

.valid.MethSimParam <- function(object) {
  # Include all .valid.MethSimParam.* functions in this vector
  msg <- c(.valid.MethSimParam.PartitionedMethylome(object),
           .valid.MethSimParam.methLevel(object),
           .valid.MethSimParam.cometh(object),
           .valid.MethSimParam.patternFreqs(object),
           .valid.MethSimParam.sampleName(object))

  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

S4Vectors::setValidity2("MethSimParam", .valid.MethSimParam)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# UP TO HERE: Fix TODOs then see if can construct a MethSimParam object for the
# E13BUF sample.
# TODO: MethSimParam() doesn't work but new("MethSimParam") does.
#' @export
MethSimParam <- function(meth_level, cometh, pattern_freqs, sample_name) {
  # TODO: Probably a good idea to check length of h, i.e., how many haplotypes
  # in pattern_freqs
  new("MethSimParam",
      partitioned_methylome = partitioned_methylome,
      methLevel = meth_level,
      cometh = cometh,
      patternFreqs = pattern_freqs,
      sampleName = sample_name)
}
