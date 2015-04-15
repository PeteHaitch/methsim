# Copyright (C) 2015 Peter Hickey
#
# This file is part of methsim.
#
# methsim is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# methsim is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with methsim  If not, see <http://www.gnu.org/licenses/>.

#' Simulate DNA methylation data.
#'
#' \pkg{methsim} is software to simulate DNA methylation sequencing data.
#' \code{methsim} can currently simulate data from bisulfite-sequencing assays,
#' such as methylC-seq or BS-seq. It uses a non-stationary, inhomogeneous Markov
#' model. This model can incorporate the strong spatial dependence of DNA
#' methylation. Simulation parameters can be estimated from data or specified
#' by the user.
#'
#' Please refer to the vignettes to see how to use the \pkg{methsim} package.
#'
#' @docType package
#' @name methsim
#' @useDynLib methsim, .registration = TRUE
#' @import methods
#' @import GenomicRanges
#' @import Rcpp
#' @import MethylationTuples
#' @import MethylSeekR
#' @import data.table
#' @import BiocParallel
#' @importFrom mhsmm hmmspec hmmfit
NULL

#' \code{\link{SimulateMethylomeParam}} object for the \emph{ADS} sample.
#'
#' A \code{\link{SimulateMethylomeParam}} object for the \emph{ADS} sample from
#' Lister \emph{et al.} (2011).
#'
#' @format A \code{\link{SimulateMethylomeParam}} object.
#'
#' @source Data originally published in Lister, R. \emph{et al.} Hotspots of
#' aberrant epigenomic reprogramming in human induced pluripotent stem cells.
#' \emph{Nature} \strong{471}, 68–73 (2011)
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/21289626} and re-analysed by Peter
#' Hickey as follows.
#'
#' In brief, aligned sequencing data were downloaded from
#' \url{http://neomorph.salk.edu/ips_methylomes/data.html} and converted to
#' \code{BAM} format using Python scripts available from
#' \url{https://github.com/PeteHaitch/Lister2BAM}.
#' Methylation patterns at m-tuples were extracted using \code{methtuple v1.4.0}
#' (\url{https://github.com/PeteHaitch/methtuple}) and processed using
#' \code{MethylationTuples}
#' (\url{https://github.com/PeteHaitch/MethylationTuples}) and \code{methsim}
#' (\url{https://github.com/PeteHaitch/methsim}). See my PhD thesis
#' for further details (\url{https://github.com/PeteHaitch/phd_thesis}), in
#' particular the \strong{Datasets} chapter.
"ADS"

#' \code{\link{SimulateMethylomeParam}} object for the \emph{ADS-adipose}
#' sample.
#'
#' A \code{\link{SimulateMethylomeParam}} object for the \emph{ADS-adipose}
#' sample from Lister \emph{et al.} (2011).
#'
#' @format A \code{\link{SimulateMethylomeParam}} object.
#'
#' @source Data originally published in Lister, R. \emph{et al.} Hotspots of
#' aberrant epigenomic reprogramming in human induced pluripotent stem cells.
#' \emph{Nature} \strong{471}, 68–73 (2011)
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/21289626} and re-analysed by Peter
#' Hickey as follows.
#'
#' In brief, aligned sequencing data were downloaded from
#' \url{http://neomorph.salk.edu/ips_methylomes/data.html} and converted to
#' \code{BAM} format using Python scripts available from
#' \url{https://github.com/PeteHaitch/Lister2BAM}.
#' Methylation patterns at m-tuples were extracted using \code{methtuple v1.4.0}
#' (\url{https://github.com/PeteHaitch/methtuple}) and processed using
#' \code{MethylationTuples}
#' (\url{https://github.com/PeteHaitch/MethylationTuples}) and \code{methsim}
#' (\url{https://github.com/PeteHaitch/methsim}). See my PhD thesis
#' for further details (\url{https://github.com/PeteHaitch/phd_thesis}), in
#' particular the \strong{Datasets} chapter.
"ADS_adipose"

#' \code{\link{SimulateMethylomeParam}} object for the \emph{ADS-iPSC} sample.
#'
#' A \code{\link{SimulateMethylomeParam}} object for the \emph{ADS-iPSC} sample
#' from Lister \emph{et al.} (2011).
#'
#' @format A \code{\link{SimulateMethylomeParam}} object.
#'
#' @source Data originally published in Lister, R. \emph{et al.} Hotspots of
#' aberrant epigenomic reprogramming in human induced pluripotent stem cells.
#' \emph{Nature} \strong{471}, 68–73 (2011)
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/21289626} and re-analysed by Peter
#' Hickey as follows.
#'
#' In brief, aligned sequencing data were downloaded from
#' \url{http://neomorph.salk.edu/ips_methylomes/data.html} and converted to
#' \code{BAM} format using Python scripts available from
#' \url{https://github.com/PeteHaitch/Lister2BAM}.
#' Methylation patterns at m-tuples were extracted using \code{methtuple v1.4.0}
#' (\url{https://github.com/PeteHaitch/methtuple}) and processed using
#' \code{MethylationTuples}
#' (\url{https://github.com/PeteHaitch/MethylationTuples}) and \code{methsim}
#' (\url{https://github.com/PeteHaitch/methsim}). See my PhD thesis
#' for further details (\url{https://github.com/PeteHaitch/phd_thesis}), in
#' particular the \strong{Datasets} chapter.
"ADS_iPSC"
