% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/methsim.R
\docType{data}
\name{ADS}
\alias{ADS}
\title{\code{\link{SimulateMethylomeParam}} object for the \emph{ADS} sample.}
\format{A \code{\link{SimulateMethylomeParam}} object.}
\source{
Data originally published in Lister, R. \emph{et al.} Hotspots of
aberrant epigenomic reprogramming in human induced pluripotent stem cells.
\emph{Nature} \strong{471}, 68–73 (2011)
\url{http://www.ncbi.nlm.nih.gov/pubmed/21289626} and re-analysed by Peter
Hickey as follows.

In brief, aligned sequencing data were downloaded from
\url{http://neomorph.salk.edu/ips_methylomes/data.html} and converted to
\code{BAM} format using Python scripts available from
\url{https://github.com/PeteHaitch/Lister2BAM}.
Methylation patterns at m-tuples were extracted using \code{methtuple v1.4.0}
(\url{https://github.com/PeteHaitch/methtuple}) and processed using
\code{MethylationTuples}
(\url{https://github.com/PeteHaitch/MethylationTuples}) and \code{methsim}
(\url{https://github.com/PeteHaitch/methsim}). See my PhD thesis
for further details (\url{https://github.com/PeteHaitch/phd_thesis}), in
particular the \strong{Datasets} chapter.
}
\usage{
ADS
}
\description{
A \code{\link{SimulateMethylomeParam}} object for the \emph{ADS} sample from
Lister \emph{et al.} (2011).
}
\keyword{datasets}

