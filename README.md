[![Project Status: Abandoned - Initial development has started, but there has not yet been a stable, usable release; the project has been abandoned and the author\(s\) do not intend on continuing development.](http://www.repostatus.org/badges/latest/abandoned.svg)]

# R package: methsim
------------------------


__This package is no longer under active development. It was developed as an experiment as part of my PhD research. The repository remains available mostly for archival purposes, but forks are of course welcome under the GPL (>= 2) license.__

`methsim` is software to simulate DNA methylation sequencing data. `methsim` 
can currently simulate data from bisulfite-sequencing assays such as 
_methylC-seq_ or _BS-seq_. It uses a non-stationary, inhomogeneous Markov model. 
This model can incorporate the strong spatial dependence of DNA methylation. 
Simulation parameters can be estimated from data or specified by the user.

`methsim` is in development and can only be installed using the development 
version of Bioconductor. Please first read 
[these instructions on installing the development version of Bioconductor](http://www.bioconductor.org/developers/how-to/useDevel/).

```R
# Install devtools if not already installed
if (suppressWarnings(!require(devtools))) {
  install.packages('devtools')
}
# Install development version of methsim
devtools::install_github("PeteHaitch/methsim")
```

## R CMD check results
Travis CI: <a href="https://travis-ci.org/PeteHaitch/methsim"><img src="https://travis-ci.org/PeteHaitch/methsim.svg?branch=master" alt="Build status"></a>

## Test coverage status
coveralls.io: [![Coverage Status](https://coveralls.io/repos/PeteHaitch/methsim/badge.svg?branch=master)](https://coveralls.io/r/PeteHaitch/methsim?branch=master)
