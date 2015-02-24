# R package: methsim

__This package is in early development. Key functionality is missing and the interface may change without notice__.

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
