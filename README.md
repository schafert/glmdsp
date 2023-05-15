
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glmdsp

<!-- badges: start -->
<!-- badges: end -->

glmdsp provides a method for trend filtering count time series with a
data adaptive Bayesian shrinkage model.

## Installation

You can install the development version of glmdsp from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("schafert/glmdsp")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(glmdsp)

# simulate data
beta <- bumps_sim(n = 200)
y <- rnbinom(n = 200, size = 5, mu = beta)
    
fit <- btf_nb(
  y = y,
  D = 1,
  nburn = 5000,
  evol0_sample = FALSE,
  verbose = FALSE,
  sigma_e = 1,
  chol0 = TRUE
)
#> 'as(<dsCMatrix>, "dgCMatrix")' is deprecated.
#> Use 'as(., "generalMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").
```
