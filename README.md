
<!-- README.md is generated from README.Rmd. Please edit that file -->

# estim <img src=man/figures/logo.png align="right" height="139" alt="logo"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/estim)](https://CRAN.R-project.org/package=estim)
[![R-CMD-check](https://github.com/thechibo/estim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/thechibo/estim/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/thechibo/estim/branch/main/graph/badge.svg)](https://app.codecov.io/gh/thechibo/estim?branch=main)
[![Codecov test
coverage](https://codecov.io/gh/thechibo/estimators/graph/badge.svg)](https://app.codecov.io/gh/thechibo/estimators)
<!-- badges: end -->

## Introduction

The `estim` R package performs parameter estimation in common
distribution families, making moment, maximum likelihood, and
state-of-the-art estimators more accessible.

### Key Features

- The common d, p, q, r function family for each distribution
  (e.g. dnorm, pnorm, qnorm, rnorm) is enriched with
  - the ll counterpart (e.g. llnorm) that calculates the log-likelihood,
  - the e counterpart (e.g. enorm) that performs parameter estimation,
  - the v counterpart (e.g. vnorm) that calculates the asymptotic
    variance-covariance matrix of an estimator.
- Moment functions (mean, median, mode, var, sd, skew, kurt) as well as
  functions that calculate the entropy and the Fisher Information are
  available for all distributions.
- An S4-class distribution system is developed. Distribution objects can
  be created (e.g. x \<- Norm(mean = 2, sd = 1) and the above functions
  can be used generically (e.g. d(x) is the density function of the
  N(2, 1) distribution).
- Distributions not included in base R are made available, such as the
  Dirichlet and the Multivariate Gamma.
- Parameter estimation is performed analytically instead of numerically
  for the estimators that can be expressed explicitly.
- Numerical optimization of the MLE (whenever required, e.g. the Beta
  and Gamma distributions) is performed with computational efficiency,
  taking advantage of the score equation system to reduce the
  dimensionality of the optimization.
- Functions to compute and plot common estimator metrics (bias,
  variance, and RMSE) are included in the package to allow the
  convenient study and comparison of the estimators.

## Installation

You can install the release version of `estim` from CRAN by running the
following line of code:

``` r
 install.packages("estim")
```

You can install the development version of `estim` from github by
running the following line of code:

``` r
 devtools::install_github("thechibo/estim")
```

More details can be found in the [estim Github
repository](https://github.com/thechibo/estim "estim Github repository").

## Documentation

Detailed documentation, along with reproducible examples, can be found
in the package vignette `vignette(topic = "estim", package = "estim")`.

## Team

The `estim` package is developed in the [Mathematics
Department](https://en.math.uoa.gr/ "Mathematics Department Homepage")
of the [University of
Athens](https://en.uoa.gr/ "University of Athens Homepage"). The package
maintainer is [Ioannis
Oikonomidis](http://users.uoa.gr/~goikon/ "Ioannis Oikonomidis Homepage"),
working under the supervision of [Prof. Samis
Trevezas](http://scholar.uoa.gr/strevezas/ "Samis Trevezas Homepage").
