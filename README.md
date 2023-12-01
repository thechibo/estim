
<!-- README.md is generated from README.Rmd. Please edit that file -->

# estimators <img src=man/figures/logo.png align="right" height="139" alt="logo"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/estimators)](https://CRAN.R-project.org/package=estimators)
[![R-CMD-check](https://github.com/thechibo/estimators/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/thechibo/estimators/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/thechibo/estimators/branch/main/graph/badge.svg)](https://app.codecov.io/gh/thechibo/estimators?branch=main)
<!-- badges: end -->

## Introduction

The `estimators` R package performs parameter estimation in common
distribution families, making moment, maximum likelihood, and
state-of-the-art estimators more accessible.

### Key Features

- The common d, p, q, r function family for each distribution
  (e.g. dnorm, pnorm, qnorm, rnorm) is enriched with
  - the ll counterpart (e.g. llnorm) that calculates the log-likelihood,
  - the e counterpart (e.g. enorm) that performs parameter estimation,
  - the v counterpart (e.g. vnorm) that calculates the asymptotic
    variance-covariance matrix of an estimators.
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
- All functions can be used with the S4-Distribution system developed by
  the `distr` package family.

## Installation

You can install the release version of `estimators` from CRAN by running
the following line of code:

``` r
 install.packages("estimators")
```

You can install the development version of `estimators` from github by
running the following line of code:

``` r
 devtools::install_github("thechibo/estimators")
```

More details can be found in the [estimators Github
repository](https://github.com/thechibo/estimators "estimators Github repository").

## Documentation

Detailed documentation, along with reproducible examples, can be found
in the package vignette
`vignette(topic = "estimators", package = "estimators")`.

## Team

The `estimators` package is developed in the [Mathematics
Department](https://en.math.uoa.gr/ "Mathematics Department Homepage")
of the [University of
Athens](https://en.uoa.gr/ "University of Athens Homepage"). The package
maintainer is [Ioannis
Oikonomidis](http://users.uoa.gr/~goikon/ "Ioannis Oikonomidis Homepage"),
working under the supervision of [Prof. Samis
Trevezas](http://scholar.uoa.gr/strevezas/ "Samis Trevezas Homepage").
