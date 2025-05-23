# estimators 0.11.3
- This package is now deprecated and will be archived.
- Please migrate to [joker](https://github.com/thechibo/joker)

# estimators 0.11.2

* Added distribution-specific documentation.
* Removed pkgdown due to incompatibility with Sweave.

# estimators 0.11.1

* Added pkgdown check.
* Attempted to fix compact-vignettes bug.

# estimators 0.11.0

* Changed plot_small_metrics, plot_large_metrics to plot
* Updated documentation pages.
* Added methods for d(D, x), p(D, x) etc.
* Added Multigam and Weibull distributions.
* Changed estim() to e().
* Changed estimators to be lists instead of atomic vectors.
* Created a more accurate idigamma() function.

# estimators 0.10.3

* Added new JSS vignette.

# estimators 0.10.2

* Tried to update test-coverage.
* Fixed bug in Nbinom test.

# estimators 0.10.1

* Removed false documentation file that caused R CMD check to fail.

# estimators 0.10.0

* Added tests for all distributions.
* Added many functions that were not available for several distributions.
* Corrected multiple bugs.
* Removed Ncchisq, Ncfisher. They will be added in later versions.

# estimators 0.9.0

* Changed package name to estim.
* Fixed Fisher df bug.
* Added NcFisher.

# estimators 0.8.5

* Removed the "v function works" Dirichlet test that failed on MKL until the reason of failure becomes clear.

# estimators 0.8.3

* Updated missing return value in calculus.Rd documentation.
* Changed default plot save path to `NULL`.
* Changed `q2` function to `qn`.

# estimators 0.8.2

* Updated missing return values in documentation.

# estimators 0.8.1

* Added support for S4 Distribution classes, therefore the distr package is no longer a dependency (nor is it imported).
* Added moment functions for all distributions.
* Added new distributions.

# estimators 0.7.3

* Removed dontrun from all examples.
* Resubmitted to CRAN.

# estimators 0.7.2

* Updated DESCRIPTION.
* Added references.
* Removed dontrun from examples.
* Resubmitted to CRAN.

# estimators 0.7.1

* Updated DESCRIPTION.
* Submitted to CRAN.

# estimators 0.7.0

* Created the estimators and avar functions that cover all distributions and estimation methods.
* Added ll<distrname> functions that calculate the log-likelihood to provide support in the default R style.
* Added e<distrname> functions that perform parameter estimation to provide support in the default R style.
* Added v<distrname> functions that calculate the estimator asymptotic variance to provide support in the default R style.
* Removed the me2 and same2 estimators. Added "dirich" argument to the MGamma me and same estimators.
* Moved package distr from Depends to Imports.
* Added testthat coverage.
* Ran rhub and devtools checks for cran release.
* Complied by goodpractice::gp().

# estimators 0.6.0

* Added MGamma ME, SAME estimators.
* Upgraded acov, plot_acov to include covariances.

# estimators 0.5.2

* Fixed a bug in dependencies.

# estimators 0.5.1

* Fixed a bug in Description.
* Enriched the documentation.
* Complied by goodpractice::gp().

# estimators 0.5.0

* Added a logo.
* Added metrics and plot_metrics.
* Added acov and plot_acov.
* Enriched the documentation.

# estimators 0.4.0

* Added GitHub and documentation files.
