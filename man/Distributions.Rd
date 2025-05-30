% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_generics.R
\name{Distributions}
\alias{Distributions}
\alias{d}
\alias{p}
\alias{q}
\alias{r}
\alias{qn}
\title{Distribution S4 Classes}
\usage{
d(distr, x, ...)

p(distr, x, ...)

qn(distr, x, ...)

r(distr, n, ...)
}
\arguments{
\item{distr}{a \code{Distribution} object.}

\item{x}{numeric. The point to evaluate the function.}

\item{...}{extra arguments.}

\item{n}{numeric. The sample size.}
}
\value{
Each type of function returns a different type of object:
\itemize{
\item Distribution Functions: When supplied with one argument (\code{distr}), the
\code{d()}, \code{p()}, \code{q()}, \code{r()}, \code{ll()} functions return the density, cumulative
probability, quantile, random sample generator, and log-likelihood functions,
respectively. When supplied with both arguments (\code{distr} and \code{x}), they
evaluate the aforementioned functions directly.
\item Moments: Returns a numeric, either vector or matrix depending on the moment
and the distribution. The \code{moments()} function returns a list with all the
available methods.
\item Estimation: Returns a list, the estimators of the unknown parameters. Note
that in distribution families like the binomial, multinomial, and negative
binomial, the size is not returned, since it is considered known.
\item Variance: Returns a named matrix. The asymptotic covariance matrix of the
estimator.
}
}
\description{
A collection of S4 classes that provide a flexible and structured way to work
with probability distributions.
}
\details{
These S4 generic methods can work both as functions and as functionals
(functions that return functions). The available distribution families are
coded as S4 classes, specifically subclasses of the \code{Distribution}
superclass. The methods can be used in two ways:

Option 1: If both the \code{distr} argument and \code{x} or \code{n} are supplied, then the
function is evaluated directly, as usual.

Option 2: If only the \code{distr} argument is supplied, the method returns a
function that takes as input the missing argument \code{x} or \code{n}, allowing the
user to work with the function object itself. See examples.

Looking for a specific distribution family?
This help page is general. Use the help page of each distribution to see the
available methods for the class, details, and examples. Check the See Also
section.
}
\section{Functions}{
\itemize{
\item \code{d()}: density function

\item \code{p()}: cumulative distribution function

\item \code{qn()}: generalized inverse distribution function

\item \code{r()}: random sample generator function

}}
\examples{
# -----------------------------------------------------
# Beta Distribution Example
# -----------------------------------------------------

# Create the distribution
a <- 3 ; b <- 5
D <- Beta(a, b)
x <- c(0.3, 0.8, 0.5)
n <- 100

# ------------------
# dpqr Functions
# ------------------

d(D, x) # density function
p(D, x) # distribution function
qn(D, 0.8) # inverse distribution function
x <- r(D, n) # random generator function

# alternative way to use the function
df <- d(D) ; df(x) # df is a function itself

# ------------------
# Moments
# ------------------

mean(D) # Expectation
var(D) # Variance
sd(D) # Standard Deviation
skew(D) # Skewness
kurt(D) # Excess Kurtosis
entro(D) # Entropy
finf(D) # Fisher Information Matrix

# List of all available moments
mom <- moments(D)
mom$mean # expectation

# ------------------
# Point Estimation
# ------------------

ll(D, x)
llbeta(x, a, b)

ebeta(x, type = "mle")
ebeta(x, type = "me")
ebeta(x, type = "same")

mle(D, x)
me(D, x)
same(D, x)
e(D, x, type = "mle")

mle("beta", x) # the distr argument can be a character

# ------------------
# As. Variance
# ------------------

vbeta(a, b, type = "mle")
vbeta(a, b, type = "me")
vbeta(a, b, type = "same")

avar_mle(D)
avar_me(D)
avar_same(D)

avar(D, type = "mle")
}
\seealso{
\link{moments}, \link{loglikelihood}, \link{estimation}, \link{Bern},
\link{Beta}, \link{Binom}, \link{Cat}, \link{Cauchy}, \link{Chisq}, \link{Dir}, \link{Exp},
\link{Fisher}, \link{Gam}, \link{Geom}, \link{Laplace}, \link{Lnorm}, \link{Multigam}, \link{Multinom},
\link{Nbinom}, \link{Norm}, \link{Pois}, \link{Stud}, \link{Unif}, \link{Weib}
}
