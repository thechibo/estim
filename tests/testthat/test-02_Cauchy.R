test_that("Cauchy distr works", {

  # Preliminaries
  D <- Cauchy(0.7)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Cauchy")

  # Errors
  expect_error(Cauchy(c(0.1, 0.2)))
  expect_error(Cauchy(2, - 1))

})

test_that("Cauchy dpqr work", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- get_params_list(D)$location
  g <- get_params_list(D)$scale

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(d(D)(m), 1 / (pi * g))
  expect_identical(p(D)(m), 0.5)
  expect_identical(qn(D)(0.5), m)
  expect_identical(qn(D)(0), -Inf)
  expect_identical(qn(D)(1), Inf)

  # 2-Way Calls
  expect_identical(d(D)(1), dcauchy(1, m, g))
  expect_identical(p(D)(1), pcauchy(1, m, g))
  expect_identical(qn(D)(1), qcauchy(1, m, g))
  expect_identical(qn(D)(0), qcauchy(0, m, g))

})

test_that("Cauchy moments work", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- get_params_list(D)$location
  g <- get_params_list(D)$scale

  # Types
  expect_warning(mean(D))
  expect_true(is.numeric(median(D)))
  expect_true(is.numeric(mode(D)))
  expect_warning(var(D))
  expect_warning(sd(D))
  expect_warning(skew(D))
  expect_warning(kurt(D))
  expect_true(is.numeric(entro(D)))
  expect_true(is.numeric(finf(D)))

})

test_that("Cauchy likelihood works", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- get_params_list(D)$location
  g <- get_params_list(D)$scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llcauchy(x, m, g)))

  # 2-Way Calls
  expect_identical(llcauchy(x, m, g), ll(x, c(m, g), D))

  # ll and lloptim convergence comparison
  method <- "L-BFGS-B"
  lower <- c(-Inf, 1e-5)
  upper <- c(Inf, Inf)

  par1 <- optim(par = me(x, D),
                fn = lloptim,
                gr = dlloptim,
                tx = x,
                distr = D,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  par2 <- optim(par = me(x, D),
                fn = function(par, x, distr) { ll(x, par, distr) },
                x = x,
                distr = D,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  expect_equal(par1, par2, tolerance = 0.01)

})

test_that("Cauchy estim works", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- get_params_list(D)$location
  g <- get_params_list(D)$scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(ecauchy(x, type = "mle")))
  expect_true(is.numeric(ecauchy(x, type = "me")))

  # 2-Way Calls
  expect_identical(ecauchy(x, type = "mle"), estim(x, D, type = "mle"))
  expect_identical(ecauchy(x, type = "me"), estim(x, D, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

})

## NEEDS CHECKING
test_that("Cauchy avar works", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- get_params_list(D)$location
  g <- get_params_list(D)$scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(vcauchy(m, g, type = "mle")))

  # 2-Way Calls
  expect_identical(vcauchy(m, g, type = "mle"), avar(D, type = "mle"))
  expect_identical(vcauchy(m, g, type = "mle"), avar_mle(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.1)

})

test_that("Cauchy small metrics work", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- get_params_list(D)$location
  g <- get_params_list(D)$scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "location",
              pos = NULL,
              val = seq(-2, 2, by = 1))

  expect_no_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me"),
                       obs = c(20, 50),
                       sam = 1e2,
                       seed = 1)
  )

  expect_no_error(
    plot_small_metrics(x,
                       save = TRUE,
                       path = tempdir())
  )

  # Types
  expect_s3_class(x, "data.frame")
  expect_true(is.numeric(x$Parameter))
  expect_s3_class(x$Observations, "factor")
  expect_s3_class(x$Estimator, "factor")
  expect_s3_class(x$Metric, "factor")
  expect_true(is.numeric(x$Value))

})

## NEEDS CHECKING
test_that("Cauchy large metrics work", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- get_params_list(D)$location
  g <- get_params_list(D)$scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "location",
              pos = NULL,
              val = seq(-2, 2, by = 1))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle"))
  )

  expect_no_error(
    plot_large_metrics(x,
                       save = TRUE,
                       path = tempdir())
  )

  # Types
  expect_s3_class(x, "data.frame")
  expect_true(is.numeric(x$Parameter))
  expect_s3_class(x$Estimator, "factor")
  expect_true(is.numeric(x$Value))

  prm <- list(name = "scale",
              pos = NULL,
              val = seq(0.5, 2, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle"))
  )

  expect_no_error(
    plot_large_metrics(x,
                       save = TRUE,
                       path = tempdir())
  )

  # Types
  expect_s3_class(x, "data.frame")
  expect_true(is.numeric(x$Parameter))
  expect_s3_class(x$Estimator, "factor")
  expect_true(is.numeric(x$Value))

})
