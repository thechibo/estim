test_that("Pois distr works", {

  # Preliminaries
  lambda <- 3
  D <- Pois(lambda)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Pois")

  # Errors
  expect_error(Pois(1:2))
  expect_error(Pois(-1))

})

test_that("Pois dpqr work", {

  # Preliminaries
  lambda <- 3
  D <- Pois(lambda)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(d(D)(-1), 0)
  expect_warning(d(D)(1.5))
  expect_identical(p(D)(-1), 0)
  expect_identical(p(D)(Inf), 1)
  expect_identical(qn(D)(1), Inf)
  expect_identical(qn(D)(0), 0)
  expect_identical(sum(r(D)(n) >= 0), n)

  # 2-Way Calls
  expect_identical(d(D)(1), dpois(1, lambda))
  expect_identical(p(D)(1), ppois(1, lambda))
  expect_equal(qn(D)(0.5), qpois(0.5, lambda), tolerance = 0.01)

})

test_that("Pois moments work", {

  # Preliminaries
  lambda <- 3
  D <- Pois(lambda)

  # Types
  expect_true(is.numeric(mean(D)))
  expect_warning(median(D))
  expect_true(is.numeric(mode(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(sd(D)))
  expect_true(is.numeric(skew(D)))
  expect_true(is.numeric(kurt(D)))
  expect_warning(entro(D))
  expect_true(is.numeric(finf(D)))

})

test_that("Pois likelihood works", {

  # Preliminaries
  lambda <- 3
  D <- Pois(lambda)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llpois(x, lambda)))

  # 2-Way Calls
  expect_identical(llpois(x, lambda), ll(x, lambda, D))

})

test_that("Pois estim works", {

  # Preliminaries
  lambda <- 3
  D <- Pois(lambda)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(epois(x, type = "mle")))
  expect_true(is.numeric(epois(x, type = "me")))

  # 2-Way Calls
  expect_identical(epois(x, type = "mle"), estim(x, D, type = "mle"))
  expect_identical(epois(x, type = "me"), estim(x, D, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)

})

test_that("Pois avar works", {

  # Preliminaries
  lambda <- 3
  D <- Pois(lambda)

  # Types
  expect_true(is.numeric(vpois(lambda, type = "mle")))
  expect_true(is.numeric(vpois(lambda, type = "me")))

  # 2-Way Calls
  expect_identical(vpois(lambda, type = "mle"), avar(D, type = "mle"))
  expect_identical(vpois(lambda, type = "me"), avar(D, type = "me"))
  expect_identical(vpois(lambda, type = "mle"), avar_mle(D))
  expect_identical(vpois(lambda, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)

})

test_that("Pois small metrics work", {

  # Preliminaries
  lambda <- 3
  D <- Pois(lambda)
  set.seed(1)

  prm <- list(name = "lambda",
              pos = NULL,
              val = seq(0.5, 5, by = 0.5))

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

test_that("Pois large metrics work", {

  # Preliminaries
  lambda <- 3
  D <- Pois(lambda)

  prm <- list(name = "lambda",
              pos = NULL,
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
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
