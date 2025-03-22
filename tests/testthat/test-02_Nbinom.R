test_that("Nbinom distr works", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Nbinom")

  # Errors
  expect_error(Nbinom(-10, 0.5))
  expect_error(Nbinom(0, 0.5))
  expect_error(Nbinom(10, 5))
  expect_error(Nbinom(3:4, 0.5))
  expect_error(Nbinom(10, c(0.5, 0.6)))

})

test_that("Nbinom dpqr work", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_equal(d(D)(0), p ^ k, tolerance = 0.01)
  expect_equal(d(D)(1), k * (1 - p) * p ^ k, tolerance = 0.01)
  expect_identical(p(D)(-1), 0)
  expect_identical(qn(D)(1), Inf)
  expect_identical(qn(D)(0), 0)
  expect_identical(sum(x >= 0), n)

  # 2-Way Calls
  expect_identical(d(D)(1), dnbinom(1, k, p))
  expect_identical(p(D)(1), pnbinom(1, k, p))
  expect_identical(qn(D)(1), qnbinom(1, k, p))
  expect_identical(qn(D)(0), qnbinom(0, k, p))

  # Special Case: Geom
  D1 <- Nbinom(1, 0.7)
  D2 <- Geom(0.7)
  expect_identical(d(D1)(3), d(D2)(3))
  expect_identical(p(D1)(3), p(D2)(3))
  expect_identical(qn(D1)(0.7), qn(D2)(0.7))

})

test_that("Nbinom moments work", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)

  # Types
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(median(D)))
  expect_true(is.numeric(mode(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(sd(D)))
  expect_true(is.numeric(skew(D)))
  expect_true(is.numeric(kurt(D)))
  expect_true(is.numeric(entro(D)))
  expect_true(is.numeric(finf(D)))

  # Values
  expect_equal(mean(D), k * (1 / p - 1), tolerance = 0.01)

})

test_that("Nbinom likelihood works", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llnbinom(x, size = k, prob = p)))

  # 2-Way Calls
  expect_identical(llnbinom(x, k, p), ll(x, c(k, p), D))

})

test_that("Nbinom estim works", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(enbinom(x, k, type = "mle")))
  expect_true(is.numeric(enbinom(x, k, type = "me")))

  # 2-Way Calls
  expect_identical(enbinom(x, k, type = "mle"), estim(x, D, type = "mle"))
  expect_identical(enbinom(x, k, type = "me"), estim(x, D, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true["prob"], d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true["prob"], d$prm_est, tolerance = 0.02)

})

test_that("Nbinom avar works", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(vnbinom(k, p, type = "mle")))
  expect_true(is.numeric(vnbinom(k, p, type = "me")))

  # 2-Way Calls
  expect_identical(vnbinom(k, p, type = "mle"), avar(D, type = "mle"))
  expect_identical(vnbinom(k, p, type = "me"), avar(D, type = "me"))
  expect_identical(vnbinom(k, p, type = "mle"), avar_mle(D))
  expect_identical(vnbinom(k, p, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.1)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.1)

})

test_that("Nbinom small metrics work", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)

  prm <- list(name = "prob",
              pos = NULL,
              val = seq(0.5, 0.8, by = 0.1))

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

test_that("Nbinom large metrics work", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)

  prm <- list(name = "prob",
              pos = NULL,
              val = seq(0.5, 0.8, by = 0.1))

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
