test_that("Exp distr works", {

  # Preliminaries
  rate <- 3
  D <- Exp(rate)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Exp")

  # Errors
  expect_error(Exp(1:2))
  expect_error(Exp(-1))

})

test_that("Exp dpqr work", {

  # Preliminaries
  rate <- 3
  D <- Exp(rate)
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
  expect_identical(p(D)(0), 0)
  expect_identical(p(D)(Inf), 1)
  expect_identical(qn(D)(1), Inf)
  expect_identical(qn(D)(0), 0)
  expect_identical(sum(r(D)(n) > 0), n)

  # 2-Way Calls
  expect_identical(d(D)(1), dexp(1, rate))
  expect_identical(p(D)(1), pexp(1, rate))
  expect_equal(qn(D)(0.5), qexp(0.5, rate), tolerance = 0.01)

})

test_that("Exp moments work", {

  # Preliminaries
  rate <- 3
  D <- Exp(rate)

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

})

test_that("Exp likelihood works", {

  # Preliminaries
  rate <- 3
  D <- Exp(rate)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llexp(x, rate)))

  # 2-Way Calls
  expect_identical(llexp(x, rate), ll(x, rate, D))

})

test_that("Exp estim works", {

  # Preliminaries
  rate <- 3
  D <- Exp(rate)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(eexp(x, type = "mle")))
  expect_true(is.numeric(eexp(x, type = "me")))

  # 2-Way Calls
  expect_identical(eexp(x, type = "mle"), estim(x, D, type = "mle"))
  expect_identical(eexp(x, type = "me"), estim(x, D, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)

})

test_that("Exp avar works", {

  # Preliminaries
  rate <- 3
  D <- Exp(rate)

  # Types
  expect_true(is.numeric(vexp(rate, type = "mle")))
  expect_true(is.numeric(vexp(rate, type = "me")))

  # 2-Way Calls
  expect_identical(vexp(rate, type = "mle"), avar(D, type = "mle"))
  expect_identical(vexp(rate, type = "me"), avar(D, type = "me"))
  expect_identical(vexp(rate, type = "mle"), avar_mle(D))
  expect_identical(vexp(rate, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)

})

test_that("Exp small metrics work", {

  # Preliminaries
  rate <- 3
  D <- Exp(rate)
  set.seed(1)

  prm <- list(name = "rate",
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

test_that("Exp large metrics work", {

  # Preliminaries
  rate <- 3
  D <- Exp(rate)
  set.seed(1)

  prm <- list(name = "rate",
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
