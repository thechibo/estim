test_that("Bern distr works", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Bern")

  # Errors
  expect_error(Bern(c(0.1, 0.2)))
  expect_error(Bern(2))
  expect_error(Bern(-1))

})

test_that("Bern dpqr work", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(d(D)(1), p)
  expect_identical(p(D)(1), 1)
  expect_identical(qn(D)(1), 1)
  expect_identical(qn(D)(0), 0)
  expect_identical(sum(x %in% c(0, 1)), n)

  # 2-Way Calls
  expect_identical(d(D)(1), dbern(1, p))
  expect_identical(p(D)(1), pbern(1, p))
  expect_identical(qn(D)(1), qbern(1, p))
  expect_identical(qn(D)(0), qbern(0, p))

})

test_that("Bern moments work", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)

  # Types
  expect_true(is.list(moments(D)))
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
  expect_identical(mean(D), p)
  expect_identical(var(D), p * (1 - p))

})

test_that("Bern likelihood works", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llbern(x, p)))

  # 2-Way Calls
  expect_identical(llbern(x, p), ll(x, p, D))

})

test_that("Bern estim works", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(ebern(x, type = "mle")))
  expect_true(is.numeric(ebern(x, type = "me")))

  # 2-Way Calls
  expect_identical(ebern(x, type = "mle"), estim(x, D, type = "mle"))
  expect_identical(ebern(x, type = "me"), estim(x, D, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true["prob"], d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true["prob"], d$prm_est, tolerance = 0.01)

})

test_that("Bern avar works", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(vbern(p, type = "mle")))
  expect_true(is.numeric(vbern(p, type = "me")))

  # 2-Way Calls
  expect_identical(vbern(p, type = "mle"), avar(D, type = "mle"))
  expect_identical(vbern(p, type = "me"), avar(D, type = "me"))
  expect_identical(vbern(p, type = "mle"), avar_mle(D))
  expect_identical(vbern(p, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est["prob"], tolerance = 0.01)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est["prob"], tolerance = 0.01)

})

test_that("Bern small metrics work", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
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

test_that("Bern large metrics work", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
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
