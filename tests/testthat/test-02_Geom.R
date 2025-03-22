test_that("Geom distr works", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Geom")

  # Errors
  expect_error(Geom(1:2))
  expect_error(Geom(-1))
  expect_error(Geom(4))

})

test_that("Geom dpqr work", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)
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
  expect_identical(d(D)(0), p)
  expect_identical(p(D)(0), p)
  expect_identical(p(D)(Inf), 1)
  expect_identical(qn(D)(1), Inf)
  expect_identical(qn(D)(0), 0)
  expect_identical(sum(r(D)(n) >= 0), n)

  # 2-Way Calls
  expect_identical(d(D)(1), dgeom(1, p))
  expect_identical(p(D)(1), pgeom(1, p))
  expect_equal(qn(D)(0.5), qgeom(0.5, p), tolerance = 0.01)

})

test_that("Geom moments work", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)

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

test_that("Geom likelihood works", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llgeom(x, p)))

  # 2-Way Calls
  expect_identical(llgeom(x, p), ll(x, p, D))

})

test_that("Geom estim works", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(eexp(x, type = "mle")))
  expect_true(is.numeric(eexp(x, type = "me")))

  # 2-Way Calls
  expect_identical(egeom(x, type = "mle"), estim(x, D, type = "mle"))
  expect_identical(egeom(x, type = "me"), estim(x, D, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)

})

test_that("Geom avar works", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)

  # Types
  expect_true(is.numeric(vgeom(p, type = "mle")))
  expect_true(is.numeric(vgeom(p, type = "me")))

  # 2-Way Calls
  expect_identical(vgeom(p, type = "mle"), avar(D, type = "mle"))
  expect_identical(vgeom(p, type = "me"), avar(D, type = "me"))
  expect_identical(vgeom(p, type = "mle"), avar_mle(D))
  expect_identical(vgeom(p, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)

})

test_that("Geom small metrics work", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)
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

test_that("Geom large metrics work", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)

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
