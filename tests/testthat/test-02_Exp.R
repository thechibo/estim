test_that("Log-Likelihood works", {

  rate <- 2
  prm <- rate
  D <- distr::Exp(rate)
  x <- rexp(100, rate)

  expect_identical(llexp(x, rate), ll(x, prm, D))

})

test_that("e functions work", {

  rate <- 2
  prm <- rate
  D <- distr::Exp(rate)
  x <- rexp(100, rate)

  expect_identical(eexp(x, "mle"), mle(x, D))
  expect_identical(eexp(x, "me"), me(x, D))

})

test_that("v functions work", {

  rate <- 2
  prm <- rate
  D <- distr::Exp(rate)

  expect_identical(vexp(rate, "mle"), avar_mle(D))
  expect_identical(vexp(rate, "me"), avar_me(D))

})

test_that("ME is consistent", {

  est <- "me"
  D0 <- distr::Exp()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("MLE is consistent", {

  est <- "mle"
  D0 <- distr::Exp()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("ME avar is correct", {

  est <- "me"
  D0 <- distr::Exp()
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 1)

})

test_that("MLE avar is correct", {

  est <- "mle"
  D0 <- distr::Exp()
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 1)

})

test_that("small_metrics works", {

  D <- distr::Exp(rate = 3)

  prm <- list(name = "rate",
              pos = NULL,
              val = seq(0.5, 2, by = 0.5))

  expect_no_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me"),
                       obs = c(20, 50),
                       sam = 1e2,
                       seed = 1)
  )
  expect_s3_class(x, "data.frame")

  expect_no_error(
    plot_small_metrics(x,
                       save = TRUE,
                       path = tempdir())
  )

})

test_that("large_metrics works", {

  D <- distr::Exp(rate = 3)

  prm <- list(name = "rate",
              pos = NULL,
              val = seq(0.5, 2, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
  )

  expect_s3_class(x, "data.frame")

  expect_no_error(
    plot_large_metrics(x,
                       save = TRUE,
                       path = tempdir())
  )

})
