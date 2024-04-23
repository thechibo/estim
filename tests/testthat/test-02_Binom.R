test_that("Log-Likelihood works", {

  set.seed(1203)
  size <- 10
  prob <- 0.7
  prm <- c(size, prob)
  D <- Binom(size = size, prob = prob)
  x <- rbinom(100, size, prob)

  expect_identical(llbinom(x, size, prob), ll(x, prm, D))

})

test_that("e functions work", {

  set.seed(1203)
  size <- 10
  prob <- 0.7
  prm <- c(size, prob)
  D <- Binom(size = size, prob = prob)
  x <- rbinom(100, size, prob)

  expect_identical(ebinom(x, "mle"), mle(x, D))
  expect_identical(ebinom(x, "me"), me(x, D))

})

test_that("v functions work", {

  size <- 10
  prob <- 0.7
  prm <- c(size, prob)
  D <- Binom(size = size, prob = prob)

  expect_identical(vbinom(size, prob, "mle"), avar_mle(D))
  expect_identical(vbinom(size, prob, "me"), avar_me(D))

})

test_that("ME is consistent", {

  set.seed(1203)
  est <- "me"
  D0 <- Binom()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true["prob"]["prob"], d$prm_est, tolerance = 0.5)

})

test_that("MLE is consistent", {

  set.seed(1203)
  est <- "mle"
  D0 <- Binom()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true["prob"], d$prm_est, tolerance = 0.5)

})

test_that("ME avar is correct", {

  set.seed(1203)
  est <- "me"
  D0 <- Binom()
  d <- test_avar(est, D0)
  expect_equal(unname(d$avar_true), d$avar_est["prob", "prob"], tolerance = 1)

})

test_that("MLE avar is correct", {

  set.seed(1203)
  est <- "mle"
  D0 <- Binom()
  d <- test_avar(est, D0)
  expect_equal(unname(d$avar_true), d$avar_est["prob", "prob"], tolerance = 1)

})

test_that("small_metrics works", {

  set.seed(1203)
  D <- Binom()

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
  expect_s3_class(x, "data.frame")

  expect_no_error(
    plot_small_metrics(x,
                       save = TRUE,
                       path = tempdir())
  )

})

test_that("large_metrics works", {

  set.seed(1203)
  D <- Binom()

  prm <- list(name = "prob",
              pos = NULL,
              val = seq(0.5, 0.8, by = 0.1))

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
