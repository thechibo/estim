test_that("Log-Likelihood works", {

  set.seed(1203)
  mean <- 1
  sd <- 2
  prm <- c(mean, sd)
  D <- Norm(mean = mean, sd = sd)
  x <- rnorm(100, mean, sd)

  expect_identical(llnorm(x, mean, sd), ll(x, prm, D))

})

test_that("e functions work", {

  set.seed(1203)
  mean <- 1
  sd <- 2
  prm <- c(mean, sd)
  D <- Norm(mean = mean, sd = sd)
  x <- rnorm(100, mean, sd)

  expect_identical(enorm(x, "mle"), mle(x, D))
  expect_identical(enorm(x, "me"), me(x, D))

})

test_that("v functions work", {

  mean <- 1
  sd <- 2
  prm <- c(mean, sd)
  D <- Norm(mean = mean, sd = sd)

  expect_identical(vnorm(mean, sd, "mle"), avar_mle(D))
  expect_identical(vnorm(mean, sd, "me"), avar_me(D))

})

test_that("ME is consistent", {

  set.seed(1203)
  est <- "me"
  D0 <- Norm()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("MLE is consistent", {

  set.seed(1203)
  est <- "mle"
  D0 <- Norm()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("ME avar is correct", {

  set.seed(1203)
  est <- "me"
  D0 <- Norm()
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("MLE avar is correct", {

  set.seed(1203)
  est <- "mle"
  D0 <- Norm(mean = 0, sd = 2)
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("small_metrics works", {

  set.seed(1203)
  D <- Norm(mean = 1, sd = 2)

  prm <- list(name = "sd",
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

  set.seed(1203)
  D <- Norm(mean = 1, sd = 2)

  prm <- list(name = "sd",
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
