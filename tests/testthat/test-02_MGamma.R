test_that("Log-Likelihood works", {

  shape <- 1:3
  scale <- 2
  prm <- c(shape, scale)
  D <- MGamma(shape, scale)
  x <- rmgamma(100, shape, scale)

  expect_identical(llmgamma(x, shape, scale), ll(x, prm, D))

})

test_that("e functions work", {

  shape <- 1:3
  scale <- 2
  prm <- c(shape, scale)
  D <- MGamma(shape, scale)
  x <- rmgamma(100, shape, scale)

  expect_identical(emgamma(x, "mle"), mle(x, D))
  expect_identical(emgamma(x, "me"), me(x, D))
  expect_identical(emgamma(x, "same"), same(x, D))

})

test_that("v functions work", {

  shape <- 1:3
  scale <- 2
  prm <- c(shape, scale)
  D <- MGamma(shape, scale)

  expect_identical(vmgamma(shape, scale, "mle"), avar_mle(D))
  expect_identical(vmgamma(shape, scale, "me"), avar_me(D))
  expect_identical(vmgamma(shape, scale, "same"), avar_same(D))

})

test_that("ME is consistent", {

  est <- "me"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("SAME is consistent", {

  est <- "same"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("Dirichlet-based ME is consistent", {

  est <- "me"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_consistency(est, D0, dirich = TRUE)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("Dirichlet-based SAME is consistent", {

  est <- "same"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_consistency(est, D0, dirich = TRUE)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("MLE is consistent", {

  est <- "mle"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("ME avar is correct", {

  est <- "me"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("SAME avar is correct", {

  est <- "same"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("Dirichlet-based ME avar is correct", {

  est <- "me"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_avar(est, D0, dirich = TRUE)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("Dirichlet-based SAME avar is correct", {

  est <- "same"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_avar(est, D0, dirich = TRUE)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("MLE avar is correct", {

  est <- "mle"
  D0 <- MGamma(shape = 1:3, scale = 4)
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("small_metrics works", {

  D <- MGamma(shape = 1:3, scale = 2)

  prm <- list(name = "shape",
              pos = 2,
              val = seq(0.5, 2, by = 0.5))

  expect_no_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me", "same"),
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

  D <- MGamma(shape = 1:3, scale = 2)

  prm <- list(name = "shape",
              pos = 2,
              val = seq(0.5, 2, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me", "same"))
  )

  expect_s3_class(x, "data.frame")

  expect_no_error(
    plot_large_metrics(x,
                       save = TRUE,
                       path = tempdir())
  )

})
