test_that("Log-Likelihood works", {

  set.seed(1203)
  shape <- 1
  scale <- 2
  prm <- c(shape, scale)
  D <- Gam(shape = shape, scale = scale)
  x <- rgamma(100, shape, scale)

  expect_identical(llgamma(x, shape, scale), ll(x, prm, D))

})

test_that("e functions work", {

  set.seed(1203)
  shape <- 1
  scale <- 2
  prm <- c(shape, scale)
  D <- Gam(shape = shape, scale = scale)
  x <- rgamma(100, shape, scale)

  expect_identical(egamma(x, "mle"), mle(x, D))
  expect_identical(egamma(x, "me"), me(x, D))
  expect_identical(egamma(x, "same"), same(x, D))

})

test_that("v functions work", {

  shape <- 1
  scale <- 2
  prm <- c(shape, scale)
  D <- Gam(shape = shape, scale = scale)

  expect_identical(vgamma(shape, scale, "mle"), avar_mle(D))
  expect_identical(vgamma(shape, scale, "me"), avar_me(D))
  expect_identical(vgamma(shape, scale, "same"), avar_same(D))

})

test_that("ME is consistent", {

  set.seed(1203)
  est <- "me"
  D0 <- Gam()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("SAME is consistent", {

  set.seed(1203)
  est <- "same"
  D0 <- Gam()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("MLE is consistent", {

  set.seed(1203)
  est <- "mle"
  D0 <- Gam()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("ME avar is correct", {

  set.seed(1203)
  est <- "me"
  D0 <- Gam()
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 1)

})

test_that("SAME avar is correct", {

  set.seed(1203)
  est <- "same"
  D0 <- Gam()
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 1)

})

test_that("MLE avar is correct", {

  set.seed(1203)
  est <- "mle"
  D0 <- Gam()
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 1)

})

test_that("small_metrics works", {

  set.seed(1203)
  D <- Gam(shape = 1, scale = 2)

  prm <- list(name = "shape",
              pos = NULL,
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

  set.seed(1203)
  D <- Gam(shape = 1, scale = 2)

  prm <- list(name = "shape",
              pos = NULL,
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
