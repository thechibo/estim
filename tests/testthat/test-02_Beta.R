test_that("Log-Likelihood works", {

  set.seed(1203)
  shape1 <- 1
  shape2 <- 2
  prm <- c(shape1, shape2)
  D <- Beta(shape1 = shape1, shape2 = shape2)
  x <- rbeta(100, shape1, shape2)

  expect_identical(llbeta(x, shape1, shape2), ll(x, prm, D))

})

test_that("e functions work", {

  set.seed(1203)
  shape1 <- 1
  shape2 <- 2
  prm <- c(shape1, shape2)
  D <- Beta(shape1 = shape1, shape2 = shape2)
  x <- rbeta(100, shape1, shape2)

  expect_identical(ebeta(x, "mle"), mle(x, D))
  expect_identical(ebeta(x, "me"), me(x, D))
  expect_identical(ebeta(x, "same"), same(x, D))

})

test_that("v functions work", {

  shape1 <- 1
  shape2 <- 2
  prm <- c(shape1, shape2)
  D <- Beta(shape1 = shape1, shape2 = shape2)

  expect_identical(vbeta(shape1, shape2, "mle"), avar_mle(D))
  expect_identical(vbeta(shape1, shape2, "me"), avar_me(D))
  expect_identical(vbeta(shape1, shape2, "same"), avar_same(D))

})

test_that("ME is consistent", {

  set.seed(1203)
  est <- "me"
  D0 <- Beta()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("SAME is consistent", {

  set.seed(1203)
  est <- "same"
  D0 <- Beta()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("MLE is consistent", {

  set.seed(1203)
  est <- "mle"
  D0 <- Beta()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("ME avar is correct", {

  set.seed(1203)
  est <- "me"
  D0 <- Beta()
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 1)

})

test_that("SAME avar is correct", {

  set.seed(1203)
  est <- "same"
  D0 <- Beta()
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 1)

})

test_that("MLE avar is correct", {

  set.seed(1203)
  est <- "mle"
  D0 <- Beta()
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 1)

})

test_that("small_metrics works", {

  set.seed(1203)
  D <- Beta(shape1 = 1, shape2 = 2)

  prm <- list(name = "shape1",
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
  D <- Beta(shape1 = 1, shape2 = 2)

  prm <- list(name = "shape1",
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
