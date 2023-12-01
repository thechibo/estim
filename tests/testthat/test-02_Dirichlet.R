test_that("Log-Likelihood works", {

  shape <- 1:4
  prm <- shape
  D <- Dirichlet(shape)
  x <- rdirich(100, shape)

  expect_identical(lldirich(x, shape), ll(x, prm, D))

})

test_that("e functions work", {

  shape <- 1:4
  prm <- shape
  D <- Dirichlet(shape)
  x <- rdirich(100, shape)

  expect_identical(edirich(x, "mle"), mle(x, D))
  expect_identical(edirich(x, "me"), me(x, D))
  expect_identical(edirich(x, "same"), same(x, D))

})

test_that("v functions work", {

  shape <- 1:4
  prm <- shape
  D <- Dirichlet(shape)

  expect_identical(vdirich(shape, "mle"), avar_mle(D))
  expect_identical(vdirich(shape, "me"), avar_me(D))
  expect_identical(vdirich(shape, "same"), avar_same(D))

})

test_that("ME is consistent", {

  est <- "me"
  D0 <- Dirichlet()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("SAME is consistent", {

  est <- "same"
  D0 <- Dirichlet()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("MLE is consistent", {

  est <- "mle"
  D0 <- Dirichlet()
  d <- test_consistency(est, D0)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.5)

})

test_that("ME avar is correct", {

  est <- "me"
  D0 <- Dirichlet(1:4)
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("SAME avar is correct", {

  est <- "same"
  D0 <- Dirichlet(1:4)
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("MLE avar is correct", {

  est <- "mle"
  D0 <- Dirichlet(1:4)
  d <- test_avar(est, D0)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.5)

})

test_that("small_metrics works", {

  D <- Dirichlet(1:4)

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

  D <- Dirichlet(1:4)

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
