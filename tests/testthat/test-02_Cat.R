test_that("Cat distr works", {

  # Preliminaries
  p <- c(0.1, 0.2, 0.3, 0.4)
  D <- Cat(p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Cat")

  # Errors
  expect_error(Cat(2))
  expect_error(Cat(-1))

})

test_that("Cat dpqr work", {

  # Preliminaries
  p <- c(0.1, 0.2, 0.3, 0.4)
  D <- Cat(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(d(D)(c(1, 0, 0, 0)), c(p[1], 0, 0, 0))
  expect_identical(sum(x %in% 1:4), n)

  # 2-Way Calls
  expect_identical(d(D)(1), dcat(1, p))
  expect_identical(d(D)(1), d(D, 1))

})

test_that("Cat moments work", {

  # Preliminaries
  p <- c(0.1, 0.2, 0.3, 0.4)
  D <- Cat(p)

  # Types
  expect_true(is.list(moments(D)))
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(mode(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(finf(D)))

  # Values
  expect_identical(mean(D), p)

})

test_that("Cat likelihood works", {

  # Preliminaries
  p <- c(0.1, 0.2, 0.3, 0.4)
  D <- Cat(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llcat(x, p)))

  # 2-Way Calls
  expect_identical(llcat(x, p), ll(D, x))
  expect_identical(ll(D)(x), ll(D, x))

})

test_that("Cat estim works", {

  # Preliminaries
  p <- c(0.1, 0.2, 0.3, 0.4)
  D <- Cat(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(ecat(x, dim = 4, type = "mle")))
  expect_true(is.list(ecat(x, dim = 4, type = "me")))

  # 2-Way Calls
  expect_identical(ecat(x, type = "mle", dim = 4), e(D, x, type = "mle"))
  expect_identical(ecat(x, type = "me", dim = 4), e(D, x, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

})

test_that("Cat avar works", {

  # Preliminaries
  p <- c(0.1, 0.2, 0.3, 0.4)
  D <- Cat(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)
  k <- length(p)

  # Types
  expect_true(is.numeric(vcat(p, type = "mle")))
  expect_true(is.numeric(vcat(p, type = "me")))

  # 2-Way Calls
  expect_identical(vcat(p, type = "mle"), avar(D, type = "mle"))
  expect_identical(vcat(p, type = "me"), avar(D, type = "me"))
  expect_identical(vcat(p, type = "mle"), avar_mle(D))
  expect_identical(vcat(p, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est[1:(k-1), 1:(k-1)], tolerance = 0.05)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est[1:(k-1), 1:(k-1)], tolerance = 0.05)

})

test_that("Cat small metrics work", {

  # Preliminaries
  p <- c(0.1, 0.2, 0.3, 0.4)
  D <- Cat(p)

  prm <- list(name = "prob",
              val = seq(0.5, 0.8, by = 0.1))

  expect_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me"),
                       obs = c(20, 50),
                       sam = 1e2,
                       seed = 1)
  )

})

test_that("Cat large metrics work", {

  # Preliminaries
  p <- c(0.1, 0.2, 0.3, 0.4)
  D <- Cat(p)

  prm <- list(name = "prob",
              val = seq(0.5, 0.8, by = 0.1))

  expect_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
  )

})
