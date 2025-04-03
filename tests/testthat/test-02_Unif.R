test_that("Unif distr works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Unif(a, b)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Unif")

  # Errors
  expect_error(Unif(c(0.1, 0.2, 0.3)))
  expect_error(Unif(-1, 2, 3))
  expect_error(Unif(1, -2))

})

test_that("Unif dpqr work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Unif(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(d(D)(a - 1), 0)
  expect_identical(d(D)((a + b) / 2), 1/ (b - a))
  expect_identical(p(D)(b), 1)
  expect_identical(p(D)(a), 0)
  expect_identical(qn(D)(1), b)
  expect_identical(qn(D)(0), a)
  expect_identical(sum(x <= b), n)
  expect_identical(sum(x >= a), n)

  # 2-Way Calls
  expect_identical(d(D)(2.5), dunif(2.5, a, b))
  expect_identical(p(D)(2.5), punif(2.5, a, b))
  expect_equal(qn(D)(0.4), qunif(0.4, a, b), tolerance = 1e-8)
  expect_identical(d(D)(2.5), d(D, 2.5))
  expect_identical(p(D)(2.5), p(D, 2.5))
  expect_equal(qn(D)(0.4), qn(D, 0.4), tolerance = 1e-8)

})

test_that("Unif moments work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Unif(a, b)

  # Types
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(median(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(sd(D)))
  expect_true(is.numeric(skew(D)))
  expect_true(is.numeric(kurt(D)))
  expect_true(is.numeric(entro(D)))

  # Warnings
  expect_warning(moments(D))
  expect_warning(mode(D))

  # Values
  expect_identical(mean(D), (a + b) / 2)
  expect_identical(median(D), (a + b) / 2)

})

test_that("Unif likelihood works", {

  # Preliminaries
  a <- 2
  b <- 5
  D <- Unif(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llunif(x, a, b)))

  # 2-Way Calls
  expect_identical(llunif(x, a, b), ll(D, x))
  expect_identical(ll(D)(x), ll(D, x))

})

test_that("Unif estim works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Unif(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(eunif(x, type = "mle")))
  expect_true(is.list(eunif(x, type = "me")))

  # 2-Way Calls
  expect_identical(eunif(x, type = "mle"), e(D, x, type = "mle"))
  expect_identical(eunif(x, type = "me"), e(D, x, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)

})

test_that("Unif small metrics work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Unif(a, b)
  set.seed(1)

  prm <- list(name = "min",
              val = seq(0.5, 2, by = 0.5))

  expect_no_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me"),
                       obs = c(20, 50),
                       sam = 1e2,
                       seed = 1)
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "SmallMetrics")

})
