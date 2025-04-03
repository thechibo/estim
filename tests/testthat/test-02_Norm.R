test_that("Norm distr works", {

  # Preliminaries
  mu <- 3
  sd <- 1
  D <- Norm(mu, sd)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Norm")

  # Errors
  expect_error(Norm(c(0, 1), 2))
  expect_error(Norm(0, -1))

})

test_that("Norm dpqr work", {

  # Preliminaries
  mu <- 3
  sd <- 1
  D <- Norm(mu, sd)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(p(D)(mu), 0.5)
  expect_identical(p(D)(Inf), 1)
  expect_identical(qn(D)(1), Inf)
  expect_identical(qn(D)(0.5), mu)

  # 2-Way Calls
  expect_identical(d(D)(1), dnorm(1, mu, sd))
  expect_identical(p(D)(1), pnorm(1, mu, sd))
  expect_equal(qn(D)(0.5), qnorm(0.5, mu, sd), tolerance = 0.01)
  expect_identical(d(D)(1), d(D, 1))
  expect_identical(p(D)(1), p(D, 1))
  expect_equal(qn(D)(0.5), qn(D, 0.5), tolerance = 0.01)

})

test_that("Norm moments work", {

  # Preliminaries
  mu <- 3
  sd <- 1
  D <- Norm(mu, sd)

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

test_that("Norm likelihood works", {

  # Preliminaries
  mu <- 3
  sd <- 1
  D <- Norm(mu, sd)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llnorm(x, mu, sd)))

  # 2-Way Calls
  expect_identical(llnorm(x, mu, sd), ll(D, x))
  expect_identical(ll(D)(x), ll(D, x))

})

test_that("Norm estim works", {

  # Preliminaries
  mu <- 3
  sd <- 1
  D <- Norm(mu, sd)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(enorm(x, type = "mle")))
  expect_true(is.list(enorm(x, type = "me")))

  # 2-Way Calls
  expect_identical(enorm(x, type = "mle"), e(D, x, type = "mle"))
  expect_identical(enorm(x, type = "me"), e(D, x, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

})

test_that("Norm avar works", {

  # Preliminaries
  mu <- 3
  sd <- 1
  D <- Norm(mu, sd)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(vnorm(mu, sd, type = "mle")))
  expect_true(is.numeric(vnorm(mu, sd, type = "me")))

  # 2-Way Calls
  expect_identical(vnorm(mu, sd, type = "mle"), avar(D, type = "mle"))
  expect_identical(vnorm(mu, sd, type = "me"), avar(D, type = "me"))
  expect_identical(vnorm(mu, sd, type = "mle"), avar_mle(D))
  expect_identical(vnorm(mu, sd, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.1)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.1)

})

test_that("Norm small metrics work", {

  # Preliminaries
  mu <- 3
  sd <- 1
  D <- Norm(mu, sd)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "mean",
              val = seq(0.5, 5, by = 0.5))

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

  prm <- list(name = "sd",
              val = seq(0.5, 5, by = 0.5))

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

test_that("Norm large metrics work", {

  # Preliminaries
  mu <- 3
  sd <- 1
  D <- Norm(mu, sd)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "mean",
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "LargeMetrics")

  prm <- list(name = "sd",
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "LargeMetrics")

})
