test_that("Lnorm distr works", {

  # Preliminaries
  meanlog <- 3
  sdlog <- 1
  D <- Lnorm(meanlog, sdlog)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Lnorm")

  # Errors
  expect_error(Lnorm(c(0, 1), 2))
  expect_error(Lnorm(0, -1))
  expect_error(Lnorm(0, c(1, 1)))

})

test_that("Lnorm dpqr work", {

  # Preliminaries
  meanlog <- 3
  sdlog <- 1
  D <- Lnorm(meanlog, sdlog)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(p(D)(0), 0)
  expect_identical(p(D)(Inf), 1)
  expect_identical(qn(D)(1), Inf)
  expect_identical(qn(D)(0), 0)

  # 2-Way Calls
  expect_identical(d(D)(1), dlnorm(1, meanlog, sdlog))
  expect_identical(p(D)(1), plnorm(1, meanlog, sdlog))
  expect_equal(qn(D)(0.5), qlnorm(0.5, meanlog, sdlog), tolerance = 0.01)
  expect_identical(d(D)(1), d(D, 1))
  expect_identical(p(D)(1), p(D, 1))
  expect_equal(qn(D)(0.5), qn(D, 0.5), tolerance = 0.01)

})

test_that("Lnorm moments work", {

  # Preliminaries
  meanlog <- 3
  sdlog <- 1
  D <- Lnorm(meanlog, sdlog)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

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

test_that("Lnorm likelihood works", {

  # Preliminaries
  meanlog <- 3
  sdlog <- 1
  D <- Lnorm(meanlog, sdlog)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(lllnorm(x, meanlog, sdlog)))

  # 2-Way Calls
  expect_identical(lllnorm(x, meanlog, sdlog), ll(D, x))
  expect_identical(ll(D)(x), ll(D, x))

})

test_that("Lnorm estim works", {

  # Preliminaries
  meanlog <- 3
  sdlog <- 1
  D <- Lnorm(meanlog, sdlog)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(eexp(x, type = "mle")))
  expect_true(is.list(eexp(x, type = "me")))

  # 2-Way Calls
  expect_identical(elnorm(x, type = "mle"), e(D, x, type = "mle"))
  expect_identical(elnorm(x, type = "me"), e(D, x, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

})

test_that("Lnorm avar works", {

  # Preliminaries
  meanlog <- 3
  sdlog <- 1
  D <- Lnorm(meanlog, sdlog)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(vlnorm(meanlog, sdlog, type = "mle")))
  expect_true(is.numeric(vlnorm(meanlog, sdlog, type = "me")))

  # 2-Way Calls
  expect_identical(vlnorm(meanlog, sdlog, type = "mle"), avar(D, type = "mle"))
  expect_identical(vlnorm(meanlog, sdlog, type = "me"), avar(D, type = "me"))
  expect_identical(vlnorm(meanlog, sdlog, type = "mle"), avar_mle(D))
  expect_identical(vlnorm(meanlog, sdlog, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.1)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.1)

})

test_that("Lnorm small metrics work", {

  # Preliminaries
  meanlog <- 3
  sdlog <- 1
  D <- Lnorm(meanlog, sdlog)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "meanlog",
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

  prm <- list(name = "sdlog",
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

test_that("Lnorm large metrics work", {

  # Preliminaries
  meanlog <- 3
  sdlog <- 1
  D <- Lnorm(meanlog, sdlog)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "meanlog",
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

  prm <- list(name = "sdlog",
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
