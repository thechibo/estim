test_that("Multinom distr works", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Multinom")

  # Errors
  expect_error(Multinom(-10, 0.5))
  expect_error(Multinom(0, 0.5))
  expect_error(Multinom(10, 5))
  expect_error(Multinom(3:4, 0.5))

})

test_that("Multinom dpqr work", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_equal(d(D)(c(N, 0, 0)), p[1] ^ N, tolerance = 0.01)
  expect_identical(sum(x %in% 0:N), length(p)*n)
  expect_identical(sum(colSums(x) == N), n)

  # 2-Way Calls
  expect_identical(d(Multinom(N, p))(c(4, 3, 3)),
                   dmultinom(c(4, 3, 3), N, p))

  # Special Case: Binomial
  expect_equal(d(Multinom(N, c(0.3, 0.7)))(c(2, N-2)),
                   dbinom(2, N, 0.3), tolerance = 0.01)

})

test_that("Multinom moments work", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)

  # Types
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(entro(D)))
  expect_true(is.numeric(finf(D)))

  # Values
  expect_identical(mean(D), N * p)
  expect_equal(var(D)[1, 1], N * p[1] * (1 - p[1]), tolerance = 0.01)

})

test_that("Multinom likelihood works", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llmultinom(x, size = N, prob = p)))

  # 2-Way Calls
  expect_identical(llmultinom(x, N, p), ll(x, c(N, p), D))

})

test_that("Multinom estim works", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(emultinom(x, type = "mle")))
  expect_true(is.numeric(emultinom(x, type = "me")))

  # 2-Way Calls
  expect_identical(emultinom(x, type = "mle"), estim(x, D, type = "mle"))
  expect_identical(emultinom(x, type = "me"), estim(x, D, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true[-1], d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true[-1], d$prm_est, tolerance = 0.01)

})

test_that("Multinom avar works", {

  # Preliminaries
  N <- 1e3
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)
  k <- length(p)

  # Types
  expect_true(is.numeric(vmultinom(N, p, type = "mle")))
  expect_true(is.numeric(vmultinom(N, p, type = "me")))

  # 2-Way Calls
  expect_identical(vmultinom(N, p, type = "mle"), avar(D, type = "mle"))
  expect_identical(vmultinom(N, p, type = "me"), avar(D, type = "me"))
  expect_identical(vmultinom(N, p, type = "mle"), avar_mle(D))
  expect_identical(vmultinom(N, p, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est[1:(k-1), 1:(k-1)], tolerance = 0.1)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est[1:(k-1), 1:(k-1)], tolerance = 0.1)

})

test_that("Multinom small metrics work", {

  # Preliminaries
  N <- 100
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)

  prm <- list(name = "prob",
              pos = NULL,
              val = seq(0.5, 0.8, by = 1))

  expect_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me"),
                       obs = c(20, 50),
                       sam = 1e2,
                       seed = 1)
  )

})

test_that("Multinom large metrics work", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)

  prm <- list(name = "prob",
              pos = NULL,
              val = seq(0.5, 0.8, by = 0.1))

  expect_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
  )

})
