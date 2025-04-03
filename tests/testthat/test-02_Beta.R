test_that("Beta distr works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Beta(a, b)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Beta")

  # Errors
  expect_error(Beta(c(0.1, 0.2, 0.3)))
  expect_error(Beta(-1, 2))
  expect_error(Beta(1, -2))

})

test_that("Beta dpqr work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Beta(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(d(D)(0), 0)
  expect_identical(d(D)(1), 0)
  expect_identical(p(D)(1), 1)
  expect_identical(p(D)(0), 0)
  expect_identical(qn(D)(1), 1)
  expect_identical(qn(D)(0), 0)
  expect_identical(sum(x <= 1), n)
  expect_identical(sum(x >= 0), n)

  # 2-Way Calls
  expect_identical(d(D)(0.4), dbeta(0.4, a, b))
  expect_identical(p(D)(0.4), pbeta(0.4, a, b))
  expect_equal(qn(D)(0.4), qbeta(0.4, a, b), tolerance = 1e-8)
  expect_identical(d(D)(0.4), d(D, 0.4))
  expect_identical(p(D)(0.4), p(D, 0.4))
  expect_equal(qn(D)(0.4), qn(D, 0.4), tolerance = 1e-8)

})

test_that("Beta moments work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Beta(a, b)

  # Types
  expect_true(is.list(moments(D)))
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(median(D)))
  expect_true(is.numeric(mode(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(sd(D)))
  expect_true(is.numeric(skew(D)))
  expect_true(is.numeric(kurt(D)))
  expect_true(is.numeric(entro(D)))
  expect_true(is.numeric(finf(D)))

  # Values
  expect_identical(mean(D), a / (a + b))

})

test_that("Beta likelihood works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Beta(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llbeta(x, a, b)))

  # 2-Way Calls
  expect_identical(llbeta(x, a, b), ll(D, x))
  expect_identical(ll(D)(x), ll(D, x))

  # ll and lloptim convergence to a0 comparison
  method <- "L-BFGS-B"
  lower <- 1e-5
  upper <- Inf
  tx  <- c(mean(log(x)), mean(log(1 - x)))

  par1 <- optim(par = sum(unlist(same(D, x))),
                fn = lloptim,
                gr = dlloptim,
                tx = tx,
                distr = D,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  par2 <- optim(par = unlist(same(D, x)),
                fn = function(par, x) { ll(Beta(par[1], par[2]), x) },
                x = x,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  expect_equal(par1, sum(par2), tolerance = 0.01)

})

test_that("Beta estim works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Beta(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(ebeta(x, type = "mle")))
  expect_true(is.list(ebeta(x, type = "me")))
  expect_true(is.list(ebeta(x, type = "same")))

  # 2-Way Calls
  expect_identical(ebeta(x, type = "mle"), e(D, x, type = "mle"))
  expect_identical(ebeta(x, type = "me"), e(D, x, type = "me"))
  expect_identical(ebeta(x, type = "same"), e(D, x, type = "same"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("same", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

})

test_that("Beta avar works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Beta(a, b)

  # Types
  expect_true(is.numeric(vbeta(a, b, type = "mle")))
  expect_true(is.numeric(vbeta(a, b, type = "me")))
  expect_true(is.numeric(vbeta(a, b, type = "same")))

  # 2-Way Calls
  expect_identical(vbeta(a, b, type = "mle"), avar(D, type = "mle"))
  expect_identical(vbeta(a, b, type = "me"), avar(D, type = "me"))
  expect_identical(vbeta(a, b, type = "same"), avar(D, type = "same"))
  expect_identical(vbeta(a, b, type = "mle"), avar_mle(D))
  expect_identical(vbeta(a, b, type = "me"), avar_me(D))
  expect_identical(vbeta(a, b, type = "same"), avar_same(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("same", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)

})

test_that("Beta small metrics work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Beta(a, b)
  set.seed(1)

  prm <- list(name = "shape1",
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me", "same"),
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

test_that("Beta large metrics work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Beta(a, b)
  set.seed(1)

  prm <- list(name = "shape1",
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me", "same"))
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "LargeMetrics")

})
