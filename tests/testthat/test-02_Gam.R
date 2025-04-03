test_that("Gam distr works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Gam(a, b)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Gam")

  # Errors
  expect_error(Gam(c(0.1, 0.2, 0.3)))
  expect_error(Gam(-1, 2))
  expect_error(Gam(1, -2))

})

test_that("Gam dpqr work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Gam(a, b)
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
  expect_identical(p(D)(Inf), 1)
  expect_identical(p(D)(0), 0)
  expect_identical(qn(D)(1), Inf)
  expect_identical(qn(D)(0), 0)
  expect_identical(sum(x >= 0), n)

  # 2-Way Calls
  expect_identical(d(D)(0.4), dgamma(0.4, shape = a, scale = b))
  expect_identical(p(D)(0.4), pgamma(0.4, shape = a, scale = b))
  expect_equal(qn(D)(0.4), qgamma(0.4, shape = a, scale = b), tolerance = 1e-8)
  expect_identical(d(D)(0.4), d(D, 0.4))
  expect_identical(p(D)(0.4), p(D, 0.4))
  expect_identical(qn(D)(0.4), qn(D, 0.4))

})

test_that("Gam moments work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Gam(a, b)

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

})

test_that("Gam likelihood works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Gam(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llgamma(x, shape = a, scale = b)))

  # 2-Way Calls
  expect_identical(llgamma(x, shape = a, scale = b), ll(D, x))
  expect_identical(ll(D)(x), ll(D, x))

  # ll and lloptim convergence to a0 comparison
  method <- "L-BFGS-B"
  lower <- 1e-5
  upper <- Inf
  tx <- c(log(mean(x)), mean(log(x)))

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
                fn = function(par, x, distr) { ll(Gam(par[1], par[2]), x) },
                x = x,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  expect_equal(c(shape = par1, scale = mean(x) / par1), par2, tolerance = 0.01)

})

test_that("Gam estim works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Gam(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(egamma(x, type = "mle")))
  expect_true(is.list(egamma(x, type = "me")))
  expect_true(is.list(egamma(x, type = "same")))

  # 2-Way Calls
  expect_identical(egamma(x, type = "mle"), e(D, x, type = "mle"))
  expect_identical(egamma(x, type = "me"), e(D, x, type = "me"))
  expect_identical(egamma(x, type = "same"), e(D, x, type = "same"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.05)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.05)
  d <- test_consistency("same", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.05)

})

test_that("Gam avar works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Gam(a, b)

  # Types
  expect_true(is.numeric(vgamma(a, b, type = "mle")))
  expect_true(is.numeric(vgamma(a, b, type = "me")))
  expect_true(is.numeric(vgamma(a, b, type = "same")))

  # 2-Way Calls
  expect_identical(vgamma(a, b, type = "mle"), avar(D, type = "mle"))
  expect_identical(vgamma(a, b, type = "me"), avar(D, type = "me"))
  expect_identical(vgamma(a, b, type = "same"), avar(D, type = "same"))
  expect_identical(vgamma(a, b, type = "mle"), avar_mle(D))
  expect_identical(vgamma(a, b, type = "me"), avar_me(D))
  expect_identical(vgamma(a, b, type = "same"), avar_same(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.07)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("same", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.07)

})

test_that("Gam small metrics work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Gam(a, b)
  set.seed(1)

  prm <- list(name = "shape",
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

test_that("Gam large metrics work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Gam(a, b)
  set.seed(1)

  prm <- list(name = "shape",
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

  prm <- list(name = "scale",
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
