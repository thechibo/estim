test_that("Chisq distr works", {

  # Preliminaries
  df <- 3
  D <- Chisq(df)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Chisq")

  # Errors
  expect_error(Chisq(1:2))
  expect_error(Chisq(-1))

})

test_that("Chisq dpqr work", {

  # Preliminaries
  df <- 3
  D <- Chisq(df)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(d(D)(-1), 0)
  expect_identical(p(D)(0), 0)
  expect_identical(p(D)(Inf), 1)
  expect_identical(qn(D)(1), Inf)
  expect_identical(qn(D)(0), 0)
  expect_identical(sum(x > 0), n)

  # 2-Way Calls
  expect_identical(d(D)(1), dchisq(1, df))
  expect_identical(p(D)(1), pchisq(1, df))
  expect_equal(qn(D)(0.5), qchisq(0.5, df), tolerance = 0.01)

})

test_that("Chisq moments work", {

  # Preliminaries
  df <- 3
  D <- Chisq(df)

  # Types
  expect_true(is.numeric(mean(D)))
  expect_warning(is.numeric(median(D)))
  expect_true(is.numeric(mode(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(sd(D)))
  expect_true(is.numeric(skew(D)))
  expect_true(is.numeric(kurt(D)))
  expect_true(is.numeric(entro(D)))
  expect_true(is.numeric(finf(D)))

})

test_that("Chisq likelihood works", {

  # Preliminaries
  df <- 3
  D <- Chisq(df)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llchisq(x, df)))

  # 2-Way Calls
  expect_identical(llchisq(x, df), ll(x, df, D))

})

test_that("Chisq estim works", {

  # Preliminaries
  df <- 3
  D <- Chisq(df)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(echisq(x, type = "mle")))
  expect_true(is.numeric(echisq(x, type = "me")))

  # 2-Way Calls
  expect_identical(echisq(x, type = "mle"), estim(x, D, type = "mle"))
  expect_identical(echisq(x, type = "me"), estim(x, D, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

})

test_that("Chisq avar works", {

  # Preliminaries
  df <- 3
  D <- Chisq(df)

  # Types
  expect_true(is.numeric(vchisq(df, type = "mle")))
  expect_true(is.numeric(vchisq(df, type = "me")))

  # 2-Way Calls
  expect_identical(vchisq(df, type = "mle"), avar(D, type = "mle"))
  expect_identical(vchisq(df, type = "me"), avar(D, type = "me"))
  expect_identical(vchisq(df, type = "mle"), avar_mle(D))
  expect_identical(vchisq(df, type = "me"), avar_me(D))

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)

})

test_that("Chisq small metrics work", {

  # Preliminaries
  df <- 3
  D <- Chisq(df)
  set.seed(1)

  prm <- list(name = "df",
              pos = NULL,
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me"),
                       obs = c(20, 50),
                       sam = 1e2,
                       seed = 1)
  )

  expect_no_error(
    plot_small_metrics(x,
                       save = TRUE,
                       path = tempdir())
  )

  # Types
  expect_s3_class(x, "data.frame")
  expect_true(is.numeric(x$Parameter))
  expect_s3_class(x$Observations, "factor")
  expect_s3_class(x$Estimator, "factor")
  expect_s3_class(x$Metric, "factor")
  expect_true(is.numeric(x$Value))

})

test_that("Chisq large metrics work", {

  # Preliminaries
  df <- 3
  D <- Chisq(df)

  prm <- list(name = "df",
              pos = NULL,
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
  )

  expect_no_error(
    plot_large_metrics(x,
                       save = TRUE,
                       path = tempdir())
  )

  # Types
  expect_s3_class(x, "data.frame")
  expect_true(is.numeric(x$Parameter))
  expect_s3_class(x$Estimator, "factor")
  expect_true(is.numeric(x$Value))

})
