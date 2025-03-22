test_that("Fisher distr works", {

  # Preliminaries
  df1 <- 2
  df2 <- 3
  D <- Fisher(df1, df2)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Fisher")

  # Errors
  expect_error(Fisher(1:2, 1))
  expect_error(Fisher(-1, 2))
  expect_error(Fisher(1, -2))

})

test_that("Fisher dpqr work", {

  # Preliminaries
  df1 <- 2
  df2 <- 3
  D <- Fisher(df1, df2)
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
  expect_equal(d(D)(1), df(1, df1, df2), tolerance = 0.01)
  expect_equal(p(D)(1), pf(1, df1, df2), tolerance = 0.01)
  expect_equal(qn(D)(0.5), qf(0.5, df1, df2), tolerance = 0.01)

})

test_that("Fisher moments work", {

  # Preliminaries
  df1 <- 12
  df2 <- 15
  D <- Fisher(df1, df2)

  # Types
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(median(D)))
  expect_true(is.numeric(mode(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(sd(D)))
  expect_true(is.numeric(skew(D)))
  expect_true(is.numeric(kurt(D)))
  expect_true(is.numeric(entro(D)))

  # Errors
  expect_error(mode(Fisher(2, 3)))
  expect_error(var(Fisher(4, 3)))
  expect_error(sd(Fisher(4, 3)))
  expect_error(skew(Fisher(6, 3)))
  expect_error(kurt(Fisher(8, 3)))

})

test_that("Fisher likelihood works", {

  # Preliminaries
  set.seed(1)
  df1 <- 3
  df2 <- 2
  D <- Fisher(df1, df2)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llf(x, df1, df2)))

  # 2-Way Calls
  expect_identical(llf(x, df1, df2), ll(x, c(df1, df2), D))

})
