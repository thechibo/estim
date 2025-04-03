test_that("Stud distr works", {

  # Preliminaries
  df <- 2
  D <- Stud(df)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Stud")

  # Errors
  expect_error(Stud(1:2))
  expect_error(Stud(-1))

})

test_that("Stud dpqr work", {

  # Preliminaries
  df <- 2
  D <- Stud(df)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(p(D)(Inf), 1)
  expect_identical(qn(D)(1), Inf)

  # 2-Way Calls
  expect_equal(d(D)(1), dt(1, df), tolerance = 0.01)
  expect_equal(p(D)(1), pt(1, df), tolerance = 0.01)
  expect_equal(qn(D)(0.5), qt(0.5, df), tolerance = 0.01)
  expect_equal(d(D)(1), d(D, 1), tolerance = 0.01)
  expect_equal(p(D)(1), p(D, 1), tolerance = 0.01)
  expect_equal(qn(D)(0.5), qn(D, 0.5), tolerance = 0.01)

})

test_that("Stud moments work", {

  # Preliminaries
  df <- 12
  D <- Stud(df)

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
  expect_error(var(Stud(2)))
  expect_error(sd(Stud(2)))
  expect_error(skew(Stud(2)))
  expect_error(kurt(Stud(2)))

  # Values
  D1 <- Stud(1e4)
  D2 <- Norm()
  expect_equal(var(D1), var(D2), tolerance = 0.01)
  expect_equal(skew(D1), skew(D2), tolerance = 0.01)
  expect_equal(kurt(D1), kurt(D2), tolerance = 0.01)
  expect_equal(entro(D1), entro(D2), tolerance = 0.01)

})

test_that("Stud likelihood works", {

  # Preliminaries
  set.seed(1)
  D <- Stud(3)
  n <- 100L
  df <- D@df
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llt(x, df)))

  # 2-Way Calls
  expect_identical(llt(x, df), ll(D, x))
  expect_identical(ll(D)(x), ll(D, x))

})
