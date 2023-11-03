test_consistency <- function(est, D0, n = 1e4, seed = 1) {

  # Random Sampling
  set.seed(seed)
  sam <- distr::r(D0)(n)

  # Estimation
  d <- max(abs(get_params(D0) - do.call(est, list(x = sam, distr = D0))))

  # Report
  cat("Maximum Absolute Distance: ", d, "\n")
  invisible(d)

}

test_acov <- function(est, D0, n = 1e4, m = 1e3, seed = 1) {

  # Preliminaries
  set.seed(seed)
  params <- get_params(D0)
  y <- matrix(nrow = m, ncol = length(params))

  # Loading bar
  pb <- loading_bar(total = m)

  for (i in 1:m) {

    # Progress Bar
    pb$tick()

    sam <- distr::r(D0)(n)
    y[i, ] <- sqrt(n) * (params - do.call(est, list(x = sam, distr = D0)))

  }

  d <- list(acov_sample = var(y),
            acov_theory = do.call(paste0("acov_", est), list(distr = D0)))
  d$mad <- max(abs(d$acov_sample - d$acov_theory))

  # Report
  cat("Maximum Absolute Distance: ", d$mad, "\n")
  invisible(d)

}
