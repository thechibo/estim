test_consistency <- function(est, D0, n = 1e4, seed = 1, ...) {

  # Random Sampling
  set.seed(seed)
  sam <- r(D0)(n)

  # Return
  list(prm_true = get_unknown_params(D0),
       prm_est = do.call(est, list(distr = D0, x = sam, ...)))

}

test_avar <- function(est, D0, n = 1e4, m = 1e3, seed = 1, bar = FALSE, ...) {

  # Preliminaries
  set.seed(seed)
  params <- get_unknown_params(D0, list = FALSE)
  y <- matrix(nrow = m, ncol = length(params))

  # Loading bar
  if (bar) {
    pb <- loading_bar(total = m)
  }

  # Estimation
  for (i in 1:m) {

    # Progress Bar
    if (bar) {
      pb$tick()
    }

    # CLT
    sam <- r(D0)(n)
    params_est <- unlist(do.call(est, list(x = sam, distr = D0, ...)))
    y[i, ] <- sqrt(n) * (params - params_est)

  }

  # Calculate avar
  avar_est <- var(y)
  prm_names <- names(get_unknown_params(D0, list = FALSE))

  if (nrow(avar_est) == 1) {
    avar_est <- as.vector(avar_est)
    names(avar_est) <- prm_names
  } else {
    dimnames(avar_est) <- list(prm_names, prm_names)
  }

  # Return
  list(avar_true = do.call(paste0("avar_", est), list(distr = D0)),
       avar_est = avar_est)

}
