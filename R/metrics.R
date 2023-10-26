# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Metrics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Calculate Metrics
#'
#' @description
#' This function performs Monte Carlo simulations to estimate the main metrics
#' (bias, variance, and RMSE) characterizing the finite sample behavior of an
#' estimator. The function evaluates the metrics as a function of a single
#' parameter, keeping the other ones constant. See Details.
#'
#' @param D0 A subclass of `Distribution`. The distribution family of interest.
#' @param prm A list containing three elements (name, pos, val). See Details.
#' @param obs numeric. The size of each sample. Can be a vector.
#' @param est character. The estimator of interest. Can be a vector.
#' @param sam numeric. The number of Monte Carlo samples used to estimate the
#' metrics.
#' @param seed numeric. Passed to set.seed() for reproducibility.
#'
#' @details
#' The distribution `D0` is used to specify an initial distribution. The list
#' `prm` contains details concerning a single parameter that is allowed to
#' change values. The quantity of interest is evaluated as a function of this
#' parameter.
#'
#' Specifically, `prm` includes three elements named "name", "pos", and "val".
#' The first two elements determine the exact parameter that changes, while the
#' third one is a numeric vector holding the values it takes. For example,
#' in the case of the Multivariate Gamma distribution,
#' `D0 <- MGamma(shape = c(1, 2), scale = 3)` and
#' `prm <- list(name = "shape", pos = 2, val = seq(1, 1.5, by = 0.1))`
#' means that the evaluation will be performed for the MGamma distributions with
#' shape parameters `(1, 1)`, `(1, 1.1)`, ..., `(1, 1.5)` and scale `3`. Notice
#' that the initial shape parameter `2` in `D0` is not utilized in the function.
#'
#' @return A data.frame with columns named "Parameter", "Observations",
#' "Estimator", "Metric", and "Value".
#'
#' @export
#'
#' @seealso [plot_metrics()]
#' @inherit mle examples
metrics <- function(D0, prm,
                    obs = c(20, 50, 100),
                    est = c("same", "me", "mle"),
                    sam = 1e4,
                    seed = 1) {

  # Preliminaries
  set.seed(seed)
  distr <- class(D0)[1]
  nmax <- max(obs)
  prm_name <- paste0(prm$name, prm$pos)

  # Create an array (prm x obs x est x sam)
  d <- list(prm = prm$val, obs = obs, est = est, sam = 1:sam)
  y <- array(dim = lengths(d), dimnames = d)

  # Loading bar
  pb <- loading_bar(total = length(prm$val) * length(obs) * length(est))

  # For each value of prm
  for (i in seq_along(prm$val)) {

    rDi <- distr::r(update_params(D0, prm, i))
    x <- replicate(sam, { rDi(nmax) })

    # For each sample size
    for(j in seq_along(obs)) {

      # For each estimator
      for (k in est) {

        # Progress Bar
        pb$tick()

        # Estimate
        y[i, j, k, ] <- apply(x[ , 1:obs[j], ],
                              MARGIN = 3,
                              FUN = k,
                              distr = distr)[prm_name, ] - prm$val[i]

      }
    }
  }

  # Calculate the metrics
  bias <- apply(y, MARGIN = c("prm", "obs", "est"), FUN = mean)
  var <- apply(y, MARGIN = c("prm", "obs", "est"), FUN = s2)
  rmse <- sqrt(bias ^ 2 + var)

  # Create the metrics data frame
  d <- append(dimnames(bias), list(metric = c("Bias", "Variance", "RMSE")))
  z <- array(c(bias, var, rmse), dim = lengths(d), dimnames = d)
  z <- array_to_df(z)

  # Data Wrangling
  names(z) <- c("Parameter", "Observations", "Estimator", "Metric", "Value")
  z$Estimator <- factor(z$Estimator)
  z$Metric <- factor(z$Metric)
  z$Observations <- factor(z$Observations)

  # Return the object
  z

}

#' @title Plot Metrics
#'
#' @description
#' This function provides an easy way to illustrate the output of `metrics()`,
#' using the `ggplot2` package. A grid of line charts is created for each metric
#' and sample size. Each estimator is plotted with a different color and
#' linetype. The plot can be saved in pdf format.
#'
#' @param x A data.frame. The result of `metrics()`.
#' @param colors character. The colors to be used in the plot.
#' @param save logical. Should the plot be saved?
#' @param path A path to the directory in which the plot will be saved.
#' @param name character. The name of the output pdf file.
#' @param width numeric. The plot width in inches.
#' @param height numeric. The plot height in inches.
#'
#' @return The plot is returned invisibly in the form of a `ggplot` object.
#'
#' @import ggplot2 ggplot geom_line aes labs facet_grid vars scale_color_manual
#' @import ggplot2 theme_minimal theme element_text unit
#' @export
#'
#' @seealso [metrics()]
#' @inherit mle examples
plot_metrics <- function(x,
                         colors,
                         save = FALSE,
                         path = getwd(),
                         name = "myplot.pdf",
                         width = 15,
                         height = 8) {

  # Global variable binding
  Parameter <- Value <- Estimator <- Observations <- Metric <- NULL

  # Save the plot
  if (save) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    pdf(file.path(path, name), width = width, height = height)
  }

  # Create the plot
  p <- ggplot2::ggplot(x) +
    ggplot2::geom_line(ggplot2::aes(x = Parameter,
                                    y = Value,
                                    col = Estimator,
                                    linetype = Estimator),
                       linewidth = 1.5) +
    ggplot2::labs(title = "Estimator Metrics Comparison",
                  y = "Value",
                  x = "Parameter Value") +
    ggplot2::facet_grid(rows = ggplot2::vars(Observations),
                        cols = ggplot2::vars(Metric),
                        switch = "y") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 25),
                   legend.key.size = ggplot2::unit(2, 'cm'),
                   plot.title = ggplot2::element_text(hjust = 0.5))

  plot(p)

  # Close the device
  if (save) grDevices::dev.off()

  # Return the plot
  invisible(p)

}
