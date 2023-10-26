# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Asymptotic Covariance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Calculate the Asymptotic Variance
#'
#' @description
#' This function calculates the asymptotic variance of an estimator as a
#' function of a single parameter, keeping the others constant. See Details.
#'
#' @param D0 A subclass of `Distribution`. The distribution family of interest.
#' @param prm A list containing three elements (name, pos, val). See Details.
#' @param est character. The estimator of interest. Can be a vector.
#'
#' @return A data.frame with columns named "Parameter", "Estimator", and
#' "Value".
#'
#' @export
#'
#' @inherit mle details examples
acov <- function(D0, prm, est = c("same", "me", "mle")) {

  # Preliminaries
  distr <- class(D0)[1]
  prm_name <- paste0(prm$name, prm$pos)

  # Create an array (prm x est)
  d <- list(prm = prm$val, est = est)
  y <- array(dim = lengths(d), dimnames = d)

  # Get the distribution
  Di <- lapply(seq_along(prm$val),
               FUN = function(i) { update_params(D0, prm, i) })

  # For each estimator
  for (j in est) {

    y[ , j] <- vapply(Di,
                      FUN = function(x) {
                        do.call(paste0("acov_", j), list(distr = x))[1, 1]
                      })

  }

  # Create the acov data frame
  z <- array_to_df(y)

  # Data Wrangling
  names(z) <- c("Parameter", "Estimator", "Value")
  z$Estimator <- factor(z$Estimator)

  # Return the object
  z

}


#' @title Plot the Asymptotic Covariance
#'
#' @description
#' This function provides an easy way to illustrate the output of `acov()`,
#' using the `ggplot2` package. A line chart is created in which each estimator
#' is plotted with a different color and line type. The plot can be saved in
#' pdf format.
#'
#' @param x A data.frame. The result of `acov()`.
#' @param colors character. The colors to be used in the plot.
#' @param save logical. Should the plot be saved?
#' @param path A path to the directory in which the plot will be saved.
#' @param name character. The name of the output pdf file.
#' @param width numeric. The plot width in inches.
#' @param height numeric. The plot height in inches.
#'
#' @return The plot is returned invisibly in the form of a `ggplot` object.
#'
#' @import ggplot2
#' @export
#'
#' @seealso [acov()]
#' @inherit mle examples
plot_acov <- function(x,
                      colors,
                      save = FALSE,
                      path = getwd(),
                      name = "myplot.pdf",
                      width = 15,
                      height = 8) {

  # Global variable binding
  Parameter <- Value <- Estimator <- NULL

  # Save the plot
  if (save) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(file.path(path, name), width = width, height = height)
  }

  # Create the plot
  p <- ggplot2::ggplot(x) +
    ggplot2::geom_line(ggplot2::aes(x = Parameter,
                                    y = Value,
                                    col = Estimator,
                                    linetype = Estimator),
                       linewidth = 1.5) +
    ggplot2::labs(title = "Estimator Asymptotic Variance Comparison",
                  y = "Value",
                  x = "Parameter Value") +
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
