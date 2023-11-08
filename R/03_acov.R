# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Asymptotic Covariance                                                     ----
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
  Row <- Col <- Parameter <- Estimator <- NULL
  distr <- class(D0)[1]
  prm_name <- paste0(prm$name, prm$pos)
  y <- list()

  # Get the distributions
  Di <- lapply(seq_along(prm$val),
               FUN = function(i) { update_params(D0, prm, i) })

  # For each estimator
  for (j in est) {

    y[[j]] <- sapply(Di, simplify = "array", FUN = paste0("acov_", j))
    dimnames(y[[j]])[3] <- list(prm = prm$val)

  }

  # Create the acov data frame
  y <- simplify2array(y)
  z <- array_to_df(y)

  # Data Wrangling
  names(z) <- c("Row", "Col", "Parameter", "Estimator", "Value")
  z$Row <- factor(z$Row)
  z$Col <- factor(z$Col)
  z$Parameter <- as.numeric(as.character(z$Parameter))
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
                      colors = NULL,
                      save = FALSE,
                      path = getwd(),
                      name = "myplot.pdf",
                      width = 15,
                      height = 8) {

  # Colors
  if (is.null(colors)) {
    colors <- c("#0073C2", "#CD534C", "#EFC000", "#868686", "#003C67",
                "#7AA6DC", "#A73030", "#8F7700", "#3B3B3B", "#4A6990")
  }

  # Global variable binding
  Row <- Col <- Parameter <- Estimator <- Value <- NULL
  x$Cov <- paste(x$Row, "-", x$Col)

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
    ggh4x::facet_grid2(rows = ggplot2::vars(Row),
                       cols = ggplot2::vars(Col),
                       scales = "free",
                       independent = "y") +
    ggplot2::scale_color_manual(values = colors[seq_along(unique(x$Estimator))]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 25),
                   legend.key.size = ggplot2::unit(2, 'cm'),
                   plot.title = ggplot2::element_text(hjust = 0.5))

  plot(p)

  # Close the device
  if (save) {
    grDevices::dev.off()
  }

  # Return the plot
  invisible(p)

}
