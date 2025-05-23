# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Metrics                                                                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculation            ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("SmallMetrics", slots = list(D = "Distribution",
                                      est = "character",
                                      df = "data.frame"))

setValidity("SmallMetrics", function(object) {
  if(!all(names(object@df) %in% c("Parameter", "Observations", "Estimator",
                                "Metric", "Value"))) {
    stop("df must have columns 'Parameter, 'Observations', 'Estimator',
                           'Metric', 'Value'")
  }
  if(!("data.frame" %in% class(object@df))) {
    stop("df has to be a data.frame")
  }
  if(!is.numeric(object@df$Parameter)) {
    stop("Column 'Parameter' has to be numeric")
  }
  if(!is.factor(object@df$Observations)) {
    stop("Column 'Observations' has to be a factor")
  }
  if(!is.factor(object@df$Estimator)) {
    stop("Column 'Estimator' has to be factor")
  }
  if(!is.factor(object@df$Metric)) {
    stop("Column 'Metric' has to be factor")
  }
  if(!is.numeric(object@df$Value)) {
    stop("Column 'Value' has to be numeric")
  }
  TRUE
})

#' @title Small Sample Metrics
#' @name SmallMetrics
#' @aliases small_metrics
#'
#' @description
#' This function performs Monte Carlo simulations to estimate the main metrics
#' (bias, variance, and RMSE) characterizing the small (finite) sample behavior
#' of an estimator. The function evaluates the metrics as a function of a single
#' parameter, keeping the other ones constant. See Details.
#'
#' @param D A subclass of `Distribution`. The distribution family of interest.
#' @param prm A list containing three elements (name, pos, val). See Details.
#' @param obs numeric. The size of each sample. Can be a vector.
#' @param est character. The estimator of interest. Can be a vector.
#' @param sam numeric. The number of Monte Carlo samples used to estimate the
#' metrics.
#' @param seed numeric. Passed to `set.seed()` for reproducibility.
#' @param df data.frame. a data.frame with columns named "Row", "Col",
#' "Parameter", "Estimator", and "Value".
#' @param ... extra arguments.
#'
#' @details
#' The distribution `D` is used to specify an initial distribution. The list
#' `prm` contains details concerning a single parameter that is allowed to
#' change values. The quantity of interest is evaluated as a function of this
#' parameter.
#'
#' The `prm` list includes two elements named "name" and "val". The first one
#' specifies the parameter that changes, and the second one is a numeric vector
#' holding the values it takes.
#'
#' In case the parameter of interest is a vector, a third element named "pos"
#' can be specified to indicate the exact paramater that changes. In the example
#' shown below, the evaluation will be performed for the Dirichlet distributions
#' with shape parameters `(0.5, 1)`, `(0.6, 1)`, ..., `(2, 1)`. Notice that the
#' initial shape parameter value (`1`) is not utilized in the function.
#'
#' @return An object of class `SmallMetrics` with slots `D`, `est`, and `df`.
#'
#' @export
#'
#' @seealso [LargeMetrics], [PlotMetrics]
#' @examples \donttest{
#' # -----------------------------------------------------
#' # Beta Distribution Example
#' # -----------------------------------------------------
#'
#' D <- Beta(shape1 = 1, shape2 = 2)
#'
#' prm <- list(name = "shape1",
#'             val = seq(0.5, 2, by = 0.1))
#'
#' x <- small_metrics(D, prm,
#'                    est = c("mle", "me", "same"),
#'                    obs = c(20, 50),
#'                    sam = 1e2,
#'                    seed = 1)
#'
#' plot(x)
#'
#' # -----------------------------------------------------
#' # Dirichlet Distribution Example
#' # -----------------------------------------------------
#'
#' D <- Dir(alpha = 1:2)
#'
#' prm <- list(name = "alpha",
#'             pos = 1,
#'             val = seq(0.5, 2, by = 0.1))
#'
#' x <- small_metrics(D, prm,
#'                    est = c("mle", "me", "same"),
#'                    obs = c(20, 50),
#'                    sam = 1e2,
#'                    seed = 1)
#'
#' plot(x)
#' }
SmallMetrics <- function(D, est, df) {
  new("SmallMetrics", D = D, est = est, df = df)
}

#' @rdname SmallMetrics
#' @export
small_metrics <- function(D,
                          prm,
                          est = c("same", "me", "mle"),
                          obs = c(20, 50, 100),
                          sam = 1e4,
                          seed = 1,
                          ...) {

  if (class(D) %in% c("Cat", "Multinom")) {
    stop("This function is not implemented for the Categorical and Multinomial
    distributions, since changing the probability value of one category forces
    the rest to change as well.")
  }

  # Preliminaries
  set.seed(seed)
  distr <- class(D)[1]
  nmax <- max(obs)
  prm_name <- paste0(prm$name, prm$pos)

  # Univariate of Multivariate
  if (distr %in% c("Dir")) {
    setk <- set1of3
    mar <- 3
  } else if (distr %in% c("Multigam")) {
    setk <- set1of3
    mar <- 3
  } else {
    setk <- set1of2
    mar <- 2
  }

  # Unidimensional of Multidimensional
  if (distr %in% c("Bern", "Binom", "Exp", "Pois", "Chisq", "Geom", "Nbinom")) {
    setd <- set1of1
  } else {
    setd <- set1of2
  }

  # Create an array (prm x obs x est x sam)
  d <- list(prm = prm$val, obs = obs, est = est, sam = 1:sam)
  y <- array(dim = lengths(d), dimnames = d)

  # Loading bar
  pb <- loading_bar(total = length(prm$val) * length(obs) * length(est))

  # For each value of prm
  for (i in seq_along(prm$val)) {

    rDi <- r(update_params(D, prm, i))
    x <- replicate(sam, { rDi(nmax) })

    # For each sample size
    for(j in seq_along(obs)) {

      # For each estimator
      for (k in est) {

        # Progress Bar
        pb$tick()

        # Estimate

        list_estim <- apply(setk(x, 1:obs[j]),
                            MARGIN = mar,
                            FUN = k,
                            distr = D,
                            ...)

        mat_estim <- do.call(cbind, lapply(list_estim, unlist))

        y[i, j, k, ] <- setd(mat_estim, prm_name) - prm$val[i]

      }
    }
  }

  # Calculate the metrics
  bias <- apply(y, MARGIN = c("prm", "obs", "est"), FUN = mean)
  var <- apply(y, MARGIN = c("prm", "obs", "est"), FUN = bsd)
  rmse <- sqrt(bias ^ 2 + var)

  # Create the metrics data frame
  d <- append(dimnames(bias), list(metric = c("Bias", "Variance", "RMSE")))
  z <- array(c(bias, var, rmse), dim = lengths(d), dimnames = d)
  z <- array_to_df(z)

  # Data Wrangling
  names(z) <- c("Parameter", "Observations", "Estimator", "Metric", "Value")
  z$Parameter <- as.numeric(as.character(z$Parameter))
  z$Estimator <- factor(z$Estimator)
  z$Metric <- factor(z$Metric)
  z$Observations <- factor(z$Observations, ordered = TRUE)

  # Return the object
  SmallMetrics(D = D, est = est, df = z)

}

setClass("LargeMetrics", slots = list(D = "Distribution",
                                      est = "character",
                                      df = "data.frame"))

#' @title Large Sample Metrics
#' @name LargeMetrics
#' @aliases large_metrics
#'
#' @description
#' This function calculates the asymptotic variance - covariance matrix
#' characterizing the large sample (asymptotic) behavior of an estimator. The
#' function evaluates the metrics as a function of a single parameter, keeping
#' the other ones constant. See Details.
#'
#' @param D A subclass of `Distribution`. The distribution family of interest.
#' @param prm A list containing three elements (name, pos, val). See Details.
#' @param est character. The estimator of interest. Can be a vector.
#' @param df data.frame. a data.frame with columns named "Row", "Col",
#' "Parameter", "Estimator", and "Value".
#' @param ... extra arguments.
#'
#' @inherit small_metrics details
#'
#' @return An object of class `LargeMetrics` with slots `D`, `est`, and `df`.
#'
#' @export
#'
#' @seealso [SmallMetrics], [PlotMetrics]
#' @examples \donttest{
#' # -----------------------------------------------------
#' # Beta Distribution Example
#' # -----------------------------------------------------
#'
#' D <- Beta(shape1 = 1, shape2 = 2)
#'
#' prm <- list(name = "shape1",
#'             val = seq(0.5, 2, by = 0.1))
#'
#' x <- large_metrics(D, prm,
#'                    est = c("mle", "me", "same"))
#'
#' plot(x)
#'
#' # -----------------------------------------------------
#' # Dirichlet Distribution Example
#' # -----------------------------------------------------
#'
#' D <- Dir(alpha = 1:2)
#'
#' prm <- list(name = "alpha",
#'             pos = 1,
#'             val = seq(0.5, 2, by = 0.1))
#'
#' x <- large_metrics(D, prm,
#'                    est = c("mle", "me", "same"))
#'
#' plot(x)
#' }
LargeMetrics <- function(D, est, df) {
  new("LargeMetrics", D = D, est = est, df = df)
}

setValidity("LargeMetrics", function(object) {
  if(!all(names(object@df) %in% c("Row", "Col", "Parameter", "Estimator", "Value"))) {
    stop("df must have columns 'Row', 'Col', 'Parameter, 'Estimator', 'Value'")
  }
  if(!("data.frame" %in% class(object@df))) {
    stop("df has to be a data.frame")
  }
  if (any(c("Row", "Col") %in% names(object@df))) {
    if(!is.factor(object@df$Row)) {
      stop("Column 'Row' has to be factor")
    }
    if(!is.factor(object@df$Col)) {
      stop("Column 'Col' has to be factor")
    }

  }
  if(!is.numeric(object@df$Parameter)) {
    stop("Column 'Parameter' has to be numeric")
  }
  if(!is.factor(object@df$Estimator)) {
    stop("Column 'Estimator' has to be factor")
  }
  if(!is.numeric(object@df$Value)) {
    stop("Column 'Value' has to be numeric")
  }
  TRUE
})

#' @rdname LargeMetrics
#' @export
large_metrics <- function(D,
                          prm,
                          est = c("same", "me", "mle"),
                          ...) {

  if (class(D) %in% c("Cat", "Multinom")) {
    stop("This function is not implemented for the Categorical and Multinomial
    distributions, since changing the probability value of one category forces
    the rest to change as well.")
  }

  # Preliminaries
  Row <- Col <- Parameter <- Estimator <- NULL
  distr <- class(D)[1]
  prm_name <- paste0(prm$name, prm$pos)
  d <- length(get_unknown_params(D, list = FALSE))
  y <- list()

  # Get the distributions
  Di <- lapply(seq_along(prm$val),
               FUN = function(i) { update_params(D, prm, i) })

  # Unidimensional or Multidimensional
  if (d > 1) {
    fvalue <- matrix(0, d, d)
  } else {
    fvalue <- 0
  }

  # For each estimator
  for (j in est) {
    y[[j]] <- vapply(Di, FUN.VALUE = fvalue,
                     FUN = paste0("avar_", j), ...)

    if (d > 1) {
      dimnames(y[[j]])[3] <- list(prm = prm$val)
    } else {
      names(y[[j]]) <- prm$val
    }

  }

  # Create the avar data frame
  y <- simplify2array(y)
  z <- array_to_df(y)

  # Data Wrangling
  if (d > 1) {
    names(z) <- c("Row", "Col", "Parameter", "Estimator", "Value")
    z$Row <- factor(z$Row)
    z$Col <- factor(z$Col)
  } else {
    names(z) <- c("Parameter", "Estimator", "Value")
  }
  z$Parameter <- as.numeric(as.character(z$Parameter))
  z$Estimator <- factor(z$Estimator)

  # Return the object
  LargeMetrics(D = D, est = est, df = z)

}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Plot Metrics
#' @aliases PlotMetrics
#'
#' @description
#' This function provides an easy way to illustrate objects of class
#' `SmallMetrics` and `LargeMetrics`, using the `ggplot2` package. See details.
#'
#' @param x An object of class `SmallMetrics` or `LargeMetrics`.
#' @param y NULL.
#' @param colors character. The colors to be used in the plot.
#' @param title character. The plot title.
#' @param save logical. Should the plot be saved?
#' @param path A path to the directory in which the plot will be saved.
#' @param name character. The name of the output pdf file.
#' @param width numeric. The plot width in inches.
#' @param height numeric. The plot height in inches.
#' @param ... extra arguments.
#'
#' @details
#' Objects of class `SmallMetrics` and `LargeMetrics` are returned by the
#' `small_metrics()` and `large_metrics()` functions, respectively.
#'
#' For the `SmallMetrics`, a grid of line charts is created for each metric and
#' sample size. For the `LargeMetrics`, a grid of line charts is created for
#' each element of the asymptotic variance - covariance matrix.
#'
#' Each estimator is plotted with a different color and line type. The plot can
#' be saved in pdf format.
#'
#' @return The plot is returned invisibly in the form of a `ggplot` object.
#'
#' @importFrom ggplot2 ggplot geom_line aes labs vars
#' @importFrom ggplot2 scale_color_manual theme_minimal theme element_text unit
#' @importFrom ggh4x facet_grid2
#' @export
#'
#' @method plot SmallMetrics,missing
#' @method plot LargeMetrics,missing
#'
#' @seealso [SmallMetrics], [LargeMetrics]
#' @inherit SmallMetrics examples
setGeneric("plot")

#' @rdname plot
setMethod("plot",
          signature(x = "SmallMetrics", y = "missing"),
          function(x,
                   y = NULL,
                   colors = NULL,
                   title = NULL,
                   save = FALSE,
                   path = NULL,
                   name = "myplot.pdf",
                   width = 15,
                   height = 8) {

  # Colors
  if (is.null(colors)) {
    colors <- c("#0073C2", "#CD534C", "#EFC000", "#868686", "#003C67",
                "#7AA6DC", "#A73030", "#8F7700", "#3B3B3B", "#4A6990")

    colors <- colors[seq_along(unique(x@est))]
  }

  # Title
  if (is.null(title)) {
    title <- "Estimator Small Sample Metrics"
  }

  # Preliminaries
  Parameter <- Value <- Estimator <- Observations <- Metric <- NULL

  # Save the plot
  if (save) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    pdf(file.path(path, name), width = width, height = height)
  }

  # Create the plot
  p <- ggplot2::ggplot(x@df) +
    ggplot2::geom_line(ggplot2::aes(x = Parameter,
                                    y = Value,
                                    col = Estimator,
                                    linetype = Estimator),
                       linewidth = 1.5) +
    ggplot2::labs(title = title,
                  y = "Value",
                  x = "Parameter Value") +
    ggh4x::facet_grid2(rows = ggplot2::vars(Observations),
                       cols = ggplot2::vars(Metric),
                       switch = "y") +
    ggplot2::scale_color_manual(values = colors) +
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

})

#' @rdname plot
setMethod("plot",
          signature(x = "LargeMetrics", y = "missing"),
          function(x,
                   y = NULL,
                   colors = NULL,
                   title = NULL,
                   save = FALSE,
                   path = NULL,
                   name = "myplot.pdf",
                   width = 15,
                   height = 8) {

  # Colors
  if (is.null(colors)) {
    colors <- c("#0073C2", "#CD534C", "#EFC000", "#868686", "#003C67",
                "#7AA6DC", "#A73030", "#8F7700", "#3B3B3B", "#4A6990")

    colors <- colors[seq_along(unique(x@est))]
  }

  # Title
  if (is.null(title)) {
    title <- "Estimator Large Sample Metrics"
  }

  # Preliminaries
  Row <- Col <- Parameter <- Estimator <- Value <- NULL

  # Save the plot
  if (save) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(file.path(path, name), width = width, height = height)
  }

  # Create the plot
  p <- ggplot2::ggplot(x@df) +
    ggplot2::geom_line(ggplot2::aes(x = Parameter,
                                    y = Value,
                                    col = Estimator,
                                    linetype = Estimator),
                       linewidth = 1.5) +
    ggplot2::labs(title = title,
                  y = "Value",
                  x = "Parameter Value") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 25),
                   legend.key.size = ggplot2::unit(2, 'cm'),
                   plot.title = ggplot2::element_text(hjust = 0.5))

  if ("Row" %in% names(x@df)) {
    p <- p + ggh4x::facet_grid2(rows = ggplot2::vars(Row),
                                cols = ggplot2::vars(Col),
                                scales = "free",
                                independent = "y")
  }

  plot(p)

  # Close the device
  if (save) {
    grDevices::dev.off()
  }

  # Return the plot
  invisible(p)

})
