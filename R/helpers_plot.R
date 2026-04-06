
#' Keep selected terms for plotting
#'
#' Internal helper function used to optionally remove the intercept term from
#' posterior sample matrices before plotting.
#'
#' @param beta Posterior sample matrix for regression coefficients.
#' @param colnames Character vector of coefficient names.
#' @param intercept Logical; should the intercept term be retained?
#'
#' @return A list with filtered posterior draws and corresponding names.
#'
#' @keywords internal
.keep_terms <- function(beta, colnames, intercept = FALSE) {
  if (is.null(colnames)) return(beta)

  if (!intercept) {
    keep <- colnames != "(Intercept)"
    beta <- beta[, keep, drop = FALSE]
    colnames <- colnames[keep]
  }

  list(beta = beta, colnames = colnames)
}

#' Summarize posterior samples for plotting
#'
#' Internal helper function that computes posterior means and credible intervals
#' for regression coefficients and returns them in a data-frame format suitable
#' for plotting.
#'
#' @param beta_draws Posterior sample matrix.
#' @param parnames Optional character vector of parameter names.
#' @param level Credible interval level.
#' @param transform Transformation applied to summaries. One of `"none"` or `"exp"`.
#'
#' @return A data frame with posterior summaries for plotting.
#'
#' @keywords internal
.posterior_summary_df <- function(beta_draws, parnames = NULL, level = 0.95,
                                  transform = c("none", "exp")) {
  transform <- match.arg(transform)

  if (!is.matrix(beta_draws)) stop("beta_draws must be a matrix.")

  alpha <- (1 - level) / 2
  probs <- c(alpha, 1 - alpha)

  est <- colMeans(beta_draws)
  ci  <- apply(beta_draws, 2, stats::quantile, probs = probs)

  # transform
  if (transform == "exp") {
    est <- exp(est)
    ci  <- exp(ci)
  }

  out <- data.frame(
    term = if (is.null(parnames)) colnames(beta_draws) else parnames,
    estimate = est,
    lower = ci[1, ],
    upper = ci[2, ],
    stringsAsFactors = FALSE
  )

  out$term <- factor(out$term, levels = rev(out$term))
  out
}

#' Internal trace plot engine for Bayesian Cox models
#'
#' Internal helper function that generates trace plots of posterior draws
#' for regression coefficients and, when present, frailty variance parameters.
#' This function is shared by both `"plcox"` and `"gplcox"` methods.
#'
#' @param x A fitted model object containing posterior samples.
#' @param intercept Logical; include the intercept term in the trace plot?
#' @param ... Further arguments, currently ignored.
#'
#' @return A plot object (typically a `ggplot2` or `patchwork` object).
#'
#' @keywords internal
.traceplot_bayescox <- function(x, intercept = FALSE, ...) {
  if (!requireNamespace("bayesplot", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Please install 'bayesplot', 'ggplot2', and 'patchwork' to use this function.")
  }

  if (is.null(x$beta) || !is.matrix(x$beta)) {
    stop("x$beta must be a posterior sample matrix.")
  }

  beta <- x$beta
  parnames <- colnames(beta)
  if (is.null(parnames)) parnames <- x$colnames

  tmp <- .keep_terms(beta, parnames, intercept = intercept)
  beta <- tmp$beta
  parnames <- tmp$colnames

  colnames(beta) <- parnames

  p1 <- bayesplot::mcmc_trace(beta) +
    ggplot2::ggtitle(expression("Trace plots: Fixed effects (" * beta * ")"))

  if (!is.null(x$sigma2)) {
    draws_sig2 <- matrix(x$sigma2, ncol = 1)
    colnames(draws_sig2) <- "sigma2"

    p2 <- bayesplot::mcmc_trace(draws_sig2) +
      ggplot2::ggtitle(expression("Trace plot: Frailty variance (" * sigma^2 * ")"))

    return(p1 / p2)
  } else {
    return(p1)
  }
}

#' Internal posterior interval plot engine for Bayesian Cox models
#'
#' Internal helper function that produces posterior interval plots for regression
#' coefficients and, when present, frailty variance parameters. This function is
#' shared by both `"plcox"` and `"gplcox"` methods.
#'
#' @param x A fitted model object containing posterior samples.
#' @param level Credible interval level.
#' @param transform Transformation applied to coefficient summaries. One of
#'   `"none"` or `"exp"`.
#' @param intercept Logical; include the intercept term in the plot?
#' @param ... Further arguments, currently ignored.
#'
#' @return A `ggplot2` object or a `patchwork` object.
#'
#' @keywords internal
.plot_bayescox <- function(x, level = 0.95,
                           transform = c("none", "exp"),
                           intercept = FALSE, ...) {
  transform <- match.arg(transform)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2' to use this function.")
  }

  if (is.null(x$beta) || !is.matrix(x$beta)) {
    stop("x$beta must be a posterior sample matrix.")
  }

  beta <- x$beta
  parnames <- colnames(beta)
  if (is.null(parnames)) parnames <- x$colnames

  tmp <- .keep_terms(beta, parnames, intercept = intercept)
  beta <- tmp$beta
  parnames <- tmp$colnames

  df_beta <- .posterior_summary_df(
    beta,
    parnames = parnames,
    level = level,
    transform = transform
  )

  xlab <- if (transform == "exp") {
    "Posterior mean and credible interval (exp scale)"
  } else {
    "Posterior mean and credible interval"
  }

  p1 <- ggplot2::ggplot(df_beta, ggplot2::aes(x = .data$estimate, y = .data$term)) +
    ggplot2::geom_vline(
      xintercept = if (transform == "exp") 1 else 0,
      linetype = 2, colour = "grey60"
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      height = 0.2
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = xlab,
      y = NULL,
      title = "Posterior interval plot: fixed effects"
    ) +
    ggplot2::theme_bw()

  if (!is.null(x$sigma2)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) return(p1)

    df_sig2 <- data.frame(
      term = factor("sigma2", levels = "sigma2"),
      estimate = mean(x$sigma2),
      lower = stats::quantile(x$sigma2, (1 - level) / 2),
      upper = stats::quantile(x$sigma2, 1 - (1 - level) / 2)
    )

    p2 <- ggplot2::ggplot(df_sig2, ggplot2::aes(x = .data$estimate, y = .data$term)) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
        height = 0.2
      ) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(
        x = "Posterior mean and credible interval",
        y = NULL,
        title = expression("Frailty variance (" * sigma^2 * ")")
      ) +
      ggplot2::theme_bw()

    return(p1 / p2)
  } else {
    return(p1)
  }
}
