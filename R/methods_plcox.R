

#' Print a fitted PL-Cox model
#'
#' Prints a concise summary of a fitted `"plcox"` object.
#'
#' @param x An object of class `"plcox"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return The fitted model object, invisibly.
#' @export
print.plcox <- function(x, ...) {
  cat("Bayesian Cox model fitted via", x$method, "\n")
  cat("Frailty:", if (x$has_frailty) "yes" else "no", "\n")
  cat("Posterior samples:", x$n_save, "\n")
  cat("Parameters:", x$n_par, "\n")
  cat("Call:\n")
  print(x$call)
  invisible(x)
}

#' Summarize a fitted PL-Cox model
#'
#' Computes posterior summaries for regression coefficients from a fitted
#' `"plcox"` object.
#'
#' @param object An object of class `"plcox"`.
#' @param level Credible interval level.
#' @param ... Further arguments, currently ignored.
#'
#' @return An object of class `"summary.plcox"`.
#' @export
summary.plcox <- function(object, level = 0.95, ...) {
  beta_draws <- object$beta
  beta_mean  <- colMeans(beta_draws)
  beta_sd    <- apply(beta_draws, 2, sd)
  beta_ci    <- apply(beta_draws, 2, quantile, probs = c(0.025, 0.975))

  if (requireNamespace("coda", quietly = TRUE)) {
    beta_ess <- coda::effectiveSize(coda::as.mcmc(beta_draws))
  } else {
    beta_ess <- rep(NA, length(beta_mean))
  }

  coef_table <- data.frame(
    Post.Mean = beta_mean,
    Post.SD   = beta_sd,
    `Lower 95%` = beta_ci[1, ],
    `Upper 95%` = beta_ci[2, ],
    `Exp(Mean)` = exp(beta_mean),
    ESS       = beta_ess,
    check.names = FALSE
  )
  rownames(coef_table) <- object$colnames

  frailty_stat <- NULL
  has_frailty  <- !is.null(object$sigma2)

  if (has_frailty) {
    sig2_draws <- object$sigma2
    frailty_stat <- list(
      sigma2 = c(
        Mean = mean(sig2_draws),
        SD   = sd(sig2_draws),
        `2.5%` = as.numeric(quantile(sig2_draws, 0.025)),
        `97.5%` = as.numeric(quantile(sig2_draws, 0.975))
      ),
      n_groups = length(object$group_levels)
    )
  }

  res <- list(
    coefficients = coef_table,
    frailty      = frailty_stat,
    has_frailty  = has_frailty
  )

  class(res) <- "summary.plcox"
  return(res)
}

#' Print a summary of a PL-Cox model
#'
#' @param x An object of class `"summary.plcox"`.
#' @param digits Number of digits to print.
#' @param ... Further arguments, currently ignored.
#'
#' @return The summary object, invisibly.
#' @export
print.summary.plcox <- function(x, digits = 3, ...) {
  cat("Fixed Effects (Coefficients):\n")
  print(round(x$coefficients, 4))

  if (x$has_frailty) {
    cat("\nFrailty Variance (Sigma^2):\n")
    print(round(x$frailty$sigma2, 4))
  }
}

#' Extract posterior mean coefficients from a PL-Cox model
#'
#' @param object An object of class `"plcox"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return A named numeric vector of posterior mean coefficients.
#' @export
coef.plcox <- function(object, ...) {
  colMeans(object$beta)
}

#' Extract posterior mean coefficients from a PL-Cox model
#'
#' @param object An object of class `"plcox"`.
#' @param parm Parameter indices or names. Currently ignored.
#' @param level Credible interval level.
#' @param ... Further arguments, currently ignored.
#'
#' @return A named numeric vector of posterior mean coefficients.
#' @export
confint.plcox <- function(object, parm, level = 0.95, ...) {
  alpha <- (1 - level) / 2
  out <- t(apply(object$beta, 2, quantile, probs = c(alpha, 1 - alpha)))
  colnames(out) <- c("lower", "upper")
  out
}

#' Trace plots for a fitted PL-Cox model
#'
#' Produces trace plots of posterior draws for regression coefficients and, when
#' present, frailty variance parameters.
#'
#' @param x An object of class `"plcox"`.
#' @param intercept Logical; include the intercept in the trace plot?
#' @param ... Further arguments, currently ignored.
#'
#' @return A plot object.
#' @export
traceplot.plcox <- function(x, intercept = FALSE, ...) {
  .traceplot_bayescox(x, intercept = intercept, ...)
}

#' Posterior interval plot for a fitted PL-Cox model
#'
#' Produces a posterior interval plot for regression coefficients from a fitted
#' `"plcox"` object.
#'
#' @param x An object of class `"plcox"`.
#' @param y Unused.
#' @param level Credible interval level.
#' @param transform Transformation applied to coefficient summaries. One of
#'   `"none"` or `"exp"`.
#' @param intercept Logical; include the intercept in the plot?
#' @param ... Further arguments, currently ignored.
#'
#' @return A `ggplot2` object or a `patchwork` object.
#' @export
plot.plcox <- function(x, y, level = 0.95,
                       transform = c("none", "exp"),
                       intercept = FALSE, ...) {
  .plot_bayescox(x, level = level, transform = transform,
                 intercept = intercept, ...)
}
