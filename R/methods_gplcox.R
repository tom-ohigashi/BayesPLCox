
#' Print a fitted GPL-Cox model
#'
#' Prints a concise summary of a fitted `"gplcox"` object.
#'
#' @param x An object of class `"gplcox"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return The fitted model object, invisibly.
#' @export
print.gplcox <- function(x, ...) {
  cat("Bayesian Cox model fitted via", x$method, "\n")
  cat("Frailty:", if (x$has_frailty) "yes" else "no", "\n")
  cat("Posterior samples:", x$n_save, "\n")
  cat("Parameters:", x$n_par, "\n")
  cat("Call:\n")
  print(x$call)
  invisible(x)
}

#' Summarize a fitted GPL-Cox model
#'
#' Computes posterior summaries for regression coefficients from a fitted
#' `"gplcox"` object.
#'
#' @param object An object of class `"gplcox"`.
#' @param level Credible interval level.
#' @param ... Further arguments, currently ignored.
#'
#' @return An object of class `"summary.gplcox"`.
#' @export
summary.gplcox <- function(object, level = 0.95, ...) {
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
  if (object$has_frailty) {
    frailty_stat <- list(
      sigma2 = c(Mean = mean(object$sigma2), SD = sd(object$sigma2),
                 quantile(object$sigma2, c(0.025, 0.975))),
      n_groups = length(object$group_levels)
    )
  }

  res <- list(coefficients = coef_table, frailty = frailty_stat,
              has_frailty = object$has_frailty, method = object$method)
  class(res) <- "summary.gplcox"
  return(res)
}

#' Print a summary of a GPL-Cox model
#'
#' @param x An object of class `"summary.gplcox"`.
#' @param digits Number of digits to print.
#' @param ... Further arguments, currently ignored.
#'
#' @return The summary object, invisibly.
#' @export
print.summary.gplcox <- function(x, digits = 3, ...) {
  cat("Fixed Effects:\n")
  print(round(x$coefficients, 4))
  if (x$has_frailty) {
    cat("\nFrailty Variance (Sigma^2):\n")
    print(round(x$frailty$sigma2, 4))
  }
}


#' Extract posterior mean coefficients from a GPL-Cox model
#'
#' @param object An object of class `"gplcox"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return A named numeric vector of posterior mean coefficients.
#' @export
coef.gplcox <- function(object, ...) {
  colMeans(object$beta)
}

#' Credible intervals for a GPL-Cox model
#'
#' Computes equal-tailed posterior credible intervals for regression coefficients.
#'
#' @param object An object of class `"gplcox"`.
#' @param parm Parameter indices or names. Currently ignored.
#' @param level Credible interval level.
#' @param ... Further arguments, currently ignored.
#'
#' @return A matrix of posterior credible intervals.
#' @export
confint.gplcox <- function(object, parm, level = 0.95, ...) {
  alpha <- (1 - level) / 2
  out <- t(apply(object$beta, 2, quantile, probs = c(alpha, 1 - alpha)))
  colnames(out) <- c("lower", "upper")
  out
}

#' Trace plots for a fitted GPL-Cox model
#'
#' Produces trace plots of posterior draws for regression coefficients and, when
#' present, frailty variance parameters.
#'
#' @param x An object of class `"gplcox"`.
#' @param intercept Logical; include the intercept in the trace plot?
#' @param ... Further arguments, currently ignored.
#'
#' @return A plot object.
#' @export
traceplot.gplcox <- function(x, intercept = FALSE, ...) {
  .traceplot_bayescox(x, intercept = intercept, ...)
}

#' Posterior interval plot for a fitted GPL-Cox model
#'
#' Produces a posterior interval plot for regression coefficients from a fitted
#' `"gplcox"` object.
#'
#' @param x An object of class `"gplcox"`.
#' @param y Unused.
#' @param level Credible interval level.
#' @param transform Transformation applied to coefficient summaries. One of
#'   `"none"` or `"exp"`.
#' @param intercept Logical; include the intercept in the plot?
#' @param ... Further arguments, currently ignored.
#'
#' @return A `ggplot2` object or a `patchwork` object.
#' @export
plot.gplcox <- function(x, y, level = 0.95,
                        transform = c("none", "exp"),
                        intercept = FALSE, ...) {
  .plot_bayescox(x, level = level, transform = transform,
                 intercept = intercept, ...)
}
