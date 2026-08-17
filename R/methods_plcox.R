

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
summary.plcox <- function(object,
                          level = 0.95,
                          ...) {

  ci_info <- credible_interval_info(level)

  beta_draws <- object$beta
  beta_names <- colnames(beta_draws)

  if (is.null(beta_names)) {
    beta_names <- object$colnames
  }

  if (is.null(beta_names)) {
    beta_names <- paste0("beta[", seq_len(ncol(beta_draws)), "]")
  }

  if (length(beta_names) != ncol(beta_draws)) {
    stop(
      "The number of coefficient names does not match ",
      "the number of columns in `object$beta`."
    )
  }

  colnames(beta_draws) <- beta_names

  # --------------------------------
  # Posterior summaries
  # --------------------------------

  beta_mean <- colMeans(beta_draws)

  beta_sd <- apply(
    beta_draws,
    2,
    stats::sd
  )

  beta_ci <- apply(
    beta_draws,
    2,
    stats::quantile,
    probs = ci_info$probs,
    names = FALSE
  )

  # --------------------------------
  # MCMC diagnostics
  # --------------------------------

  beta_diag <- compute_mcmc_diagnostics(
    draws = beta_draws,
    chain_id = object$chain_id,
    variable_names = colnames(beta_draws)
  )

  # --------------------------------
  # Coefficient table
  # --------------------------------

  coef_table <- data.frame(
    Post.Mean = beta_mean,
    Post.SD = beta_sd,
    Lower = beta_ci[1, ],
    Upper = beta_ci[2, ],
    `Exp(Mean)` = exp(beta_mean),
    Rhat = beta_diag$Rhat,
    `Bulk ESS` = beta_diag$`Bulk ESS`,
    `Tail ESS` = beta_diag$`Tail ESS`,
    MCSE = beta_diag$MCSE,
    check.names = FALSE
  )

  names(coef_table)[3:4] <- c(
    ci_info$lower_name,
    ci_info$upper_name
  )

  rownames(coef_table) <- colnames(beta_draws)

  frailty_stat <- NULL
  has_frailty <- object$has_frailty

  if (has_frailty) {

    sig2_draws <- object$sigma2

    sig2_ci <- stats::quantile(
      sig2_draws,
      probs = ci_info$probs,
      names = FALSE
    )

    sig2_diag <- compute_mcmc_diagnostics(
      draws = sig2_draws,
      chain_id = object$chain_id,
      variable_names = "frailty_var"
    )

    sigma2_summary <- data.frame(
      Post.Mean = mean(sig2_draws),
      Post.SD = stats::sd(sig2_draws),
      Lower = sig2_ci[1],
      Upper = sig2_ci[2],
      Rhat = sig2_diag$Rhat,
      `Bulk ESS` = sig2_diag$`Bulk ESS`,
      `Tail ESS` = sig2_diag$`Tail ESS`,
      MCSE = sig2_diag$MCSE,
      check.names = FALSE
    )

    names(sigma2_summary)[3:4] <- c(
      ci_info$lower_name,
      ci_info$upper_name
    )

    rownames(sigma2_summary) <- "frailty_var"

    frailty_stat <- list(
      sigma2 = sigma2_summary,
      n_groups = length(object$group_levels)
    )
  }

  res <- list(
    coefficients = coef_table,
    frailty = frailty_stat,
    has_frailty = has_frailty,
    level = level,
    n_chains = object$n_chains
  )

  class(res) <- "summary.plcox"

  res
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
#' Produces chain-specific trace plots of posterior draws for regression
#' coefficients and, when present, frailty variance parameters and selected
#' individual frailty effects.
#'
#' For models fitted with multiple MCMC chains, each chain is displayed
#' separately in the trace plot.
#'
#' @param x An object of class `"plcox"`.
#' @param intercept Logical; include the intercept in the trace plot when
#'   \code{pars = NULL}?
#' @param pars Regression coefficients to plot. Can be \code{NULL},
#'   numeric indices, or character parameter names. If \code{NULL},
#'   all regression coefficients are plotted, subject to
#'   \code{intercept}.
#' @param include_frailty_var Logical; include the frailty variance
#'   parameter when a frailty model is fitted?
#' @param frailty Optional individual frailty effects to plot. Can be
#'   \code{NULL}, numeric indices, or character names.
#' @param ncol Number of columns used for parameter facets.
#' @param ... Further arguments reserved for future use.
#'
#' @return A \code{ggplot2} object or, when multiple plot components are
#'   produced, a \code{patchwork} object.
#'
#' @export
traceplot.plcox <- function(
    x,
    intercept = FALSE,
    pars = NULL,
    include_frailty_var = FALSE,
    frailty = NULL,
    ncol = 1L,
    ...) {

  .traceplot_bayescox(
    x = x,
    intercept = intercept,
    pars = pars,
    include_frailty_var = include_frailty_var,
    frailty = frailty,
    ncol = ncol,
    ...
  )
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
