
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

  if (is.null(colnames)) {
    return(
      list(
        beta = beta,
        colnames = colnames
      )
    )
  }

  if (!intercept) {
    keep <- colnames != "(Intercept)"
    beta <- beta[, keep, drop = FALSE]
    colnames <- colnames[keep]
  }

  list(
    beta = beta,
    colnames = colnames
  )
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

#' Resolve Regression Parameters for Trace Plots
#'
#' Internal helper function used to determine which regression
#' coefficients should be included in a trace plot.
#'
#' @param x A fitted BayesPLCox model object.
#' @param pars Regression coefficients to plot. Can be \code{NULL},
#'   numeric indices, or character parameter names.
#' @param intercept Logical; should the intercept be included when
#'   \code{pars = NULL}?
#'
#' @return A list containing selected column indices and parameter names.
#'
#' @keywords internal
.resolve_trace_pars <- function(x,
                                pars = NULL,
                                intercept = FALSE) {

  if (is.null(x$beta) || !is.matrix(x$beta)) {
    stop("`x$beta` must be a posterior sample matrix.")
  }

  p <- ncol(x$beta)

  parnames <- colnames(x$beta)

  if (is.null(parnames)) {
    parnames <- x$colnames
  }

  if (is.null(parnames)) {
    parnames <- paste0(
      "beta[",
      seq_len(p),
      "]"
    )
  }

  if (length(parnames) != p) {
    stop(
      "The number of parameter names does not match ",
      "the number of columns in `x$beta`."
    )
  }

  # Default: retain all regression coefficients,
  # optionally excluding the intercept.
  if (is.null(pars)) {

    idx <- seq_len(p)

    if (!isTRUE(intercept)) {
      idx <- idx[
        parnames[idx] != "(Intercept)"
      ]
    }

  } else if (is.numeric(pars)) {

    if (length(pars) == 0L ||
        any(!is.finite(pars)) ||
        any(pars != floor(pars)) ||
        any(pars < 1L) ||
        any(pars > p)) {
      stop(
        "Numeric `pars` must contain valid coefficient indices."
      )
    }

    idx <- as.integer(pars)

  } else if (is.character(pars)) {

    unknown <- setdiff(
      pars,
      parnames
    )

    if (length(unknown) > 0L) {
      stop(
        "Unknown parameter(s) in `pars`: ",
        paste(unknown, collapse = ", "),
        "."
      )
    }

    idx <- match(
      pars,
      parnames
    )

  } else {

    stop(
      "`pars` must be NULL, numeric indices, ",
      "or character parameter names."
    )
  }

  if (length(idx) == 0L) {
    stop(
      "No regression coefficients were selected for plotting."
    )
  }

  list(
    index = idx,
    names = parnames[idx]
  )
}


#' Resolve Frailty Parameters for Trace Plots
#'
#' Internal helper function used to determine which individual
#' frailty effects should be included in a trace plot.
#'
#' @param x A fitted BayesPLCox frailty model object.
#' @param frailty Frailty effects to plot. Can be \code{NULL},
#'   numeric indices, or character names.
#'
#' @return A list containing selected column indices and frailty names.
#'   If \code{frailty = NULL}, an empty selection is returned.
#'
#' @keywords internal
.resolve_trace_frailty <- function(x,
                                   frailty = NULL) {

  if (is.null(frailty)) {
    return(
      list(
        index = integer(0),
        names = character(0)
      )
    )
  }

  if (!isTRUE(x$has_frailty) ||
      is.null(x$u) ||
      !is.matrix(x$u)) {
    stop(
      "`frailty` can only be specified for a fitted frailty model."
    )
  }

  G <- ncol(x$u)

  frailty_names <- colnames(x$u)

  if (is.null(frailty_names) &&
      !is.null(x$group_levels) &&
      length(x$group_levels) == G) {

    frailty_names <- as.character(
      x$group_levels
    )
  }

  if (is.null(frailty_names)) {
    frailty_names <- paste0(
      "frailty[",
      seq_len(G),
      "]"
    )
  }

  if (is.numeric(frailty)) {

    if (length(frailty) == 0L ||
        any(!is.finite(frailty)) ||
        any(frailty != floor(frailty)) ||
        any(frailty < 1L) ||
        any(frailty > G)) {
      stop(
        "Numeric `frailty` must contain valid frailty indices."
      )
    }

    idx <- as.integer(frailty)

  } else if (is.character(frailty)) {

    unknown <- setdiff(
      frailty,
      frailty_names
    )

    if (length(unknown) > 0L) {
      stop(
        "Unknown frailty effect(s): ",
        paste(unknown, collapse = ", "),
        "."
      )
    }

    idx <- match(
      frailty,
      frailty_names
    )

  } else {

    stop(
      "`frailty` must be NULL, numeric indices, ",
      "or character frailty names."
    )
  }

  list(
    index = idx,
    names = frailty_names[idx]
  )
}


#' Internal trace plot engine for Bayesian Cox models
#'
#' Internal helper function that generates chain-specific trace plots
#' of posterior draws for regression coefficients and, when present,
#' frailty variance parameters and selected individual frailty effects.
#' This function is shared by both `"plcox"` and `"gplcox"` methods.
#'
#' @param x A fitted model object containing posterior samples.
#' @param intercept Logical; include the intercept term in the trace plot
#'   when \code{pars = NULL}?
#' @param pars Regression coefficients to plot. Can be \code{NULL},
#'   numeric indices, or character parameter names.
#' @param include_frailty_var Logical; include the frailty variance
#'   parameter when a frailty model is fitted?
#' @param frailty Optional individual frailty effects to plot. Can be
#'   \code{NULL}, numeric indices, or character names.
#' @param ncol Number of columns used for parameter facets.
#' @param ... Further arguments, currently ignored.
#'
#' @return A \code{ggplot2} object or a \code{patchwork} object.
#'
#' @keywords internal
.traceplot_bayescox <- function(
    x,
    intercept = FALSE,
    pars = NULL,
    include_frailty_var = TRUE,
    frailty = NULL,
    ncol = 1L,
    ...) {

  if (!requireNamespace("bayesplot", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Please install 'bayesplot' and 'ggplot2' ",
      "to use this function."
    )
  }

  if (is.null(x$beta) || !is.matrix(x$beta)) {
    stop(
      "`x$beta` must be a posterior sample matrix."
    )
  }

  # --------------------------------------------------
  # Chain information
  # --------------------------------------------------

  if (is.null(x$chain_id)) {

    # Backward compatibility for fitted objects created
    # before multiple-chain support was introduced.
    chain_id <- rep(
      1L,
      nrow(x$beta)
    )

  } else {

    chain_id <- x$chain_id
  }

  if (length(chain_id) != nrow(x$beta)) {
    stop(
      "The length of `x$chain_id` must match ",
      "the number of rows in `x$beta`."
    )
  }

  if (!is.numeric(ncol) ||
      length(ncol) != 1L ||
      !is.finite(ncol) ||
      ncol < 1 ||
      ncol != floor(ncol)) {
    stop(
      "`ncol` must be a positive integer."
    )
  }

  ncol <- as.integer(ncol)

  # --------------------------------------------------
  # Regression coefficients
  # --------------------------------------------------

  beta_sel <- .resolve_trace_pars(
    x = x,
    pars = pars,
    intercept = intercept
  )

  beta <- x$beta[
    ,
    beta_sel$index,
    drop = FALSE
  ]

  colnames(beta) <- beta_sel$names

  beta_array <- as_chain_draws_array(
    draws = beta,
    chain_id = chain_id,
    variable_names = beta_sel$names
  )

  p_beta <- bayesplot::mcmc_trace(
    beta_array,
    facet_args = list(
      ncol = ncol
    )
  ) +
    ggplot2::ggtitle(
      expression(
        "Trace plots: Fixed effects (" * beta * ")"
      )
    )

  plots <- list(
    p_beta
  )

  # --------------------------------------------------
  # Frailty variance
  # --------------------------------------------------

  if (isTRUE(include_frailty_var) &&
      isTRUE(x$has_frailty) &&
      !is.null(x$sigma2)) {

    if (length(x$sigma2) != length(chain_id)) {
      stop(
        "The number of frailty-variance draws is inconsistent ",
        "with `x$chain_id`."
      )
    }

    sig2_array <- as_chain_draws_array(
      draws = x$sigma2,
      chain_id = chain_id,
      variable_names = "sigma2"
    )

    p_sig2 <- bayesplot::mcmc_trace(
      sig2_array
    ) +
      ggplot2::ggtitle(
        expression(
          "Trace plot: Frailty variance (" * sigma^2 * ")"
        )
      )

    plots[[length(plots) + 1L]] <- p_sig2
  }

  # --------------------------------------------------
  # Individual frailty effects
  # --------------------------------------------------

  frailty_sel <- .resolve_trace_frailty(
    x = x,
    frailty = frailty
  )

  if (length(frailty_sel$index) > 0L) {

    if (nrow(x$u) != length(chain_id)) {
      stop(
        "The number of frailty draws is inconsistent ",
        "with `x$chain_id`."
      )
    }

    u_draws <- x$u[
      ,
      frailty_sel$index,
      drop = FALSE
    ]

    colnames(u_draws) <- frailty_sel$names

    u_array <- as_chain_draws_array(
      draws = u_draws,
      chain_id = chain_id,
      variable_names = frailty_sel$names
    )

    p_u <- bayesplot::mcmc_trace(
      u_array,
      facet_args = list(
        ncol = ncol
      )
    ) +
      ggplot2::ggtitle(
        "Trace plots: Frailty effects"
      )

    plots[[length(plots) + 1L]] <- p_u
  }

  # --------------------------------------------------
  # Return plot
  # --------------------------------------------------

  if (length(plots) == 1L) {
    return(
      plots[[1L]]
    )
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Please install 'patchwork' to combine multiple trace plots."
    )
  }

  patchwork::wrap_plots(
    plots,
    ncol = 1L
  )
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
