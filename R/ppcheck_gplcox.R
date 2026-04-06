
#' Posterior predictive check for tie block sizes in GPL-Cox
#'
#' Generates posterior predictive summaries of tie block sizes across distinct
#' event times under a fitted GPL-Cox model.
#'
#' @param object An object of class `"gplcox"`.
#' @param ndraw Number of posterior draws used for the predictive check. If
#'   `NULL`, all stored draws are used.
#' @param level Predictive interval level.
#' @param seed Optional random seed for reproducibility.
#' @param ... Further arguments, currently ignored.
#'
#' @return An object of class `"pp_check_tie_gplcox"` and `"data.frame"`.
#'
#' @examples
#' \dontrun{
#' fit <- gplcox(Surv(time, status) ~ x1 + x2, data = dat)
#' pp <- pp_check_tie.gplcox(fit)
#' print(pp)
#' plot(pp)
#' }
#'
#' @export
pp_check_tie.gplcox <- function(object,
                                ndraw = NULL,
                                level = 0.95,
                                seed = NULL,
                                ...) {
  if (!inherits(object, "gplcox")) {
    stop("object must be of class 'gplcox'.")
  }

  needed <- c("beta", "X", "R_sets", "d_obs")
  miss <- needed[!needed %in% names(object)]
  if (length(miss) > 0) {
    stop("The object does not contain required components: ",
         paste(miss, collapse = ", "),
         ". Refit the model with stored GPL risk-set information.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  beta_draws <- object$beta
  n_save <- nrow(beta_draws)

  if (is.null(ndraw)) {
    idx_draw <- seq_len(n_save)
  } else {
    ndraw <- min(ndraw, n_save)
    idx_draw <- sort(sample.int(n_save, ndraw))
  }

  R <- length(object$R_sets)
  d_rep <- matrix(NA_integer_, nrow = length(idx_draw), ncol = R)

  has_frailty <- isTRUE(object$has_frailty) && !is.null(object$u)

  for (m in seq_along(idx_draw)) {
    s <- idx_draw[m]

    beta_m <- beta_draws[s, ]

    if (has_frailty) {
      u_m <- object$u[s, ]
      d_rep[m, ] <- .sim_tie_block_size_once(
        beta = beta_m,
        X = object$X,
        R_sets = object$R_sets,
        u = u_m,
        group_id_ord = object$group_id_ord
      )
    } else {
      d_rep[m, ] <- .sim_tie_block_size_once(
        beta = beta_m,
        X = object$X,
        R_sets = object$R_sets
      )
    }
  }

  alpha <- (1 - level) / 2

  out <- data.frame(
    event_index  = seq_len(R),
    observed     = object$d_obs,
    pred_mean    = colMeans(d_rep),
    pred_median  = apply(d_rep, 2, stats::median),
    pred_lower   = apply(d_rep, 2, stats::quantile, probs = alpha),
    pred_upper   = apply(d_rep, 2, stats::quantile, probs = 1 - alpha),
    pp_p_ge_obs  = colMeans(sweep(d_rep, 2, object$d_obs, FUN = ">=")),
    pp_p_le_obs  = colMeans(sweep(d_rep, 2, object$d_obs, FUN = "<="))
  )

  attr(out, "replicates") <- d_rep
  class(out) <- c("pp_check_tie_gplcox", "data.frame")
  out
}

#' Print a posterior predictive tie-block check
#'
#' @param x An object of class `"pp_check_tie_gplcox"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return The object, invisibly.
#' @export
print.pp_check_tie_gplcox <- function(x, ...) {
  cat("Posterior predictive check for tie block sizes\n")
  cat("Number of event times:", nrow(x), "\n\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' Plot posterior predictive tie-block summaries
#'
#' Produces a posterior predictive check plot comparing observed tie block sizes
#' with posterior predictive medians and intervals.
#'
#' @param x An object of class `"pp_check_tie_gplcox"`.
#' @param y Unused.
#' @param show_mean Logical; add the posterior predictive mean?
#' @param connect Logical; connect observed and predictive median values by lines?
#' @param subtitle Logical; include a subtitle explaining the plot elements?
#' @param ... Further arguments, currently ignored.
#'
#' @return A `ggplot2` object.
#' @export
plot.pp_check_tie_gplcox <- function(x, y,
                                     show_mean = FALSE,
                                     connect = FALSE,
                                     subtitle = FALSE,
                                     ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2' to use this function.")
  }

  df <- as.data.frame(x)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$event_index)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$pred_lower, ymax = .data$pred_upper),
      width = 0.2,
      colour = "#2C7FB8"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$pred_median),
      colour = "#2C7FB8",
      size = 2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$observed),
      colour = "black",
      size = 2
    ) +
    ggplot2::labs(
      x = "Index of distinct event time",
      y = "Tie block size",
      title = "Posterior predictive check for tie block sizes"
    ) +
    ggplot2::theme_bw()

  if (subtitle) {
    p <- p + ggplot2::labs(
      subtitle = "Black = observed, blue = posterior predictive median and 95% interval"
    )
  }

  if (connect) {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$observed),
        colour = "black",
        linewidth = 0.4
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$pred_median),
        colour = "#2C7FB8",
        linewidth = 0.4
      )
  }

  if (show_mean) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes(y = .data$pred_mean),
      colour = "#2C7FB8",
      linetype = 2,
      linewidth = 0.5
    )
  }

  return(p)
}
