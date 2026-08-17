
#' MCMC Diagnostics
#'
#' Computes convergence and sampling-efficiency diagnostics for a
#' fitted BayesPLCox model.
#'
#' @param object A fitted model object.
#' @param include_frailty Logical; whether individual frailty effects
#'   should also be included.
#' @param ... Further arguments.
#'
#' @return A list of MCMC diagnostic tables.
#'
#' @export
mcmc_diagnostics <- function(object,
                             include_frailty = FALSE,
                             ...) {
  UseMethod("mcmc_diagnostics")
}


#' MCMC Diagnostics for BayesPLCox Models
#'
#' @param object A fitted object inheriting from class
#'   \code{"bayescoxrank"}.
#' @param include_frailty Logical; whether diagnostics for individual
#'   frailty effects should be included.
#' @param ... Further arguments.
#'
#' @return A list containing MCMC diagnostics for regression
#'   coefficients and, when applicable, frailty components.
#'
#' @export
mcmc_diagnostics.bayescoxrank <- function(
    object,
    include_frailty = FALSE,
    ...) {

  out <- list(
    coefficients = NULL,
    frailty_var = NULL,
    frailty = NULL
  )

  # Regression coefficients
  beta_names <- colnames(object$beta)

  if (is.null(beta_names)) {
    beta_names <- object$colnames
  }

  if (is.null(beta_names)) {
    beta_names <- paste0(
      "beta[",
      seq_len(ncol(object$beta)),
      "]"
    )
  }

  out$coefficients <- compute_mcmc_diagnostics(
    draws = object$beta,
    chain_id = object$chain_id,
    variable_names = beta_names
  )

  if (object$has_frailty) {

    # Frailty variance
    out$frailty_var <- compute_mcmc_diagnostics(
      draws = object$sigma2,
      chain_id = object$chain_id,
      variable_names = "frailty_var"
    )

    # Individual frailty effects
    if (include_frailty) {

      frailty_names <- colnames(object$u)

      if (is.null(frailty_names)) {
        frailty_names <- paste0(
          "frailty[",
          seq_len(ncol(object$u)),
          "]"
        )
      }

      out$frailty <- compute_mcmc_diagnostics(
        draws = object$u,
        chain_id = object$chain_id,
        variable_names = frailty_names
      )
    }
  }

  out
}
