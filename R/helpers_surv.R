#' Construct risk-set information for PL-Cox
#'
#' Internal helper function that constructs risk sets and related indexing
#' objects from survival data for the PL-Cox model.
#'
#' @param time Numeric vector of observed times.
#' @param status Integer or numeric event indicator (`1` = event, `0` = censored).
#'
#' @return A list containing event times, risk sets, membership indices,
#'   tie counts, and event indicators.
#'
#' @keywords internal
build_pl_from_surv <- function(time, status) {
  n <- length(time)
  event_times <- sort(unique(time[status == 1]))
  num_events <- length(event_times)

  # Individual indices for risk set R_t at each event time t
  # and the set of risk sets containing individual i (t's such that i \in R_t)
  risk_sets <- vector("list", num_events)
  event_sets <- vector("list", num_events)
  member_of <- vector("list", n)

  # Number of concurrent events |E_r| per event time (handling ties)
  e_counts <- numeric(num_events)

  for (t in seq_along(event_times)) {
    time_t <- event_times[t]
    r_idx <- which(time >= time_t)
    e_idx <- which(time == time_t & status == 1)
    risk_sets[[t]] <- r_idx
    event_sets[[t]] <- e_idx
    e_counts[t] <- sum(time == time_t & status == 1)

    for (idx in r_idx) {
      member_of[[idx]] <- c(member_of[[idx]], t)
    }
  }

  return(list(
    event_times = event_times,
    risk_sets = risk_sets,
    event_sets = event_sets,
    member_of = member_of,
    e_counts = e_counts,
    d = status # Event indicator for each individual (0 or 1)
  ))
}

#' Construct risk-set and tied-event information for GPL-Cox
#'
#' Internal helper function that constructs ordered risk sets and concurrent
#' event sets from survival data for the GPL-Cox model.
#'
#' @param time Numeric vector of observed times.
#' @param status Integer or numeric event indicator (`1` = event, `0` = censored).
#' @param X Numeric design matrix.
#'
#' @return A list containing ordered survival data, risk sets, event sets,
#'   event indicators, and the ordering index.
#'
#' @keywords internal
build_gpl_from_surv <- function(time, status, X) {
  stopifnot(length(time) == length(status), nrow(X) == length(time))
  n <- length(time)
  ord <- order(time, decreasing = FALSE)       # Sort by time (ascending)
  time  <- time[ord]; status <- status[ord]; X <- X[ord, , drop=FALSE]

  # Unique time points where events occurred
  event_times <- sort(unique(time[status == 1]))

  R_sets <- vector("list", length(event_times))  # Risk sets R_t
  E_sets <- vector("list", length(event_times))  # Concurrent event sets E_t

  for (t in seq_along(event_times)) {
    time_t <- event_times[t]
    # Risk set: Individuals who have not dropped out (time >= t_t)
    R_t <- which(time >= time_t)
    # Concurrent events: Individuals who experienced the event at t_t
    E_t <- which(time == time_t & status == 1)
    R_sets[[t]] <- R_t
    E_sets[[t]] <- E_t
  }

  # Selection count d_i for each individual (= event count; usually 0 or 1)
  d <- integer(n)
  for (t in seq_along(E_sets)) d[E_sets[[t]]] <- d[E_sets[[t]]] + 1L

  list(
    time  = time,
    status= status,
    X     = X,
    R_sets= R_sets,
    E_sets= E_sets,
    d     = d,
    ord   = ord            # Used for reverting to original order
  )
}

#' Sample latent geometric variables for GPL-Cox
#'
#' Internal helper function that samples latent geometric variables associated
#' with each risk set under the GPL-Cox model.
#'
#' @param theta Vector of success probabilities.
#' @param R_sets A list of risk sets.
#'
#' @return An integer vector of latent geometric variables.
#'
#' @keywords internal
sample_Z_geom <- function(theta, R_sets) {
  Z <- integer(length(R_sets))
  for (t in seq_along(R_sets)) {
    Rt <- R_sets[[t]]
    p_t <- 1 - prod(1 - theta[Rt])   # Probability of at least one success
    # Numerical stabilization (prevent p_r from reaching 0 or 1)
    p_t <- pmin(pmax(p_t, 1e-10), 1 - 1e-10)
    Z[t] <- 1L + rgeom(1L, p_t)      # rgeom returns # of failures; +1 for # of trials
  }
  Z
}

#' Accumulate latent geometric counts over risk sets
#'
#' Internal helper function that computes the cumulative latent count for each
#' subject by summing latent geometric variables over all risk sets containing
#' that subject.
#'
#' @param Z Integer vector of latent geometric variables.
#' @param R_sets A list of risk sets.
#' @param n Number of subjects.
#'
#' @return An integer vector of cumulative latent counts.
#'
#' @keywords internal
accumulate_zeta_i <- function(Z, R_sets, n) {
  zeta <- integer(n)
  for (t in seq_along(R_sets)) {
    zeta[R_sets[[t]]] <- zeta[R_sets[[t]]] + Z[t]
  }
  zeta
}

#' Simulate posterior predictive tie block sizes once
#'
#' Internal helper function that simulates one posterior predictive replicate
#' of tie block sizes across event times under the GPL-Cox model.
#'
#' @param beta Numeric vector of regression coefficients.
#' @param X Numeric design matrix ordered consistently with `R_sets`.
#' @param R_sets A list of risk sets.
#' @param u Optional numeric vector of frailty terms.
#' @param group_id_ord Optional integer vector of ordered group memberships.
#'
#' @return An integer vector of posterior predictive tie block sizes.
#'
#' @keywords internal
.sim_tie_block_size_once <- function(beta, X, R_sets,
                                     u = NULL, group_id_ord = NULL) {
  eta <- as.numeric(X %*% beta)
  if (!is.null(u)) {
    eta <- eta + u[group_id_ord]
  }
  theta <- stats::plogis(eta)
  out <- integer(length(R_sets))
  for (r in seq_along(R_sets)) {
    idx <- R_sets[[r]]
    p_r <- theta[idx]

    p_r <- pmin(pmax(p_r, 1e-10), 1 - 1e-10)

    w_r <- stats::rgeom(length(idx), prob = p_r) + 1L
    out[r] <- sum(w_r == min(w_r))
  }
  out
}

#' Log-likelihood for the PL-Cox model
#'
#' Computes the log partial likelihood of the Cox model under the
#' Plackett--Luce (PL) formulation given posterior draws of regression
#' coefficients.
#'
#' The likelihood is evaluated based on risk sets and event sets derived
#' from survival data. Tied events are handled via multiplicity in
#' \code{event_sets}.
#'
#' @param beta Numeric vector of regression coefficients.
#' @param X Numeric design matrix.
#' @param risk_sets A list of integer vectors. Each element contains the indices
#'   of individuals in the risk set at a given event time.
#' @param event_sets A list of integer vectors. Each element contains the indices
#'   of individuals experiencing the event at that time.
#'
#' @return A numeric value representing the log partial likelihood.
#'
#' @details
#' The log-likelihood is computed as
#' \deqn{
#' \sum_r \left( \sum_{i \in E_r} \log \lambda_i
#' - d_r \log \sum_{j \in R_r} \lambda_j \right),
#' }
#' where \eqn{\lambda_i = \exp(x_i^\top \beta)} and \eqn{d_r = |E_r|}.
#'
#' @keywords internal
loglik_plcox <- function(beta, X, risk_sets, event_sets) {
  eta <- as.vector(X %*% beta)
  lambda <- exp(eta)
  out <- 0
  R <- length(risk_sets)
  for (r in seq_len(R)) {
    Er <- event_sets[[r]]
    Rr <- risk_sets[[r]]
    dr <- length(Er)
    out <- out + sum(log(lambda[Er])) - dr * log(sum(lambda[Rr]))
  }
  out
}

#' Log-likelihood for the GPL-Cox model
#'
#' Computes the log-likelihood of the Cox model under the generalized
#' Plackett--Luce (GPL) formulation.
#'
#' The GPL model introduces a geometric latent variable representation,
#' leading to a likelihood that accounts for tied events through success
#' probabilities.
#'
#' @param beta Numeric vector of regression coefficients.
#' @param X Numeric design matrix.
#' @param risk_sets A list of integer vectors representing risk sets.
#' @param event_sets A list of integer vectors representing event sets.
#'
#' @return A numeric value representing the GPL-Cox log-likelihood.
#'
#' @details
#' Let \eqn{\theta_i = \mathrm{logit}^{-1}(x_i^\top \beta)}.
#' The log-likelihood is given by
#' \deqn{
#' \sum_r \left(
#' \sum_{i \in E_r} \log \theta_i
#' + \sum_{i \in R_r \setminus E_r} \log (1 - \theta_i)
#' - \log \left( 1 - \prod_{j \in R_r} (1 - \theta_j) \right)
#' \right).
#' }
#'
#' Numerical stabilization is achieved using \code{log1p}.
#'
#' @keywords internal
loglik_gplcox <- function(beta, X, risk_sets, event_sets) {
  eta <- as.vector(X %*% beta)
  theta <- plogis(eta)
  out <- 0
  R <- length(risk_sets)
  for (r in seq_len(R)) {
    Er <- event_sets[[r]]
    Rr <- risk_sets[[r]]
    nonEr <- setdiff(Rr, Er)

    num1 <- sum(log(theta[Er]))
    num2 <- if (length(nonEr) > 0) sum(base::log1p(-theta[nonEr])) else 0

    s <- sum(base::log1p(-theta[Rr]))
    den <- base::log1p(-exp(s))

    out <- out + num1 + num2 - den
  }
  out
}

#' Compute the deviance information criterion from posterior draws
#'
#' Computes the deviance information criterion (DIC) from posterior draws
#' of regression coefficients for either the PL-Cox or GPL-Cox model.
#'
#' Risk-set and event-set information are constructed internally from the
#' observed survival outcome, so users only need to provide posterior draws,
#' the survival outcome, and the design matrix.
#'
#' @param beta_draws A numeric matrix of posterior draws, with rows
#'   corresponding to MCMC samples and columns corresponding to model
#'   parameters.
#' @param time A numeric vector of observed times.
#' @param event An integer or numeric event indicator
#'   (\code{1} = event, \code{0} = censored).
#' @param X A numeric design matrix.
#' @param loglik Character string specifying the likelihood to use:
#'   \code{"pl"} (default) for the PL-Cox likelihood or
#'   \code{"gpl"} for the GPL-Cox likelihood.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Dbar}{Posterior mean deviance.}
#'   \item{Dhat}{Deviance evaluated at the posterior mean of \code{beta}.}
#'   \item{pD}{Effective number of parameters, defined as
#'   \code{Dbar - Dhat}.}
#'   \item{DIC}{Deviance information criterion, defined as
#'   \code{Dbar + pD}.}
#'   \item{deviance}{Numeric vector of deviance values for all posterior draws.}
#' }
#'
#' @details
#' The deviance is defined as
#' \deqn{D(\beta) = -2 \log L(\beta),}
#' where \eqn{L(\beta)} is the model-specific likelihood.
#'
#' The DIC is computed as
#' \deqn{\mathrm{DIC} = \bar{D} + p_D,}
#' where \eqn{\bar{D}} is the posterior mean deviance and
#' \eqn{p_D = \bar{D} - D(\bar{\beta})}.
#'
#' For the GPL-Cox likelihood, the survival data and design matrix are
#' internally reordered by increasing observed time so that the event-set
#' indexing is consistent with the likelihood evaluation.
#'
#' @examples
#' \dontrun{
#' dic_from_draws(
#'   beta_draws = fit_pl$beta,
#'   time = readmission$time,
#'   event = readmission$event,
#'   X = X,
#'   loglik = "pl"
#' )
#'
#' dic_from_draws(
#'   beta_draws = fit_gpl$beta,
#'   time = readmission$time,
#'   event = readmission$event,
#'   X = X,
#'   loglik = "gpl"
#' )
#' }
#'
#' @export
dic_from_draws <- function(beta_draws, time, event, X,
                           loglik = c("pl", "gpl")) {
  loglik <- match.arg(loglik)

  if (!is.matrix(beta_draws)) {
    beta_draws <- as.matrix(beta_draws)
  }

  if (!is.numeric(time)) {
    stop("'time' must be numeric.")
  }
  if (!is.numeric(event) && !is.integer(event)) {
    stop("'event' must be numeric or integer.")
  }
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (length(time) != length(event)) {
    stop("'time' and 'event' must have the same length.")
  }
  if (nrow(X) != length(time)) {
    stop("nrow(X) must equal length(time).")
  }

  if (loglik == "pl") {
    dat <- build_pl_from_surv(time = time, status = event)
    X_use <- X
    risk_sets <- dat$risk_sets
    event_sets <- dat$event_sets
    loglik_fun <- loglik_plcox
  } else {
    dat <- build_gpl_from_surv(time = time, status = event, X = X)
    X_use <- dat$X
    risk_sets <- dat$R_sets
    event_sets <- dat$E_sets
    loglik_fun <- loglik_gplcox
  }

  dev <- apply(beta_draws, 1, function(b) {
    -2 * loglik_fun(b, X_use, risk_sets, event_sets)
  })

  beta_bar <- colMeans(beta_draws)
  d_bar <- mean(dev)
  d_at_bar <- -2 * loglik_fun(beta_bar, X_use, risk_sets, event_sets)
  p_d <- d_bar - d_at_bar
  dic <- d_bar + p_d

  list(
    Dbar = d_bar,
    Dhat = d_at_bar,
    pD = p_d,
    DIC = dic,
    deviance = dev
  )
}
