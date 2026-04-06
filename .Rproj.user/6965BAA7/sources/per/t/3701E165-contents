#' Fit the GPL-Cox model without frailty
#'
#' Internal fitting routine for the GPL-Cox model with fixed effects only.
#'
#' @param formula A formula with a `Surv()` response.
#' @param data A data frame.
#' @param prior_mean Prior mean vector for regression coefficients.
#' @param prior_var Prior variance for regression coefficients.
#' @param n_iter Total number of MCMC iterations.
#' @param burn Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param standardize Logical; standardize covariates except the intercept.
#' @param verbose Logical; print progress messages?
#'
#' @return A list containing posterior samples and supporting information.
#'
#' @keywords internal
gplcox_fixed_fit <- function(formula, data,
                             prior_mean = NULL, prior_var = 100,
                             n_iter = 4000, burn = 2000, thin = 1,
                             standardize = FALSE,
                             verbose = FALSE) {

  mf <- model.frame(formula, data = data)
  Y <- model.response(mf)
  if (!inherits(Y, "Surv")) stop("The response variable must be a Surv(time, status) object.")
  X <- model.matrix(formula, data = mf)

  # --- Standardization (excluding intercept) ---
  if (standardize) {
    Xc <- X
    if (colnames(Xc)[1] %in% c("(Intercept)", "(Intercept)")) {
      j0 <- 2:ncol(Xc)
    } else {
      j0 <- 1:ncol(Xc)
    }
    mu <- colMeans(Xc[, j0, drop=FALSE])
    sd <- apply(Xc[, j0, drop=FALSE], 2, sd)
    sd[sd == 0] <- 1
    Xc[, j0] <- scale(Xc[, j0, drop=FALSE], center = mu, scale = sd)
    X <- Xc
  }

  dat   <- build_gpl_from_surv(Y[,1], Y[,2], X)
  X     <- dat$X
  n     <- nrow(X)
  p     <- ncol(X)
  Rsets <- dat$R_sets
  d     <- dat$d

  if (is.null(prior_mean)) b0 <- rep(0, p)
  V0_inv <- as.matrix(Matrix::Diagonal(x = rep(1 / prior_var, p)))

  beta  <- rep(0, p)
  eta   <- as.numeric(X %*% beta)
  theta <- plogis(eta)

  # Initialize Z, zeta, and kappa
  Z     <- sample_Z_geom(theta, Rsets)
  zeta  <- accumulate_zeta_i(Z, Rsets, n)
  kappa <- d - 0.5 * zeta

  omega <- numeric(n)
  idx   <- zeta > 0
  if (any(idx)) {
    omega[idx] <- BayesLogit::rpg(sum(idx), h = as.numeric(zeta[idx]), z = as.numeric(eta[idx]))
    # --- Robustness: Set NA/Inf from PG to 0 ---
    bad <- !is.finite(omega[idx])
    if (any(bad)) omega[idx][bad] <- 0
  }

  n_save <- floor((n_iter - burn)/thin)
  beta_draws <- matrix(NA_real_, n_save, p)
  colnames(beta_draws) <- colnames(X)

  for (it in 1:n_iter) {
    # -- Z | theta
    eta   <- as.numeric(X %*% beta)
    theta <- plogis(eta)
    Z     <- sample_Z_geom(theta, Rsets)
    zeta  <- accumulate_zeta_i(Z, Rsets, n)

    # Guard against NA/Inf in kappa
    kappa <- d - 0.5 * zeta
    kappa[!is.finite(kappa)] <- 0

    # -- omega | beta, Z (Polya-Gamma; skip if h=0)
    omega[] <- 0
    idx <- zeta > 0
    if (any(idx)) {
      pg <- BayesLogit::rpg(sum(idx), h = as.numeric(zeta[idx]), z = as.numeric(eta[idx]))
      # Robustness: handle non-finite values
      pg[!is.finite(pg)] <- 0
      omega[idx] <- pg
    }

    # -- beta | omega, Z (Gaussian: Symmetrization)
    Omega   <- Matrix::Diagonal(x = omega)
    XtOmega <- t(X) %*% Omega
    A       <- XtOmega %*% X + V0_inv
    A <- as.matrix(A)
    # Symmetrize to remove numerical asymmetry
    A <- 0.5 * (A + t(A))

    b <- t(X) %*% kappa + V0_inv %*% b0

    cholA <- chol(A)
    V_post <- chol2inv(cholA)
    m_post <- V_post %*% b

    beta   <- as.numeric(MASS::mvrnorm(1, mu = as.numeric(m_post), Sigma = V_post))

    if (verbose && (round(10*it/n_iter)==(10*it/n_iter))) {
      message(sprintf("MCMC %3d%% completed", round(100 * it / n_iter)))
    }
    if (it > burn && ((it - burn) %% thin == 0)) {
      beta_draws[(it - burn) / thin, ] <- beta
    }
  }

  if (standardize) {
    for (j in 2:ncol(X)) {
      beta_draws[, j] <- beta_draws[, j] / sd[j-1]
    }
    intercept_adjustment <- rowSums(sweep(beta_draws[, -1, drop=FALSE], 2, mu, "*"))
    beta_draws[, 1] <- beta_draws[, 1] - intercept_adjustment
  }

  return(list(
    beta         = beta_draws,
    u            = NULL,
    sigma2       = NULL,
    colnames     = colnames(X),
    group_levels = NULL,
    time         = dat$time,
    status       = dat$status,
    X            = dat$X,
    ord          = dat$ord,
    R_sets       = dat$R_sets,
    E_sets       = dat$E_sets,
    d_obs        = sapply(dat$E_sets, length),
    group_id_ord = NULL
  ))
}

#' Fit the GPL-Cox model with frailty
#'
#' Internal fitting routine for the GPL-Cox model with group-specific frailty.
#'
#' @inheritParams gplcox_fixed_fit
#' @param a_sig,b_sig Hyperparameters for the frailty variance prior.
#'
#' @return A list containing posterior samples and supporting information.
#'
#' @keywords internal
gplcox_frailty_fit <- function(formula, data,
                               prior_mean = NULL, prior_var = 100,
                               a_sig = 0.01, b_sig = 0.01,
                               n_iter = 4000, burn = 2000, thin = 1,
                               standardize = FALSE,
                               verbose = FALSE) {

  # --- 1. Parse Formula and Extract Data ---
  fb <- lme4::findbars(formula)
  if (length(fb) == 0) stop("No random effects found in formula.")

  fixed_formula <- lme4::nobars(formula)
  mf <- model.frame(fixed_formula, data)
  Y  <- model.response(mf)
  if (!inherits(Y, "Surv")) stop("The response variable must be a Surv(time, status) object.")
  X  <- model.matrix(fixed_formula, data)

  # --- Standardization (excluding intercept) ---
  if (standardize) {
    Xc <- X
    if (colnames(Xc)[1] %in% c("(Intercept)", "(Intercept)")) {
      j0 <- 2:ncol(Xc)
    } else {
      j0 <- 1:ncol(Xc)
    }
    mu <- colMeans(Xc[, j0, drop=FALSE])
    sd <- apply(Xc[, j0, drop=FALSE], 2, sd)
    sd[sd == 0] <- 1
    Xc[, j0] <- scale(Xc[, j0, drop=FALSE], center = mu, scale = sd)
    X <- Xc
  }

  # Extract Grouping Factor
  group_var <- as.character(fb[[1]][[3]])
  group_factor <- as.factor(data[[group_var]])
  group_id <- as.integer(group_factor)
  n_groups <- max(group_id)
  group_list <- split(1:nrow(data), group_id)

  # Build GPL Risk Set Data
  dat <- build_gpl_from_surv(Y[,1], Y[,2], X)
  X_ord <- dat$X
  n <- nrow(X_ord); p <- ncol(X_ord)
  Rsets <- dat$R_sets
  d     <- dat$d
  group_id_ord <- group_id[dat$ord]
  group_list_ord <- split(seq_len(length(group_id_ord)), group_id_ord)

  # --- 2. Setup Priors and Initial Values ---
  b0 <- if(is.null(prior_mean)) rep(0, p) else prior_mean
  V0_inv <- as.matrix(Matrix::Diagonal(x = rep(1 / prior_var, p)))

  beta   <- rep(0, p)
  u      <- rep(0, n_groups)
  sigma2 <- 1.0
  omega  <- numeric(n)

  # Storage
  n_save <- floor((n_iter - burn) / thin)
  beta_draws   <- matrix(0, n_save, p)
  u_draws      <- matrix(0, n_save, n_groups)
  sigma2_draws <- numeric(n_save)

  # --- 3. Gibbs Sampler ---
  for (it in 1:n_iter) {

    # -- Step A: Update theta and sample Z (Geometric)
    # eta includes frailty: eta = Xb + u_g
    eta <- as.numeric(X_ord %*% beta) + u[group_id_ord]
    theta <- plogis(eta)
    Z     <- sample_Z_geom(theta, Rsets)
    zeta  <- accumulate_zeta_i(Z, Rsets, n)
    kappa <- d - 0.5 * zeta

    # -- Step B: Update Polya-Gamma omega_i
    omega <- numeric(n)
    idx <- zeta > 0
    if (any(idx)) {
      pg <- BayesLogit::rpg(sum(idx), h = as.numeric(zeta[idx]), z = as.numeric(eta[idx]))
      pg[!is.finite(pg)] <- 0
      omega[idx] <- pg
    }

    # -- B: Construct joint precision matrix for (beta, u) --
    A_fixed <- as.matrix(t(X_ord) %*% (omega * X_ord) + V0_inv)

    # Lower-right: random effects precision (Diagonal matrix)
    A_rand_diag <- sapply(group_list_ord, function(idx) sum(omega[idx])) + (1 / sigma2)

    # Off-diagonal: cross-precision between fixed effects and frailties
    A_off <- matrix(0, nrow = p, ncol = n_groups)
    for (j in seq_len(n_groups)) {
      idx_j <- group_list_ord[[j]]
      if (length(idx_j) > 0) {
        A_off[, j] <- colSums(omega[idx_j] * X_ord[idx_j, , drop = FALSE])
      }
    }

    A_joint <- Matrix::rbind2(
      Matrix::cbind2(A_fixed, A_off),
      Matrix::cbind2(t(A_off), Matrix::Diagonal(x = A_rand_diag))
    )

    # -- C: Construct Joint RHS Vector b --
    b_fixed <- as.numeric(t(X_ord) %*% kappa + V0_inv %*% b0)
    b_rand  <- sapply(group_list_ord, function(idx) sum(kappa[idx]))
    b_joint <- c(b_fixed, b_rand)

    # -- D: Joint Sampling --
    cholA <- Matrix::chol(A_joint)
    alpha <- Matrix::solve(cholA, rnorm(p + n_groups) + Matrix::solve(Matrix::t(cholA), b_joint, system = "L"), system = "Lt")

    beta <- as.numeric(alpha[1:p])
    u    <- as.numeric(alpha[(p+1):(p + n_groups)])

    # -- Step E: Update sigma2 (Inverse-Gamma)
    shape_post <- a_sig + (n_groups / 2)
    rate_post  <- b_sig + (sum(u^2) / 2)
    sigma2     <- 1 / rgamma(1, shape = shape_post, rate = rate_post)

    # --- Progress and Storage ---
    if (verbose && (round(10*it/n_iter) == (10*it/n_iter))) {
      message(sprintf("MCMC %3d%% completed", round(100 * it / n_iter)))
    }

    if (it > burn && ((it - burn) %% thin == 0)) {
      idx_save <- (it - burn) / thin
      beta_draws[idx_save, ]   <- beta
      u_draws[idx_save, ]      <- u
      sigma2_draws[idx_save]   <- sigma2
    }
  }

  if (standardize) {
    for (j in 2:ncol(X)) {
      beta_draws[, j] <- beta_draws[, j] / sd[j-1]
    }
    intercept_adjustment <- rowSums(sweep(beta_draws[, -1, drop=FALSE], 2, mu, "*"))
    beta_draws[, 1] <- beta_draws[, 1] - intercept_adjustment
  }

  return(list(
    beta         = beta_draws,
    u            = u_draws,
    sigma2       = sigma2_draws,
    colnames     = colnames(X_ord),
    group_levels = levels(group_factor),
    time         = dat$time,
    status       = dat$status,
    X            = dat$X,
    ord          = dat$ord,
    R_sets       = dat$R_sets,
    E_sets       = dat$E_sets,
    d_obs        = sapply(dat$E_sets, length),
    group_id_ord = group_id_ord
  ))
}

#' Fit a Bayesian Cox model via the GPL-Cox likelihood
#'
#' Fits a Bayesian Cox regression model based on the generalized Plackett-Luce
#' (GPL) likelihood using Gibbs sampling. The GPL-Cox model provides a generative
#' representation for tied event sets and can therefore be used for posterior
#' predictive checks of tie structures.
#'
#' @param formula A formula with a `Surv()` response. A frailty term of the form
#'   `(1 | group)` may be included.
#' @param data A data frame.
#' @param prior_mean Prior mean vector for regression coefficients.
#' @param prior_var Prior variance for regression coefficients.
#' @param a_sig,b_sig Hyperparameters for the frailty variance prior.
#' @param n_iter Total number of MCMC iterations.
#' @param burn Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param standardize Logical; standardize covariates except the intercept.
#' @param verbose Logical; print progress messages?
#'
#' @return An object of class `"gplcox"`.
#'
#' @examples
#' \dontrun{
#' fit <- gplcox(Surv(time, status) ~ x1 + x2, data = dat)
#' summary(fit)
#' plot(fit)
#' traceplot(fit)
#' pp <- pp_check_tie.gplcox(fit)
#' plot(pp)
#' }
#'
#' @export
gplcox <- function(formula, data,
                   n_iter = 4000, burn = 2000, thin = 1,
                   prior_mean = NULL, prior_var = 100,
                   a_sig = 0.01, b_sig = 0.01,
                   standardize = FALSE,
                   verbose = FALSE) {
  has_frailty <- any(grepl("\\|", format(formula)))

  if (has_frailty) {
    res <- gplcox_frailty_fit(formula, data,
                              prior_mean, prior_var, a_sig, b_sig,
                              n_iter, burn, thin, standardize,
                              verbose)
  } else {
    res <- gplcox_fixed_fit(formula, data,
                            prior_mean, prior_var,
                            n_iter, burn, thin, standardize,
                            verbose)
  }

  res$method <- "GPL-Cox"
  res$has_frailty <- has_frailty
  res$call <- match.call()
  res$formula <- formula
  res$n_iter <- n_iter
  res$burn <- burn
  res$thin <- thin
  res$standardize <- standardize
  res$n_save <- nrow(res$beta)
  res$n_par  <- ncol(res$beta)

  class(res) <- c("gplcox", "bayescoxrank")
  return(res)
}
