#' Fit the PL-Cox model without frailty
#'
#' Internal fitting routine for the PL-Cox model with fixed effects only.
#'
#' @param formula A formula with a `Surv()` response.
#' @param data A data frame.
#' @param delta Precision parameter for the negative-binomial approximation.
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
plcox_fixed_fit <- function(formula, data,
                            delta = 10.0,
                            prior_mean = NULL, prior_var = 100,
                            n_iter = 2000, burn = 1000, thin = 1,
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

  n <- nrow(X); p <- ncol(X)
  dat <- build_pl_from_surv(Y[,1], Y[,2])
  num_t <- length(dat$event_times)

  # Set prior distribution
  b0 <- if(is.null(prior_mean)) rep(0, p) else prior_mean
  V0_inv <- diag(1/prior_var, p)

  # Initial values
  beta <- rep(0, p)
  Z_t  <- rep(0.1, num_t) # Latent variables for exponential expansion
  omega <- rep(1, n)      # Polya-Gamma latent variables

  # Storage
  beta_draws <- matrix(0, nrow = floor((n_iter - burn) / thin), ncol = p)
  colnames(beta_draws) <- colnames(X)

  for (it in 1:n_iter) {
    # -- 1. Sample Z_t
    # Z_t | \beta \sim Gamma(|E_t|, \sum_{j \in R_t} \exp(x_j \beta))
    eta_all <- as.numeric(X %*% beta)
    exp_eta <- exp(eta_all)

    for (t in 1:num_t) {
      denom <- sum(exp_eta[dat$risk_sets[[t]]])
      Z_t[t] <- rgamma(1, shape = dat$e_counts[t], rate = denom)
    }

    # Cumulative latent variable \zeta_i = \sum_{t: i \in R_t} Z_t for each individual
    zeta <- sapply(dat$member_of, function(rs) sum(Z_t[rs]))

    # -- 2. Sample \omega_i (Hamura et al. 2025)
    # \omega_i | \beta, \zeta_i \sim PG(w_i + \delta, x_i \beta + \log(\zeta_i) - \log(\delta))
    psi <- eta_all + log(pmax(zeta, 1e-12)) - log(delta)
    omega <- BayesLogit::rpg(n, h = dat$d + delta, z = psi)

    # -- 3. Sample \beta
    kappa <- (dat$d - delta) / 2

    # Precision matrix A and mean vector b
    Omega <- Matrix::Diagonal(x = omega)
    A <- t(X) %*% Omega %*% X + V0_inv
    z_adj <- log(pmax(zeta, 1e-12)) - log(delta)
    b_vec <- t(X) %*% (kappa - omega * z_adj) + V0_inv %*% b0

    cholA <- chol(A)
    beta <- backsolve(cholA, rnorm(p) + backsolve(cholA, b_vec, transpose = TRUE))

    # Storage
    if (it > burn && (it - burn) %% thin == 0) {
      beta_draws[(it - burn) / thin, ] <- as.numeric(beta)
    }

    # --- Progress Monitoring ---
    if (verbose) {
      # Calculate percentage based on 10% increments
      if (round(10 * it / n_iter) == (10 * it / n_iter)) {
        prog <- round(100 * it / n_iter)
        message(paste0("MCMC ", prog, "% completed"))
      }
    }
  }

  if (standardize) {
    for (j in 2:ncol(X)) {
      beta_draws[, j] <- beta_draws[, j] / sd[j-1]
    }
    intercept_adjustment <- rowSums(sweep(beta_draws[, -1, drop=FALSE], 2, mu, "*"))
    beta_draws[, 1] <- beta_draws[, 1] - intercept_adjustment
  }

  return(list(beta         = beta_draws,
              u            = NULL,
              sigma2       = NULL,
              colnames     = colnames(X),
              group_levels = NULL
  ))
}

#' Fit the PL-Cox model with log-normal frailty
#'
#' Internal fitting routine for the PL-Cox model with group-specific frailty.
#'
#' @inheritParams plcox_fixed_fit
#' @param a_sig,b_sig Hyperparameters for the frailty variance prior.
#'
#' @return A list containing posterior samples and supporting information.
#'
#' @keywords internal
plcox_frailty_fit <- function(formula, data,
                              delta = 10.0,
                              prior_mean = NULL, prior_var = 100,
                              a_sig = 0.01, b_sig = 0.01,
                              n_iter = 2000, burn = 1000, thin = 1,
                              standardize = FALSE,
                              verbose = FALSE) {

  # --- 1. Parse Formula and Extract Data ---
  fb <- lme4::findbars(formula)
  if (length(fb) == 0) stop("No random effects found in formula.")

  # Get fixed effects part
  fixed_formula <- lme4::nobars(formula)
  mf <- model.frame(fixed_formula, data)
  Y <- model.response(mf)
  if (!inherits(Y, "Surv")) stop("The response variable must be a Surv(time, status) object.")
  X <- model.matrix(fixed_formula, data)

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

  # Extract Grouping Factor for Frailty
  group_var <- as.character(fb[[1]][[3]])
  group_factor <- as.factor(data[[group_var]])
  group_id <- as.integer(group_factor)
  n_groups <- max(group_id)
  group_list <- split(1:nrow(data), group_id) # Pre-index for efficiency

  # Build Risk Set Data
  n <- nrow(X); p <- ncol(X)
  dat <- build_pl_from_surv(Y[,1], Y[,2])
  num_t <- length(dat$event_times)

  # --- 2. Setup Priors and Initial Values ---
  b0 <- if(is.null(prior_mean)) rep(0, p) else prior_mean
  V0_inv <- as.matrix(Matrix::Diagonal(x = rep(1 / prior_var, p)))

  beta   <- rep(0, p)
  u      <- rep(0, n_groups)
  sigma2 <- 1.0
  Z_t    <- rep(1, num_t) # Latent Gamma variables

  # Storage for MCMC samples
  n_save <- floor((n_iter - burn) / thin)
  beta_draws   <- matrix(0, nrow = n_save, ncol = p)
  u_draws      <- matrix(0, nrow = n_save, ncol = n_groups)
  sigma2_draws <- numeric(n_save)

  # --- 3. Gibbs Sampler ---
  for (it in 1:n_iter) {

    # -- Step A: Sample Z_t (Gamma)
    # eta includes both fixed and random effects
    eta_all <- as.numeric(X %*% beta) + u[group_id]
    exp_eta <- exp(eta_all)

    for (t in seq_along(dat$event_times)) {
      denom <- sum(exp_eta[dat$risk_sets[[t]]])
      Z_t[t] <- rgamma(1, shape = dat$e_counts[t], rate = denom)
    }

    # Cumulative latent variable zeta_i = sum of Z_t where individual i is at risk
    zeta <- sapply(dat$member_of, function(rs) sum(Z_t[rs]))

    # -- Step B: Sample Polya-Gamma latent variables (omega_i)
    # psi_i = eta_i + log(zeta_i) - log(delta)
    psi   <- eta_all + log(pmax(zeta, 1e-12)) - log(delta)
    omega <- BayesLogit::rpg(n, h = dat$d + delta, z = psi)

    # # -- Step C: Sample Random Effects (u_j)
    # 1. Construct the ingredients
    # Upper-left: fixed effects precision
    A_fixed <- as.matrix(t(X) %*% (omega * X) + V0_inv)

    # Lower-right: random effects precision (Diagonal matrix)
    # Each diagonal element is the sum of omega in that group + 1/sigma2
    A_rand_diag <- sapply(group_list, function(idx) sum(omega[idx])) + (1 / sigma2)

    # Off-diagonal: cross-precision between fixed and random effects
    # Each column j is sum(omega_i * x_i) for individuals in group j
    A_off <- matrix(0, nrow = p, ncol = n_groups)
    for(j in 1:n_groups) {
      idx_j <- group_list[[j]]
      if(length(idx_j) > 0) {
        A_off[, j] <- colSums(omega[idx_j] * X[idx_j, , drop=FALSE])
      }
    }

    # 2. Assemble the Joint Precision Matrix (Full size: p + n_groups)
    # Using Matrix package for sparse efficiency if n_groups is large
    A_joint <- Matrix::rbind2(
      Matrix::cbind2(A_fixed, A_off),
      Matrix::cbind2(t(A_off), Matrix::Diagonal(x = A_rand_diag))
    )

    # 3. Construct the Joint RHS Vector b
    kappa_adj <- (dat$d - delta)/2 - omega * (log(pmax(zeta, 1e-12)) - log(delta))
    b_fixed <- as.numeric(t(X) %*% kappa_adj + V0_inv %*% b0)
    b_rand  <- sapply(group_list, function(idx) sum(kappa_adj[idx]))
    b_joint <- c(b_fixed, b_rand)

    # 4. Jointly Sample (beta, u)
    # Using Cholesky decomposition on the joint matrix
    cholA_j <- Matrix::chol(A_joint)
    # Sample from Multivariate Normal: A_joint^{-1} * b_joint + cholA_j^{-1} * epsilon
    alpha_joint <- Matrix::solve(cholA_j, rnorm(p + n_groups) +
                                   Matrix::solve(Matrix::t(cholA_j), b_joint, system = "L"),
                                 system = "Lt")

    # 5. Split back into beta and u
    beta <- as.numeric(alpha_joint[1:p])
    u    <- as.numeric(alpha_joint[(p+1):(p + n_groups)])


    # -- Step E: Sample Frailty Variance (sigma2)
    shape_post <- a_sig + (n_groups / 2)
    rate_post  <- b_sig + (sum(u^2) / 2)
    sigma2     <- 1 / rgamma(1, shape = shape_post, rate = rate_post)

    # -- Step F: Store draws after burn-in and thinning
    if (it > burn && (it - burn) %% thin == 0) {
      idx_save <- (it - burn) / thin
      beta_draws[idx_save, ]   <- beta
      u_draws[idx_save, ]      <- u
      sigma2_draws[idx_save]   <- sigma2
    }

    # --- Progress Monitoring ---
    if (verbose) {
      # Calculate percentage based on 10% increments
      if (round(10 * it / n_iter) == (10 * it / n_iter)) {
        prog <- round(100 * it / n_iter)
        message(paste0("MCMC ", prog, "% completed"))
      }
    }
  }

  if (standardize) {
    for (j in 2:ncol(X)) {
      beta_draws[, j] <- beta_draws[, j] / sd[j-1]
    }
    intercept_adjustment <- rowSums(sweep(beta_draws[, -1, drop=FALSE], 2, mu, "*"))
    beta_draws[, 1] <- beta_draws[, 1] - intercept_adjustment
  }

  # --- 4. Return results ---
  return(list(
    beta         = beta_draws,
    u            = u_draws,
    sigma2       = sigma2_draws,
    colnames     = colnames(X),
    group_levels = levels(group_factor)
  ))
}

#' Fit a Bayesian Cox model via the PL-Cox likelihood
#'
#' Fits a Bayesian Cox regression model based on the PL-Cox likelihood using
#' Gibbs sampling. The acronym "PL" reflects both the Cox partial likelihood
#' and the Plackett-Luce model for rank-ordered data.
#'
#' @param formula A formula with a `Surv()` response. A frailty term of the form
#'   `(1 | group)` may be included.
#' @param data A data frame.
#' @param delta Precision parameter for the negative-binomial approximation.
#' @param prior_mean Prior mean vector for regression coefficients.
#' @param prior_var Prior variance for regression coefficients.
#' @param a_sig,b_sig Hyperparameters for the frailty variance prior.
#' @param n_iter Total number of MCMC iterations.
#' @param burn Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param standardize Logical; standardize covariates except the intercept.
#' @param verbose Logical; print progress messages?
#'
#' @return An object of class `"plcox"`.
#'
#' @examples
#' \dontrun{
#' fit <- plcox(Surv(time, status) ~ x1 + x2, data = dat)
#' summary(fit)
#' plot(fit)
#' traceplot(fit)
#' }
#'
#' @export
plcox <- function(formula, data,
                  delta = 10.0,
                  prior_mean = NULL, prior_var = 100,
                  a_sig = 0.01, b_sig = 0.01,
                  n_iter = 2000, burn = 1000, thin = 1,
                  standardize = FALSE,
                  verbose = FALSE) {
  has_frailty <- any(grepl("\\|", format(formula)))

  if (has_frailty) {
    res <- plcox_frailty_fit(formula, data,
                             delta, prior_mean, prior_var, a_sig, b_sig,
                             n_iter, burn, thin, standardize, verbose)
  } else {
    res <- plcox_fixed_fit(formula, data,
                           delta, prior_mean, prior_var,
                           n_iter, burn, thin, standardize, verbose)
  }

  res$method <- "PL-Cox"
  res$has_frailty <- has_frailty
  res$call <- match.call()
  res$formula <- formula
  res$n_iter <- n_iter
  res$burn <- burn
  res$thin <- thin
  res$standardize <- standardize
  res$n_save <- nrow(res$beta)
  res$n_par  <- ncol(res$beta)

  class(res) <- c("plcox", "bayescoxrank")
  return(res)
}

