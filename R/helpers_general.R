
#' Standardize a Design Matrix
#'
#' Standardizes the non-intercept columns of a design matrix.
#'
#' If an intercept is present, covariates are centered and scaled.
#' If no intercept is present, covariates are scaled without centering
#' so that the no-intercept model specification is preserved.
#'
#' @param X A numeric design matrix.
#' @param standardize Logical; whether to standardize the covariates.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{X}: The transformed design matrix.
#'   \item \code{standardized}: Whether standardization was performed.
#'   \item \code{has_intercept}: Whether an intercept is present.
#'   \item \code{intercept_idx}: Column index of the intercept, if present.
#'   \item \code{coef_idx}: Column indices of the non-intercept covariates.
#'   \item \code{x_center}: Centering values used for the covariates.
#'   \item \code{x_scale}: Scaling values used for the covariates.
#' }
#'
#' @details
#' When an intercept is included, each covariate is transformed as
#' \deqn{
#' x_j^\ast = \frac{x_j - \bar{x}_j}{s_j}.
#' }
#' When no intercept is included, centering would introduce a constant
#' shift in the linear predictor. Therefore, only scaling is performed:
#' \deqn{
#' x_j^\ast = \frac{x_j}{s_j}.
#' }
#'
#' @keywords internal
standardize_design <- function(X, standardize = FALSE) {

  has_intercept <- "(Intercept)" %in% colnames(X)

  intercept_idx <- if (has_intercept) {
    which(colnames(X) == "(Intercept)")
  } else {
    integer(0)
  }

  coef_idx <- which(colnames(X) != "(Intercept)")

  # Standardization is actually performed only when requested
  # and at least one non-intercept covariate exists.
  standardized <- isTRUE(standardize) && length(coef_idx) > 0L

  out <- list(
    X = X,
    standardized = standardized,
    has_intercept = has_intercept,
    intercept_idx = intercept_idx,
    coef_idx = coef_idx,
    x_center = NULL,
    x_scale = NULL
  )

  if (!standardized) {
    return(out)
  }

  Xc <- X

  # Scaling factors
  x_scale <- apply(
    Xc[, coef_idx, drop = FALSE],
    2,
    stats::sd
  )

  # Avoid division by zero or non-finite scaling factors
  x_scale[!is.finite(x_scale) | x_scale == 0] <- 1

  if (has_intercept) {

    # Center and scale when an intercept is present
    x_center <- colMeans(
      Xc[, coef_idx, drop = FALSE]
    )

    Xc[, coef_idx] <- sweep(
      Xc[, coef_idx, drop = FALSE],
      2,
      x_center,
      "-"
    )

    Xc[, coef_idx] <- sweep(
      Xc[, coef_idx, drop = FALSE],
      2,
      x_scale,
      "/"
    )

  } else {

    # Scale only when there is no intercept
    x_center <- rep(0, length(coef_idx))
    names(x_center) <- colnames(Xc)[coef_idx]

    Xc[, coef_idx] <- sweep(
      Xc[, coef_idx, drop = FALSE],
      2,
      x_scale,
      "/"
    )
  }

  out$X <- Xc
  out$x_center <- x_center
  out$x_scale <- x_scale

  out
}


#' Back-transform Regression Coefficient Draws
#'
#' Transforms posterior draws of regression coefficients from the
#' standardized covariate scale back to the original covariate scale.
#'
#' @param beta_draws A numeric matrix of posterior draws for regression
#'   coefficients, with parameters in columns.
#' @param standardization A list returned by
#'   \code{standardize_design()}.
#'
#' @return A numeric matrix of regression coefficient draws on the
#'   original covariate scale.
#'
#' @details
#' For non-intercept coefficients,
#' \deqn{
#' \beta_j = \frac{\beta_j^\ast}{s_j}.
#' }
#' If an intercept is present, it is additionally transformed as
#' \deqn{
#' \alpha
#' =
#' \alpha^\ast
#' -
#' \sum_j \bar{x}_j \beta_j.
#' }
#' If no intercept is present, only the scaling transformation is
#' required because the covariates were not centered.
#'
#' @keywords internal
backtransform_beta <- function(beta_draws, standardization) {

  if (!standardization$standardized) {
    return(beta_draws)
  }

  coef_idx <- standardization$coef_idx
  x_center <- standardization$x_center
  x_scale <- standardization$x_scale

  # Back-transform non-intercept coefficients
  beta_draws[, coef_idx] <- sweep(
    beta_draws[, coef_idx, drop = FALSE],
    2,
    x_scale,
    "/"
  )

  # Adjust intercept when present
  if (standardization$has_intercept) {

    intercept_idx <- standardization$intercept_idx

    intercept_adjustment <- rowSums(
      sweep(
        beta_draws[, coef_idx, drop = FALSE],
        2,
        x_center,
        "*"
      )
    )

    beta_draws[, intercept_idx] <-
      beta_draws[, intercept_idx] - intercept_adjustment
  }

  beta_draws
}


#' Transform Initial Regression Coefficients
#'
#' Transforms initial regression coefficients from the original
#' covariate scale to the scale used internally by the MCMC sampler.
#'
#' @param beta Numeric vector of initial regression coefficients on the
#'   original covariate scale.
#' @param standardization A list returned by
#'   \code{standardize_design()}.
#'
#' @return A numeric vector of regression coefficients on the internal
#'   covariate scale.
#'
#' @details
#' For non-intercept coefficients,
#' \deqn{
#' \beta_j^\ast = s_j \beta_j.
#' }
#' If an intercept is present, it is additionally transformed as
#' \deqn{
#' \alpha^\ast
#' =
#' \alpha
#' +
#' \sum_j \bar{x}_j \beta_j.
#' }
#' If no intercept is present, only the scaling transformation is
#' required because the covariates are not centered.
#'
#' @keywords internal
transform_initial_beta <- function(beta, standardization) {

  beta <- as.numeric(beta)

  if (!standardization$standardized) {
    return(beta)
  }

  coef_idx <- standardization$coef_idx
  x_center <- standardization$x_center
  x_scale <- standardization$x_scale

  beta_original <- beta

  # Transform intercept before scaling slope coefficients
  if (standardization$has_intercept) {

    intercept_idx <- standardization$intercept_idx

    beta[intercept_idx] <-
      beta_original[intercept_idx] +
      sum(
        x_center * beta_original[coef_idx]
      )
  }

  # Transform non-intercept coefficients
  beta[coef_idx] <-
    beta_original[coef_idx] * x_scale

  beta
}


#' Process a Credible Interval Level
#'
#' Validates a credible interval level and returns the corresponding
#' lower and upper tail probabilities and display labels.
#'
#' @param level Numeric scalar specifying the credible interval level.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{level}: The credible interval level.
#'   \item \code{probs}: Lower and upper quantile probabilities.
#'   \item \code{lower_name}: Display name for the lower interval bound.
#'   \item \code{upper_name}: Display name for the upper interval bound.
#' }
#'
#' @keywords internal
credible_interval_info <- function(level = 0.95) {

  if (!is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop("`level` must be a single numeric value between 0 and 1.")
  }

  alpha <- (1 - level) / 2
  probs <- c(alpha, 1 - alpha)

  level_percent <- format(
    100 * level,
    trim = TRUE,
    scientific = FALSE
  )

  lower_name <- paste0("Lower ", level_percent, "%")
  upper_name <- paste0("Upper ", level_percent, "%")

  list(
    level = level,
    probs = probs,
    lower_name = lower_name,
    upper_name = upper_name
  )
}

#' Validate Initial Values
#'
#' Validates user-supplied initial values for PL-Cox and GPL-Cox models.
#'
#' @param init A named list of initial values, or \code{NULL}.
#' @param p Number of regression coefficients.
#' @param has_frailty Logical; whether the model contains frailty.
#' @param n_groups Number of frailty groups.
#'
#' @return The validated initialization list.
#'
#' @keywords internal
validate_init <- function(init,
                          p,
                          has_frailty = FALSE,
                          n_groups = NULL) {

  if (is.null(init)) {
    return(NULL)
  }

  if (!is.list(init) || is.null(names(init))) {
    stop("`init` must be NULL or a named list.")
  }

  allowed <- if (has_frailty) {
    c("beta", "frailty", "frailty_var")
  } else {
    "beta"
  }

  unknown <- setdiff(names(init), allowed)

  if (length(unknown) > 0L) {
    stop(
      "Unknown element(s) in `init`: ",
      paste(unknown, collapse = ", "),
      "."
    )
  }

  # beta
  if (!is.null(init$beta)) {

    if (!is.numeric(init$beta) ||
        length(init$beta) != p ||
        any(!is.finite(init$beta))) {
      stop(
        "`init$beta` must be a finite numeric vector of length ",
        p,
        "."
      )
    }
  }

  if (has_frailty) {

    # frailty
    if (!is.null(init$frailty)) {

      if (!is.numeric(init$frailty) ||
          length(init$frailty) != n_groups ||
          any(!is.finite(init$frailty))) {
        stop(
          "`init$frailty` must be a finite numeric vector of length ",
          n_groups,
          "."
        )
      }
    }

    # frailty variance
    if (!is.null(init$frailty_var)) {

      if (!is.numeric(init$frailty_var) ||
          length(init$frailty_var) != 1L ||
          !is.finite(init$frailty_var) ||
          init$frailty_var <= 0) {
        stop(
          "`init$frailty_var` must be a single positive finite value."
        )
      }
    }
  }

  init
}


#' Validate MCMC Execution Settings
#'
#' Validates settings used for multiple MCMC chains.
#'
#' @param chains Number of MCMC chains.
#' @param parallel Logical; whether chains should be run in parallel.
#' @param cores Number of CPU cores requested for parallel execution.
#' @param seed Optional master random seed.
#'
#' @return A list containing validated execution settings.
#'
#' @keywords internal
validate_mcmc_execution <- function(chains = 1L,
                                    parallel = FALSE,
                                    cores = 1L,
                                    seed = NULL) {

  if (!is.numeric(chains) ||
      length(chains) != 1L ||
      !is.finite(chains) ||
      chains < 1 ||
      chains != floor(chains)) {
    stop("`chains` must be a positive integer.")
  }

  if (!is.logical(parallel) ||
      length(parallel) != 1L ||
      is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.")
  }

  if (!is.numeric(cores) ||
      length(cores) != 1L ||
      !is.finite(cores) ||
      cores < 1 ||
      cores != floor(cores)) {
    stop("`cores` must be a positive integer.")
  }

  if (!is.null(seed)) {

    if (!is.numeric(seed) ||
        length(seed) != 1L ||
        !is.finite(seed) ||
        seed < 0 ||
        seed > .Machine$integer.max ||
        seed != floor(seed)) {
      stop(
        "`seed` must be NULL or a single non-negative integer ",
        "not exceeding .Machine$integer.max."
      )
    }

    seed <- as.integer(seed)
  }

  chains <- as.integer(chains)
  cores <- as.integer(cores)

  list(
    chains = chains,
    parallel = parallel,
    cores = cores,
    seed = seed
  )
}

#' Prepare Initial Values for Multiple Chains
#'
#' Expands user-supplied initial values into one initialization
#' specification for each MCMC chain.
#'
#' @param init Initial values supplied by the user.
#' @param chains Number of MCMC chains.
#'
#' @return A list of length \code{chains}.
#'
#' @keywords internal
prepare_chain_inits <- function(init, chains) {

  # No user-supplied initialization
  if (is.null(init)) {
    return(rep(list(NULL), chains))
  }

  init_names <- c(
    "beta",
    "frailty",
    "frailty_var"
  )

  # One named initialization specification:
  # use the same initial values for every chain.
  if (is.list(init) &&
      !is.null(names(init)) &&
      any(names(init) %in% init_names)) {

    return(rep(list(init), chains))
  }

  # Separate initialization specification for each chain.
  if (is.list(init) &&
      length(init) == chains &&
      all(
        vapply(
          init,
          function(x) {
            is.null(x) ||
              (is.list(x) && !is.null(names(x)))
          },
          logical(1)
        )
      )) {

    return(init)
  }

  stop(
    "`init` must be NULL, a named initialization list, ",
    "or a list of initialization lists of length `chains`."
  )
}

#' Generate RNG States for Multiple MCMC Chains
#'
#' Generates independent L'Ecuyer-CMRG random-number streams for
#' multiple MCMC chains.
#'
#' @param chains Number of MCMC chains.
#' @param seed Optional master random seed.
#'
#' @return A list containing the master seed and one RNG state
#'   for each chain.
#'
#' @keywords internal
make_chain_rng_states <- function(chains, seed = NULL) {

  # If no seed is supplied, draw one from the current RNG state.
  if (is.null(seed)) {
    seed_used <- sample.int(
      .Machine$integer.max,
      size = 1L
    )
  } else {
    seed_used <- as.integer(seed)
  }

  # Save the current RNG state after generating seed_used.
  old_kind <- RNGkind()

  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (had_seed) {
    old_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }

  on.exit({
    do.call(
      RNGkind,
      as.list(old_kind)
    )
    if (had_seed) {
      assign(
        ".Random.seed",
        old_seed,
        envir = .GlobalEnv
      )
    } else if (exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )) {
      rm(
        ".Random.seed",
        envir = .GlobalEnv
      )
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed_used)

  states <- vector(
    "list",
    chains
  )

  states[[1L]] <- .Random.seed

  if (chains > 1L) {
    for (k in 2:chains) {
      states[[k]] <- parallel::nextRNGStream(
        states[[k - 1L]]
      )
    }
  }

  list(
    seed = seed_used,
    kind = "L'Ecuyer-CMRG",
    states = states
  )
}

#' Run a Single MCMC Chain
#'
#' Runs one MCMC chain using a predetermined random-number state.
#'
#' @param job A list containing \code{chain_id}, \code{init}, and
#'   \code{rng_state}.
#' @param fit_fun Single-chain fitting function.
#' @param fit_args Arguments passed to the fitting function.
#'
#' @return A fitted single-chain object.
#'
#' @keywords internal
run_one_mcmc_chain <- function(job,
                               fit_fun,
                               fit_args) {

  RNGkind("L'Ecuyer-CMRG")

  assign(
    ".Random.seed",
    job$rng_state,
    envir = .GlobalEnv
  )

  start_time <- proc.time()[["elapsed"]]

  fit <- do.call(
    fit_fun,
    c(
      fit_args,
      list(init = job$init)
    )
  )

  elapsed <-
    proc.time()[["elapsed"]] - start_time

  fit$chain_id <- job$chain_id
  fit$chain_elapsed <- elapsed

  fit
}

#' Run Multiple MCMC Chains
#'
#' Runs independent MCMC chains sequentially or in parallel.
#'
#' @param fit_fun_name Single-chain fitting function.
#' @param fit_args Arguments passed to \code{fit_fun}.
#' @param chains Number of MCMC chains.
#' @param parallel Logical; whether chains should be run in parallel.
#' @param cores Number of CPU cores.
#' @param init Initial values.
#' @param seed Optional master random seed.
#'
#' @return A list containing fitted chain objects and execution metadata.
#'
#' @keywords internal
run_mcmc_chains <- function(fit_fun_name,
                            fit_args,
                            chains = 1L,
                            parallel = FALSE,
                            cores = 1L,
                            init = NULL,
                            seed = NULL) {

  ns <- asNamespace("BayesPLCox")

  fit_fun <- get(
    fit_fun_name,
    envir = ns,
    inherits = FALSE
  )

  run_fun <- get(
    "run_one_mcmc_chain",
    envir = ns,
    inherits = FALSE
  )

  control <- validate_mcmc_execution(
    chains = chains,
    parallel = parallel,
    cores = cores,
    seed = seed
  )

  chain_inits <- prepare_chain_inits(
    init = init,
    chains = control$chains
  )

  rng_info <- make_chain_rng_states(
    chains = control$chains,
    seed = control$seed
  )

  # Preserve the main-process RNG state.
  main_kind <- RNGkind()

  had_main_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (had_main_seed) {
    main_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }

  on.exit({
    do.call(
      RNGkind,
      as.list(main_kind)
    )
    if (had_main_seed) {
      assign(
        ".Random.seed",
        main_seed,
        envir = .GlobalEnv
      )
    } else if (exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )) {
      rm(
        ".Random.seed",
        envir = .GlobalEnv
      )
    }
  }, add = TRUE)

  chain_ids <- seq_len(
    control$chains
  )

  jobs <- lapply(
    chain_ids,
    function(k) {
      list(
        chain_id = k,
        init = chain_inits[[k]],
        rng_state = rng_info$states[[k]]
      )
    }
  )

  cores_used <- min(
    control$cores,
    control$chains
  )

  parallel_used <-
    control$parallel &&
    control$chains > 1L &&
    cores_used > 1L

  start_time <- proc.time()[["elapsed"]]

  if (!parallel_used) {
    fits <- lapply(
      jobs,
      run_fun,
      fit_fun = fit_fun,
      fit_args = fit_args
    )
    cores_used <- 1L
  } else {
    cl <- parallel::makePSOCKcluster(
      cores_used
    )

    on.exit(
      parallel::stopCluster(cl),
      add = TRUE
    )

    # Load the installed BayesPLCox namespace on each worker
    parallel::clusterCall(
      cl,
      function() {
        loadNamespace("BayesPLCox")
        NULL
      }
    )

    worker_run <- function(job,
                           fit_fun_name,
                           fit_args) {

      ns <- asNamespace("BayesPLCox")

      fit_fun <- get(
        fit_fun_name,
        envir = ns,
        inherits = FALSE
      )

      run_fun <- get(
        "run_one_mcmc_chain",
        envir = ns,
        inherits = FALSE
      )

      run_fun(
        job = job,
        fit_fun = fit_fun,
        fit_args = fit_args
      )
    }

    fits <- parallel::parLapply(
      cl,
      jobs,
      worker_run,
      fit_fun_name = fit_fun_name,
      fit_args = fit_args
    )
  }

  elapsed <- proc.time()[["elapsed"]] - start_time

  list(
    fits = fits,
    seed = rng_info$seed,
    rng_kind = rng_info$kind,
    rng_states = rng_info$states,
    chains = control$chains,
    parallel = parallel_used,
    cores = cores_used,
    elapsed = elapsed
  )
}

#' Combine Multiple MCMC Chain Fits
#'
#' Combines posterior draws and metadata from independent
#' single-chain fits.
#'
#' @param chain_run Output from \code{run_mcmc_chains()}.
#' @param has_frailty Logical; whether frailty terms are present.
#'
#' @return A combined fitted object.
#'
#' @keywords internal
combine_mcmc_chains <- function(chain_run,
                                has_frailty = FALSE) {

  chain_fits <- chain_run$fits
  chains <- length(chain_fits)

  res <- chain_fits[[1L]]

  # Regression coefficients
  res$beta <- do.call(
    rbind,
    lapply(chain_fits, function(x) x$beta)
  )

  # Number of saved draws
  n_save_per_chain <- vapply(
    chain_fits,
    function(x) nrow(x$beta),
    integer(1)
  )

  names(n_save_per_chain) <-
    paste0("chain", seq_len(chains))

  # Chain ID
  res$chain_id <- rep(
    seq_len(chains),
    times = n_save_per_chain
  )

  # Initial values
  res$initial_values <- lapply(
    chain_fits,
    function(x) x$initial_values
  )

  # Frailty components
  if (has_frailty) {

    res$u <- do.call(
      rbind,
      lapply(chain_fits, function(x) x$u)
    )

    res$sigma2 <- unlist(
      lapply(
        chain_fits,
        function(x) x$sigma2
      ),
      use.names = FALSE
    )
  }

  # Timing
  res$chain_elapsed <- vapply(
    chain_fits,
    function(x) x$chain_elapsed,
    numeric(1)
  )

  names(res$chain_elapsed) <-
    paste0("chain", seq_len(chains))

  res$elapsed <- chain_run$elapsed

  # Execution information
  res$n_chains <- chains
  res$seed <- chain_run$seed
  res$rng_kind <- chain_run$rng_kind
  res$parallel <- chain_run$parallel
  res$cores <- chain_run$cores

  # Posterior dimensions
  res$n_save <- nrow(res$beta)
  res$n_save_per_chain <- n_save_per_chain
  res$n_par <- ncol(res$beta)

  res
}

#' Convert Combined MCMC Draws to a Draws Array
#'
#' Converts posterior draws stored as a combined matrix together with
#' chain identifiers into a posterior draws_array object.
#'
#' @param draws A numeric matrix or vector of posterior draws.
#' @param chain_id Integer vector identifying the chain of each draw.
#' @param variable_names Optional character vector of variable names.
#'
#' @return A \code{posterior::draws_array} object.
#'
#' @keywords internal
as_chain_draws_array <- function(draws,
                                 chain_id,
                                 variable_names = NULL) {

  if (is.null(dim(draws))) {
    draws <- matrix(draws, ncol = 1L)
  }

  if (!is.matrix(draws)) {
    draws <- as.matrix(draws)
  }

  if (length(chain_id) != nrow(draws)) {
    stop(
      "`chain_id` must have the same length as the number of draws."
    )
  }

  chains <- sort(unique(chain_id))
  n_chains <- length(chains)

  n_per_chain <- vapply(
    chains,
    function(k) sum(chain_id == k),
    integer(1)
  )

  if (length(unique(n_per_chain)) != 1L) {
    stop(
      "All chains must contain the same number of saved draws."
    )
  }

  n_iter <- n_per_chain[1L]
  n_var <- ncol(draws)

  if (is.null(variable_names)) {

    variable_names <- colnames(draws)

    if (is.null(variable_names)) {
      variable_names <- paste0("V", seq_len(n_var))
    }
  }

  if (length(variable_names) != n_var) {
    stop(
      "`variable_names` must have one name for each variable."
    )
  }

  out <- array(
    NA_real_,
    dim = c(
      n_iter,
      n_chains,
      n_var
    ),
    dimnames = list(
      iteration = as.character(seq_len(n_iter)),
      chain = as.character(chains),
      variable = variable_names
    )
  )

  for (j in seq_along(chains)) {

    idx <- chain_id == chains[j]

    out[, j, ] <- draws[
      idx,
      ,
      drop = FALSE
    ]
  }

  posterior::as_draws_array(out)
}

#' Compute MCMC Diagnostics
#'
#' Computes convergence and sampling-efficiency diagnostics from
#' multiple MCMC chains.
#'
#' @param draws A numeric matrix or vector of posterior draws.
#' @param chain_id Integer vector identifying the chain of each draw.
#' @param variable_names Optional variable names.
#'
#' @return A data frame containing R-hat, bulk ESS, tail ESS,
#'   and MCSE for the posterior mean.
#'
#' @keywords internal
compute_mcmc_diagnostics <- function(draws,
                                     chain_id,
                                     variable_names = NULL) {

  draws_array <- as_chain_draws_array(
    draws = draws,
    chain_id = chain_id,
    variable_names = variable_names
  )

  diag <- posterior::summarise_draws(
    draws_array,
    Rhat = posterior::rhat,
    `Bulk ESS` = posterior::ess_bulk,
    `Tail ESS` = posterior::ess_tail,
    MCSE = posterior::mcse_mean
  )

  diag <- as.data.frame(
    diag,
    stringsAsFactors = FALSE
  )

  if (length(unique(chain_id)) < 2L) {
    diag$Rhat <- NA_real_
  }

  diag
}
