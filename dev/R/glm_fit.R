#' @rdname feglm
#' @export
feglm.fit <- function(y, X, fixef_df, family = "gaussian", vcov, offset, split,
                      fsplit, split.keep, split.drop, cluster, se, ssc, weights, subset, start = NULL,
                      etastart = NULL, mustart = NULL, fixef.rm = "perfect",
                      fixef.tol = 1e-6, fixef.iter = 10000, collin.tol = 1e-10,
                      glm.iter = 25, glm.tol = 1e-8, nthreads = getFixest_nthreads(),
                      lean = FALSE, warn = TRUE, notes = getFixest_notes(), mem.clean = FALSE,
                      verbose = 0, only.env = FALSE, only.coef = FALSE, env, ...) {
  dots <- list(...)

  lean_internal <- isTRUE(dots$lean_internal)
  means <- 1
  if (!missing(env)) {
    # This is an internal call from the function feglm
    # no need to further check the arguments
    # we extract them from the env

    if ((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
      stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
    }

    # main variables
    if (missing(y)) y <- get("lhs", env)
    if (missing(X)) X <- get("linear.mat", env)
    if (!missing(fixef_df) && is.null(fixef_df)) {
      assign("isFixef", FALSE, env)
    }

    if (missing(offset)) offset <- get("offset.value", env)
    if (missing(weights)) weights <- get("weights.value", env)

    # other params
    if (missing(fixef.tol)) fixef.tol <- get("fixef.tol", env)
    if (missing(fixef.iter)) fixef.iter <- get("fixef.iter", env)
    if (missing(collin.tol)) collin.tol <- get("collin.tol", env)
    if (missing(glm.iter)) glm.iter <- get("glm.iter", env)
    if (missing(glm.tol)) glm.tol <- get("glm.tol", env)
    if (missing(warn)) warn <- get("warn", env)
    if (missing(verbose)) verbose <- get("verbose", env)
    if (missing(only.coef)) only.coef <- get("only.coef", env)

    # starting point of the fixed-effects
    if (!is.null(dots$means)) means <- dots$means

    # init
    init.type <- get("init.type", env)
    starting_values <- get("starting_values", env)
    if (lean_internal) {
      # Call within here => either null model or fe only
      init.type <- "default"
      if (!is.null(etastart)) {
        init.type <- "eta"
        starting_values <- etastart
      }
    }
  } else {
    if (missing(weights)) weights <- NULL

    time_start <- proc.time()

    set_defaults("fixest_estimation")
    call_env <- new.env(parent = parent.frame())

    env <- try(fixest_env(
      y = y, X = X, fixef_df = fixef_df, family = family,
      nthreads = nthreads, lean = lean, offset = offset,
      weights = weights, subset = subset, split = split,
      fsplit = fsplit,
      split.keep = split.keep, split.drop = split.drop,
      vcov = vcov, cluster = cluster,
      se = se, ssc = ssc, linear.start = start,
      etastart = etastart, mustart = mustart, fixef.rm = fixef.rm,
      fixef.tol = fixef.tol, fixef.iter = fixef.iter,
      collin.tol = collin.tol, glm.iter = glm.iter,
      glm.tol = glm.tol, notes = notes, mem.clean = mem.clean,
      only.coef = only.coef,
      warn = warn, verbose = verbose, origin = "feglm.fit",
      mc_origin = match.call(), call_env = call_env, ...
    ), silent = TRUE)

    if ("try-error" %in% class(env)) {
      stop(format_error_msg(env, "feglm.fit"))
    }

    check_arg(only.env, "logical scalar")
    if (only.env) {
      return(env)
    }

    verbose <- get("verbose", env)
    if (verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep = "")

    # y/X
    y <- get("lhs", env)
    X <- get("linear.mat", env)

    # offset
    offset <- get("offset.value", env)

    # weights
    weights <- get("weights.value", env)

    # warn flag
    warn <- get("warn", env)

    # init
    init.type <- get("init.type", env)
    starting_values <- get("starting_values", env)
  }

  IN_MULTI <- get("IN_MULTI", env)
  is_multi_root <- get("is_multi_root", env)
  if (is_multi_root) {
    on.exit(release_multi_notes())
    assign("is_multi_root", FALSE, env)
  }

  # Split ----

  do_split <- get("do_split", env)
  if (do_split) {
    res <- multi_split(env, feglm.fit)

    return(res)
  }

  # Multi fixef ----

  do_multi_fixef <- get("do_multi_fixef", env)
  if (do_multi_fixef) {
    res <- multi_fixef(env, feglm.fit)

    return(res)
  }

  # Multi LHS and RHS ----

  do_multi_lhs <- get("do_multi_lhs", env)
  do_multi_rhs <- get("do_multi_rhs", env)
  if (do_multi_lhs || do_multi_rhs) {
    res <- multi_LHS_RHS(env, feglm.fit)

    return(res)
  }

  # Regular estimation -----

  # Setup:
  family <- get("family_funs", env)
  isFixef <- get("isFixef", env)
  nthreads <- get("nthreads", env)
  isWeight <- length(weights) > 1
  isOffset <- length(offset) > 1
  nobs <- length(y)
  onlyFixef <- length(X) == 1

  # the preformatted results
  res <- get("res", env)

  #
  # glm functions:
  #

  variance <- family$variance
  linkfun <- family$linkfun
  linkinv <- family$linkinv

  dev.resids <- family$dev.resids
  sum_dev.resids <- function(y, mu, eta, wt) sum(dev.resids(y, mu, wt))

  fun_mu.eta <- family$mu.eta
  mu.eta <- function(mu, eta) fun_mu.eta(eta)

  valideta <- family$valideta
  if (is.null(valideta)) valideta <- function(...) TRUE

  validmu <- family$validmu
  if (is.null(validmu)) validmu <- function(mu) TRUE

  family_equiv <- family$family_equiv

  # Optimizing poisson/logit
  if (family_equiv == "poisson") {
    linkfun <- function(mu) cpppar_log(mu, nthreads)
    linkinv <- function(eta) cpppar_poisson_linkinv(eta, nthreads)

    y_pos <- y[y > 0]
    qui_pos <- y > 0
    if (isWeight) {
      constant <- sum(weights[qui_pos] * y_pos * cpppar_log(y_pos, nthreads) - weights[qui_pos] * y_pos)
      sum_dev.resids <- function(y, mu, eta, wt) 2 * (constant - sum(wt[qui_pos] * y_pos * eta[qui_pos]) + sum(wt * mu))
    } else {
      constant <- sum(y_pos * cpppar_log(y_pos, nthreads) - y_pos)
      sum_dev.resids <- function(y, mu, eta, wt) 2 * (constant - sum(y_pos * eta[qui_pos]) + sum(mu))
    }

    mu.eta <- function(mu, eta) mu
    validmu <- function(mu) cpppar_poisson_validmu(mu, nthreads)
  } else if (family_equiv == "logit") {
    # To avoid annoying warnings
    family$initialize <- quasibinomial()$initialize

    linkfun <- function(mu) cpppar_logit_linkfun(mu, nthreads)
    linkinv <- function(eta) cpppar_logit_linkinv(eta, nthreads)
    if (isWeight) {
      sum_dev.resids <- function(y, mu, eta, wt) sum(cpppar_logit_devresids(y, mu, wt, nthreads))
    } else {
      sum_dev.resids <- function(y, mu, eta, wt) sum(cpppar_logit_devresids(y, mu, 1, nthreads))
    }
    sum_dev.resids <- sum_dev.resids

    mu.eta <- function(mu, eta) cpppar_logit_mueta(eta, nthreads)
  }

  #
  # Init
  #

  if (mem.clean) {
    gc()
  }

  if (init.type == "mu") {
    mu <- starting_values

    if (!valideta(mu)) {
      stop("In 'mustart' the values provided are not valid.")
    }

    eta <- linkfun(mu)
  } else if (init.type == "eta") {
    eta <- starting_values

    if (!valideta(eta)) {
      stop("In 'etastart' the values provided are not valid.")
    }

    mu <- linkinv(eta)
  } else if (init.type == "coef") {
    # If there are fixed-effects we MUST first compute the FE model with starting values as offset
    #   otherwise we are too far away from the solution and starting values may lead to divergence
    #   (hence step halving would be required)
    # This means that initializing with coefficients incurs large computational costs
    #   with fixed-effects

    start <- get("start", env)
    offset_fe <- offset + cpppar_xbeta(X, start, nthreads)

    if (isFixef) {
      mustart <- 0
      eval(family$initialize)
      eta <- linkfun(mustart)

      # just a rough estimate (=> high tol values) [no benefit in high precision]
      model_fe <- try(feglm.fit(X = 0, etastart = eta, offset = offset_fe, glm.tol = 1e-2, fixef.tol = 1e-2, env = env, lean_internal = TRUE))

      if ("try-error" %in% class(model_fe)) {
        stop("Estimation failed during initialization when getting the fixed-effects, maybe change the values of 'start'? \n", model_fe)
      }

      eta <- model_fe$linear.predictors
      mu <- model_fe$fitted.values
      devold <- model_fe$deviance
    } else {
      eta <- offset_fe
      mu <- linkinv(eta)
      devold <- sum_dev.resids(y, mu, eta, wt = weights)
    }

    wols_old <- list(fitted.values = eta - offset)
  } else {
    mustart <- 0
    eval(family$initialize)
    eta <- linkfun(mustart)
    mu <- linkinv(eta)

    # NOTA: FE only => ADDS LOTS OF COMPUTATIONAL COSTS without convergence benefit
  }

  if (mem.clean) {
    gc()
  }

  if (init.type != "coef") {
    # starting deviance with constant equal to 1e-5
    # this is important for getting in step halving early (when deviance goes awry right from the start)
    devold <- sum_dev.resids(y, rep(linkinv(1e-5), nobs), rep(1e-5, nobs), wt = weights)
    wols_old <- list(fitted.values = rep(1e-5, nobs))
  }

  if (!validmu(mu) || !valideta(eta)) {
    stop("Current starting values are not valid.")
  }

  assign("nb_sh", 0, env)
  on.exit(warn_step_halving(env, stack_multi = IN_MULTI))

  if ((init.type == "coef" && verbose >= 1) || verbose >= 4) {
    cat("Deviance at initializat.  = ", numberFormatNormal(devold), "\n", sep = "")
  }

  #
  # The main loop
  #

  wols_means <- 1
  conv <- FALSE
  warning_msg <- div_message <- ""
  for (iter in 1:glm.iter) {
    if (mem.clean) {
      gc()
    }

    mu.eta.val <- mu.eta(mu, eta)
    var_mu <- variance(mu)

    # controls
    any_pblm_mu <- cpp_any_na_null(var_mu)
    if (any_pblm_mu) {
      if (anyNA(var_mu)) {
        stop("NAs in V(mu), at iteration ", iter, ".")
      } else if (any(var_mu == 0)) {
        stop("0s in V(mu), at iteration ", iter, ".")
      }
    }

    if (anyNA(mu.eta.val)) {
      stop("NAs in d(mu)/d(eta), at iteration ", iter, ".")
    }

    if (isOffset) {
      z <- (eta - offset) + (y - mu) / mu.eta.val
    } else {
      z <- eta + (y - mu) / mu.eta.val
    }

    w <- as.vector(weights * mu.eta.val**2 / var_mu)

    is_0w <- w == 0
    any_0w <- any(is_0w)
    if (any_0w && all(is_0w)) {
      warning_msg <- paste0("No informative observation at iteration ", iter, ".")
      div_message <- "No informative observation."
      break
    }

    if (mem.clean && iter > 1) {
      rm(wols)
      gc()
    }

    wols <- feols(
      y = z, X = X, weights = w, means = wols_means,
      correct_0w = any_0w, env = env, fixef.tol = fixef.tol,
      fixef.iter = fixef.iter, collin.tol = collin.tol, nthreads = nthreads,
      mem.clean = mem.clean, warn = warn,
      verbose = FALSE, fromGLM = TRUE
    )

    if (isTRUE(wols$NA_model)) {
      return(wols)
    }

    # In theory OLS estimation is guaranteed to exist
    # yet, NA coef may happen with non-infinite very large values of z/w (e.g. values > 1e100)
    if (anyNA(wols$coefficients)) {
      if (iter == 1) {
        stop("Weighted-OLS returns NA coefficients at first iteration, step halving cannot be performed. Try other starting values?")
      }

      warning_msg <- paste0("Divergence at iteration ", iter, ": Weighted-OLS returns NA coefficients. Last evaluated coefficients with finite deviance are returned for information purposes.")
      div_message <- "Weighted-OLS returned NA coefficients."
      wols <- wols_old
      break
    } else {
      wols_means <- wols$means
    }

    eta <- wols$fitted.values

    if (isOffset) {
      eta <- eta + offset
    }

    if (mem.clean) {
      gc()
    }

    mu <- linkinv(eta)

    dev <- sum_dev.resids(y, mu, eta, wt = weights)
    dev_evol <- dev - devold

    if (verbose >= 1) cat("Iteration: ", sprintf("%02i", iter), " -- Deviance = ", numberFormatNormal(dev), " -- Evol. = ", dev_evol, "\n", sep = "")

    #
    # STEP HALVING
    #

    no_SH <- is.finite(dev) && (abs(dev_evol) < glm.tol || abs(dev_evol) / (0.1 + abs(dev)) < glm.tol)
    if (no_SH == FALSE && (!is.finite(dev) || dev_evol > 0 || !valideta(eta) || !validmu(mu))) {
      if (!is.finite(dev)) {
        # we report step-halving but only for non-finite deviances
        # other situations are OK (it just happens)
        nb_sh <- get("nb_sh", env)
        assign("nb_sh", nb_sh + 1, env)
      }

      eta_new <- wols$fitted.values
      eta_old <- wols_old$fitted.values

      iter_sh <- 0
      do_exit <- FALSE
      while (!is.finite(dev) || dev_evol > 0 || !valideta(eta_new) || !validmu(mu)) {
        if (iter == 1 && (is.finite(dev) && valideta(eta_new) && validmu(mu)) && iter_sh >= 2) {
          # BEWARE FIRST ITERATION:
          # at first iteration, the deviance can be higher than the init, and SH may not help
          # we need to make sure we get out of SH before it's messed up
          break
        } else if (iter_sh == glm.iter) {
          # if first iteration => means algo did not find viable solution
          if (iter == 1) {
            stop("Algorithm failed at first iteration. Step-halving could not find a valid set of parameters.")
          }

          # Problem only if the deviance is non-finite or eta/mu not valid
          # Otherwise, it means that we're at a maximum

          if (!is.finite(dev) || !valideta(eta_new) || !validmu(mu)) {
            # message
            msg <- ifelse(!is.finite(dev), "non-finite deviance", "no valid eta/mu")

            warning_msg <- paste0("Divergence at iteration ", iter, ": ", msg, ". Step halving: no valid correction found. Last evaluated coefficients with finite deviance are returned for information purposes.")
            div_message <- paste0(msg, " despite step-halving")
            wols <- wols_old
            do_exit <- TRUE
          }

          break
        }

        iter_sh <- iter_sh + 1
        eta_new <- (eta_old + eta_new) / 2

        if (mem.clean) {
          gc()
        }

        mu <- linkinv(eta_new + offset)
        dev <- sum_dev.resids(y, mu, eta_new + offset, wt = weights)
        dev_evol <- dev - devold

        if (verbose >= 3) cat("Step-halving: iter =", iter_sh, "-- dev:", numberFormatNormal(dev), "-- evol:", numberFormatNormal(dev_evol), "\n")
      }

      if (do_exit) break

      # it worked: update
      eta <- eta_new + offset
      wols$fitted.values <- eta_new
      # NOTA: we must NOT end with a step halving => we need a proper weighted-ols estimation
      # we force the algorithm to continue
      dev_evol <- Inf

      if (verbose >= 2) {
        cat("Step-halving: new deviance = ", numberFormatNormal(dev), "\n", sep = "")
      }
    }

    if (abs(dev_evol) / (0.1 + abs(dev)) < glm.tol) {
      conv <- TRUE
      break
    } else {
      devold <- dev
      wols_old <- wols
    }
  }

  # Convergence flag
  if (!conv) {
    if (iter == glm.iter) {
      warning_msg <- paste0("Absence of convergence: Maximum number of iterations reached (", glm.iter, "). Final deviance: ", numberFormatNormal(dev), ".")
      div_message <- "no convergence: Maximum number of iterations reached"
    }

    res$convStatus <- FALSE
    res$message <- div_message
  } else {
    res$convStatus <- TRUE
  }

  #
  # post processing
  #

  if (only.coef) {
    if (wols$multicol) {
      beta <- rep(NA_real_, length(wols$is_excluded))
      beta[!wols$is_excluded] <- wols$coefficients
    } else {
      beta <- wols$coefficients
    }

    names(beta) <- colnames(X)
    return(beta)
  }

  # Collinearity message
  collin.adj <- 0
  if (wols$multicol) {
    var_collinear <- colnames(X)[wols$is_excluded]
    if (notes) {
      msg <- dsb("w!The variable.[*s_, q, ' has'V, 3KO, C?var_collinear] been
                       removed because of collinearity (see $collin.var).")
      if (IN_MULTI) {
        stack_multi_notes(msg)
      } else {
        message(msg)
      }
    }

    res$collin.var <- var_collinear

    # full set of coeffficients with NAs
    collin.coef <- setNames(rep(NA, ncol(X)), colnames(X))
    collin.coef[!wols$is_excluded] <- wols$coefficients
    res$collin.coef <- collin.coef

    wols$X_demean <- wols$X_demean[, !wols$is_excluded, drop = FALSE]
    X <- X[, !wols$is_excluded, drop = FALSE]

    collin.adj <- sum(wols$is_excluded)
  }

  res$nparams <- res$nparams - collin.adj

  res$irls_weights <- w # weights from the iteratively reweighted least square

  res$coefficients <- coef <- wols$coefficients
  res$collin.min_norm <- wols$collin.min_norm

  if (!is.null(wols$warn_varying_slope)) {
    msg <- wols$warn_varying_slope
    if (IN_MULTI) {
      stack_multi_notes(msg)
    } else {
      warning(msg)
    }
  }

  res$linear.predictors <- wols$fitted.values
  if (isOffset) {
    res$linear.predictors <- res$linear.predictors + offset
  }

  res$fitted.values <- linkinv(res$linear.predictors)
  res$residuals <- y - res$fitted.values

  if (onlyFixef) res$onlyFixef <- onlyFixef

  # dispersion + scores
  if (family$family %in% c("poisson", "binomial")) {
    res$dispersion <- 1
  } else {
    weighted_resids <- wols$residuals * res$irls_weights
    # res$dispersion = sum(weighted_resids ** 2) / sum(res$irls_weights)
    # I use the second line to fit GLM's
    res$dispersion <- sum(weighted_resids * wols$residuals) / max(res$nobs - res$nparams, 1)
  }

  res$working_residuals <- wols$residuals

  if (!onlyFixef && !lean_internal) {
    # score + hessian + vcov

    if (mem.clean) {
      gc()
    }

    # dispersion + scores
    if (family$family %in% c("poisson", "binomial")) {
      res$scores <- (wols$residuals * res$irls_weights) * wols$X_demean
      res$hessian <- cpppar_crossprod(wols$X_demean, res$irls_weights, nthreads)
    } else {
      res$scores <- (weighted_resids / res$dispersion) * wols$X_demean
      res$hessian <- cpppar_crossprod(wols$X_demean, res$irls_weights, nthreads) / res$dispersion
    }

    if (any(diag(res$hessian) < 0)) {
      # This should not occur, but I prefer to be safe
      # In fact it's the opposite of the Hessian
      stop("Negative values in the diagonal of the Hessian found after the weighted-OLS stage. (If possible, could you send a replicable example to fixest's author? He's curious about when that actually happens, since in theory it should never happen.)")
    }

    # I put tol = 0, otherwise we may remove everything mistakenly
    # when VAR(Y) >>> VAR(X) // that is especially TRUE for missspecified
    # Poisson models
    info_inv <- cpp_cholesky(res$hessian, tol = 0, nthreads = nthreads)

    if (!is.null(info_inv$all_removed)) {
      # This should not occur, but I prefer to be safe
      stop("Not a single variable with a minimum of explanatory power found after the weighted-OLS stage. (If possible, could you send a replicable example to fixest's author? He's curious about when that actually happens, since in theory it should never happen.)")
    }

    var <- info_inv$XtX_inv
    is_excluded <- info_inv$id_excl

    if (any(is_excluded)) {
      # There should be no remaining collinearity
      warning_msg <- paste(warning_msg, "Residual collinearity was found after the weighted-OLS stage. The covariance is not defined. (This should not happen. If possible, could you send a replicable example to fixest's author? He's curious about when that actually happens.)")
      var <- matrix(NA, length(is_excluded), length(is_excluded))
    }
    res$cov.iid <- var

    rownames(res$cov.iid) <- colnames(res$cov.iid) <- names(coef)

    # for compatibility with conleyreg
    res$cov.unscaled <- res$cov.iid

    # se
    se <- diag(res$cov.iid)
    se[se < 0] <- NA
    se <- sqrt(se)

    # coeftable
    zvalue <- coef / se
    use_t <- !family$family %in% c("poisson", "binomial")
    if (use_t) {
      pvalue <- 2 * pt(-abs(zvalue), max(res$nobs - res$nparams, 1))
      ctable_names <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    } else {
      pvalue <- 2 * pnorm(-abs(zvalue))
      ctable_names <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    }

    coeftable <- data.frame("Estimate" = coef, "Std. Error" = se, "z value" = zvalue, "Pr(>|z|)" = pvalue)
    names(coeftable) <- ctable_names
    row.names(coeftable) <- names(coef)

    attr(se, "type") <- attr(coeftable, "type") <- "IID"
    res$coeftable <- coeftable
    res$se <- se
  }

  if (nchar(warning_msg) > 0) {
    if (warn) {
      if (IN_MULTI) {
        stack_multi_notes(warning_msg)
      } else {
        warning(warning_msg, call. = FALSE)
        options("fixest_last_warning" = proc.time())
      }
    }
  }

  n <- length(y)
  res$nobs <- n

  df_k <- res$nparams

  # r2s
  if (!cpp_isConstant(res$fitted.values)) {
    res$sq.cor <- tryCatch(stats::cor(y, res$fitted.values), warning = function(x) NA_real_)**2
  } else {
    res$sq.cor <- NA
  }

  # deviance
  res$deviance <- dev

  # simpler form for poisson
  if (family_equiv == "poisson") {
    if (isWeight) {
      if (mem.clean) {
        gc()
      }

      res$loglik <- sum((y * eta - mu - cpppar_lgamma(y + 1, nthreads)) * weights)
    } else {
      # lfact is later used in model0 and is costly to compute
      lfact <- sum(rpar_lgamma(y + 1, env))
      assign("lfactorial", lfact, env)
      res$loglik <- sum(y * eta - mu) - lfact
    }
  } else {
    res$loglik <- family$aic(y = y, n = rep.int(1, n), mu = res$fitted.values, wt = weights, dev = dev) / -2
  }

  if (lean_internal) {
    return(res)
  }

  # The pseudo_r2
  if (family_equiv %in% c("poisson", "logit")) {
    model0 <- get_model_null(env, theta.init = NULL)
    ll_null <- model0$loglik
    fitted_null <- linkinv(model0$constant)
  } else {
    if (verbose >= 1) cat("Null model:\n")

    if (mem.clean) {
      gc()
    }

    model_null <- feglm.fit(X = matrix(1, nrow = n, ncol = 1), fixef_df = NULL, env = env, lean_internal = TRUE)
    ll_null <- model_null$loglik
    fitted_null <- model_null$fitted.values
  }
  res$ll_null <- ll_null
  res$pseudo_r2 <- 1 - (res$loglik - df_k) / (ll_null - 1)

  # fixef info
  if (isFixef) {
    if (onlyFixef) {
      res$sumFE <- res$linear.predictors
    } else {
      res$sumFE <- res$linear.predictors - cpppar_xbeta(X, res$coefficients, nthreads)
    }

    if (isOffset) {
      res$sumFE <- res$sumFE - offset
    }
  }

  # other
  res$iterations <- iter
  class(res) <- "fixest"

  do_summary <- get("do_summary", env)
  if (do_summary) {
    vcov <- get("vcov", env)
    lean <- get("lean", env)
    ssc <- get("ssc", env)
    agg <- get("agg", env)
    summary_flags <- get("summary_flags", env)

    # To compute the RMSE and lean = TRUE
    if (lean) res$ssr <- cpp_ssq(res$residuals, weights)

    res <- summary(res, vcov = vcov, agg = agg, ssc = ssc, lean = lean, summary_flags = summary_flags)
  }

  return(res)
}
