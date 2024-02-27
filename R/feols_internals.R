# OLS internals ----

ols_fit <- function(y, X, w, correct_0w = FALSE, collin.tol, nthreads, xwx = NULL,
                    xwy = NULL, only.coef = FALSE) {
  # No control here -- done before

  if (is.null(xwx)) {
    info_products <- cpp_sparse_products(X, w, y, correct_0w, nthreads)
    xwx <- info_products$XtX
    xwy <- info_products$Xty
  }

  multicol <- FALSE
  info_inv <- cpp_cholesky(xwx, collin.tol, nthreads)

  if (!is.null(info_inv$all_removed)) {
    # Means all variables are collinear! => can happen when using FEs
    return(list(all_removed = TRUE))
  }

  xwx_inv <- info_inv$XtX_inv
  is_excluded <- info_inv$id_excl

  multicol <- any(is_excluded)

  if (only.coef) {
    if (multicol) {
      beta <- rep(NA_real_, length(is_excluded))
      beta[!is_excluded] <- as.vector(xwx_inv %*% xwy[!is_excluded])
    } else {
      beta <- as.vector(xwx_inv %*% xwy)
    }

    return(beta)
  }

  if (multicol) {
    beta <- as.vector(xwx_inv %*% xwy[!is_excluded])
    fitted.values <- cpp_xbeta(X[, !is_excluded, drop = FALSE], beta, nthreads)
  } else {
    # avoids copies
    beta <- as.vector(xwx_inv %*% xwy)
    fitted.values <- cpp_xbeta(X, beta, nthreads)
  }


  residuals <- y - fitted.values

  res <- list(
    xwx = xwx, coefficients = beta, fitted.values = fitted.values,
    xwx_inv = xwx_inv, multicol = multicol, residuals = residuals,
    is_excluded = is_excluded, collin.min_norm = info_inv$min_norm
  )

  res
}

check_conv <- function(y, X, fixef_id_list, slope_flag, slope_vars, weights, full = FALSE, fixef_names = NULL) {
  # VERY SLOW!!!!
  # IF THIS FUNCTION LASTS => TO BE PORTED TO C++

  # y, X => variables that were demeaned

  # For each variable: we compute the optimal FE coefficient
  # it should be 0 if the algorithm converged

  Q <- length(slope_flag)

  use_y <- TRUE
  if (length(X) == 1) {
    K <- 1
    nobs <- length(y)
  } else {
    nobs <- NROW(X)

    if (length(y) <= 1) {
      use_y <- FALSE
    }

    K <- NCOL(X) + use_y
  }

  info_max <- list()
  info_full <- vector("list", K)

  for (k in 1:K) {
    if (k == 1 && use_y) {
      x <- y
    } else if (use_y) {
      x <- X[, k - use_y]
    } else {
      x <- X[, k]
    }

    info_max_tmp <- c()
    info_full_tmp <- list()
    index_slope <- 1
    for (q in 1:Q) {
      fixef_id <- fixef_id_list[[q]]

      if (slope_flag[q] >= 0) {
        x_means <- tapply(weights * x, fixef_id, mean)
        info_max_tmp <- c(info_max_tmp, max(abs(x_means)))

        if (full) {
          info_full_tmp[[length(info_full_tmp) + 1]] <- x_means
        }
      }

      n_slopes <- abs(slope_flag[q])
      if (n_slopes > 0) {
        for (i in 1:n_slopes) {
          var <- slope_vars[[index_slope]]

          num <- tapply(weights * x * var, fixef_id, sum)
          denom <- tapply(weights * var^2, fixef_id, sum)

          x_means <- num / denom
          info_max_tmp <- c(info_max_tmp, max(abs(x_means)))

          if (full) {
            info_full_tmp[[length(info_full_tmp) + 1]] <- x_means
          }

          index_slope <- index_slope + 1
        }
      }
    }

    if (!is.null(fixef_names)) {
      names(info_full_tmp) <- fixef_names
    }

    info_max[[k]] <- info_max_tmp

    if (full) {
      info_full[[k]] <- info_full_tmp
    }
  }

  info_max <- do.call("rbind", info_max)

  if (full) {
    res <- info_full
  } else {
    res <- info_max
  }

  res
}


#' @rdname feols
#' @export
feols.fit <- function(y, X, fixef_df, vcov, offset, split, fsplit, split.keep, split.drop,
                      cluster, se, ssc, weights,
                      subset, fixef.rm = "perfect", fixef.tol = 1e-6, fixef.iter = 10000,
                      fixef.algo = NULL, collin.tol = 1e-10,
                      nthreads = getFixest_nthreads(), lean = FALSE,
                      warn = TRUE, notes = getFixest_notes(), mem.clean = FALSE, verbose = 0,
                      only.env = FALSE, only.coef = FALSE, env, ...) {
  if (missing(weights)) weights <- NULL

  time_start <- proc.time()

  if (missing(env)) {
    set_defaults("fixest_estimation")
    call_env <- new.env(parent = parent.frame())

    env <- try(fixest_env(
      y = y, X = X, fixef_df = fixef_df, vcov = vcov, offset = offset,
      split = split, fsplit = fsplit,
      split.keep = split.keep, split.drop = split.drop,
      cluster = cluster, se = se, ssc = ssc,
      weights = weights, subset = subset, fixef.rm = fixef.rm,
      fixef.tol = fixef.tol, fixef.iter = fixef.iter,
      fixef.algo = fixef.algo, collin.tol = collin.tol,
      nthreads = nthreads, lean = lean, warn = warn, notes = notes,
      mem.clean = mem.clean, verbose = verbose, only.coef = only.coef,
      origin = "feols.fit",
      mc_origin = match.call(), call_env = call_env, ...
    ), silent = TRUE)
  } else if ((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
    stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
  }

  if ("try-error" %in% class(env)) {
    mc <- match.call()
    origin <- ifelse(is.null(mc$origin), "feols.fit", mc$origin)
    stop(format_error_msg(env, origin))
  }

  check_arg(only.env, "logical scalar")
  if (only.env) {
    return(env)
  }

  verbose <- get("verbose", env)
  if (verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep = "")

  # workhorse is feols (OK if error msg leads to feols [clear enough])
  res <- feols(env = env)

  res
}

#' Estimates a `fixest` estimation from a `fixest` environment
#'
#' This is a function advanced users which allows to estimate any `fixest` estimation from a
#' `fixest` environment obtained with `only.env = TRUE` in a `fixest` estimation.
#'
#' @param env An environment obtained from a `fixest` estimation with `only.env = TRUE`. This is
#'  intended for advanced users so there is no error handling: any other kind of input will
#' fail with a poor error message.
#' @param y A vector representing the dependent variable. Should be of the same length
#' as the number of observations in the initial estimation.
#' @param X A matrix representing the independent variables. Should be of the same dimension
#' as in the initial estimation.
#' @param weights A vector of weights (i.e. with only positive values). Should be of
#' the same length as the number of observations in the initial estimation. If identical
#' to the scalar 1, this will mean that no weights will be used in the estimation.
#' @param endo A matrix representing the endogenous regressors in IV estimations. It should
#' be of the same dimension as the original endogenous regressors.
#' @param inst A matrix representing the instruments in IV estimations. It should be of
#' the same dimension as the original instruments.
#'
#' @return
#'
#' It returns the results of a `fixest` estimation: the one that was summoned when
#' obtaining the environment.
#'
#' @details
#'
#' This function has been created for advanced users, mostly to avoid overheads
#' when making simulations with `fixest`.
#'
#' How can it help you make simulations? First make a core estimation with `only.env = TRUE`,
#' and usually with `only.coef = TRUE` (to avoid having extra things that take time to compute).
#' Then loop while modifying the appropriate things directly in the environment. Beware that
#' if you make a mistake here (typically giving stuff of the wrong length),
#' then you can make the R session crash because there is no more error-handling!
#' Finally estimate with `est_env(env = core_env)` and store the results.
#'
#' Instead of `est_env`, you could use directly `fixest` estimations too, like `feols`,
#' since they accept the `env` argument. The function `est_env` is only here to add a
#' bit of generality to avoid the trouble to the user to write conditions
#' (look at the source, it's just a one liner).
#'
#' Objects of main interest in the environment are:
#' \describe{
#' \item{lhs}{The left hand side, or dependent variable.}
#' \item{linear.mat}{The matrix of the right-hand-side, or explanatory variables.}
#' \item{iv_lhs}{The matrix of the endogenous variables in IV regressions.}
#' \item{iv.mat}{The matrix of the instruments in IV regressions.}
#' \item{weights.value}{The vector of weights.}
#' }
#'
#' I strongly discourage changing the dimension of any of these elements, or else crash can occur.
#' However, you can change their values at will (given the dimension stay the same).
#' The only exception is the weights, which tolerates changing its dimension: it can
#' be identical to the scalar `1` (meaning no weights), or to something of the length the
#' number of observations.
#'
#' I also discourage changing anything in the fixed-effects (even their value)
#' since this will almost surely lead to a crash.
#'
#' Note that this function is mostly useful when the overheads/estimation ratio is high.
#' This means that OLS will benefit the most from this function. For GLM/Max.Lik. estimations,
#' the ratio is small since the overheads is only a tiny portion of the total estimation time.
#' Hence this function will be less useful for these models.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Let's make a short simulation
#' # Inspired from Grant McDermott bboot function
#' # See https://twitter.com/grant_mcdermott/status/1487528757418102787
#'
#' # Simple function that computes a Bayesian bootstrap
#' bboot <- function(x, n_sim = 100) {
#'   # We bootstrap on the weights
#'   # Works with fixed-effects/IVs
#'   #  and with any fixest function that accepts weights
#'
#'   core_env <- update(x, only.coef = TRUE, only.env = TRUE)
#'   n_obs <- x$nobs
#'
#'   res_all <- vector("list", n_sim)
#'   for (i in 1:n_sim) {
#'     ## begin: NOT RUN
#'     ## We could directly assign in the environment:
#'     # assign("weights.value", rexp(n_obs, rate = 1), core_env)
#'     # res_all[[i]] = est_env(env = core_env)
#'     ##   end: NOT RUN
#'
#'     ## Instead we can use the argument weights, which does the same
#'     res_all[[i]] <- est_env(env = core_env, weights = rexp(n_obs, rate = 1))
#'   }
#'
#'   do.call(rbind, res_all)
#' }
#'
#'
#' est <- feols(mpg ~ wt + hp, mtcars)
#'
#' boot_res <- bboot(est)
#' coef <- colMeans(boot_res)
#' std_err <- apply(boot_res, 2, sd)
#'
#' # Comparing the results with the main estimation
#' coeftable(est)
#' cbind(coef, std_err)
#'
#' @export
est_env <- function(env, y, X, weights, endo, inst) {
  # No check whatsoever: for advanced users

  if (!missnull(y)) {
    assign("lhs", y, env)
  }

  if (!missnull(X)) {
    assign("linear.mat", X, env)
  }

  if (!missnull(weights)) {
    assign("weights.value", weights, env)
  }

  if (!missnull(endo)) {
    assign("is_lhs", endo, env)
  }

  if (!missnull(inst)) {
    assign("iv.mat", inst, env)
  }

  switch(env$res$method_type,
    "feols"  = feols(env = env),
    "feglm"  = feglm.fit(env = env),
    "feNmlm" = feNmlm(env = env)
  )
}
