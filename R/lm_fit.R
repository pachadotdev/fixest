#' @rdname feols
#' @export
feols.fit <- function(y, X, fixef_df, vcov, offset, split, fsplit, split.keep, split.drop,
                      cluster, se, ssc, weights,
                      subset, fixef.rm = "perfect", fixef.tol = 1e-6, fixef.iter = 10000,
                      collin.tol = 1e-10, nthreads = getFixest_nthreads(), lean = FALSE,
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
      fixef.tol = fixef.tol, fixef.iter = fixef.iter, collin.tol = collin.tol,
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
