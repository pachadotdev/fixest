# GLM Wrappers ----

#' @rdname femlm
#' @export
fenegbin <- function(fml, data, vcov, theta.init, start = 0, fixef, fixef.rm = "perfect",
                     offset, subset, split, fsplit, split.keep, split.drop,
                     cluster, se, ssc, panel.id,
                     fixef.tol = 1e-5, fixef.iter = 10000, nthreads = getFixest_nthreads(),
                     lean = FALSE, verbose = 0, warn = TRUE, notes = getFixest_notes(),
                     combine.quick, mem.clean = FALSE,
                     only.env = FALSE, only.coef = FALSE, data.save = FALSE, env, ...) {
  # We control for the problematic argument family
  if ("family" %in% names(match.call())) {
    stop("Function fenegbin does not accept the argument 'family'.")
  }

  # This is just an alias
  set_defaults("fixest_estimation")
  call_env_bis <- new.env(parent = parent.frame())

  res <- try(feNmlm(
    fml = fml, data = data, family = "negbin", theta.init = theta.init,
    start = start, fixef = fixef, fixef.rm = fixef.rm, offset = offset,
    subset = subset, split = split, fsplit = fsplit,
    split.keep = split.keep, split.drop = split.drop,
    vcov = vcov, cluster = cluster,
    se = se, ssc = ssc, panel.id = panel.id, fixef.tol = fixef.tol,
    fixef.iter = fixef.iter, nthreads = nthreads, lean = lean,
    verbose = verbose, warn = warn, notes = notes, combine.quick = combine.quick,
    mem.clean = mem.clean, only.env = only.env,
    only.coef = only.coef, data.save = data.save,
    origin = "fenegbin", mc_origin_bis = match.call(),
    call_env_bis = call_env_bis, env = env, ...
  ), silent = TRUE)

  if ("try-error" %in% class(res)) {
    stop(format_error_msg(res, "fenegbin"))
  }

  return(res)
}

#' @rdname feglm
#' @export
fepois <- function(fml, data, vcov, offset, weights, subset, split, fsplit,
                   split.keep, split.drop,
                   cluster, se, ssc, panel.id, start = NULL, etastart = NULL,
                   mustart = NULL, fixef, fixef.rm = "perfect", fixef.tol = 1e-6,
                   fixef.iter = 10000, fixef.algo = NULL,
                   collin.tol = 1e-10, glm.iter = 25, glm.tol = 1e-8,
                   nthreads = getFixest_nthreads(), lean = FALSE,
                   warn = TRUE, notes = getFixest_notes(),
                   verbose = 0, combine.quick, mem.clean = FALSE, only.env = FALSE,
                   only.coef = FALSE, data.save = FALSE, env, ...) {
  # We control for the problematic argument family
  if ("family" %in% names(match.call())) {
    stop("Function fepois does not accept the argument 'family'.")
  }

  # This is just an alias
  set_defaults("fixest_estimation")
  call_env_bis <- new.env(parent = parent.frame())

  res <- try(feglm(
    fml = fml, data = data, family = "poisson", offset = offset,
    weights = weights, subset = subset, split = split, fsplit = fsplit,
    split.keep = split.keep, split.drop = split.drop,
    vcov = vcov, cluster = cluster, se = se, ssc = ssc, panel.id = panel.id,
    start = start, etastart = etastart, mustart = mustart, fixef = fixef,
    fixef.rm = fixef.rm, fixef.tol = fixef.tol, fixef.iter = fixef.iter,
    fixef.algo = fixef.algo,
    collin.tol = collin.tol, glm.iter = glm.iter, glm.tol = glm.tol,
    nthreads = nthreads, lean = lean, warn = warn, notes = notes,
    verbose = verbose, combine.quick = combine.quick, mem.clean = mem.clean,
    only.env = only.env, only.coef = only.coef, data.save = data.save,
    origin_bis = "fepois", mc_origin_bis = match.call(),
    call_env_bis = call_env_bis, env = env, ...
  ), silent = TRUE)

  if ("try-error" %in% class(res)) {
    stop(format_error_msg(res, "fepois"))
  }

  return(res)
}
