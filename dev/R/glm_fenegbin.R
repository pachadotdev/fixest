#' @rdname femlm
#' @export
fenegbin <- function(fml, data, vcov, theta.init, start = 0, fixef, fixef.rm = "perfect",
                     offset, subset, split, fsplit, split.keep, split.drop,
                     cluster, se, ssc, panel.id,
                     fixef.tol = 1e-5, fixef.iter = 10000, nthreads = getFixest_nthreads(),
                     lean = FALSE, verbose = 0, warn = TRUE, notes = getFixest_notes(),
                     combine.quick, mem.clean = FALSE, only.env = FALSE, only.coef = FALSE, env, ...) {
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
    mem.clean = mem.clean, only.env = only.env, only.coef = only.coef,
    origin = "fenegbin", mc_origin_bis = match.call(),
    call_env_bis = call_env_bis, env = env, ...
  ), silent = TRUE)

  if ("try-error" %in% class(res)) {
    stop(format_error_msg(res, "fenegbin"))
  }

  return(res)
}
