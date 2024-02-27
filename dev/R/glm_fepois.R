#' @rdname feglm
#' @export
fepois <- function(fml, data, vcov, offset, weights, subset, split, fsplit, split.keep, split.drop,
                   cluster, se, ssc, panel.id, start = NULL, etastart = NULL,
                   mustart = NULL, fixef, fixef.rm = "perfect", fixef.tol = 1e-6,
                   fixef.iter = 10000, collin.tol = 1e-10, glm.iter = 25, glm.tol = 1e-8,
                   nthreads = getFixest_nthreads(), lean = FALSE, warn = TRUE, notes = getFixest_notes(),
                   verbose = 0, combine.quick, mem.clean = FALSE, only.env = FALSE,
                   only.coef = FALSE, env, ...) {
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
    collin.tol = collin.tol, glm.iter = glm.iter, glm.tol = glm.tol,
    nthreads = nthreads, lean = lean, warn = warn, notes = notes,
    verbose = verbose, combine.quick = combine.quick, mem.clean = mem.clean,
    only.env = only.env, only.coef = only.coef,
    origin_bis = "fepois", mc_origin_bis = match.call(),
    call_env_bis = call_env_bis, env = env, ...
  ), silent = TRUE)

  if ("try-error" %in% class(res)) {
    stop(format_error_msg(res, "fepois"))
  }

  return(res)
}
