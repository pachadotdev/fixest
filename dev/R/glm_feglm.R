#' Fixed-effects GLM estimations
#'
#' Estimates GLM models with any number of fixed-effects.
#'
#' @inheritParams femlm
#' @inheritParams feols
#' @inheritSection feols Combining the fixed-effects
#' @inheritSection feols Varying slopes
#' @inheritSection feols Lagging variables
#' @inheritSection feols Interactions
#' @inheritSection feols On standard-errors
#' @inheritSection feols Multiple estimations
#' @inheritSection feols Argument sliding
#' @inheritSection feols Piping
#' @inheritSection feols Tricks to estimate multiple LHS
#' @inheritSection xpd Dot square bracket operator in formulas
#'
#' @param family Family to be used for the estimation. Defaults to `gaussian()`. See [`family`] for details of family functions.
#' @param start Starting values for the coefficients. Can be: i) a numeric of length 1 (e.g. `start = 0`), ii) a numeric vector of the exact same length as the number of variables, or iii) a named vector of any length (the names will be used to initialize the appropriate coefficients). Default is missing.
#' @param etastart Numeric vector of the same length as the data. Starting values for the linear predictor. Default is missing.
#' @param mustart Numeric vector of the same length as the data. Starting values for the vector of means. Default is missing.
#' @param fixef.tol Precision used to obtain the fixed-effects. Defaults to `1e-6`. It corresponds to the maximum absolute difference allowed between two coefficients of successive iterations.
#' @param glm.iter Number of iterations of the glm algorithm. Default is 25.
#' @param glm.tol Tolerance level for the glm algorithm. Default is `1e-8`.
#' @param verbose Integer. Higher values give more information. In particular, it can detail the number of iterations in the demeaning algoritmh (the first number is the left-hand-side, the other numbers are the right-hand-side variables). It can also detail the step-halving algorithm.
#' @param notes Logical. By default, three notes are displayed: when NAs are removed, when some fixed-effects are removed because of only 0 (or 0/1) outcomes, or when a variable is dropped because of collinearity. To avoid displaying these messages, you can set `notes = FALSE`. You can remove these messages permanently by using `setFixest_notes(FALSE)`.
#'
#' @details
#' The core of the GLM are the weighted OLS estimations. These estimations are performed with [`feols`]. The method used to demean each variable along the fixed-effects is based on Berge (2018), since this is the same problem to solve as for the Gaussian case in a ML setup.
#'
#' @return
#' A `fixest` object. Note that `fixest` objects contain many elements and most of them are for internal use, they are presented here only for information. To access them, it is safer to use the user-level methods (e.g. [`vcov.fixest`], [`resid.fixest`], etc) or functions (like for instance [`fitstat`] to access any fit statistic).
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{fml_all}{A list containing different parts of the formula. Always contain the linear formula. Then, if relevant: `fixef`: the fixed-effects.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{y}{(When relevant.) The dependent variable (used to compute the within-R2 when fixed-effects are present).}
#' \item{convStatus}{Logical, convergence status of the IRWLS algorithm.}
#' \item{irls_weights}{The weights of the last iteration of the IRWLS algorithm.}
#' \item{obs_selection}{(When relevant.) List containing vectors of integers. It represents the sequential selection of observation vis a vis the original data set.}
#' \item{fixef_removed}{(When relevant.) In the case there were fixed-effects and some observations were removed because of only 0/1 outcome within a fixed-effect, it gives the list (for each fixed-effect dimension) of the fixed-effect identifiers that were removed.}
#' \item{coefficients}{The named vector of estimated coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The loglikelihood.}
#' \item{deviance}{Deviance of the fitted model.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
#' \item{ssr_null}{Sum of the squared residuals of the null model (containing only with the intercept).}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{fitted.values}{The fitted values are the expected value of the dependent variable for the fitted model: that is \eqn{E(Y|X)}.}
#' \item{linear.predictors}{The linear predictors.}
#' \item{residuals}{The residuals (y minus the fitted values).}
#' \item{sq.cor}{Squared correlation between the dependent variable and the expected predictor (i.e. fitted.values) obtained by the estimation.}
#' \item{hessian}{The Hessian of the parameters.}
#' \item{cov.iid}{The variance-covariance matrix of the parameters.}
#' \item{se}{The standard-error of the parameters.}
#' \item{scores}{The matrix of the scores (first derivative for each observation).}
#' \item{residuals}{The difference between the dependent variable and the expected predictor.}
#' \item{sumFE}{The sum of the fixed-effects coefficients for each observation.}
#' \item{offset}{(When relevant.) The offset formula.}
#' \item{weights}{(When relevant.) The weights formula.}
#' \item{collin.var}{(When relevant.) Vector containing the variables removed because of collinearity.}
#' \item{collin.coef}{(When relevant.) Vector of coefficients, where the values of the variables removed because of collinearity are NA.}
#'
#'
#'
#'
#' @seealso
#' See also [`summary.fixest`] to see the results with the appropriate standard-errors, [`fixef.fixest`] to extract the fixed-effects coefficients.
#' And other estimation methods: [`feols`], [`femlm`], [`fenegbin`], [`feNmlm`].
#'
#' @author
#' Laurent Berge
#'
#' @references
#'
#' Berge, Laurent, 2018, "Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm." CREA Discussion Papers, 13 ([](https://wwwen.uni.lu/content/download/110162/1299525/file/2018_13)).
#'
#' For models with multiple fixed-effects:
#'
#' Gaure, Simen, 2013, "OLS with multiple high dimensional category variables", Computational Statistics & Data Analysis 66 pp. 8--18
#'
#'
#' @examples
#' # Single estimation ----
#'
#' # Poisson estimation
#' res <- feglm(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris, "poisson")
#'
#' # You could also use fepois
#' res_pois <- fepois(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris)
#'
#' # With the fit method:
#' res_fit <- feglm.fit(iris$Sepal.Length, iris[, 2:3], iris$Species, "poisson")
#'
#' # Note that you have many more examples in feols
#'
#' # Multiple estimations ----
#'
#' # 6 estimations
#' est_mult <- fepois(c(Ozone, Solar.R) ~ Wind + Temp + csw0(Wind:Temp, Day), airquality)
#'
#' # Combining with split
#' est_split <- fepois(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
#'   airquality,
#'   split = ~Month
#' )
#'
#' # You can display everything at once with the print method
#' est_split
#'
#' # Different way of displaying the results with "compact"
#' summary(est_split, "compact")
#'
#' # You can still select which sample/LHS/RHS to display
#' est_split[sample = 1:2, lhs = 1, rhs = 1]
#'
#' @export
feglm <- function(fml, data, family = "gaussian", vcov, offset, weights, subset, split,
                  fsplit, split.keep, split.drop, cluster, se, ssc, panel.id, start = NULL,
                  etastart = NULL, mustart = NULL, fixef, fixef.rm = "perfect",
                  fixef.tol = 1e-6, fixef.iter = 10000, collin.tol = 1e-10,
                  glm.iter = 25, glm.tol = 1e-8, nthreads = getFixest_nthreads(),
                  lean = FALSE, warn = TRUE, notes = getFixest_notes(), verbose = 0,
                  only.coef = FALSE,
                  combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...) {
  if (missing(weights)) weights <- NULL

  time_start <- proc.time()

  if (missing(env)) {
    set_defaults("fixest_estimation")
    call_env <- new.env(parent = parent.frame())

    env <- try(fixest_env(
      fml = fml, data = data, family = family,
      offset = offset, weights = weights, subset = subset,
      split = split, fsplit = fsplit,
      split.keep = split.keep, split.drop = split.drop,
      vcov = vcov,
      cluster = cluster, se = se, ssc = ssc,
      panel.id = panel.id, linear.start = start,
      etastart = etastart, mustart = mustart, fixef = fixef,
      fixef.rm = fixef.rm, fixef.tol = fixef.tol,
      fixef.iter = fixef.iter, collin.tol = collin.tol,
      glm.iter = glm.iter, glm.tol = glm.tol,
      nthreads = nthreads, lean = lean, warn = warn, notes = notes,
      verbose = verbose, combine.quick = combine.quick,
      mem.clean = mem.clean, only.coef = only.coef, origin = "feglm",
      mc_origin = match.call(), call_env = call_env, ...
    ), silent = TRUE)
  } else if ((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
    stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
  }

  if ("try-error" %in% class(env)) {
    mc <- match.call()
    origin <- ifelse(is.null(mc$origin), "feglm", mc$origin)
    stop(format_error_msg(env, origin))
  }

  check_arg(only.env, "logical scalar")
  if (only.env) {
    return(env)
  }

  verbose <- get("verbose", env)
  if (verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep = "")

  # workhorse is feglm.fit (OK if error msg leads to feglm.fit [clear enough])
  res <- feglm.fit(env = env)

  res
}
