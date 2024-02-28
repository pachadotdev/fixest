# FE MLM ----

#' Fixed-effects maximum likelihood models
#'
#' This function estimates maximum likelihood models with any number of fixed-effects.
#'
#' @inheritParams feNmlm
#' @inherit feNmlm return details
#' @inheritSection feols Combining the fixed-effects
#' @inheritSection feols Lagging variables
#' @inheritSection feols Interactions
#' @inheritSection feols On standard-errors
#' @inheritSection feols Multiple estimations
#' @inheritSection feols Argument sliding
#' @inheritSection feols Piping
#' @inheritSection feols Tricks to estimate multiple LHS
#' @inheritSection xpd Dot square bracket operator in formulas
#'
#' @param fml A formula representing the relation to be estimated. For example: `fml = z~x+y`.
#' To include fixed-effects, insert them in this formula using a pipe: e.g.
#' `fml = z~x+y|fixef_1+fixef_2`. Multiple estimations can be performed at once:
#' for multiple dep. vars, wrap them in `c()`: ex `c(y1, y2)`. For multiple indep.
#' vars, use the stepwise functions: ex `x1 + csw(x2, x3)`.
#' The formula `fml = c(y1, y2) ~ x1 + cw0(x2, x3)` leads to 6 estimation, see details.
#' Square brackets starting with a dot can be used to call global variables:
#' `y.[i] ~ x.[1:2]` will lead to `y3 ~ x1 + x2` if `i` is equal to 3 in
#' the current environment (see details in [`xpd`]).
#' @param start Starting values for the coefficients. Can be: i) a numeric of length 1
#' (e.g. `start = 0`, the default), ii) a numeric vector of the exact same length as the
#' number of variables, or iii) a named vector of any length (the names will be
#' used to initialize the appropriate coefficients).
#'
#' @details
#' Note that the functions [`feglm`] and [`femlm`] provide the same results when using
#' the same families but differ in that the latter is a direct maximum likelihood
#' optimization (so the two can really have different convergence rates).
#'
#' @return
#' A `fixest2` object. Note that `fixest2` objects contain many elements and most of
#' them are for internal use, they are presented here only for information.
#' To access them, it is safer to use the user-level methods
#' (e.g. [`vcov.fixest2`], [`resid.fixest2`], etc) or functions (like for instance
#' [`fitstat`] to access any fit statistic).
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{data}{The original data set used when calling the function. Only available when
#' the estimation was called with `data.save = TRUE`}
#' \item{fml_all}{A list containing different parts of the formula. Always contain the
#' linear formula. Then, if relevant: `fixef`: the fixed-effects;
#' `NL`: the non linear part of the formula.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the
#' fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique
#' identifier for each fixed-effect dimension).}
#' \item{convStatus}{Logical, convergence status.}
#' \item{message}{The convergence message from the optimization procedures.}
#' \item{obs_selection}{(When relevant.) List containing vectors of integers. It represents
#' the sequential selection of observation vis a vis the original data set.}
#' \item{fixef_removed}{(When relevant.) In the case there were fixed-effects and some
#' observations were removed because of only 0/1 outcome within a fixed-effect, it gives the
#' list (for each fixed-effect dimension) of the fixed-effect identifiers that were removed.}
#' \item{coefficients}{The named vector of estimated coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values
#' and p-values.}
#' \item{loglik}{The log-likelihood.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
#' \item{ll_fe_only}{Log-likelihood of the model with only the fixed-effects.}
#' \item{ssr_null}{Sum of the squared residuals of the null model (containing only with
#' the intercept).}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{fitted.values}{The fitted values are the expected value of the dependent variable
#' for the fitted model: that is \eqn{E(Y|X)}.}
#' \item{residuals}{The residuals (y minus the fitted values).}
#' \item{sq.cor}{Squared correlation between the dependent variable and the
#' expected predictor (i.e. fitted.values) obtained by the estimation.}
#' \item{hessian}{The Hessian of the parameters.}
#' \item{cov.iid}{The variance-covariance matrix of the parameters.}
#' \item{se}{The standard-error of the parameters.}
#' \item{scores}{The matrix of the scores (first derivative for each observation).}
#' \item{residuals}{The difference between the dependent variable and the expected predictor.}
#' \item{sumFE}{The sum of the fixed-effects coefficients for each observation.}
#' \item{offset}{(When relevant.) The offset formula.}
#' \item{weights}{(When relevant.) The weights formula.}
#'
#'
#' @seealso
#' See also [`summary.fixest2`] to see the results with the appropriate standard-errors,
#' [`fixef.fixest2`] to extract the fixed-effects coefficients, and the function
#' [`etable`] to visualize the results of multiple estimations.
#' And other estimation methods: [`feols`], [`feglm`], [`fepois`], [`feNmlm`].
#'
#' @author
#' Laurent Berge
#'
#' @references
#'
#' Berge, Laurent, 2018, "Efficient estimation of maximum likelihood models with
#' multiple fixed-effects: the R package FENmlm." CREA Discussion Papers,
#' 13 ([](https://github.com/lrberge/fixest2/blob/master/_DOCS/FENmlm_paper.pdf)).
#'
#' For models with multiple fixed-effects:
#'
#' Gaure, Simen, 2013, "OLS with multiple high dimensional category variables",
#' Computational Statistics & Data Analysis 66 pp. 8--18
#'
#' On the unconditionnal Negative Binomial model:
#'
#' Allison, Paul D and Waterman, Richard P, 2002, "Fixed-Effects Negative
#' Binomial Regression Models", Sociological Methodology 32(1) pp. 247--265
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade => we account for 3 fixed-effects
#' # 1) Poisson estimation
#' est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination + Product, trade)
#'
#' # 2) Log-Log Gaussian estimation (with same FEs)
#' est_gaus <- update(est_pois, log(Euros + 1) ~ ., family = "gaussian")
#'
#' # Comparison of the results using the function etable
#' etable(est_pois, est_gaus)
#' # Now using two way clustered standard-errors
#' etable(est_pois, est_gaus, se = "twoway")
#'
#' # Comparing different types of standard errors
#' sum_hetero <- summary(est_pois, se = "hetero")
#' sum_oneway <- summary(est_pois, se = "cluster")
#' sum_twoway <- summary(est_pois, se = "twoway")
#' sum_threeway <- summary(est_pois, se = "threeway")
#'
#' etable(sum_hetero, sum_oneway, sum_twoway, sum_threeway)
#'
#'
#' #
#' # Multiple estimations:
#' #
#'
#' # 6 estimations
#' est_mult <- femlm(c(Ozone, Solar.R) ~ Wind + Temp + csw0(Wind:Temp, Day), airquality)
#'
#' # We can display the results for the first lhs:
#' etable(est_mult[lhs = 1])
#'
#' # And now the second (access can be made by name)
#' etable(est_mult[lhs = "Solar.R"])
#'
#' # Now we focus on the two last right hand sides
#' # (note that .N can be used to specify the last item)
#' etable(est_mult[rhs = 2:.N])
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
femlm <- function(fml, data, family = c("poisson", "negbin", "logit", "gaussian"), vcov,
                  start = 0, fixef, fixef.rm = "perfect", offset, subset,
                  split, fsplit, split.keep, split.drop,
                  cluster, se, ssc, panel.id, fixef.tol = 1e-5, fixef.iter = 10000,
                  nthreads = getFixest_nthreads(), lean = FALSE, verbose = 0, warn = TRUE,
                  notes = getFixest_notes(), theta.init, combine.quick, mem.clean = FALSE,
                  only.env = FALSE, only.coef = FALSE, data.save = FALSE, env, ...) {
  # This is just an alias

  set_defaults("fixest2_estimation")
  call_env_bis <- new.env(parent = parent.frame())

  res <- try(feNmlm(
    fml = fml, data = data, family = family, fixef = fixef,
    fixef.rm = fixef.rm, offset = offset, subset = subset,
    split = split, fsplit = fsplit,
    split.keep = split.keep, split.drop = split.drop,
    vcov = vcov, cluster = cluster,
    se = se, ssc = ssc, panel.id = panel.id, start = start,
    fixef.tol = fixef.tol, fixef.iter = fixef.iter, nthreads = nthreads,
    lean = lean, verbose = verbose, warn = warn, notes = notes,
    theta.init = theta.init, combine.quick = combine.quick,
    mem.clean = mem.clean, origin = "femlm", mc_origin_bis = match.call(),
    call_env_bis = call_env_bis, only.env = only.env,
    only.coef = only.coef, data.save = data.save,
    env = env, ...
  ), silent = TRUE)

  if ("try-error" %in% class(res)) {
    stop(format_error_msg(res, "femlm"))
  }

  return(res)
}
