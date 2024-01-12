#' Fixed effects nonlinear maximum likelihood models
#'
#' This function estimates maximum likelihood models (e.g., Poisson or Logit) with non-linear in parameters right-hand-sides and is efficient to handle any number of fixed effects. If you do not use non-linear in parameters right-hand-side, use [`femlm`] or [`feglm`] instead (their design is simpler).
#'
#' @inheritParams summary.fixest
#' @inheritParams panel
#' @inheritSection feols Lagging variables
#' @inheritSection feols Interactions
#' @inheritSection feols On standard-errors
#' @inheritSection feols Multiple estimations
#' @inheritSection feols Argument sliding
#' @inheritSection feols Piping
#' @inheritSection feols Tricks to estimate multiple LHS
#' @inheritSection xpd Dot square bracket operator in formulas
#'
#' @param fml A formula. This formula gives the linear formula to be estimated (it is similar to a `lm` formula), for example: `fml = z~x+y`. To include fixed-effects variables, insert them in this formula using a pipe (e.g. `fml = z~x+y|fixef_1+fixef_2`). To include a non-linear in parameters element, you must use the argment `NL.fml`. Multiple estimations can be performed at once: for multiple dep. vars, wrap them in `c()`: ex `c(y1, y2)`. For multiple indep. vars, use the stepwise functions: ex `x1 + csw(x2, x3)`. This leads to 6 estimation `fml = c(y1, y2) ~ x1 + cw0(x2, x3)`. See details. Square brackets starting with a dot can be used to call global variables: `y.[i] ~ x.[1:2]` will lead to `y3 ~ x1 + x2` if `i` is equal to 3 in the current environment (see details in [`xpd`]).
#' @param start Starting values for the coefficients in the linear part (for the non-linear part, use NL.start). Can be: i) a numeric of length 1 (e.g. `start = 0`, the default), ii) a numeric vector of the exact same length as the number of variables, or iii) a named vector of any length (the names will be used to initialize the appropriate coefficients).
#' @param NL.fml A formula. If provided, this formula represents the non-linear part of the right hand side (RHS). Note that contrary to the `fml` argument, the coefficients must explicitly appear in this formula. For instance, it can be `~a*log(b*x + c*x^3)`, where `a`, `b`, and `c` are the coefficients to be estimated. Note that only the RHS of the formula is to be provided, and NOT the left hand side.
#' @param split A one sided formula representing a variable (eg `split = ~var`) or a vector. If provided, the sample is split according to the variable and one estimation is performed for each value of that variable. If you also want to include the estimation for the full sample, use the argument `fsplit` instead. You can use the special operators `%keep%` and `%drop%` to select only a subset of values for which to split the sample. E.g. `split = ~var %keep% c("v1", "v2")` will split the sample only according to the values `v1` and `v2` of the variable `var`; it is equivalent to supplying the argument `split.keep = c("v1", "v2")`. By default there is partial matching on each value, you can trigger a regular expression evaluation by adding a `'@'` first, as in: `~var %drop% "@^v[12]"` which will drop values starting with `"v1"` or `"v2"` (of course you need to know regexes!).
#' @param fsplit A one sided formula representing a variable (eg `split = ~var`) or a vector. If provided, the sample is split according to the variable and one estimation is performed for each value of that variable. This argument is the same as split but also includes the full sample as the first estimation. You can use the special operators `%keep%` and `%drop%` to select only a subset of values for which to split the sample. E.g. `split = ~var %keep% c("v1", "v2")` will split the sample only according to the values `v1` and `v2` of the variable `var`; it is equivalent to supplying the argument `split.keep = c("v1", "v2")`. By default there is partial matching on each value, you can trigger a regular expression evaluation by adding an `'@'` first, as in: `~var %drop% "@^v[12]"` which will drop values starting with `"v1"` or `"v2"` (of course you need to know regexes!).
#' @param split.keep A character vector. Only used when `split`, or `fsplit`, is supplied. If provided, then the sample will be split only on the values of `split.keep`. The values in `split.keep` will be partially matched to the values of `split`. To enable regular expressions, you need to add an `'@'` first. For example `split.keep = c("v1", "@other|var")` will keep only the value in `split` partially matched by `"v1"` or the values containing `"other"` or `"var"`.
#' @param split.drop A character vector. Only used when `split`, or `fsplit`, is supplied. If provided, then the sample will be split only on the values that are not in `split.drop`. The values in `split.drop` will be partially matched to the values of `split`. To enable regular expressions, you need to add an `'@'` first. For example `split.drop = c("v1", "@other|var")` will drop only the value in `split` partially matched by `"v1"` or the values containing `"other"` or `"var"`.
#' @param data A data.frame containing the necessary variables to run the model. The variables of the non-linear right hand side of the formula are identified with this `data.frame` names. Can also be a matrix.
#' @param family Character scalar. It should provide the family. The possible values are "poisson" (Poisson model with log-link, the default), "negbin" (Negative Binomial model with log-link), "logit" (LOGIT model with log-link), "gaussian" (Gaussian model).
#' @param fixef Character vector. The names of variables to be used as fixed-effects. These variables should contain the identifier of each observation (e.g., think of it as a panel identifier). Note that the recommended way to include fixed-effects is to insert them directly in the formula.
#' @param subset A vector (logical or numeric) or a one-sided formula. If provided, then the estimation will be performed only on the observations defined by this argument.
#' @param NL.start (For NL models only) A list of starting values for the non-linear parameters. ALL the parameters are to be named and given a staring value. Example: `NL.start=list(a=1,b=5,c=0)`. Though, there is an exception: if all parameters are to be given the same starting value, you can use a numeric scalar.
#' @param lower (For NL models only) A list. The lower bound for each of the non-linear parameters that requires one. Example: `lower=list(b=0,c=0)`. Beware, if the estimated parameter is at his lower bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param upper (For NL models only) A list. The upper bound for each of the non-linear parameters that requires one. Example: `upper=list(a=10,c=50)`. Beware, if the estimated parameter is at his upper bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param NL.start.init (For NL models only) Numeric scalar. If the argument `NL.start` is not provided, or only partially filled (i.e. there remain non-linear parameters with no starting value), then the starting value of all remaining non-linear parameters is set to `NL.start.init`.
#' @param offset A formula or a numeric vector. An offset can be added to the estimation. If equal to a formula, it should be of the form (for example) `~0.5*x**2`. This offset is linearly added to the elements of the main formula 'fml'.
#' @param jacobian.method (For NL models only) Character scalar. Provides the method used to numerically compute the Jacobian of the non-linear part. Can be either `"simple"` or `"Richardson"`. Default is `"simple"`. See the help of [`jacobian`] for more information.
#' @param useHessian Logical. Should the Hessian be computed in the optimization stage? Default is `TRUE`.
#' @param hessian.args List of arguments to be passed to function [`genD`]. Defaults is missing. Only used with the presence of `NL.fml`.
#' @param opt.control List of elements to be passed to the optimization method [`nlminb`]. See the help page of [`nlminb`] for more information.
#' @param nthreads The number of threads. Can be: a) an integer lower than, or equal to, the maximum number of threads; b) 0: meaning all available threads will be used; c) a number strictly between 0 and 1 which represents the fraction of all threads to use. The default is to use 50% of all threads. You can set permanently the number of threads used within this package using the function [`setFixest_nthreads`].
#' @param verbose Integer, default is 0. It represents the level of information that should be reported during the optimisation process. If `verbose=0`: nothing is reported. If `verbose=1`: the value of the coefficients and the likelihood are reported. If `verbose=2`: `1` + information on the computing time of the null model, the fixed-effects coefficients and the hessian are reported.
#' @param theta.init Positive numeric scalar. The starting value of the dispersion parameter if `family="negbin"`. By default, the algorithm uses as a starting value the theta obtained from the model with only the intercept.
#' @param fixef.rm Can be equal to "perfect" (default), "singleton", "both" or "none". Controls which observations are to be removed. If "perfect", then observations having a fixed-effect with perfect fit (e.g. only 0 outcomes in Poisson estimations) will be removed. If "singleton", all observations for which a fixed-effect appears only once will be removed. The meaning of "both" and "none" is direct.
#' @param fixef.tol Precision used to obtain the fixed-effects. Defaults to `1e-5`. It corresponds to the maximum absolute difference allowed between two coefficients of successive iterations. Argument `fixef.tol` cannot be lower than `10000*.Machine$double.eps`. Note that this parameter is dynamically controlled by the algorithm.
#' @param fixef.iter Maximum number of iterations in fixed-effects algorithm (only in use for 2+ fixed-effects). Default is 10000.
#' @param deriv.iter Maximum number of iterations in the algorithm to obtain the derivative of the fixed-effects (only in use for 2+ fixed-effects). Default is 1000.
#' @param deriv.tol Precision used to obtain the fixed-effects derivatives. Defaults to `1e-4`. It corresponds to the maximum absolute difference allowed between two coefficients of successive iterations. Argument `deriv.tol` cannot be lower than `10000*.Machine$double.eps`.
#' @param warn Logical, default is `TRUE`. Whether warnings should be displayed (concerns warnings relating to convergence state).
#' @param notes Logical. By default, two notes are displayed: when NAs are removed (to show additional information) and when some observations are removed because of only 0 (or 0/1) outcomes in a fixed-effect setup (in Poisson/Neg. Bin./Logit models). To avoid displaying these messages, you can set `notes = FALSE`. You can remove these messages permanently by using `setFixest_notes(FALSE)`.
#' @param combine.quick Logical. When you combine different variables to transform them into a single fixed-effects you can do e.g. `y ~ x | paste(var1, var2)`. The algorithm provides a shorthand to do the same operation: `y ~ x | var1^var2`. Because pasting variables is a costly operation, the internal algorithm may use a numerical trick to hasten the process. The cost of doing so is that you lose the labels. If you are interested in getting the value of the fixed-effects coefficients after the estimation, you should use `combine.quick = FALSE`. By default it is equal to `FALSE` if the number of observations is lower than 50,000, and to `TRUE` otherwise.
#' @param only.env (Advanced users.) Logical, default is `FALSE`. If `TRUE`, then only the environment used to make the estimation is returned.
#' @param mem.clean Logical, default is `FALSE`. Only to be used if the data set is large compared to the available RAM. If `TRUE` then intermediary objects are removed as much as possible and [`gc`] is run before each substantial C++ section in the internal code to avoid memory issues.
#' @param lean Logical, default is `FALSE`. If `TRUE` then all large objects are removed from the returned result: this will save memory but will block the possibility to use many methods. It is recommended to use the arguments `se` or `cluster` to obtain the appropriate standard-errors at estimation time, since obtaining different SEs won't be possible afterwards.
#' @param env (Advanced users.) A `fixest` environment created by a `fixest` estimation with `only.env = TRUE`. Default is missing. If provided, the data from this environment will be used to perform the estimation.
#' @param only.coef Logical, default is `FALSE`. If `TRUE`, then only the estimated coefficients are returned. Note that the length of the vector returned is always the length of the number of coefficients to be estimated: this means that the variables found to be collinear are returned with an NA value.
#' @param ... Not currently used.
#'
#' @details
#' This function estimates maximum likelihood models where the conditional expectations are as follows:
#'
#' Gaussian likelihood:
#' \deqn{E(Y|X)=X\beta}{E(Y|X) = X*beta}
#' Poisson and Negative Binomial likelihoods:
#' \deqn{E(Y|X)=\exp(X\beta)}{E(Y|X) = exp(X*beta)}
#' where in the Negative Binomial there is the parameter \eqn{\theta}{theta} used to model the variance as \eqn{\mu+\mu^2/\theta}{mu+mu^2/theta}, with \eqn{\mu}{mu} the conditional expectation.
#' Logit likelihood:
#' \deqn{E(Y|X)=\frac{\exp(X\beta)}{1+\exp(X\beta)}}{E(Y|X) = exp(X*beta) / (1 + exp(X*beta))}
#'
#' When there are one or more fixed-effects, the conditional expectation can be written as:
#' \deqn{E(Y|X) = h(X\beta+\sum_{k}\sum_{m}\gamma_{m}^{k}\times C_{im}^{k}),}
#' where \eqn{h(.)} is the function corresponding to the likelihood function as shown before. \eqn{C^k} is the matrix associated to fixed-effect dimension \eqn{k} such that \eqn{C^k_{im}} is equal to 1 if observation \eqn{i} is of category \eqn{m} in the fixed-effect dimension \eqn{k} and 0 otherwise.
#'
#' When there are non linear in parameters functions, we can schematically split the set of regressors in two:
#' \deqn{f(X,\beta)=X^1\beta^1 + g(X^2,\beta^2)}
#' with first a linear term and then a non linear part expressed by the function g. That is, we add a non-linear term to the linear terms (which are \eqn{X*beta} and the fixed-effects coefficients). It is always better (more efficient) to put into the argument `NL.fml` only the non-linear in parameter terms, and add all linear terms in the `fml` argument.
#'
#' To estimate only a non-linear formula without even the intercept, you must exclude the intercept from the linear formula by using, e.g., `fml = z~0`.
#'
#' The over-dispersion parameter of the Negative Binomial family, theta, is capped at 10,000. If theta reaches this high value, it means that there is no overdispersion.
#'
#' @return
#' A `fixest` object. Note that `fixest` objects contain many elements and most of them are for internal use, they are presented here only for information. To access them, it is safer to use the user-level methods (e.g. [`vcov.fixest`], [`resid.fixest`], etc) or functions (like for instance [`fitstat`] to access any fit statistic).
#' \item{coefficients}{The named vector of coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The loglikelihood.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{nobs}{The number of observations.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{call}{The call.}
#' \item{fml}{The linear formula of the call.}
#' \item{fml_all}{A list containing different parts of the formula. Always contain the linear formula. Then, if relevant: `fixef`: the fixed-effects; `NL`: the non linear part of the formula.}
#' \item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{message}{The convergence message from the optimization procedures.}
#' \item{sq.cor}{Squared correlation between the dependent variable and the expected predictor (i.e. fitted.values) obtained by the estimation.}
#' \item{hessian}{The Hessian of the parameters.}
#' \item{fitted.values}{The fitted values are the expected value of the dependent variable for the fitted model: that is \eqn{E(Y|X)}.}
#' \item{cov.iid}{The variance-covariance matrix of the parameters.}
#' \item{se}{The standard-error of the parameters.}
#' \item{scores}{The matrix of the scores (first derivative for each observation).}
#' \item{family}{The ML family that was used for the estimation.}
#' \item{residuals}{The difference between the dependent variable and the expected predictor.}
#' \item{sumFE}{The sum of the fixed-effects for each observation.}
#' \item{offset}{The offset formula.}
#' \item{NL.fml}{The nonlinear formula of the call.}
#' \item{bounds}{Whether the coefficients were upper or lower bounded. -- This can only be the case when a non-linear formula is included and the arguments 'lower' or 'upper' are provided.}
#' \item{isBounded}{The logical vector that gives for each coefficient whether it was bounded or not. This can only be the case when a non-linear formula is included and the arguments 'lower' or 'upper' are provided.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{obs_selection}{(When relevant.) List containing vectors of integers. It represents the sequential selection of observation vis a vis the original data set.}
#' \item{fixef_removed}{In the case there were fixed-effects and some observations were removed because of only 0/1 outcome within a fixed-effect, it gives the list (for each fixed-effect dimension) of the fixed-effect identifiers that were removed.}
#' \item{theta}{In the case of a negative binomial estimation: the overdispersion parameter.}
#'
#'  @seealso
#' See also [`summary.fixest`] to see the results with the appropriate standard-errors, [`fixef.fixest`] to extract the fixed-effects coefficients.
#'
#' And other estimation methods: [`feols`], [`femlm`], [`feglm`], [`fepois`][fixest2::feglm], [`fenegbin`][fixest2::femlm].
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
#' On the unconditionnal Negative Binomial model:
#'
#' Allison, Paul D and Waterman, Richard P, 2002, "Fixed-Effects Negative Binomial Regression Models", Sociological Methodology 32(1) pp. 247--265
#'
#' @examples
#'
#' # This section covers only non-linear in parameters examples
#' # For linear relationships: use femlm or feglm instead
#'
#' # Generating data for a simple example
#' set.seed(1)
#' n <- 100
#' x <- rnorm(n, 1, 5)**2
#' y <- rnorm(n, -1, 5)**2
#' z1 <- rpois(n, x * y) + rpois(n, 2)
#' base <- data.frame(x, y, z1)
#'
#' # Estimating a 'linear' relation:
#' est1_L <- femlm(z1 ~ log(x) + log(y), base)
#' # Estimating the same 'linear' relation using a 'non-linear' call
#' est1_NL <- feNmlm(z1 ~ 1, base, NL.fml = ~ a * log(x) + b * log(y), NL.start = list(a = 0, b = 0))
#'
#' # Now generating a non-linear relation (E(z2) = x + y + 1):
#' z2 <- rpois(n, x + y) + rpois(n, 1)
#' base$z2 <- z2
#'
#' # Estimation using this non-linear form
#' est2_NL <- feNmlm(z2 ~ 0, base,
#'   NL.fml = ~ log(a * x + b * y),
#'   NL.start = 2, lower = list(a = 0, b = 0)
#' )
#' # we can't estimate this relation linearily
#' # => closest we can do:
#' est2_L <- femlm(z2 ~ log(x) + log(y), base)
#'
#' @export
feNmlm <- function(fml, data, family = c("poisson", "negbin", "logit", "gaussian"), NL.fml, vcov,
                   fixef, fixef.rm = "perfect", NL.start, lower, upper, NL.start.init,
                   offset, subset, split, fsplit, split.keep, split.drop,
                   cluster, se, ssc, panel.id,
                   start = 0, jacobian.method = "simple", useHessian = TRUE,
                   hessian.args = NULL, opt.control = list(), nthreads = getFixest_nthreads(),
                   lean = FALSE, verbose = 0, theta.init, fixef.tol = 1e-5, fixef.iter = 10000,
                   deriv.tol = 1e-4, deriv.iter = 1000, warn = TRUE, notes = getFixest_notes(),
                   combine.quick, mem.clean = FALSE, only.env = FALSE, only.coef = FALSE, env, ...) {
  time_start <- proc.time()

  if (missing(env)) {
    set_defaults("fixest_estimation")
    call_env <- new.env(parent = parent.frame())

    env <- try(fixest_env(
      fml = fml, data = data, family = family, NL.fml = NL.fml,
      fixef = fixef, fixef.rm = fixef.rm, NL.start = NL.start,
      lower = lower, upper = upper, NL.start.init = NL.start.init,
      offset = offset, subset = subset, split = split, fsplit = fsplit,
      split.keep = split.keep, split.drop = split.drop,
      vcov = vcov, cluster = cluster, se = se, ssc = ssc,
      panel.id = panel.id, linear.start = start, jacobian.method = jacobian.method,
      useHessian = useHessian, opt.control = opt.control, nthreads = nthreads,
      lean = lean, verbose = verbose, theta.init = theta.init,
      fixef.tol = fixef.tol, fixef.iter = fixef.iter, deriv.iter = deriv.iter,
      warn = warn, notes = notes, combine.quick = combine.quick,
      mem.clean = mem.clean, only.coef = only.coef, mc_origin = match.call(),
      call_env = call_env, computeModel0 = TRUE, ...
    ), silent = TRUE)
  } else if ((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
    stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
  }

  check_arg(only.env, "logical scalar")
  if (only.env) {
    return(env)
  }

  if ("try-error" %in% class(env)) {
    mc <- match.call()
    origin <- ifelse(is.null(mc$origin), "feNmlm", mc$origin)
    stop(format_error_msg(env, origin))
  }

  verbose <- get("verbose", env)
  if (verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep = "")

  IN_MULTI <- get("IN_MULTI", env)
  is_multi_root <- get("is_multi_root", env)
  if (is_multi_root) {
    on.exit(release_multi_notes())
    assign("is_multi_root", FALSE, env)
  }

  # Split ----

  do_split <- get("do_split", env)
  if (do_split) {
    res <- multi_split(env, feNmlm)

    return(res)
  }

  # Multi fixef ----

  do_multi_fixef <- get("do_multi_fixef", env)
  if (do_multi_fixef) {
    res <- multi_fixef(env, feNmlm)

    return(res)
  }

  # Multi LHS and RHS ----

  do_multi_lhs <- get("do_multi_lhs", env)
  do_multi_rhs <- get("do_multi_rhs", env)
  if (do_multi_lhs || do_multi_rhs) {
    res <- multi_LHS_RHS(env, feNmlm)

    return(res)
  }

  # Regular estimation ----


  # Objects needed for optimization + misc
  start <- get("start", env)
  lower <- get("lower", env)
  upper <- get("upper", env)
  gradient <- get("gradient", env)
  hessian <- get("hessian", env)
  family <- get("family", env)
  isLinear <- get("isLinear", env)
  isNonLinear <- get("isNL", env)
  opt.control <- get("opt.control", env)
  only.coef <- get("only.coef", env)

  lhs <- get("lhs", env)

  family <- get("family", env)
  famFuns <- get("famFuns", env)
  params <- get("params", env)

  isFixef <- get("isFixef", env)
  onlyFixef <- !isLinear && !isNonLinear && isFixef

  #
  # Model 0 + theta init
  #

  theta.init <- get("theta.init", env)
  model0 <- get_model_null(env, theta.init)

  # For the negative binomial:
  if (family == "negbin") {
    theta.init <- get("theta.init", env)
    if (is.null(theta.init)) {
      theta.init <- model0$theta
    }

    params <- c(params, ".theta")
    start <- c(start, theta.init)
    names(start) <- params
    upper <- c(upper, 10000)
    lower <- c(lower, 1e-3)

    assign("params", params, env)
  }

  assign("model0", model0, env)

  # the result
  res <- get("res", env)

  # NO VARIABLE -- ONLY FIXED-EFFECTS
  if (onlyFixef) {
    if (family == "negbin") {
      stop("To estimate the negative binomial model, you need at least one variable. (The estimation of the model with only the fixed-effects is not implemented.)")
    }

    res <- femlm_only_clusters(env)
    res$onlyFixef <- TRUE

    return(res)
  }

  # warnings => to avoid accumulation, but should appear even if the user stops the algorithm
  on.exit(warn_fixef_iter(env, stack_multi = IN_MULTI))

  #
  # Maximizing the likelihood
  #

  opt <- try(stats::nlminb(start = start, objective = femlm_ll, env = env, lower = lower, upper = upper, gradient = gradient, hessian = hessian, control = opt.control), silent = TRUE)

  if ("try-error" %in% class(opt)) {
    # We return the coefficients (can be interesting for debugging)
    iter <- get("iter", env)
    origin <- get("origin", env)
    warning_msg <- paste0("[", origin, "] Optimization failed at iteration ", iter, ". Reason: ", gsub("^[^\n]+\n *(.+\n)", "\\1", opt))
    if (IN_MULTI) {
      stack_multi_notes(warning_msg)
      return(fixest_NA_results(env))
    } else if (!"coef_evaluated" %in% names(env)) {
      # big problem right from the start
      stop(warning_msg)
    } else {
      coef <- get("coef_evaluated", env)
      warning(warning_msg, " Last evaluated coefficients returned.", call. = FALSE)
      return(coef)
    }
  } else {
    convStatus <- TRUE
    warning_msg <- ""
    if (!opt$message %in% c("X-convergence (3)", "relative convergence (4)", "both X-convergence and relative convergence (5)")) {
      warning_msg <- " The optimization algorithm did not converge, the results are not reliable."
      convStatus <- FALSE
    }

    coef <- opt$par
  }

  if (only.coef) {
    if (convStatus == FALSE) {
      if (IN_MULTI) {
        stack_multi_notes(warning_msg)
      } else {
        warning(warning_msg)
      }
    }

    names(coef) <- params
    return(coef)
  }



  # The Hessian
  hessian <- femlm_hessian(coef, env = env)
  # we add the names of the non linear variables in the hessian
  if (isNonLinear || family == "negbin") {
    dimnames(hessian) <- list(params, params)
  }

  # we create the Hessian without the bounded parameters
  hessian_noBounded <- hessian

  # Handling the bounds
  if (!isNonLinear) {
    NL.fml <- NULL
    bounds <- NULL
    isBounded <- NULL
  } else {
    nonlinear.params <- get("nonlinear.params", env)
    # we report the bounds & if the estimated parameters are bounded
    upper_bound <- upper[nonlinear.params]
    lower_bound <- lower[nonlinear.params]

    # 1: are the estimated parameters at their bounds?
    coef_NL <- coef[nonlinear.params]
    isBounded <- rep(FALSE, length(params))
    isBounded[seq_along(coef_NL)] <- (coef_NL == lower_bound) | (coef_NL == upper_bound)

    # 2: we save the bounds
    upper_bound_small <- upper_bound[is.finite(upper_bound)]
    lower_bound_small <- lower_bound[is.finite(lower_bound)]
    bounds <- list()
    if (length(upper_bound_small) > 0) bounds$upper <- upper_bound_small
    if (length(lower_bound_small) > 0) bounds$lower <- lower_bound_small
    if (length(bounds) == 0) {
      bounds <- NULL
    }

    # 3: we update the Hessian (basically, we drop the bounded element)
    if (any(isBounded)) {
      hessian_noBounded <- hessian[-which(isBounded), -which(isBounded), drop = FALSE]

      boundText <- ifelse(coef_NL == upper_bound, "Upper bounded", "Lower bounded")[isBounded]

      attr(isBounded, "type") <- boundText
    }
  }

  # Variance

  var <- NULL
  try(var <- solve(hessian_noBounded), silent = TRUE)
  if (is.null(var)) {
    warning_msg <- paste(warning_msg, "The information matrix is singular: presence of collinearity.")
    var <- hessian_noBounded * NA
    se <- diag(var)
  } else {
    se <- diag(var)
    se[se < 0] <- NA
    se <- sqrt(se)
  }

  # Warning message
  if (nchar(warning_msg) > 0) {
    if (warn) {
      if (IN_MULTI) {
        stack_multi_notes(warning_msg)
      } else {
        warning("[femlm]:", warning_msg, call. = FALSE)
        options("fixest_last_warning" = proc.time())
      }
    }
  }

  # To handle the bounded coefficient, we set its SE to NA
  if (any(isBounded)) {
    se <- se[params]
    names(se) <- params
  }

  zvalue <- coef / se
  pvalue <- 2 * pnorm(-abs(zvalue))

  # We add the information on the bound for the se & update the var to drop the bounded vars
  se_format <- se
  if (any(isBounded)) {
    se_format[!isBounded] <- decimalFormat(se_format[!isBounded])
    se_format[isBounded] <- boundText
  }

  coeftable <- data.frame("Estimate" = coef, "Std. Error" = se_format, "z value" = zvalue, "Pr(>|z|)" = pvalue, stringsAsFactors = FALSE)
  names(coeftable) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  row.names(coeftable) <- params

  attr(se, "type") <- attr(coeftable, "type") <- "IID"

  mu_both <- get_mu(coef, env, final = TRUE)
  mu <- mu_both$mu
  exp_mu <- mu_both$exp_mu

  # calcul pseudo r2
  loglik <- -opt$objective # moins car la fonction minimise
  ll_null <- model0$loglik

  # dummies are constrained, they don't have full dof (cause you need to take one value off for unicity)
  # this is an approximation, in some cases there can be more than one ref. But good approx.
  nparams <- res$nparams
  pseudo_r2 <- 1 - (loglik - nparams + 1) / ll_null

  # Calcul residus
  expected.predictor <- famFuns$expected.predictor(mu, exp_mu, env)
  residuals <- lhs - expected.predictor

  # calcul squared corr
  if (cpp_isConstant(expected.predictor)) {
    sq.cor <- NA
  } else {
    sq.cor <- stats::cor(lhs, expected.predictor)**2
  }
  ssr_null <- cpp_ssr_null(lhs)

  # The scores
  scores <- femlm_scores(coef, env)
  if (isNonLinear) {
    # we add the names of the non linear params in the score
    colnames(scores) <- params
  }

  n <- length(lhs)

  # Saving
  res$coefficients <- coef
  res$coeftable <- coeftable
  res$loglik <- loglik
  res$iterations <- opt$iterations
  res$ll_null <- ll_null
  res$ssr_null <- ssr_null
  res$pseudo_r2 <- pseudo_r2
  res$message <- opt$message
  res$convStatus <- convStatus
  res$sq.cor <- sq.cor
  res$fitted.values <- expected.predictor
  res$hessian <- hessian

  dimnames(var) <- list(params, params)
  res$cov.iid <- var
  # for compatibility with conleyreg
  res$cov.unscaled <- res$cov.iid

  res$se <- se
  res$scores <- scores
  res$family <- family
  res$residuals <- residuals

  # The value of mu (if cannot be recovered from fitted())
  if (family == "logit") {
    qui_01 <- expected.predictor %in% c(0, 1)
    if (any(qui_01)) {
      res$mu <- mu
    }
  } else if (family %in% c("poisson", "negbin")) {
    qui_0 <- expected.predictor == 0
    if (any(qui_0)) {
      res$mu <- mu
    }
  }


  if (!is.null(bounds)) {
    res$bounds <- bounds
    res$isBounded <- isBounded
  }

  # Fixed-effects
  if (isFixef) {
    useExp_fixefCoef <- family %in% c("poisson")

    sumFE <- attr(mu, "sumFE")
    if (useExp_fixefCoef) {
      sumFE <- rpar_log(sumFE, env)
    }

    res$sumFE <- sumFE

    # The LL and SSR with FE only
    if ("ll_fe_only" %in% names(env)) {
      res$ll_fe_only <- get("ll_fe_only", env)
      res$ssr_fe_only <- get("ssr_fe_only", env)
    } else {
      # we need to compute it

      # indicator of whether we compute the exp(mu)
      useExp <- family %in% c("poisson", "logit", "negbin")

      # mu, using the offset
      if (!is.null(res$offset)) {
        mu_noDum <- res$offset
      } else {
        mu_noDum <- 0
      }

      if (length(mu_noDum) == 1) mu_noDum <- rep(mu_noDum, n)

      exp_mu_noDum <- NULL
      if (useExp_fixefCoef) {
        exp_mu_noDum <- rpar_exp(mu_noDum, env)
      }

      assign("fixef.tol", 1e-4, env) # no need of supa precision
      dummies <- getDummies(mu_noDum, exp_mu_noDum, env, coef)

      exp_mu <- NULL
      if (useExp_fixefCoef) {
        # despite being called mu, it is in fact exp(mu)!!!
        exp_mu <- exp_mu_noDum * dummies
        mu <- rpar_log(exp_mu, env)
      } else {
        mu <- mu_noDum + dummies
        if (useExp) {
          exp_mu <- rpar_exp(mu, env)
        }
      }

      res$ll_fe_only <- famFuns$ll(lhs, mu, exp_mu, env, coef)
      ep <- famFuns$expected.predictor(mu, exp_mu, env)
      res$ssr_fe_only <- cpp_ssq(lhs - ep)
    }
  }

  if (family == "negbin") {
    theta <- coef[".theta"]
    res$theta <- theta

    if (notes && theta > 1000) {
      msg <- paste0("Very high value of theta (", theta, "). There is no sign of overdispersion, you may consider a Poisson model.")
      if (IN_MULTI) {
        stack_multi_notes(msg)
      } else {
        message(msg)
      }
    }
  }

  class(res) <- "fixest"

  if (verbose > 0) {
    cat("\n")
  }

  do_summary <- get("do_summary", env)
  if (do_summary) {
    vcov <- get("vcov", env)
    lean <- get("lean", env)
    ssc <- get("ssc", env)
    agg <- get("agg", env)
    summary_flags <- get("summary_flags", env)

    # To compute the RMSE and lean = TRUE
    if (lean) res$ssr <- cpp_ssq(res$residuals)

    res <- summary(res, vcov = vcov, ssc = ssc, agg = agg, lean = lean, summary_flags = summary_flags)
  }

  return(res)
}
