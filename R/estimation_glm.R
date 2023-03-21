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
#' See also [`summary.fixest`] to see the results with the appropriate standard-errors, [`fixef.fixest`] to extract the fixed-effects coefficients, and the function [`etable`] to visualize the results of multiple estimations.
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
#'
#' # Poisson estimation
#' res = feglm(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris, "poisson")
#'
#' # You could also use fepois
#' res_pois = fepois(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris)
#'
#' # With the fit method:
#' res_fit = feglm.fit(iris$Sepal.Length, iris[, 2:3], iris$Species, "poisson")
#'
#' # All results are identical:
#' etable(res, res_pois, res_fit)
#'
#' # Note that you have many more examples in feols
#'
#' #
#' # Multiple estimations:
#' #
#'
#' # 6 estimations
#' est_mult = fepois(c(Ozone, Solar.R) ~ Wind + Temp + csw0(Wind:Temp, Day), airquality)
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
#' est_split = fepois(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
#'                   airquality, split = ~ Month)
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
feglm = function(fml, data, family = "gaussian", vcov, offset, weights, subset, split,
                 fsplit, split.keep, split.drop, cluster, se, ssc, panel.id, start = NULL,
                 etastart = NULL, mustart = NULL, fixef, fixef.rm = "perfect",
                 fixef.tol = 1e-6, fixef.iter = 10000, collin.tol = 1e-10,
                 glm.iter = 25, glm.tol = 1e-8, nthreads = getFixest_nthreads(),
                 lean = FALSE, warn = TRUE, notes = getFixest_notes(), verbose = 0,
                 only.coef = FALSE,
                 combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

    if(missing(weights)) weights = NULL

    time_start = proc.time()

    if(missing(env)){
        set_defaults("fixest_estimation")
        call_env = new.env(parent = parent.frame())

        env = try(fixest_env(fml=fml, data=data, family = family,
                             offset = offset, weights = weights, subset = subset,
                             split = split, fsplit = fsplit,
                             split.keep = split.keep, split.drop = split.drop,
                             vcov = vcov,
                             cluster = cluster, se = se, ssc = ssc,
                             panel.id = panel.id, linear.start = start,
                             etastart=etastart, mustart = mustart, fixef = fixef,
                             fixef.rm = fixef.rm, fixef.tol = fixef.tol,
                             fixef.iter=fixef.iter, collin.tol = collin.tol,
                             glm.iter = glm.iter, glm.tol = glm.tol,
                             nthreads = nthreads, lean = lean, warn = warn, notes = notes,
                             verbose = verbose, combine.quick = combine.quick,
                             mem.clean = mem.clean, only.coef = only.coef, origin = "feglm",
                             mc_origin = match.call(), call_env = call_env, ...), silent = TRUE)

    } else if((r <- !is.environment(env)) || !isTRUE(env$fixest_env)){
        stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
    }

    if("try-error" %in% class(env)){
        mc = match.call()
        origin = ifelse(is.null(mc$origin), "feglm", mc$origin)
        stop(format_error_msg(env, origin))
    }

    check_arg(only.env, "logical scalar")
    if(only.env){
        return(env)
    }

    verbose = get("verbose", env)
    if(verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep="")

    # workhorse is feglm.fit (OK if error msg leads to feglm.fit [clear enough])
    res = feglm.fit(env = env)

    res
}



#' @rdname feglm
#' @export
feglm.fit = function(y, X, fixef_df, family = "gaussian", vcov, offset, split,
                     fsplit, split.keep, split.drop, cluster, se, ssc, weights, subset, start = NULL,
                     etastart = NULL, mustart = NULL, fixef.rm = "perfect",
                     fixef.tol = 1e-6, fixef.iter = 10000, collin.tol = 1e-10,
                     glm.iter = 25, glm.tol = 1e-8, nthreads = getFixest_nthreads(),
                     lean = FALSE, warn = TRUE, notes = getFixest_notes(), mem.clean = FALSE,
                     verbose = 0, only.env = FALSE, only.coef = FALSE, env, ...){

    dots = list(...)

    lean_internal = isTRUE(dots$lean_internal)
    means = 1
    if(!missing(env)){
        # This is an internal call from the function feglm
        # no need to further check the arguments
        # we extract them from the env

        if((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
            stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
        }

        # main variables
        if(missing(y)) y = get("lhs", env)
        if(missing(X)) X = get("linear.mat", env)
        if(!missing(fixef_df) && is.null(fixef_df)){
            assign("isFixef", FALSE, env)
        }

        if(missing(offset)) offset = get("offset.value", env)
        if(missing(weights)) weights = get("weights.value", env)

        # other params
        if(missing(fixef.tol)) fixef.tol = get("fixef.tol", env)
        if(missing(fixef.iter)) fixef.iter = get("fixef.iter", env)
        if(missing(collin.tol)) collin.tol = get("collin.tol", env)
        if(missing(glm.iter)) glm.iter = get("glm.iter", env)
        if(missing(glm.tol)) glm.tol = get("glm.tol", env)
        if(missing(warn)) warn = get("warn", env)
        if(missing(verbose)) verbose = get("verbose", env)
        if(missing(only.coef)) only.coef = get("only.coef", env)

        # starting point of the fixed-effects
        if(!is.null(dots$means)) means = dots$means

        # init
        init.type = get("init.type", env)
        starting_values = get("starting_values", env)
        if(lean_internal){
            # Call within here => either null model or fe only
            init.type = "default"
            if(!is.null(etastart)){
                init.type = "eta"
                starting_values = etastart
            }
        }

    } else {

        if(missing(weights)) weights = NULL

        time_start = proc.time()

        set_defaults("fixest_estimation")
        call_env = new.env(parent = parent.frame())

        env = try(fixest_env(y = y, X = X, fixef_df = fixef_df, family = family,
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
                             mc_origin = match.call(), call_env = call_env, ...), silent = TRUE)

        if("try-error" %in% class(env)){
            stop(format_error_msg(env, "feglm.fit"))
        }

        check_arg(only.env, "logical scalar")
        if(only.env){
            return(env)
        }

        verbose = get("verbose", env)
        if(verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep="")

        # y/X
        y = get("lhs", env)
        X = get("linear.mat", env)

        # offset
        offset = get("offset.value", env)

        # weights
        weights = get("weights.value", env)

        # warn flag
        warn = get("warn", env)

        # init
        init.type = get("init.type", env)
        starting_values = get("starting_values", env)
    }

    IN_MULTI = get("IN_MULTI", env)
    is_multi_root = get("is_multi_root", env)
    if(is_multi_root){
        on.exit(release_multi_notes())
        assign("is_multi_root", FALSE, env)
    }

    # Split ----

    do_split = get("do_split", env)
    if(do_split){

        res = multi_split(env, feglm.fit)

        return(res)
    }

    # Multi fixef ----

    do_multi_fixef = get("do_multi_fixef", env)
    if(do_multi_fixef){

        res = multi_fixef(env, feglm.fit)

        return(res)
    }

    # Multi LHS and RHS ----

    do_multi_lhs = get("do_multi_lhs", env)
    do_multi_rhs = get("do_multi_rhs", env)
    if(do_multi_lhs || do_multi_rhs){

        res = multi_LHS_RHS(env, feglm.fit)

        return(res)
    }

    # Regular estimation -----

    # Setup:
    family = get("family_funs", env)
    isFixef = get("isFixef", env)
    nthreads = get("nthreads", env)
    isWeight = length(weights) > 1
    isOffset = length(offset) > 1
    nobs = length(y)
    onlyFixef = length(X) == 1

    # the preformatted results
    res = get("res", env)

    #
    # glm functions:
    #

    variance = family$variance
    linkfun = family$linkfun
    linkinv = family$linkinv

    dev.resids = family$dev.resids
    sum_dev.resids = function(y, mu, eta, wt) sum(dev.resids(y, mu, wt))

    fun_mu.eta = family$mu.eta
    mu.eta = function(mu, eta) fun_mu.eta(eta)

    valideta = family$valideta
    if(is.null(valideta)) valideta = function(...) TRUE

    validmu = family$validmu
    if(is.null(validmu)) validmu = function(mu) TRUE

    family_equiv = family$family_equiv

    # Optimizing poisson/logit
    if(family_equiv == "poisson"){
        linkfun = function(mu) cpppar_log(mu, nthreads)
        linkinv = function(eta) cpppar_poisson_linkinv(eta, nthreads)

        y_pos = y[y > 0]
        qui_pos = y > 0
        if(isWeight){
            constant = sum(weights[qui_pos] * y_pos * cpppar_log(y_pos, nthreads) - weights[qui_pos] * y_pos)
            sum_dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(wt[qui_pos] * y_pos * eta[qui_pos]) + sum(wt * mu))
        } else {
            constant = sum(y_pos * cpppar_log(y_pos, nthreads) - y_pos)
            sum_dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(y_pos * eta[qui_pos]) + sum(mu))
        }

        mu.eta = function(mu, eta) mu
        validmu = function(mu) cpppar_poisson_validmu(mu, nthreads)

    } else if(family_equiv == "logit"){
        # To avoid annoying warnings
        family$initialize = quasibinomial()$initialize

        linkfun = function(mu) cpppar_logit_linkfun(mu, nthreads)
        linkinv = function(eta) cpppar_logit_linkinv(eta, nthreads)
        if(isWeight){
            sum_dev.resids = function(y, mu, eta, wt) sum(cpppar_logit_devresids(y, mu, wt, nthreads))
        } else {
            sum_dev.resids = function(y, mu, eta, wt) sum(cpppar_logit_devresids(y, mu, 1, nthreads))
        }
        sum_dev.resids = sum_dev.resids

        mu.eta = function(mu, eta) cpppar_logit_mueta(eta, nthreads)
    }

    #
    # Init
    #

    if(mem.clean){
        gc()
    }

    if(init.type == "mu"){

        mu = starting_values

        if(!valideta(mu)){
            stop("In 'mustart' the values provided are not valid.")
        }

        eta = linkfun(mu)

    } else if(init.type == "eta"){

        eta = starting_values

        if(!valideta(eta)){
            stop("In 'etastart' the values provided are not valid.")
        }

        mu = linkinv(eta)

    } else if(init.type == "coef"){
        # If there are fixed-effects we MUST first compute the FE model with starting values as offset
        #   otherwise we are too far away from the solution and starting values may lead to divergence
        #   (hence step halving would be required)
        # This means that initializing with coefficients incurs large computational costs
        #   with fixed-effects

        start = get("start", env)
        offset_fe = offset + cpppar_xbeta(X, start, nthreads)

        if(isFixef){
            mustart = 0
            eval(family$initialize)
            eta = linkfun(mustart)

            # just a rough estimate (=> high tol values) [no benefit in high precision]
            model_fe = try(feglm.fit(X = 0, etastart = eta, offset = offset_fe, glm.tol = 1e-2, fixef.tol = 1e-2, env = env, lean_internal = TRUE))

            if("try-error" %in% class(model_fe)){
                stop("Estimation failed during initialization when getting the fixed-effects, maybe change the values of 'start'? \n", model_fe)
            }

            eta = model_fe$linear.predictors
            mu = model_fe$fitted.values
            devold = model_fe$deviance
        } else {
            eta = offset_fe
            mu = linkinv(eta)
            devold = sum_dev.resids(y, mu, eta, wt = weights)
        }

        wols_old = list(fitted.values = eta - offset)

    } else {
        mustart = 0
        eval(family$initialize)
        eta = linkfun(mustart)
        mu = linkinv(eta)

        # NOTA: FE only => ADDS LOTS OF COMPUTATIONAL COSTS without convergence benefit
    }

    if(mem.clean){
        gc()
    }

    if(init.type != "coef"){
        # starting deviance with constant equal to 1e-5
        # this is important for getting in step halving early (when deviance goes awry right from the start)
        devold = sum_dev.resids(y, rep(linkinv(1e-5), nobs), rep(1e-5, nobs), wt = weights)
        wols_old = list(fitted.values = rep(1e-5, nobs))
    }

    if(!validmu(mu) || !valideta(eta)){
        stop("Current starting values are not valid.")
    }

    assign("nb_sh", 0, env)
    on.exit(warn_step_halving(env, stack_multi = IN_MULTI))

    if((init.type == "coef" && verbose >= 1) || verbose >= 4) {
        cat("Deviance at initializat.  = ", numberFormatNormal(devold), "\n", sep = "")
    }

    #
    # The main loop
    #

    wols_means = 1
    conv = FALSE
    warning_msg = div_message = ""
    for (iter in 1:glm.iter) {

        if(mem.clean){
            gc()
        }

        mu.eta.val = mu.eta(mu, eta)
        var_mu = variance(mu)

        # controls
        any_pblm_mu = cpp_any_na_null(var_mu)
        if(any_pblm_mu){
            if (anyNA(var_mu)){
                stop("NAs in V(mu), at iteration ", iter, ".")
            } else if (any(var_mu == 0)){
                stop("0s in V(mu), at iteration ", iter, ".")
            }
        }

        if(anyNA(mu.eta.val)){
            stop("NAs in d(mu)/d(eta), at iteration ", iter, ".")
        }

        if(isOffset){
            z = (eta - offset) + (y - mu)/mu.eta.val
        } else {
            z = eta + (y - mu)/mu.eta.val
        }

        w = as.vector(weights * mu.eta.val**2 / var_mu)

        is_0w = w == 0
        any_0w = any(is_0w)
        if(any_0w && all(is_0w)){
            warning_msg = paste0("No informative observation at iteration ", iter, ".")
            div_message = "No informative observation."
            break
        }

        if(mem.clean && iter > 1){
            rm(wols)
            gc()
        }

        wols = feols(y = z, X = X, weights = w, means = wols_means,
                     correct_0w = any_0w, env = env, fixef.tol = fixef.tol,
                     fixef.iter = fixef.iter, collin.tol = collin.tol, nthreads = nthreads,
                     mem.clean = mem.clean, warn = warn,
                     verbose = FALSE, fromGLM = TRUE)

        if(isTRUE(wols$NA_model)){
            return(wols)
        }

        # In theory OLS estimation is guaranteed to exist
        # yet, NA coef may happen with non-infinite very large values of z/w (e.g. values > 1e100)
        if(anyNA(wols$coefficients)){
            if(iter == 1){
                stop("Weighted-OLS returns NA coefficients at first iteration, step halving cannot be performed. Try other starting values?")
            }

            warning_msg = paste0("Divergence at iteration ", iter, ": Weighted-OLS returns NA coefficients. Last evaluated coefficients with finite deviance are returned for information purposes.")
            div_message = "Weighted-OLS returned NA coefficients."
            wols = wols_old
            break
        } else {
            wols_means = wols$means
        }

        eta = wols$fitted.values

        if(isOffset){
            eta = eta + offset
        }

        if(mem.clean){
            gc()
        }

        mu = linkinv(eta)

        dev = sum_dev.resids(y, mu, eta, wt = weights)
        dev_evol = dev - devold

        if(verbose >= 1) cat("Iteration: ", sprintf("%02i", iter), " -- Deviance = ", numberFormatNormal(dev), " -- Evol. = ", dev_evol, "\n", sep = "")

        #
        # STEP HALVING
        #

        no_SH = is.finite(dev) && (abs(dev_evol) < glm.tol || abs(dev_evol)/(0.1 + abs(dev)) < glm.tol)
        if(no_SH == FALSE && (!is.finite(dev) || dev_evol > 0 || !valideta(eta) || !validmu(mu))){

            if(!is.finite(dev)){
                # we report step-halving but only for non-finite deviances
                # other situations are OK (it just happens)
                nb_sh = get("nb_sh", env)
                assign("nb_sh", nb_sh + 1, env)
            }

            eta_new = wols$fitted.values
            eta_old = wols_old$fitted.values

            iter_sh = 0
            do_exit = FALSE
            while(!is.finite(dev) || dev_evol > 0 || !valideta(eta_new) || !validmu(mu)){

                if(iter == 1 && (is.finite(dev) && valideta(eta_new) && validmu(mu)) && iter_sh >= 2){
                    # BEWARE FIRST ITERATION:
                    # at first iteration, the deviance can be higher than the init, and SH may not help
                    # we need to make sure we get out of SH before it's messed up
                    break
                } else if(iter_sh == glm.iter){

                    # if first iteration => means algo did not find viable solution
                    if(iter == 1){
                        stop("Algorithm failed at first iteration. Step-halving could not find a valid set of parameters.")
                    }

                    # Problem only if the deviance is non-finite or eta/mu not valid
                    # Otherwise, it means that we're at a maximum

                    if(!is.finite(dev) || !valideta(eta_new) || !validmu(mu)){
                        # message
                        msg = ifelse(!is.finite(dev), "non-finite deviance", "no valid eta/mu")

                        warning_msg = paste0("Divergence at iteration ", iter, ": ", msg, ". Step halving: no valid correction found. Last evaluated coefficients with finite deviance are returned for information purposes.")
                        div_message = paste0(msg, " despite step-halving")
                        wols = wols_old
                        do_exit = TRUE
                    }

                    break
                }

                iter_sh = iter_sh + 1
                eta_new = (eta_old + eta_new) / 2

                if(mem.clean){
                    gc()
                }

                mu = linkinv(eta_new + offset)
                dev = sum_dev.resids(y, mu, eta_new + offset, wt = weights)
                dev_evol = dev - devold

                if(verbose >= 3) cat("Step-halving: iter =", iter_sh, "-- dev:", numberFormatNormal(dev), "-- evol:", numberFormatNormal(dev_evol), "\n")
            }

            if(do_exit) break

            # it worked: update
            eta = eta_new + offset
            wols$fitted.values = eta_new
            # NOTA: we must NOT end with a step halving => we need a proper weighted-ols estimation
            # we force the algorithm to continue
            dev_evol = Inf

            if(verbose >= 2){
                cat("Step-halving: new deviance = ", numberFormatNormal(dev), "\n", sep = "")
            }

        }

        if(abs(dev_evol)/(0.1 + abs(dev)) < glm.tol){
            conv = TRUE
            break
        } else {
            devold = dev
            wols_old = wols
        }
    }

    # Convergence flag
    if(!conv){
        if(iter == glm.iter){
            warning_msg = paste0("Absence of convergence: Maximum number of iterations reached (", glm.iter, "). Final deviance: ", numberFormatNormal(dev), ".")
            div_message = "no convergence: Maximum number of iterations reached"
        }

        res$convStatus = FALSE
        res$message = div_message
    } else {
        res$convStatus = TRUE
    }

    #
    # post processing
    #

    if(only.coef){
        if(wols$multicol){
            beta = rep(NA_real_, length(wols$is_excluded))
            beta[!wols$is_excluded] = wols$coefficients
        } else {
            beta = wols$coefficients
        }

        names(beta) = colnames(X)
        return(beta)
    }

    # Collinearity message
    collin.adj = 0
    if(wols$multicol){
        var_collinear = colnames(X)[wols$is_excluded]
        if(notes){
            msg = dsb("w!The variable.[*s_, q, ' has'V, 3KO, C?var_collinear] been
                       removed because of collinearity (see $collin.var).")
            if(IN_MULTI){
                stack_multi_notes(msg)
            } else {
                message(msg)
            }
        }

        res$collin.var = var_collinear

        # full set of coeffficients with NAs
        collin.coef = setNames(rep(NA, ncol(X)), colnames(X))
        collin.coef[!wols$is_excluded] = wols$coefficients
        res$collin.coef = collin.coef

        wols$X_demean = wols$X_demean[, !wols$is_excluded, drop = FALSE]
        X = X[, !wols$is_excluded, drop = FALSE]

        collin.adj = sum(wols$is_excluded)
    }

    res$nparams = res$nparams - collin.adj

    res$irls_weights = w # weights from the iteratively reweighted least square

    res$coefficients = coef = wols$coefficients
    res$collin.min_norm = wols$collin.min_norm

    if(!is.null(wols$warn_varying_slope)){
        msg = wols$warn_varying_slope
        if(IN_MULTI){
            stack_multi_notes(msg)
        } else {
            warning(msg)
        }
    }

    res$linear.predictors = wols$fitted.values
    if(isOffset){
        res$linear.predictors = res$linear.predictors + offset
    }

    res$fitted.values = linkinv(res$linear.predictors)
    res$residuals = y - res$fitted.values

    if(onlyFixef) res$onlyFixef = onlyFixef

    # dispersion + scores
    if(family$family %in% c("poisson", "binomial")){
        res$dispersion = 1
    } else {
        weighted_resids = wols$residuals * res$irls_weights
        # res$dispersion = sum(weighted_resids ** 2) / sum(res$irls_weights)
        # I use the second line to fit GLM's
        res$dispersion = sum(weighted_resids * wols$residuals) / max(res$nobs - res$nparams, 1)
    }

    res$working_residuals = wols$residuals

    if(!onlyFixef && !lean_internal){
        # score + hessian + vcov

        if(mem.clean){
            gc()
        }

        # dispersion + scores
        if(family$family %in% c("poisson", "binomial")){
            res$scores = (wols$residuals * res$irls_weights) * wols$X_demean
            res$hessian = cpppar_crossprod(wols$X_demean, res$irls_weights, nthreads)
        } else {
            res$scores = (weighted_resids / res$dispersion) * wols$X_demean
            res$hessian = cpppar_crossprod(wols$X_demean, res$irls_weights, nthreads) / res$dispersion
        }

        if(any(diag(res$hessian) < 0)){
            # This should not occur, but I prefer to be safe
            # In fact it's the opposite of the Hessian
            stop("Negative values in the diagonal of the Hessian found after the weighted-OLS stage. (If possible, could you send a replicable example to fixest's author? He's curious about when that actually happens, since in theory it should never happen.)")
        }

        # I put tol = 0, otherwise we may remove everything mistakenly
        # when VAR(Y) >>> VAR(X) // that is especially TRUE for missspecified
        # Poisson models
        info_inv = cpp_cholesky(res$hessian, tol = 0, nthreads = nthreads)

        if(!is.null(info_inv$all_removed)){
            # This should not occur, but I prefer to be safe
            stop("Not a single variable with a minimum of explanatory power found after the weighted-OLS stage. (If possible, could you send a replicable example to fixest's author? He's curious about when that actually happens, since in theory it should never happen.)")
        }

        var = info_inv$XtX_inv
        is_excluded = info_inv$id_excl

        if(any(is_excluded)){
            # There should be no remaining collinearity
            warning_msg = paste(warning_msg, "Residual collinearity was found after the weighted-OLS stage. The covariance is not defined. (This should not happen. If possible, could you send a replicable example to fixest's author? He's curious about when that actually happens.)")
            var = matrix(NA, length(is_excluded), length(is_excluded))
        }
        res$cov.iid = var

        rownames(res$cov.iid) = colnames(res$cov.iid) = names(coef)

        # for compatibility with conleyreg
        res$cov.unscaled = res$cov.iid

        # se
        se = diag(res$cov.iid)
        se[se < 0] = NA
        se = sqrt(se)

        # coeftable
        zvalue = coef/se
        use_t = !family$family %in% c("poisson", "binomial")
        if(use_t){
            pvalue = 2*pt(-abs(zvalue), max(res$nobs - res$nparams, 1))
            ctable_names = c("Estimate", "Std. Error", "t value",  "Pr(>|t|)")
        } else {
            pvalue = 2*pnorm(-abs(zvalue))
            ctable_names = c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
        }

        coeftable = data.frame("Estimate"=coef, "Std. Error"=se, "z value"=zvalue, "Pr(>|z|)"=pvalue)
        names(coeftable) = ctable_names
        row.names(coeftable) = names(coef)

        attr(se, "type") = attr(coeftable, "type") = "IID"
        res$coeftable = coeftable
        res$se = se
    }

    if(nchar(warning_msg) > 0){
        if(warn){
            if(IN_MULTI){
                stack_multi_notes(warning_msg)
            } else {
                warning(warning_msg, call. = FALSE)
                options("fixest_last_warning" = proc.time())
            }
        }
    }

    n = length(y)
    res$nobs = n

    df_k = res$nparams

    # r2s
    if(!cpp_isConstant(res$fitted.values)){
        res$sq.cor = tryCatch(stats::cor(y, res$fitted.values), warning = function(x) NA_real_)**2
    } else {
        res$sq.cor = NA
    }

    # deviance
    res$deviance = dev

    # simpler form for poisson
    if(family_equiv == "poisson"){
        if(isWeight){

            if(mem.clean){
                gc()
            }

            res$loglik = sum( (y * eta - mu - cpppar_lgamma(y + 1, nthreads)) * weights)
        } else {
            # lfact is later used in model0 and is costly to compute
            lfact = sum(rpar_lgamma(y + 1, env))
            assign("lfactorial", lfact, env)
            res$loglik = sum(y * eta - mu) - lfact
        }
    } else {
        res$loglik = family$aic(y = y, n = rep.int(1, n), mu = res$fitted.values, wt = weights, dev = dev) / -2
    }

    if(lean_internal){
        return(res)
    }

    # The pseudo_r2
    if(family_equiv %in% c("poisson", "logit")){
        model0 = get_model_null(env, theta.init = NULL)
        ll_null = model0$loglik
        fitted_null = linkinv(model0$constant)

    } else {
        if(verbose >= 1) cat("Null model:\n")

        if(mem.clean){
            gc()
        }

        model_null = feglm.fit(X = matrix(1, nrow = n, ncol = 1), fixef_df = NULL, env = env, lean_internal = TRUE)
        ll_null = model_null$loglik
        fitted_null = model_null$fitted.values
    }
    res$ll_null = ll_null
    res$pseudo_r2 = 1 - (res$loglik - df_k)/(ll_null - 1)

    # fixef info
    if(isFixef){
        if(onlyFixef){
            res$sumFE = res$linear.predictors
        } else {
            res$sumFE = res$linear.predictors - cpppar_xbeta(X, res$coefficients, nthreads)
        }

        if(isOffset){
            res$sumFE = res$sumFE - offset
        }

    }

    # other
    res$iterations = iter
    class(res) = "fixest"

    do_summary = get("do_summary", env)
    if(do_summary){
        vcov = get("vcov", env)
        lean = get("lean", env)
        ssc = get("ssc", env)
        agg = get("agg", env)
        summary_flags = get("summary_flags", env)

        # To compute the RMSE and lean = TRUE
        if(lean) res$ssr = cpp_ssq(res$residuals, weights)

        res = summary(res, vcov = vcov, agg = agg, ssc = ssc, lean = lean, summary_flags = summary_flags)
    }

    return(res)
}


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
#' @param fml A formula representing the relation to be estimated. For example: `fml = z~x+y`. To include fixed-effects, insert them in this formula using a pipe: e.g. `fml = z~x+y|fixef_1+fixef_2`. Multiple estimations can be performed at once: for multiple dep. vars, wrap them in `c()`: ex `c(y1, y2)`. For multiple indep. vars, use the stepwise functions: ex `x1 + csw(x2, x3)`. The formula `fml = c(y1, y2) ~ x1 + cw0(x2, x3)` leads to 6 estimation, see details. Square brackets starting with a dot can be used to call global variables: `y.[i] ~ x.[1:2]` will lead to `y3 ~ x1 + x2` if `i` is equal to 3 in the current environment (see details in [`xpd`]).
#' @param start Starting values for the coefficients. Can be: i) a numeric of length 1 (e.g. `start = 0`, the default), ii) a numeric vector of the exact same length as the number of variables, or iii) a named vector of any length (the names will be used to initialize the appropriate coefficients).
#'
#' @details
#' Note that the functions [`feglm`] and [`femlm`] provide the same results when using the same families but differ in that the latter is a direct maximum likelihood optimization (so the two can really have different convergence rates).
#'
#' @return
#' A `fixest` object. Note that `fixest` objects contain many elements and most of them are for internal use, they are presented here only for information. To access them, it is safer to use the user-level methods (e.g. [`vcov.fixest`], [`resid.fixest`], etc) or functions (like for instance [`fitstat`] to access any fit statistic).
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{fml_all}{A list containing different parts of the formula. Always contain the linear formula. Then, if relevant: `fixef`: the fixed-effects; `NL`: the non linear part of the formula.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{convStatus}{Logical, convergence status.}
#' \item{message}{The convergence message from the optimization procedures.}
#' \item{obs_selection}{(When relevant.) List containing vectors of integers. It represents the sequential selection of observation vis a vis the original data set.}
#' \item{fixef_removed}{(When relevant.) In the case there were fixed-effects and some observations were removed because of only 0/1 outcome within a fixed-effect, it gives the list (for each fixed-effect dimension) of the fixed-effect identifiers that were removed.}
#' \item{coefficients}{The named vector of estimated coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The log-likelihood.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
#' \item{ll_fe_only}{Log-likelihood of the model with only the fixed-effects.}
#' \item{ssr_null}{Sum of the squared residuals of the null model (containing only with the intercept).}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{fitted.values}{The fitted values are the expected value of the dependent variable for the fitted model: that is \eqn{E(Y|X)}.}
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
#'
#'
#' @seealso
#' See also [`summary.fixest`] to see the results with the appropriate standard-errors, [`fixef.fixest`] to extract the fixed-effects coefficients, and the function [`etable`] to visualize the results of multiple estimations.
#' And other estimation methods: [`feols`], [`feglm`], [`fepois`], [`feNmlm`].
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
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade => we account for 3 fixed-effects
#' # 1) Poisson estimation
#' est_pois = femlm(Euros ~ log(dist_km) | Origin + Destination + Product, trade)
#'
#' # 2) Log-Log Gaussian estimation (with same FEs)
#' est_gaus = update(est_pois, log(Euros+1) ~ ., family = "gaussian")
#'
#' # Comparison of the results using the function etable
#' etable(est_pois, est_gaus)
#' # Now using two way clustered standard-errors
#' etable(est_pois, est_gaus, se = "twoway")
#'
#' # Comparing different types of standard errors
#' sum_hetero   = summary(est_pois, se = "hetero")
#' sum_oneway   = summary(est_pois, se = "cluster")
#' sum_twoway   = summary(est_pois, se = "twoway")
#' sum_threeway = summary(est_pois, se = "threeway")
#'
#' etable(sum_hetero, sum_oneway, sum_twoway, sum_threeway)
#'
#'
#' #
#' # Multiple estimations:
#' #
#'
#' # 6 estimations
#' est_mult = femlm(c(Ozone, Solar.R) ~ Wind + Temp + csw0(Wind:Temp, Day), airquality)
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
#' est_split = fepois(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
#'                   airquality, split = ~ Month)
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
femlm = function(fml, data, family = c("poisson", "negbin", "logit", "gaussian"), vcov,
                 start = 0, fixef, fixef.rm = "perfect", offset, subset,
                 split, fsplit, split.keep, split.drop,
                 cluster, se, ssc, panel.id, fixef.tol = 1e-5, fixef.iter = 10000,
                 nthreads = getFixest_nthreads(), lean = FALSE, verbose = 0, warn = TRUE,
                 notes = getFixest_notes(), theta.init, combine.quick, mem.clean = FALSE,
                 only.env = FALSE, only.coef = FALSE, env, ...){

    # This is just an alias

    set_defaults("fixest_estimation")
    call_env_bis = new.env(parent = parent.frame())

    res = try(feNmlm(fml = fml, data = data, family = family, fixef = fixef,
                     fixef.rm = fixef.rm, offset = offset, subset = subset,
                     split = split, fsplit = fsplit,
                     split.keep = split.keep, split.drop = split.drop,
                     vcov = vcov, cluster = cluster,
                     se = se, ssc = ssc, panel.id = panel.id, start = start,
                     fixef.tol=fixef.tol, fixef.iter=fixef.iter, nthreads=nthreads,
                     lean = lean, verbose=verbose, warn=warn, notes=notes,
                     theta.init = theta.init, combine.quick = combine.quick,
                     mem.clean = mem.clean, origin = "femlm", mc_origin_bis = match.call(),
                     call_env_bis = call_env_bis, only.env = only.env, only.coef = only.coef,
                     env = env, ...), silent = TRUE)

    if("try-error" %in% class(res)){
        stop(format_error_msg(res, "femlm"))
    }

    return(res)
}

#' @rdname femlm
#' @export
fenegbin = function(fml, data, vcov, theta.init, start = 0, fixef, fixef.rm = "perfect",
                    offset, subset, split, fsplit, split.keep, split.drop,
                    cluster, se, ssc, panel.id,
                    fixef.tol = 1e-5, fixef.iter = 10000, nthreads = getFixest_nthreads(),
                    lean = FALSE, verbose = 0, warn = TRUE, notes = getFixest_notes(),
                    combine.quick, mem.clean = FALSE, only.env = FALSE, only.coef = FALSE, env, ...){

    # We control for the problematic argument family
    if("family" %in% names(match.call())){
        stop("Function fenegbin does not accept the argument 'family'.")
    }

    # This is just an alias
    set_defaults("fixest_estimation")
    call_env_bis = new.env(parent = parent.frame())

    res = try(feNmlm(fml = fml, data=data, family = "negbin", theta.init = theta.init,
                     start = start, fixef = fixef, fixef.rm = fixef.rm, offset = offset,
                     subset = subset, split = split, fsplit = fsplit,
                     split.keep = split.keep, split.drop = split.drop,
                     vcov = vcov, cluster = cluster,
                     se = se, ssc = ssc, panel.id = panel.id, fixef.tol = fixef.tol,
                     fixef.iter = fixef.iter, nthreads = nthreads, lean = lean,
                     verbose = verbose, warn = warn, notes = notes, combine.quick = combine.quick,
                     mem.clean = mem.clean, only.env = only.env, only.coef = only.coef,
                     origin = "fenegbin", mc_origin_bis = match.call(),
                     call_env_bis = call_env_bis, env = env, ...), silent = TRUE)

    if("try-error" %in% class(res)){
        stop(format_error_msg(res, "fenegbin"))
    }

    return(res)
}

#' @rdname feglm
#' @export
fepois = function(fml, data, vcov, offset, weights, subset, split, fsplit, split.keep, split.drop,
                  cluster, se, ssc, panel.id, start = NULL, etastart = NULL,
                  mustart = NULL, fixef, fixef.rm = "perfect", fixef.tol = 1e-6,
                  fixef.iter = 10000, collin.tol = 1e-10, glm.iter = 25, glm.tol = 1e-8,
                  nthreads = getFixest_nthreads(), lean = FALSE, warn = TRUE, notes = getFixest_notes(),
                  verbose = 0, combine.quick, mem.clean = FALSE, only.env = FALSE,
                  only.coef = FALSE, env, ...){

    # We control for the problematic argument family
    if("family" %in% names(match.call())){
        stop("Function fepois does not accept the argument 'family'.")
    }

    # This is just an alias
    set_defaults("fixest_estimation")
    call_env_bis = new.env(parent = parent.frame())

    res = try(feglm(fml = fml, data = data, family = "poisson", offset = offset,
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
                    call_env_bis = call_env_bis, env = env, ...), silent = TRUE)

    if("try-error" %in% class(res)){
        stop(format_error_msg(res, "fepois"))
    }

    return(res)
}



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
#' See also [`summary.fixest`] to see the results with the appropriate standard-errors, [`fixef.fixest`] to extract the fixed-effects coefficients, and the function [`etable`] to visualize the results of multiple estimations.
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
#' n = 100
#' x = rnorm(n, 1, 5)**2
#' y = rnorm(n, -1, 5)**2
#' z1 = rpois(n, x*y) + rpois(n, 2)
#' base = data.frame(x, y, z1)
#'
#' # Estimating a 'linear' relation:
#' est1_L = femlm(z1 ~ log(x) + log(y), base)
#' # Estimating the same 'linear' relation using a 'non-linear' call
#' est1_NL = feNmlm(z1 ~ 1, base, NL.fml = ~a*log(x)+b*log(y), NL.start = list(a=0, b=0))
#' # we compare the estimates with the function esttable (they are identical)
#' etable(est1_L, est1_NL)
#'
#' # Now generating a non-linear relation (E(z2) = x + y + 1):
#' z2 = rpois(n, x + y) + rpois(n, 1)
#' base$z2 = z2
#'
#' # Estimation using this non-linear form
#' est2_NL = feNmlm(z2 ~ 0, base, NL.fml = ~log(a*x + b*y),
#'                NL.start = 2, lower = list(a=0, b=0))
#' # we can't estimate this relation linearily
#' # => closest we can do:
#' est2_L = femlm(z2 ~ log(x) + log(y), base)
#'
#' # Difference between the two models:
#' etable(est2_L, est2_NL)
#'
#' # Plotting the fits:
#' plot(x, z2, pch = 18)
#' points(x, fitted(est2_L), col = 2, pch = 1)
#' points(x, fitted(est2_NL), col = 4, pch = 2)
#'
#' @export
feNmlm = function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), NL.fml, vcov,
                  fixef, fixef.rm = "perfect", NL.start, lower, upper, NL.start.init,
                  offset, subset, split, fsplit, split.keep, split.drop,
                  cluster, se, ssc, panel.id,
                  start = 0, jacobian.method="simple", useHessian = TRUE,
                  hessian.args = NULL, opt.control = list(), nthreads = getFixest_nthreads(),
                  lean = FALSE, verbose = 0, theta.init, fixef.tol = 1e-5, fixef.iter = 10000,
                  deriv.tol = 1e-4, deriv.iter = 1000, warn = TRUE, notes = getFixest_notes(),
                  combine.quick, mem.clean = FALSE, only.env = FALSE, only.coef = FALSE, env, ...){

    time_start = proc.time()

    if(missing(env)){

        set_defaults("fixest_estimation")
        call_env = new.env(parent = parent.frame())

        env = try(fixest_env(fml = fml, data = data, family = family, NL.fml = NL.fml,
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
                             call_env = call_env, computeModel0 = TRUE, ...), silent = TRUE)

    } else if((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
        stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
    }

    check_arg(only.env, "logical scalar")
    if(only.env){
        return(env)
    }

    if("try-error" %in% class(env)){
        mc = match.call()
        origin = ifelse(is.null(mc$origin), "feNmlm", mc$origin)
        stop(format_error_msg(env, origin))
    }

    verbose = get("verbose", env)
    if(verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep="")

    IN_MULTI = get("IN_MULTI", env)
    is_multi_root = get("is_multi_root", env)
    if(is_multi_root){
        on.exit(release_multi_notes())
        assign("is_multi_root", FALSE, env)
    }

    # Split ----

    do_split = get("do_split", env)
    if(do_split){

        res = multi_split(env, feNmlm)

        return(res)
    }

    # Multi fixef ----

    do_multi_fixef = get("do_multi_fixef", env)
    if(do_multi_fixef){

        res = multi_fixef(env, feNmlm)

        return(res)
    }

    # Multi LHS and RHS ----

    do_multi_lhs = get("do_multi_lhs", env)
    do_multi_rhs = get("do_multi_rhs", env)
    if(do_multi_lhs || do_multi_rhs){

        res = multi_LHS_RHS(env, feNmlm)

        return(res)

    }

    # Regular estimation ----


    # Objects needed for optimization + misc
    start = get("start", env)
    lower = get("lower", env)
    upper = get("upper", env)
    gradient = get("gradient", env)
    hessian = get("hessian", env)
    family = get("family", env)
    isLinear = get("isLinear", env)
    isNonLinear = get("isNL", env)
    opt.control = get("opt.control", env)
    only.coef = get("only.coef", env)

    lhs = get("lhs", env)

    family = get("family", env)
    famFuns = get("famFuns", env)
    params = get("params", env)

    isFixef = get("isFixef", env)
    onlyFixef = !isLinear && !isNonLinear && isFixef

    #
    # Model 0 + theta init
    #

    theta.init = get("theta.init", env)
    model0 = get_model_null(env, theta.init)

    # For the negative binomial:
    if(family == "negbin"){
        theta.init = get("theta.init", env)
        if(is.null(theta.init)){
            theta.init = model0$theta
        }

        params = c(params, ".theta")
        start = c(start, theta.init)
        names(start) = params
        upper = c(upper, 10000)
        lower = c(lower, 1e-3)

        assign("params", params, env)
    }

    assign("model0", model0, env)

    # the result
    res = get("res", env)

    # NO VARIABLE -- ONLY FIXED-EFFECTS
    if(onlyFixef){
        if(family == "negbin"){
            stop("To estimate the negative binomial model, you need at least one variable. (The estimation of the model with only the fixed-effects is not implemented.)")
        }

        res = femlm_only_clusters(env)
        res$onlyFixef = TRUE

        return(res)
    }

    # warnings => to avoid accumulation, but should appear even if the user stops the algorithm
    on.exit(warn_fixef_iter(env, stack_multi = IN_MULTI))

    #
    # Maximizing the likelihood
    #

    opt = try(stats::nlminb(start=start, objective=femlm_ll, env=env, lower=lower, upper=upper, gradient=gradient, hessian=hessian, control=opt.control), silent = TRUE)

    if("try-error" %in% class(opt)){
        # We return the coefficients (can be interesting for debugging)
        iter = get("iter", env)
        origin = get("origin", env)
        warning_msg = paste0("[", origin, "] Optimization failed at iteration ", iter, ". Reason: ", gsub("^[^\n]+\n *(.+\n)", "\\1", opt))
        if(IN_MULTI){
            stack_multi_notes(warning_msg)
            return(fixest_NA_results(env))
        } else if(!"coef_evaluated" %in% names(env)){
            # big problem right from the start
            stop(warning_msg)
        } else {
            coef = get("coef_evaluated", env)
            warning(warning_msg, " Last evaluated coefficients returned.", call. = FALSE)
            return(coef)
        }

    } else {
        convStatus = TRUE
        warning_msg = ""
        if(!opt$message %in% c("X-convergence (3)", "relative convergence (4)", "both X-convergence and relative convergence (5)")){
            warning_msg = " The optimization algorithm did not converge, the results are not reliable."
            convStatus = FALSE
        }

        coef = opt$par
    }

    if(only.coef){
        if(convStatus == FALSE){
            if(IN_MULTI){
                stack_multi_notes(warning_msg)
            } else {
                warning(warning_msg)
            }
        }

        names(coef) = params
        return(coef)
    }



    # The Hessian
    hessian = femlm_hessian(coef, env = env)
    # we add the names of the non linear variables in the hessian
    if(isNonLinear || family == "negbin"){
        dimnames(hessian) = list(params, params)
    }

    # we create the Hessian without the bounded parameters
    hessian_noBounded = hessian

    # Handling the bounds
    if(!isNonLinear){
        NL.fml = NULL
        bounds = NULL
        isBounded = NULL
    } else {
        nonlinear.params = get("nonlinear.params", env)
        # we report the bounds & if the estimated parameters are bounded
        upper_bound = upper[nonlinear.params]
        lower_bound = lower[nonlinear.params]

        # 1: are the estimated parameters at their bounds?
        coef_NL = coef[nonlinear.params]
        isBounded = rep(FALSE, length(params))
        isBounded[1:length(coef_NL)] = (coef_NL == lower_bound) | (coef_NL == upper_bound)

        # 2: we save the bounds
        upper_bound_small = upper_bound[is.finite(upper_bound)]
        lower_bound_small = lower_bound[is.finite(lower_bound)]
        bounds = list()
        if(length(upper_bound_small) > 0) bounds$upper = upper_bound_small
        if(length(lower_bound_small) > 0) bounds$lower = lower_bound_small
        if(length(bounds) == 0){
            bounds = NULL
        }

        # 3: we update the Hessian (basically, we drop the bounded element)
        if(any(isBounded)){
            hessian_noBounded = hessian[-which(isBounded), -which(isBounded), drop = FALSE]

            boundText = ifelse(coef_NL == upper_bound, "Upper bounded", "Lower bounded")[isBounded]

            attr(isBounded, "type") = boundText
        }

    }

    # Variance

    var = NULL
    try(var <- solve(hessian_noBounded), silent = TRUE)
    if(is.null(var)){
        warning_msg = paste(warning_msg, "The information matrix is singular: presence of collinearity.")
        var = hessian_noBounded * NA
        se = diag(var)
    } else {
        se = diag(var)
        se[se < 0] = NA
        se = sqrt(se)
    }

    # Warning message
    if(nchar(warning_msg) > 0){
        if(warn){

            if(IN_MULTI){
                stack_multi_notes(warning_msg)
            } else {
                warning("[femlm]:", warning_msg, call. = FALSE)
                options("fixest_last_warning" = proc.time())
            }
        }
    }

    # To handle the bounded coefficient, we set its SE to NA
    if(any(isBounded)){
        se = se[params]
        names(se) = params
    }

    zvalue = coef/se
    pvalue = 2*pnorm(-abs(zvalue))

    # We add the information on the bound for the se & update the var to drop the bounded vars
    se_format = se
    if(any(isBounded)){
        se_format[!isBounded] = decimalFormat(se_format[!isBounded])
        se_format[isBounded] = boundText
    }

    coeftable = data.frame("Estimate"=coef, "Std. Error"=se_format, "z value"=zvalue, "Pr(>|z|)"=pvalue, stringsAsFactors = FALSE)
    names(coeftable) = c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
    row.names(coeftable) = params

    attr(se, "type") = attr(coeftable, "type") = "IID"

    mu_both = get_mu(coef, env, final = TRUE)
    mu = mu_both$mu
    exp_mu = mu_both$exp_mu

    # calcul pseudo r2
    loglik = -opt$objective # moins car la fonction minimise
    ll_null = model0$loglik

    # dummies are constrained, they don't have full dof (cause you need to take one value off for unicity)
    # this is an approximation, in some cases there can be more than one ref. But good approx.
    nparams = res$nparams
    pseudo_r2 = 1 - (loglik - nparams + 1) / ll_null

    # Calcul residus
    expected.predictor = famFuns$expected.predictor(mu, exp_mu, env)
    residuals = lhs - expected.predictor

    # calcul squared corr
    if(cpp_isConstant(expected.predictor)){
        sq.cor = NA
    } else {
        sq.cor = stats::cor(lhs, expected.predictor)**2
    }
    ssr_null = cpp_ssr_null(lhs)

    # The scores
    scores = femlm_scores(coef, env)
    if(isNonLinear){
        # we add the names of the non linear params in the score
        colnames(scores) = params
    }

    n = length(lhs)

    # Saving
    res$coefficients = coef
    res$coeftable = coeftable
    res$loglik = loglik
    res$iterations = opt$iterations
    res$ll_null = ll_null
    res$ssr_null = ssr_null
    res$pseudo_r2 = pseudo_r2
    res$message = opt$message
    res$convStatus = convStatus
    res$sq.cor = sq.cor
    res$fitted.values = expected.predictor
    res$hessian = hessian

    dimnames(var) = list(params, params)
    res$cov.iid = var
    # for compatibility with conleyreg
    res$cov.unscaled = res$cov.iid

    res$se = se
    res$scores = scores
    res$family = family
    res$residuals = residuals

    # The value of mu (if cannot be recovered from fitted())
    if(family == "logit"){
        qui_01 = expected.predictor %in% c(0, 1)
        if(any(qui_01)){
            res$mu = mu
        }
    } else if(family %in% c("poisson", "negbin")){
        qui_0 = expected.predictor == 0
        if(any(qui_0)){
            res$mu = mu
        }
    }


    if(!is.null(bounds)){
        res$bounds = bounds
        res$isBounded = isBounded
    }

    # Fixed-effects
    if(isFixef){

        useExp_fixefCoef = family %in% c("poisson")

        sumFE = attr(mu, "sumFE")
        if(useExp_fixefCoef){
            sumFE = rpar_log(sumFE, env)
        }

        res$sumFE = sumFE

        # The LL and SSR with FE only
        if("ll_fe_only" %in% names(env)){
            res$ll_fe_only = get("ll_fe_only", env)
            res$ssr_fe_only = get("ssr_fe_only", env)
        } else {
            # we need to compute it

            # indicator of whether we compute the exp(mu)
            useExp = family %in% c("poisson", "logit", "negbin")

            # mu, using the offset
            if(!is.null(res$offset)){
                mu_noDum = res$offset
            } else {
                mu_noDum = 0
            }

            if(length(mu_noDum) == 1) mu_noDum = rep(mu_noDum, n)

            exp_mu_noDum = NULL
            if(useExp_fixefCoef){
                exp_mu_noDum = rpar_exp(mu_noDum, env)
            }

            assign("fixef.tol", 1e-4, env) # no need of supa precision
            dummies = getDummies(mu_noDum, exp_mu_noDum, env, coef)

            exp_mu = NULL
            if(useExp_fixefCoef){
                # despite being called mu, it is in fact exp(mu)!!!
                exp_mu = exp_mu_noDum*dummies
                mu = rpar_log(exp_mu, env)
            } else {
                mu = mu_noDum + dummies
                if(useExp){
                    exp_mu = rpar_exp(mu, env)
                }
            }

            res$ll_fe_only = famFuns$ll(lhs, mu, exp_mu, env, coef)
            ep = famFuns$expected.predictor(mu, exp_mu, env)
            res$ssr_fe_only = cpp_ssq(lhs - ep)

        }
    }

    if(family == "negbin"){
        theta = coef[".theta"]
        res$theta = theta

        if(notes && theta > 1000){
            msg = paste0("Very high value of theta (", theta, "). There is no sign of overdispersion, you may consider a Poisson model.")
            if(IN_MULTI){
                stack_multi_notes(msg)
            } else {
                message(msg)
            }
        }

    }

    class(res) = "fixest"

    if(verbose > 0){
        cat("\n")
    }

    do_summary = get("do_summary", env)
    if(do_summary){
        vcov = get("vcov", env)
        lean = get("lean", env)
        ssc = get("ssc", env)
        agg = get("agg", env)
        summary_flags = get("summary_flags", env)

        # To compute the RMSE and lean = TRUE
        if(lean) res$ssr = cpp_ssq(res$residuals)

        res = summary(res, vcov = vcov, ssc = ssc, agg = agg, lean = lean, summary_flags = summary_flags)
    }

    return(res)
}


#' Estimates a `fixest` estimation from a `fixest` environment
#'
#' This is a function advanced users which allows to estimate any `fixest` estimation from a `fixest` environment obtained with `only.env = TRUE` in a `fixest` estimation.
#'
#' @param env An environment obtained from a `fixest` estimation with `only.env = TRUE`. This is intended for advanced users so there is no error handling: any other kind of input will fail with a poor error message.
#' @param y A vector representing the dependent variable. Should be of the same length as the number of observations in the initial estimation.
#' @param X A matrix representing the independent variables. Should be of the same dimension as in the initial estimation.
#' @param weights A vector of weights (i.e. with only positive values). Should be of the same length as the number of observations in the initial estimation. If identical to the scalar 1, this will mean that no weights will be used in the estimation.
#' @param endo A matrix representing the endogenous regressors in IV estimations. It should be of the same dimension as the original endogenous regressors.
#' @param inst A matrix representing the instruments in IV estimations. It should be of the same dimension as the original instruments.
#'
#' @return
#'
#' It returns the results of a `fixest` estimation: the one that was summoned when obtaining the environment.
#'
#' @details
#'
#' This function has been created for advanced users, mostly to avoid overheads when making simulations with `fixest`.
#'
#' How can it help you make simulations? First make a core estimation with `only.env = TRUE`, and usually with `only.coef = TRUE` (to avoid having extra things that take time to compute). Then loop while modifying the appropriate things directly in the environment. Beware that if you make a mistake here (typically giving stuff of the wrong length), then you can make the R session crash because there is no more error-handling! Finally estimate with `est_env(env = core_env)` and store the results.
#'
#' Instead of `est_env`, you could use directly `fixest` estimations too, like `feols`, since they accept the `env` argument. The function `est_env` is only here to add a bit of generality to avoid the trouble to the user to write conditions (look at the source, it's just a one liner).
#'
#' Objects of main interest in the environment are:
#' \itemize{
#' \item{lhs}{The left hand side, or dependent variable.}
#' \item{linear.mat}{The matrix of the right-hand-side, or explanatory variables.}
#' \item{iv_lhs}{The matrix of the endogenous variables in IV regressions.}
#' \item{iv.mat}{The matrix of the instruments in IV regressions.}
#' \item{weights.value}{The vector of weights.}
#' }
#'
#' I strongly discourage changing the dimension of any of these elements, or else crash can occur. However, you can change their values at will (given the dimension stay the same). The only exception is the weights, which tolerates changing its dimension: it can be identical to the scalar `1` (meaning no weights), or to something of the length the number of observations.
#'
#' I also discourage changing anything in the fixed-effects (even their value) since this will almost surely lead to a crash.
#'
#' Note that this function is mostly useful when the overheads/estimation ratio is high. This means that OLS will benefit the most from this function. For GLM/Max.Lik. estimations, the ratio is small since the overheads is only a tiny portion of the total estimation time. Hence this function will be less useful for these models.
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
#' bboot = function(x, n_sim = 100){
#'   # We bootstrap on the weights
#'   # Works with fixed-effects/IVs
#'   #  and with any fixest function that accepts weights
#'
#'   core_env = update(x, only.coef = TRUE, only.env = TRUE)
#'   n_obs = x$nobs
#'
#'   res_all = vector("list", n_sim)
#'   for(i in 1:n_sim){
#'     # # begin: NOT RUN
#'     # # We could directly assign in the environment:
#'     # assign("weights.value", rexp(n_obs, rate = 1), core_env)
#'     # res_all[[i]] = est_env(env = core_env)
#'     #   end: NOT RUN
#'
#'     # Instead we can use the argument weights, which does the same
#'     res_all[[i]] = est_env(env = core_env, weights = rexp(n_obs, rate = 1))
#'   }
#'
#'   do.call(rbind, res_all)
#' }
#'
#'
#' est = feols(mpg ~ wt + hp, mtcars)
#'
#' boot_res = bboot(est)
#' coef = colMeans(boot_res)
#' std_err = apply(boot_res, 2, sd)
#'
#' # Comparing the results with the main estimation
#' coeftable(est)
#' cbind(coef, std_err)
#'
#' @export
est_env = function(env, y, X, weights, endo, inst){
    # No check whatsoever: for advanced users

    if(!missnull(y)){
        assign("lhs", y, env)
    }

    if(!missnull(X)){
        assign("linear.mat", X, env)
    }

    if(!missnull(weights)){
        assign("weights.value", weights, env)
    }

    if(!missnull(endo)){
        assign("is_lhs", endo, env)
    }

    if(!missnull(inst)){
        assign("iv.mat", inst, env)
    }

    switch(env$res$method_type,
           "feols"  = feols(env = env),
           "feglm"  = feglm.fit(env = env),
           "feNmlm" = feNmlm(env = env))
}

#### Delayed warnings and notes ----

setup_multi_notes = function(){
    # We reset all the notes
    options("fixest_multi_notes" = NULL)
}

stack_multi_notes = function(note){
    all_notes = getOption("fixest_multi_notes")
    options("fixest_multi_notes" = c(all_notes, note))
}

release_multi_notes = function(){
    notes = getOption("fixest_multi_notes")

    if(length(notes) > 0){
        dict = c("\\$collin\\.var" = "Some variables have been removed because of collinearity (see $collin.var).",
                 "collinear with the fixed-effects" = "All the variables are collinear with the fixed-effects (the estimation is void).",
                 "virtually constant and equal to 0" = "All the variables are virtually constant and equal to 0 (the estimation is void).",
                 "Very high value of theta" = "Very high value of theta, there is no sign of overdispersion, you may consider a Poisson model.")

        for(i in seq_along(dict)){
            d_value = dict[i]
            d_name = names(dict)[i]

            x_i = grepl(d_name, notes)
            x = notes[x_i]
            # we prefer specific messages
            if(length(unique(x)) == 1){
                next
            } else {
                notes[x_i] = d_value
            }
        }

        tn = table(notes)
        new_notes = paste0("[x ", tn, "] ", names(tn))

        message("Notes from the estimations:\n", dsb("'\n'c?new_notes"))
    }

    options("fixest_multi_notes" = NULL)
}


warn_fixef_iter = function(env, stack_multi = FALSE){
    # Show warnings related to the nber of times the maximum of iterations was reached

    # For fixed-effect
    fixef.iter = get("fixef.iter", env)
    fixef.iter.limit_reached = get("fixef.iter.limit_reached", env)
    origin = get("origin", env)
    warn = get("warn", env)

    if(!warn) return(invisible(NULL))

    goWarning = FALSE
    warning_msg = ""
    if(fixef.iter.limit_reached > 0){
        goWarning = TRUE
        warning_msg = paste0(origin, ": [Getting the fixed-effects] iteration limit reached (", fixef.iter, ").", ifelse(fixef.iter.limit_reached > 1, paste0(" (", fixef.iter.limit_reached, " times.)"), " (Once.)"))
    }

    # For the fixed-effect derivatives
    deriv.iter = get("deriv.iter", env)
    deriv.iter.limit_reached = get("deriv.iter.limit_reached", env)

    if(deriv.iter.limit_reached > 0){
        prefix = ifelse(goWarning, paste0("\n", sprintf("% *s", nchar(origin) + 2, " ")), paste0(origin, ": "))
        warning_msg = paste0(warning_msg, prefix, "[Getting fixed-effects derivatives] iteration limit reached (", deriv.iter, ").", ifelse(deriv.iter.limit_reached > 1, paste0(" (", deriv.iter.limit_reached, " times.)"), " (Once.)"))
        goWarning = TRUE
    }

    if(goWarning){
        if(stack_multi){
            stack_multi_notes(warning_msg)
        } else {
            warning(warning_msg, call. = FALSE, immediate. = TRUE)
        }
    }

}


warn_step_halving = function(env, stack_multi = FALSE){

    nb_sh = get("nb_sh", env)
    warn = get("warn", env)

    if(!warn) return(invisible(NULL))

    if(nb_sh > 0){
        msg = paste0("feglm: Step halving due to non-finite deviance (", ifelse(nb_sh > 1, paste0(nb_sh, " times"), "once"), ").")
        if(stack_multi){
            stack_multi_notes(msg)
        } else {
            warning(msg, call. = FALSE, immediate. = TRUE)
        }
    }

}


format_error_msg = function(x, origin){
    # Simple formatting of the error msg

    # LATER:
    # - for object not found: provide a better error msg by calling the name of the missing
    #   argument => likely I'll need a match.call argument

    x = gsub("\n+$", "", x)

    if(grepl("^Error (in|:|: in) (fe|fixest|fun|fml_split)[^\n]+\n", x)){
        res = gsub("^Error (in|:|: in) (fe|fixest|fun|fml_split)[^\n]+\n *(.+)", "\\3", x)
    } else if(grepl("[Oo]bject '.+' not found", x) || grepl("memory|cannot allocate", x)) {
        res = x
    } else {
        res = paste0(x, "\nThis error was unforeseen by the author of the function ", origin, ". If you think your call to the function is legitimate, could you report?")
    }
    res
}

#### Multiple estimation tools ----

multi_split = function(env, fun){
    split = get("split", env)
    split.items = get("split.items", env)
    split.id = get("split.id", env)
    split.name = get("split.name", env)

    assign("do_split", FALSE, env)

    n_split = length(split.id)
    res_all = vector("list", n_split)
    I = 1
    for(i in split.id){
        if(i == 0){
            my_env = reshape_env(env)
            my_res = fun(env = my_env)
        } else {
            my_res = fun(env = reshape_env(env, obs2keep = which(split == i)))
        }

        res_all[[I]] = my_res
        I = I + 1
    }

    values = list(sample = split.items)

    # result
    res_multi = setup_multi(res_all, values, var = split.name)

    return(res_multi)
}



multi_LHS_RHS = function(env, fun){
    do_multi_lhs = get("do_multi_lhs", env)
    do_multi_rhs = get("do_multi_rhs", env)

    assign("do_multi_lhs", FALSE, env)
    assign("do_multi_rhs", FALSE, env)

    nthreads = get("nthreads", env)

    # IMPORTANT NOTE:
    # contrary to feols, the preprocessing is only a small fraction of the
    # computing time in ML models
    # Therefore we don't need to optimize processing as hard as in FEOLS
    # because the gains are only marginal

    fml = get("fml", env)

    # LHS
    lhs_names = get("lhs_names", env)
    lhs = get("lhs", env)

    if(do_multi_lhs == FALSE){
        lhs = list(lhs)
    }

    # RHS
    if(do_multi_rhs){
        rhs_info_stepwise = get("rhs_info_stepwise", env)
        multi_rhs_fml_full = rhs_info_stepwise$fml_all_full
        multi_rhs_fml_sw = rhs_info_stepwise$fml_all_sw
        multi_rhs_cumul = rhs_info_stepwise$is_cumul

        linear_core = get("linear_core", env)
        rhs_sw = get("rhs_sw", env)

    } else {
        multi_rhs_fml_full = list(.xpd(rhs = fml[[3]]))
        multi_rhs_cumul = FALSE
        linear.mat = get("linear.mat", env)
        linear_core = list(left = linear.mat, right = 1)
        rhs_sw = list(1)
    }

    isLinear_left = length(linear_core$left) > 1
    isLinear_right = length(linear_core$right) > 1

    n_lhs = length(lhs)
    n_rhs = length(rhs_sw)
    res = vector("list", n_lhs * n_rhs)

    rhs_names = sapply(multi_rhs_fml_full, function(x) as.character(x)[[2]])

    for(i in seq_along(lhs)){
        for(j in seq_along(rhs_sw)){
            # reshaping the env => taking care of the NAs

            # Forming the RHS
            my_rhs = linear_core[1]

            if(multi_rhs_cumul){
                if(j > 1){
                    # if j == 1, the RHS is already in the data
                    my_rhs[1 + (2:j - 1)] = rhs_sw[2:j]
                }
            } else {
                my_rhs[2] = rhs_sw[j]
            }

            if(isLinear_right){
                my_rhs[[length(my_rhs) + 1]] = linear_core$right
            }

            n_all = lengths(my_rhs)
            if(any(n_all == 1)){
                my_rhs = my_rhs[n_all > 1]
            }

            if(length(my_rhs) == 0){
                my_rhs = 1
            } else {
                my_rhs = do.call("cbind", my_rhs)
            }

            if(length(my_rhs) == 1){
                is_na_current = !is.finite(lhs[[i]])
            } else {
                is_na_current = !is.finite(lhs[[i]]) | cpppar_which_na_inf_mat(my_rhs, nthreads)$is_na_inf
            }

            my_fml = .xpd(lhs = lhs_names[i], rhs = multi_rhs_fml_full[[j]])

            if(any(is_na_current)){
                my_env = reshape_env(env, which(!is_na_current), lhs = lhs[[i]], rhs = my_rhs, fml_linear = my_fml)
            } else {
                # We still need to check the RHS (only 0/1)
                my_env = reshape_env(env, lhs = lhs[[i]], rhs = my_rhs, fml_linear = my_fml, check_lhs = TRUE)
            }

            my_res = fun(env = my_env)

            res[[index_2D_to_1D(i, j, n_rhs)]] = my_res
        }
    }

    # Meta information for fixest_multi
    values = list(lhs = rep(lhs_names, each = n_rhs),
                  rhs = rep(rhs_names, n_lhs))

    # result
    res_multi = setup_multi(res, values)

    return(res_multi)
}


multi_fixef = function(env, estfun){
    # Honestly had I known it was so painful, I wouldn't have done it...
    assign("do_multi_fixef", FALSE, env)

    multi_fixef_fml_full = get("multi_fixef_fml_full", env)
    combine.quick = get("combine.quick", env)
    fixef.rm = get("fixef.rm", env)
    family = get("family", env)
    origin_type = get("origin_type", env)
    nthreads = get("nthreads", env)

    data = get("data", env)

    n_fixef = length(multi_fixef_fml_full)

    data_results = vector("list", n_fixef)
    for(i in 1:n_fixef){

        fml_fixef = multi_fixef_fml_full[[i]]

        if(length(all.vars(fml_fixef)) > 0){

            #
            # Evaluation of the fixed-effects
            #

            fixef_terms_full = fixef_terms(fml_fixef)

            # fixef_terms_full computed in the formula section
            fixef_terms = fixef_terms_full$fml_terms

            # FEs
            fixef_df = error_sender(prepare_df(fixef_terms_full$fe_vars, data, combine.quick),
                                    "Problem evaluating the fixed-effects part of the formula:\n")

            fixef_vars = names(fixef_df)

            # Slopes
            isSlope = any(fixef_terms_full$slope_flag != 0)
            slope_vars_list = list(0)
            if(isSlope){

                slope_df = error_sender(prepare_df(fixef_terms_full$slope_vars, data),
                                        "Problem evaluating the variables with varying slopes in the fixed-effects part of the formula:\n")

                slope_flag = fixef_terms_full$slope_flag
                slope_vars = fixef_terms_full$slope_vars
                slope_vars_list = fixef_terms_full$slope_vars_list

                # Further controls
                not_numeric = !sapply(slope_df, is.numeric)
                if(any(not_numeric)){
                    stop("In the fixed-effects part of the formula (i.e. in ", as.character(fml_fixef[2]), "), variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope_df)[not_numeric], "s.is.quote"), " not.")
                }

                # slope_flag: 0: no Varying slope // > 0: varying slope AND fixed-effect // < 0: varying slope WITHOUT fixed-effect
                onlySlope = all(slope_flag < 0)

            }

            # fml update
            fml_fixef = .xpd(rhs = fixef_terms)

            #
            # NA
            #

            for(j in seq_along(fixef_df)){
                if(!is.numeric(fixef_df[[j]]) && !is.character(fixef_df[[j]])){
                    fixef_df[[j]] = as.character(fixef_df[[j]])
                }
            }

            is_NA = !complete.cases(fixef_df)

            if(isSlope){
                # Convert to double
                who_not_double = which(sapply(slope_df, is.integer))
                for(j in who_not_double){
                    slope_df[[j]] = as.numeric(slope_df[[j]])
                }

                info = cpppar_which_na_inf_df(slope_df, nthreads)
                if(info$any_na_inf){
                    is_NA = is_NA | info$is_na_inf
                }
            }

            if(any(is_NA)){
                # Remember that isFixef is FALSE so far => so we only change the reg vars
                my_env = reshape_env(env = env, obs2keep = which(!is_NA))

                # NA removal in fixef
                fixef_df = fixef_df[!is_NA, , drop = FALSE]

                if(isSlope){
                    slope_df = slope_df[!is_NA, , drop = FALSE]
                }
            } else {
                my_env = new.env(parent = env)
            }

            # We remove the linear part if needed


            if(get("do_multi_rhs", env)){
                linear_core = get("linear_core", my_env)
                if("(Intercept)" %in% colnames(linear_core$left)){
                    int_col = which("(Intercept)" %in% colnames(linear_core$left))
                    if(ncol(linear_core$left) == 1){
                        linear_core$left = 1
                    } else {
                        linear_core$left = linear_core$left[, -int_col, drop = FALSE]
                    }
                    assign("linear_core", linear_core, my_env)
                }
            } else {
                linear.mat = get("linear.mat", my_env)
                if("(Intercept)" %in% colnames(linear.mat)){
                    int_col = which("(Intercept)" %in% colnames(linear.mat))
                    if(ncol(linear.mat) == 1){
                        assign("linear.mat", 1, my_env)
                    } else {
                        assign("linear.mat", linear.mat[, -int_col, drop = FALSE], my_env)
                    }
                }
            }

            # We assign the fixed-effects
            lhs = get("lhs", my_env)

            # We delay the computation by using isSplit = TRUE and split.full = FALSE
            # Real QUF will be done in the last reshape env
            info_fe = setup_fixef(fixef_df = fixef_df, lhs = lhs, fixef_vars = fixef_vars, fixef.rm = fixef.rm, family = family, isSplit = TRUE, split.full = FALSE, origin_type = origin_type, isSlope = isSlope, slope_flag = slope_flag, slope_df = slope_df, slope_vars_list = slope_vars_list, nthreads = nthreads)

            fixef_id        = info_fe$fixef_id
            fixef_names     = info_fe$fixef_names
            fixef_sizes     = info_fe$fixef_sizes
            fixef_table     = info_fe$fixef_table
            sum_y_all       = info_fe$sum_y_all
            lhs             = info_fe$lhs

            obs2remove      = info_fe$obs2remove
            fixef_removed   = info_fe$fixef_removed
            message_fixef   = info_fe$message_fixef

            slope_variables = info_fe$slope_variables
            slope_flag      = info_fe$slope_flag

            fixef_id_res    = info_fe$fixef_id_res
            fixef_sizes_res = info_fe$fixef_sizes_res
            new_order       = info_fe$new_order

            assign("isFixef", TRUE, my_env)
            assign("new_order_original", new_order, my_env)
            assign("fixef_names", fixef_names, my_env)
            assign("fixef_vars", fixef_vars, my_env)

            assign_fixef_env(env, family, origin_type, fixef_id, fixef_sizes, fixef_table, sum_y_all, slope_flag, slope_variables, slope_vars_list)

            #
            # Formatting the fixef stuff from res
            #

            # fml & fixef_vars => other stuff will be taken care of in reshape
            res = get("res", my_env)
            res$fml_all$fixef = fml_fixef
            res$fixef_vars = fixef_vars
            if(isSlope){
                res$fixef_terms = fixef_terms
            }
            assign("res", res, my_env)

            #
            # Last reshape
            #

            my_env_est = reshape_env(my_env, assign_fixef = TRUE)

        } else {
            # No fixed-effect // new.env is indispensable => otherwise multi RHS/LHS not possible
            my_env_est = reshape_env(env)
        }

        data_results[[i]] = estfun(env = my_env_est)
    }

    index = list(fixef = n_fixef)
    fixef_names = sapply(multi_fixef_fml_full, function(x) as.character(x)[[2]])
    values = list(fixef = fixef_names)

    res_multi = setup_multi(data_results, values)

    if("lhs" %in% names(attr(res_multi, "meta")$index)){
        res_multi = res_multi[lhs = TRUE]
    }

    return(res_multi)

}
