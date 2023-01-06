# additional methods ----

# Here we add common statistical functions

#' Extracts the number of observations form a \code{fixest} object
#'
#' This function simply extracts the number of observations form a \code{fixest} object, obtained using the functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams summary.fixest
#'
#' @param ... Not currently used.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @return
#' It returns an interger.
#'
#' @examples
#'
#' # simple estimation on iris data with "Species" fixed-effects
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' nobs(res)
#' logLik(res)
#'
#'
nobs.fixest = function(object, ...){
    object$nobs
}

#' Aikake's an information criterion
#'
#' This function computes the AIC (Aikake's, an information criterion) from a \code{fixest} estimation.
#'
#' @inheritParams nobs.fixest
#'
#' @param ... Optionally, more fitted objects.
#' @param k A numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC (i.e. \code{AIC=-2*LL+k*nparams}).
#'
#' @details
#' The AIC is computed as:
#' \deqn{AIC = -2\times LogLikelihood + k\times nbParams}
#' with k the penalty parameter.
#'
#' You can have more information on this criterion on \code{\link[stats]{AIC}}.
#'
#' @return
#' It return a numeric vector, with length the same as the number of objects taken as arguments.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Other statictics methods: \code{\link[fixest]{BIC.fixest}}, \code{\link[fixest]{logLik.fixest}}, \code{\link[fixest]{nobs.fixest}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'              Petal.Width | Species, iris)
#' res2 = femlm(Sepal.Length ~ Petal.Width | Species, iris)
#'
#' AIC(res1, res2)
#' BIC(res1, res2)
#'
#'
AIC.fixest = function(object, ..., k = 2){

    dots = list(...)
    if(length(dots) > 0){
        # we check consistency with observations
        nobs_all = c(nobs(object), sapply(dots, nobs))

        if(any(diff(nobs_all) != 0)){
            warning("Models are not all fitted to the same number of observations.")
        }

        otherAIC = sapply(dots, AIC)
    } else {
        otherAIC = c()
    }

    all_AIC = c(-2*logLik(object) + k*object$nparams, otherAIC)

    all_AIC
}

#' Bayesian information criterion
#'
#' This function computes the BIC (Bayesian information criterion) from a \code{fixest} estimation.
#'
#'
#' @inheritParams nobs.fixest
#'
#' @param ... Optionally, more fitted objects.
#'
#' @details
#' The BIC is computed as follows:
#' \deqn{BIC = -2\times LogLikelihood + \log(nobs)\times nbParams}
#' with k the penalty parameter.
#'
#' You can have more information on this criterion on \code{\link[stats]{AIC}}.
#'
#' @return
#' It return a numeric vector, with length the same as the number of objects taken as arguments.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Other statistics functions: \code{\link[fixest]{AIC.fixest}}, \code{\link[fixest]{logLik.fixest}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#' res2 = femlm(Sepal.Length ~ Petal.Width | Species, iris)
#'
#' AIC(res1, res2)
#' BIC(res1, res2)
#'
BIC.fixest = function(object, ...){

    dots = list(...)
    if(length(dots) > 0){
        # we check consistency with observations
        nobs_all = c(nobs(object), sapply(dots, nobs))

        if(any(diff(nobs_all) != 0)){
            warning("Models are not all fitted to the same number of observations.")
        }

        otherBIC = sapply(dots, BIC)
    } else {
        otherBIC = c()
    }

    all_BIC = c(-2*logLik(object) + object$nparams*log(nobs(object)), otherBIC)

    all_BIC
}

#' Extracts the log-likelihood
#'
#' This function extracts the log-likelihood from a \code{fixest} estimation.
#'
#' @inheritParams nobs.fixest
#'
#' @param ... Not currently used.
#'
#' @details
#' This function extracts the log-likelihood based on the model fit. You can have more information on the likelihoods in the details of the function \code{\link[fixest]{femlm}}.
#'
#' @return
#' It returns a numeric scalar.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Other statistics functions: \code{\link[fixest]{AIC.fixest}}, \code{\link[fixest]{BIC.fixest}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # simple estimation on iris data with "Species" fixed-effects
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' nobs(res)
#' logLik(res)
#'
#'
logLik.fixest = function(object, ...){

    if(object$method_type == "feols"){
        # if the summary is 'lean', then no way we can compute that
        resid = object$residuals
        if(is.null(resid)) resid = NA

        sigma = sqrt(mean(resid^2))
        n = length(resid)
        ll = -1/2/sigma^2 * sum(resid^2) - n * log(sigma) - n * log(2*pi)/2
    } else {
        ll = object$loglik
    }

    ll
}

#' Extracts the coefficients from a \code{fixest} estimation
#'
#' This function extracts the coefficients obtained from a model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#' @inheritParams etable
#'
#' @param agg Logical scalar, default is \code{TRUE}. If the coefficients of the estimation have been aggregated, whether to report the aggregated coefficients. If \code{FALSE}, the raw coefficients will be returned.
#' @param ... Not currently used.
#'
#' @details
#' The coefficients are the ones that have been found to maximize the log-likelihood of the specified model. More information can be found on the models from the estimations help pages: \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' Note that if the model has been estimated with fixed-effects, to obtain the fixed-effect coefficients, you need to use the function \code{\link[fixest]{fixef.fixest}}.
#'
#' @return
#' This function returns a named numeric vector.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{confint.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{etable}}, \code{\link[fixest]{fixef.fixest}}.
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # the coefficients of the variables:
#' coef(res)
#'
#' # the fixed-effects coefficients:
#' fixef(res)
#'
#'
coef.fixest = coefficients.fixest = function(object, keep, drop, order, agg = TRUE, ...){

    check_arg(keep, drop, order, "NULL character vector no na")
    check_arg(agg, "logical scalar")

    if(isTRUE(object$is_agg) && agg){
        res = object$coeftable[, 1]
        names(res) = rownames(object$coeftable)
    } else {
        res = object$coefficients
    }

    if(!missnull(keep) || !missnull(drop) || !missnull(order)){
        cnames = names(res)
        cnames = keep_apply(cnames, keep)
        cnames = drop_apply(cnames, drop)
        cnames = order_apply(cnames, order)

        if(length(cnames) == 0){
            return(numeric(0))
        }

        res = res[cnames]
    }

    if(identical(object$family, "negbin")){
        res = res[-length(res)]
    }


    # deltaMethod tweak
    if(is_calling_fun("deltaMethod")){
        sysOrigin = sys.parent()
        mc_DM = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin))

        if("parameterNames" %in% names(mc_DM)){
            PN = eval(mc_DM$parameterNames, parent.frame())

            check_value(PN, "character vector no na len(data)", .data = res,
                        .arg_name = "parameterNames", .up = 1)

            names(res) = PN
        }
    }

    res
}

#' @rdname coef.fixest
coefficients.fixest <- coef.fixest


#' Extracts fitted values from a \code{fixest} fit
#'
#' This function extracts the fitted values from a model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. The fitted values that are returned are the \emph{expected predictor}.
#'
#' @inheritParams nobs.fixest
#'
#' @param type Character either equal to \code{"response"} (default) or \code{"link"}. If \code{type="response"}, then the output is at the level of the response variable, i.e. it is the expected predictor \eqn{E(Y|X)}. If \code{"link"}, then the output is at the level of the explanatory variables, i.e. the linear predictor \eqn{X\cdot \beta}.
#' @param na.rm Logical, default is \code{TRUE}. If \code{FALSE} the number of observation returned will be the number of observations in the original data set, otherwise it will be the number of observations used in the estimation.
#' @param ... Not currently used.
#'
#' @details
#' This function returns the \emph{expected predictor} of a \code{fixest} fit. The likelihood functions are detailed in \code{\link[fixest]{femlm}} help page.
#'
#' @return
#' It returns a numeric vector of length the number of observations used to estimate the model.
#'
#' If \code{type = "response"}, the value returned is the expected predictor, i.e. the expected value of the dependent variable for the fitted model: \eqn{E(Y|X)}.
#' If \code{type = "link"}, the value returned is the linear predictor of the fitted model, that is \eqn{X\cdot \beta} (remind that \eqn{E(Y|X) = f(X\cdot \beta)}).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{resid.fixest}}, \code{\link[fixest]{predict.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{fixef.fixest}}.
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res_poisson = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris)
#'
#' # we extract the fitted values
#' y_fitted_poisson = fitted(res_poisson)
#'
#' # Same estimation but in OLS (Gaussian family)
#' res_gaussian = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris, family = "gaussian")
#'
#' y_fitted_gaussian = fitted(res_gaussian)
#'
#' # comparison of the fit for the two families
#' plot(iris$Sepal.Length, y_fitted_poisson)
#' points(iris$Sepal.Length, y_fitted_gaussian, col = 2, pch = 2)
#'
#'
fitted.fixest = fitted.values.fixest = function(object, type = c("response", "link"), na.rm = TRUE, ...){

    # Checking the arguments
    if(is_user_level_call()){
        validate_dots(suggest_args = "type")
    }

    type = match.arg(type)

    fit = predict(object)

    if(type == "response" || object$method_type == "feols"){
        res = fit
    } else if(!is.null(object$mu)){
        res = object$mu
    } else if(object$method == "femlm"){
        family = object$family
        famFuns = switch(family,
                         poisson = ml_poisson(),
                         negbin = ml_negbin(),
                         logit = ml_logit(),
                         gaussian = ml_gaussian())

        res = famFuns$linearFromExpected(fit)
    } else {
        res = object$family$linkfun(fit)
    }

    # Nota: obs can be removed: either because of NA, either because perfect fit
    # Shall I put perfect fit as NA since they're out of the estimation???
    # Still pondering...
    # Actually adding them means a lot of work to ensure consistency (also in predict...)
    if(!na.rm) res = fill_with_na(res, object)

    res
}

#' @rdname fitted.fixest
#' @method fitted.values fixest
fitted.values.fixest <- fitted.fixest

#' Extracts residuals from a \code{fixest} object
#'
#' This function extracts residuals from a fitted model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#'
#' @param type A character scalar, either \code{"response"} (default), \code{"deviance"}, \code{"pearson"}, or \code{"working"}. Note that the \code{"working"} corresponds to the residuals from the weighted least square and only applies to \code{\link[fixest]{feglm}} models.
#' @param na.rm Logical, default is \code{TRUE}. Whether to remove the observations with NAs from the original data set. If \code{FALSE}, then the vector returned is always of the same length as the original data set.
#' @param ... Not currently used.
#'
#'
#' @return
#' It returns a numeric vector of the length the number of observations used for the estimation (if \code{na.rm = TRUE}) or of the length of the original data set (if \code{na.rm = FALSE}).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{fitted.fixest}}, \code{\link[fixest]{predict.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{fixef.fixest}}.
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res_poisson = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris)
#'
#' # we plot the residuals
#' plot(resid(res_poisson))
#'
resid.fixest = residuals.fixest = function(object, type = c("response", "deviance", "pearson", "working"), na.rm = TRUE, ...){

    check_arg_plus(type, "match")
    check_arg_plus(na.rm, "logical scalar")

    method = object$method
    family = object$family

    r = object$residuals
    w = object[["weights"]]

    if(isTRUE(object$lean)){
        stop("The method 'resid.fixest' cannot be applied to a 'lean' fixest object. Please apply reestimate with 'lean = FALSE'.")
    }

    if(method %in% c("feols", "feols.fit") || (method %in% c("feNmlm", "femlm") && family == "gaussian")){

        if(type == "working") stop("Type 'working' only applies to models fitted via feglm (thus is not valid for feols).")

        if(type %in% c("deviance", "pearson") && !is.null(w)){
            res = r * sqrt(w)
        } else {
            res = r
        }

    } else if(method %in% c("fepois", "feglm")){

        if(type == "response"){
            res = r

        } else if(type == "working"){
            res = object$working_residuals

        } else {
            mu = object$fitted.values
            if(is.null(w)) w = rep(1, length(r))

            if(type == "deviance"){
                y = r + mu

                res = sqrt(pmax((object$family$dev.resids)(y, mu, w), 0))
                qui = y < mu
                res[qui] = -res[qui]

            } else if(type == "pearson"){
                res = r * sqrt(w)/sqrt(object$family$variance(object$fitted.values))

            }
        }


    } else {

        if(type == "working") stop("Type 'working' only applies to models fitted via feglm (thus is not valid for ", method, ").")

        if(type == "response"){
            res = r

        } else {
            # deviance or pearson
            mu = object$fitted.values
            if(is.null(w)) w = rep(1, length(r))

            theta = ifelse(family == "negbin", object$theta, 1)

            if(type == "deviance"){

                # dev.resids function
                if(family == "poisson"){
                    dev.resids = poisson()$dev.resids

                } else if(family == "logit"){
                    dev.resids = binomial()$dev.resids

                } else if(family == "negbin"){
                    dev.resids = function(y, mu, wt) 2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))

                }

                y = object$residuals + mu

                res = sqrt(pmax(dev.resids(y, mu, w), 0))
                qui = y < mu
                res[qui] = -res[qui]

            } else if(type == "pearson"){

                # variance function
                if(family == "poisson"){
                    variance = poisson()$variance

                } else if(family == "logit"){
                    variance = binomial()$variance

                } else if(family == "negbin"){
                    variance = function(mu) mu + mu^2/theta

                }

                res = r * sqrt(w)/sqrt(variance(mu))

            }
        }


    }

    if(!na.rm){
        res = fill_with_na(res, object)
    }

    res
}

#' @rdname resid.fixest
residuals.fixest <- resid.fixest

#' Predict method for \code{fixest} fits
#'
#' This function obtains prediction from a fitted model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#' @inheritParams fitted.fixest
#' @inheritParams summary.fixest
#'
#' @param newdata A data.frame containing the variables used to make the prediction. If not provided, the fitted expected (or linear if \code{type = "link"}) predictors are returned.
#' @param sample Either "estimation" (default) or "original". This argument is only used when arg. 'newdata' is missing, and is ignored otherwise. If equal to "estimation", the vector returned matches the sample used for the estimation. If equal to "original", it matches the original data set (the observations not used for the estimation being filled with NAs).
#' @param se.fit Logical, default is \code{FALSE}. If \code{TRUE}, the standard-error of the predicted value is computed and returned in a column named \code{se.fit}. This feature is only available for OLS models not containing fixed-effects.
#' @param interval Either "none" (default), "confidence" or "prediction". What type of confidence interval to compute. Note that this feature is only available for OLS models not containing fixed-effects (GLM/ML models are not covered).
#' @param level A numeric scalar in between 0.5 and 1, defaults to 0.95. Only used when the argument 'interval' is requested, it corresponds to the width of the confidence interval.
#' @param fixef Logical scalar, default is \code{FALSE}. If \code{TRUE}, a data.frame is returned, with each column representing the fixed-effects coefficients for each observation in \code{newdata} -- with as many columns as fixed-effects. Note that when there are variables with varying slopes, the slope coefficients are returned (i.e. they are not multiplied by the variable).
#' @param vs.coef Logical scalar, default is \code{FALSE}. Only used when \code{fixef = TRUE} and when variables with varying slopes are present. If \code{TRUE}, the coefficients of the variables with varying slopes are returned instead of the coefficient multiplied by the value of the variables (default).
#' @param ... Not currently used.
#'
#'
#' @return
#' It returns a numeric vector of length equal to the number of observations in argument \code{newdata}.
#' If \code{newdata} is missing, it returns a vector of the same length as the estimation sample, except if \code{sample = "original"}, in which case the length of the vector will match the one of the original data set (which can, but also cannot, be the estimation sample).
#' If \code{fixef = TRUE}, a \code{data.frame} is returned.
#' If \code{se.fit = TRUE} or \code{interval != "none"}, the object returned is a data.frame with the following columns: \code{fit}, \code{se.fit}, and, if CIs are requested, \code{ci_low} and \code{ci_high}.
#'
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{update.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{fixef.fixest}}.
#'
#' @examples
#'
#' # Estimation on iris data
#' res = fepois(Sepal.Length ~ Petal.Length | Species, iris)
#'
#' # what would be the prediction if the data was all setosa?
#' newdata = data.frame(Petal.Length = iris$Petal.Length, Species = "setosa")
#' pred_setosa = predict(res, newdata = newdata)
#'
#' # Let's look at it graphically
#' plot(c(1, 7), c(3, 11), type = "n", xlab = "Petal.Length",
#'      ylab = "Sepal.Length")
#'
#' newdata = iris[order(iris$Petal.Length), ]
#' newdata$Species = "setosa"
#' lines(newdata$Petal.Length, predict(res, newdata))
#'
#' # versicolor
#' newdata$Species = "versicolor"
#' lines(newdata$Petal.Length, predict(res, newdata), col=2)
#'
#' # virginica
#' newdata$Species = "virginica"
#' lines(newdata$Petal.Length, predict(res, newdata), col=3)
#'
#' # The original data
#' points(iris$Petal.Length, iris$Sepal.Length, col = iris$Species, pch = 18)
#' legend("topleft", lty = 1, col = 1:3, legend = levels(iris$Species))
#'
#'
#' #
#' # Getting the fixed-effect coefficients for each obs.
#' #
#'
#' data(trade)
#' est_trade = fepois(Euros ~ log(dist_km) | Destination^Product +
#'                                            Origin^Product + Year, trade)
#' obs_fe = predict(est_trade, fixef = TRUE)
#' head(obs_fe)
#'
#' # can we check we get the right sum of fixed-effects
#' head(cbind(rowSums(obs_fe), est_trade$sumFE))
#'
#'
#' #
#' # Standard-error of the prediction
#' #
#'
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#'
#' est = feols(y ~ x1 + species, base)
#'
#' head(predict(est, se.fit = TRUE))
#'
#' # regular confidence interval
#' head(predict(est, interval = "conf"))
#'
#' # adding the residual to the CI
#' head(predict(est, interval = "predi"))
#'
#' # You can change the type of SE on the fly
#' head(predict(est, interval = "conf", vcov = ~species))
#'
#'
#'
predict.fixest = function(object, newdata, type = c("response", "link"), se.fit = FALSE,
                          interval = "none", level = 0.95, fixef = FALSE,
                          vs.coef = FALSE, sample = c("estimation", "original"),
                          vcov = NULL, ssc = NULL, ...){

    # Checking the arguments
    if(is_user_level_call()){
        validate_dots(suggest_args = c("newdata", "type"))
    }

    # Controls
    check_arg_plus(type, sample, "match")
    check_arg(fixef, vs.coef, "logical scalar")
    check_arg(se.fit, "logical scalar")
    check_arg(level, "numeric scalar GT{.50} LT{1}")
    check_arg_plus(interval, "match(none, confidence, prediction)")
    if(!se.fit && interval != "none"){
        se.fit = TRUE
    }

    if(se.fit && object$method_type != "feols"){
        stop("The standard-error of the prediction is currently only available for OLS models, sorry.")
    }

    # renaming to clarify
    fixef.return = fixef
    do_anyway = fixef.return || se.fit

    # if newdata is missing
    is_original_data = FALSE
    if(missing(newdata)){

        if(do_anyway || isTRUE(object$lean)){
            newdata = fetch_data(object, "In 'predict', ")
            is_original_data = TRUE
        } else {
            if(type == "response" || object$method_type == "feols"){
                res = object$fitted.values
            } else if(object$method == "femlm") {
                if("mu" %in% names(object)){
                    res = object$mu
                } else {
                    family = object$family
                    famFuns = switch(family,
                                     poisson = ml_poisson(),
                                     negbin = ml_negbin(),
                                     logit = ml_logit(),
                                     gaussian = ml_gaussian())

                    res = famFuns$linearFromExpected(object$fitted.values)
                }
            } else {
                res = object$family$linkfun(object$fitted.values)
            }

            if(sample == "original") res = fill_with_na(res, object)

            return(res)
        }

    }

    if(!is.matrix(newdata) && !"data.frame" %in% class(newdata)){
        stop("Argument 'newdata' must be a data.frame.")
    }

    # we ensure it really is a clean data.frame
    newdata = as.data.frame(newdata)

    mc = match.call()
    if(fixef.return){

        if(is.null(object$fixef_vars)){
            stop("The argument 'fixef=TRUE' cannot work since the estimation did not contain fixed-effects.")
        }

        if("type" %in% names(mc)){
            warning("Argument 'type' is ignored when fixef = TRUE.")
        }
    }

    # We deconstruct it in four steps:
    # 1) cluster
    # 2) linear
    # 3) non-linear
    # 4) offset

    # + step 0: panel setup

    n = nrow(newdata)

    # NOTA 2019-11-26: I'm pondering whether to include NA-related messages
    # (would it be useful???)


    # STEP 0: panel setup

    fml = object$fml
    panel__meta__info = set_panel_meta_info(object, newdata)

    #
    # 1) Fixed-effects
    #

    # init fixed-effect values
    value_fixef = 0

    fixef_vars = object$fixef_vars
    if(!is.null(fixef_vars)){

        n_fe = length(fixef_vars)

        # Extraction of the FEs
        id_fixef = list()
        for(i in 1:n_fe){
            # checking if the variable is in the newdata
            fe_var = fixef_vars[i]
            variable = all.vars(str2lang(fe_var))
            isNotHere = !variable %in% names(newdata)
            if(any(isNotHere)){
                stop("The variable ", variable[isNotHere][1], " is absent from the 'newdata' but is needed for prediction (it is a fixed-effect variable).")
            }

            # The values taken by the FE variable
            fixef_values_possible = attr(object$fixef_id[[i]], "fixef_names")

            # Checking if ^ is present
            if(grepl("\\^", fe_var)){
                # If fastCombine was used => we're screwed, impossible to recover
                if(!all(grepl("_", fixef_values_possible, fixed = TRUE))){
                    stop("You cannot use predict() based on the initial regression since the fixed-effect '", fe_var, "' was combined using an algorithm dropping the FE values (but fast). Please re-run the regression using the argument 'combine.quick=FALSE'.")
                }

                fe_var = fml_combine(fe_var, fastCombine = FALSE, vars = TRUE)
            }

            # Obtaining the vector of fixed-effect
            fixef_current = eval(str2lang(fe_var), newdata)

            fixef_current_num = unclass(factor(fixef_current, levels = fixef_values_possible))
            id_fixef[[i]] = fixef_current_num
        }

        names(id_fixef) = fixef_vars

        # Value of the fixef coefficients // we don't show the notes, it's inelegant
        fixef_coef = fixef(object, sorted = FALSE, notes = FALSE)

        # We create the DF to be returned
        if(fixef.return){
            fixef_df = list()
        }

        # Adding the FEs and Slopes
        if(!is.null(object$fixef_terms)){

            terms_full = extract_fe_slope(object$fixef_terms)
            fixef_vars = terms_full$fixef_vars
            slope_fe = terms_full$slope_fe
            slope_vars = terms_full$slope_vars
            slope_terms = terms_full$slope_terms

            # We extract the slope variables
            slope_vars_unik = unique(slope_vars)

            slope_var_list = list()
            for(i in 1:length(slope_vars_unik)){
                variable = all.vars(str2lang(slope_vars_unik[i]))
                isNotHere = !variable %in% names(newdata)
                if(any(isNotHere)){
                    stop("The variable ", variable[isNotHere][1], " is absent from the 'newdata' but is needed for prediction (it is a variable with varying slope).")
                }

                slope_var_list[[slope_vars_unik[i]]] = eval(str2lang(slope_vars_unik[i]), newdata)
            }

            # Adding the FE values
            for(var in fixef_vars){
                fixef_current_num = id_fixef[[var]]
                fixef_coef_current = fixef_coef[[var]]

                if(fixef.return){
                    fixef_df[[var]] = fixef_coef_current[fixef_current_num]

                } else {
                    value_fixef = value_fixef + fixef_coef_current[fixef_current_num]
                }

            }

            # Adding the slopes
            for(i in seq_along(slope_vars)){

                fixef_current_num = id_fixef[[slope_fe[i]]]
                fixef_coef_current = fixef_coef[[slope_terms[i]]]

                if(fixef.return){
                    vname = slope_terms[i]

                    # We return only the coefs OR the coef * the variable
                    if(vs.coef){
                        fixef_df[[vname]] = fixef_coef_current[fixef_current_num]

                    } else {
                        fixef_df[[vname]] = fixef_coef_current[fixef_current_num] * slope_var_list[[slope_vars[i]]]
                    }


                } else {
                    value_fixef = value_fixef + fixef_coef_current[fixef_current_num] * slope_var_list[[slope_vars[i]]]
                }
            }


        } else {
            # Adding only FEs
            for(i in 1:n_fe){
                fixef_current_num = id_fixef[[i]]
                fixef_coef_current = fixef_coef[[i]]

                if(fixef.return){
                    fixef_df[[fixef_vars[i]]] = fixef_coef_current[fixef_current_num]

                } else {
                    value_fixef = value_fixef + fixef_coef_current[fixef_current_num]
                }

            }
        }

        if(fixef.return){

            # putting the results into a DF
            res = fixef_df
            attr(res, "row.names") = .set_row_names(length(res[[1L]]))
            oldClass(res) = "data.frame"

            if(is_original_data && sample == "estimation"){
                # here we want the same nber of obs
                # as in the estimation sample
                for(i in seq_along(object$obs_selection)){
                    res = res[object$obs_selection[[i]], , drop = FALSE]
                }
            }

            return(res)
        }

        # dropping names
        value_fixef = as.vector(value_fixef)
    }

    #
    # 2) Linear values
    #

    coef = object$coefficients

    value_linear = 0
    var_keep = NULL
    rhs_fml = fml_split(fml, 1)
    linear.varnames = all_vars_with_i_prefix(rhs_fml[[3]])

    if(length(linear.varnames) > 0){
        # Checking all variables are there

        if(isTRUE(object$iv) && object$iv_stage == 2){
            names(coef) = gsub("^fit_", "", names(coef))
            linear.varnames = c(linear.varnames, all_vars_with_i_prefix(object$fml_all$iv[[2]]))
            iv_fml = object$fml_all$iv
            rhs_fml = .xpd(..lhs ~ ..endo + ..rhs, ..lhs = rhs_fml[[2]], ..endo = iv_fml[[2]], ..rhs = rhs_fml[[3]])
        }

        varNotHere = setdiff(linear.varnames, names(newdata))
        if(length(varNotHere) > 0){
            stop("The variable", enumerate_items(varNotHere, "s.quote"),
                 " used to estimate the model (in fml) ", ifsingle(varNotHere, "is", "are"),
                 " missing in the data.frame given by the argument 'newdata'.")
        }

        # we create the matrix
        matrix_linear = error_sender(fixest_model_matrix_extra(object = object, newdata = newdata, original_data = FALSE, fml = rhs_fml, i_noref = TRUE), "Error when creating the linear matrix: ")

        # Checking the levels created with i()
        mm_info_new = attr(matrix_linear, "model_matrix_info")
        if(!is.null(mm_info_new)){
            mm_info = object$model_matrix_info
            # The order of creation is exactly the same (same fun used),
            # so the two mm_info are identical in structure
            for(i in seq_along(mm_info)){
                mm_new_i = mm_info_new[[i]]
                mm_i = mm_info[[i]]
                if("coef_names_full" %in% names(mm_i)){
                    pblm = setdiff(mm_new_i$coef_names_full, mm_i$coef_names_full)
                    if(length(pblm) > 0){
                        stop(dsb("In i(), predictions cannot be done for values that were not present at estimation time.",
                                 " It concerns the value.[*s_, 3KO, C?pblm]."))
                    }
                }
            }
        }

        var_keep = intersect(names(coef), colnames(matrix_linear))
        value_linear = value_linear + as.vector(matrix_linear[, var_keep, drop = FALSE] %*% coef[var_keep])
    }

    #
    # 3) Non linear terms
    #

    value_NL = 0
    NL_fml = object$NL.fml
    if(!is.null(NL_fml)){
        # controlling that we can evaluate that
        NL_vars = all.vars(NL_fml)
        varNotHere = setdiff(NL_vars, c(names(coef), names(newdata)))
        if(length(varNotHere) > 0){
            stop("Some variables used to estimate the model (in the non-linear formula) are missing from argument 'newdata': ", enumerate_items(varNotHere), ".")
        }

        var2send = intersect(NL_vars, names(newdata))
        env = new.env()
        for(var in var2send){
            assign(var, newdata[[var]], env)
        }

        coef2send = setdiff(NL_vars, names(newdata))
        for(iter_coef in coef2send){
            assign(iter_coef, coef[iter_coef], env)
        }

        # Evaluation of the NL part
        value_NL = eval(NL_fml[[2]], env)
    }

    #
    # 4) offset value
    #

    value_offset = 0
    offset = object$call$offset
    if(!is.null(offset)){
        # evaluation of the offset

        if(is.numeric(offset)){
            # simple numeric offset
            value_offset = offset

        } else {
            # offset valid only if formula
            offset_char = as.character(offset)

            if(length(offset_char) == 2 && offset_char[1] == "~"){
                offset_fml = eval(offset)
                varNotHere = setdiff(all.vars(offset_fml), names(newdata))
                if(length(varNotHere) > 0){
                    stop("In the offset, the variable", enumerate_items(varNotHere, "s.is"), " not present in 'newdata'.")
                }

                value_offset = eval(offset_fml[[length(offset_fml)]], newdata)
            } else {
                stop("Predict can't be applied to this estimation because the offset (", deparse_long(offset), ") cannot be evaluated for the new data. Use a formula for the offset in the first estimation to avoid this.")
            }

        }

    }

    value_predicted = value_fixef + value_linear + value_NL + value_offset

    if(type == "link" || object$method_type == "feols"){
        res = value_predicted
    } else if(object$method == "femlm") {
        # Now the expected predictor
        family = object$family
        famFuns = switch(family,
                         poisson = ml_poisson(),
                         negbin = ml_negbin(),
                         logit = ml_logit(),
                         gaussian = ml_gaussian())

        if(family == "gaussian"){
            exp_value = 0
        } else {
            exp_value = exp(value_predicted)
        }

        res = famFuns$expected.predictor(value_predicted, exp_value)
    } else {
        res = object$family$linkinv(value_predicted)
    }


    #
    # se.fit
    #

    if(se.fit){

        if(!is.null(object$fixef_vars)){
            stop("The standard-errors (SEs) of the prediction cannot be computed in the presence of fixed-effects. To obtain the SEs, you would need to include the FEs as standard factors in the model.")
        }

        if(!is.null(NL_fml)){
            stop("The standard-errors (SEs) of the prediction cannot be computed in models containing non-linear in parameter elements.")
        }

        # The matrix has been created already

        V_raw = vcov(object, attr = TRUE, vcov = vcov, ssc = ssc)
        V = V_raw[var_keep, var_keep, drop = FALSE]
        X = matrix_linear[, var_keep, drop = FALSE]

        var.fit = rowSums((X %*% V) * X)
        se.fit = sqrt(var.fit)

        res = data.frame(fit = res, se.fit = se.fit)

        if(interval != "none"){
            fact = fixest_CI_factor(object, level, V_raw)

            if(interval == "prediction"){
                w = object$weights

                if(!is.null(w)){
                    var.u = cpp_ssq(resid(object), w) / w
                } else {
                    var.u = cpp_ssq(resid(object))
                }
                var.u = var.u / degrees_freedom(object, "resid")

                se_obs = sqrt(var.fit + var.u)

            } else {
                se_obs = se.fit
            }

            res$ci_low  = res$fit + fact[1] * se_obs
            res$ci_high = res$fit + fact[2] * se_obs
        }

    }

    res
}


#' Confidence interval for parameters estimated with \code{fixest}
#'
#' This function computes the confidence interval of parameter estimates obtained from a model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#' @inheritParams vcov.fixest
#'
#' @param parm The parameters for which to compute the confidence interval (either an integer vector OR a character vector with the parameter name). If missing, all parameters are used.
#' @param level The confidence level. Default is 0.95.
#'
#' @return
#' Returns a data.frame with two columns giving respectively the lower and upper bound of the confidence interval. There is as many rows as parameters.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade (with 3 fixed-effects)
#' est_pois = femlm(Euros ~ log(dist_km) + log(Year) | Origin + Destination +
#'                  Product, trade)
#'
#' # confidence interval with "normal" VCOV
#' confint(est_pois)
#'
#' # confidence interval with "clustered" VCOV (w.r.t. the Origin factor)
#' confint(est_pois, se = "cluster")
#'
#'
confint.fixest = function(object, parm, level = 0.95, vcov, se, cluster, ssc = NULL, ...){

    # Checking the arguments
    if(is_user_level_call()){
        validate_dots(suggest_args = c("parm", "level", "se", "cluster"),
                      valid_args = c("forceCovariance", "keepBounded"))
    }

    # Control
    if(!is.numeric(level) || !length(level) == 1 || level >= 1 || level <= .50){
        stop("The argument 'level' must be a numeric scalar greater than 0.50 and strictly lower than 1.")
    }

    # The proper SE
    sum_object = summary(object, vcov = vcov, se = se, cluster = cluster, ssc = ssc, ...)

    se_all = sum_object$coeftable[, 2]
    coef_all = sum_object$coeftable[, 1]

    # the parameters for which we should compute the confint
    all_params = names(coef_all)

    if(missing(parm)){
        parm_use = all_params
    } else if(is.numeric(parm)){
        if(any(parm %% 1 != 0)){
            stop("If the argument 'parm' is numeric, it must be integers.")
        }

        parm_use = unique(na.omit(all_params[parm]))
        if(length(parm_use) == 0){
            stop("There are ", length(all_params), " coefficients, the argument 'parm' does not correspond to any of them.")
        }
    } else if(is.character(parm)){
        parm_pblm = setdiff(parm, all_params)
        if(length(parm_pblm) > 0){
            stop("some parameters of 'parm' have no estimated coefficient: ", paste0(parm_pblm, collapse=", "), ".")
        }

        parm_use = intersect(parm, all_params)
    }

    # multiplicative factor
    fact = fixest_CI_factor(object, level, sum_object$cov.scaled)

    # The confints
    # Note that for glm models, there is no profiling
    lower_bound = coef_all[parm_use] + fact[1] * se_all[parm_use]
    upper_bound = coef_all[parm_use] + fact[2] * se_all[parm_use]

    res = data.frame(lower_bound, upper_bound, row.names = parm_use)

    val = (1 - level) / 2
    names(res) = paste0(round(100*c(val, 1-val), 1), " %")

    attr(res, "type") = attr(se_all, "type")

    res
}

#' Updates a \code{fixest} estimation
#'
#' Updates and re-estimates a \code{fixest} model (estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). This function updates the formulas and use previous starting values to estimate a new \code{fixest} model. The data is obtained from the original \code{call}.
#'
#' @method update fixest
#'
#' @inheritParams nobs.fixest
#'
#' @param fml.update Changes to be made to the original argument \code{fml}. See more information on \code{\link[stats]{update.formula}}. You can add/withdraw both variables and fixed-effects. E.g. \code{. ~ . + x2 | . + z2} would add the variable \code{x2} and the cluster \code{z2} to the former estimation.
#' @param nframes (Advanced users.) Defaults to 1. Number of frames up the stack where to perform the evaluation of the updated call. By default, this is the parent frame.
#' @param evaluate Logical, default is \code{TRUE}. If \code{FALSE}, only the updated call is returned.
#' @param ... Other arguments to be passed to the functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @return
#' It returns a \code{fixest} object (see details in \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}).
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{predict.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{fixef.fixest}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Example using trade data
#' data(trade)
#'
#' # main estimation
#' est_pois = fepois(Euros ~ log(dist_km) | Origin + Destination, trade)
#'
#' # we add the variable log(Year)
#' est_2 = update(est_pois, . ~ . + log(Year))
#'
#' # we add another fixed-effect: "Product"
#' est_3 = update(est_2, . ~ . | . + Product)
#'
#' # we remove the fixed-effect "Origin" and the variable log(dist_km)
#' est_4 = update(est_3, . ~ . - log(dist_km) | . - Origin)
#'
#' # Quick look at the 4 estimations
#' etable(est_pois, est_2, est_3, est_4)
#'
update.fixest = function(object, fml.update, nframes = 1, evaluate = TRUE, ...){
    # Update method
    # fml.update: update the formula
    # If 1) SAME DATA and 2) SAME dep.var, then we make initialisation


    if(missing(fml.update)){
        fml.update = . ~ .
    } else {
        check_arg(fml.update, "formula")
    }

    check_arg(evaluate, "logical scalar")

    if(isTRUE(object$is_fit)){
        stop("update method not available for fixest estimations obtained from fit methods.")
    }

    if(!isScalar(nframes) || nframes < 1 || nframes %% 1 != 0){
        stop("Argument 'nframes' must be a single integer greater than, or equal to, 1.")
    }

    call_new = match.call()
    dots = list(...)

    dot_names = names(dots)
    if("fixef" %in% dot_names){
        stop("Argument 'fixef' is not accepted in the 'update.fixest' method. Please make modifications to fixed-effects directly in the argument 'fml.update'. (E.g. .~.|.+v5 to add variable v5 as a fixed-effect.)")
    }

    if(any(dot_names == "")){
        call_new_names = names(call_new)
        problems = call_new[call_new_names == ""][-1]
        stop("In 'update.fixest' the arguments of '...' are passed to the function ", object$method, ", and must be named. Currently there are un-named arguments (e.g. '", deparse_long(problems[[1]]), "').")
    }

    #
    # I) Linear formula update
    #

    fml_old = object$fml
    fml_linear = update(fml_old, fml_split(fml.update, 1))

    # Family information
    if(!is.null(dots$family)){
        if(object$method_type == "feols"){
            stop("'family' is not an argument of function feols().")
        } else if(object$method %in% c("femlm", "feNmlm", "fepois", "fenegbin")){
            family_new = match.arg(dots$family, c("poisson", "negbin", "gaussian", "logit"))
        }
    }

    #
    # II) fixed-effects updates
    #

    fml_fixef = NULL

    updt_fml_parts = fml_split(fml.update, raw = TRUE)
    n_parts = length(updt_fml_parts)

    if(n_parts > 2 + (object$method_type == "feols")){
        stop("The update formula cannot have more than ", 2 + (object$method_type == "feols"), " parts for the method ", object$method, ".")
    }

    is_fe = n_parts > 1 && !is_fml_inside(updt_fml_parts[[2]])

    fixef_vars = object$fixef_vars

    if(is_fe){

        fixef_old = object$fml_all$fixef

        # I use it as text to catch the var1^var2 FEs (update does not work)
        if(is.null(fixef_old)){
            fixef_old_text = "~ 1"
        } else {
            fixef_old_text = deparse_long(fixef_old)
        }

        fixef_new_fml = fml_maker(updt_fml_parts[[2]])
        fixef_new_text = deparse_long(fixef_new_fml)

        if(fixef_new_text == "~."){
            # nothing happens
            fixef_new = fixef_old

        } else if(fixef_new_text %in% c("~0", "~1")){
            fixef_new = ~1

        } else if(grepl("\\^", fixef_old_text) || grepl("\\^", fixef_new_text)){
            # we update manually.... dammmit
            # Note that what follows does not work ONLY when you have number^var or number^number
            # and both cases don't make much sense -- I need not control for them
            fml_text_old = gsub("\\^", "__666__", fixef_old_text)
            fml_text_new = gsub("\\^", "__666__", fixef_new_text)

            fixef_new_wip = update(as.formula(fml_text_old), as.formula(fml_text_new))

            fixef_new = as.formula(gsub("__666__", "^", fixef_new_wip))
        } else {
            fixef_new = update(fixef_old, fixef_new_fml)
        }

        if(length(all.vars(fixef_new)) > 0){
            # means there is a fixed-effect
            fml_fixef = fixef_new
        }

    } else if(!is.null(fixef_vars)){
        # the formula updated:
        fml_fixef = object$fml_all$fixef

    }

    #
    # III) IV updates
    #

    if(n_parts > 2 || (n_parts == 2 && !is_fe)){

        iv_new_fml = fml_maker(updt_fml_parts[[n_parts]])

        if(!is_fml_inside(iv_new_fml)){
            stop("The third part of the update formula in 'feols' must be a formula.")
        }

        iv_old = object$fml_all$iv

        if(is.null(iv_old)){
            fml_iv = iv_new_fml

        } else {
            fml_iv = update(iv_old, iv_new_fml)
        }

    } else {
        fml_iv = object$fml_all$iv
    }


    fml_new = merge_fml(fml_linear, fml_fixef, fml_iv)


    #
    # The call
    #

    call_old = object$call

    # we drop the argument fixef from old call (now it's in the fml_new)
    call_old$fixef = NULL

    # We also drop the arguments for multiple estimations:
    call_old$split = call_old$fsplit = NULL

    # new call: call_clear
    call_clear = call_old
    for(arg in setdiff(names(call_new)[-1], c("fml.update", "nframes", "evaluate", "object"))){
        call_clear[[arg]] = call_new[[arg]]
    }

    call_clear$fml = as.call(fml_new)

    if(!evaluate) return(call_clear)

    res = eval(call_clear, parent.frame(nframes))

    res
}


#' Extract the formula of a \code{fixest} fit
#'
#' This function extracts the formula from a \code{fixest} estimation (obtained with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). If the estimation was done with fixed-effects, they are added in the formula after a pipe (\dQuote{|}). If the estimation was done with a non linear in parameters part, then this will be added in the formula in between \code{I()}.
#'
#'
#' @param x An object of class \code{fixest}. Typically the result of a \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}} estimation.
#' @param type A character scalar. Default is \code{type = "full"} which gives back a formula containing the linear part of the model along with the fixed-effects (if any) and the IV part (if any). If \code{type = "linear"} then only the linear formula is returned. If \code{type = "NL"} then only the non linear in parameters part is returned.
#' @param ... Not currently used.
#'
#' @return
#' It returns a formula.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{model.matrix.fixest}}, \code{\link[fixest]{update.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # formula with the fixed-effect variable
#' formula(res)
#'
#' # linear part without the fixed-effects
#' formula(res, "linear")
#'
#'
formula.fixest = function(x, type = c("full", "linear", "iv", "NL"), ...){
    # Extract the formula from the object
    # we add the clusters in the formula if needed

    # Checking the arguments
    if(is_user_level_call()){
        validate_dots(suggest_args = "type")
    }

    if(isTRUE(x$is_fit)){
        stop("formula method not available for fixest estimations obtained from fit methods.")
    }

    check_arg_plus(type, "match")

    if(type == "linear"){
        return(x$fml)

    } else if(type == "NL"){

        if(!x$method == "feNmlm"){
            stop("type = 'NL' is not valid for a ", x$method, " estimation.")
        }

        NL.fml = x$NL.fml
        if(is.null(NL.fml)){
            stop("There was no nonlinear part estimated, option type = 'NL' cannot be used.")
        }

        return(NL.fml)

    } else if(type == "iv"){
        if(is.null(x$fml_all$iv)){
            stop("type = 'iv' is only available for feols estimations with IV.")
        }
    }

    # Shall I add LHS ~ RHS + NL(NL fml) | fe | iv ???
    res = merge_fml(x$fml_all$linear, x$fml_all$fixef, x$fml_all$iv)

    res
}

#' Extract the terms
#'
#' This function extracts the terms of a \code{fixest} estimation, excluding the fixed-effects part.
#'
#' @param x A \code{fixest} object. Obtained using the functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param ... Not currently used.
#'
#' @return
#' An object of class \code{c("terms", "formula")} which contains the terms representation of a symbolic model.
#'
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res = feols(Sepal.Length ~ Sepal.Width*Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # Terms of the linear part
#' terms(res)
#'
#'
terms.fixest = function(x, ...){
    terms(formula(x, type = "linear"))
}
