#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Tue Apr 23 16:41:47 2019
# Purpose: All estimation functions
#----------------------------------------------#



#' Fixed-effects OLS estimation
#'
#' Estimates OLS with any number of fixed-effects.
#'
#' @inheritParams femlm
#' @inheritSection xpd Dot square bracket operator in formulas
#'
#' @param fml A formula representing the relation to be estimated. For example: \code{fml = z~x+y}. To include fixed-effects, insert them in this formula using a pipe: e.g. \code{fml = z~x+y | fe_1+fe_2}. You can combine two fixed-effects with \code{^}: e.g. \code{fml = z~x+y|fe_1^fe_2}, see details. You can also use variables with varying slopes using square brackets: e.g. in \code{fml = z~y|fe_1[x] + fe_2}, see details. To add IVs, insert the endogenous vars./instruments after a pipe, like in \code{y ~ x | c(x_endo1, x_endo2) ~ x_inst1 + x_inst2}. Note that it should always be the last element, see details. Multiple estimations can be performed at once: for multiple dep. vars, wrap them in \code{c()}: ex \code{c(y1, y2)}. For multiple indep. vars, use the stepwise functions: ex \code{x1 + csw(x2, x3)}. The formula \code{fml = c(y1, y2) ~ x1 + cw0(x2, x3)} leads to 6 estimation, see details. Square brackets starting with a dot can be used to call global variables: \code{y.[i] ~ x.[1:2]} will lead to \code{y3 ~ x1 + x2} if \code{i} is equal to 3 in the current environment (see details in \code{\link[fixest2]{xpd}}).
#' @param weights A formula or a numeric vector. Each observation can be weighted, the weights must be greater than 0. If equal to a formula, it should be one-sided: for example \code{~ var_weight}.
#' @param verbose Integer. Higher values give more information. In particular, it can detail the number of iterations in the demeaning algorithm (the first number is the left-hand-side, the other numbers are the right-hand-side variables).
#' @param demeaned Logical, default is \code{FALSE}. Only used in the presence of fixed-effects: should the centered variables be returned? If \code{TRUE}, it creates the items \code{y_demeaned} and \code{X_demeaned}.
#' @param notes Logical. By default, two notes are displayed: when NAs are removed (to show additional information) and when some observations are removed because of collinearity. To avoid displaying these messages, you can set \code{notes = FALSE}. You can remove these messages permanently by using \code{setFixest_notes(FALSE)}.
#' @param collin.tol Numeric scalar, default is \code{1e-10}. Threshold deciding when variables should be considered collinear and subsequently removed from the estimation. Higher values means more variables will be removed (if there is presence of collinearity). One signal of presence of collinearity is t-stats that are extremely low (for instance when t-stats < 1e-3).
#' @param y Numeric vector/matrix/data.frame of the dependent variable(s). Multiple dependent variables will return a \code{fixest_multi} object.
#' @param X Numeric matrix of the regressors.
#' @param fixef_df Matrix/data.frame of the fixed-effects.
#'
#' @details
#' The method used to demean each variable along the fixed-effects is based on Berge (2018), since this is the same problem to solve as for the Gaussian case in a ML setup.
#'
#' @section Combining the fixed-effects:
#' You can combine two variables to make it a new fixed-effect using \code{^}. The syntax is as follows: \code{fe_1^fe_2}. Here you created a new variable which is the combination of the two variables fe_1 and fe_2. This is identical to doing \code{paste0(fe_1, "_", fe_2)} but more convenient.
#'
#' Note that pasting is a costly operation, especially for large data sets. Thus, the internal algorithm uses a numerical trick which is fast, but the drawback is that the identity of each observation is lost (i.e. they are now equal to a meaningless number instead of being equal to \code{paste0(fe_1, "_", fe_2)}). These \dQuote{identities} are useful only if you're interested in the value of the fixed-effects (that you can extract with \code{\link[fixest2]{fixef.fixest}}). If you're only interested in coefficients of the variables, it doesn't matter. Anyway, you can use \code{combine.quick = FALSE} to tell the internal algorithm to use \code{paste} instead of the numerical trick. By default, the numerical trick is performed only for large data sets.
#'
#' @section Varying slopes:
#' You can add variables with varying slopes in the fixed-effect part of the formula. The syntax is as follows: fixef_var[var1, var2]. Here the variables var1 and var2 will be with varying slopes (one slope per value in fixef_var) and the fixed-effect fixef_var will also be added.
#'
#' To add only the variables with varying slopes and not the fixed-effect, use double square brackets: fixef_var[[var1, var2]].
#'
#' In other words:
#' \itemize{
#'   \item fixef_var[var1, var2] is equivalent to fixef_var + fixef_var[[var1]] + fixef_var[[var2]]
#'   \item fixef_var[[var1, var2]] is equivalent to fixef_var[[var1]] + fixef_var[[var2]]
#' }
#'
#' In general, for convergence reasons, it is recommended to always add the fixed-effect and avoid using only the variable with varying slope (i.e. use single square brackets).
#'
#' @section Lagging variables:
#'
#' To use leads/lags of variables in the estimation, you can: i) either provide the argument \code{panel.id}, ii) either set your data set as a panel with the function \code{\link[fixest2]{panel}}. Doing either of the two will give you acceess to the lagging functions \code{\link[fixest2]{l}},  \code{\link[fixest2:l]{f}} and \code{\link[fixest2:l]{d}}.
#'
#' You can provide several leads/lags/differences at once: e.g. if your formula is equal to \code{f(y) ~ l(x, -1:1)}, it means that the dependent variable is equal to the lead of \code{y}, and you will have as explanatory variables the lead of \code{x1}, \code{x1} and the lag of \code{x1}. See the examples in function \code{\link[fixest2]{l}} for more details.
#'
#' @section Interactions:
#'
#' You can interact a numeric variable with a "factor-like" variable by using \code{i(factor_var, continuous_var, ref)}, where \code{continuous_var} will be interacted with each value of \code{factor_var} and the argument \code{ref} is a value of \code{factor_var} taken as a reference (optional).
#'
#' Using this specific way to create interactions leads to a different display of the interacted values in \code{\link[fixest2]{etable}} and offers a special representation of the interacted coefficients in the function \code{\link[fixest2]{coefplot}}. See examples.
#'
#'  It is important to note that *if you do not care about the standard-errors of the interactions*, then you can add interactions in the fixed-effects part of the formula, it will be incomparably faster (using the syntax \code{factor_var[continuous_var]}, as explained in the section \dQuote{Varying slopes}).
#'
#' The function \code{\link[fixest2:i]{i}} has in fact more arguments, please see details in its associated help page.
#'
#' @section On standard-errors:
#'
#' Standard-errors can be computed in different ways, you can use the arguments \code{se} and \code{ssc} in \code{\link[fixest2]{summary.fixest}} to define how to compute them. By default, in the presence of fixed-effects, standard-errors are automatically clustered.
#'
#' The following vignette: \href{https://lrberge.github.io/fixest/articles/standard_errors.html}{On standard-errors} describes in details how the standard-errors are computed in \code{fixest} and how you can replicate standard-errors from other software.
#'
#' You can use the functions \code{\link[fixest2]{setFixest_vcov}} and \code{\link[fixest2:ssc]{setFixest_ssc}} to permanently set the way the standard-errors are computed.
#'
#' @section Instrumental variables:
#'
#' To estimate two stage least square regressions, insert the relationship between the endogenous regressor(s) and the instruments in a formula, after a pipe.
#'
#' For example, \code{fml = y ~ x1 | x_endo ~ x_inst} will use the variables \code{x1} and \code{x_inst} in the first stage to explain \code{x_endo}. Then will use the fitted value of \code{x_endo} (which will be named \code{fit_x_endo}) and \code{x1} to explain \code{y}.
#' To include several endogenous regressors, just use "+", like in: \code{fml = y ~ x1 | x_endo1 + x_end2 ~ x_inst1 + x_inst2}.
#'
#' Of course you can still add the fixed-effects, but the IV formula must always come last, like in \code{fml = y ~ x1 | fe1 + fe2 | x_endo ~ x_inst}.
#'
#' If you want to estimate a model without exogenous variables, use \code{"1"} as a placeholder: e.g. \code{fml = y ~ 1 | x_endo + x_inst}.
#'
#' By default, the second stage regression is returned. You can access the first stage(s) regressions either directly in the slot \code{iv_first_stage} (not recommended), or using the argument \code{stage = 1} from the function \code{\link[fixest2]{summary.fixest}}. For example \code{summary(iv_est, stage = 1)} will give the first stage(s). Note that using summary you can display both the second and first stages at the same time using, e.g., \code{stage = 1:2} (using \code{2:1} would reverse the order).
#'
#'
#' @section Multiple estimations:
#'
#' Multiple estimations can be performed at once, they just have to be specified in the formula. Multiple estimations yield a \code{fixest_multi} object which is \sQuote{kind of} a list of all the results but includes specific methods to access the results in a handy way. Please have a look at the dedicated vignette: \href{https://lrberge.github.io/fixest/articles/multiple_estimations.html}{Multiple estimations}.
#'
#' To include multiple dependent variables, wrap them in \code{c()} (\code{list()} also works). For instance \code{fml = c(y1, y2) ~ x1} would estimate the model \code{fml = y1 ~ x1} and then the model \code{fml = y2 ~ x1}.
#'
#' To include multiple independent variables, you need to use the stepwise functions. There are 4 stepwise functions: \code{sw}, \code{sw0}, \code{csw}, \code{csw0}, and \code{mvsw}. Of course \code{sw} stands for stepwise, and \code{csw} for cumulative stepwise. Finally \code{mvsw} is a bit special, it stands for multiverse stepwise. Let's explain that.
#' Assume you have the following formula: \code{fml = y ~ x1 + sw(x2, x3)}. The stepwise function \code{sw} will estimate the following two models: \code{y ~ x1 + x2} and \code{y ~ x1 + x3}. That is, each element in \code{sw()} is sequentially, and separately, added to the formula. Would have you used \code{sw0} in lieu of \code{sw}, then the model \code{y ~ x1} would also have been estimated. The \code{0} in the name means that the model without any stepwise element also needs to be estimated.
#' The prefix \code{c} means cumulative: each stepwise element is added to the next. That is, \code{fml = y ~ x1 + csw(x2, x3)} would lead to the following models \code{y ~ x1 + x2} and \code{y ~ x1 + x2 + x3}. The \code{0} has the same meaning and would also lead to the model without the stepwise elements to be estimated: in other words, \code{fml = y ~ x1 + csw0(x2, x3)} leads to the following three models: \code{y ~ x1}, \code{y ~ x1 + x2} and \code{y ~ x1 + x2 + x3}.
#' Finally \code{mvsw} will add, in a stepwise fashion all possible combinations of the variables in its arguments. For example \code{mvsw(x1, x2, x3)} is equivalent to \code{sw0(x1, x2, x3, x1 + x2, x1 + x3, x2 + x3, x1 + x2 + x3)}. The number of models to estimate grows at a factorial rate: so be cautious!
#'
#' Multiple independent variables can be combined with multiple dependent variables, as in \code{fml = c(y1, y2) ~ cw(x1, x2, x3)} which would lead to 6 estimations. Multiple estimations can also be combined to split samples (with the arguments \code{split}, \code{fsplit}).
#'
#' You can also add fixed-effects in a stepwise fashion. Note that you cannot perform stepwise estimations on the IV part of the formula (\code{feols} only).
#'
#' If NAs are present in the sample, to avoid too many messages, only NA removal concerning the variables common to all estimations is reported.
#'
#' A note on performance. The feature of multiple estimations has been highly optimized for \code{feols}, in particular in the presence of fixed-effects. It is faster to estimate multiple models using the formula rather than with a loop. For non-\code{feols} models using the formula is roughly similar to using a loop performance-wise.
#'
#' @section Tricks to estimate multiple LHS:
#'
#' To use multiple dependent variables in \code{fixest} estimations, you need to include them in a vector: like in \code{c(y1, y2, y3)}.
#'
#' First, if names are stored in a vector, they can readily be inserted in a formula to perform multiple estimations using the dot square bracket operator. For instance if \code{my_lhs = c("y1", "y2")}, calling \code{fixest} with, say \code{feols(.[my_lhs] ~ x1, etc)} is equivalent to using \code{feols(c(y1, y2) ~ x1, etc)}. Beware that this is a special feature unique to the \emph{left-hand-side} of \code{fixest} estimations (the default behavior of the DSB operator is to aggregate with sums, see \code{\link[fixest2]{xpd}}).
#'
#' Second, you can use a regular expression to grep the left-hand-sides on the fly. When the \code{..("regex")} feature is used naked on the LHS, the variables grepped are inserted into \code{c()}. For example \code{..("Pe") ~ Sepal.Length, iris} is equivalent to \code{c(Petal.Length, Petal.Width) ~ Sepal.Length, iris}. Beware that this is a special feature unique to the \emph{left-hand-side} of \code{fixest} estimations (the default behavior of \code{..("regex")} is to aggregate with sums, see \code{\link[fixest2]{xpd}}).
#'
#' @section Argument sliding:
#'
#' When the data set has been set up globally using \code{\link[fixest2]{setFixest_estimation}}\code{(data = data_set)}, the argument \code{vcov} can be used implicitly. This means that calls such as \code{feols(y ~ x, "HC1")}, or \code{feols(y ~ x, ~id)}, are valid: i) the data is automatically deduced from the global settings, and ii) the \code{vcov} is deduced to be the second argument.
#'
#' @section Piping:
#'
#' Although the argument 'data' is placed in second position, the data can be piped to the estimation functions. For example, with R >= 4.1, \code{mtcars |> feols(mpg ~ cyl)} works as \code{feols(mpg ~ cyl, mtcars)}.
#'
#'
#' @return
#' A \code{fixest} object. Note that \code{fixest} objects contain many elements and most of them are for internal use, they are presented here only for information. To access them, it is safer to use the user-level methods (e.g. \code{\link[fixest2]{vcov.fixest}}, \code{\link[fixest2]{resid.fixest}}, etc) or functions (like for instance \code{\link[fixest2]{fitstat}} to access any fit statistic).
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{fml_all}{A list containing different parts of the formula. Always contain the linear formula. Then depending on the cases: \code{fixef}: the fixed-effects, \code{iv}: the IV part of the formula.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{coefficients}{The named vector of estimated coefficients.}
#' \item{multicol}{Logical, if multicollinearity was found.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The loglikelihood.}
#' \item{ssr_null}{Sum of the squared residuals of the null model (containing only with the intercept).}
#' \item{ssr_fe_only}{Sum of the squared residuals of the model estimated with fixed-effects only.}
#' \item{ll_null}{The log-likelihood of the null model (containing only with the intercept).}
#' \item{ll_fe_only}{The log-likelihood of the model estimated with fixed-effects only.}
#' \item{fitted.values}{The fitted values.}
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
#' \item{obs_selection}{(When relevant.) List containing vectors of integers. It represents the sequential selection of observation vis a vis the original data set.}
#' \item{collin.var}{(When relevant.) Vector containing the variables removed because of collinearity.}
#' \item{collin.coef}{(When relevant.) Vector of coefficients, where the values of the variables removed because of collinearity are NA.}
#' \item{collin.min_norm}{The minimal diagonal value of the Cholesky decomposition. Small values indicate possible presence collinearity.}
#' \item{y_demeaned}{Only when \code{demeaned = TRUE}: the centered dependent variable.}
#' \item{X_demeaned}{Only when \code{demeaned = TRUE}: the centered explanatory variable.}
#'
#'
#' @seealso
#' See also \code{\link[fixest2]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest2]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest2]{etable}} to visualize the results of multiple estimations. For plotting coefficients: see \code{\link[fixest2]{coefplot}}.
#'
#' And other estimation methods: \code{\link[fixest2]{femlm}}, \code{\link[fixest2]{feglm}}, \code{\link[fixest2:feglm]{fepois}}, \code{\link[fixest2:femlm]{fenegbin}}, \code{\link[fixest2]{feNmlm}}.
#'
#' @author
#' Laurent Berge
#'
#' @references
#'
#' Berge, Laurent, 2018, "Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm." CREA Discussion Papers, 13 (\url{https://wwwen.uni.lu/content/download/110162/1299525/file/2018_13}).
#'
#' For models with multiple fixed-effects:
#'
#' Gaure, Simen, 2013, "OLS with multiple high dimensional category variables", Computational Statistics & Data Analysis 66 pp. 8--18
#'
#' @examples
#'
#' #
#' # Basic estimation
#' #
#'
#' res <- feols(Sepal.Length ~ Sepal.Width + Petal.Length, iris)
#' # You can specify clustered standard-errors in summary:
#' summary(res, cluster = ~Species)
#'
#' #
#' # Just one set of fixed-effects:
#' #
#'
#' res <- feols(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris)
#' # By default, the SEs are clustered according to the first fixed-effect
#' summary(res)
#'
#' #
#' # Varying slopes:
#' #
#'
#' res <- feols(Sepal.Length ~ Petal.Length | Species[Sepal.Width], iris)
#' summary(res)
#'
#' #
#' # Combining the FEs:
#' #
#'
#' base <- iris
#' base$fe_2 <- rep(1:10, 15)
#' res_comb <- feols(Sepal.Length ~ Petal.Length | Species^fe_2, base)
#' summary(res_comb)
#' fixef(res_comb)[[1]]
#'
#' #
#' # Using leads/lags:
#' #
#'
#' data(base_did)
#' # We need to set up the panel with the arg. panel.id
#' est1 <- feols(y ~ l(x1, 0:1), base_did, panel.id = ~ id + period)
#' est2 <- feols(f(y) ~ l(x1, -1:1), base_did, panel.id = ~ id + period)
#' etable(est1, est2, order = "f", drop = "Int")
#'
#' #
#' # Using interactions:
#' #
#'
#' data(base_did)
#' # We interact the variable 'period' with the variable 'treat'
#' est_did <- feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)
#'
#' # Now we can plot the result of the interaction with coefplot
#' coefplot(est_did)
#' # You have many more example in coefplot help
#'
#' #
#' # Instrumental variables
#' #
#'
#' # To estimate Two stage least squares,
#' # insert a formula describing the endo. vars./instr. relation after a pipe:
#'
#' base <- iris
#' names(base) <- c("y", "x1", "x2", "x3", "fe1")
#' base$x_inst1 <- 0.2 * base$x1 + 0.7 * base$x2 + rpois(150, 2)
#' base$x_inst2 <- 0.2 * base$x2 + 0.7 * base$x3 + rpois(150, 3)
#' base$x_endo1 <- 0.5 * base$y + 0.5 * base$x3 + rnorm(150, sd = 2)
#' base$x_endo2 <- 1.5 * base$y + 0.5 * base$x3 + 3 * base$x_inst1 + rnorm(150, sd = 5)
#'
#' # Using 2 controls, 1 endogenous var. and 1 instrument
#' res_iv <- feols(y ~ x1 + x2 | x_endo1 ~ x_inst1, base)
#'
#' # The second stage is the default
#' summary(res_iv)
#'
#' # To show the first stage:
#' summary(res_iv, stage = 1)
#'
#' # To show both the first and second stages:
#' summary(res_iv, stage = 1:2)
#'
#' # Adding a fixed-effect => IV formula always last!
#' res_iv_fe <- feols(y ~ x1 + x2 | fe1 | x_endo1 ~ x_inst1, base)
#'
#' # With two endogenous regressors
#' res_iv2 <- feols(y ~ x1 + x2 | x_endo1 + x_endo2 ~ x_inst1 + x_inst2, base)
#'
#' # Now there's two first stages => a fixest_multi object is returned
#' sum_res_iv2 <- summary(res_iv2, stage = 1)
#'
#' # You can navigate through it by subsetting:
#' sum_res_iv2[iv = 1]
#'
#' # The stage argument also works in etable:
#' etable(res_iv, res_iv_fe, res_iv2, order = "endo")
#'
#' etable(res_iv, res_iv_fe, res_iv2,
#'   stage = 1:2, order = c("endo", "inst"),
#'   group = list(control = "!endo|inst")
#' )
#'
#' #
#' # Multiple estimations:
#' #
#'
#' # 6 estimations
#' est_mult <- feols(c(Ozone, Solar.R) ~ Wind + Temp + csw0(Wind:Temp, Day), airquality)
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
#' est_split <- feols(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
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
#'
#' #
#' # Argument sliding
#' #
#'
#' # When the data set is set up globally, you can use the vcov argument implicitly
#'
#' base <- iris
#' names(base) <- c("y", "x1", "x2", "x3", "species")
#'
#' no_sliding <- feols(y ~ x1 + x2, base, ~species)
#'
#' # With sliding
#' setFixest_estimation(data = base)
#'
#' # ~species is implicitly deduced to be equal to 'vcov'
#' sliding <- feols(y ~ x1 + x2, ~species)
#'
#' etable(no_sliding, sliding)
#'
#' # Resetting the global options
#' setFixest_estimation(data = NULL)
#'
#'
#' #
#' # Formula expansions
#' #
#'
#' # By default, the features of the xpd function are enabled in
#' # all fixest estimations
#' # Here's a few examples
#'
#' base <- setNames(iris, c("y", "x1", "x2", "x3", "species"))
#'
#' # dot square bracket operator
#' feols(y ~ x.[1:3], base)
#'
#' # fetching variables via regular expressions: ..("regex")
#' feols(y ~ ..("1|2"), base)
#'
#' # NOTA: it also works for multiple LHS
#' mult1 <- feols(x.[1:2] ~ y + species, base)
#' mult2 <- feols(..("y|3") ~ x.[1:2] + species, base)
#' etable(mult1, mult2)
#'
#'
#' # Use .[, stuff] to include variables in functions:
#' feols(y ~ csw(x.[, 1:3]), base)
#'
#' # Same for ..(, "regex")
#' feols(y ~ csw(..(, "x")), base)
#' @export
feols <- function(fml, data, vcov, weights, offset, subset, split, fsplit, cluster, se,
                  ssc, panel.id, fixef, fixef.rm = "none", fixef.tol = 1e-6,
                  fixef.iter = 10000, collin.tol = 1e-10, nthreads = getFixest_nthreads(),
                  lean = FALSE, verbose = 0, warn = TRUE, notes = getFixest_notes(),
                  only.coef = FALSE,
                  combine.quick, demeaned = FALSE, mem.clean = FALSE, only.env = FALSE, env, ...) {
  dots <- list(...)

  # 1st: is the call coming from feglm? ----
  fromGLM <- FALSE
  skip_fixef <- FALSE
  if ("fromGLM" %in% names(dots)) {
    fromGLM <- TRUE
    # env is provided by feglm
    X <- dots$X
    y <- as.vector(dots$y)
    init <- dots$means
    correct_0w <- dots$correct_0w
    only.coef <- FALSE

    if (verbose) {
      # I can't really mutualize these three lines of code since the verbose
      # needs to be checked before using it, and here it's an internal call
      time_start <- proc.time()
      gt <- function(x, nl = TRUE) cat(sfill(x, 20), ": ", -(t0 - (t0 <<- proc.time()))[3], "s", ifelse(nl, "\n", ""), sep = "")
      t0 <- proc.time()
    }
  } else {
    time_start <- proc.time()
    gt <- function(x, nl = TRUE) cat(sfill(x, 20), ": ", -(t0 - (t0 <<- proc.time()))[3], "s", ifelse(nl, "\n", ""), sep = "")
    t0 <- proc.time()

    # we use fixest_env for appropriate controls and data handling

    if (missing(env)) {
      set_defaults("fixest_estimation")
      call_env <- new.env(parent = parent.frame())

      env <- try(fixest_env(
        fml = fml, data = data, weights = weights, offset = offset,
        subset = subset, split = split, fsplit = fsplit,
        vcov = vcov, cluster = cluster, se = se, ssc = ssc,
        panel.id = panel.id, fixef = fixef, fixef.rm = fixef.rm,
        fixef.tol = fixef.tol, fixef.iter = fixef.iter, collin.tol = collin.tol,
        nthreads = nthreads, lean = lean, verbose = verbose, warn = warn,
        notes = notes, only.coef = only.coef, combine.quick = combine.quick, demeaned = demeaned,
        mem.clean = mem.clean, origin = "feols", mc_origin = match.call(),
        call_env = call_env, ...
      ), silent = TRUE)
    } else if ((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
      stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
    }

    if ("try-error" %in% class(env)) {
      stop(format_error_msg(env, "feols"))
    }

    check_arg(only.env, "logical scalar")
    if (only.env) {
      return(env)
    }

    y <- get("lhs", env)
    X <- get("linear.mat", env)
    nthreads <- get("nthreads", env)
    init <- 0

    # demeaned variables
    if (!is.null(dots$X_demean)) {
      skip_fixef <- TRUE
      X_demean <- dots$X_demean
      y_demean <- dots$y_demean
    }

    # offset
    offset <- get("offset.value", env)
    isOffset <- length(offset) > 1
    if (isOffset) {
      y <- y - offset
    }

    # weights
    weights <- get("weights.value", env)
    isWeight <- length(weights) > 1
    correct_0w <- FALSE

    mem.clean <- get("mem.clean", env)

    demeaned <- get("demeaned", env)

    warn <- get("warn", env)

    only.coef <- get("only.coef", env)

    verbose <- get("verbose", env)
    if (verbose >= 2) gt("Setup")
  }

  isFixef <- get("isFixef", env)

  # Used to solve with the reduced model ----
  xwx <- dots$xwx
  xwy <- dots$xwy

  # Split ----

  do_split <- get("do_split", env)
  if (do_split) {
    res <- multi_split(env, feols)

    return(res)
  }

  # Multi fixef ----

  do_multi_fixef <- get("do_multi_fixef", env)
  if (do_multi_fixef) {
    res <- multi_fixef(env, feols)

    return(res)
  }

  # Multi LHS and RHS ----

  do_multi_lhs <- get("do_multi_lhs", env)
  do_multi_rhs <- get("do_multi_rhs", env)
  if (do_multi_lhs || do_multi_rhs) {
    assign("do_multi_lhs", FALSE, env)
    assign("do_multi_rhs", FALSE, env)
    do_iv <- get("do_iv", env)

    fml <- get("fml", env)
    lhs_names <- get("lhs_names", env)
    lhs <- y

    if (do_multi_lhs) {
      # We find out which LHS have the same NA patterns => saves a lot of computation

      n_lhs <- length(lhs)
      lhs_group_is_na <- list()
      lhs_group_id <- c()
      lhs_group_n_na <- c()
      for (i in 1:n_lhs) {
        is_na_current <- !is.finite(lhs[[i]])
        n_na_current <- sum(is_na_current)

        if (i == 1) {
          lhs_group_id <- 1
          lhs_group_is_na[[1]] <- is_na_current
          lhs_group_n_na[1] <- n_na_current
        } else {
          qui <- which(lhs_group_n_na == n_na_current)

          if (length(qui) > 0) {
            if (n_na_current == 0) {
              # no need to check the pattern
              lhs_group_id[i] <- lhs_group_id[qui[1]]
              next
            }

            for (j in qui) {
              if (all(is_na_current == lhs_group_is_na[[j]])) {
                lhs_group_id[i] <- lhs_group_id[j]
                next
              }
            }
          }

          # if here => new group because couldn't be matched
          id <- max(lhs_group_id) + 1
          lhs_group_id[i] <- id
          lhs_group_is_na[[id]] <- is_na_current
          lhs_group_n_na[id] <- n_na_current
        }
      }

      # we make groups
      lhs_group <- list()
      for (i in 1:max(lhs_group_id)) {
        lhs_group[[i]] <- which(lhs_group_id == i)
      }
    } else if (do_multi_lhs == FALSE) {
      lhs_group_is_na <- list(FALSE)
      lhs_group_n_na <- 0
      lhs_group <- list(1)
      lhs <- list(lhs) # I really abuse R shallow copy system...
      names(lhs) <- deparse_long(fml[[2]])
    }

    if (do_multi_rhs) {
      rhs_info_stepwise <- get("rhs_info_stepwise", env)
      multi_rhs_fml_full <- rhs_info_stepwise$fml_all_full
      multi_rhs_fml_sw <- rhs_info_stepwise$fml_all_sw
      multi_rhs_cumul <- rhs_info_stepwise$is_cumul

      linear_core <- get("linear_core", env)
      rhs <- get("rhs_sw", env)

      # Two schemes:
      #  - if cumulative: we take advantage of it => both in demeaning and in estimation
      #  - if regular stepwise => only in demeaning
      # => of course this is dependent on the pattern of NAs
      #

      n_core_left <- ifelse(length(linear_core$left) == 1, 0, ncol(linear_core$left))
      n_core_right <- ifelse(length(linear_core$right) == 1, 0, ncol(linear_core$right))

      # rnc: running number of columns
      rnc <- n_core_left
      if (rnc == 0) {
        col_start <- integer(0)
      } else {
        col_start <- 1:rnc
      }

      rhs_group_is_na <- list()
      rhs_group_id <- c()
      rhs_group_n_na <- c()
      rhs_n_vars <- c()
      rhs_col_id <- list()
      any_na_rhs <- FALSE
      for (i in seq_along(multi_rhs_fml_sw)) {
        # We evaluate the extra data and check the NA pattern

        my_fml <- multi_rhs_fml_sw[[i]]

        if (i == 1 && (multi_rhs_cumul || identical(my_fml[[3]], 1))) {
          # That case is already in the main linear.mat => no NA
          rhs_group_id <- 1
          rhs_group_is_na[[1]] <- FALSE
          rhs_group_n_na[1] <- 0
          rhs_n_vars[1] <- 0
          rhs[[1]] <- 0
          if (rnc == 0) {
            rhs_col_id[[1]] <- integer(0)
          } else {
            rhs_col_id[[1]] <- 1:rnc
          }

          next
        }

        rhs_current <- rhs[[i]]
        rhs_n_vars[i] <- ncol(rhs_current)
        info <- cpppar_which_na_inf_mat(rhs_current, nthreads)

        is_na_current <- info$is_na_inf
        if (multi_rhs_cumul && any_na_rhs) {
          # we cumulate the NAs
          is_na_current <- is_na_current | rhs_group_is_na[[rhs_group_id[i - 1]]]
          info$any_na_inf <- any(is_na_current)
        }

        n_na_current <- 0
        if (info$any_na_inf) {
          any_na_rhs <- TRUE
          n_na_current <- sum(is_na_current)
        } else {
          # NULL would lead to problems down the road
          is_na_current <- FALSE
        }

        if (i == 1) {
          rhs_group_id <- 1
          rhs_group_is_na[[1]] <- is_na_current
          rhs_group_n_na[1] <- n_na_current
        } else {
          qui <- which(rhs_group_n_na == n_na_current)

          if (length(qui) > 0) {
            if (n_na_current == 0) {
              # no need to check the pattern
              rhs_group_id[i] <- rhs_group_id[qui[1]]
              next
            }

            go_next <- FALSE
            for (j in qui) {
              if (all(is_na_current == rhs_group_is_na[[j]])) {
                rhs_group_id[i] <- j
                go_next <- TRUE
                break
              }
            }

            if (go_next) next
          }

          # if here => new group because couldn't be matched
          id <- max(rhs_group_id) + 1
          rhs_group_id[i] <- id
          rhs_group_is_na[[id]] <- is_na_current
          rhs_group_n_na[id] <- n_na_current
        }
      }


      # we make groups
      rhs_group <- list()
      for (i in 1:max(rhs_group_id)) {
        rhs_group[[i]] <- which(rhs_group_id == i)
      }

      # Finding the right column IDs to select

      rhs_group_n_vars <- rep(0, length(rhs_group)) # To get the total nber of cols per group

      for (i in seq_along(multi_rhs_fml_sw)) {
        if (multi_rhs_cumul) {
          rnc <- rnc + rhs_n_vars[i]
          if (rnc == 0) {
            rhs_col_id[[i]] <- integer(0)
          } else {
            rhs_col_id[[i]] <- 1:rnc
          }
        } else {
          id <- rhs_group_id[i]
          rhs_col_id[[i]] <- c(col_start, seq(rnc + rhs_group_n_vars[id] + 1, length.out = rhs_n_vars[i]))
          rhs_group_n_vars[id] <- rhs_group_n_vars[id] + rhs_n_vars[i]
        }
      }

      if (n_core_right > 0) {
        # We adjust
        if (multi_rhs_cumul) {
          for (i in seq_along(multi_rhs_fml_sw)) {
            id <- rhs_group_id[i]
            gmax <- max(rhs_group[[id]])
            rhs_col_id[[i]] <- c(rhs_col_id[[i]], n_core_left + sum(rhs_n_vars[1:gmax]) + 1:n_core_right)
          }
        } else {
          for (i in seq_along(multi_rhs_fml_sw)) {
            id <- rhs_group_id[i]
            rhs_col_id[[i]] <- c(rhs_col_id[[i]], n_core_left + rhs_group_n_vars[id] + 1:n_core_right)
          }
        }
      }
    } else if (do_multi_rhs == FALSE) {
      multi_rhs_fml_full <- list(.xpd(rhs = fml[[3]]))
      multi_rhs_cumul <- FALSE
      rhs_group_is_na <- list(FALSE)
      rhs_group_n_na <- 0
      rhs_n_vars <- 0
      rhs_group <- list(1)
      rhs <- list(0)
      rhs_col_id <- list(1:NCOL(X))
      linear_core <- list(left = X, right = 1)
    }

    isLinear_right <- length(linear_core$right) > 1

    isLinear <- length(linear_core$left) > 1 || isLinear_right

    n_lhs <- length(lhs)
    n_rhs <- length(rhs)
    res <- vector("list", n_lhs * n_rhs)

    rhs_names <- sapply(multi_rhs_fml_full, function(x) as.character(x)[[2]])

    for (i in seq_along(lhs_group)) {
      for (j in seq_along(rhs_group)) {
        # NA removal
        no_na <- FALSE
        if (lhs_group_n_na[i] > 0) {
          if (rhs_group_n_na[j] > 0) {
            is_na_current <- lhs_group_is_na[[i]] | rhs_group_is_na[[j]]
          } else {
            is_na_current <- lhs_group_is_na[[i]]
          }
        } else if (rhs_group_n_na[j] > 0) {
          is_na_current <- rhs_group_is_na[[j]]
        } else {
          no_na <- TRUE
        }

        # Here it depends on whether there are FEs or not, whether it's cumul or not
        my_lhs <- lhs[lhs_group[[i]]]
        if (isLinear) {
          my_rhs <- linear_core[1]
          if (multi_rhs_cumul) {
            gmax <- max(rhs_group[[j]])
            my_rhs[1 + (1:gmax)] <- rhs[1:gmax]
          } else {
            for (u in rhs_group[[j]]) {
              if (length(rhs[[u]]) > 1) {
                my_rhs[[length(my_rhs) + 1]] <- rhs[[u]]
              }
            }
          }

          if (isLinear_right) {
            my_rhs[[length(my_rhs) + 1]] <- linear_core$right
          }
        } else {
          rhs_len <- lengths(rhs)
          if (multi_rhs_cumul) {
            gmax <- max(rhs_group[[j]])
            my_rhs <- rhs[rhs_len > 1 & seq_along(rhs) <= gmax]
          } else {
            my_rhs <- rhs[rhs_len > 1 & seq_along(rhs) %in% rhs_group[[j]]]
          }

          if (isLinear_right) {
            my_rhs[[length(my_rhs) + 1]] <- linear_core$right
          }
        }

        len_all <- lengths(my_rhs)
        if (any(len_all == 1)) {
          my_rhs <- my_rhs[len_all > 1]
        }

        if (!no_na) {
          # NA removal
          for (u in seq_along(my_lhs)) {
            my_lhs[[u]] <- my_lhs[[u]][!is_na_current]
          }

          for (u in seq_along(my_rhs)) {
            if (length(my_rhs[[u]]) > 1) my_rhs[[u]] <- my_rhs[[u]][!is_na_current, , drop = FALSE]
          }

          my_env <- reshape_env(env, obs2keep = which(!is_na_current), assign_lhs = FALSE, assign_rhs = FALSE)
        } else {
          my_env <- reshape_env(env)
        }

        weights <- get("weights.value", my_env)

        all_varnames <- NULL
        isLinear_current <- TRUE
        if (length(my_rhs) == 0) {
          X_all <- 0
          isLinear_current <- FALSE
        } else {
          # We try to avoid repeating variables
          # => can happen in stepwise estimations (not in csw)

          all_varnames <- unlist(sapply(my_rhs, colnames))
          all_varnames_unik <- unique(all_varnames)
          all_varnames_done <- rep(FALSE, length(all_varnames_unik))
          all_varnames_done <- all_varnames_unik %in% colnames(my_rhs[[1]])

          dict_vars <- 1:length(all_varnames_unik)
          names(dict_vars) <- all_varnames_unik

          n_rhs_current <- length(my_rhs)
          args_cbind <- vector("list", n_rhs_current)
          args_cbind[[1]] <- my_rhs[[1]]

          id_rhs <- 2
          while (!all(all_varnames_done) && id_rhs <= n_rhs_current) {
            rhs_current <- my_rhs[[id_rhs]]
            qui <- !colnames(rhs_current) %in% all_varnames_unik[all_varnames_done]
            if (any(qui)) {
              args_cbind[[id_rhs]] <- rhs_current[, qui, drop = FALSE]
              all_varnames_done <- all_varnames_done | all_varnames_unik %in% colnames(rhs_current)
            }

            id_rhs <- id_rhs + 1
          }

          X_all <- do.call("cbind", args_cbind)
        }

        if (do_iv) {
          # We need to GET them => they have been modified in my_env
          iv_lhs <- get("iv_lhs", my_env)
          iv.mat <- get("iv.mat", my_env)
          n_inst <- ncol(iv.mat)
        }

        if (isFixef) {
          # We batch demean

          n_vars_X <- ifelse(is.null(ncol(X_all)), 0, ncol(X_all))

          # fixef information
          fixef_sizes <- get("fixef_sizes", my_env)
          fixef_table_vector <- get("fixef_table_vector", my_env)
          fixef_id_list <- get("fixef_id_list", my_env)

          slope_flag <- get("slope_flag", my_env)
          slope_vars <- get("slope_variables", my_env)

          if (mem.clean) gc()

          vars_demean <- cpp_demean(my_lhs, X_all, weights,
            iterMax = fixef.iter,
            diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
            fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
            slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
            r_init = init, nthreads = nthreads, save_fixef = FALSE
          )

          X_demean <- vars_demean$X_demean
          y_demean <- vars_demean$y_demean

          if (do_iv) {
            iv_vars_demean <- cpp_demean(iv_lhs, iv.mat, weights,
              iterMax = fixef.iter,
              diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
              fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
              slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
              r_init = init, nthreads = nthreads, save_fixef = FALSE
            )

            iv.mat_demean <- iv_vars_demean$X_demean
            iv_lhs_demean <- iv_vars_demean$y_demean
          }
        }

        # We precompute the solution
        if (do_iv) {
          if (isFixef) {
            iv_products <- cpp_iv_products(
              X = X_demean, y = y_demean,
              Z = iv.mat_demean, u = iv_lhs_demean,
              w = weights, nthreads = nthreads
            )
          } else {
            iv_products <- cpp_iv_products(
              X = X_all, y = my_lhs, Z = iv.mat,
              u = iv_lhs, w = weights, nthreads = nthreads
            )
          }
        } else {
          if (isFixef) {
            my_products <- cpp_sparse_products(X_demean, weights, y_demean, nthreads = nthreads)
          } else {
            my_products <- cpp_sparse_products(X_all, weights, my_lhs, nthreads = nthreads)
          }

          xwx <- my_products$XtX
          xwy <- my_products$Xty
        }

        for (ii in seq_along(my_lhs)) {
          i_lhs <- lhs_group[[i]][ii]

          for (jj in rhs_group[[j]]) {
            # linking the unique variables to the variables
            qui_X <- rhs_col_id[[jj]]
            if (!is.null(all_varnames)) {
              qui_X <- dict_vars[all_varnames[qui_X]]
            }

            if (isLinear_current) {
              my_X <- X_all[, qui_X, drop = FALSE]
            } else {
              my_X <- 0
            }

            my_fml <- .xpd(lhs = lhs_names[i_lhs], rhs = multi_rhs_fml_full[[jj]])
            current_env <- reshape_env(my_env, lhs = my_lhs[[ii]], rhs = my_X, fml_linear = my_fml)

            if (do_iv) {
              if (isLinear_current) {
                qui_iv <- c(1:n_inst, n_inst + qui_X)

                XtX <- iv_products$XtX[qui_X, qui_X, drop = FALSE]
                Xty <- iv_products$Xty[[ii]][qui_X]
              } else {
                qui_iv <- 1:n_inst

                XtX <- matrix(0, 1, 1)
                Xty <- matrix(0, 1, 1)
              }

              my_iv_products <- list(
                XtX = XtX,
                Xty = Xty,
                ZXtZX = iv_products$ZXtZX[qui_iv, qui_iv, drop = FALSE],
                ZXtu = lapply(iv_products$ZXtu, function(x) x[qui_iv])
              )

              if (isFixef) {
                my_res <- feols(
                  env = current_env, iv_products = my_iv_products,
                  X_demean = X_demean[, qui_X, drop = FALSE],
                  y_demean = y_demean[[ii]],
                  iv.mat_demean = iv.mat_demean, iv_lhs_demean = iv_lhs_demean
                )
              } else {
                my_res <- feols(env = current_env, iv_products = my_iv_products)
              }
            } else {
              if (isFixef) {
                my_res <- feols(
                  env = current_env, xwx = xwx[qui_X, qui_X, drop = FALSE], xwy = xwy[[ii]][qui_X],
                  X_demean = X_demean[, qui_X, drop = FALSE],
                  y_demean = y_demean[[ii]]
                )
              } else {
                my_res <- feols(env = current_env, xwx = xwx[qui_X, qui_X, drop = FALSE], xwy = xwy[[ii]][qui_X])
              }
            }


            res[[index_2D_to_1D(i_lhs, jj, n_rhs)]] <- my_res
          }
        }
      }
    }

    # Meta information for fixest_multi

    index <- list(lhs = n_lhs, rhs = n_rhs)
    all_names <- list(lhs = lhs_names, rhs = rhs_names)

    # result
    res_multi <- setup_multi(index, all_names, res)

    return(res_multi)
  }


  #
  # IV ####
  #

  do_iv <- get("do_iv", env)
  if (do_iv) {
    assign("do_iv", FALSE, env)
    assign("verbose", 0, env)

    # Loaded already
    # y: lhs
    # X: linear.mat

    iv_lhs <- get("iv_lhs", env)
    iv_lhs_names <- get("iv_lhs_names", env)
    iv.mat <- get("iv.mat", env) # we enforce (before) at least one variable in iv.mat
    K <- ncol(iv.mat)
    n_endo <- length(iv_lhs)
    lean <- get("lean", env)

    # Simple check that the function is not misused
    pblm <- intersect(iv_lhs_names, colnames(X))
    if (length(pblm) > 0) {
      any_exo <- length(setdiff(colnames(X), iv_lhs_names)) > 0
      msg <- if (any_exo) "" else " If there is no exogenous variable, just use '1' in the first part of the formula."
      stop("Endogenous variables should not be used as exogenous regressors. The variable", enumerate_items(pblm, "s.quote.were"), " found in the first part of the multipart formula: ", ifsingle(pblm, "it", "they"), " should not be there.", msg)
    }

    if (isFixef) {
      # we batch demean first

      n_vars_X <- ifelse(is.null(ncol(X)), 0, ncol(X))

      if (mem.clean) gc()

      if (!is.null(dots$iv_products)) {
        # means this is a call from multiple LHS/RHS

        X_demean <- dots$X_demean
        y_demean <- dots$y_demean

        iv.mat_demean <- dots$iv.mat_demean
        iv_lhs_demean <- dots$iv_lhs_demean

        iv_products <- dots$iv_products
      } else {
        # fixef information
        fixef_sizes <- get("fixef_sizes", env)
        fixef_table_vector <- get("fixef_table_vector", env)
        fixef_id_list <- get("fixef_id_list", env)

        slope_flag <- get("slope_flag", env)
        slope_vars <- get("slope_variables", env)

        vars_demean <- cpp_demean(y, X, weights,
          iterMax = fixef.iter,
          diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
          fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
          slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
          r_init = init, nthreads = nthreads, save_fixef = FALSE
        )

        iv_vars_demean <- cpp_demean(iv_lhs, iv.mat, weights,
          iterMax = fixef.iter,
          diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
          fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
          slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
          r_init = init, nthreads = nthreads, save_fixef = FALSE
        )

        X_demean <- vars_demean$X_demean
        y_demean <- vars_demean$y_demean

        iv.mat_demean <- iv_vars_demean$X_demean
        iv_lhs_demean <- iv_vars_demean$y_demean

        # We precompute the solution
        iv_products <- cpp_iv_products(
          X = X_demean, y = y_demean, Z = iv.mat_demean,
          u = iv_lhs_demean, w = weights, nthreads = nthreads
        )
      }

      if (n_vars_X == 0) {
        ZX_demean <- iv.mat_demean
        ZX <- iv.mat
      } else {
        ZX_demean <- cbind(iv.mat_demean, X_demean)
        ZX <- cbind(iv.mat, X)
      }

      # First stage(s)

      ZXtZX <- iv_products$ZXtZX
      ZXtu <- iv_products$ZXtu

      res_first_stage <- list()

      for (i in 1:n_endo) {
        current_env <- reshape_env(env, lhs = iv_lhs[[i]], rhs = ZX, fml_iv_endo = iv_lhs_names[i])

        my_res <- feols(
          env = current_env, xwx = ZXtZX, xwy = ZXtu[[i]],
          X_demean = ZX_demean, y_demean = iv_lhs_demean[[i]],
          add_fitted_demean = TRUE, iv_call = TRUE, notes = FALSE
        )

        # For the F-stats
        if (n_vars_X == 0) {
          my_res$ssr_no_inst <- cpp_ssq(iv_lhs_demean[[i]], weights)
        } else {
          fit_no_inst <- ols_fit(iv_lhs_demean[[i]], X_demean,
            w = weights, correct_0w = FALSE,
            collin.tol = collin.tol, nthreads = nthreads,
            xwx = iv_products$XtX, xwy = ZXtu[[i]][-(1:K)]
          )
          my_res$ssr_no_inst <- cpp_ssq(fit_no_inst$residuals, weights)
        }

        my_res$iv_stage <- 1
        my_res$iv_inst_names_xpd <- colnames(iv.mat)

        res_first_stage[[iv_lhs_names[i]]] <- my_res
      }

      if (verbose >= 2) gt("1st stage(s)")

      # Second stage

      if (n_endo == 1) {
        res_FS <- res_first_stage[[1]]
        U <- as.matrix(res_FS$fitted.values)
        U_demean <- as.matrix(res_FS$fitted.values_demean)
      } else {
        U_list <- list()
        U_dm_list <- list()
        for (i in 1:n_endo) {
          res_FS <- res_first_stage[[i]]
          U_list[[i]] <- res_FS$fitted.values
          U_dm_list[[i]] <- res_FS$fitted.values_demean
        }

        U <- do.call("cbind", U_list)
        U_demean <- do.call("cbind", U_dm_list)
      }

      colnames(U) <- colnames(U_demean) <- paste0("fit_", iv_lhs_names)

      if (n_vars_X == 0) {
        UX <- as.matrix(U)
        UX_demean <- as.matrix(U_demean)
      } else {
        UX <- cbind(U, X)
        UX_demean <- cbind(U_demean, X_demean)
      }

      XtX <- iv_products$XtX
      Xty <- iv_products$Xty
      iv_prod_second <- cpp_iv_product_completion(
        XtX = XtX, Xty = Xty, X = X_demean,
        y = y_demean, U = U_demean, w = weights, nthreads = nthreads
      )

      UXtUX <- iv_prod_second$UXtUX
      UXty <- iv_prod_second$UXty

      resid_s1 <- lapply(res_first_stage, function(x) x$residuals)

      current_env <- reshape_env(env, rhs = UX)
      res_second_stage <- feols(
        env = current_env, xwx = UXtUX, xwy = UXty,
        X_demean = UX_demean, y_demean = y_demean,
        resid_1st_stage = resid_s1, iv_call = TRUE, notes = FALSE
      )

      # For the F-stats
      if (n_vars_X == 0) {
        res_second_stage$ssr_no_endo <- cpp_ssq(y_demean, weights)
      } else {
        fit_no_endo <- ols_fit(y_demean, X_demean,
          w = weights, correct_0w = FALSE,
          collin.tol = collin.tol, nthreads = nthreads,
          xwx = XtX, xwy = Xty
        )
        res_second_stage$ssr_no_endo <- cpp_ssq(fit_no_endo$residuals, weights)
      }
    } else {
      # fixef == FALSE

      # We precompute the solution
      if (!is.null(dots$iv_products)) {
        # means this is a call from multiple LHS/RHS
        iv_products <- dots$iv_products
      } else {
        iv_products <- cpp_iv_products(
          X = X, y = y, Z = iv.mat,
          u = iv_lhs, w = weights, nthreads = nthreads
        )
      }

      if (verbose >= 2) gt("IV products")

      ZX <- cbind(iv.mat, X)

      # First stage(s)

      ZXtZX <- iv_products$ZXtZX
      ZXtu <- iv_products$ZXtu

      # Let's put the intercept first => I know it's not really elegant, but that's life
      is_int <- "(Intercept)" %in% colnames(X)
      if (is_int) {
        nz <- ncol(iv.mat)
        nzx <- ncol(ZX)
        qui <- c(nz + 1, (1:nzx)[-(nz + 1)])

        ZX <- ZX[, qui, drop = FALSE]
        ZXtZX <- ZXtZX[qui, qui, drop = FALSE]
        for (i in seq_along(ZXtu)) {
          ZXtu[[i]] <- ZXtu[[i]][qui]
        }
      }

      res_first_stage <- list()

      for (i in 1:n_endo) {
        current_env <- reshape_env(env, lhs = iv_lhs[[i]], rhs = ZX, fml_iv_endo = iv_lhs_names[i])
        my_res <- feols(
          env = current_env, xwx = ZXtZX, xwy = ZXtu[[i]],
          iv_call = TRUE, notes = FALSE
        )

        # For the F-stats
        fit_no_inst <- ols_fit(iv_lhs[[i]], X,
          w = weights, correct_0w = FALSE,
          collin.tol = collin.tol, nthreads = nthreads,
          xwx = ZXtZX[-(1:K + is_int), -(1:K + is_int), drop = FALSE],
          xwy = ZXtu[[i]][-(1:K + is_int)]
        )
        my_res$ssr_no_inst <- cpp_ssq(fit_no_inst$residuals, weights)

        my_res$iv_stage <- 1
        my_res$iv_inst_names_xpd <- colnames(iv.mat)

        res_first_stage[[iv_lhs_names[i]]] <- my_res
      }

      if (verbose >= 2) gt("1st stage(s)")

      # Second stage

      if (n_endo == 1) {
        res_FS <- res_first_stage[[1]]
        U <- as.matrix(res_FS$fitted.values)
      } else {
        U_list <- list()
        U_dm_list <- list()
        for (i in 1:n_endo) {
          res_FS <- res_first_stage[[i]]
          U_list[[i]] <- res_FS$fitted.values
        }

        U <- do.call("cbind", U_list)
      }

      colnames(U) <- paste0("fit_", iv_lhs_names)

      UX <- cbind(U, X)

      XtX <- iv_products$XtX
      Xty <- iv_products$Xty
      iv_prod_second <- cpp_iv_product_completion(
        XtX = XtX, Xty = Xty, X = X,
        y = y, U = U, w = weights, nthreads = nthreads
      )

      UXtUX <- iv_prod_second$UXtUX
      UXty <- iv_prod_second$UXty

      if (is_int) {
        nu <- ncol(U)
        nux <- ncol(UX)
        qui <- c(nu + 1, (1:nux)[-(nu + 1)])

        UX <- UX[, qui, drop = FALSE]
        UXtUX <- UXtUX[qui, qui, drop = FALSE]
        UXty <- UXty[qui]
      }

      resid_s1 <- lapply(res_first_stage, function(x) x$residuals)

      current_env <- reshape_env(env, rhs = UX)
      res_second_stage <- feols(
        env = current_env, xwx = UXtUX, xwy = UXty,
        resid_1st_stage = resid_s1, iv_call = TRUE, notes = FALSE
      )

      # For the F-stats
      fit_no_endo <- ols_fit(y, X,
        w = weights, correct_0w = FALSE,
        collin.tol = collin.tol, nthreads = nthreads,
        xwx = XtX, xwy = Xty
      )
      res_second_stage$ssr_no_endo <- cpp_ssq(fit_no_endo$residuals, weights)
    }

    if (verbose >= 2) gt("2nd stage")

    #
    # Wu-Hausman endogeneity test
    #

    # Current limitation => only standard vcov => later add argument (which would yield the full est.)?
    # The problem of the full est. is that it takes memory very likely needlessly

    if (isFixef) {
      ENDO_demean <- do.call(cbind, iv_lhs_demean)
      iv_prod_wh <- cpp_iv_product_completion(
        XtX = UXtUX, Xty = UXty,
        X = UX_demean, y = y_demean, U = ENDO_demean,
        w = weights, nthreads = nthreads
      )

      RHS_wh <- cbind(ENDO_demean, UX_demean)
      fit_wh <- ols_fit(y_demean, RHS_wh,
        w = weights, correct_0w = FALSE, collin.tol = collin.tol,
        nthreads = nthreads, xwx = iv_prod_wh$UXtUX, xwy = iv_prod_wh$UXty
      )
    } else {
      ENDO <- do.call(cbind, iv_lhs)
      iv_prod_wh <- cpp_iv_product_completion(
        XtX = UXtUX, Xty = UXty,
        X = UX, y = y, U = ENDO,
        w = weights, nthreads = nthreads
      )

      RHS_wh <- cbind(ENDO, UX)
      fit_wh <- ols_fit(y, RHS_wh,
        w = weights, correct_0w = FALSE, collin.tol = collin.tol,
        nthreads = nthreads, xwx = iv_prod_wh$UXtUX, xwy = iv_prod_wh$UXty
      )
    }

    df1 <- n_endo
    df2 <- length(y) - (res_second_stage$nparams + df1)
    if (any(fit_wh$is_excluded)) {
      stat <- p <- NA
    } else {
      qui <- df1 + 1:df1 + ("(Intercept)" %in% names(res_second_stage$coefficients))
      my_coef <- fit_wh$coefficients[qui]
      vcov_wh <- fit_wh$xwx_inv[qui, qui] * cpp_ssq(fit_wh$residuals, weights) / df2
      stat <- drop(my_coef %*% solve(vcov_wh) %*% my_coef) / df1
      p <- pf(stat, df1, df2, lower.tail = FALSE)
    }

    res_second_stage$iv_wh <- list(stat = stat, p = p, df1 = df1, df2 = df2)

    #
    # Sargan
    #

    if (n_endo < ncol(iv.mat)) {
      df <- ncol(iv.mat) - n_endo
      resid_2nd <- res_second_stage$residuals

      if (isFixef) {
        xwy <- cpppar_xwy(ZX_demean, resid_2nd, weights, nthreads)
        fit_sargan <- ols_fit(resid_2nd, ZX_demean,
          w = weights, correct_0w = FALSE, collin.tol = collin.tol,
          nthreads = nthreads, xwx = ZXtZX, xwy = xwy
        )
      } else {
        xwy <- cpppar_xwy(ZX, resid_2nd, weights, nthreads)
        fit_sargan <- ols_fit(resid_2nd, ZX,
          w = weights, correct_0w = FALSE, collin.tol = collin.tol,
          nthreads = nthreads, xwx = ZXtZX, xwy = xwy
        )
      }

      r <- fit_sargan$residuals
      stat <- length(r) * (1 - cpp_ssq(r, weights) / cpp_ssr_null(resid_2nd))
      p <- pchisq(stat, df, lower.tail = FALSE)
      res_second_stage$iv_sargan <- list(stat = stat, p = p, df = df)
    }

    # extra information
    res_second_stage$iv_inst_names_xpd <- res_first_stage[[1]]$iv_inst_names_xpd
    res_second_stage$iv_endo_names_fit <- paste0("fit_", res_second_stage$iv_endo_names)

    # Collinearity message

    collin.vars <- c(res_second_stage$collin.var, res_first_stage[[1]]$collin.var)
    res_second_stage$collin.var <- unique(collin.vars)
    if (notes && length(collin.vars) > 0) {
      coll.endo <- intersect(collin.vars, res_second_stage$iv_endo_names_fit)
      coll.inst <- intersect(collin.vars, res_second_stage$iv_inst_names_xpd)
      coll.exo <- setdiff(collin.vars, c(coll.endo, coll.inst))

      n_c <- length(collin.vars)
      n_c_endo <- length(coll.endo)
      n_c_inst <- length(coll.inst)
      n_c_exo <- length(coll.exo)

      msg_endo <- msg_exo <- msg_inst <- NULL
      if (n_c_endo > 0) {
        msg_endo <- paste0("The endogenous regressor", plural(n_c_endo), " ", enumerate_items(coll.endo, "quote", nmax = 3))
      } else if (n_c_inst > 0) {
        msg_inst <- paste0("the instrument", plural(n_c_inst), " ", enumerate_items(coll.inst, "quote", nmax = 3))
      } else if (n_c_exo > 0) {
        msg_exo <- paste0("the exogenous variable", plural(n_c_exo), " ", enumerate_items(coll.exo, "quote", nmax = 3))
      }

      msg <- enumerate_items(c(msg_endo, msg_inst, msg_exo))
      msg <- gsub("^t", "T", msg)
      message(msg, " ", plural(n_c, "has"), " been removed because of collinearity (see $collin.var).")
    }

    # if lean = TRUE: we clean the IV residuals (which were needed so far)
    if (lean) {
      for (i in 1:n_endo) {
        res_first_stage[[i]]$residuals <- NULL
        res_first_stage[[i]]$fitted.values <- NULL
        res_first_stage[[i]]$fitted.values_demean <- NULL
      }

      res_second_stage$residuals <- NULL
      res_second_stage$fitted.values <- NULL
      res_second_stage$fitted.values_demean <- NULL
    }

    res_second_stage$iv_first_stage <- res_first_stage

    # meta info
    res_second_stage$iv_stage <- 2

    return(res_second_stage)
  }


  #
  # Regular estimation ####
  #

  onlyFixef <- length(X) == 1 || ncol(X) == 0

  if (fromGLM) {
    res <- list(coefficients = NA)
  } else {
    res <- get("res", env)
  }

  if (skip_fixef) {
    # Variables were already demeaned
  } else if (!isFixef) {
    # No Fixed-effects
    y_demean <- y
    X_demean <- X
    res$means <- 0
  } else {
    time_demean <- proc.time()

    # Number of nthreads
    n_vars_X <- ifelse(is.null(ncol(X)), 0, ncol(X))

    # fixef information
    fixef_sizes <- get("fixef_sizes", env)
    fixef_table_vector <- get("fixef_table_vector", env)
    fixef_id_list <- get("fixef_id_list", env)

    slope_flag <- get("slope_flag", env)
    slope_vars <- get("slope_variables", env)

    if (mem.clean) {
      # we can't really rm many variables... but gc can be enough
      # cpp_demean is the most mem intensive bit
      gc()
    }

    vars_demean <- cpp_demean(y, X, weights,
      iterMax = fixef.iter,
      diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
      fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
      slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
      r_init = init, nthreads = nthreads, save_fixef = FALSE
    )

    y_demean <- vars_demean$y_demean
    if (onlyFixef) {
      X_demean <- matrix(1, nrow = length(y_demean))
    } else {
      X_demean <- vars_demean$X_demean
    }

    res$iterations <- vars_demean$iterations
    if (fromGLM) {
      res$means <- vars_demean$means
    }

    if (mem.clean) {
      rm(vars_demean)
    }

    if (any(abs(slope_flag) > 0) && any(res$iterations > 300)) {
      # Maybe we have a convergence problem
      # This is poorly coded, but it's a temporary fix
      opt_fe <- check_conv(y_demean, X_demean, fixef_id_list, slope_flag, slope_vars, weights)

      # This is a bit too rough a check but it should catch the most problematic cases
      if (anyNA(opt_fe) || any(opt_fe > 1e-4)) {
        msg <- "There seems to be a convergence problem due to the presence of variables with varying slopes. The precision of the estimates may not be great."
        if (any(slope_flag < 0)) {
          sugg <- "This convergence problem mostly arises when there are varying slopes without their associated fixed-effect, as is the case in your estimation. Why not try to include the fixed-effect (i.e. use '[' instead of '[[')?"
        } else {
          sugg <- "As a workaround, and if there are not too many slopes, you can use the variables with varying slopes as regular variables using the function i (see ?i)."
        }

        sugg_tol <- "Or use a lower 'fixef.tol' (sometimes it works)?"

        msg <- paste(msg, sugg, sugg_tol)

        res$convStatus <- FALSE

        res$message <- paste0("tol: ", signif_plus(fixef.tol), ", iter: ", max(res$iterations))

        if (fromGLM) {
          res$warn_varying_slope <- msg
        } else if (warn) {
          warning(msg)
        }
      }
    } else if (any(res$iterations >= fixef.iter)) {
      msg <- paste0("Demeaning algorithm: Absence of convergence after reaching the maximum number of iterations (", fixef.iter, ").")

      res$convStatus <- FALSE
      res$message <- paste0("Maximum of ", fixef.iter, " iterations reached.")

      if (fromGLM) {
        res$warn_varying_slope <- msg
      } else {
        warning(msg)
      }
    }

    if (verbose >= 1) {
      if (length(fixef_sizes) > 1) {
        gt("Demeaning", FALSE)
        cat(" (iter: ", paste0(c(tail(res$iterations, 1), res$iterations[-length(res$iterations)]), collapse = ", "), ")\n", sep = "")
      } else {
        gt("Demeaning")
      }
    }
  }

  #
  # Estimation
  #

  if (mem.clean) {
    gc()
  }

  if (!onlyFixef) {
    est <- ols_fit(y_demean, X_demean, weights, correct_0w, collin.tol, nthreads,
      xwx, xwy,
      only.coef = only.coef
    )

    if (only.coef) {
      names(est) <- colnames(X)
      return(est)
    }

    if (mem.clean) {
      gc()
    }

    # Corner case: not any relevant variable
    if (!is.null(est$all_removed)) {
      all_vars <- colnames(X)

      IN_MULTI <- get("IN_MULTI", env)

      if (isFixef) {
        msg <- paste0(ifsingle(all_vars, "The only variable ", "All variables"), enumerate_items(all_vars, "quote.is", nmax = 3), " collinear with the fixed effects. In such circumstances, the estimation is void.")
      } else {
        msg <- paste0(ifsingle(all_vars, "The only variable ", "All variables"), enumerate_items(all_vars, "quote.is", nmax = 3), " virtually constant and equal to 0. In such circumstances, the estimation is void.")
      }

      if (IN_MULTI || !warn) {
        if (warn) warning(msg)

        return(fixest_NA_results(env))
      } else {
        stop_up(msg, up = fromGLM)
      }
    }

    # Formatting the result
    coef <- est$coefficients
    names(coef) <- colnames(X)[!est$is_excluded]
    res$coefficients <- coef
    # Additional stuff
    res$residuals <- est$residuals
    res$multicol <- est$multicol
    res$collin.min_norm <- est$collin.min_norm
    if (fromGLM) res$is_excluded <- est$is_excluded

    if (demeaned) {
      res$y_demeaned <- y_demean
      res$X_demeaned <- X_demean
      colnames(res$X_demeaned) <- colnames(X)
    }
  } else {
    res$residuals <- y_demean
    res$coefficients <- coef <- NULL
    res$onlyFixef <- TRUE
    res$multicol <- FALSE

    if (demeaned) {
      res$y_demeaned <- y_demean
    }
  }

  time_post <- proc.time()
  if (verbose >= 1) {
    gt("Estimation")
  }

  if (mem.clean) {
    gc()
  }

  if (fromGLM) {
    res$fitted.values <- y - res$residuals
    if (!onlyFixef) {
      res$X_demean <- X_demean
    }

    return(res)
  }

  #
  # Post processing
  #

  # Collinearity message
  collin.adj <- 0
  if (res$multicol) {
    var_collinear <- colnames(X)[est$is_excluded]
    if (notes) {
      message(ifsingle(var_collinear, "The variable ", "Variables "), enumerate_items(var_collinear, "quote.has", nmax = 3), " been removed because of collinearity (see $collin.var).")
    }

    res$collin.var <- var_collinear

    # full set of coefficients with NAs
    collin.coef <- setNames(rep(NA, ncol(X)), colnames(X))
    collin.coef[!est$is_excluded] <- res$coefficients
    res$collin.coef <- collin.coef

    if (isFixef) {
      X <- X[, !est$is_excluded, drop = FALSE]
    }
    X_demean <- X_demean[, !est$is_excluded, drop = FALSE]

    collin.adj <- sum(est$is_excluded)
  }

  n <- length(y)
  res$nparams <- res$nparams - collin.adj
  df_k <- res$nparams
  res$nobs <- n

  if (isWeight) res$weights <- weights

  #
  # IV correction
  #

  resid_origin <- NULL
  if (!is.null(dots$resid_1st_stage)) {
    # We correct the residual
    is_int <- "(Intercept)" %in% names(res$coefficients)
    resid_origin <- res$residuals
    resid_new <- cpp_iv_resid(res$residuals, res$coefficients, dots$resid_1st_stage, is_int, nthreads)
    res$iv_residuals <- res$residuals
    res$residuals <- resid_new
  }

  #
  # Hessian, score, etc
  #

  if (onlyFixef) {
    res$fitted.values <- res$sumFE <- y - res$residuals
  } else {
    if (mem.clean) {
      gc()
    }

    # X_beta / fitted / sumFE
    if (isFixef) {
      x_beta <- cpppar_xbeta(X, coef, nthreads)

      if (!is.null(resid_origin)) {
        res$sumFE <- y - x_beta - resid_origin
      } else {
        res$sumFE <- y - x_beta - res$residuals
      }

      res$fitted.values <- x_beta + res$sumFE
      if (isTRUE(dots$add_fitted_demean)) {
        res$fitted.values_demean <- est$fitted.values
      }
    } else {
      res$fitted.values <- est$fitted.values
    }

    if (isOffset) {
      res$fitted.values <- res$fitted.values + offset
    }

    #
    # score + hessian + vcov
    if (isWeight) {
      res$scores <- (res$residuals * weights) * X_demean
    } else {
      res$scores <- res$residuals * X_demean
    }

    res$hessian <- est$xwx

    if (mem.clean) {
      gc()
    }

    res$sigma2 <- cpp_ssq(res$residuals, weights) / (length(y) - df_k)

    res$cov.iid <- est$xwx_inv * res$sigma2

    rownames(res$cov.iid) <- colnames(res$cov.iid) <- names(coef)

    # for compatibility with conleyreg
    res$cov.unscaled <- res$cov.iid

    # se
    se <- diag(res$cov.iid)
    se[se < 0] <- NA
    se <- sqrt(se)

    # coeftable
    zvalue <- coef / se
    pvalue <- 2 * pt(-abs(zvalue), max(n - df_k, 1))

    coeftable <- data.frame("Estimate" = coef, "Std. Error" = se, "t value" = zvalue, "Pr(>|t|)" = pvalue)
    names(coeftable) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    row.names(coeftable) <- names(coef)

    attr(se, "type") <- attr(coeftable, "type") <- "IID"
    res$coeftable <- coeftable
    res$se <- se
  }

  # fit stats
  if (!cpp_isConstant(res$fitted.values)) {
    res$sq.cor <- stats::cor(y, res$fitted.values)**2
  } else {
    res$sq.cor <- NA
  }

  if (mem.clean) {
    gc()
  }

  res$ssr_null <- cpp_ssr_null(y, weights)
  res$ssr <- cpp_ssq(res$residuals, weights)
  sigma_null <- sqrt(res$ssr_null / ifelse(isWeight, sum(weights), n))
  res$ll_null <- -1 / 2 / sigma_null^2 * res$ssr_null - (log(sigma_null) + log(2 * pi) / 2) * ifelse(isWeight, sum(weights), n)

  # fixef info
  if (isFixef) {
    # For the within R2
    if (!onlyFixef) {
      res$ssr_fe_only <- cpp_ssq(y_demean, weights)
      sigma <- sqrt(res$ssr_fe_only / ifelse(isWeight, sum(weights), n))
      res$ll_fe_only <- -1 / 2 / sigma^2 * res$ssr_fe_only - (log(sigma) + log(2 * pi) / 2) * ifelse(isWeight, sum(weights), n)
    }
  }

  if (verbose >= 3) gt("Post-processing")

  class(res) <- "fixest"

  do_summary <- get("do_summary", env)
  if (do_summary) {
    vcov <- get("vcov", env)
    lean <- get("lean", env)
    ssc <- get("ssc", env)
    agg <- get("agg", env)
    summary_flags <- get("summary_flags", env)

    # If lean = TRUE, 1st stage residuals are still needed for the 2nd stage
    if (isTRUE(dots$iv_call) && lean) {
      r <- res$residuals
      fv <- res$fitted.values
      fvd <- res$fitted.values_demean
    }

    res <- summary(res, vcov = vcov, agg = agg, ssc = ssc, lean = lean, summary_flags = summary_flags)

    if (isTRUE(dots$iv_call) && lean) {
      res$residuals <- r
      res$fitted.values <- fv
      res$fitted.values_demean <- fvd
    }
  }

  res
}

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
    fitted.values <- cpppar_xbeta(X[, !is_excluded, drop = FALSE], beta, nthreads)
  } else {
    # avoids copies
    beta <- as.vector(xwx_inv %*% xwy)
    fitted.values <- cpppar_xbeta(X, beta, nthreads)
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
feols.fit <- function(y, X, fixef_df, vcov, offset, split, fsplit, cluster, se, ssc, weights,
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
      split = split, fsplit = fsplit, cluster = cluster, se = se, ssc = ssc,
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
