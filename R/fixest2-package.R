#' Fast and User-Friendly Fixed-Effects Estimations
#'
#' The package \pkg{fixest} provides a family of functions to perform estimations
#' with multiple fixed-effects. Standard-errors can be easily and intuitively clustered.
#' It also includes tools to seamlessly export the results of various estimations.
#' * To get started, look at the [introduction](https://cran.r-project.org/package=fixest/vignettes/fixest_walkthrough.html).
#'
#' The main features are:
#' * Estimation. The core functions are: [`feols`], [`feglm`] and [`femlm`] to
#' estimate, respectively, linear models, generalized linear models and maximum likelihood
#' models with multiple fixed-effects. The function [`feNmlm`] allows the inclusion of
#' non-linear in parameters right hand sides. Finally [`fepois`][fixest::feglm]
#' and [`fenegbin`][fixest::femlm] are shorthands to estimate Poisson and Negative
#' Binomial models.
#' * Multiple estimations: You can perform multiple estimations at once with
#' the [`stepwise`] functions. It's then very easy to manipulate multiple results
#' with the associated methods. See an introduction in the dedicated vignette:
#' [Multiple estimations](https://cran.r-project.org/package=fixest/vignettes/multiple_estimations.html)
#' * Easy and flexible clustering of standard-errors. By using the arguments `vcov`
#' and `ssc` (see [`summary.fixest`]). To have a sense of how the standard errors are computed,
#' see the vignette [On standard-errors](https://lrberge.github.io/fixest/articles/standard_errors.html).
#' * Visualization and exportation of results. You can visualize the results of
#' multiple estimations in R, or export them in Latex using the function [`etable`].
#' This vignette details how to customize the Latex tables: [Exporting estimation tables](https://lrberge.github.io/fixest/articles/exporting_tables.html).
#' * Plot multiple results. You can plot the coefficients and confidence intervals of
#' estimations easily with the function [`coefplot`]. This function also offers a specific
#' layout for interactions.
#'
#'
#' @references
#' Berge, Laurent, 2018, "Efficient estimation of maximum likelihood models
#' with multiple fixed-effects: the R package FENmlm." CREA Discussion Papers,
#' 13 ([](https://github.com/lrberge/fixest/blob/master/_DOCS/FENmlm_paper.pdf)).
#'
#' @keywords internal
#'
"_PACKAGE"

#' @useDynLib fixest2, .registration = TRUE
NULL

#' @importFrom dreamerr check_arg check_arg_plus check_value check_value_plus
#'  set_up stop_up validate_dots sfill plural plural_len ifunit warn_up n_th
#'  check_arg fit_screen ifsingle enumerate_items check_set_arg check_set_value
NULL

#' @importFrom stats .getXlevels AIC BIC aggregate as.formula binomial cancor
#'  coef complete.cases confint deviance formula logLik median model.frame
#'  model.matrix na.omit na.pass nlminb nobs pchisq pf pnorm poisson predict
#'  pt qnorm qt quantile resid rnorm sd setNames terms update var
NULL

#' @importFrom utils as.roman capture.output combn head object.size
#'  packageVersion read.csv tail write.csv
NULL

#' @importFrom nlme fixef
#' @export
NULL

#' @importFrom sandwich estfun bread
#' @export
NULL

#' @importFrom stringmagic string_any
NULL

#' Trade data sample
#'
#' This data reports trade information between countries of the European Union (EU15).
#'
#' @usage
#' data(trade)
#'
#' @format
#' \code{trade} is a data frame with 38,325 observations and 6 variables named \code{Destination}, \code{Origin}, \code{Product}, \code{Year}, \code{dist_km} and \code{Euros}.
#'
#' \itemize{
#' \item{Origin: 2-digits codes of the countries of origin of the trade flow.}
#' \item{Destination: 2-digits codes of the countries of destination of the trade flow.}
#' \item{Products: Number representing the product categories (from 1 to 20).}
#' \item{Year: Years from 2007 to 2016}
#' \item{dist_km: Geographic distance in km between the centers of the countries of origin and destination.}
#' \item{Euros: The total amount in euros of the trade flow for the specific year/product category/origin-destination country pair.}
#'
#' }
#'
#' @source
#' This data has been extrated from Eurostat on October 2017.
"trade"

#' @importFrom stringmagic string_magic_alias cat_magic_alias message_magic_alias
NULL
