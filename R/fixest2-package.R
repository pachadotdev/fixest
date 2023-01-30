#' @useDynLib fixest2, .registration = TRUE
NULL

#' @importFrom dreamerr check_arg check_arg_plus check_value check_value_plus
#'  set_up
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
#'
#'
#'
"trade"
