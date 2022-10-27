#' Obtain various statistics from an estimation
#'
#' Set of functions to directly extract some commonly used statistics, like the p-value or the table of coefficients, from estimations. This was first implemented for \code{fixest} estimations, but has some support for other models.
#'
#' @inheritParams etable
#'
#' @param object An estimation. For example obtained from \code{\link[fixest]{feols}}.
#' @param se [Fixest specific.] Character scalar. Which kind of standard error should be computed: \dQuote{iid}, \dQuote{hetero}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}? By default if there are fixed-effects in the estimation: \code{se = "cluster"}, otherwise \code{se = "iid"}. Note that this argument is not needed if the argument \code{cluster} is present.
#' @param cluster [Fixest specific.] Tells how to cluster the standard-errors (if clustering is requested). Can be either a list of vectors, a character vector of variable names, a formula or an integer vector. Assume we want to perform 2-way clustering over \code{var1} and \code{var2} contained in the data.frame \code{base} used for the estimation. All the following \code{cluster} arguments are valid and do the same thing: \code{cluster = base[, c("var1, "var2")]}, \code{cluster = c("var1, "var2")}, \code{cluster = ~var1+var2}. If the two variables were used as clusters in the estimation, you could further use \code{cluster = 1:2} or leave it blank with \code{se = "twoway"} (assuming \code{var1} [resp. \code{var2}] was the 1st [res. 2nd] cluster).
#' @param ... Other arguments to be passed to \code{summary}.
#'
#' @details
#' This set of functions is primarily constructed for \code{fixest} estimations. Although it can work for non-\code{fixest} estimations, support for exotic estimation procedures that do not report standardized coefficient tables is highly limited.
#'
#' @return
#' Returns a table of coefficients, with in rows the variables and four columns: the estimate, the standard-error, the t-statistic and the p-value.
#'
#' @examples
#'
#' # Some data and estimation
#' data(trade)
#' est = fepois(Euros ~ log(dist_km) | Origin^Product + Year, trade)
#'
#' #
#' # Coeftable/se/tstat/pvalue
#' #
#'
#' # Default is clustering along Origin^Product
#' coeftable(est)
#' se(est)
#' tstat(est)
#' pvalue(est)
#'
#' # Now with two-way clustered standard-errors
#' #  and using coeftable()
#'
#' coeftable(est, cluster = ~Origin + Product)
#' se(est, cluster = ~Origin + Product)
#' pvalue(est, cluster = ~Origin + Product)
#' tstat(est, cluster = ~Origin + Product)
#'
#' # Or you can cluster only once:
#' est_sum = summary(est, cluster = ~Origin + Product)
#' coeftable(est_sum)
#' se(est_sum)
#' tstat(est_sum)
#' pvalue(est_sum)
#'
#' # You can use the arguments keep, drop, order
#' # to rearrange the results
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' est_iv = feols(y ~ x1 | x2 ~ x3, base)
#'
#' tstat(est_iv, keep = "x1")
#' coeftable(est_iv, keep = "x1|Int")
#'
#' coeftable(est_iv, order = "!Int")
#'
#'
#'
coeftable = function(object, vcov = NULL, ssc = NULL, cluster = NULL, keep, drop, order, ...){
    # We don't explicitly refer to the other arguments

    check_arg(keep, drop, order, "NULL character vector no na")

    IS_FIXEST = "fixest" %in% class(object)
    # IS_FIXEST_MULI = "fixest_multi" %in% class(object)
    # IS_FIXEST_LIST = "fixest_multi" %in% class(object)

    if(IS_FIXEST){
        if(!isTRUE(object$summary) || !(all_missing(vcov, ssc, cluster) && ...length() > 0)){
            object = summary(object, vcov = vcov, ssc = ssc, cluster = cluster, ...)
        }
    } else {
        if(!any(grepl("summary", class(object)))){
            # Known issue:
            # - if vcov/.../order are also args of summary, that will be a problem

            args = list(object = object)

            if(...length() > 0){
                args = append(args, list(...))
            }

            object = do.call("summary", args)
        }
    }

    # Let's find out the coefficients table
    if(IS_FIXEST){
        res = object$coeftable
    } else {
        list_mat = object[sapply(object, is.matrix)]

        ok = FALSE
        for(i in seq_along(list_mat)){
            mat = list_mat[[i]]
            if(!is.null(colnames(mat)) && any(grepl("(?i)(estimate|value|Pr\\()", colnames(mat)))){
                ok = TRUE
                res = mat
            }
        }

        if(ok == FALSE){
            stop("Sorry, the coeffficients table could not be extracted.")
        }

    }

    if(!missnull(keep) || !missnull(drop) || !missnull(order)){
        r_names = rownames(res)
        r_names = keep_apply(r_names, keep)
        r_names = drop_apply(r_names, drop)
        r_names = order_apply(r_names, order)

        if(length(r_names) == 0){
            return(NULL)
        }

        res = res[r_names, , drop = FALSE]
    }

    res
}

#' @describeIn coeftable Extracts the p-value of an estimation
pvalue = function(object, vcov = NULL, ssc = NULL, cluster = NULL, keep, drop, order, ...){

    check_arg(keep, drop, order, "NULL character vector no na")

    mat = coeftable(object, vcov = vcov, ssc = ssc, cluster = cluster, ...)

    if(ncol(mat) != 4){
        stop("No appropriate coefficient table found (number of columns is ", ncol(mat), " instead of 4), sorry.")
    }

    res = mat[, 4]
    if(is.null(names(res))) {
        names(res) = rownames(mat)
    }

    if(!missnull(keep) || !missnull(drop) || !missnull(order)){
        r_names = names(res)
        r_names = keep_apply(r_names, keep)
        r_names = drop_apply(r_names, drop)
        r_names = order_apply(r_names, order)

        if(length(r_names) == 0){
            return(numeric(0))
        }

        res = res[r_names]
    }

    res
}

#' @describeIn coeftable Extracts the t-statistics of an estimation
tstat = function(object, vcov = NULL, ssc = NULL, cluster = NULL, keep, drop, order, ...){

    check_arg(keep, drop, order, "NULL character vector no na")

    mat = coeftable(object, vcov = vcov, ssc = ssc, cluster = cluster, ...)

    if(ncol(mat) != 4){
        stop("No appropriate coefficient table found (number of columns is ", ncol(mat), " instead of 4), sorry.")
    }

    res = mat[, 3]
    if(is.null(names(res))) {
        names(res) = rownames(mat)
    }

    if(!missnull(keep) || !missnull(drop) || !missnull(order)){
        r_names = names(res)
        r_names = keep_apply(r_names, keep)
        r_names = drop_apply(r_names, drop)
        r_names = order_apply(r_names, order)

        if(length(r_names) == 0){
            return(numeric(0))
        }

        res = res[r_names]
    }

    res
}

#' @describeIn coeftable Extracts the standard-error of an estimation
se = function(object, vcov = NULL, ssc = NULL, cluster = NULL, keep, drop, order, ...){

    check_arg(keep, drop, order, "NULL character vector no na")

    if(is.matrix(object) && nrow(object) == ncol(object)){
        # special case => object is a VCOV matrix and NOT an estimation
        res = sqrt(diag(object))

    } else {

        mat = coeftable(object, vcov = vcov, ssc = ssc, cluster = cluster, ...)

        if(ncol(mat) != 4){
            stop("No appropriate coefficient table found (number of columns is ", ncol(mat), " instead of 4), sorry.")
        }

        res = mat[, 2]
        if(is.null(names(res))) {
            names(res) = rownames(mat)
        }
    }

    if(!missnull(keep) || !missnull(drop) || !missnull(order)){
        r_names = names(res)
        r_names = keep_apply(r_names, keep)
        r_names = drop_apply(r_names, drop)
        r_names = order_apply(r_names, order)

        if(length(r_names) == 0){
            return(numeric(0))
        }

        res = res[r_names]
    }

    res
}
