#' Design matrix of a \code{fixest} object
#'
#' This function creates the left-hand-side or the right-hand-side(s) of a \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}} estimation.
#'
#' @method model.matrix fixest
#'
#' @inheritParams nobs.fixest
#'
#' @param data If missing (default) then the original data is obtained by evaluating the \code{call}. Otherwise, it should be a \code{data.frame}.
#' @param type Character vector or one sided formula, default is "rhs". Contains the type of matrix/data.frame to be returned. Possible values are: "lhs", "rhs", "fixef", "iv.rhs1" (1st stage RHS), "iv.rhs2" (2nd stage RHS), "iv.endo" (endogenous vars.), "iv.exo" (exogenous vars), "iv.inst" (instruments).
#' @param na.rm Default is \code{TRUE}. Should observations with NAs be removed from the matrix?
#' @param subset Logical or character vector. Default is \code{FALSE}. If \code{TRUE}, then the matrix created will be restricted only to the variables contained in the argument \code{data}, which can then contain a subset of the variables used in the estimation. If a character vector, then only the variables matching the elements of the vector via regular expressions will be created.
#' @param as.matrix Logical scalar, default is \code{FALSE}. Whether to coerce the result to a matrix.
#' @param as.df Logical scalar, default is \code{FALSE}. Whether to coerce the result to a data.frame.
#' @param collin.rm Logical scalar, default is \code{TRUE}. Whether to remove variables that were found to be collinear during the estimation. Beware: it does not perform a collinearity check.
#' @param ... Not currently used.
#'
#' @return
#' It returns either a vector, a matrix or a data.frame. It returns a vector for the dependent variable ("lhs"), a data.frame for the fixed-effects ("fixef") and a matrix for any other type.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{formula.fixest}}, \code{\link[fixest]{update.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}.
#'
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' est = feols(y ~ poly(x1, 2) + x2, base)
#' head(model.matrix(est))
#'
#' # Illustration of subset
#'
#' # subset => character vector
#' head(model.matrix(est, subset = "x1"))
#'
#' # subset => TRUE, only works with data argument!!
#' head(model.matrix(est, data = base[, "x1", drop = FALSE], subset = TRUE))
#'
#'
#'
model.matrix.fixest = function(object, data, type = "rhs", na.rm = TRUE, subset = FALSE,
                               as.matrix = FALSE, as.df = FALSE, collin.rm = TRUE, ...){
    # We evaluate the formula with the past call
    # type: lhs, rhs, fixef, iv.endo, iv.inst, iv.rhs1, iv.rhs2
    # if fixef => return a DF

    # Checking the arguments
    if(is_user_level_call()){
        validate_dots(suggest_args = c("data", "type"))
    }

    # We allow type to be used in the location of data if data is missing
    if(!missing(data) && missing(type)){
        sc = sys.call()
        if(!"data" %in% names(sc)){
            if(!is.null(data) && (is.character(data) || "formula" %in% class(data))){
                # data is in fact the type
                type = data
                data = NULL
            }
        }
    }


    type = check_set_types(type, c("lhs", "rhs", "fixef", "iv.endo", "iv.inst", "iv.exo", "iv.rhs1", "iv.rhs2"))

    if(isTRUE(object$is_fit)){
        stop("model.matrix method not available for fixest estimations obtained from fit methods.")
    }

    if(any(grepl("^iv", type)) && !isTRUE(object$iv)){
        stop("The type", enumerate_items(grep("^iv", type, value = TRUE), "s.is"), " only valid for IV estimations.")
    }

    check_arg(subset, "logical scalar | character vector no na")

    check_arg_plus(as.matrix, as.df, collin.rm, "logical scalar")

    # The formulas
    fml_full = formula(object, type = "full")
    fml_linear = formula(object, type = "linear")

    # Evaluation with the data
    original_data = FALSE
    if(missnull(data)){
        original_data = TRUE

        data = fetch_data(object, "To apply 'model.matrix.fixest', ")

    }

    # control of the data
    if(is.matrix(data)){
        if(is.null(colnames(data))){
            stop("If argument 'data' is to be a matrix, its columns must be named.")
        }
        data = as.data.frame(data)
    }

    if(!"data.frame" %in% class(data)){
        stop("The argument 'data' must be a data.frame or a matrix.")
    }

    data = as.data.frame(data)

    # Panel setup
    panel__meta__info = set_panel_meta_info(object, data)

    res = list()

    if("lhs" %in% type){
        lhs = list()

        namesLHS = all.vars(fml_linear[[2]])
        if(length(pblm <- setdiff(namesLHS, names(data)))){
            stop("In 'model.matrix', to create the LHS, the variable", enumerate_items(pblm, "s.is.quote"), " not in the data set.")
        }

        lhs_text = deparse_long(fml_linear[[2]])
        lhs[[lhs_text]] = eval(fml_linear[[2]], data)

        res[["lhs"]] = as.data.frame(lhs)
    }

    if("rhs" %in% type && !isTRUE(object$onlyFixef)){
        # we kick out the intercept if there is presence of fixed-effects
        fake_intercept = !is.null(object$fixef_vars) && !(!is.null(object$slope_flag) && all(object$slope_flag < 0))

        fml = fml_linear
        if(isTRUE(object$iv)){
            fml_iv = object$fml_all$iv
            fml = .xpd(..lhs ~ ..endo + ..rhs, ..lhs = fml[[2]], ..endo = fml_iv[[2]], ..rhs = fml[[3]])
        }

        linear.mat = error_sender(fixest_model_matrix_extra(
            object = object, newdata = data, original_data = original_data,
            fml = fml, fake_intercept = fake_intercept,
            subset = subset),
            "In 'model.matrix', the RHS could not be evaluated: ")

        if(collin.rm){
            qui = which(colnames(linear.mat) %in% object$collin.var)
            if(length(qui) == ncol(linear.mat)){
                linear.mat = NULL
            } else if(length(qui) > 0){
                linear.mat =  linear.mat[, -qui, drop = FALSE]
            }

            coefs = object$coefficients
            if(length(coefs) == ncol(linear.mat) && any(colnames(linear.mat) != names(coefs))){
                # we reorder the matrix
                # This can happen in multiple estimations, where we respect the
                # order of the user

                if(all(names(coefs) %in% colnames(linear.mat))){
                    linear.mat = linear.mat[, names(coefs), drop = FALSE]
                }
            }
        }

        res[["rhs"]] = linear.mat
    }

    if("fixef" %in% type){

        if(is.null(object$fixef_vars)){
            stop("In model.matrix, the type 'fixef' is only valid for models with fixed-effects. This estimation does not contain fixed-effects.")
        }

        fixef_terms_full = fixef_terms(object$fml_all$fixef)
        fixef_terms = fixef_terms_full$fml_terms

        fixef_df = error_sender(prepare_df(fixef_terms_full$fe_vars, data, fastCombine = FALSE),
                                "In 'model.matrix', problem evaluating the fixed-effects part of the formula:\n")

        isSlope = any(fixef_terms_full$slope_flag != 0)
        if(isSlope){
            slope_df = error_sender(prepare_df(fixef_terms_full$slope_vars, data),
                                    "In 'model.matrix', problem evaluating the variables with varying slopes in the fixed-effects part of the formula:\n")

            fixef_df = cbind(fixef_df, slope_df)
        }

        res[["fixef"]] = fixef_df
    }

    if("iv.endo" %in% type){
        fml = object$iv_endo_fml

        endo.mat = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = TRUE), "In 'model.matrix', the endogenous variables could not be evaluated: ")

        if(collin.rm){
            qui = which(colnames(endo.mat) %in% object$collin.var)
            if(length(qui) == ncol(endo.mat)){
                endo.mat = NULL
            } else if(length(qui) > 0){
                endo.mat =  endo.mat[, -qui, drop = FALSE]
            }
        }

        res[["iv.endo"]] = endo.mat
    }

    if("iv.inst" %in% type){
        fml = object$fml_all$iv

        inst.mat = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = TRUE), "In 'model.matrix', the instruments could not be evaluated: ")

        if(collin.rm){
            qui = which(colnames(inst.mat) %in% object$collin.var)
            if(length(qui) == ncol(inst.mat)){
                inst.mat = NULL
            } else if(length(qui) > 0){
                inst.mat =  inst.mat[, -qui, drop = FALSE]
            }
        }

        res[["iv.inst"]] = inst.mat
    }

    if("iv.exo" %in% type){

        fake_intercept = !is.null(object$fixef_vars) && !(!is.null(object$slope_flag) && all(object$slope_flag < 0))
        fml = object$fml_all$linear

        exo.mat = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = fake_intercept), "In 'model.matrix', the instruments could not be evaluated: ")

        if(is.atomic(exo.mat) && length(exo.mat) == 1){
            # This is the intercept only
            # Two cases:
            is_int = attr(terms(fml), "intercept")
            if(is_int && is.null(object$fixef_vars)){
                # Valid intercept
                exo.mat = matrix(1, nrow(data))
            } else {
                # should be NULL
                exo.mat = NULL
            }
        } else if(collin.rm){
            qui = which(colnames(exo.mat) %in% object$collin.var)
            if(length(qui) == ncol(exo.mat)){
                exo.mat = NULL
            } else if(length(qui) > 0){
                exo.mat =  exo.mat[, -qui, drop = FALSE]
            }
        }

        res[["iv.exo"]] = exo.mat
    }

    if("iv.rhs1" %in% type){
        # First stage

        if(!isTRUE(object$iv)){
            stop("In model.matrix, the type 'iv.rhs1' is only valid for IV models. This estimation is no IV.")
        }

        fml = object$fml
        if(object$iv_stage == 2){
            fml_iv = object$fml_all$iv
            fml = .xpd(..lhs ~ ..inst + ..rhs, ..lhs = fml[[2]], ..inst = fml_iv[[3]], ..rhs = fml[[3]])
        }

        fake_intercept = !is.null(object$fixef_vars) && !(!is.null(object$slope_flag) && all(object$slope_flag < 0))
        # iv_rhs1 = error_sender(fixest_model_matrix(fml, data, fake_intercept = fake_intercept),
        #                        "In 'model.matrix', the RHS of the 1st stage could not be evaluated: ")
        iv_rhs1 = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = fake_intercept, subset = subset), "In 'model.matrix', the RHS of the 1st stage could not be evaluated: ")

        if(collin.rm){
            qui = which(colnames(iv_rhs1) %in% object$collin.var)
            if(length(qui) == ncol(iv_rhs1)){
                iv_rhs1 = NULL
            } else if(length(qui) > 0){
                iv_rhs1 =  iv_rhs1[, -qui, drop = FALSE]
            }
        }

        res[["iv.rhs1"]] = iv_rhs1
    }

    if("iv.rhs2" %in% type){
        # Second stage

        if(!isTRUE(object$iv)){
            stop("In model.matrix, the type 'iv.rhs2' is only valid for second stage IV models. This estimation is not even IV.")
        }

        if(!object$iv_stage == 2){
            stop("In model.matrix, the type 'iv.rhs2' is only valid for second stage IV models. This estimation is the first stage.")
        }

        # I) we get the fit
        stage_1 = object$iv_first_stage

        fit_vars = c()
        for(i in seq_along(stage_1)){
            fit_vars[i] = v = paste0("fit_", names(stage_1)[i])
            data[[v]] = predict(stage_1[[i]], newdata = data, sample = "original")
        }

        # II) we create the variables

        fml = object$fml
        fml = .xpd(..lhs ~ ..fit + ..rhs, ..lhs = fml[[2]], ..fit = fit_vars, ..rhs = fml[[3]])

        fake_intercept = !is.null(object$fixef_vars) && !(!is.null(object$slope_flag) && all(object$slope_flag < 0))
        # iv_rhs2 = error_sender(fixest_model_matrix(fml, data, fake_intercept = fake_intercept),
        #                        "In 'model.matrix', the RHS of the 2nd stage could not be evaluated: ")
        iv_rhs2 = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = fake_intercept, subset = subset), "In 'model.matrix', the RHS of the 2nd stage could not be evaluated: ")

        if(collin.rm){
            qui = which(colnames(iv_rhs2) %in% object$collin.var)
            if(length(qui) == ncol(iv_rhs2)){
                iv_rhs2 = NULL
            } else if(length(qui) > 0){
                iv_rhs2 =  iv_rhs2[, -qui, drop = FALSE]
            }
        }

        res[["iv.rhs2"]] = iv_rhs2
    }

    # Formatting res
    if(length(res) == 0){
        return(NULL)
    } else if(length(type) > 1){
        res = res[type]
        res = do.call(cbind, unname(res))
    } else {
        res = res[[1]]
    }

    #
    # Removing obs if needed
    #

    check_0 = FALSE
    if(original_data){

        if(na.rm == FALSE){
            # We do nothing. Or shall I add NA values for obs not
            # included in the estimation?
            if(FALSE && length(object$obs_selection) > 0){

                # we reconstruct the full vector of obs
                # and we fill with NA
                obs_id = 1:nrow(data)
                for(i in seq_along(object$obs_selection)){
                    obs_id = select_obs(obs_id, object$obs_selection[[i]])
                }

                res[!1:nrow(res) %in% obs_id, ] = NA

            }

        } else {
            for(i in seq_along(object$obs_selection)){
                check_0 = TRUE
                res = select_obs(res, object$obs_selection[[i]])
            }
        }



        na.rm = FALSE
    }

    if(na.rm){

        if(is.numeric(res) || all(sapply(res, is.numeric))){
            info = cpp_which_na_inf(res, nthreads = 1)
        } else {
            info = list(any_na_inf = anyNA(res))
            if(info$any_na_inf) info$is_na_inf = !complete.cases(res)
        }

        if(info$any_na_inf){
            check_0 = TRUE
            isNA_L = info$is_na_inf

            if(sum(isNA_L) == nrow(res)){
                warning("All observations contain NA values.")
                return(res[-which(isNA_L), , drop = FALSE])
            }

            res = select_obs(res, -which(isNA_L))
        }
    }


    if(as.matrix){
        res = as.matrix(res)
    } else if(as.df){
        res = as.data.frame(res)
    } else if(identical(type, "lhs")){
        res = res[[1]]
    }

    if(check_0 && !"fixef" %in% type){
        only_0 = cpppar_check_only_0(base::as.matrix(res), nthreads = 1)
        if(all(only_0 == 1)){
            stop("After removing NAs, not a single explanatory variable is different from 0.")

        } else if(any(only_0 == 1)){
            # At that point it must be either a matrix or a DF
            # (can't be a vector)
            res = res[, only_0 == 0, drop = FALSE]
        }
    }

    res
}
