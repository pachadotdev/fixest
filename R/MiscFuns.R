#' Extract the Fixed-Effects from a \code{fixest} estimation.
#'
#' This function retrieves the fixed effects from a \code{fixest} estimation. It is useful only when there are one or more fixed-effect dimensions.
#'
#' @inheritParams feNmlm
#'
#' @param object A \code{fixest} estimation (e.g. obtained using \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}).
#' @param notes Logical. Whether to display a note when the fixed-effects coefficients are not regular.
#' @param sorted Logical, default is \code{TRUE}. Whether to order the fixed-effects by their names. If \code{FALSE}, then the order used in the demeaning algorithm is used.
#'
#' @details
#' If the fixed-effect coefficients not regular, then several reference points need to be set, leading to the coefficients to be NOT interpretable. If this is the case, then a warning is raised.
#'
#' @return
#' A list containing the vectors of the fixed effects.
#'
#' If there is more than 1 fixed-effect, then the attribute \dQuote{references} is created. This is a vector of length the number of fixed-effects, each element contains the number of coefficients set as references. By construction, the elements of the first fixed-effect dimension are never set as references. In the presence of regular fixed-effects, there should be Q-1 references (with Q the number of fixed-effects).
#'
#' @seealso
#' \code{\link[fixest]{plot.fixest.fixef}}. See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effect coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' data(trade)
#'
#' # We estimate the effect of distance on trade => we account for 3 fixed-effects
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # Obtaining the fixed-effects coefficients:
#' fe_trade = fixef(est_pois)
#'
#' # The fixed-effects of the first fixed-effect dimension:
#' head(fe_trade$Origin)
#'
#' # Summary information:
#' summary(fe_trade)
#'
#' # Plotting them:
#' plot(fe_trade)
#'
fixef.fixest = function(object, notes = getFixest_notes(), sorted = TRUE, nthreads = getFixest_nthreads(),
                        fixef.tol = 1e-5, fixef.iter = 10000, ...){

    # object is a fixest object
    # This function retrieves the dummies

    check_arg(notes, sorted, "logical scalar")

    # Checking the arguments
    if(is_user_level_call()){
        validate_dots()
    }

    check_value(fixef.tol, "numeric scalar GT{0} LT{1}")
    check_value(fixef.iter, "strict integer scalar GT{0}")

    if(isTRUE(object$lean)){
        # LATER: recompute the FEs by extracting them from the data
        stop("Fixed-effects from 'lean' fixest objects cannot be extracted. Please re-estimate with 'lean = FALSE'.")
    }

    # Preliminary stuff
    S = object$sumFE

    if(is.null(S)){
        stop("The estimation was done without fixed-effects (FE). The FE coefficients cannot be retrieved.")
    }

    family = object$family
    fixef_names = object$fixef_vars

    fixef_id = object$fixef_id

    Q = length(fixef_id)
    N = length(S)

    # either (we need to clean its attributes for unlist to be efficient)
    id_dummies_vect = list()
    for(i in 1:Q) id_dummies_vect[[i]] = as.vector(fixef_id[[i]])

    is_ref_approx = FALSE
    isSlope = FALSE
    if(!is.null(object$fixef_terms)){
        isSlope = TRUE
        # This is an estimation with slopes
        # we apply another method => we use the demeaning function

        slope_variables = object$slope_variables_reordered
        slope_flag = object$slope_flag_reordered

        new_order = object$fe.reorder
        fixef_vars = object$fixef_vars[new_order]
        fixef_sizes = as.integer(object$fixef_sizes[new_order])

        # We reconstruct the terms
        fixef_terms = c()
        start = c(0, cumsum(abs(slope_flag)))
        for(i in seq_along(slope_flag)){
            sf = slope_flag[i]
            if(sf >= 0){
                fixef_terms = c(fixef_terms, fixef_vars[i])
            }

            if(abs(sf) > 0){
                fixef_terms = c(fixef_terms, paste0(fixef_vars[i], "[[", names(slope_variables)[start[i] + 1:abs(sf)], "]]"))
            }
        }

        fe_id_list = object$fixef_id[new_order]

        #
        # STEP 2: demeaning
        #

        nthreads = check_set_nthreads(nthreads)

        table_id_I = as.integer(unlist(lapply(fe_id_list, table), use.names = FALSE))

        S_demean = cpp_demean(y = S, X_raw = 0, r_weights = 0, iterMax = as.integer(fixef.iter),
                              diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
                              fe_id_list = fe_id_list, table_id_I = table_id_I,
                              slope_flag_Q = slope_flag, slope_vars_list = slope_variables,
                              r_init = 0, nthreads = nthreads, save_fixef = TRUE)

        fixef_coef = S_demean$fixef_coef

        names(fixef_sizes) = fixef_vars

        fe_all = c()
        for(i in seq_along(slope_flag)){
            fe_all = c(fe_all, rep(fixef_vars[i], 1 + abs(slope_flag[i]) - (slope_flag[i] < 0)))
        }

        start = 1
        i = 1
        fixef_values = list()
        for(q in seq_along(slope_flag)){
            sf = slope_flag[q]
            if(sf == 0){
                fixef_values[[i]] = fixef_coef[seq(start, length.out = fixef_sizes[q])]
                i = i + 1
                start = start + fixef_sizes[q]
            } else {
                nb = abs(sf) + (sf > 0)

                adj = 0
                if(sf > 0){
                    # The fixed-effects is in the last position
                    j_fe = nb - 1
                    fixef_values[[i]] = fixef_coef[seq(start + j_fe, by = nb, length.out = fixef_sizes[q])]
                    adj = 1
                }

                for(j in 0:(nb - 1 - adj)){
                    fixef_values[[i + j + adj]] = fixef_coef[seq(start + j, by = nb, length.out = fixef_sizes[q])]
                }
                i = i + nb
                start = start + fixef_sizes[q] * nb
            }

        }

        #
        # Now the referenes
        #

        nb_ref = integer(length(fixef_terms))

        # FE references
        who_fe = slope_flag >= 0
        Q_fe = sum(who_fe)
        if(Q_fe >= 2){

            my_dum = fe_id_list[who_fe]

            dumMat = matrix(unlist(my_dum, use.names = FALSE), N, Q_fe) - 1
            orderCluster = matrix(unlist(lapply(my_dum, order), use.names = FALSE), N, Q_fe) - 1

            nbCluster = sapply(my_dum, max)

            fixef_values_tmp = cpp_get_fe_gnl(Q_fe, N, rep(1, N), dumMat, nbCluster, orderCluster)

            # the information on the references
            nb_ref_fe = fixef_values_tmp[[Q_fe+1]]
        } else {
            nb_ref_fe = integer(Q_fe)
        }

        # Slope references (if associated FE + constant)

        names(slope_flag) = fixef_vars

        Q_slope = sum(abs(slope_flag))
        nb_ref_slope = integer(Q_slope)
        i_noVS = 1
        for(i in seq_along(fixef_terms)){

            ft = fixef_terms[i]

            if(!grepl("[[", ft, fixed = TRUE)){
                # No slope => already computed
                nb_ref[i] = nb_ref_fe[i_noVS]
                i_noVS = i_noVS + 1

            } else {
                # Slope
                fe_name = gsub("\\[.+", "", ft)
                my_dum = fe_id_list[[fe_name]]

                my_order = order(my_dum)
                var_sorted = slope_variables[[gsub(".+\\[|\\]+", "", ft)]][my_order]

                # if no associated FE => we check only 0 values
                if(slope_flag[fe_name] < 0){
                    nb_ref[i] = cpp_constant_dum(fixef_sizes[fe_name], var_sorted, my_dum[my_order], only_0 = TRUE)
                } else {
                    nb_ref[i] = cpp_constant_dum(fixef_sizes[fe_name], var_sorted, my_dum[my_order])
                }
            }
        }

        # we recreate that to avoid conditioning on isSlope later
        fixef_id = fixef_id[fe_all]
        fixef_names = fixef_terms

    } else if(Q == 1){
        # This is the simplest case
        id = id_dummies_vect[[1]]

        myOrder = order(id)
        myDiff = c(1, diff(id[myOrder]))

        select = myOrder[myDiff == 1]

        fixef_values = list(S[select])

        # There are no references => no need to set nb_ref
    } else {
        # We apply a Rcpp script to handle complicated cases (and we don't know beforehand if the input is one)

        dumMat = matrix(unlist(id_dummies_vect), N, Q) - 1
        orderCluster = matrix(unlist(lapply(id_dummies_vect, order)), N, Q) - 1

        nbCluster = sapply(fixef_id, max)

        fixef_values = cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)

        # The algorithm is fast but may fail on some instance. We need to check
        if(any(fixef_values[[Q + 1]] > 1) && Q >= 3){
            # we re-compute the "sumFE"
            sum_FE = fixef_values[[1]][1 + dumMat[, 1]]
            for(i in 2:Q){
                sum_FE = sum_FE + fixef_values[[i]][1 + dumMat[, i]]
            }

            if(max(abs(sum_FE - S)) > 1e-1){
                # divergence => we need to correct
                # we recompute the FEs

                is_ref_approx = TRUE

                fixef_sizes = as.integer(object$fixef_sizes)
                new_order = order(object$fixef_sizes, decreasing = TRUE)
                fixef_sizes = fixef_sizes[new_order]

                fe_id_list = object$fixef_id[new_order]
                table_id_I = as.integer(unlist(lapply(fe_id_list, table), use.names = FALSE))

                slope_flag = rep(0L, Q)
                slope_variables = list()

                S_demean = cpp_demean(y = S, X_raw = 0, r_weights = 0, iterMax = as.integer(fixef.iter),
                                      diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
                                      fe_id_list = fe_id_list, table_id_I = table_id_I,
                                      slope_flag_Q = slope_flag, slope_vars_list = slope_variables,
                                      r_init = 0, nthreads = nthreads, save_fixef = TRUE)

                fixef_coef = S_demean$fixef_coef

                end = cumsum(fixef_sizes)
                start = c(1, end + 1)
                for(i in 1:Q){
                    fixef_values[[new_order[i]]] = fixef_coef[start[i]:end[i]]
                }
            }
        }

        # the information on the references
        nb_ref = fixef_values[[Q + 1]]
        fixef_values[[Q + 1]] = NULL
    }

    # now saving & adding information
    all_clust = list()
    Q_all = ifelse(isSlope, length(fixef_terms), Q)
    for(i in 1:Q_all){
        # We put it in the right order, if requested
        fn = attr(fixef_id[[i]], "fixef_names")

        if(sorted){
            if(all(!grepl("[^[:digit:]]", fn))) fn = as.numeric(fn)
            my_order = order(fn)

            cv = fixef_values[[i]][my_order]
            names(cv) = fn[my_order]
            all_clust[[fixef_names[i]]] = cv
        } else {
            cv = fixef_values[[i]]
            names(cv) = fn
            all_clust[[fixef_names[i]]] = cv
        }

    }

    class(all_clust) = c("fixest.fixef", "list")

    # Dealing with the references
    if(Q_all > 1){
        names(nb_ref) = fixef_names
        attr(all_clust, "references") = nb_ref

        if(!isSlope) slope_flag = rep(FALSE, Q)

        # warning if unbalanced
        if(notes && sum(nb_ref[!slope_flag]) > Q-1){
            if(is_ref_approx){
                message("NOTE: The fixed-effects are not regular, they cannot be straightforwardly interpreted. The number of references is only approximate.")
            } else {
                message("NOTE: The fixed-effects are not regular, they cannot be straightforwardly interpreted.")
            }

        }
    }

    # Family information
    attr(all_clust, "exponential") = FALSE
    if(object$method_type == "feNmlm" && object$family %in% c("poisson", "negbin")){
        attr(all_clust, "exponential") = TRUE
    } else if(object$method_type == "feglm" && object$family$link == "log"){
        attr(all_clust, "exponential") = TRUE
    }

    return(all_clust)
}

#' Functions exported from \pkg{nlme} to implement \pkg{fixest} methods
#'
#' The package \pkg{fixest} uses the \code{fixef} method from \pkg{nlme}. Unfortunately, re-exporting this method is required in order not to attach package \pkg{nlme}.
#'
#' \itemize{
#' \item Here is the help from package \pkg{nlme}: \code{\link[nlme:fixed.effects]{fixef}}. The help from package \pkg{fixest} is here: \code{\link[fixest]{fixef.fixest}}.
#' }
#'
#' @note
#' I could find this workaround thanks to the package \pkg{plm}.
#'
#' @name fixef_reexported
#' @keywords internal
NULL

#' @rdname fixef_reexported
#' @name fixef
NULL




#' Displaying the most notable fixed-effects
#'
#' This function plots the 5 fixed-effects with the highest and lowest values, for each of the fixed-effect dimension. It takes as an argument the fixed-effects obtained from the function \code{\link{fixef.fixest}} after an estimation using \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @method plot fixest.fixef
#'
#' @param x An object obtained from the function \code{\link{fixef.fixest}}.
#' @param n The number of fixed-effects to be drawn. Defaults to 5.
#' @param ... Not currently used.
#'
#' Note that the fixed-effect coefficients might NOT be interpretable. This function is useful only for fully regular panels.
#'
#' If the data are not regular in the fixed-effect coefficients, this means that several \sQuote{reference points} are set to obtain the fixed-effects, thereby impeding their interpretation. In this case a warning is raised.
#'
#' @seealso
#' \code{\link[fixest]{fixef.fixest}} to extract clouster coefficients. See also the main estimation function \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' data(trade)
#'
#' # We estimate the effect of distance on trade
#' # => we account for 3 fixed-effects
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # obtaining the fixed-effects coefficients
#' fe_trade = fixef(est_pois)
#'
#' # plotting them
#' plot(fe_trade)
#'
#'
plot.fixest.fixef = function(x, n = 5, ...){

    # Checking the arguments
    if(is_user_level_call()){
        validate_dots(suggest_args = "n")
    }

    Q = length(x)

    mfrow = as.character(c(11, 12, 22, 22, 32, 32, 33, 33))

    fixef_names = names(x)
    slope_flag = grepl("\\[", fixef_names)

    if(Q > 1 && sum(attr(x, "references")[!slope_flag]) > sum(!slope_flag)-1){
        warning("The fixed-effects are not regular, they cannot be straightforwardly interpreted.", call. = FALSE)
    }

    # modification par:
    opar = par(no.readonly =TRUE)
    on.exit(par(opar))

    par(mfrow = as.numeric(strsplit(mfrow[Q], "")[[1]]), mar = c(3, 3, 2.5, 3))

    addExp = attr(x, "exponential")
    for(i in 1:Q){
        plot_single_cluster(x[[i]], n = n, addExp = addExp, fe_name = fixef_names[i])
    }

}


#' Collinearity diagnostics for \code{fixest} objects
#'
#' In some occasions, the optimization algorithm of \code{\link[fixest]{femlm}} may fail to converge, or the variance-covariance matrix may not be available. The most common reason of why this happens is colllinearity among variables. This function helps to find out which set of variables is problematic.
#'
#'
#' @param x A \code{fixest} object obtained from, e.g. functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param verbose An integer. If higher than or equal to 1, then a note is prompted at each step of the algorithm. By default \code{verbose = 0} for small problems and to 1 for large problems.
#'
#' @details
#' This function tests: 1) collinearity with the fixed-effect variables, 2) perfect multi-collinearity between the variables, 4) perfect multi-collinearity between several variables and the fixed-effects, and 4) identification issues when there are non-linear in parameters parts.
#'
#' @return
#' It returns a text message with the identified diagnostics.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Creating an example data base:
#' set.seed(1)
#' fe_1 = sample(3, 100, TRUE)
#' fe_2 = sample(20, 100, TRUE)
#' x = rnorm(100, fe_1)**2
#' y = rnorm(100, fe_2)**2
#' z = rnorm(100, 3)**2
#' dep = rpois(100, x*y*z)
#' base = data.frame(fe_1, fe_2, x, y, z, dep)
#'
#' # creating collinearity problems:
#' base$v1 = base$v2 = base$v3 = base$v4 = 0
#' base$v1[base$fe_1 == 1] = 1
#' base$v2[base$fe_1 == 2] = 1
#' base$v3[base$fe_1 == 3] = 1
#' base$v4[base$fe_2 == 1] = 1
#'
#' # Estimations:
#'
#' # Collinearity with the fixed-effects:
#' res_1 = femlm(dep ~ log(x) + v1 + v2 + v4 | fe_1 + fe_2, base)
#' collinearity(res_1)
#'
#' # => collinearity with the first fixed-effect identified, we drop v1 and v2
#' res_1bis = femlm(dep ~ log(x) + v4 | fe_1 + fe_2, base)
#' collinearity(res_1bis)
#'
#' # Multi-Collinearity:
#' res_2 =  femlm(dep ~ log(x) + v1 + v2 + v3 + v4, base)
#' collinearity(res_2)
#'
#'
collinearity = function(x, verbose){
    # x: fixest estimation

    # stop("Sorry, it does not work. A new version will hopefully come soon.")

    if(!inherits(x, "fixest")){
        stop("Argument 'x' must be a fixest object.")
    }

    # I) (linear) collinearity with fixed-effects
    # II) (linear) multi collinearity
    # III) (non-linear) overidentification

    # flags
    isFixef = !is.null(x$fixef_vars)
    isFE = FALSE
    isSlope = FALSE
    if(isFixef){
        isSlope = !is.null(x$fixef_terms)
        if(isSlope){
            fixef_terms = x$fixef_terms
            slope_flag = grepl("\\[", fixef_terms)
            isFE = any(!slope_flag)
        } else {
            isFE = TRUE
        }
    }

    linear_fml = fml_split(formula(x, "linear"), 2, split.lhs = TRUE)

    rhs_fml = fml_split(formula(x, "linear"), 1)
    linear.varnames = all.vars(rhs_fml[[3]])

    isLinear = length(linear.varnames) > 0

    NL_fml = x$NL.fml
    isNL = !is.null(NL_fml)
    coef = x$coefficients

    # Getting the data
    data = fetch_data(x, "To apply function 'collinearity', ")

    if(is.matrix(data)){
        data = as.data.frame(data)
    } else {
        class(data) = "data.frame"
    }

    if(isFE){
        linear_fml = update(linear_fml, ~ . + 1)
    }

    # Panel setup
    panel__meta__info = set_panel_meta_info(x, data)

    if(isLinear || isFixef || "(Intercept)" %in% names(coef)){
        # linear.matrix = model.matrix(linear_fml, data)
        linear.matrix = fixest_model_matrix(rhs_fml, data)
    }

    for(i in seq_along(x$obs_selection)){
        linear.matrix = linear.matrix[x$obs_selection[[i]], , drop = FALSE]
    }

    if(isLinear){
        # We do that to drop interaction variables that should not be there any more
        # if factors with only NA values
        varkeep = intersect(names(x$coefficients), colnames(linear.matrix))
        if(length(varkeep) < ncol(linear.matrix)){
            linear.matrix = linear.matrix[, varkeep, drop = FALSE]
        }
    }


    if(isLinear){
        # We center and scale to have something comparable across data sets
        first_row = linear.matrix[1, ]
        linear.matrix = scale(linear.matrix, center = FALSE, scale = TRUE)
        # The constants => we set it to 1
        constant_id = apply(linear.matrix, 2, anyNA)
        constant_value = first_row[constant_id]
        linear.matrix[is.na(linear.matrix)] = 1
        linear_mat_noIntercept = linear.matrix[, -1, drop = FALSE]
    }

    n_obs = nrow(data)
    Q = length(x$fixef_id)

    # information display
    ccat = function(...) NULL
    if(missing(verbose)){
        verbose = FALSE
        # verbose only if task is long
        if(100 * nrow(linear.matrix) * (ncol(linear.matrix)**2 * Q**2) >= 1e9){
            verbose = TRUE
        }
    }
    if(verbose) ccat = cat

    #
    # 0) Data preparation (FE/slopes)
    #

    # Variables for FE / Slopes
    if(isSlope){
        fixef_terms = x$fixef_terms
        terms_full = extract_fe_slope(fixef_terms)
        fixef_vars = terms_full$fixef_vars
        slope_fe = terms_full$slope_fe
        fe_all = terms_full$fe_all
        slope_vars = terms_full$slope_vars
        slope_terms = terms_full$slope_terms

        # dictionary mapping fe var names to the ids of id_dummies_vect
        dict_fe = 1:Q
        names(dict_fe) = x$fixef_vars
        slope_vars_unik = unique(slope_vars)

        if(any(!slope_vars_unik %in% names(data))){
            var_pblm = setdiff(slope_vars_unik, names(data))
            stop("To check collinearity, we need to fetch some variables in the original dataset (", deparse_long(x$call$data), "). However, the variable", enumerate_items(var_pblm, "s.is"), " not present in the original dataset any more.")
        }

        slope_var_list = list()
        for(i in 1:length(slope_vars_unik)){
            variable = all.vars(str2lang(slope_vars_unik[i]))

            # as.numeric => we'll use cpp so required
            svar = as.numeric(as.vector(eval(str2lang(slope_vars_unik[i]), data)))
            if(length(svar) == 1) svar = rep(svar, n_obs)

            slope_var_list[[slope_vars_unik[i]]] = svar
        }

        Q_slope = length(slope_fe)
        slope_fe_id = x$fixef_id[slope_fe]
        slope_var_all = slope_var_list[slope_vars]

        Q_fe = length(fixef_vars)
        fe_id = x$fixef_id[fixef_vars]

    } else if(isFE){
        Q_fe = length(x$fixef_id)
        fe_id = x$fixef_id
    }


    #
    # I) collinearity with clusters
    #

    ccat("Checking Collinearity: ")

    #
    # Variables are constant or equal to 0
    #

    if(isLinear){

        if(isFE){
            if(any(constant_id[-1])){
                # var_problem = colnames(linear_mat_noIntercept)[is_const]
                var_problem = colnames(linear_mat_noIntercept)[constant_id[-1]]
                message = paste0("Variable", enumerate_items(var_problem, "s.is.quote"), " constant, thus collinear with the fixed-effects.")
                print(message)
                return(invisible(message))
            }

        } else {

            isIntercept = attr(terms(x$fml),"intercept")

            if(isIntercept && any(constant_id[-1])){
                var_problem = colnames(linear_mat_noIntercept)[constant_id[-1]]
                message = paste0("Variable", enumerate_items(var_problem, "s.is.quote"), " constant, thus collinear with the intercept.")
                print(message)
                return(invisible(message))
            }

            if(any(constant_value == 0)){
                # constant and equal to 0
                var_problem = colnames(linear.matrix)[first_row == 0 & constant_id]
                message = paste0("Variable", enumerate_items(var_problem, "s.is.quote"), " constant and equal to 0.")

                print(message)
                return(invisible(message))
            }
        }
    }

    if(isFE && isLinear){
        ccat("simple with fixed-effects:")
        # We project each variable onto the cluster subspace

        cluster = fe_id
        Q_fe = length(cluster)
        for(q in 1:Q_fe){
            ccat(".")
            dum = cluster[[q]]
            k = max(dum)

            value = cpp_tapply_sum(Q = k, x = linear_mat_noIntercept, dum = dum)
            nb_per_cluster = cpp_table(Q = k, dum = dum)

            # residuals of the linear projection on the cluster space
            residuals = linear_mat_noIntercept - (value/nb_per_cluster)[dum, ]

            max_residuals = apply(abs(residuals), 2, max)

            if(any(max_residuals < 1e-6)){
                ccat("\n")
                varnames = colnames(linear_mat_noIntercept)
                collin_var = varnames[max_residuals < 1e-6]
                message = paste0("Variable", enumerate_items(collin_var, "s.is.quote"), " collinear with fixed-effects '", names(cluster)[q], "'.")

                print(message)
                return(invisible(message))

            }

        }
        ccat("OK")
    }

    if(isSlope && isLinear){
        ccat("simple with variables with varying slopes:")

        for(q in 1:Q_slope){
            ccat(".")
            dum = fe_id[[q]]
            my_var = slope_var_all[[q]]
            k = max(dum)
            value = cpp_tapply_sum(Q = k, x = linear_mat_noIntercept*my_var, dum = dum)

            denom = as.vector(tapply(my_var**2, dum, sum))

            # residuals of the linear projection on the cluster space
            residuals = linear_mat_noIntercept - (value/denom)[dum, ]*my_var

            max_residuals = apply(abs(residuals), 2, max)

            if(any(max_residuals < 1e-6)){
                ccat("\n")
                varnames = colnames(linear_mat_noIntercept)
                collin_var = varnames[max_residuals < 1e-6]

                message = paste0("Variable", enumerate_items(collin_var, "s.is.quote"), " collinear with variable with varying slope '", slope_vars[q], "' (on '", slope_fe[q], "').")

                print(message)
                return(invisible(message))

            }

        }
        ccat("OK")
    }

    #
    # II) perfect multicollinearity
    #

    name2change = grepl("[^[:alnum:]\\._]", colnames(linear.matrix))
    dict_name = colnames(linear.matrix)[!grepl("(Intercept)", colnames(linear.matrix))]
    if(any(name2change)){
        linbase = as.data.frame(linear.matrix)
        names(linbase)[name2change] = gsub("[^[:alnum:]\\._]", "_", names(linbase)[name2change])
    } else {
        linbase = as.data.frame(linear.matrix)
    }
    linearVars = names(linbase)[!grepl("Intercept", names(linbase))]
    names(dict_name) = linearVars
    dict_name["_Intercept_"] = "(Intercept)"

    # for multicol without cluster
    mat_base = as.matrix(linbase)

    # we add possibly missing variables
    varmiss = setdiff(all.vars(linear_fml), names(linbase))
    for(v in varmiss) linbase[[v]] = data[[v]]

    if(isLinear && length(linearVars) >= 2){
        ccat(ifelse(Q >= 1, ", ", " "), "multiple:")
        for(v in linearVars){
            ccat(".")

            i = which(colnames(mat_base) == v)
            res = ols_fit(y = mat_base[, i], X = mat_base[, -i, drop = FALSE], w = 1, collin.tol = 1e-10,  nthreads = 1)

            max_residuals = max(abs(res$residuals))

            if(max_residuals < 1e-4){
                ccat("\n")
                coef_lm = res$coefficients
                vars = colnames(mat_base)[-i]
                collin_var = vars[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
                message = paste0("Variable '", dict_name[v], "' is collinear with variable", enumerate_items(dict_name[collin_var], "s.quote"), ".")

                print(message)
                return(invisible(message))
            }
        }
        ccat("OK")
    }

    #
    # II.b) perfect multicollinearity + fixed-effects
    #

    # linearVars = setdiff(colnames(linear.matrix), "(Intercept)")
    if(isLinear && length(linearVars) >= 2 && isFixef){
        ccat(", multiple with cluster")

        dum_names = names(x$fixef_id)
        n_clust = length(dum_names)
        new_dum_names = paste0("__DUM_", 1:n_clust)
        for(i in 1:n_clust){
            linbase[[paste0("__DUM_", i)]] = x$fixef_id[[i]]
        }

        for(v in linearVars){
            ccat(".")
            fml2estimate = as.formula(paste0(v, "~", paste0(setdiff(linearVars, v), collapse = "+")))

            for(id_cluster in 1:n_clust){

                res = feols(fml2estimate, linbase, fixef = new_dum_names[id_cluster], warn = FALSE)

                max_residuals = max(abs(res$residuals))

                if(max_residuals < 1e-4){
                    ccat("\n")
                    coef_lm = coef(res)
                    collin_var = names(coef_lm)[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
                    message = paste0("Variable '", dict_name[v], "' is collinear with variable", enumerate_items(dict_name[collin_var], "s.quote"), ", together with the fixed-effects ", dum_names[id_cluster], ".")

                    print(message)
                    return(invisible(message))
                }
            }
        }
        ccat("OK")
    }

    #
    # III) NL problem
    #

    if(isNL){
        NL_vars = all.vars(NL_fml)
        varNotHere = setdiff(NL_vars, c(names(coef), names(data)))
        if(length(varNotHere) > 0){
            stop("Some variables used to estimate the model (in the non-linear formula) are missing from the original data: ", paste0(varNotHere, collapse = ", "), ".")
        }

        var2send = intersect(NL_vars, names(data))
        env = new.env()
        for(var in var2send){
            assign(var, data[[var]], env)
        }
    }

    if(isNL && (length(coef) >= 2 || isFixef)){
        ccat(", in non-linear term:")
        if(isFixef){
            # we add the constant
            coef["CONSTANT"] = 1
            data$CONSTANT = 1
            linear.matrix = model.matrix(update(linear_fml, ~.-1+CONSTANT), data)
        }

        coef_name = names(coef)

        NL_coef = setdiff(all.vars(NL_fml), names(data))
        # L_coef = setdiff(coef_name, NL_coef)
        L_coef = colnames(linear.matrix)

        #
        # We compute mu:
        #

        mu = 0
        if(length(L_coef) > 0){
            mu = linear.matrix %*% coef[L_coef]
        }

        for(iter_coef in NL_coef){
            assign(iter_coef, coef[iter_coef], env)
        }

        # Evaluation of the NL part
        value_NL = eval(NL_fml[[2]], env)
        mu = mu + value_NL
        data$mu = mu

        if(var(mu) == 0){
            message = "Variance of the NL part is 0."
            print(message)
            return(invisible(message))
        }

        #
        # The loop
        #

        for(var in coef_name){
            ccat(".")
            # we modify the coef of var, then we fix it:
            # if we can find a set of other parameters that give the same fit,
            # then there is an identification issue

            if(var %in% L_coef){
                # if linear coef: we modify the fml and add an offset
                if(length(L_coef) == 1){
                    fml = mu ~ 0
                } else {
                    fml = as.formula(paste0("mu ~ 0+", paste0(setdiff(L_coef, var), collapse = "+")))
                }

                # offset with the new value
                offset = as.formula(paste0("~", coef[var] + 1, "*", var))

                res = feNmlm(fml, data = data, family = "gaussian", NL.fml = NL_fml, NL.start = as.list(coef), offset = offset)

            } else {
                # NL case:
                # we modify both fml and NL.fml

                if(isFixef){
                    fml = update(x$fml, mu ~ . - 1 + CONSTANT)
                } else {
                    fml = update(x$fml, mu ~ .)
                }

                # the new formula with the new variable
                NL_fml_char = as.character(NL_fml)[2]
                NL_fml_new = as.formula(paste0("~", gsub(paste0("\\b", var, "\\b"), 1 + coef[var], NL_fml_char)))

                if(length(NL_coef) == 1){
                    # there is no more parameter to estimate => offset
                    res = femlm(fml, data = data, family = "gaussian", offset = NL_fml_new)
                } else {
                    res = feNmlm(fml, data = data, family = "gaussian", NL.fml = NL_fml_new, NL.start = as.list(coef))
                }

            }

            sum_resids = sum(abs(res$residuals))
            if(sum_resids < 1e-4){
                coef_esti = coef(res)
                coef_diff = abs(coef_esti - coef[names(coef_esti)])
                collin_var = names(coef_diff)[coef_diff > 1e-3]
                message = paste0("Coefficients ", show_vars_limited_width(c(var, collin_var)), " are not uniquely identifed.")

                print(message)
                return(invisible(message))
            }

        }
        ccat("OK")
    }

    ccat("\n")
    message = "No visible collinearity problem. (Doesn't mean there's none!)"
    print(message)
    return(invisible(message))

}



#' Treated and control sample descriptives
#'
#' This function shows the means and standard-deviations of several variables conditional on whether they are from the treated or the control group. The groups can further be split according to a pre/post variable. Results can be seamlessly be exported to Latex.
#'
#'
#' @param fml Either a formula of the type \code{var1 + ... + varN ~ treat} or \code{var1 + ... + varN ~ treat | post}. Either a data.frame/matrix containing all the variables for which the means are to be computed (they must be numeric of course). Both the treatment and the post variables must contain only exactly two values. You can use a point to select all the variables of the data set: \code{. ~ treat}.
#' @param base A data base containing all the variables in the formula \code{fml}.
#' @param treat_var Only if argument \code{fml} is *not* a formula. The vector identifying the treated and the control observations (the vector can be of any type but must contain only two possible values). Must be of the same length as the data.
#' @param post_var Only if argument \code{fml} is *not* a formula. The vector identifying the periods (pre/post) of the observations (the vector can be of any type but must contain only two possible values). The first value (in the sorted sense) of the vector is taken as the pre period. Must be of the same length as the data.
#' @param treat_first Which value of the 'treatment' vector should appear on the left? By default the max value appears first (e.g. if the treatment variable is a 0/1 vector, 1 appears first).
#' @param tex Should the result be displayed in Latex? Default is \code{FALSE}. Automatically set to \code{TRUE} if the table is to be saved in a file using the argument \code{file}.
#' @param treat_dict A character vector of length two. What are the names of the treated and the control? This should be a dictionary: e.g. \code{c("1"="Treated", "0" = "Control")}.
#' @param dict A named character vector. A dictionary between the variables names and an alias. For instance \code{dict=c("x"="Inflation Rate")} would replace the variable name \code{x} by \dQuote{Inflation Rate}.
#' @param file A file path. If given, the table is written in Latex into this file.
#' @param replace Default is \code{TRUE}, which means that when the table is exported, the existing file is not erased.
#' @param title Character string giving the Latex title of the table. (Only if exported.)
#' @param label Character string giving the Latex label of the table. (Only if exported.)
#' @param raw Logical, default is \code{FALSE}. If \code{TRUE}, it returns the information without formatting.
#' @param indiv Either the variable name of individual identifiers, a one sided formula, or a vector. If the data is that of a panel, this can be used to track the number of individuals per group.
#' @param prepostnames Only if there is a 'post' variable. The names of the pre and post periods to be displayed in Latex. Default is \code{c("Before", "After")}.
#' @param diff.inv Logical, default to \code{FALSE}. Whether to inverse the difference.
#'
#' @details
#' By default, when the user tries to apply this function to nun-numeric variables, an error is raised. The exception is when the all variables are selected with the dot (like in \code{. ~ treat}. In this case, non-numeric variables are automatically omitted (with a message).
#'
#' NAs are removed automatically: if the data contains NAs an information message will be prompted. First all observations containing NAs relating to the treatment or post variables are removed. Then if there are still NAs for the variables, they are excluded separately for each variable, and a new message detailing the NA breakup is prompted.
#'
#' @return
#' It returns a data.frame or a Latex table with the conditional means and statistical differences between the groups.
#'
#' @examples
#'
#' # Playing around with the DiD data
#' data(base_did)
#'
#' # means of treat/control
#' did_means(y+x1+period~treat, base_did)
#'
#' # same but inverting the difference
#' did_means(y+x1+period~treat, base_did, diff.inv = TRUE)
#'
#' # now treat/control, before/after
#' did_means(y+x1+period~treat|post, base_did)
#'
#' # same but with a new line giving the number of unique "indiv" for each case
#' did_means(y+x1+period~treat|post, base_did, indiv = "id")
#'
#' # same but with the treat case "0" coming first
#' did_means(y+x1+period~treat|post, base_did, indiv = ~id, treat_first = 0)
#'
#' # Selecting all the variables with "."
#' did_means(.~treat|post, base_did, indiv = "id")
#'
#'
did_means = function(fml, base, treat_var, post_var, tex = FALSE, treat_dict,
                     dict = getFixest_dict(), file, replace = FALSE, title,
                     label, raw = FALSE, indiv, treat_first, prepostnames = c("Before", "After"),
                     diff.inv = FALSE){
    # x is a data.frame
    # treat_vector is a list of IDs
    # treat_first: the treated-value to appear first

    # Ajouter protection quand un des groupes vaut 0!!!
    # ie un des treat-post n'a pas d'observation

    fml_in = fml

    #
    # Nber of units of observations (I do it before the main data)
    #

    IS_INDIV = FALSE
    if(!missing(indiv)){
        IS_INDIV = TRUE
        if("formula" %in% class(indiv)){

            check_arg(indiv, "os formula")
            indiv_varname = attr(terms(indiv), "term.labels")

            if(missing(base)){
                stop("To use 'indiv' as a formula, you must provide the argument 'base'.")
            }

            if(length(indiv_varname) > 1){
                stop("The argument 'indiv' must refer to only one variable.")
            }

            if(!all(all.vars(indiv) %in% names(base))){
                pblm = setdiff(all.vars(indiv), names(base))
                stop("In argument 'indiv': the variable", enumerate_items(pblm, "s.is")," not in the data set.")
            }

            indiv_var = try(eval(str2lang(indiv_varname), base), silent = TRUE)
            if("try-error" %in% class(indiv_var)){
                stop("Evaluation of 'indiv' raises and error:\n", indiv_var)
            }
        } else if(length(indiv) == 1 && is.character(indiv)){
            indiv_varname = indiv

            if(missing(base) || !indiv %in% names(base)){
                stop("To use 'indiv' as a 'character string' you must provide the argument 'base'.")
            } else {
                indiv_var = base[[indiv]]
            }
        } else {
            if(!missing(base) && length(indiv) != NROW(base)){
                stop("The length of 'indiv' must be the same as the data.")
            } else if(missing(base) && length(indiv) != NROW(fml)){
                stop("The length of 'indiv' must be the same as the data.")
            }
            indiv_var = indiv
        }
    } else {
        indiv_var = NULL
    }

    #
    # Extracting the information
    #

    usePost = FALSE
    if("formula" %in% class(fml_in)){
        if(missing(base) || !is.data.frame(base)){
            stop("If you provide a formula, a data.frame must be given in argument 'base'.")
        }

        # Creation of x and the condition
        if(!length(fml_in) == 3){
            stop("The formula must be of the type 'var1 + ... + varN ~ treat' or 'var1 + ... + varN ~ treat | post'.")
        }

        fml_parts = fml_split(fml_in)

        fml = fml_parts[[1]]
        pipe = NULL
        if(length(fml_parts) > 1){
            pipe = fml_parts[[2]]
        }

        #
        # Treat & post
        #

        if(!(all(all.vars(fml[[3]]) %in% names(base)))){
            pblm = setdiff(all.vars(fml[[3]]), names(base))
            stop("In the evaluation of the treatment variable: ", enumerate_items(pblm, "is"), " not in the data set.")
        }

        treat_var = try(eval(fml[[3]], base), silent = TRUE)
        if("try-error" %in% class(treat_var)){
            stop("Evaluation of the 'treatment' variable raises and error: \n", treat_var)
        }

        if(!is.null(pipe)){
            if(!(all(all.vars(pipe) %in% names(base)))){
                pblm = setdiff(all.vars(pipe), names(base))
                stop("In the evaluation of the 'post' variable: ", enumerate_items(pblm, "is"), " not in the data set.")
            }

            post_var = try(eval(pipe[[2]], base), silent = TRUE)
            if("try-error" %in% class(post_var)){
                stop("Evaluation of the 'post' variable raises and error: \n", treat_var)
            }

            usePost = TRUE
        }

        #
        # The variables
        #

        # special case: the point
        if(deparse(fml[[2]])[1] == "."){
            all_vars = setdiff(names(base), deparse(fml[[3]])[1])
            if(usePost) all_vars = setdiff(all_vars, as.character(pipe))
            if(IS_INDIV) all_vars = setdiff(all_vars, indiv_varname)

            if("data.table" %in% class(base)){
                mat_vars = as.data.frame(base[, all_vars, with = FALSE])
            } else {
                mat_vars = base[, all_vars, FALSE]
            }

            # we drop non-numeric info + note
            base_small = head(mat_vars, 10)
            is_num = sapply(base_small, function(x) is.numeric(x) || is.logical(x))
            pblm = all_vars[!is_num]

            if(all(!is_num)){
                stop("Function cannot be performed: not any numeric variable.")
            } else if(length(pblm) > 0){
                message("NOTE: The variable", enumerate_items(pblm, "s.is"), " removed because they are non-numeric.")
            }

            mat_vars = as.matrix(mat_vars[, is_num, FALSE])
        } else {
            # we swap the formula to use model.frame
            fml_x = as.formula(paste0("1~", deparse_long(fml[[2]]), "-1"))

            # Variables control:
            all_vars = all.vars(fml_x)
            pblm = setdiff(all_vars, names(base))
            if(length(pblm) > 0){
                stop("The variable", enumerate_items(pblm, "s.is"), " not in the data set (", deparse(substitute(base)), ").")
            }

            # Evaluation control
            base_small = head(base, 10)
            var2eval = attr(terms(fml_x), "term.labels")
            var2eval = gsub(":", "*", var2eval)
            for(i in seq_along(var2eval)){
                var = var2eval[i]
                x_small = try(eval(parse(text=var), base_small), silent = TRUE)
                if("try-error" %in% class(x_small)){
                    stop("Evaluation of the variable '", var, "' raises and error:\n", x_small)
                }

                if(!is.numeric(x_small) && !is.logical(x_small)){
                    stop("The variable ", var, " is not numeric, please provide only numeric variables.")
                }
            }

            mat_vars = prepare_matrix(fml_x, base)
        }

        if(is.logical(mat_vars)){
            mat_vars = as.numeric(mat_vars)
        }

        # other info
        x_name = colnames(mat_vars)
    } else {
        if(missing(treat_var)){
            stop("If argument 'fml' is not a formula, you must provide the argument 'treat_var'.")
        } else {
            mat_vars = fml_in

            if(NROW(mat_vars) != length(treat_var)){
                stop("The arguments 'x' and 'treat_var' must be of the same length.")
            }

            if(!missing(post_var) && !is.null(post_var)){
                if(NROW(mat_vars) != length(post_var)){
                    stop("If provided, the argument 'post_var' must be of the same length of 'x'.")
                }
                usePost = TRUE
            }

            # We exclude non numeric variables
            if("data.frame" %in% class(mat_vars)){
                is_num = sapply(mat_vars, function(x) is.numeric(x) || is.logical(x))
                if(any(!is_num)){
                    pblm = names(mat_vars)[!is_num]
                    stop("This function works only for numeric variables. Variable", enumerate_items(pblm, "s.is"), " not numeric.")
                }
                mat_vars = as.matrix(mat_vars)
            } else if(!is.matrix(mat_vars)){
                stop("If not a formula, argument 'fml' must be a data.frame or a matrix. Currently it is of class ", class(mat_vars)[1], ".")
            } else if(!is.numeric(mat_vars) && !is.logical(mat_vars)){
                stop("If not a formula, argument 'fml' must be a data.frame or a matrix with numeric variables. Currenlty its is not numeric.")
            }

            if(is.logical(mat_vars)) mat_vars = as.numeric(mat_vars)

            # other info
            x_name = colnames(mat_vars)
        }
    }

    #
    # NA control
    #

    # Treat & post
    ANY_NA = FALSE
    if(anyNA(treat_var) || (usePost && anyNA(post_var))){
        ANY_NA = TRUE
        qui_NA = is.na(treat_var)
        if(usePost) qui_NA = qui_NA | is.na(post_var)

        if(all(qui_NA)){
            msg = "treatment variable."
            if(usePost && anyNA(post_var)){
                if(anyNA(treat_var)){
                    msg = "'treatment' and 'post' variables."
                } else {
                    msg = "'post' variable."
                }
            }
            stop("All observation contain NA values for the ", msg)
        }

        mat_vars = mat_vars[!qui_NA, ]
        treat_var = treat_var[!qui_NA]
        if(usePost) post_var = post_var[!qui_NA]

        if(IS_INDIV) indiv_var = indiv_var[!qui_NA]

        msg = ifelse(usePost, "/post variables.", " variable.")
        message("NOTE: ", sum(qui_NA), " observations removed because of NA values in the treatment", msg)

    }



    # Treatment / post controls
    if(length(unique(treat_var)) != 2){
        n = length(unique(treat_var))
        msg = ifelse(n == 1, "only one value.", paste0(n, " values."))
        stop("This function supports only 2 conditional values for the treatment variable. Currently, it contains ", msg)
    }

    if(usePost && length(unique(post_var)) != 2){
        n = length(unique(post_var))
        msg = ifelse(n == 1, "only one value.", paste0(n, " values."))
        stop("This function supports only 2 conditional values for the 'post' variable. Currently, it contains ", msg)
    }

    if(usePost){
        check_arg(prepostnames, "character vector len(2) no na")
    }

    treat_cases = sort(unique(treat_var), decreasing = TRUE)
    if(!missing(treat_dict) && !is.null(treat_dict)){

        if(!isVector(treat_dict) || is.null(names(treat_dict))){
            stop("The argument 'treat_dict' must be a named character vector.")
        }

        pblm = setdiff(treat_cases, names(treat_dict))
        if(length(pblm) > 0){
            stop("The value", enumerate_items(pblm, "s.is"), " not in 'treat_dict'.")
        }
    } else {
        treat_dict = paste0("cond: ", treat_cases)
        names(treat_dict) = treat_cases
    }

    #
    # The main function
    #

    if(!missing(treat_first) && !is.null(treat_first)){
        if(!treat_first %in% treat_cases){
            stop("Argument 'treat_first' must be an element of the treated variable.")
        } else if(treat_first != treat_cases[1]){
            treat_cases = rev(treat_cases)
        }
    }

    if(!usePost){
        # the simple case
        res = .meanDiffTable(mat_vars, treat_var, treat_cases, indiv_var, raw = raw, diff.inv = diff.inv)

        if(any(attr(res, "na"))){
            nb_na = attr(res, "na")
            var_with_na = x_name[nb_na > 0]
            message("NA values were omitted for the variable", enumerate_items(paste0(var_with_na, " (", nb_na[nb_na > 0], ")"), "s"), ".")
        }

        attr(res, "na") = NULL


    } else {
        # the more complex case
        post_values = sort(unique(post_var))

        qui_pre = which(post_var == post_values[1])
        res_pre = .meanDiffTable(mat_vars[qui_pre, ], treat_var[qui_pre], treat_cases, indiv_var[qui_pre], raw = raw, diff.inv = diff.inv)

        qui_post = which(post_var == post_values[2])
        res_post = .meanDiffTable(mat_vars[qui_post, ], treat_var[qui_post], treat_cases, indiv_var[qui_post], raw = raw, diff.inv = diff.inv)

        nb_na = attr(res_pre, "na") + attr(res_post, "na")

        res = cbind(res_pre, res_post[, -1])

        if(any(nb_na)){
            var_with_na = x_name[nb_na > 0]
            message("NA values were omitted for the variable", enumerate_items(paste0(var_with_na, " (", nb_na[nb_na > 0], ")"), "s"), ".")
        }

        attr(res, "na") = NULL
    }

    #
    # Tex exportation
    #

    if(tex || !missing(file)){
        # we sink!
        if(!missing(file)){
            sink(file = file, append = (replace == FALSE))
            on.exit(sink())
        }

        if(!missing(title)){
            # if there is the title, then we create a "table" environment
            myTitle = title
            if(!missing(label)) myTitle = paste0("\\label{", label, "} ", myTitle)
            cat("\\begin{table}[htbp]\\centering\n\\caption{",  myTitle, "}\n")
        }

        #
        # we write the tex code
        #

        treat_text = paste0(treat_dict[as.character(treat_cases)], collapse = " & ")

        if(!usePost){
            # The simple case
            cat("\\begin{tabular}{lcccc}\n")
            cat(" &  &  &  & \\tabularnewline\n\\hline\n\\hline\n")

            cat(" & ", treat_text,  "\\tabularnewline\n")

            # the variables "line"
            cat("Variables & mean (s.e.) & mean (s.e.) & Diff. & t-stat \\tabularnewline\n")
            cat("\\hline\n")
        } else {
            # The case with pre and post
            cat("\\begin{tabular}{lcccccccc}\n")
            cat(" &  &  &  &  &  &  &  & \\tabularnewline\n\\hline\n\\hline\n")

            # the treat post line
            cat(" & ", "\\multicolumn{4}{c}{", prepostnames[1], "} & \\multicolumn{4}{c}{", prepostnames[2], "}\\tabularnewline\n")

            # the 'treated' line
            cat(" & ", treat_text, " &  &  & ", treat_text, "\\tabularnewline\n")

            # the variables "line"
            cat("Variables & mean (s.e.) & mean (s.e.) & Diff. & t-stat & mean (s.e.) & mean (s.e.) & Diff. & t-stat \\tabularnewline\n")
            cat("\\hline\n")
        }


        # Changing the names of the variables
        if(!missnull(dict)){
            qui = which(res$vars %in% names(dict))
            for(i in qui) res$vars[i] = escape_latex(dict[res$vars[i]])
        }

        # The data
        for(i in 1:nrow(res)){
            cat(paste0(res[i,], collapse = " & "), "\\tabularnewline\n")
        }

        cat("\\hline\n\\hline\n &  &  &  & \\tabularnewline\n")
        cat("\\end{tabular}\n")

        if(!missing(title)) cat("\\end{table}\n")

        return(invisible(res))
    }

    res
}

.meanDiffTable = function(mat_vars, treat_var, treat_cases, indiv_var, raw, diff.inv = FALSE){
    # This function is the workhorse of did_means

    x_name = colnames(mat_vars)

    treat_01 = 1 * (treat_var == treat_cases[2])

    nthreads = getFixest_nthreads()
    info = cpppar_cond_means(mat_vars, treat_01, nthreads)

    means = info$means
    sds = info$sd
    ns = info$n
    n_01 = info$n_01

    name_1 = as.character(treat_cases[1])
    name_2 = as.character(treat_cases[2])

    means_1 = means[, 1]
    means_2 = means[, 2]

    sds_1 = sds[, 1]
    sds_2 = sds[, 2]

    n_1 = ns[, 1]
    n_2 = ns[, 2]

    format_1 = paste0(mysignif(means_1, 2), " (", mysignif(sds_1, 2), ")")
    format_2 = paste0(mysignif(means_2, 2), " (", mysignif(sds_2, 2), ")")

    if(diff.inv){
        D = (means_2 - means_1)
    } else {
        D = (means_1 - means_2)
    }

    sd_D = sqrt(sds_1**2/n_1 + sds_2**2/n_2)

    t_stat = signif(D/sd_D, 3)

    #
    # Formatting
    #

    if(raw){
        res = data.frame(vars = x_name, stringsAsFactors = FALSE)
        res[[paste0("Mean: ", treat_cases[1])]] = c(means_1, recursive = TRUE)
        res[[paste0("Mean: ", treat_cases[2])]] = c(means_2, recursive = TRUE)
        res[[paste0("se: ", treat_cases[1])]] = c(sds_1, recursive = TRUE)
        res[[paste0("se: ", treat_cases[2])]] = c(sds_2, recursive = TRUE)
        res[["Difference"]] = c(D, recursive=TRUE)
        res[["t-stat"]] = c(as.character(t_stat), recursive=TRUE)
    } else {
        res = data.frame(vars = x_name, stringsAsFactors = FALSE)
        res[[paste0("cond: ", treat_cases[1])]] = format_1
        res[[paste0("cond: ", treat_cases[2])]] = format_2
        res[["Difference"]] = c(mysignif(D, 3), recursive = TRUE)
        res[["t-stat"]] = c(as.character(t_stat), recursive = TRUE)

        res[nrow(res) + 1, ] = c("Observations", addCommas(n_01[1]), addCommas(n_01[2]), "", "")
        # Le nombre d'individus
        if(!missing(indiv_var) && !is.null(indiv_var)){
            nb_indiv =  tapply(indiv_var, treat_var, function(x) length(unique(x)))
            res[nrow(res) + 1, ] = c("# Individuals", addCommas(nb_indiv[name_1]), addCommas(nb_indiv[name_2]), "", "")
        }
    }

    attr(res, "na") = info$na

    res
}



#' Create, or interact variables with, factors
#'
#' Treat a variable as a factor, or interacts a variable with a factor. Values to be dropped/kept from the factor can be easily set. Note that to interact fixed-effects, this function should not be used: instead use directly the syntax \code{fe1^fe2}.
#'
#'
#' @inheritParams bin
#'
#' @param factor_var  A vector (of any type) that will be treated as a factor. You can set references (i.e. exclude values for which to create dummies) with the \code{ref} argument.
#' @param var A variable of the same length as \code{factor_var}. This variable will be interacted with the factor in \code{factor_var}. It can be numeric or factor-like. To force a numeric variable to be treated as a factor, you can add the \code{i.} prefix to a variable name. For instance take a numeric variable \code{x_num}: \code{i(x_fact, x_num)} will treat \code{x_num} as numeric while \code{i(x_fact, i.x_num)} will treat \code{x_num} as a factor (it's a shortcut to \code{as.factor(x_num)}).
#' @param ref A vector of values to be taken as references from \code{factor_var}. Can also be a logical: if \code{TRUE}, then the first value of \code{factor_var} will be removed. If \code{ref} is a character vector, partial matching is applied to values; use "@" as the first character to enable regular expression matching. See examples.
#' @param keep A vector of values to be kept from \code{factor_var} (all others are dropped). By default they should be values from \code{factor_var} and if \code{keep} is a character vector partial matching is applied. Use "@" as the first character to enable regular expression matching instead.
#' @param ref2 A vector of values to be dropped from \code{var}. By default they should be values from \code{var} and if \code{ref2} is a character vector partial matching is applied. Use "@" as the first character to enable regular expression matching instead.
#' @param keep2 A vector of values to be kept from \code{var} (all others are dropped). By default they should be values from \code{var} and if \code{keep2} is a character vector partial matching is applied. Use "@" as the first character to enable regular expression matching instead.
#' @param bin2 A list or vector defining the binning of the second variable. See help for the argument \code{bin} for details (or look at the help of the function \code{\link[fixest]{bin}}). You can use \code{.()} for \code{list()}.
#' @param ... Not currently used.
#'
#' @details
#' To interact fixed-effects, this function should not be used: instead use directly the syntax \code{fe1^fe2} in the fixed-effects part of the formula. Please see the details and examples in the help page of \code{\link[fixest]{feols}}.
#'
#' @return
#' It returns a matrix with number of rows the length of \code{factor_var}. If there is no interacted variable or it is interacted with a numeric variable, the number of columns is equal to the number of cases contained in \code{factor_var} minus the reference(s). If the interacted variable is a factor, the number of columns is the number of combined cases between \code{factor_var} and \code{var}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest:coefplot]{iplot}} to plot interactions or factors created with \code{i()}, \code{\link[fixest]{feols}} for OLS estimation with multiple fixed-effects.
#'
#' See the function \code{\link[fixest]{bin}} for binning variables.
#'
#' @examples
#'
#' #
#' # Simple illustration
#' #
#'
#' x = rep(letters[1:4], 3)[1:10]
#' y = rep(1:4, c(1, 2, 3, 4))
#'
#' # interaction
#' data.frame(x, y, i(x, y, ref = TRUE))
#'
#' # without interaction
#' data.frame(x, i(x, "b"))
#'
#' # you can interact factors too
#' z = rep(c("e", "f", "g"), c(5, 3, 2))
#' data.frame(x, z, i(x, z))
#'
#' # to force a numeric variable to be treated as a factor: use i.
#' data.frame(x, y, i(x, i.y))
#'
#' # Binning
#' data.frame(x, i(x, bin = list(ab = c("a", "b"))))
#'
#' # Same as before but using .() for list() and a regular expression
#' # note that to trigger a regex, you need to use an @ first
#' data.frame(x, i(x, bin = .(ab = "@a|b")))
#'
#' #
#' # In fixest estimations
#' #
#'
#' data(base_did)
#' # We interact the variable 'period' with the variable 'treat'
#' est_did = feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)
#'
#' # => plot only interactions with iplot
#' iplot(est_did)
#'
#' # Using i() for factors
#' est_bis = feols(y ~ x1 + i(period, keep = 3:6) + i(period, treat, 5) | id, base_did)
#'
#' # we plot the second set of variables created with i()
#' # => we need to use keep (otherwise only the first one is represented)
#' coefplot(est_bis, keep = "trea")
#'
#' # => special treatment in etable
#' etable(est_bis, dict = c("6" = "six"))
#'
#' #
#' # Interact two factors
#' #
#'
#' # We use the i. prefix to consider week as a factor
#' data(airquality)
#' aq = airquality
#' aq$week = aq$Day %/% 7 + 1
#'
#' # Interacting Month and week:
#' res_2F = feols(Ozone ~ Solar.R + i(Month, i.week), aq)
#'
#' # Same but dropping the 5th Month and 1st week
#' res_2F_bis = feols(Ozone ~ Solar.R + i(Month, i.week, ref = 5, ref2 = 1), aq)
#'
#' etable(res_2F, res_2F_bis)
#'
#' #
#' # Binning
#' #
#'
#' data(airquality)
#'
#' feols(Ozone ~ i(Month, bin = "bin::2"), airquality)
#'
#' feols(Ozone ~ i(Month, bin = list(summer = 7:9)), airquality)
#'
#'
#'
i = function(factor_var, var, ref, keep, bin, ref2, keep2, bin2, ...){
    # Used to create interactions

    # Later: binning (bin = 1:3 // bin = list("a" = "[abc]")). Default name is bin name (eg "1:3")

    # gt = function(x) cat(sfill(x, 20), ": ", -(t0 - (t0<<-proc.time()))[3], "s\n", sep = "")
    # t0 = proc.time()

    validate_dots(valid_args = c("f2", "f_name", "ref_special", "sparse"))

    dots = list(...)
    is_sparse = isTRUE(dots$sparse)

    mc = match.call()

    # Finding if it's a call from fixest
    FROM_FIXEST = is_fixest_call()

    # General checks
    check_arg(factor_var, "mbt vector")

    # NOTA:
    # the user can use the prefix "i." to tell the algorithm to consider the
    # variable var as a factor. This requires a non standard evaluation
    #

    var_name = NULL
    if(!is.null(mc$f2)){
        # should be only used internally

        var_name = deparse_long(mc$f2)
        var = dots$f2
        check_value(var, "vector", .arg_name = "f2")
    }

    # General information on the factor variable
    if(!is.null(dots$f_name)){
        f_name = dots$f_name
    } else {
        f_name = deparse_long(mc$factor_var)
    }
    f = factor_var # renaming for clarity

    # checking var
    IS_INTER_NUMERIC = IS_INTER_FACTOR = FALSE
    if(!missing(var)){
        # Evaluation of var (with possibly, user prepending with i.)
        if(is.null(var_name)){
            var_name = deparse_long(mc$var)
            if(grepl("^i\\.", var_name)){
                var_name = gsub("^i\\.", "", var_name)
                var = str2lang(var_name)
                check_value_plus(var, "evalset vector", .data = parent.frame(), .arg_name = "var")
                IS_INTER_FACTOR = TRUE
            } else {
                check_value(var, "vector", .arg_name = "var")
            }
        } else {
            # f2 was used
            IS_INTER_FACTOR = TRUE
        }

        # is it an interaction with a factor or with a continuous variable?
        if(length(var) == length(f)){
            # Logical => numeric by default
            # otherwise, just setting the flags

            if(IS_INTER_FACTOR == FALSE){
                # If previous condition == TRUE, means was requested by the user

                if(is.logical(var)){
                    var = as.numeric(var)
                    IS_INTER_NUMERIC = TRUE
                } else if(is.numeric(var)){
                    IS_INTER_NUMERIC = TRUE
                } else {
                    IS_INTER_FACTOR = TRUE
                }
            }

        } else if(length(var) <= 2 && !"ref" %in% names(mc)){
            # ref is implicitly called via the location of var
            ref = var
            dots$ref_special = FALSE
        } else {

            if(grepl("^[[:alpha:]\\.][[:alnum:]\\._]*:[[:alpha:]\\.][[:alnum:]\\._]*$", var_name)){
                info = strsplit(var_name, ":")[[1]]
                stop("In i(): When 'var' is equal to a product, please use I(", info[1], "*", info[2], ") instead of ", var_name, ".")
            } else {
                stop("The arguments 'var' and 'f' must be of the same length (currently ", length(var), " vs ", length(f), ").")
            }
        }
    }

    IS_INTER = IS_INTER_NUMERIC || IS_INTER_FACTOR

    if(isTRUE(dots$ref_special) && (!IS_INTER || IS_INTER_FACTOR)){
        ref = TRUE
    }

    #
    # QUFing + NA
    #

    if(IS_INTER){
        is_na_all = is.na(f) | is.na(var)
    } else {
        is_na_all = is.na(f)
    }

    if(!missing(bin)){
        bin = error_sender(eval_dot(bin), arg_name = "bin")

        if(!is.null(bin)){
            f = bin_factor(bin, f, f_name)
        }
    }

    if(IS_INTER_FACTOR && !MISSNULL(bin2)){
        bin2 = error_sender(eval_dot(bin2), arg_name = "bin2")

        if(!is.null(bin2)){
            var = bin_factor(bin2, var, f_name)
        }
    }

    if(IS_INTER_FACTOR){
        info = to_integer(f, var, add_items = TRUE, items.list = TRUE, sorted = TRUE, multi.join = "__%%__")
    } else {
        info = to_integer(f, add_items = TRUE, items.list = TRUE, sorted = TRUE)
    }

    fe_num = info$x
    items = info$items

    if(!IS_INTER_NUMERIC){
        # neutral var in C code
        var = 1
    }

    if(IS_INTER_FACTOR){
        items_split = strsplit(items, "__%%__", fixed = TRUE)

        f_items = sapply(items_split, `[`, 1)
        var_items = sapply(items_split, `[`, 2)
    } else {
        f_items = items
    }

    check_arg(ref, "logical scalar | vector no na")

    check_arg(ref2, keep, keep2, "vector no na")

    no_rm = TRUE
    id_drop = c()
    if(!missing(ref)){
        if(isTRUE(ref)){
            # We always delete the first value
            # Que ce soit items ici est normal (et pas f_items)
            id_drop = which(items == items[1])
        } else {
            id_drop = items_to_drop(f_items, ref, "factor_var")
        }
        ref_id = id_drop
    }


    if(!missing(keep)){
        id_drop = c(id_drop, items_to_drop(f_items, keep, "factor_var", keep = TRUE))
    }

    if(IS_INTER_FACTOR){
        if(!missing(ref2)){
            id_drop = c(id_drop, items_to_drop(var_items, ref2, "var"))
        }

        if(!missing(keep2)){
            id_drop = c(id_drop, items_to_drop(var_items, keep2, "var", keep = TRUE))
        }
    }

    if(length(id_drop) > 0){
        id_drop = unique(sort(id_drop))
        if(length(id_drop) == length(items)) stop("All items from the interaction have been removed.")
        who_is_dropped = id_drop
        no_rm = FALSE
    } else {
        # -1 is neutral
        who_is_dropped = -1
    }

    # The column names

    if(length(id_drop) > 0){
        items_name = items[-id_drop]
        f_items = f_items[-id_drop]
        if(IS_INTER_FACTOR){
            var_items = var_items[-id_drop]
        }
    } else {
        items_name = items
    }

    if(FROM_FIXEST){
        # Pour avoir des jolis noms c'est un vrai gloubiboulga,
        # mais j'ai pas trouve plus simple...
        if(IS_INTER_FACTOR){
            col_names = paste0("__CLEAN__", f_name, "::", f_items, ":", var_name, "::", var_items)
            col_names_full = NULL

        } else if(IS_INTER_NUMERIC){
            col_names = paste0("__CLEAN__", f_name, "::", f_items, ":", var_name)
            col_names_full = paste0(f_name, "::", items, ":", var_name)

        } else {
            col_names = paste0("__CLEAN__", f_name, "::", items_name)
            col_names_full = paste0(f_name, "::", items)
        }
    } else {

        if(IS_INTER_FACTOR){
            items_name = gsub("__%%__", ":", items_name, fixed = TRUE)
        }

        col_names = items_name
    }

    if(is_sparse){
        # Internal call: we return the row ids + the values + the indexes + the names

        if(length(who_is_dropped) > 0){
            valid_row = !is_na_all & !fe_num %in% who_is_dropped
        } else {
            valid_row = !is_na_all
        }

        # we need to ensure the IDs go from 1 to # Unique
        fe_colid = to_integer(fe_num[valid_row], sorted = TRUE)

        values = if(length(var) == 1) rep(1, length(valid_row)) else var
        res = list(rowid = which(valid_row), values = values,
                   colid = fe_colid, col_names = col_names)

        class(res) = "i_sparse"

        return(res)
    }

    res = cpp_factor_matrix(fe_num, is_na_all, who_is_dropped, var, col_names)
    # res => matrix with...
    #  - NAs where appropriate
    #  - appropriate number of columns
    #  - interacted if needed
    #


    # We send the information on the reference
    if(FROM_FIXEST){
        is_GLOBAL = FALSE
        for(where in 1:min(8, sys.nframe())){
            if(exists("GLOBAL_fixest_mm_info", parent.frame(where))){
                GLOBAL_fixest_mm_info = get("GLOBAL_fixest_mm_info", parent.frame(where))
                is_GLOBAL = TRUE
                break
            }
        }

        if(is_GLOBAL && !IS_INTER_FACTOR){
            # NOTA:
            # if you change stuff here => change them also in sunab()

            info = list()
            info$coef_names_full = col_names_full
            info$items = items
            if(!missing(ref)){
                info$ref_id = ref_id
                info$ref = items[ref_id]
            }
            info$is_num = is.numeric(items)
            info$f_name = f_name
            info$is_inter_num = IS_INTER_NUMERIC
            if(IS_INTER_NUMERIC){
                info$var_name = var_name
            }
            info$is_inter_fact = IS_INTER_FACTOR

            GLOBAL_fixest_mm_info[[length(GLOBAL_fixest_mm_info) + 1]] = info
            # re assignment
            assign("GLOBAL_fixest_mm_info", GLOBAL_fixest_mm_info, parent.frame(where))
        }

    }

    res
}


i_ref = function(factor_var, var, ref, bin, keep, ref2, keep2, bin2){
    # To automatically add references when i(x) is used

    mc = match.call()

    mc[[1]] = as.name("i")

    if(!any(c("ref", "keep", "ref2", "keep2") %in% names(mc))){
        mc$ref_special = TRUE
    }

    return(deparse_long(mc))
}

i_noref = function(factor_var, var, ref, bin, keep, ref2, keep2, bin2){
    # Used only in predict => to create data without restriction

    mc = match.call()

    mc[[1]] = as.name("i")

    mc$ref = mc$keep = mc$ref2 = mc$keep2 = NULL

    return(deparse_long(mc))
}


#' Bins the values of a variable (typically a factor)
#'
#' Tool to easily group the values of a given variable.
#'
#' @param x A vector whose values have to be grouped. Can be of any type but must be atomic.
#' @param bin A list of values to be grouped, a vector, a formula, or the special values \code{"bin::digit"} or \code{"cut::values"}. To create a new value from old values, use \code{bin = list("new_value"=old_values)} with \code{old_values} a vector of existing values. You can use \code{.()} for \code{list()}.
#' It accepts regular expressions, but they must start with an \code{"@"}, like in \code{bin="@Aug|Dec"}. It accepts one-sided formulas which must contain the variable \code{x}, e.g. \code{bin=list("<2" = ~x < 2)}.
#' The names of the list are the new names. If the new name is missing, the first value matched becomes the new name. In the name, adding \code{"@d"}, with \code{d} a digit, will relocate the value in position \code{d}: useful to change the position of factors. Use \code{"@"} as first item to make subsequent items be located first in the factor.
#' Feeding in a vector is like using a list without name and only a single element. If the vector is numeric, you can use the special value \code{"bin::digit"} to group every \code{digit} element.
#' For example if \code{x} represents years, using \code{bin="bin::2"} creates bins of two years.
#' With any data, using \code{"!bin::digit"} groups every digit consecutive values starting from the first value.
#' Using \code{"!!bin::digit"} is the same but starting from the last value.
#' With numeric vectors you can: a) use \code{"cut::n"} to cut the vector into \code{n} equal parts, b) use \code{"cut::a]b["} to create the following bins: \code{[min, a]}, \code{]a, b[}, \code{[b, max]}.
#' The latter syntax is a sequence of number/quartile (q0 to q4)/percentile (p0 to p100) followed by an open or closed square bracket. You can add custom bin names by adding them in the character vector after \code{'cut::values'}. See details and examples. Dot square bracket expansion (see \code{\link[fixest]{dsb}}) is enabled.
#'
#' @section "Cutting" a numeric vector:
#'
#' Numeric vectors can be cut easily into: a) equal parts, b) user-specified bins.
#'
#' Use \code{"cut::n"} to cut the vector into \code{n} (roughly) equal parts. Percentiles are used to partition the data, hence some data distributions can lead to create less than \code{n} parts (for example if P0 is the same as P50).
#'
#' The user can specify custom bins with the following syntax: \code{"cut::a]b]c]"etc}. Here the numbers \code{a}, \code{b}, \code{c}, etc, are a sequence of increasing numbers, each followed by an open or closed square bracket. The numbers can be specified as either plain numbers (e.g. \code{"cut::5]12[32["}), quartiles (e.g. \code{"cut::q1]q3["}), or percentiles (e.g. \code{"cut::p10]p15]p90]"}). Values of different types can be mixed: \code{"cut::5]q2[p80["} is valid provided the median (\code{q2}) is indeed greater than \code{5}, otherwise an error is thrown.
#'
#' The square bracket right of each number tells whether the numbers should be included or excluded from the current bin. For example, say \code{x} ranges from 0 to 100, then \code{"cut::5]"} will create two  bins: one from 0 to 5 and a second from 6 to 100. With \code{"cut::5["} the bins would have been 0-4 and 5-100.
#'
#' A factor is returned. The labels report the min and max values in each bin.
#'
#' To have user-specified bin labels, just add them in the character vector following \code{'cut::values'}. You don't need to provide all of them, and \code{NA} values fall back to the default label. For example, \code{bin = c("cut::4", "Q1", NA, "Q3")} will modify only the first and third label that will be displayed as \code{"Q1"} and \code{"Q3"}.
#'
#' @section \code{bin} vs \code{ref}:
#'
#' The functions \code{\link[fixest]{bin}} and \code{\link[fixest]{ref}} are able to do the same thing, then why use one instead of the other? Here are the differences:
#'
#' \itemize{
#' \item{}{\code{ref} always returns a factor. This is in contrast with \code{bin} which returns, when possible, a vector of the same type as the vector in input.}
#' \item{}{\code{ref} always places the values modified in the first place of the factor levels. On the other hand, \code{bin} tries to not modify the ordering of the levels. It is possible to make \code{bin} mimic the behavior of \code{ref} by adding an \code{"@"} as the first element of the list in the argument \code{bin}.}
#'  \item{}{when a vector (and not a list) is given in input, \code{ref} will place each element of the vector in the first place of the factor levels. The behavior of \code{bin} is totally different, \code{bin} will transform all the values in the vector into a single value in \code{x} (i.e. it's binning).}
#' }
#'
#' @return
#' It returns a vector of the same length as \code{x}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' To re-factor variables: \code{\link[fixest]{ref}}.
#'
#' @examples
#'
#' data(airquality)
#' month_num = airquality$Month
#' table(month_num)
#'
#' # Grouping the first two values
#' table(bin(month_num, 5:6))
#'
#' # ... plus changing the name to '10'
#' table(bin(month_num, list("10" = 5:6)))
#'
#' # ... and grouping 7 to 9
#' table(bin(month_num, list("g1" = 5:6, "g2" = 7:9)))
#'
#' # Grouping every two months
#' table(bin(month_num, "bin::2"))
#'
#' # ... every 2 consecutive elements
#' table(bin(month_num, "!bin::2"))
#'
#' # ... idem starting from the last one
#' table(bin(month_num, "!!bin::2"))
#'
#' # Using .() for list():
#' table(bin(month_num, .("g1" = 5:6)))
#'
#'
#' #
#' # with non numeric data
#' #
#'
#' month_lab = c("may", "june", "july", "august", "september")
#' month_fact = factor(month_num, labels = month_lab)
#'
#' # Grouping the first two elements
#' table(bin(month_fact, c("may", "jun")))
#'
#' # ... using regex
#' table(bin(month_fact, "@may|jun"))
#'
#' # ...changing the name
#' table(bin(month_fact, list("spring" = "@may|jun")))
#'
#' # Grouping every 2 consecutive months
#' table(bin(month_fact, "!bin::2"))
#'
#' # ...idem but starting from the last
#' table(bin(month_fact, "!!bin::2"))
#'
#' # Relocating the months using "@d" in the name
#' table(bin(month_fact, .("@5" = "may", "@1 summer" = "@aug|jul")))
#'
#' # Putting "@" as first item means subsequent items will be placed first
#' table(bin(month_fact, .("@", "aug", "july")))
#'
#' #
#' # "Cutting" numeric data
#' #
#'
#' data(iris)
#' plen = iris$Petal.Length
#'
#' # 3 parts of (roughly) equal size
#' table(bin(plen, "cut::3"))
#'
#' # Three custom bins
#' table(bin(plen, "cut::2]5]"))
#'
#' # .. same, excluding 5 in the 2nd bin
#' table(bin(plen, "cut::2]5["))
#'
#' # Using quartiles
#' table(bin(plen, "cut::q1]q2]q3]"))
#'
#' # Using percentiles
#' table(bin(plen, "cut::p20]p50]p70]p90]"))
#'
#' # Mixing all
#' table(bin(plen, "cut::2[q2]p90]"))
#'
#' # NOTA:
#' # -> the labels always contain the min/max values in each bin
#'
#' # Custom labels can be provided, just give them in the char. vector
#' # NA values lead to the default label
#' table(bin(plen, c("cut::2[q2]p90]", "<2", "]2; Q2]", NA, ">90%")))
#'
#'
#'
#' #
#' # With a formula
#' #
#'
#' data(iris)
#' plen = iris$Petal.Length
#'
#' # We need to use "x"
#' table(bin(plen, list("< 2" = ~x < 2, ">= 2" = ~x >= 2)))
#'
#'
bin = function(x, bin){
    check_arg(x, "vector mbt")

    bin = error_sender(eval_dot(bin), arg_name = "bin")

    check_arg(bin, "list | vector mbt")

    varname = deparse(substitute(x))[1]
    bin_factor(bin, x, varname)
}

#' Refactors a variable
#'
#' Takes a variables of any types, transforms it into a factors, and modifies the values of the factors. Useful in estimations when you want to set some value of a vector as a reference.
#'
#' @inheritSection bin \code{bin} vs \code{ref}
#'
#' @param x A vector of any type (must be atomic though).
#' @param ref A vector or a list, or special binning values (explained later). If a vector, it must correspond to (partially matched) values of the vector \code{x}. The vector \code{x} which will be transformed into a factor and these values will be placed first in the levels. That's the main usage of this function. You can also bin on-the-fly the values of \code{x}, using the same syntax as the function \code{\link[fixest]{bin}}. Here's a description of what bin does: To create a new value from old values, use \code{bin = list("new_value"=old_values)} with \code{old_values} a vector of existing values. You can use \code{.()} for \code{list()}.
#' It accepts regular expressions, but they must start with an \code{"@"}, like in \code{bin="@Aug|Dec"}. It accepts one-sided formulas which must contain the variable \code{x}, e.g. \code{bin=list("<2" = ~x < 2)}.
#' The names of the list are the new names. If the new name is missing, the first value matched becomes the new name. In the name, adding \code{"@d"}, with \code{d} a digit, will relocate the value in position \code{d}: useful to change the position of factors.
#' If the vector \code{x} is numeric, you can use the special value \code{"bin::digit"} to group every \code{digit} element.
#' For example if \code{x} represents years, using \code{bin="bin::2"} creates bins of two years.
#' With any data, using \code{"!bin::digit"} groups every digit consecutive values starting from the first value.
#' Using \code{"!!bin::digit"} is the same but starting from the last value.
#' With numeric vectors you can: a) use \code{"cut::n"} to cut the vector into \code{n} equal parts, b) use \code{"cut::a]b["} to create the following bins: \code{[min, a]}, \code{]a, b[}, \code{[b, max]}.
#' The latter syntax is a sequence of number/quartile (q0 to q4)/percentile (p0 to p100) followed by an open or closed square bracket. You can add custom bin names by adding them in the character vector after \code{'cut::values'}. See details and examples. Dot square bracket expansion (see \code{\link[fixest]{dsb}}) is enabled.
#'
#' @return
#' It returns a factor of the same length as \code{x}, where levels have been modified according to the argument \code{ref}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' To bin the values of a vectors: \code{\link[fixest]{bin}}.
#'
#' @examples
#'
#' data(airquality)
#'
#' # A vector of months
#' month_num = airquality$Month
#' month_lab = c("may", "june", "july", "august", "september")
#' month_fact = factor(month_num, labels = month_lab)
#' table(month_num)
#' table(month_fact)
#'
#' #
#' # Main use
#' #
#'
#' # Without argument: equivalent to as.factor
#' ref(month_num)
#'
#' # Main usage: to set a level first:
#' # (Note that partial matching is enabled.)
#' table(ref(month_fact, "aug"))
#'
#' # You can rename the level on-the-fly
#' # (Northern hemisphere specific!)
#' table(ref(month_fact, .("Hot month"="aug",
#'                         "Late summer" = "sept")))
#'
#'
#' # Main use is in estimations:
#' a = feols(Petal.Width ~ Petal.Length + Species, iris)
#'
#' # We change the reference
#' b = feols(Petal.Width ~ Petal.Length + ref(Species, "vers"), iris)
#'
#' etable(a, b)
#'
#'
#' #
#' # Binning
#' #
#'
#' # You can also bin factor values on the fly
#' # Using @ first means a regular expression will be used to match the values.
#' # Note that the value created is placed first.
#' # To avoid that behavior => use the function "bin"
#' table(ref(month_fact, .(summer = "@jul|aug|sep")))
#'
#' # Please refer to the example in the bin help page for more example.
#' # The syntax is the same.
#'
#'
#' #
#' # Precise relocation
#' #
#'
#' # You can place a factor at the location you want
#' #  by adding "@digit" in the name first:
#' table(ref(month_num, .("@5"=5)))
#'
#' # Same with renaming
#' table(ref(month_num, .("@5 five"=5)))
#'
#'
ref = function(x, ref){
    check_arg(x, "vector mbt")

    if(missing(ref)){
        return(as.factor(x))
    }

    ref = error_sender(eval_dot(ref), arg_name = "ref")
    check_arg(ref, "list | vector mbt")

    varname = deparse(substitute(x))[1]

    IS_SPECIAL = FALSE
    if(!is.list(ref)){
        if(is.character(ref[1]) && grepl("^(cut|bin)", ref[1])){
            IS_SPECIAL = TRUE
            if(!is.numeric(x)){
                stop(.dsb("To use the special binning '.[ref[1]]' the variable ",
                          "'.[varname]' must be numeric. Currently this is not the case ",
                          "(it is of class .[3KO, C?class(x)] instead)."))
            }
        } else {
            ref = as.list(ref)
        }
    }

    if(!IS_SPECIAL && !is.factor(x)){
        x = as.factor(x)
    }

    if(!IS_SPECIAL && ref[[1]] != "@"){
        ref_new = list("@")
        index = 1:length(ref) + 1
        ref_new[index] = ref
        names(ref_new)[index] = names(ref)
        ref = ref_new
    }

    bin_factor(ref, x, varname)
}



#' Expands formula macros
#'
#' Create macros within formulas and expand them with character vectors or other formulas.
#'
#' @inheritParams setFixest_fml
#'
#' @param fml A formula containing macros variables. Each macro variable must start with two dots. The macro variables can be set globally using \code{setFixest_fml}, or can be defined in \code{...}. Special macros of the form \code{..("regex")} can be used to fetch, through a regular expression, variables directly in a character vector (or in column names) given in the argument \code{data}. square brackets have a special meaning: Values in them are evaluated and parsed accordingly. Example: \code{y~x.[1:2] + z.[i]} will lead to \code{y~x1+x2+z3} if \code{i==3}. See examples.
#' @param lhs If present then a formula will be constructed with \code{lhs} as the full left-hand-side. The value of \code{lhs} can be a one-sided formula, a call, or a character vector. Note that the macro variables wont be applied. You can use it in combination with the argument \code{rhs}. Note that if \code{fml} is not missing, its LHS will be replaced by \code{lhs}.
#' @param rhs If present, then a formula will be constructed with \code{rhs} as the full right-hand-side. The value of \code{rhs} can be a one-sided formula, a call, or a character vector. Note that the macro variables wont be applied. You can use it in combination with the argument \code{lhs}. Note that if \code{fml} is not missing, its RHS will be replaced by \code{rhs}.
#' @param data Either a character vector or a data.frame. This argument will only be used if a macro of the type \code{..("regex")} is used in the formula of the argument \code{fml}. If so, any variable name from \code{data} that matches the regular expression will be added to the formula.
#'
#' @details
#' In \code{xpd}, the default macro variables are taken from \code{getFixest_fml}. Any value in the \code{...} argument of \code{xpd} will replace these default values.
#'
#' The definitions of the macro variables will replace in verbatim the macro variables. Therefore, you can include multi-part formulas if you wish but then beware of the order of the macros variable in the formula. For example, using the \code{airquality} data, say you want to set as controls the variable \code{Temp} and \code{Day} fixed-effects, you can do \code{setFixest_fml(..ctrl = ~Temp | Day)}, but then \code{feols(Ozone ~ Wind + ..ctrl, airquality)} will be quite different from \code{feols(Ozone ~ ..ctrl + Wind, airquality)}, so beware!
#'
#' @section Dot square bracket operator in formulas:
#'
#' In a formula, the dot square bracket (DSB) operator can: i) create manifold variables at once, or ii) capture values from the current environment and put them verbatim in the formula.
#'
#' Say you want to include the variables \code{x1} to \code{x3} in your formula. You can use \code{xpd(y ~ x.[1:3])} and you'll get \code{y ~ x1 + x2 + x3}.
#'
#' To summon values from the environment, simply put the variable in square brackets. For example: \code{for(i in 1:3) xpd(y.[i] ~ x)} will create the formulas \code{y1 ~ x} to \code{y3 ~ x} depending on the value of \code{i}.
#'
#' You can include a full variable from the environment in the same way: \code{for(y in c("a", "b")) xpd(.[y] ~ x)} will create the two formulas \code{a ~ x} and \code{b ~ x}.
#'
#' The DSB can even be used within variable names, but then the variable must be nested in character form. For example \code{y ~ .["x.[1:2]_sq"]} will create \code{y ~ x1_sq + x2_sq}. Using the character form is important to avoid a formula parsing error. Double quotes must be used. Note that the character string that is nested will be parsed with the function \code{\link[fixest]{dsb}}, and thus it will return a vector.
#'
#' By default, the DSB operator expands vectors into sums. You can add a comma, like in \code{.[, x]}, to expand with commas--the content can then be used within functions. For instance: \code{c(x.[, 1:2])} will create \code{c(x1, x2)} (and \emph{not} \code{c(x1 + x2)}).
#'
#' In all \code{fixest} estimations, this special parsing is enabled, so you don't need to use \code{xpd}.
#'
#' You can even use multiple square brackets within a single variable, but then the use of nesting is required. For example, the following \code{xpd(y ~ .[".[letters[1:2]]_.[1:2]"])} will create \code{y ~ a_1 + b_2}. Remember that the nested character string is parsed with \code{\link[fixest]{dsb}}, which explains this behavior.
#'
#' @section Regular expressions:
#'
#' You can catch several variable names at once by using regular expressions. To use regular expressions, you need to enclose it in the dot-dot or the regex function: \code{..("regex")} or \code{regex("regex")}. For example, \code{regex("Sepal")} will catch both the variables \code{Sepal.Length} and \code{Sepal.Width} from the \code{iris} data set. In a \code{fixest} estimation, the variables names from which the regex will be applied come from the data set. If you use \code{xpd}, you need to provide either a data set or a vector of names in the argument \code{data}.
#'
#' By default the variables are aggregated with a sum. For example in a data set with the variables x1 to x10, \code{regex("x(1|2)"} will yield \code{x1 + x2 + x10}. You can instead ask for "comma" aggregation by using a comma first, just before the regular expression: \code{y ~ sw(regex(,"x(1|2)"))} would lead to \code{y ~ sw(x1, x2, x10)}.
#'
#' Note that the dot square bracket operator (DSB, see before) is applied before the regular expression is evaluated. This means that \code{regex("x.[3:4]_sq")} will lead, after evaluation of the DSB, to \code{regex("x3_sq|x4_sq")}. It is a handy way to insert range of numbers in a regular expression.
#'
#'
#' @return
#' It returns a formula where all macros have been expanded.
#'
#' @author
#' Laurent Berge
#'
#'
#' @seealso
#' \code{\link[fixest]{setFixest_fml}} to set formula macros, and \code{\link[fixest]{dsb}} to modify character strings with the DSB operator.
#'
#' @examples
#'
#' # Small examples with airquality data
#' data(airquality)
#' # we set two macro variables
#' setFixest_fml(..ctrl = ~ Temp + Day,
#'               ..ctrl_long = ~ poly(Temp, 2) + poly(Day, 2))
#'
#' # Using the macro in lm with xpd:
#' lm(xpd(Ozone ~ Wind + ..ctrl), airquality)
#' lm(xpd(Ozone ~ Wind + ..ctrl_long), airquality)
#'
#' # You can use the macros without xpd() in fixest estimations
#' a = feols(Ozone ~ Wind + ..ctrl, airquality)
#' b = feols(Ozone ~ Wind + ..ctrl_long, airquality)
#' etable(a, b, keep = "Int|Win")
#'
#'
#' # Using .[]
#'
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' i = 2:3
#' z = "species"
#' lm(xpd(y ~ x.[2:3] + .[z]), base)
#'
#' # No xpd() needed in feols
#' feols(y ~ x.[2:3] + .[z], base)
#'
#' #
#' # You can use xpd for stepwise estimations
#' #
#'
#' # Note that for stepwise estimations in fixest, you can use
#' # the stepwise functions: sw, sw0, csw, csw0
#' # -> see help in feols or in the dedicated vignette
#'
#' # we want to look at the effect of x1 on y
#' # controlling for different variables
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # We first create a matrix with all possible combinations of variables
#' my_args = lapply(names(base)[-(1:2)], function(x) c("", x))
#' (all_combs = as.matrix(do.call("expand.grid", my_args)))
#'
#' res_all = list()
#' for(i in 1:nrow(all_combs)){
#'   res_all[[i]] = feols(xpd(y ~ x1 + ..v, ..v = all_combs[i, ]), base)
#' }
#'
#' etable(res_all)
#' coefplot(res_all, group = list(Species = "^^species"))
#'
#' #
#' # You can use macros to grep variables in your data set
#' #
#'
#' # Example 1: setting a macro variable globally
#'
#' data(longley)
#' setFixest_fml(..many_vars = grep("GNP|ployed", names(longley), value = TRUE))
#' feols(Armed.Forces ~ Population + ..many_vars, longley)
#'
#' # Example 2: using ..("regex") or regex("regex") to grep the variables "live"
#'
#' feols(Armed.Forces ~ Population + ..("GNP|ployed"), longley)
#'
#' # Example 3: same as Ex.2 but without using a fixest estimation
#'
#' # Here we need to use xpd():
#' lm(xpd(Armed.Forces ~ Population + regex("GNP|ployed"), data = longley), longley)
#'
#' # Stepwise estimation with regex: use a comma after the parenthesis
#' feols(Armed.Forces ~ Population + sw(regex(,"GNP|ployed")), longley)
#'
#' # Multiple LHS
#' etable(feols(..("GNP|ployed") ~ Population, longley))
#'
#'
#' #
#' # lhs and rhs arguments
#' #
#'
#' # to create a one sided formula from a character vector
#' vars = letters[1:5]
#' xpd(rhs = vars)
#'
#' # Alternatively, to replace the RHS
#' xpd(y ~ 1, rhs = vars)
#'
#' # To create a two sided formula
#' xpd(lhs = "y", rhs = vars)
#'
#'
#' #
#' # Dot square bracket operator
#' #
#'
#' # You can create multiple variables at once
#' xpd(y ~ x.[1:5] + z.[2:3])
#'
#' # You can summon variables from the environment
#' var = "a"
#' xpd(y ~ x.[var])
#'
#' # ... the variables can be multiple
#' vars = LETTERS[1:3]
#' xpd(y ~ x.[vars])
#'
#' # You can have "complex" variable names but they must be nested in character form
#' xpd(y ~ .["x.[vars]_sq"])
#'
#' # DSB can be used within regular expressions
#' re = c("GNP", "Pop")
#' xpd(Unemployed ~ regex(".[re]"), data = longley)
#'
#' # => equivalent to regex("GNP|Pop")
#'
#' # Use .[,var] (NOTE THE COMMA!) to expand with commas
#' # !! can break the formula if missused
#' vars = c("wage", "unemp")
#' xpd(c(y.[,1:3]) ~ csw(.[,vars]))
#'
#'
#' # Example of use of .[] within a loop
#' res_all = list()
#' for(p in 1:3){
#'   res_all[[p]] = feols(Ozone ~ Wind + poly(Temp, .[p]), airquality)
#' }
#'
#' etable(res_all)
#'
#' # The former can be compactly estimated with:
#' res_compact = feols(Ozone ~ Wind + sw(.[, "poly(Temp, .[1:3])"]), airquality)
#'
#' etable(res_compact)
#'
#' # How does it work?
#' # 1)  .[, stuff] evaluates stuff and, if a vector, aggregates it with commas
#' #     Comma aggregation is done thanks to the comma placed after the square bracket
#' #     If .[stuff], then aggregation is with sums.
#' # 2) stuff is evaluated, and if it is a character string, it is evaluated with
#' # the function dsb which expands values in .[]
#' #
#' # Wrapping up:
#' # 2) evaluation of dsb("poly(Temp, .[1:3])") leads to the vector:
#' #    c("poly(Temp, 1)", "poly(Temp, 2)", "poly(Temp, 3)")
#' # 1) .[, c("poly(Temp, 1)", "poly(Temp, 2)", "poly(Temp, 3)")] leads to
#' #    poly(Temp, 1), poly(Temp, 2), poly(Temp, 3)
#' #
#' # Hence sw(.[, "poly(Temp, .[1:3])"]) becomes:
#' #       sw(poly(Temp, 1), poly(Temp, 2), poly(Temp, 3))
#'
#'
xpd = function(fml, ..., lhs, rhs, data = NULL){
    .xpd(fml = fml, ..., lhs = lhs, rhs = rhs, data = data, check = TRUE, macro = TRUE, frame = parent.frame())
}

.xpd = function(fml, ..., lhs, rhs, data = NULL, check = FALSE, macro = FALSE, frame = NULL){

    is_lhs = !missing(lhs)
    is_rhs = !missing(rhs)
    if(is_lhs || is_rhs){

        if(check) check_arg(fml, .type = "formula", .up = 1)

        # Direct formula creation
        if(missing(fml)){
            res = if(is_lhs) 1 ~ 1 else ~ 1
        } else {
            if(is_lhs && length(fml)[1] == 2){
                # Means: RHS formula and we add the LHS
                res = 1 ~ 1
                res[[3]] = fml[[2]]

            } else {
                res = fml
            }
        }

        if(is_lhs){
            res[[2]] = value2stringCall(lhs, call = TRUE, check = check)
        }

        if(is_rhs){
            res[[length(res)]] = value2stringCall(rhs, call = TRUE, check = check)
        }

        fml = res

        if(!macro) return(fml)

        # NOTA:
        # if we allow for macro implementation ex post:
        # This entails a 50% performance drop in terms of speed.
        # Now, without macro variables, speed is at 30us while it was 20us before
        # so in .xpd => macro argument

    } else if(check){
        check_arg(fml, .type = "formula mbt", .up = 1)
    }

    macros = parse_macros(..., from_xpd = TRUE, check = check, frame = frame)

    is_macro = length(macros) != 0
    is_data = !missnull(data)
    fml_funs = all.vars(fml, functions = TRUE)
    is_brackets = "[" %in% fml_funs
    is_regex = any(c("regex", "..") %in% fml_funs)
    if(!(is_macro || is_data || is_brackets || is_regex)) return(fml)

    fml_dp = NULL

    if(is_macro){
        # We allow recursion + auto-blocking at 5th depth
        # to allow crazy things from the user side (but not too much)
        max_depth = 5
        depth  = 0
        any_done = TRUE
        qui = which(names(macros) %in% all.vars(fml))
        while(length(qui) > 0 && any_done && depth < max_depth){
            depth = depth + 1
            fml_dp_origin = fml_dp = deparse_long(fml)
            # We need to add a lookafter assertion: otherwise if we have ..ctrl + ..ctrl_long, there will be a replacement in ..ctrl_long
            for(i in qui){
                fml_dp = gsub(paste0(escape_regex(names(macros)[i]), "(?=$|[^[:alnum:]_\\.])"), macros[[i]], fml_dp, perl = TRUE)
            }

            fml = as.formula(fml_dp, frame)

            qui = which(names(macros) %in% all.vars(fml))
            any_done = fml_dp_origin != fml_dp
        }

        if(depth == max_depth && length(qui) > 0){
            warning(dsb("In xpd, max recursivity has been reached. ",
                        "Please revise the recursive definition of your variables. ",
                        "It concerns the variable.[*s_, q, C?names(macros)[qui]]."))
        }
    }

    if(is_regex){
        # We expand only if data is provided (it means the user wants us to check)
        # if .[]: we expand inside the ..(".[1:3]") => ..("1|2|3")

        if(is.null(fml_dp)) fml_dp = deparse_long(fml)
        is_brackets = grepl(".[", fml_dp, fixed = TRUE)

        if(is_data || is_brackets){

            if(is.null(fml_dp)) fml_dp = deparse_long(fml)

            if(is_data){
                check_arg(data, "character vector no na | matrix | data.frame")

                if(is.matrix(data)){
                    data = colnames(data)
                    if(is.null(data)){
                        stop("The argument 'data' must contain variable names. It is currently a matrix without column names.")
                    }
                } else if(is.data.frame(data)){
                    data = names(data)
                }
            }

            fml_dp_split = strsplit(fml_dp, '(?!<[[:alnum:]._])(regex|\\.\\.)\\((?=[,"])',
                                    perl = TRUE)[[1]]

            res = fml_dp_split
            for(i in 2:length(res)){

                re = sub('"\\).*', "", res[i])

                re_width = nchar(re)

                is_comma = grepl("^,", re)
                re = gsub("^,? ?\"", "", re)

                re = dot_square_bracket(re, frame, regex = TRUE)

                if(is_data){
                    vars = grep(re, data, value = TRUE, perl = TRUE)
                    if(length(vars) == 0){
                        vars = "1"
                    }

                    coll = if(is_comma) ", " else " + "

                    res[i] = paste0(paste(vars, collapse = coll), substr(res[i], re_width + 3, nchar(res[i])))
                } else {
                    res[i] = paste0("..(\"", re, "\")", substr(res[i], re_width + 3, nchar(res[i])))
                }

            }

            fml_txt = paste(res, collapse = "")
            fml = error_sender(as.formula(fml_txt, frame),
                               "Expansion of variables in ..(\"regex\"): coercion of the following text to a formula led to an error.\n",
                               fit_screen(paste0("     TEXT: ", fml_txt)),
                               "\n  PROBLEM: see below", up = 1)
        }
    }

    if("[" %in% all.vars(fml, functions = TRUE)){
        fml_txt = deparse_long(fml)
        if(grepl(".[", fml_txt, fixed = TRUE)){
            fml_txt = dot_square_bracket(fml_txt, frame)

            # level 1 recursivity
            if(grepl(".[", fml_txt, fixed = TRUE)){
                fml_txt = dot_square_bracket(fml_txt, frame)
            }

            fml = error_sender(as.formula(fml_txt, frame),
                               "Dot square bracket operator: coercion of the following text to a formula led to an error.\n",
                               fit_screen(paste0("     TEXT: ", fml_txt)),
                               "\n  PROBLEM: see below", up = 1)
        }

    }

    fml
}


#' Fast transform of any type of vector(s) into an integer vector
#'
#' Tool to transform any type of vector, or even combination of vectors, into an integer vector ranging from 1 to the number of unique values. This actually creates an unique identifier vector.
#'
#' @param ... Vectors of any type, to be transformed in integer.
#' @param sorted Logical, default is \code{FALSE}. Whether the integer vector should make reference to sorted values?
#' @param add_items Logical, default is \code{FALSE}. Whether to add the unique values of the original vector(s). If requested, an attribute \code{items} is created containing the values (alternatively, they can appear in a list if \code{items.list=TRUE}).
#' @param items.list Logical, default is \code{FALSE}. Only used if \code{add_items=TRUE}. If \code{TRUE}, then a list of length 2 is returned with \code{x} the integer vector and \code{items} the vector of items.
#' @param multi.df Logical, default is \code{FALSE}. If \code{TRUE} then a data.frame listing the unique elements is returned in the form of a data.frame. Ignored if \code{add_items = FALSE}.
#' @param multi.join Character scalar used to join the items of multiple vectors. The default is \code{"_"}. Ignored if \code{add_items = FALSE}.
#' @param internal Logical, default is \code{FALSE}. For programming only. If this function is used within another function, setting \code{internal = TRUE} is needed to make the evaluation of \code{...} valid. End users of \code{to_integer} should not care.
#'
#'
#' @return
#' Reruns a vector of the same length as the input vectors.
#' If \code{add_items=TRUE} and \code{items.list=TRUE}, a list of two elements is returned: \code{x} being the integer vector and \code{items} being the unique values to which the values in \code{x} make reference.
#'
#' @examples
#'
#' x1 = iris$Species
#' x2 = as.integer(iris$Sepal.Length)
#'
#' # transforms the species vector into integers
#' to_integer(x1)
#'
#' # To obtain the "items":
#' to_integer(x1, add_items = TRUE)
#' # same but in list form
#' to_integer(x1, add_items = TRUE, items.list = TRUE)
#'
#' # transforms x2 into an integer vector from 1 to 4
#' to_integer(x2, add_items = TRUE)
#'
#' # To have the sorted items:
#' to_integer(x2, add_items = TRUE, sorted = TRUE)
#'
#' # The result can safely be used as an index
#' res = to_integer(x2, add_items = TRUE, sorted = TRUE, items.list = TRUE)
#' all(res$items[res$x] == x2)
#'
#'
#' #
#' # Multiple vectors
#' #
#'
#' to_integer(x1, x2, add_items = TRUE)
#'
#' # You can use multi.join to handle the join of the items:
#' to_integer(x1, x2, add_items = TRUE, multi.join = "; ")
#'
to_integer = function(..., sorted = FALSE, add_items = FALSE, items.list = FALSE,
                      multi.df = FALSE, multi.join = "_", internal = FALSE){

    if(!internal) check_arg(..., "vector mbt")
    check_arg(sorted, add_items, items.list, "logical scalar")
    check_arg(multi.join, "character scalar")

    dots = list(...)

    # Removing NAs
    Q = length(dots)
    n_all = lengths(dots)
    n = n_all[1]

    if(length(unique(n_all)) != 1) stop("All elements in '...' should be of the same length (current lenghts are ", enumerate_items(n_all), ").")

    is_na = is.na(dots[[1]])
    for(q in seq(from = 2, length.out = Q - 1)){
        is_na = is_na | is.na(dots[[q]])
    }

    ANY_NA = FALSE
    if(any(is_na)){
        ANY_NA = TRUE

        if(all(is_na)){
            message("NOTE: All values are NA.")
            res = rep(NA, n)
            if(add_items){
                if(items.list){
                    res = list(x = res, items = NA)
                } else {
                    attr(res, "items") = NA
                }
            }

            return(res)
        }

        for(q in 1:Q) dots[[q]] = dots[[q]][!is_na]
    }

    #
    # Creating the ID
    #

    if(Q == 1){
        if(sorted && !is.numeric(dots[[1]]) && !is.character(dots[[1]])){
            # general way => works for any type with a sort method
            f = dots[[1]]
            res_raw = quickUnclassFactor(f, addItem = TRUE, sorted = FALSE)
            obs_1st = cpp_get_first_item(res_raw$x, length(res_raw$items))
            f_unik = f[obs_1st]
            f_order = order(f_unik)
            x_new = order(f_order)[res_raw$x]
            if(add_items){
                items_new = f_unik[f_order]
                res = list(x = x_new, items = items_new)
            } else {
                res = x_new
            }

        } else {
            res = quickUnclassFactor(dots[[1]], addItem = add_items, sorted = sorted)
        }

    } else {

        QUF_raw = list()
        for(q in 1:Q){
            QUF_raw[[q]] = quickUnclassFactor(dots[[q]], sorted = FALSE, addItem = TRUE)
        }

        # Then we combine
        power = floor(1 + log10(sapply(QUF_raw, function(x) length(x$items))))

        is_large = sum(power) > 14
        if(is_large){
            # 15 Aug 2021, finally found a solution. It was so obvious with hindsight...
            QUF_raw_value = lapply(QUF_raw, `[[`, 1)
            order_index = do.call(order, QUF_raw_value)
            index = cpp_combine_clusters(QUF_raw_value, order_index)
        } else {
            # quicker, but limited by the precision of doubles
            index = QUF_raw[[1]]$x
            for(q in 2:Q){
                index = index + QUF_raw[[q]]$x*10**sum(power[1:(q-1)])
            }
        }

        res = quickUnclassFactor(index, addItem = add_items || sorted, sorted = sorted)

        if(add_items || sorted){
            # we re order appropriately
            # f prefix means factor

            obs_1st = cpp_get_first_item(res$x, length(res$items))

            f_all = list()
            for(q in 1:Q){
                f_all[[q]] = dots[[q]][obs_1st]
            }

            f_order = do.call("order", f_all)

            x_new = order(f_order)[res$x]

            if(multi.df){
                # Putting into a DF => we take care of names
                mc_dots = match.call(expand.dots = FALSE)[["..."]]
                n_dots = length(mc_dots)
                mc_dots_names = names(mc_dots)
                if(is.null(mc_dots_names)) mc_dots_names = character(n_dots)

                my_names = character(n_dots)
                for(q in 1:n_dots){
                    if(nchar(mc_dots_names[q]) > 0){
                        my_names[q] = mc_dots_names[q]
                    } else {
                        my_names[q] = deparse_long(mc_dots[[q]])
                    }
                }

                names(f_all) = my_names

                f_df = as.data.frame(f_all)
                items_new = f_df[f_order, , drop = FALSE]
                row.names(items_new) = 1:nrow(items_new)
            } else {
                # we "paste" them
                arg_list = f_all
                arg_list$sep = multi.join
                f_char = do.call("paste", arg_list)
                items_new = f_char[f_order]
            }

            if(add_items){
                res = list(x = x_new, items = items_new)
            } else {
                res = x_new
            }
        }
    }

    if(ANY_NA){
        if(is.list(res)){
            x = res$x
        } else {
            x = res
        }

        x_na = rep(NA, n)
        x_na[!is_na] = x

        if(is.list(res)){
            res$x = x_na
        } else {
            res = x_na
        }

    }

    if(add_items && isFALSE(items.list)){
        res_tmp = res$x
        attr(res_tmp, "items") = res$items
        res = res_tmp
    }

    res
}



#' Centers a set of variables around a set of factors
#'
#' User-level access to internal demeaning algorithm of \code{fixest}.
#'
#' @inheritSection feols Varying slopes
#'
#' @param X A matrix, vector, data.frame or a list OR a formula OR a \code{\link[fixest]{feols}} estimation. If equal to a formula, then the argument \code{data} is required, and it must be of the type: \code{x1 + x2 ~ f1 + fe2} with on the LHS the variables to be centered, and on the RHS the factors used for centering. Note that you can use variables with varying slopes with the syntax \code{fe[v1, v2]} (see details in \code{\link[fixest]{feols}}). If a \code{feols} estimation, all variables (LHS+RHS) are demeaned and then returned (only if it was estimated with fixed-effects). Otherwise, it must represent the data to be centered. Of course the number of observations of that data must be the same as the factors used for centering (argument \code{f}).
#' @param f A matrix, vector, data.frame or list. The factors used to center the variables in argument \code{X}. Matrices will be coerced using \code{as.data.frame}.
#' @param slope.vars A vector, matrix or list representing the variables with varying slopes. Matrices will be coerced using \code{as.data.frame}. Note that if this argument is used it MUST be in conjunction with the argument \code{slope.flag} that maps the factors to which the varying slopes are attached. See examples.
#' @param slope.flag An integer vector of the same length as the number of variables in \code{f} (the factors used for centering). It indicates for each factor the number of variables with varying slopes to which it is associated. Positive values mean that the raw factor should also be included in the centering, negative values that it should be excluded. Sorry it's complicated... but see the examples it may get clearer.
#' @param data A data.frame containing all variables in the argument \code{X}. Only used if \code{X} is a formula, in which case \code{data} is mandatory.
#' @param weights Vector, can be missing or NULL. If present, it must contain the same number of observations as in \code{X}.
#' @param nthreads Number of threads to be used. By default it is equal to \code{getFixest_nthreads()}.
#' @param notes Logical, whether to display a message when NA values are removed. By default it is equal to \code{getFixest_notes()}.
#' @param iter Number of iterations, default is 2000.
#' @param tol Stopping criterion of the algorithm. Default is \code{1e-6}. The algorithm stops when the maximum absolute increase in the coefficients values is lower than \code{tol}.
#' @param na.rm Logical, default is \code{TRUE}. If \code{TRUE} and the input data contains any NA value, then any observation with NA will be discarded leading to an output with less observations than the input. If \code{FALSE}, if NAs are present the output will also be filled with NAs for each NA observation in input.
#' @param as.matrix Logical, if \code{TRUE} a matrix is returned, if \code{FALSE} it will be a data.frame. The default depends on the input, if atomic then a matrix will be returned.
#' @param im_confident Logical, default is \code{FALSE}. FOR EXPERT USERS ONLY! This argument allows to skip some of the preprocessing of the arguments given in input. If \code{TRUE}, then \code{X} MUST be a numeric vector/matrix/list (not a formula!), \code{f} MUST be a list, \code{slope.vars} MUST be a list, \code{slope.vars} MUST be consistent with \code{slope.flag}, and \code{weights}, if given, MUST be numeric (not integer!). Further there MUST be not any NA value, and the number of observations of each element MUST be consistent. Non compliance to these rules may simply lead your R session to break.
#' @param fixef.reorder Logical, default is \code{TRUE}. Whether to reorder the fixed-effects by frequencies before feeding them into the algorithm. If \code{FALSE}, the original fixed-effects order provided by the user is maintained. In general, reordering leads to faster and more precise performance.
#' @param ... Not currently used.
#'
#' @return
#' It returns a data.frame of the same number of columns as the number of variables to be centered.
#'
#' If \code{na.rm = TRUE}, then the number of rows is equal to the number of rows in input minus the number of NA values (contained in \code{X}, \code{f}, \code{slope.vars} or \code{weights}). The default is to have an output of the same number of observations as the input (filled with NAs where appropriate).
#'
#' A matrix can be returned if \code{as.matrix = TRUE}.
#'
#' @examples
#'
#' # Illustration of the FWL theorem
#' data(trade)
#'
#' base = trade
#' base$ln_dist = log(base$dist_km)
#' base$ln_euros = log(base$Euros)
#'
#' # We center the two variables ln_dist and ln_euros
#' #  on the factors Origin and Destination
#' X_demean = demean(X = base[, c("ln_dist", "ln_euros")],
#'                   f = base[, c("Origin", "Destination")])
#' base[, c("ln_dist_dm", "ln_euros_dm")] = X_demean
#'
#' est = feols(ln_euros_dm ~ ln_dist_dm, base)
#' est_fe = feols(ln_euros ~ ln_dist | Origin + Destination, base)
#'
#' # The results are the same as if we used the two factors
#' # as fixed-effects
#' etable(est, est_fe, se = "st")
#'
#' #
#' # Variables with varying slopes
#' #
#'
#' # You can center on factors but also on variables with varying slopes
#'
#' # Let's have an illustration
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' #
#' # We center y and x1 on species and x2 * species
#'
#' # using a formula
#' base_dm = demean(y + x1 ~ species[x2], data = base)
#'
#' # using vectors
#' base_dm_bis = demean(X = base[, c("y", "x1")], f = base$species,
#'                      slope.vars = base$x2, slope.flag = 1)
#'
#' # Let's look at the equivalences
#' res_vs_1 = feols(y ~ x1 + species + x2:species, base)
#' res_vs_2 = feols(y ~ x1, base_dm)
#' res_vs_3 = feols(y ~ x1, base_dm_bis)
#'
#' # only the small sample adj. differ in the SEs
#' etable(res_vs_1, res_vs_2, res_vs_3, keep = "x1")
#'
#' #
#' # center on x2 * species and on another FE
#'
#' base$fe = rep(1:5, 10)
#'
#' # using a formula => double square brackets!
#' base_dm = demean(y + x1 ~ fe + species[[x2]], data = base)
#'
#' # using vectors => note slope.flag!
#' base_dm_bis = demean(X = base[, c("y", "x1")], f = base[, c("fe", "species")],
#'                      slope.vars = base$x2, slope.flag = c(0, -1))
#'
#' # Explanations slope.flag = c(0, -1):
#' # - the first 0: the first factor (fe) is associated to no variable
#' # - the "-1":
#' #    * |-1| = 1: the second factor (species) is associated to ONE variable
#' #    *   -1 < 0: the second factor should not be included as such
#'
#' # Let's look at the equivalences
#' res_vs_1 = feols(y ~ x1 + i(fe) + x2:species, base)
#' res_vs_2 = feols(y ~ x1, base_dm)
#' res_vs_3 = feols(y ~ x1, base_dm_bis)
#'
#' # only the small sample adj. differ in the SEs
#' etable(res_vs_1, res_vs_2, res_vs_3, keep = "x1")
#'
#'
#'
#'
demean = function(X, f, slope.vars, slope.flag, data, weights,
                  nthreads = getFixest_nthreads(), notes = getFixest_notes(),
                  iter = 2000, tol = 1e-6, fixef.reorder = TRUE, na.rm = TRUE,
                  as.matrix = is.atomic(X),
                  im_confident = FALSE, ...) {


    ANY_NA = FALSE
    # SK: to reassign class if X is data.frame, data.table or tibble. Optimally you would preserve all attributes,
    # but using attributes<- is slow on data.frames. What I did in collapse is export the SET_ATTRIB and DUPLICATE_ATTRIB
    # functions from the C-API to use them internally in R -> copy attributes without checks at 0 cost, even for large data.frames.

    # LB: next line is needed if input data is matrix and as.matrix is set to FALSE
    clx = NULL
    is_fixest = inherits(X, "fixest")
    if(lX <- is.list(X) && !is_fixest) {
        clx <- oldClass(X)
        # SK: oldClass is faster and returns NULL when the list is plain. class returns the implicit class "list".
        # This is the key to fast R code -> all data.frame methods are super slow and should duly be avoided in internal code.
        oldClass(X) <- NULL
    }
    # LB: The reassignment of attributes to data.frames is actually a very good idea, thanks Seb!

    # LB: To avoid delayed evaluation problems (due to new default is.atomic(X))
    as_matrix = as.matrix

    # internal argument
    fe_info = FALSE

    # Step 1: formatting the input
    if(!im_confident){

        check_arg(X, "numeric vmatrix | list | formula | class(fixest) mbt")
        check_arg(iter, "integer scalar GE{1}")
        check_arg(tol, "numeric scalar GT{0}")
        check_arg(notes, "logical scalar")

        validate_dots(valid_args = "fe_info", stop = TRUE)
        dots = list(...)
        fe_info = isTRUE(dots$fe_info)

        data_mbt = TRUE
        if(inherits(X, "fixest")){
            data_mbt = FALSE
            if(!identical(X$method, "feols")){
                stop("This function only works for 'feols' estimations (not for ", X$method, ").")
            }

            if(!fe_info && !is.null(X$y_demeaned)){
                return(cbind(X$y_demeaned, X$X_demeaned))
            }

            if(!"fixef_vars" %in% names(X)){
                stop("To apply demeaning, the estimation must contain fixed-effects, this is currently not the case.")
            }

            # Otherwise => we create data and X: formula
            data = fetch_data(X)
            fml_linear = X$fml_all$linear
            fml_fixef = X$fml_all$fixef

            fml_vars = .xpd(~ ..lhs + ..rhs,
                            ..lhs = fml_linear[[2]],
                            ..rhs = fml_linear[[3]])

            if(!is.null(X$fml_all$iv)){
                fml_iv = X$fml_all$iv
                fml_vars = .xpd(~ ..vars + ..iv_lhs + ..iv_rhs,
                                ..vars = fml_vars,
                                ..iv_lhs = fml_iv[[2]],
                                ..iv_rhs = fml_iv[[3]])
            }

            fml = .xpd(lhs = fml_vars, rhs = fml_fixef)

            mc = match.call()
            if(!"tol" %in% names(mc)) tol = X$fixef.tol
            if(!"iter" %in% names(mc)) iter = X$fixef.iter

            weights = X[["weights"]]

            X = fml

            as_matrix = TRUE
        }

        #
        # X
        #

        fe_done = FALSE
        if(is.call(X)) {

            if(data_mbt) check_arg(data, "data.frame mbt")
            check_arg(X, "ts formula var(data)", .data = data)

            # Extracting the information
            terms_fixef = fixef_terms(.xpd(rhs = X[[3L]]))
            # We add the intercept only for only slope models, otherwise this would be void since the result would be always 0
            X = fixest_model_matrix(.xpd(lhs = quote(y), rhs = X[[2L]]), data, fake_intercept = any(terms_fixef$slope_flag >= 0))
            var_names = dimnames(X)[[2]]

            lX = FALSE # Needed for the rest of the code to work

            # FE
            fe_done = TRUE
            f = unclass(error_sender(prepare_df(terms_fixef$fe_vars, data),
                                     "Error when evaluating the fixed-effects variables: "))

            isSlope = any(terms_fixef$slope_flag != 0)

            if(isSlope) {

                slope.vars = unclass(error_sender(prepare_df(terms_fixef$slope_vars, data),
                                                  "Error when evaluating the variable with varying slopes: "))

                if(anyDuplicated(terms_fixef$slope_vars)){
                    # we need it!!!!
                    slope.vars = slope.vars[terms_fixef$slope_vars]
                }

                slope.flag = terms_fixef$slope_flag

                is_numeric = vapply(`attributes<-`(slope.vars, NULL), is.numeric, TRUE)
                if(!all(is_numeric)){
                    stop("In the fixed-effects part of the formula, variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope.vars)[!is_numeric], "s.is.quote"), " not.")
                }
            } else {
                slope.vars = list(0)
                slope.flag = rep(0L, length(f))
            }

        } else if(lX) {
            var_names = names(X)
            if(!any(clx == "data.frame")) {
                n_all = lengths(X)
                if(!all(n_all == n_all[1L])) {
                    n_unik = unique(n_all)
                    stop("In argument 'X' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
                }
            }

            # SK: If your using Rcpp and the input is NumericVector or NumericMatrix, you should get a C++ error for wrongly typed data, so I don't see the need for this
            # LB: I disagree, such error messages are poor and hard to debug for end users.

            is_numeric = vapply(`attributes<-`(X, NULL), is.numeric, TRUE)
            if(!all(is_numeric)){
                # Faster than any(!is_numeric) => LB: yet a few nano seconds saved... ;-)
                n_non_num = sum(!is_numeric)
                stop("All values in 'X' must be numeric. This is not the case for ", n_non_num, " variable", plural(n_non_num), ".")
            }

        } else {
            # Here: numeric matrix or vector
            var_names = dimnames(X)[[2L]]
        }

        # nobs: useful for checks
        nobs = if(is.list(X)) length(.subset2(X, 1L)) else NROW(X)

        #
        # f + slope.vars
        #

        if(!fe_done){
            check_arg(f, "vmatrix | list mbt")
            check_arg(slope.vars, "numeric vmatrix | list")
            check_arg(slope.flag, "integer vector no na")

            if(is.list(f)) {
                if(!is.data.frame(f)) {
                    n_all = lengths(f)
                    if(!all(n_all == n_all[1L])) {
                        n_unik = unique(n_all)
                        stop("In argument 'f' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
                    }
                    # f = as.data.frame(f)
                    # SK: as.data.frame is super slow, especially on large data. You already checked the lengths, so it's ok
                    # LB: true
                } else {
                    oldClass(f) <- NULL
                }

            } else {
                # SK: as.data.frame on matrix is also very slow, you could do in C++ as in at collapse::mctl.. but I suppose most will pass lists of factors anyway..
                # LB: In the long run I'll deal with that but that's really low priority
                f = if(!is.array(f)) list(f) else unclass(as.data.frame(f))
            }

            is_pblm = vapply(`attributes<-`(f, NULL), function(x) !(is.numeric(x) || is.character(x)), TRUE)
            if(any(is_pblm)) {
                for(i in which(is_pblm)) f[[i]] = as.character(f[[i]])
            }

            if(length(f[[1L]]) != nobs){
                stop("The number of observations in 'X' and in 'f' don't match (", if(lX) length(X[[1L]]) else NROW(X), " vs ", length(f[[1L]]), ").")
            }

            # Now the slopes
            isSlope = FALSE
            if(!missnull(slope.vars) || !missnull(slope.flag)) {

                isSlope = TRUE

                if(missnull(slope.vars)) stop("If argument 'slope.flag' is provided, you must also provide the argument 'slope.vars'.")
                if(missnull(slope.flag)) stop("If argument 'slope.vars' is provided, you must also provide the argument 'slope.flag'.")

                if(is.list(slope.vars)){
                    if(!is.data.frame(slope.vars)){
                        n_all = lengths(slope.vars)
                        if(!all(n_all == n_all[1L])) {
                            n_unik = unique(n_all)
                            stop("In argument 'slope.vars' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
                        }

                    } else {
                        oldClass(slope.vars) <- NULL
                    }

                    is_numeric = vapply(`attributes<-`(slope.vars, NULL), is.numeric, TRUE)
                    # SK: Much faster than !sapply(slope.vars, is.numeric)
                    # LB: 10us diff, only impedes readability IMO
                    if(!all(is_numeric)) stop("In the argument 'slope.vars', the variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope.vars)[!is_numeric], "s.is.quote"), " not.")

                } else {
                    # SK: Same as above.. Necessary to convert ??
                    # LB: The new C++ code accepts lists of vector/matrices or vector/matrices directly // So if we're here (meaning vector or matrix), that's fine
                }

                if(length(slope.flag) != length(f)) stop("The argument 'slope.flag' must be a vector representing, for each fixed-effect, the number of variables with varying slopes associated to it. Problem: the lengths of slope.flag and the fixed-effect differ: ", length(slope.flag), " vs ", length(f), ".")


                # SK: if(sum(abs(slope.flag)) != length(slope.vars))  : This means that slope flag can only be 0, 1 or -1. In the documentation you only talk about positive and negative values. The change I made reflects this.
                # LB: Cmon Sebastian why change sum(abs(slope.flag)) != length(slope.vars) into sum(slope.flag != 0L) != length(slope.vars)?  The time gains lie in a handful of nano second and you've just introduced a bug!

                n_sv = if(is.list(slope.vars)) length(slope.vars) else NCOL(slope.vars)
                if(sum(abs(slope.flag)) != n_sv) stop("The argument 'slope.flag' must be a vector representing, for each fixed-effect, the number of variables with varying slopes associated to it. Problem: currently the number of variables in 'slope.flag' differ from the number of variables in 'slope.vars': ", sum(abs(slope.flag)), " vs ", n_sv, ".")

            } else {
                slope.vars = list(0)
                slope.flag = rep(0L, length(f))
            }
        }


        ## weights
        check_arg_plus(weights, "NULL numeric conv vector len(value) GE{0}", .value = nobs)
        is_weight = !missing(weights) && !is.null(weights)

        ## nthreads (Also seems unnecessarily costly: 250 microseconds)
        # LB: I know... but it has to be checked
        if(!missing(nthreads)) nthreads = check_set_nthreads(nthreads)

        #
        # FORMATTING
        #

        # NAs
        is_NA = FALSE
        info_X = cpp_which_na_inf(X, nthreads)

        if(info_X$any_na_inf){
            is_NA = info_X$is_na_inf
            n_na_x = 1L
        } else n_na_x = 0L


        if(anyNA(f)){
            is_na_fe = !complete.cases(f)
            n_na_fe = 1L
            is_NA = is_NA | is_na_fe
        } else n_na_fe = 0L


        is_na_slope = 0L
        if(isSlope){
            info_slopes = cpp_which_na_inf(slope.vars, nthreads)

            if(info_slopes$any_na_inf){
                is_na_slope = info_slopes$is_na_inf
                is_NA = is_NA | is_na_slope
            }
        }


        is_na_weights = 0L
        if(is_weight){
            info_weights = cpp_which_na_inf(weights, nthreads)

            if(info_weights$any_na_inf){
                is_na_weights = info_weights$is_na_inf
                is_NA = is_NA | is_na_weights
            }
        }

        n_na = sum(is_NA)
        if(n_na && notes) {
            # Here is all the error message stuff now: Only evaluated if needed
            if(n_na_x) n_na_x = sum(info_X$is_na_inf)
            if(n_na_fe) n_na_fe = sum(is_na_fe)
            slopes_msg = if(isSlope) paste0(", slope vars: ", sum(is_na_slope)) else ""
            weight_msg = if(is_weight) paste0(", weights: ", sum(is_na_weights)) else ""

            if(n_na == length(is_NA)) stop("All observations contain NA values (Breakup: X: ", n_na_x, ", f: ", n_na_fe, slopes_msg, weight_msg, ").")
            message("NOTE: ", fsignif(n_na), " observation", plural(n_na), " removed because of NA values (Breakup: X: ", n_na_x, ", f: ", n_na_fe, slopes_msg, weight_msg, ").")
        }

        if(n_na) {
            ANY_NA = TRUE
            # SK: subsetting integer is faster than logical !!
            # LB: Indeed, just checked, quite a diff. Good to know!
            cc <- which(!is_NA)
            X = select_obs(X, cc)

            f = lapply(f, `[`, cc)
            if(isSlope) slope.vars = lapply(slope.vars, `[`, cc)
            if(is_weight) weights = weights[cc]
        }

        if(!is_weight) weights = 1

    } else {
        # Need this here, otherwise error:
        var_names = if(lX) names(X) else dimnames(X)[[2L]]
        if(missing(weights) || is.null(weights)) weights = 1
        if(missnull(slope.vars) || missnull(slope.flag)){
            slope.vars = list(0)
            slope.flag = rep(0L, length(unclass(f)))
            # SK: unclass gives extra speed if f is data.frame, and no cost if f is list.
        }
    }

    #
    # Unclassing fes
    #

    quf_info_all = quf_table_sum(x = f, y = 0, do_sum_y = FALSE, rm_0 = FALSE,
                                 rm_1 = FALSE, rm_single = FALSE, do_refactor = FALSE,
                                 r_x_sizes = 0L, obs2keep = 0L, only_slope = slope.flag < 0L,
                                 nthreads = as.integer(nthreads))

    # table/sum_y/sizes
    fixef_table = quf_info_all$table
    fixef_sizes = lengths(fixef_table)

    if(fixef.reorder){
        new_order = order(fixef_sizes, decreasing = TRUE)
        if(is.unsorted(new_order)){
            fixef_table = fixef_table[new_order]
            fixef_sizes = fixef_sizes[new_order]
            quf_info_all$quf = quf_info_all$quf[new_order]

            # We reorder slope.vars only if needed (because it is a pain)
            if(sum(slope.flag != 0) >= 2){

                # here slope.vars have 2 or more variables:
                # either matrix or data.frame/list
                if(!is.list(slope.vars)){
                    slope.vars = as.data.frame(slope.vars)
                }
                # we ensure it's a plain list
                slope.vars = unclass(slope.vars)

                n_fixef = length(slope.flag)
                slope_vars_list = vector("list", n_fixef)
                id_current = 0
                for(i in 1:n_fixef){
                    n_slopes = abs(slope.flag[i])

                    if(n_slopes == 0){
                        slope_vars_list[[i]] = 0
                    } else {
                        slope_vars_list[[i]] = slope.vars[1:n_slopes + id_current]
                        id_current = id_current + n_slopes
                    }
                }

                # Reordering => currently it's quite clumsy (but I HATE VSs!)
                slope.vars = vector("list", sum(abs(slope.flag)))
                slope_vars_list_reordered = slope_vars_list[new_order]
                id_current = 0
                for(i in 1:n_fixef){
                    vs = slope_vars_list_reordered[[i]]
                    if(is.list(vs)){
                        n_vs = length(vs)
                        for(j in 1:n_vs){
                            slope.vars[[id_current + j]] = vs[[j]]
                        }
                        id_current = id_current + n_vs
                    }
                }

                # slope.flag must be reordered afterwards
                slope.flag = slope.flag[new_order]
            }
        }
    }

    fixef_table_vector = unlist(fixef_table)
    if(!is.integer(slope.flag)) slope.flag = as.integer(slope.flag)

    #
    # The algorithm
    #

    # y => list, X => matrix
    if(as_matrix){
        y = 0
    } else {
        # Quirk => y returns a list only if NCOL(y) > 1 or is.list(y)
        y = if(lX) X else list(X)
        # SK: This does the same, a bit more frugal
        # LB: I tend to prefer late evaluations instead of lX which requires bookkeeping
        X = 0
    }

    vars_demean = cpp_demean(y, X, weights, iterMax = iter,
                             diffMax = tol, r_nb_id_Q = fixef_sizes,
                             fe_id_list = quf_info_all$quf, table_id_I = fixef_table_vector,
                             slope_flag_Q = slope.flag, slope_vars_list = slope.vars,
                             r_init = 0, nthreads = nthreads)

    # Internal call
    if(fe_info){

        res = list(y = vars_demean$y_demean, X = vars_demean$X_demean,
                   weights = weights, fixef_id_list = quf_info_all$quf,
                   slope_vars = slope.vars, slope_flag = slope.flag, varnames = colnames(X))

        return(res)

    }


    if(as_matrix) {
        if(ANY_NA && na.rm == FALSE) {
            K = NCOL(vars_demean$X_demean)
            res = matrix(NA_real_, length(is_NA), K)
            res[cc, ] = vars_demean$X_demean

        } else {
            res = vars_demean$X_demean
        }

        dimnames(res)[[2L]] = var_names

    } else {
        if(ANY_NA && na.rm == FALSE) {
            n = length(is_NA)
            # SK: One-line efficient solution to the task:
            res = lapply(vars_demean$y_demean, function(x) `[<-`(rep(NA_real_, n), cc, value = x))
            # LB: nice compactness

        } else {
            res = vars_demean$y_demean
        }

        names(res) = var_names
        # SK: This makes your function class-agnostic to lists, data.frame's, data.table's and tibbles (with the exception for any custom data.frame row.names which are not preserved, but the as.data.frame solution also did not do that)
        if(!lX || length(clx)) attr(res, "row.names") <- .set_row_names(length(res[[1L]]))
        # SK: Here !lX || length(clx) means add row.names if either X was a matrix or is classed that is a data.frame
        oldClass(res) <- if(lX) clx else "data.frame"
        # SK: If X is a plain list, and since oldClass returns NULL, this will assign a null class again.
    }
    res
}

#' Extracts the observations used for the estimation
#'
#' This function extracts the observations used in \code{fixest} estimation.
#'
#' @param x A \code{fixest} object.
#'
#' @return
#' It returns a simple vector of integers.
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#' base$y[1:5] = NA
#'
#' # Split sample estimations
#' est_split = feols(y ~ x1, base, split = ~species)
#' (obs_setosa = obs(est_split$setosa))
#' (obs_versi = obs(est_split$versicolor))
#'
#' est_versi = feols(y ~ x1, base, subset = obs_versi)
#'
#' etable(est_split, est_versi)
#'
#'
#'
#'
obs = function(x){
    check_arg(x, "class(fixest)")

    if(isTRUE(x$lean)){
        stop("obs() does not work with models estimated with 'lean = TRUE'.")
    }

    id = 1:x$nobs_origin

    for(i in seq_along(x$obs_selection)){
        id = id[x$obs_selection[[i]]]
    }

    return(id)
}




#' Check the fixed-effects convergence of a \code{feols} estimation
#'
#' Checks the convergence of a \code{feols} estimation by computing the first-order conditions of all fixed-effects (all should be close to 0)
#'
#' @param x A \code{\link[fixest]{feols}} estimation that should contain fixed-effects.
#' @param object An object returned by \code{check_conv_feols}.
#' @param type Either "short" (default) or "detail". If "short", only the maximum absolute FOC are displayed, otherwise the 2 smallest and the 2 largest FOC are reported for each fixed-effect and each variable.
#' @param ... Not currently used.
#'
#' Note that this function first re-demeans the variables, thus possibly incurring some extra computation time.
#'
#' @return
#' It returns a list of \code{N} elements, \code{N} being the number of variables in the estimation (dependent variable + explanatory variables +, if IV, endogenous variables and instruments). For each variable, all the first-order conditions for each fixed-effect are returned.
#'
#' @examples
#'
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' base$FE = rep(1:30, 5)
#'
#' # one estimation with fixed-effects + varying slopes
#' est = feols(y ~ x1 | species[x2] + FE[x3], base)
#'
#' # Checking the convergence
#' conv = check_conv_feols(est)
#'
#' # We can check that al values are close to 0
#' summary(conv)
#'
#' summary(conv, "detail")
#'
#'
#'
check_conv_feols = function(x){

    check_arg(x, "class(fixest) mbt")
    if(!identical(x$method, "feols")){
        stop("This function only works for 'feols' estimations (not for ", x$method, ").")
    }

    if(!"fixef_vars" %in% names(x)){
        message("This function only works with fixed-effects (which are currently not present).")
        return(NULL)
    }

    # fixef names for information
    new_order = x$fe.reorder
    fixef_vars = x$fixef_vars[new_order]
    if(!is.null(x$fixef_terms)){

        slope_variables = x$slope_variables_reordered
        slope_flag = x$slope_flag_reordered

        # We reconstruct the terms
        fixef_names = c()
        start = c(0, cumsum(abs(slope_flag)))
        for(i in seq_along(slope_flag)){
            sf = slope_flag[i]
            if(sf >= 0){
                fixef_names = c(fixef_names, fixef_vars[i])
            }

            if(abs(sf) > 0){
                fixef_names = c(fixef_names, paste0(fixef_vars[i], "[[", names(slope_variables)[start[i] + 1:abs(sf)], "]]"))
            }
        }
    } else {
        fixef_names = fixef_vars
    }


    info = demean(x, fe_info = TRUE)

    res = check_conv(y = info$y, X = info$X, fixef_id_list = info$fixef_id_list,
                     slope_flag = info$slope_flag, slope_vars = info$slope_vars,
                     weights = info$weights, full = TRUE, fixef_names = fixef_names)

    names(res) = info$varnames

    class(res) = "fixest_check_conv"

    res
}

#' Prints the number of unique elements in a data set
#'
#' This utility tool displays the number of unique elements in one or multiple data.frames as well as their number of NA values.
#'
#' @param x A formula, with data set names on the LHS and variables on the RHS, like \code{data1 + data2 ~ var1 + var2}. The following special variables are admitted: \code{"."} to get default values, \code{".N"} for the number of observations, \code{".U"} for the number of unique rows, \code{".NA"} for the number of rows with at least one NA. Variables can be combined with \code{"^"}, e.g. \code{df~id^period}; use \code{id\%^\%period} to also include the terms on both sides. Note that using \code{:} and \code{*} is equivalent to \code{^} and \code{\%^\%}. Sub select with \code{id[cond]}, when doing so \code{id} is automatically included. Conditions can be chained, as in \code{id[cond1, cond2]}. Use \code{NA(x, y)} in conditions instead of \code{is.na(x) | is.na(y)}. Use the \code{!!} operator to have both a condition and its opposite. To compare the keys in two data sets, use \code{data1:data2}. If not a formula, \code{x} can be: a vector (displays the # of unique values); a \code{data.frame} (default values are displayed), or a "sum" of data sets like in \code{x = data1 + data2}, in that case it is equivalent to \code{data1 + data2 ~ .}.
#' @param ... Not currently used.
#'
#' @section Special values and functions:
#'
#' In the formula, you can use the following special values: \code{"."}, \code{".N"}, \code{".U"}, and \code{".NA"}.
#'
#' \itemize{
#'
#' \item{\code{"."}}{Accesses the default values. If there is only one data set and the data set is \emph{not} a \code{data.table}, then the default is to display the number of observations and the number of unique rows. If the data is a \code{data.table}, the number of unique items in the key(s) is displayed instead of the number of unique rows (if the table has keys of course). If there are two or more data sets, then the default is to display the unique items for: a) the variables common across all data sets, if there's less than 4, and b) if no variable is shown in a), the number of variables common across at least two data sets, provided there are less than 5. If the data sets are data tables, the keys are also displayed on top of the common variables. In any case, the number of observations is always displayed.}
#'
#' \item{\code{".N"}}{Displays the number of observations.}
#'
#' \item{\code{".U"}}{Displays the number of unique rows.}
#'
#' \item{\code{".NA"}}{Displays the number of rows with at least one NA.}
#'
#' }
#'
#' @section The \code{NA} function:
#'
#' The special function \code{NA} is an equivalent to \code{is.na} but can handle several variables. For instance, \code{NA(x, y)} is equivalent to \code{is.na(x) | is.na(y)}. You can add as many variables as you want as arguments. If no argument is provided, as in \code{NA()}, it is identical to having all the variables of the data set as argument.
#'
#' @section Combining variables:
#'
#' Use the "hat", \code{"^"}, operator to combine several variables. For example \code{id^period} will display the number of unique values of id x period combinations.
#'
#' Use the "super hat", \code{"\%^\%"}, operator to also include the terms on both sides. For example, instead of writing \code{id + period + id^period}, you can simply write \code{id\%^\%period}.
#'
#' Alternatively, you can use \code{:} for \code{^} and \code{*} for \code{\%^\%}.
#'
#' @section Sub-selections:
#'
#' To show the number of unique values for sub samples, simply use \code{[]}. For example, \code{id[x > 10]} will display the number of unique \code{id} for which \code{x > 10}.
#'
#' Simple square brackets lead to the inclusion of both the variable and its subset. For example \code{id[x > 10]} is equivalent to \code{id + id[x > 10]}. To include only the sub selection, use double square brackets, as in \code{id[[x > 10]]}.
#'
#' You can add multiple sub selections at once, only separate them with a comma. For example \code{id[x > 10, NA(y)]} is equivalent to \code{id[x > 10] + id[NA(y)]}.
#'
#' Use the double negative operator, i.e. \code{!!}, to include both a condition and its opposite at once. For example \code{id[!!x > 10]} is equivalent to \code{id[x > 10, !x > 10]}. Double negative operators can be chained, like in \code{id[!!cond1 & !!cond2]}, then the cardinal product of all double negatived conditions is returned.
#'
#' @return
#' It returns a vector containing the number of unique values per element. If several data sets were provided, a list is returned, as long as the number of data sets, each element being a vector of unique values.
#'
#' @examples
#'
#' data = base_did
#' data$x1.L1 = round(lag(x1~id+period, 1, data))
#'
#' # By default, just the formatted number of observations
#' n_unik(data)
#'
#' # Or the nber of unique elements of a vector
#' n_unik(data$id)
#'
#' # number of unique id values and id x period pairs
#' n_unik(data ~.N + id + id^period)
#'
#' # use the %^% operator to include the terms on the two sides at once
#' # => same as id*period
#' n_unik(data ~.N + id %^% period)
#'
#' # using sub selection with []
#' n_unik(data ~.N + period[!NA(x1.L1)])
#'
#' # to show only the sub selection: [[]]
#' n_unik(data ~.N + period[[!NA(x1.L1)]])
#'
#' # you can have multiple values in [],
#' # just separate them with a comma
#' n_unik(data ~.N + period[!NA(x1.L1), x1 > 7])
#'
#' # to have both a condition and its opposite,
#' # use the !! operator
#' n_unik(data ~.N[!!NA(x1.L1)])
#'
#' # the !! operator works within condition chains
#' n_unik(data ~.N[!!NA(x1.L1) & !!x1 > 7])
#'
#' # Conditions can be distributed
#' n_unik(data ~ (id + period)[x1 > 7])
#'
#' #
#' # Several data sets
#' #
#'
#' # Typical use case: merging
#' # Let's create two data sets and merge them
#'
#' data(base_did)
#' base_main = base_did
#' base_extra = sample_df(base_main[, c("id", "period")], 100)
#' base_extra$id[1:10] = 111:120
#' base_extra$period[11:20] = 11:20
#' base_extra$z = rnorm(100)
#'
#' # You can use db1:db2 to compare the common keys in two data sets
#'  n_unik(base_main:base_extra)
#'
#' tmp = merge(base_main, base_extra, all.x = TRUE, by = c("id", "period"))
#'
#' # You can show unique values for any variable, as before
#' n_unik(tmp + base_main + base_extra ~ id[!!NA(z)] + id^period)
#'
#'
#'
n_unik = function(x){
    # returns a vector with the nber of unique values
    # attr("na.info") => nber of NA values, vector

    if(missing(x)){
        stop("Argument 'x' must be provided. Problem: it is missing.")
    }

    # Non standard-evaluation
    x_dp = deparse_long(substitute(x))
    if(!grepl("~", x_dp, fixed = TRUE)){
        if(grepl("[+:*]", x_dp)){
            # Non standard-evaluation
            x = as.formula(paste0(x_dp, "~ ."))
        } else {
            check_arg(x, "data.frame | vector l0 | ts formula")
        }
    } else {
        check_arg(x, "ts formula")
    }

    nthreads = getFixest_nthreads()

    # If vector
    comp_pairs = list()
    if(is.vector(x)){

        x_name = gsub("^[[:alpha:]\\.][[:alnum:]\\._]*\\$", "", x_dp)

        na.info = 0
        if(anyNA(x)){
            who_NA = is.na(x)
            na.info = sum(who_NA)
            x = x[!who_NA]
        }

        res = len_unique(x, nthreads)
        names(res) = x_name

        attr(res, "na.info") = na.info

        class(res) = "vec_n_unik"
        return(res)

    } else if(is.data.frame(x)){

        n_x = 1
        x_all = list(x)
        fml = ~.

    } else if(inherits(x, "formula")){

        x_terms = terms(x[1:2])

        fact_mat = attr(x_terms, "factors")

        x_all_names = rownames(fact_mat)

        # We get the comparison pairs
        vars = colnames(fact_mat)
        for(pair in grep(":", vars, fixed = TRUE, value = TRUE)){
            dict = setNames(1:length(x_all_names), x_all_names)
            pair_split = strsplit(pair, ":", fixed = TRUE)[[1]]
            comb = combn(pair_split, 2)
            for(i in 1:ncol(comb)){
                comp_pairs[[length(comp_pairs) + 1]] = unname(dict[comb[, i]])
            }
        }

        if(length(comp_pairs) > 0){
            comp_pairs = unique(comp_pairs)
        }

        # Construction of the list  + sanity check
        n_x = length(x_all_names)
        x_all = vector("list", n_x)
        for(i in 1:n_x){
            x_all[[i]] = error_sender(eval(str2lang(x_all_names[i]), parent.frame()),
                                      "The left-hand-side of the formula must contain valid data.frames. Problem in the evaluation of '", x_all_names[i], "':")
            check_value(x_all[[i]], "data.frame",
                        .message = paste0("The value '", x_all_names[i], "' (in the LHS of the formula) must be a data.frame."))
        }

        fml = x[c(1, 3)]

    }

    # If DF + formula
    # ex: fml = ~ id^year + author_id[sw0(is.na(author_name), is.na(affil_name), year == min_year)]
    # fml = ~ id^sw0(fe)

    # check variable names
    naked_vars = origin_vars = all.vars(fml)
    valid_vars = c(".N", ".U", ".NA", ".", unique(unlist(lapply(x_all, names))))

    check_value_plus(naked_vars, "multi match", .choices = valid_vars,
                     .message = paste0("The formula must only use variables in the data set", plural_len(x_all), "."))

    if("." %in% naked_vars){
        # The default values

        IS_DT = requireNamespace("data.table", quietly = TRUE)

        dot_default = ".N"
        x_keys = c()
        x_common_vars = c()

        # keys
        if(IS_DT){
            for(I in 1:n_x){
                if(inherits(x_all[[I]], "data.table")){
                    x_keys = c(x_keys, data.table::key(x_all[[I]]))
                }
            }
        }

        # common variables
        if(n_x > 1){

            # Common to all
            x_names_current = names(x_all[[1]])
            for(i in 2:n_x){
                x_names_current = intersect(x_names_current, names(x_all[[i]]))
                if(length(x_names_current) == 0) break
            }

            if(length(x_names_current) %in% 1:4){
                # No more than 4 common variables by default
                x_common_vars = x_names_current

                if(length(comp_pairs) > 0 && length(x_names_current) <= 3){
                    x_common_vars = paste0(x_names_current, collapse = "*")
                }

            } else if(n_x > 2){
                # Common to at least 2 data sets
                for(i in 1:(n_x - 1)){
                    x_names_current = names(x_all[[i]])
                    for(j in (i + 1):n_x){
                        qui_common = x_names_current %in% names(x_all[[j]])
                        if(any(qui_common)){
                            x_common_vars = c(x_common_vars, x_names_current[qui_common])
                        }
                    }
                }

                x_common_vars = unique(x_common_vars)

                if(length(x_common_vars) > 5){
                    #  we keep max 5
                    x_common_vars = c()
                }
            }
        }

        if(length(x_keys) > 0 || length(x_common_vars) > 0){
            dot_default = unique(c(dot_default, x_keys, x_common_vars))
        }

        if(length(dot_default) == 1){
            dot_default = c(".N", ".U")
        }

        rhs_txt = as.character(fml)[[2]]
        rhs_txt = gsub("(^| )\\.( |$)", dsb(".['+'c? dot_default]"), rhs_txt)
        fml = .xpd(rhs = rhs_txt)
    }

    if(any(naked_vars != origin_vars)){
        # we complete the partial matching
        fml_txt = deparse_long(fml)

        var_diff = which(naked_vars != origin_vars)

        for(i in var_diff){
            re = paste0("(?!<[[:alnum:]._])", origin_vars[i], "(?!=[[:alnum:]._])")
            fml_txt = gsub(re, naked_vars[i], fml_txt, perl = TRUE)
        }

        fml = as.formula(fml_txt)
    }

    tm = terms_hat(fml, n_unik = TRUE)
    all_vars = attr(tm, "term.labels")
    var_final = c()

    # stepwise function
    for(var in all_vars){
        # var = all_vars[1]

        if(grepl("combine_clusters(_fast)?", var)){
            var = gsub("combine_clusters(_fast)?", "to_integer", var)
        }

        IS_HAT_SW = grepl("sw0?\\)\\(", var)
        if(IS_HAT_SW){
            var = paste0(gsub("(sw0?)\\)\\(", "\\1(", var), ")")
        }

        if(grepl("^[[:alpha:].][[:alnum:]._]*\\[", var, perl = TRUE) && !grepl("sw0?\\(", var)){
            # We use sw to make id[cond1, cond2] work
            # That's what's called path dependency!

            var_new = gsub("^([[:alpha:].][[:alnum:]._]*)\\[", "\\1[sw0(", var)
            if(grepl("sw0([", var_new, fixed = TRUE)){
                var_new = sub("sw0([", "sw(", str_trim(var_new, -1), fixed = TRUE)
            }
            var_new = sub("\\]$", ")]", var_new)

            var = var_new
        }

        if(is_fun_in_char(var, "sw0?")){

            info = extract_fun(var, "sw0?", err_msg = "The stepwise function can be used only once per variable.")

            sw_value = eval(str2lang(info$fun))

            qui_double = grepl("!!", trimws(sw_value))
            while(any(qui_double)){
                i = which(qui_double)[1]

                value = sw_value[i]
                value_split = strsplit(value, "!!")[[1]]

                n_split = length(value_split)
                n_s = 2 ** (n_split - 1)
                new_values = rep(value_split[1], n_s)

                prefixes = do.call(expand.grid, rep(list(c("", "!")), n_split - 1))

                for(j in 1:ncol(prefixes)){
                    new_values = paste0(new_values, prefixes[, j], value_split[j + 1])
                }

                sw_value = insert(sw_value[-i], new_values, i)
                qui_double = grepl("!!", sw_value)
            }

            var_char_new = paste0(info$before, sw_value, info$after)

            if(IS_HAT_SW){
                var_char_new[1] = gsub(", (,|\\))", "\\1", var_char_new[1])

                qui_nested = grepl("^(to_integer\\(.+){2,}", var_char_new)
                if(any(qui_nested)){
                    # later => add check problem if another function is used
                    var_char_new[qui_nested] = paste0("to_integer(", gsub("to_integer\\(|\\)", "", var_char_new[qui_nested]), ")")
                }
            }

            if(grepl("[]", var_char_new[1], fixed = TRUE)){
                var_char_new[1] = gsub("[]", "", var_char_new[1], fixed = TRUE)
            }

            var_final = c(var_final, var_char_new)
        } else {
            var_final = c(var_final, var)
        }
    }

    var_final = unique(var_final)

    var_final_names = var_final
    # we take care of id^year cases
    if(any(who_combine <- is_fun_in_char(var_final, "to_integer"))){

        for(i in which(who_combine)){
            info = extract_fun(var_final[i], "to_integer")

            # we use sw() to get the vars
            sw_char = gsub("^[^\\(]+\\(", "sw(", info$fun)
            sw_value = eval(str2lang(sw_char))

            hat_value = paste(sw_value, collapse = "^")

            var_final_names[i] = paste0(info$before, hat_value, info$after)
        }
    }

    # NA counting + unique counting

    res_all = list()
    for(I in 1:n_x){

        x = x_all[[I]]
        x_list = unclass(x)
        x_list[["NA_fun"]] = function(...) NA_fun(..., df = x)

        vars_legit = c(".N", ".U", ".NA", names(x))

        res = c()
        na.info = c()

        for(i in seq_along(var_final)){

            vf = var_final[i]
            vf_name = var_final_names[i]
            na_i = 0

            vf = gsub("((?<=[^[:alnum:]_\\.])|^)NA\\(", "NA_fun(", vf, perl = TRUE)

            vf_call = str2lang(vf)

            if(!all(all.vars(vf_call) %in% vars_legit)){
                res_i = NA_real_
                vf_name = ""

            } else if(grepl("^\\.(N|U)($|\\[)", vf)){


                if(vf %in% c(".N", ".N[]")){
                    res_i = nrow(x)
                    vf_name = "# Observations"

                } else if(vf %in% c(".U", ".U[]")){
                    res_i = nrow(unique(x))
                    vf_name = "# Unique rows"

                } else {
                    # Other methods are faster but are less general

                    vf_new = gsub("(^\\.(N|U)\\[)|(\\]$)", "", vf)
                    do_unik = grepl("^\\.U", vf)

                    val = eval(str2lang(vf_new), x_list)

                    # we want to drop the NAs for indices
                    if(anyNA(val)){
                        val = val[!is.na(val)]
                    }

                    if(length(val) == 0){
                        res_i = 0

                    } else if(do_unik){
                        res_i = NROW(unique(x[val, ]))

                    } else {
                        res_i = if(is.logical(val)) sum(val) else length(val)
                    }

                    msg = if(do_unik) "# Unique rows" else "# Obs."
                    vf_new = gsub("NA_fun", "NA", vf_new, fixed = TRUE)
                    vf_name = paste0(msg, " with ", vf_new)
                }

            } else if(grepl("^(\\.NA|NA_fun\\(\\))", vf)) {

                if(vf %in% c(".NA", ".NA[]", "NA_fun()", "NA_fun()[]")){
                    res_i = sum(!complete.cases(x))
                    vf_name = "# Rows with NAs"

                } else {

                    subselect = NULL
                    sub_text = ""
                    if(grepl("\\]$", vf)){
                        sub_text_split = strsplit(vf, "[", fixed = TRUE)[[1]]
                        sub_text = paste0(sub_text_split[-1], collapse = "[")
                        sub_text = gsub("\\]$", "", sub_text)

                        subselect = eval(str2lang(sub_text), x_list)
                    }

                    if(length(subselect) == 0){
                        res_i = 0
                    } else {
                        res_i = sum(!complete.cases(x[subselect, ]))
                    }

                    vf_name = dsb("# Rows with NAs, with .[sub_text]")
                }
            } else {
                val = eval(vf_call, x_list)

                if(anyNA(val)){
                    who_NA = is.na(val)
                    na_i = sum(who_NA)
                    val = val[!who_NA]
                }

                if(grepl("^NA_fun\\(", vf)){
                    res_i = sum(val)

                    # now the message
                    leftover = sub("^[^\\)]+\\)", "", vf)
                    begin = substr(vf, 1, nchar(vf) - nchar(leftover))

                    my_fun = list("NA_fun" = function(x) enumerate_items(sapply(sys.call()[-1], deparse_long)))
                    na_vars = eval(str2lang(begin), my_fun)
                    msg_start = paste0("# NAs in ", na_vars)

                    msg_end = ""
                    if(grepl("[", vf, fixed = TRUE)){
                        msg_end = paste0(" with ", gsub("^\\[|\\]$", "", leftover))
                    }

                    vf_name = paste0(msg_start, msg_end)

                } else {
                    res_i = len_unique(val, nthreads)
                }

            }

            res[vf_name] = res_i
            na.info[i] = na_i
        }

        attr(res, "na.info") = na.info
        class(res) = "vec_n_unik"

        if(n_x == 1){
            return(res)
        }

        res_all[[x_all_names[I]]] = res
    }


    # Specific data set comparisons
    info_pairs = list()
    for(pair in comp_pairs){

        i_x = pair[1]
        i_y = pair[2]

        x = x_all[[i_x]]
        x_list = unclass(x)
        x_list[["NA_fun"]] = function(...) NA_fun(..., df = x)

        y = x_all[[i_y]]
        y_list = unclass(y)
        y_list[["NA_fun"]] = function(...) NA_fun(..., df = y)

        vars_legit = intersect(names(x), names(y))

        # two last elements: id and common
        all_rows_id_common = c()
        row_temp = rep(NA_real_, n_x + 2)

        for(i in seq_along(var_final)){

            vf = var_final[i]
            vf = gsub("((?<=[^[:alnum:]_\\.])|^)NA\\(", "NA_fun(", vf, perl = TRUE)
            vf_call = str2lang(vf)

            if(!all(all.vars(vf_call) %in% vars_legit) || grepl("^NA_fun", vf)){
                next

            } else {

                if(grepl("to_integer", vf, fixed = TRUE)){
                    # specific scheme for combinations

                    if(grepl("NA_fun()", vf, fixed = TRUE)){
                        # we don't allow that
                        next
                    }

                    vars_keep = all.vars(vf_call)

                    xy_list = list()
                    for(v in vars_keep){
                        xy_list[[v]] = c(x_list[[v]], y_list[[v]])
                    }

                    xy_list[["NA_fun"]] = function(...) NA_fun(..., df = stop("Internal error."))


                    val_xy = eval(vf_call, xy_list)

                    nrow_x = length(x_list[[1]])
                    nrow_y = length(y_list[[1]])
                    val_x = val_xy[1:nrow_x]
                    val_y = val_xy[(nrow_x + 1):(nrow_x + nrow_y)]

                } else {
                    val_x = eval(vf_call, x_list)
                    val_y = eval(vf_call, y_list)
                }


                if(anyNA(val_x)){
                    val_x = val_x[!is.na(val_x)]
                }

                if(anyNA(val_y)){
                    val_y = val_y[!is.na(val_y)]
                }

                # first col is ID
                val_x_unik = unique(val_x)
                val_y_unik = unique(val_y)

                n_x_not_in_y = sum(!val_x_unik %in% val_y_unik)
                n_y_not_in_x = sum(!val_y_unik %in% val_x_unik)
                n_common = sum(val_x_unik %in% val_y_unik)

                row = row_temp
                row[i_x] = n_x_not_in_y
                row[i_y] = n_y_not_in_x

                row[n_x + 1] = i
                row[n_x + 2] = n_common

                all_rows_id_common = rbind(all_rows_id_common, row)

            }
        }

        if(length(all_rows_id_common) > 0){
            attr(all_rows_id_common, "data_id") = c(i_x, i_y)
            info_pairs[[length(info_pairs) + 1]] = all_rows_id_common
        }
    }

    if(length(info_pairs) > 0){
        attr(res_all, "info_pairs") = info_pairs
    }

    class(res_all) = "list_n_unik"

    return(res_all)
}

#' Formatted object size
#'
#' Tools that returns a formatted object size, where the appropriate unit is automatically chosen.
#'
#' @param x Any R object.
#' @param ... Not currently used.
#'
#' @return
#' Returns a character scalar.
#'
#' @examples
#'
#' osize(iris)
#'
#' data(trade)
#' osize(trade)
#'
#'
osize = function(x){

    size = as.numeric(utils::object.size(x))
    n = log10(size)

    if (n < 3) {
        res = paste0(size, " Octets.")
    } else if (n < 6) {
        res = paste0(fsignif(size/1000), " Ko.")
    } else {
        res = paste0(fsignif(size/1e+06), " Mo.")
    }

    class(res) = "osize"
    res
}

#' Randomly draws observations from a data set
#'
#' This function is useful to check a data set. It gives a random number of rows of the input data set.
#'
#' @param x A data set: either a vector, a matrix or a data frame.
#' @param n The number of random rows/elements to sample randomly.
#' @param previous Logical scalar. Whether the results of the previous draw should be returned.
#'
#' @return
#' A data base (resp vector) with \code{n} rows (resp elements).
#'
#' @examples
#'
#' sample_df(iris)
#'
#' sample_df(iris, previous = TRUE)
#'
sample_df = function(x, n = 10, previous = FALSE){

    check_arg(n, "integer scalar")
    check_arg(previous, "logical scalar")

    if(MISSNULL(x)){
        if(missing(x)) stop("The argument 'x' cannot be missing, please provide it.")
        if(is.null(x)) stop("The argument 'x' must be a vector, matrix or data.frame. Currently it is NULL.")
    }

    all_draws = getOption("fixest_sample_df")
    x_dp = deparse(substitute(x))

    make_draw = TRUE
    if(previous){
        if(x_dp %in% names(all_draws)){
            draw__ = all_draws[[x_dp]]
            make_draw = FALSE
        } else {
            warning("No previous draw found for this data. Making a new draw.")
        }
    }

    is_unidim = is.null(dim(x)) || inherits(x, "table")
    n_data = if(is_unidim) length(x) else nrow(x)


    n = min(n, n_data)

    if(make_draw){
        # complicated var name to avoid issue with data.table variables
        draw__ = sample(n_data, n)
    }

    # saving
    all_draws[[x_dp]] = draw__
    options(fixest_sample_df = all_draws)

    # returning
    if(is_unidim){
        return(x[draw__])
    } else {
        return(x[draw__, ])
    }
}

len_unique = function(x, nthreads = getFixest_nthreads()){

    if(!is.vector(x)){
        return(length(unique(x)))
    }

    ANY_NA = FALSE
    if(is.numeric(x)){
        info_na = cpp_which_na_inf(x, nthreads)
        if(info_na$any_na_inf){
            ANY_NA = TRUE
            x = x[!info_na$is_na_inf]
        }
    } else {
        if(anyNA(x)){
            ANY_NA = TRUE
            x = x[!is.na(x)]
        }
    }

    if(length(x) == 0){
        return(0)
    }

    x_quf = quickUnclassFactor(x, addItem = TRUE, sorted = FALSE)

    length(x_quf$items) + ANY_NA
}

#' Replicates \code{fixest} objects
#'
#' Simple function that replicates \code{fixest} objects while (optionally) computing different standard-errors. Useful mostly in combination with \code{\link[fixest]{etable}} or \code{\link[fixest]{coefplot}}.
#'
#' @param x Either a \code{fixest} object, either a list of \code{fixest} objects created with \code{.l()}.
#' @param times Integer vector giving the number of repetitions of the vector of elements. By default \code{times = 1}. It must be either of length 1, either of the same length as the argument \code{x}.
#' @param each Integer scalar indicating the repetition of each element. Default is 1.
#' @param vcov A list containing the types of standard-error to be computed, default is missing. If not missing, it must be of the same length as \code{times}, \code{each}, or the final vector. Note that if the arguments \code{times} and \code{each} are missing, then \code{times} becomes equal to the length of \code{vcov}. To see how to summon a VCOV, see the dedicated section in the \href{https://lrberge.github.io/fixest/articles/fixest_walkthrough.html#the-vcov-argument-1}{vignette}.
#' @param ... In \code{.l()}: \code{fixest} objects. In \code{rep()}: not currently used.
#'
#' @details
#' To apply \code{rep.fixest} on a list of \code{fixest} objects, it is absolutely necessary to use \code{.l()} and not \code{list()}.
#'
#' @return
#' Returns a list of the appropriate length. Each element of the list is a \code{fixest} object.
#'
#' @examples
#'
#' # Let's show results with different standard-errors
#'
#' est = feols(Ozone ~ Solar.R + Wind + Temp, data = airquality)
#'
#' my_vcov = list(~ Month, ~ Day, ~ Day + Month)
#'
#' etable(rep(est, vcov = my_vcov))
#'
#' coefplot(rep(est, vcov = my_vcov), drop = "Int")
#'
#' #
#' # To rep multiple objects, you need to use .l()
#' #
#'
#' est_bis = feols(Ozone ~ Solar.R + Wind + Temp | Month, airquality)
#'
#' etable(rep(.l(est, est_bis), vcov = my_vcov))
#'
#' # using each
#' etable(rep(.l(est, est_bis), each = 3, vcov = my_vcov))
#'
#'
rep.fixest = function(x, times = 1, each = 1, vcov, ...){
    # each is applied first, then times
    # x can be either a list of fixest objects, either a fixest object

    check_arg(x, "class(fixest, fixest_list) mbt")
    check_arg(times, "integer scalar GE{1} | integer vector no na GE{0}")
    check_arg(each, "integer scalar GE{1} | logical scalar")
    check_arg(vcov, "class(list)")

    if(is_user_level_call()){
        validate_dots(suggest_args = c("times", "each"), stop = TRUE)
    }

    # Checking the arguments
    IS_LIST = FALSE
    if("fixest_list" %in% class(x)){
        IS_LIST = TRUE
        class(x) = "list"

        n = length(x)

    } else {
        n = 1
    }

    if(is.logical(each) && each == FALSE){
        stop("Argument 'each' cannot be equal to FALSE.")
    }

    IS_MULTI_VCOV = !missing(vcov)
    if(IS_MULTI_VCOV){
        n_vcov = length(vcov)

        if(times == 1 && each == 1){
            if(isTRUE(each)){
                each = n_vcov
            } else {
                times = n_vcov
            }

        }
    }

    res_int = rep(1:n, times = times, each = each)
    n_res = length(res_int)

    if(IS_MULTI_VCOV){
        # Checking and expanding

        vcov_mapping = 1:n_res
        if(times == 1){
            if(n_vcov != each && n_vcov != n_res){
                stop("In rep, the argument 'vcov' (currently of length ", n_vcov, ") must be a list either of length ", each, " or of length ", n_res, ".")
            }

            if(n_vcov == each) vcov_mapping = rep(1:each, times = n)

        } else if(each == 1){
            if(n_vcov != times && n_vcov != n_res){
                stop("In rep, the argument 'vcov' (currently of length ", n_vcov, ") must be a list either of length ", times, " or of length ", n_res, ".")
            }

            if(n_vcov == times) vcov_mapping = rep(1:n_vcov, each = n)

        } else {
            if(n_vcov != n_res){
                stop("In rep, the argument 'vcov' (currently of length ", n_vcov, ") must be a list either of length ", n_res, ".")
            }
        }
    }

    res = vector("list", length(res_int))

    if(IS_MULTI_VCOV){
        for(i in 1:n_res){
            if(IS_LIST){
                res[[i]] = summary(x[[res_int[i]]], vcov = vcov[[vcov_mapping[i]]])
            } else {
                res[[i]] = summary(x, vcov = vcov[[vcov_mapping[i]]])
            }
        }

    } else {
        for(i in unique(res_int)){
            if(IS_LIST){
                res[res_int == i] = x[[i]]
            } else {
                res[res_int == i] = list(x)
            }
        }
    }

    for(i in 1:n_res){
        res[[i]]$model_id = res_int[i]
    }

    res
}

#' @rdname rep.fixest
rep.fixest_list = function(x, times = 1, each = 1, vcov, ...){
    rep.fixest(x, times = times, each = each, vcov = vcov, ...)
}

#' @rdname rep.fixest
.l = function(...){

    check_arg(..., "mbt class(fixest) | list")

    dots = list(...)
    if(all(sapply(dots, function(x) "fixest" %in% class(x)))){
        class(dots) = "fixest_list"

        return(dots)
    }

    if(length(dots) == 1){

        if("fixest_multi" %in% class(dots[[1]])){
            res = attr(dots[[1]], "data")
            class(res) = "fixest_list"
            return(res)
        }

        if(all(sapply(dots[[1]], function(x) "fixest" %in% class(x)))){
            res = dots[[1]]
            class(res) = "fixest_list"
            return(res)
        }
    }

    res = list()
    for(i in seq_along(dots)){
        if("fixest" %in% class(dots[[i]])){
            res[[length(res) + 1]] = dots[[i]]
        } else {
            obj = dots[[i]]

            if(class(obj) %in% "fixest_multi"){
                data = attr(obj, "data")
                for(j in seq_along(data)){
                    res[[length(res) + 1]] = data[[j]]
                }

            } else {
                for(j in seq_along(obj)){
                    if(!"fixest" %in% class(obj[[j]])){
                        stop("In .l(...), each argument must be either a fixest object, or a list of fixest objects. Problem: The ", n_th(j), " element of the ", n_th(i), " argument (the latter being a list) is not a fixest object.")
                    }

                    res[[length(res) + 1]] = obj[[j]]
                }
            }
        }
    }

    class(res) = "fixest_list"
    res
}
