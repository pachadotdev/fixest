#' Extract the Fixed-Effects from a \code{fixest} estimation.
#'
#' This function retrieves the fixed effects from a \code{fixest} estimation. It is useful only when there are one or more fixed-effect dimensions.
#'
#' @inheritParams feNmlm
#'
#' @param object A \code{fixest} estimation (e.g. obtained using \code{\link[fixest2]{feols}} or \code{\link[fixest2]{feglm}}).
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
#' \code{\link[fixest2]{plot.fixest.fixef}}. See also the main estimation functions \code{\link[fixest2]{femlm}}, \code{\link[fixest2]{feols}} or \code{\link[fixest2]{feglm}}. Use \code{\link[fixest2]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest2]{fixef.fixest}} to extract the fixed-effect coefficients, and the function \code{\link[fixest2]{etable}} to visualize the results of multiple estimations.
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
#' est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination + Product, trade)
#'
#' # Obtaining the fixed-effects coefficients:
#' fe_trade <- fixef(est_pois)
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
fixef.fixest <- function(object, notes = getFixest_notes(), sorted = TRUE, nthreads = getFixest_nthreads(),
                         fixef.tol = 1e-5, fixef.iter = 10000, ...) {
  # object is a fixest object
  # This function retrieves the dummies

  check_arg(notes, sorted, "logical scalar")

  # Checking the arguments
  if (is_user_level_call()) {
    validate_dots()
  }

  check_value(fixef.tol, "numeric scalar GT{0} LT{1}")
  check_value(fixef.iter, "strict integer scalar GT{0}")

  if (isTRUE(object$lean)) {
    # LATER: recompute the FEs by extracting them from the data
    stop("Fixed-effects from 'lean' fixest objects cannot be extracted. Please re-estimate with 'lean = FALSE'.")
  }

  # Preliminary stuff
  S <- object$sumFE

  if (is.null(S)) {
    stop("The estimation was done without fixed-effects (FE). The FE coefficients cannot be retrieved.")
  }

  family <- object$family
  fixef_names <- object$fixef_vars

  fixef_id <- object$fixef_id

  Q <- length(fixef_id)
  N <- length(S)

  # either (we need to clean its attributes for unlist to be efficient)
  id_dummies_vect <- list()
  for (i in 1:Q) id_dummies_vect[[i]] <- as.vector(fixef_id[[i]])

  is_ref_approx <- FALSE
  isSlope <- FALSE
  if (!is.null(object$fixef_terms)) {
    isSlope <- TRUE
    # This is an estimation with slopes
    # we apply another method => we use the demeaning function

    slope_variables <- object$slope_variables_reordered
    slope_flag <- object$slope_flag_reordered

    new_order <- object$fe.reorder
    fixef_vars <- object$fixef_vars[new_order]
    fixef_sizes <- as.integer(object$fixef_sizes[new_order])

    # We reconstruct the terms
    fixef_terms <- c()
    start <- c(0, cumsum(abs(slope_flag)))
    for (i in seq_along(slope_flag)) {
      sf <- slope_flag[i]
      if (sf >= 0) {
        fixef_terms <- c(fixef_terms, fixef_vars[i])
      }

      if (abs(sf) > 0) {
        fixef_terms <- c(fixef_terms, paste0(fixef_vars[i], "[[", names(slope_variables)[start[i] + 1:abs(sf)], "]]"))
      }
    }

    fe_id_list <- object$fixef_id[new_order]

    #
    # STEP 2: demeaning
    #

    nthreads <- check_set_nthreads(nthreads)

    table_id_I <- as.integer(unlist(lapply(fe_id_list, table), use.names = FALSE))

    S_demean <- cpp_demean(
      y = S, X_raw = 0, r_weights = 0, iterMax = as.integer(fixef.iter),
      diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
      fe_id_list = fe_id_list, table_id_I = table_id_I,
      slope_flag_Q = slope_flag, slope_vars_list = slope_variables,
      r_init = 0, nthreads = nthreads, save_fixef = TRUE
    )

    fixef_coef <- S_demean$fixef_coef

    names(fixef_sizes) <- fixef_vars

    fe_all <- c()
    for (i in seq_along(slope_flag)) {
      fe_all <- c(fe_all, rep(fixef_vars[i], 1 + abs(slope_flag[i]) - (slope_flag[i] < 0)))
    }

    start <- 1
    i <- 1
    fixef_values <- list()
    for (q in seq_along(slope_flag)) {
      sf <- slope_flag[q]
      if (sf == 0) {
        fixef_values[[i]] <- fixef_coef[seq(start, length.out = fixef_sizes[q])]
        i <- i + 1
        start <- start + fixef_sizes[q]
      } else {
        nb <- abs(sf) + (sf > 0)

        adj <- 0
        if (sf > 0) {
          # The fixed-effects is in the last position
          j_fe <- nb - 1
          fixef_values[[i]] <- fixef_coef[seq(start + j_fe, by = nb, length.out = fixef_sizes[q])]
          adj <- 1
        }

        for (j in 0:(nb - 1 - adj)) {
          fixef_values[[i + j + adj]] <- fixef_coef[seq(start + j, by = nb, length.out = fixef_sizes[q])]
        }
        i <- i + nb
        start <- start + fixef_sizes[q] * nb
      }
    }

    #
    # Now the referenes
    #

    nb_ref <- integer(length(fixef_terms))

    # FE references
    who_fe <- slope_flag >= 0
    Q_fe <- sum(who_fe)
    if (Q_fe >= 2) {
      my_dum <- fe_id_list[who_fe]

      dumMat <- matrix(unlist(my_dum, use.names = FALSE), N, Q_fe) - 1
      orderCluster <- matrix(unlist(lapply(my_dum, order), use.names = FALSE), N, Q_fe) - 1

      nbCluster <- sapply(my_dum, max)

      fixef_values_tmp <- cpp_get_fe_gnl(Q_fe, N, rep(1, N), dumMat, nbCluster, orderCluster)

      # the information on the references
      nb_ref_fe <- fixef_values_tmp[[Q_fe + 1]]
    } else {
      nb_ref_fe <- integer(Q_fe)
    }

    # Slope references (if associated FE + constant)

    names(slope_flag) <- fixef_vars

    Q_slope <- sum(abs(slope_flag))
    nb_ref_slope <- integer(Q_slope)
    i_noVS <- 1
    for (i in seq_along(fixef_terms)) {
      ft <- fixef_terms[i]

      if (!grepl("[[", ft, fixed = TRUE)) {
        # No slope => already computed
        nb_ref[i] <- nb_ref_fe[i_noVS]
        i_noVS <- i_noVS + 1
      } else {
        # Slope
        fe_name <- gsub("\\[.+", "", ft)
        my_dum <- fe_id_list[[fe_name]]

        my_order <- order(my_dum)
        var_sorted <- slope_variables[[gsub(".+\\[|\\]+", "", ft)]][my_order]

        # if no associated FE => we check only 0 values
        if (slope_flag[fe_name] < 0) {
          nb_ref[i] <- cpp_constant_dum(fixef_sizes[fe_name], var_sorted, my_dum[my_order], only_0 = TRUE)
        } else {
          nb_ref[i] <- cpp_constant_dum(fixef_sizes[fe_name], var_sorted, my_dum[my_order])
        }
      }
    }

    # we recreate that to avoid conditioning on isSlope later
    fixef_id <- fixef_id[fe_all]
    fixef_names <- fixef_terms
  } else if (Q == 1) {
    # This is the simplest case
    id <- id_dummies_vect[[1]]

    myOrder <- order(id)
    myDiff <- c(1, diff(id[myOrder]))

    select <- myOrder[myDiff == 1]

    fixef_values <- list(S[select])

    # There are no references => no need to set nb_ref
  } else {
    # We apply a Rcpp script to handle complicated cases (and we don't know beforehand if the input is one)

    dumMat <- matrix(unlist(id_dummies_vect), N, Q) - 1
    orderCluster <- matrix(unlist(lapply(id_dummies_vect, order)), N, Q) - 1

    nbCluster <- sapply(fixef_id, max)

    fixef_values <- cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)

    # The algorithm is fast but may fail on some instance. We need to check
    if (any(fixef_values[[Q + 1]] > 1) && Q >= 3) {
      # we re-compute the "sumFE"
      sum_FE <- fixef_values[[1]][1 + dumMat[, 1]]
      for (i in 2:Q) {
        sum_FE <- sum_FE + fixef_values[[i]][1 + dumMat[, i]]
      }

      if (max(abs(sum_FE - S)) > 1e-1) {
        # divergence => we need to correct
        # we recompute the FEs

        is_ref_approx <- TRUE

        fixef_sizes <- as.integer(object$fixef_sizes)
        new_order <- order(object$fixef_sizes, decreasing = TRUE)
        fixef_sizes <- fixef_sizes[new_order]

        fe_id_list <- object$fixef_id[new_order]
        table_id_I <- as.integer(unlist(lapply(fe_id_list, table), use.names = FALSE))

        slope_flag <- rep(0L, Q)
        slope_variables <- list()

        S_demean <- cpp_demean(
          y = S, X_raw = 0, r_weights = 0, iterMax = as.integer(fixef.iter),
          diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
          fe_id_list = fe_id_list, table_id_I = table_id_I,
          slope_flag_Q = slope_flag, slope_vars_list = slope_variables,
          r_init = 0, nthreads = nthreads, save_fixef = TRUE
        )

        fixef_coef <- S_demean$fixef_coef

        end <- cumsum(fixef_sizes)
        start <- c(1, end + 1)
        for (i in 1:Q) {
          fixef_values[[new_order[i]]] <- fixef_coef[start[i]:end[i]]
        }
      }
    }

    # the information on the references
    nb_ref <- fixef_values[[Q + 1]]
    fixef_values[[Q + 1]] <- NULL
  }

  # now saving & adding information
  all_clust <- list()
  Q_all <- ifelse(isSlope, length(fixef_terms), Q)
  for (i in 1:Q_all) {
    # We put it in the right order, if requested
    fn <- attr(fixef_id[[i]], "fixef_names")

    if (sorted) {
      if (all(!grepl("[^[:digit:]]", fn))) fn <- as.numeric(fn)
      my_order <- order(fn)

      cv <- fixef_values[[i]][my_order]
      names(cv) <- fn[my_order]
      all_clust[[fixef_names[i]]] <- cv
    } else {
      cv <- fixef_values[[i]]
      names(cv) <- fn
      all_clust[[fixef_names[i]]] <- cv
    }
  }

  class(all_clust) <- c("fixest.fixef", "list")

  # Dealing with the references
  if (Q_all > 1) {
    names(nb_ref) <- fixef_names
    attr(all_clust, "references") <- nb_ref

    if (!isSlope) slope_flag <- rep(FALSE, Q)

    # warning if unbalanced
    if (notes && sum(nb_ref[!slope_flag]) > Q - 1) {
      if (is_ref_approx) {
        message("NOTE: The fixed-effects are not regular, they cannot be straightforwardly interpreted. The number of references is only approximate.")
      } else {
        message("NOTE: The fixed-effects are not regular, they cannot be straightforwardly interpreted.")
      }
    }
  }

  # Family information
  attr(all_clust, "exponential") <- FALSE
  if (object$method_type == "feNmlm" && object$family %in% c("poisson", "negbin")) {
    attr(all_clust, "exponential") <- TRUE
  } else if (object$method_type == "feglm" && object$family$link == "log") {
    attr(all_clust, "exponential") <- TRUE
  }

  return(all_clust)
}

#' Functions exported from \pkg{nlme} to implement \pkg{fixest} methods
#'
#' The package \pkg{fixest} uses the \code{fixef} method from \pkg{nlme}. Unfortunately, re-exporting this method is required in order not to attach package \pkg{nlme}.
#'
#' \itemize{
#' \item Here is the help from package \pkg{nlme}: \code{\link[nlme:fixed.effects]{fixef}}. The help from package \pkg{fixest} is here: \code{\link[fixest2]{fixef.fixest}}.
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

#' Collinearity diagnostics for \code{fixest} objects
#'
#' In some occasions, the optimization algorithm of \code{\link[fixest2]{femlm}} may fail to converge, or the variance-covariance matrix may not be available. The most common reason of why this happens is colllinearity among variables. This function helps to find out which set of variables is problematic.
#'
#'
#' @param x A \code{fixest} object obtained from, e.g. functions \code{\link[fixest2]{femlm}}, \code{\link[fixest2]{feols}} or \code{\link[fixest2]{feglm}}.
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
#' fe_1 <- sample(3, 100, TRUE)
#' fe_2 <- sample(20, 100, TRUE)
#' x <- rnorm(100, fe_1)**2
#' y <- rnorm(100, fe_2)**2
#' z <- rnorm(100, 3)**2
#' dep <- rpois(100, x * y * z)
#' base <- data.frame(fe_1, fe_2, x, y, z, dep)
#'
#' # creating collinearity problems:
#' base$v1 <- base$v2 <- base$v3 <- base$v4 <- 0
#' base$v1[base$fe_1 == 1] <- 1
#' base$v2[base$fe_1 == 2] <- 1
#' base$v3[base$fe_1 == 3] <- 1
#' base$v4[base$fe_2 == 1] <- 1
#'
#' # Estimations:
#'
#' # Collinearity with the fixed-effects:
#' res_1 <- femlm(dep ~ log(x) + v1 + v2 + v4 | fe_1 + fe_2, base)
#' collinearity(res_1)
#'
#' # => collinearity with the first fixed-effect identified, we drop v1 and v2
#' res_1bis <- femlm(dep ~ log(x) + v4 | fe_1 + fe_2, base)
#' collinearity(res_1bis)
#'
#' # Multi-Collinearity:
#' res_2 <- femlm(dep ~ log(x) + v1 + v2 + v3 + v4, base)
#' collinearity(res_2)
#'
collinearity <- function(x, verbose) {
  # x: fixest estimation

  # stop("Sorry, it does not work. A new version will hopefully come soon.")

  if (!inherits(x, "fixest")) {
    stop("Argument 'x' must be a fixest object.")
  }

  # I) (linear) collinearity with fixed-effects
  # II) (linear) multi collinearity
  # III) (non-linear) overidentification

  # flags
  isFixef <- !is.null(x$fixef_vars)
  isFE <- FALSE
  isSlope <- FALSE
  if (isFixef) {
    isSlope <- !is.null(x$fixef_terms)
    if (isSlope) {
      fixef_terms <- x$fixef_terms
      slope_flag <- grepl("\\[", fixef_terms)
      isFE <- any(!slope_flag)
    } else {
      isFE <- TRUE
    }
  }

  linear_fml <- fml_split(formula(x, "linear"), 2, split.lhs = TRUE)

  rhs_fml <- fml_split(formula(x, "linear"), 1)
  linear.varnames <- all.vars(rhs_fml[[3]])

  isLinear <- length(linear.varnames) > 0

  NL_fml <- x$NL.fml
  isNL <- !is.null(NL_fml)
  coef <- x$coefficients

  # Getting the data
  data <- fetch_data(x, "To apply function 'collinearity', ")

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  } else {
    class(data) <- "data.frame"
  }

  if (isFE) {
    linear_fml <- update(linear_fml, ~ . + 1)
  }

  # Panel setup
  panel__meta__info <- set_panel_meta_info(x, data)

  if (isLinear || isFixef || "(Intercept)" %in% names(coef)) {
    # linear.matrix = model.matrix(linear_fml, data)
    linear.matrix <- fixest_model_matrix(rhs_fml, data)
  }

  for (i in seq_along(x$obs_selection)) {
    linear.matrix <- linear.matrix[x$obs_selection[[i]], , drop = FALSE]
  }

  if (isLinear) {
    # We do that to drop interaction variables that should not be there any more
    # if factors with only NA values
    varkeep <- intersect(names(x$coefficients), colnames(linear.matrix))
    if (length(varkeep) < ncol(linear.matrix)) {
      linear.matrix <- linear.matrix[, varkeep, drop = FALSE]
    }
  }


  if (isLinear) {
    # We center and scale to have something comparable across data sets
    first_row <- linear.matrix[1, ]
    linear.matrix <- scale(linear.matrix, center = FALSE, scale = TRUE)
    # The constants => we set it to 1
    constant_id <- apply(linear.matrix, 2, anyNA)
    constant_value <- first_row[constant_id]
    linear.matrix[is.na(linear.matrix)] <- 1
    linear_mat_noIntercept <- linear.matrix[, -1, drop = FALSE]
  }

  n_obs <- nrow(data)
  Q <- length(x$fixef_id)

  # information display
  ccat <- function(...) NULL
  if (missing(verbose)) {
    verbose <- FALSE
    # verbose only if task is long
    if (100 * nrow(linear.matrix) * (ncol(linear.matrix)**2 * Q**2) >= 1e9) {
      verbose <- TRUE
    }
  }
  if (verbose) ccat <- cat

  #
  # 0) Data preparation (FE/slopes)
  #

  # Variables for FE / Slopes
  if (isSlope) {
    fixef_terms <- x$fixef_terms
    terms_full <- extract_fe_slope(fixef_terms)
    fixef_vars <- terms_full$fixef_vars
    slope_fe <- terms_full$slope_fe
    fe_all <- terms_full$fe_all
    slope_vars <- terms_full$slope_vars
    slope_terms <- terms_full$slope_terms

    # dictionary mapping fe var names to the ids of id_dummies_vect
    dict_fe <- 1:Q
    names(dict_fe) <- x$fixef_vars
    slope_vars_unik <- unique(slope_vars)

    if (any(!slope_vars_unik %in% names(data))) {
      var_pblm <- setdiff(slope_vars_unik, names(data))
      stop("To check collinearity, we need to fetch some variables in the original dataset (", deparse_long(x$call$data), "). However, the variable", enumerate_items(var_pblm, "s.is"), " not present in the original dataset any more.")
    }

    slope_var_list <- list()
    for (i in 1:length(slope_vars_unik)) {
      variable <- all.vars(str2lang(slope_vars_unik[i]))

      # as.numeric => we'll use cpp so required
      svar <- as.numeric(as.vector(eval(str2lang(slope_vars_unik[i]), data)))
      if (length(svar) == 1) svar <- rep(svar, n_obs)

      slope_var_list[[slope_vars_unik[i]]] <- svar
    }

    Q_slope <- length(slope_fe)
    slope_fe_id <- x$fixef_id[slope_fe]
    slope_var_all <- slope_var_list[slope_vars]

    Q_fe <- length(fixef_vars)
    fe_id <- x$fixef_id[fixef_vars]
  } else if (isFE) {
    Q_fe <- length(x$fixef_id)
    fe_id <- x$fixef_id
  }


  #
  # I) collinearity with clusters
  #

  ccat("Checking Collinearity: ")

  #
  # Variables are constant or equal to 0
  #

  if (isLinear) {
    if (isFE) {
      if (any(constant_id[-1])) {
        # var_problem = colnames(linear_mat_noIntercept)[is_const]
        var_problem <- colnames(linear_mat_noIntercept)[constant_id[-1]]
        message <- paste0("Variable", enumerate_items(var_problem, "s.is.quote"), " constant, thus collinear with the fixed-effects.")
        print(message)
        return(invisible(message))
      }
    } else {
      isIntercept <- attr(terms(x$fml), "intercept")

      if (isIntercept && any(constant_id[-1])) {
        var_problem <- colnames(linear_mat_noIntercept)[constant_id[-1]]
        message <- paste0("Variable", enumerate_items(var_problem, "s.is.quote"), " constant, thus collinear with the intercept.")
        print(message)
        return(invisible(message))
      }

      if (any(constant_value == 0)) {
        # constant and equal to 0
        var_problem <- colnames(linear.matrix)[first_row == 0 & constant_id]
        message <- paste0("Variable", enumerate_items(var_problem, "s.is.quote"), " constant and equal to 0.")

        print(message)
        return(invisible(message))
      }
    }
  }

  if (isFE && isLinear) {
    ccat("simple with fixed-effects:")
    # We project each variable onto the cluster subspace

    cluster <- fe_id
    Q_fe <- length(cluster)
    for (q in 1:Q_fe) {
      ccat(".")
      dum <- cluster[[q]]
      k <- max(dum)

      value <- cpp_tapply_sum(Q = k, x = linear_mat_noIntercept, dum = dum)
      nb_per_cluster <- cpp_table(Q = k, dum = dum)

      # residuals of the linear projection on the cluster space
      residuals <- linear_mat_noIntercept - (value / nb_per_cluster)[dum, ]

      max_residuals <- apply(abs(residuals), 2, max)

      if (any(max_residuals < 1e-6)) {
        ccat("\n")
        varnames <- colnames(linear_mat_noIntercept)
        collin_var <- varnames[max_residuals < 1e-6]
        message <- paste0("Variable", enumerate_items(collin_var, "s.is.quote"), " collinear with fixed-effects '", names(cluster)[q], "'.")

        print(message)
        return(invisible(message))
      }
    }
    ccat("OK")
  }

  if (isSlope && isLinear) {
    ccat("simple with variables with varying slopes:")

    for (q in 1:Q_slope) {
      ccat(".")
      dum <- fe_id[[q]]
      my_var <- slope_var_all[[q]]
      k <- max(dum)
      value <- cpp_tapply_sum(Q = k, x = linear_mat_noIntercept * my_var, dum = dum)

      denom <- as.vector(tapply(my_var**2, dum, sum))

      # residuals of the linear projection on the cluster space
      residuals <- linear_mat_noIntercept - (value / denom)[dum, ] * my_var

      max_residuals <- apply(abs(residuals), 2, max)

      if (any(max_residuals < 1e-6)) {
        ccat("\n")
        varnames <- colnames(linear_mat_noIntercept)
        collin_var <- varnames[max_residuals < 1e-6]

        message <- paste0("Variable", enumerate_items(collin_var, "s.is.quote"), " collinear with variable with varying slope '", slope_vars[q], "' (on '", slope_fe[q], "').")

        print(message)
        return(invisible(message))
      }
    }
    ccat("OK")
  }

  #
  # II) perfect multicollinearity
  #

  name2change <- grepl("[^[:alnum:]\\._]", colnames(linear.matrix))
  dict_name <- colnames(linear.matrix)[!grepl("(Intercept)", colnames(linear.matrix))]
  if (any(name2change)) {
    linbase <- as.data.frame(linear.matrix)
    names(linbase)[name2change] <- gsub("[^[:alnum:]\\._]", "_", names(linbase)[name2change])
  } else {
    linbase <- as.data.frame(linear.matrix)
  }
  linearVars <- names(linbase)[!grepl("Intercept", names(linbase))]
  names(dict_name) <- linearVars
  dict_name["_Intercept_"] <- "(Intercept)"

  # for multicol without cluster
  mat_base <- as.matrix(linbase)

  # we add possibly missing variables
  varmiss <- setdiff(all.vars(linear_fml), names(linbase))
  for (v in varmiss) linbase[[v]] <- data[[v]]

  if (isLinear && length(linearVars) >= 2) {
    ccat(ifelse(Q >= 1, ", ", " "), "multiple:")
    for (v in linearVars) {
      ccat(".")

      i <- which(colnames(mat_base) == v)
      res <- ols_fit(y = mat_base[, i], X = mat_base[, -i, drop = FALSE], w = 1, collin.tol = 1e-10, nthreads = 1)

      max_residuals <- max(abs(res$residuals))

      if (max_residuals < 1e-4) {
        ccat("\n")
        coef_lm <- res$coefficients
        vars <- colnames(mat_base)[-i]
        collin_var <- vars[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
        message <- paste0("Variable '", dict_name[v], "' is collinear with variable", enumerate_items(dict_name[collin_var], "s.quote"), ".")

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
  if (isLinear && length(linearVars) >= 2 && isFixef) {
    ccat(", multiple with cluster")

    dum_names <- names(x$fixef_id)
    n_clust <- length(dum_names)
    new_dum_names <- paste0("__DUM_", 1:n_clust)
    for (i in 1:n_clust) {
      linbase[[paste0("__DUM_", i)]] <- x$fixef_id[[i]]
    }

    for (v in linearVars) {
      ccat(".")
      fml2estimate <- as.formula(paste0(v, "~", paste0(setdiff(linearVars, v), collapse = "+")))

      for (id_cluster in 1:n_clust) {
        res <- feols(fml2estimate, linbase, fixef = new_dum_names[id_cluster], warn = FALSE)

        max_residuals <- max(abs(res$residuals))

        if (max_residuals < 1e-4) {
          ccat("\n")
          coef_lm <- coef(res)
          collin_var <- names(coef_lm)[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
          message <- paste0("Variable '", dict_name[v], "' is collinear with variable", enumerate_items(dict_name[collin_var], "s.quote"), ", together with the fixed-effects ", dum_names[id_cluster], ".")

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

  if (isNL) {
    NL_vars <- all.vars(NL_fml)
    varNotHere <- setdiff(NL_vars, c(names(coef), names(data)))
    if (length(varNotHere) > 0) {
      stop("Some variables used to estimate the model (in the non-linear formula) are missing from the original data: ", paste0(varNotHere, collapse = ", "), ".")
    }

    var2send <- intersect(NL_vars, names(data))
    env <- new.env()
    for (var in var2send) {
      assign(var, data[[var]], env)
    }
  }

  if (isNL && (length(coef) >= 2 || isFixef)) {
    ccat(", in non-linear term:")
    if (isFixef) {
      # we add the constant
      coef["CONSTANT"] <- 1
      data$CONSTANT <- 1
      linear.matrix <- model.matrix(update(linear_fml, ~ . - 1 + CONSTANT), data)
    }

    coef_name <- names(coef)

    NL_coef <- setdiff(all.vars(NL_fml), names(data))
    # L_coef = setdiff(coef_name, NL_coef)
    L_coef <- colnames(linear.matrix)

    #
    # We compute mu:
    #

    mu <- 0
    if (length(L_coef) > 0) {
      mu <- linear.matrix %*% coef[L_coef]
    }

    for (iter_coef in NL_coef) {
      assign(iter_coef, coef[iter_coef], env)
    }

    # Evaluation of the NL part
    value_NL <- eval(NL_fml[[2]], env)
    mu <- mu + value_NL
    data$mu <- mu

    if (var(mu) == 0) {
      message <- "Variance of the NL part is 0."
      print(message)
      return(invisible(message))
    }

    #
    # The loop
    #

    for (var in coef_name) {
      ccat(".")
      # we modify the coef of var, then we fix it:
      # if we can find a set of other parameters that give the same fit,
      # then there is an identification issue

      if (var %in% L_coef) {
        # if linear coef: we modify the fml and add an offset
        if (length(L_coef) == 1) {
          fml <- mu ~ 0
        } else {
          fml <- as.formula(paste0("mu ~ 0+", paste0(setdiff(L_coef, var), collapse = "+")))
        }

        # offset with the new value
        offset <- as.formula(paste0("~", coef[var] + 1, "*", var))

        res <- feNmlm(fml, data = data, family = "gaussian", NL.fml = NL_fml, NL.start = as.list(coef), offset = offset)
      } else {
        # NL case:
        # we modify both fml and NL.fml

        if (isFixef) {
          fml <- update(x$fml, mu ~ . - 1 + CONSTANT)
        } else {
          fml <- update(x$fml, mu ~ .)
        }

        # the new formula with the new variable
        NL_fml_char <- as.character(NL_fml)[2]
        NL_fml_new <- as.formula(paste0("~", gsub(paste0("\\b", var, "\\b"), 1 + coef[var], NL_fml_char)))

        if (length(NL_coef) == 1) {
          # there is no more parameter to estimate => offset
          res <- femlm(fml, data = data, family = "gaussian", offset = NL_fml_new)
        } else {
          res <- feNmlm(fml, data = data, family = "gaussian", NL.fml = NL_fml_new, NL.start = as.list(coef))
        }
      }

      sum_resids <- sum(abs(res$residuals))
      if (sum_resids < 1e-4) {
        coef_esti <- coef(res)
        coef_diff <- abs(coef_esti - coef[names(coef_esti)])
        collin_var <- names(coef_diff)[coef_diff > 1e-3]
        message <- paste0("Coefficients ", show_vars_limited_width(c(var, collin_var)), " are not uniquely identifed.")

        print(message)
        return(invisible(message))
      }
    }
    ccat("OK")
  }

  ccat("\n")
  message <- "No visible collinearity problem. (Doesn't mean there's none!)"
  print(message)
  return(invisible(message))
}

.meanDiffTable <- function(mat_vars, treat_var, treat_cases, indiv_var, raw, diff.inv = FALSE) {
  # This function is the workhorse of did_means

  x_name <- colnames(mat_vars)

  treat_01 <- 1 * (treat_var == treat_cases[2])

  nthreads <- getFixest_nthreads()
  info <- cpppar_cond_means(mat_vars, treat_01, nthreads)

  means <- info$means
  sds <- info$sd
  ns <- info$n
  n_01 <- info$n_01

  name_1 <- as.character(treat_cases[1])
  name_2 <- as.character(treat_cases[2])

  means_1 <- means[, 1]
  means_2 <- means[, 2]

  sds_1 <- sds[, 1]
  sds_2 <- sds[, 2]

  n_1 <- ns[, 1]
  n_2 <- ns[, 2]

  format_1 <- paste0(mysignif(means_1, 2), " (", mysignif(sds_1, 2), ")")
  format_2 <- paste0(mysignif(means_2, 2), " (", mysignif(sds_2, 2), ")")

  if (diff.inv) {
    D <- (means_2 - means_1)
  } else {
    D <- (means_1 - means_2)
  }

  sd_D <- sqrt(sds_1**2 / n_1 + sds_2**2 / n_2)

  t_stat <- signif(D / sd_D, 3)

  #
  # Formatting
  #

  if (raw) {
    res <- data.frame(vars = x_name, stringsAsFactors = FALSE)
    res[[paste0("Mean: ", treat_cases[1])]] <- c(means_1, recursive = TRUE)
    res[[paste0("Mean: ", treat_cases[2])]] <- c(means_2, recursive = TRUE)
    res[[paste0("se: ", treat_cases[1])]] <- c(sds_1, recursive = TRUE)
    res[[paste0("se: ", treat_cases[2])]] <- c(sds_2, recursive = TRUE)
    res[["Difference"]] <- c(D, recursive = TRUE)
    res[["t-stat"]] <- c(as.character(t_stat), recursive = TRUE)
  } else {
    res <- data.frame(vars = x_name, stringsAsFactors = FALSE)
    res[[paste0("cond: ", treat_cases[1])]] <- format_1
    res[[paste0("cond: ", treat_cases[2])]] <- format_2
    res[["Difference"]] <- c(mysignif(D, 3), recursive = TRUE)
    res[["t-stat"]] <- c(as.character(t_stat), recursive = TRUE)

    res[nrow(res) + 1, ] <- c("Observations", addCommas(n_01[1]), addCommas(n_01[2]), "", "")
    # Le nombre d'individus
    if (!missing(indiv_var) && !is.null(indiv_var)) {
      nb_indiv <- tapply(indiv_var, treat_var, function(x) length(unique(x)))
      res[nrow(res) + 1, ] <- c("# Individuals", addCommas(nb_indiv[name_1]), addCommas(nb_indiv[name_2]), "", "")
    }
  }

  attr(res, "na") <- info$na

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
#' @param bin2 A list or vector defining the binning of the second variable. See help for the argument \code{bin} for details (or look at the help of the function \code{\link[fixest2]{bin}}). You can use \code{.()} for \code{list()}.
#' @param ... Not currently used.
#'
#' @details
#' To interact fixed-effects, this function should not be used: instead use directly the syntax \code{fe1^fe2} in the fixed-effects part of the formula. Please see the details and examples in the help page of \code{\link[fixest2]{feols}}.
#'
#' @return
#' It returns a matrix with number of rows the length of \code{factor_var}. If there is no interacted variable or it is interacted with a numeric variable, the number of columns is equal to the number of cases contained in \code{factor_var} minus the reference(s). If the interacted variable is a factor, the number of columns is the number of combined cases between \code{factor_var} and \code{var}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest2:coefplot]{iplot}} to plot interactions or factors created with \code{i()}, \code{\link[fixest2]{feols}} for OLS estimation with multiple fixed-effects.
#'
#' See the function \code{\link[fixest2]{bin}} for binning variables.
#'
#' @examples
#'
#' #
#' # Simple illustration
#' #
#'
#' x <- rep(letters[1:4], 3)[1:10]
#' y <- rep(1:4, c(1, 2, 3, 4))
#'
#' # interaction
#' data.frame(x, y, i(x, y, ref = TRUE))
#'
#' # without interaction
#' data.frame(x, i(x, "b"))
#'
#' # you can interact factors too
#' z <- rep(c("e", "f", "g"), c(5, 3, 2))
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
#' est_did <- feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)
#'
#' # => plot only interactions with iplot
#' iplot(est_did)
#'
#' # Using i() for factors
#' est_bis <- feols(y ~ x1 + i(period, keep = 3:6) + i(period, treat, 5) | id, base_did)
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
#' aq <- airquality
#' aq$week <- aq$Day %/% 7 + 1
#'
#' # Interacting Month and week:
#' res_2F <- feols(Ozone ~ Solar.R + i(Month, i.week), aq)
#'
#' # Same but dropping the 5th Month and 1st week
#' res_2F_bis <- feols(Ozone ~ Solar.R + i(Month, i.week, ref = 5, ref2 = 1), aq)
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
i <- function(factor_var, var, ref, keep, bin, ref2, keep2, bin2, ...) {
  # Used to create interactions

  # Later: binning (bin = 1:3 // bin = list("a" = "[abc]")). Default name is bin name (eg "1:3")

  # gt = function(x) cat(sfill(x, 20), ": ", -(t0 - (t0<<-proc.time()))[3], "s\n", sep = "")
  # t0 = proc.time()

  validate_dots(valid_args = c("f2", "f_name", "ref_special", "sparse"))

  dots <- list(...)
  is_sparse <- isTRUE(dots$sparse)

  mc <- match.call()

  # Finding if it's a call from fixest
  FROM_FIXEST <- is_fixest_call()

  # General checks
  check_arg(factor_var, "mbt vector")

  # NOTA:
  # the user can use the prefix "i." to tell the algorithm to consider the
  # variable var as a factor. This requires a non standard evaluation
  #

  var_name <- NULL
  if (!is.null(mc$f2)) {
    # should be only used internally

    var_name <- deparse_long(mc$f2)
    var <- dots$f2
    check_value(var, "vector", .arg_name = "f2")
  }

  # General information on the factor variable
  if (!is.null(dots$f_name)) {
    f_name <- dots$f_name
  } else {
    f_name <- deparse_long(mc$factor_var)
  }
  f <- factor_var # renaming for clarity

  # checking var
  IS_INTER_NUMERIC <- IS_INTER_FACTOR <- FALSE
  if (!missing(var)) {
    # Evaluation of var (with possibly, user prepending with i.)
    if (is.null(var_name)) {
      var_name <- deparse_long(mc$var)
      if (grepl("^i\\.", var_name)) {
        var_name <- gsub("^i\\.", "", var_name)
        var <- str2lang(var_name)
        check_value_plus(var, "evalset vector", .data = parent.frame(), .arg_name = "var")
        IS_INTER_FACTOR <- TRUE
      } else {
        check_value(var, "vector", .arg_name = "var")
      }
    } else {
      # f2 was used
      IS_INTER_FACTOR <- TRUE
    }

    # is it an interaction with a factor or with a continuous variable?
    if (length(var) == length(f)) {
      # Logical => numeric by default
      # otherwise, just setting the flags

      if (IS_INTER_FACTOR == FALSE) {
        # If previous condition == TRUE, means was requested by the user

        if (is.logical(var)) {
          var <- as.numeric(var)
          IS_INTER_NUMERIC <- TRUE
        } else if (is.numeric(var)) {
          IS_INTER_NUMERIC <- TRUE
        } else {
          IS_INTER_FACTOR <- TRUE
        }
      }
    } else if (length(var) <= 2 && !"ref" %in% names(mc)) {
      # ref is implicitly called via the location of var
      ref <- var
      dots$ref_special <- FALSE
    } else {
      if (grepl("^[[:alpha:]\\.][[:alnum:]\\._]*:[[:alpha:]\\.][[:alnum:]\\._]*$", var_name)) {
        info <- strsplit(var_name, ":")[[1]]
        stop("In i(): When 'var' is equal to a product, please use I(", info[1], "*", info[2], ") instead of ", var_name, ".")
      } else {
        stop("The arguments 'var' and 'f' must be of the same length (currently ", length(var), " vs ", length(f), ").")
      }
    }
  }

  IS_INTER <- IS_INTER_NUMERIC || IS_INTER_FACTOR

  if (isTRUE(dots$ref_special) && (!IS_INTER || IS_INTER_FACTOR)) {
    ref <- TRUE
  }

  #
  # QUFing + NA
  #

  if (IS_INTER) {
    is_na_all <- is.na(f) | is.na(var)
  } else {
    is_na_all <- is.na(f)
  }

  if (!missing(bin)) {
    bin <- error_sender(eval_dot(bin), arg_name = "bin")

    if (!is.null(bin)) {
      f <- bin_factor(bin, f, f_name)
    }
  }

  if (IS_INTER_FACTOR && !MISSNULL(bin2)) {
    bin2 <- error_sender(eval_dot(bin2), arg_name = "bin2")

    if (!is.null(bin2)) {
      var <- bin_factor(bin2, var, f_name)
    }
  }

  if (IS_INTER_FACTOR) {
    info <- to_integer(f, var, add_items = TRUE, items.list = TRUE, sorted = TRUE, multi.join = "__%%__")
  } else {
    info <- to_integer(f, add_items = TRUE, items.list = TRUE, sorted = TRUE)
  }

  fe_num <- info$x
  items <- info$items

  if (!IS_INTER_NUMERIC) {
    # neutral var in C code
    var <- 1
  }

  if (IS_INTER_FACTOR) {
    items_split <- strsplit(items, "__%%__", fixed = TRUE)

    f_items <- sapply(items_split, `[`, 1)
    var_items <- sapply(items_split, `[`, 2)
  } else {
    f_items <- items
  }

  check_arg(ref, "logical scalar | vector no na")

  check_arg(ref2, keep, keep2, "vector no na")

  no_rm <- TRUE
  id_drop <- c()
  if (!missing(ref)) {
    if (isTRUE(ref)) {
      # We always delete the first value
      # Que ce soit items ici est normal (et pas f_items)
      id_drop <- which(items == items[1])
    } else {
      id_drop <- items_to_drop(f_items, ref, "factor_var")
    }
    ref_id <- id_drop
  }


  if (!missing(keep)) {
    id_drop <- c(id_drop, items_to_drop(f_items, keep, "factor_var", keep = TRUE))
  }

  if (IS_INTER_FACTOR) {
    if (!missing(ref2)) {
      id_drop <- c(id_drop, items_to_drop(var_items, ref2, "var"))
    }

    if (!missing(keep2)) {
      id_drop <- c(id_drop, items_to_drop(var_items, keep2, "var", keep = TRUE))
    }
  }

  if (length(id_drop) > 0) {
    id_drop <- unique(sort(id_drop))
    if (length(id_drop) == length(items)) stop("All items from the interaction have been removed.")
    who_is_dropped <- id_drop
    no_rm <- FALSE
  } else {
    # -1 is neutral
    who_is_dropped <- -1
  }

  # The column names

  if (length(id_drop) > 0) {
    items_name <- items[-id_drop]
    f_items <- f_items[-id_drop]
    if (IS_INTER_FACTOR) {
      var_items <- var_items[-id_drop]
    }
  } else {
    items_name <- items
  }

  if (FROM_FIXEST) {
    # Pour avoir des jolis noms c'est un vrai gloubiboulga,
    # mais j'ai pas trouve plus simple...
    if (IS_INTER_FACTOR) {
      col_names <- paste0("__CLEAN__", f_name, "::", f_items, ":", var_name, "::", var_items)
      col_names_full <- NULL
    } else if (IS_INTER_NUMERIC) {
      col_names <- paste0("__CLEAN__", f_name, "::", f_items, ":", var_name)
      col_names_full <- paste0(f_name, "::", items, ":", var_name)
    } else {
      col_names <- paste0("__CLEAN__", f_name, "::", items_name)
      col_names_full <- paste0(f_name, "::", items)
    }
  } else {
    if (IS_INTER_FACTOR) {
      items_name <- gsub("__%%__", ":", items_name, fixed = TRUE)
    }

    col_names <- items_name
  }

  if (is_sparse) {
    # Internal call: we return the row ids + the values + the indexes + the names

    if (length(who_is_dropped) > 0) {
      valid_row <- !is_na_all & !fe_num %in% who_is_dropped
    } else {
      valid_row <- !is_na_all
    }

    # we need to ensure the IDs go from 1 to # Unique
    fe_colid <- to_integer(fe_num[valid_row], sorted = TRUE)

    values <- if (length(var) == 1) rep(1, length(valid_row)) else var
    res <- list(
      rowid = which(valid_row), values = values,
      colid = fe_colid, col_names = col_names
    )

    class(res) <- "i_sparse"

    return(res)
  }

  res <- cpp_factor_matrix(fe_num, is_na_all, who_is_dropped, var, col_names)
  # res => matrix with...
  #  - NAs where appropriate
  #  - appropriate number of columns
  #  - interacted if needed
  #


  # We send the information on the reference
  if (FROM_FIXEST) {
    is_GLOBAL <- FALSE
    for (where in 1:min(8, sys.nframe())) {
      if (exists("GLOBAL_fixest_mm_info", parent.frame(where))) {
        GLOBAL_fixest_mm_info <- get("GLOBAL_fixest_mm_info", parent.frame(where))
        is_GLOBAL <- TRUE
        break
      }
    }

    if (is_GLOBAL && !IS_INTER_FACTOR) {
      # NOTA:
      # if you change stuff here => change them also in sunab()

      info <- list()
      info$coef_names_full <- col_names_full
      info$items <- items
      if (!missing(ref)) {
        info$ref_id <- ref_id
        info$ref <- items[ref_id]
      }
      info$is_num <- is.numeric(items)
      info$f_name <- f_name
      info$is_inter_num <- IS_INTER_NUMERIC
      if (IS_INTER_NUMERIC) {
        info$var_name <- var_name
      }
      info$is_inter_fact <- IS_INTER_FACTOR

      GLOBAL_fixest_mm_info[[length(GLOBAL_fixest_mm_info) + 1]] <- info
      # re assignment
      assign("GLOBAL_fixest_mm_info", GLOBAL_fixest_mm_info, parent.frame(where))
    }
  }

  res
}


i_ref <- function(factor_var, var, ref, bin, keep, ref2, keep2, bin2) {
  # To automatically add references when i(x) is used

  mc <- match.call()

  mc[[1]] <- as.name("i")

  if (!any(c("ref", "keep", "ref2", "keep2") %in% names(mc))) {
    mc$ref_special <- TRUE
  }

  return(deparse_long(mc))
}

i_noref <- function(factor_var, var, ref, bin, keep, ref2, keep2, bin2) {
  # Used only in predict => to create data without restriction

  mc <- match.call()

  mc[[1]] <- as.name("i")

  mc$ref <- mc$keep <- mc$ref2 <- mc$keep2 <- NULL

  return(deparse_long(mc))
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
#' The DSB can even be used within variable names, but then the variable must be nested in character form. For example \code{y ~ .["x.[1:2]_sq"]} will create \code{y ~ x1_sq + x2_sq}. Using the character form is important to avoid a formula parsing error. Double quotes must be used. Note that the character string that is nested will be parsed with the function \code{\link[fixest2]{dsb}}, and thus it will return a vector.
#'
#' By default, the DSB operator expands vectors into sums. You can add a comma, like in \code{.[, x]}, to expand with commas--the content can then be used within functions. For instance: \code{c(x.[, 1:2])} will create \code{c(x1, x2)} (and \emph{not} \code{c(x1 + x2)}).
#'
#' In all \code{fixest} estimations, this special parsing is enabled, so you don't need to use \code{xpd}.
#'
#' You can even use multiple square brackets within a single variable, but then the use of nesting is required. For example, the following \code{xpd(y ~ .[".[letters[1:2]]_.[1:2]"])} will create \code{y ~ a_1 + b_2}. Remember that the nested character string is parsed with \code{\link[fixest2]{dsb}}, which explains this behavior.
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
#' \code{\link[fixest2]{setFixest_fml}} to set formula macros, and \code{\link[fixest2]{dsb}} to modify character strings with the DSB operator.
#'
#' @examples
#'
#' # Small examples with airquality data
#' data(airquality)
#' # we set two macro variables
#' setFixest_fml(
#'   ..ctrl = ~ Temp + Day,
#'   ..ctrl_long = ~ poly(Temp, 2) + poly(Day, 2)
#' )
#'
#' # Using the macro in lm with xpd:
#' lm(xpd(Ozone ~ Wind + ..ctrl), airquality)
#' lm(xpd(Ozone ~ Wind + ..ctrl_long), airquality)
#'
#' # You can use the macros without xpd() in fixest estimations
#' a <- feols(Ozone ~ Wind + ..ctrl, airquality)
#' b <- feols(Ozone ~ Wind + ..ctrl_long, airquality)
#' etable(a, b, keep = "Int|Win")
#'
#'
#' # Using .[]
#'
#' base <- setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' i <- 2:3
#' z <- "species"
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
#' base <- iris
#' names(base) <- c("y", "x1", "x2", "x3", "species")
#'
#' # We first create a matrix with all possible combinations of variables
#' my_args <- lapply(names(base)[-(1:2)], function(x) c("", x))
#' (all_combs <- as.matrix(do.call("expand.grid", my_args)))
#'
#' res_all <- list()
#' for (i in 1:nrow(all_combs)) {
#'   res_all[[i]] <- feols(xpd(y ~ x1 + ..v, ..v = all_combs[i, ]), base)
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
#' feols(Armed.Forces ~ Population + sw(regex(, "GNP|ployed")), longley)
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
#' vars <- letters[1:5]
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
#' var <- "a"
#' xpd(y ~ x.[var])
#'
#' # ... the variables can be multiple
#' vars <- LETTERS[1:3]
#' xpd(y ~ x.[vars])
#'
#' # You can have "complex" variable names but they must be nested in character form
#' xpd(y ~ .["x.[vars]_sq"])
#'
#' # DSB can be used within regular expressions
#' re <- c("GNP", "Pop")
#' xpd(Unemployed ~ regex(".[re]"), data = longley)
#'
#' # => equivalent to regex("GNP|Pop")
#'
#' # Use .[,var] (NOTE THE COMMA!) to expand with commas
#' # !! can break the formula if missused
#' vars <- c("wage", "unemp")
#' xpd(c(y.[, 1:3]) ~ csw(.[, vars]))
#'
#'
#' # Example of use of .[] within a loop
#' res_all <- list()
#' for (p in 1:3) {
#'   res_all[[p]] <- feols(Ozone ~ Wind + poly(Temp, .[p]), airquality)
#' }
#'
#' etable(res_all)
#'
#' # The former can be compactly estimated with:
#' res_compact <- feols(Ozone ~ Wind + sw(.[, "poly(Temp, .[1:3])"]), airquality)
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
xpd <- function(fml, ..., lhs, rhs, data = NULL) {
  .xpd(fml = fml, ..., lhs = lhs, rhs = rhs, data = data, check = TRUE, macro = TRUE, frame = parent.frame())
}

.xpd <- function(fml, ..., lhs, rhs, data = NULL, check = FALSE, macro = FALSE, frame = NULL) {
  is_lhs <- !missing(lhs)
  is_rhs <- !missing(rhs)
  if (is_lhs || is_rhs) {
    if (check) check_arg(fml, .type = "formula", .up = 1)

    # Direct formula creation
    if (missing(fml)) {
      res <- if (is_lhs) 1 ~ 1 else ~1
    } else {
      if (is_lhs && length(fml)[1] == 2) {
        # Means: RHS formula and we add the LHS
        res <- 1 ~ 1
        res[[3]] <- fml[[2]]
      } else {
        res <- fml
      }
    }

    if (is_lhs) {
      res[[2]] <- value2stringCall(lhs, call = TRUE, check = check)
    }

    if (is_rhs) {
      res[[length(res)]] <- value2stringCall(rhs, call = TRUE, check = check)
    }

    fml <- res

    if (!macro) {
      return(fml)
    }

    # NOTA:
    # if we allow for macro implementation ex post:
    # This entails a 50% performance drop in terms of speed.
    # Now, without macro variables, speed is at 30us while it was 20us before
    # so in .xpd => macro argument
  } else if (check) {
    check_arg(fml, .type = "formula mbt", .up = 1)
  }

  macros <- parse_macros(..., from_xpd = TRUE, check = check, frame = frame)

  is_macro <- length(macros) != 0
  is_data <- !missnull(data)
  fml_funs <- all.vars(fml, functions = TRUE)
  is_brackets <- "[" %in% fml_funs
  is_regex <- any(c("regex", "..") %in% fml_funs)
  if (!(is_macro || is_data || is_brackets || is_regex)) {
    return(fml)
  }

  fml_dp <- NULL

  if (is_macro) {
    # We allow recursion + auto-blocking at 5th depth
    # to allow crazy things from the user side (but not too much)
    max_depth <- 5
    depth <- 0
    any_done <- TRUE
    qui <- which(names(macros) %in% all.vars(fml))
    while (length(qui) > 0 && any_done && depth < max_depth) {
      depth <- depth + 1
      fml_dp_origin <- fml_dp <- deparse_long(fml)
      # We need to add a lookafter assertion: otherwise if we have ..ctrl + ..ctrl_long, there will be a replacement in ..ctrl_long
      for (i in qui) {
        fml_dp <- gsub(paste0(escape_regex(names(macros)[i]), "(?=$|[^[:alnum:]_\\.])"), macros[[i]], fml_dp, perl = TRUE)
      }

      fml <- as.formula(fml_dp, frame)

      qui <- which(names(macros) %in% all.vars(fml))
      any_done <- fml_dp_origin != fml_dp
    }

    if (depth == max_depth && length(qui) > 0) {
      warning(dsb(
        "In xpd, max recursivity has been reached. ",
        "Please revise the recursive definition of your variables. ",
        "It concerns the variable.[*s_, q, C?names(macros)[qui]]."
      ))
    }
  }

  if (is_regex) {
    # We expand only if data is provided (it means the user wants us to check)
    # if .[]: we expand inside the ..(".[1:3]") => ..("1|2|3")

    if (is.null(fml_dp)) fml_dp <- deparse_long(fml)
    is_brackets <- grepl(".[", fml_dp, fixed = TRUE)

    if (is_data || is_brackets) {
      if (is.null(fml_dp)) fml_dp <- deparse_long(fml)

      if (is_data) {
        check_arg(data, "character vector no na | matrix | data.frame")

        if (is.matrix(data)) {
          data <- colnames(data)
          if (is.null(data)) {
            stop("The argument 'data' must contain variable names. It is currently a matrix without column names.")
          }
        } else if (is.data.frame(data)) {
          data <- names(data)
        }
      }

      fml_dp_split <- strsplit(fml_dp, '(?!<[[:alnum:]._])(regex|\\.\\.)\\((?=[,"])',
        perl = TRUE
      )[[1]]

      res <- fml_dp_split
      for (i in 2:length(res)) {
        re <- sub('"\\).*', "", res[i])

        re_width <- nchar(re)

        is_comma <- grepl("^,", re)
        re <- gsub("^,? ?\"", "", re)

        re <- dot_square_bracket(re, frame, regex = TRUE)

        if (is_data) {
          vars <- grep(re, data, value = TRUE, perl = TRUE)
          if (length(vars) == 0) {
            vars <- "1"
          }

          coll <- if (is_comma) ", " else " + "

          res[i] <- paste0(paste(vars, collapse = coll), substr(res[i], re_width + 3, nchar(res[i])))
        } else {
          res[i] <- paste0("..(\"", re, "\")", substr(res[i], re_width + 3, nchar(res[i])))
        }
      }

      fml_txt <- paste(res, collapse = "")
      fml <- error_sender(as.formula(fml_txt, frame),
        "Expansion of variables in ..(\"regex\"): coercion of the following text to a formula led to an error.\n",
        fit_screen(paste0("     TEXT: ", fml_txt)),
        "\n  PROBLEM: see below",
        up = 1
      )
    }
  }

  if ("[" %in% all.vars(fml, functions = TRUE)) {
    fml_txt <- deparse_long(fml)
    if (grepl(".[", fml_txt, fixed = TRUE)) {
      fml_txt <- dot_square_bracket(fml_txt, frame)

      # level 1 recursivity
      if (grepl(".[", fml_txt, fixed = TRUE)) {
        fml_txt <- dot_square_bracket(fml_txt, frame)
      }

      fml <- error_sender(as.formula(fml_txt, frame),
        "Dot square bracket operator: coercion of the following text to a formula led to an error.\n",
        fit_screen(paste0("     TEXT: ", fml_txt)),
        "\n  PROBLEM: see below",
        up = 1
      )
    }
  }

  fml
}

#' Centers a set of variables around a set of factors
#'
#' User-level access to internal demeaning algorithm of \code{fixest}.
#'
#' @inheritSection feols Varying slopes
#'
#' @param X A matrix, vector, data.frame or a list OR a formula OR a \code{\link[fixest2]{feols}} estimation. If equal to a formula, then the argument \code{data} is required, and it must be of the type: \code{x1 + x2 ~ f1 + fe2} with on the LHS the variables to be centered, and on the RHS the factors used for centering. Note that you can use variables with varying slopes with the syntax \code{fe[v1, v2]} (see details in \code{\link[fixest2]{feols}}). If a \code{feols} estimation, all variables (LHS+RHS) are demeaned and then returned (only if it was estimated with fixed-effects). Otherwise, it must represent the data to be centered. Of course the number of observations of that data must be the same as the factors used for centering (argument \code{f}).
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
#' base <- trade
#' base$ln_dist <- log(base$dist_km)
#' base$ln_euros <- log(base$Euros)
#'
#' # We center the two variables ln_dist and ln_euros
#' #  on the factors Origin and Destination
#' X_demean <- demean(
#'   X = base[, c("ln_dist", "ln_euros")],
#'   f = base[, c("Origin", "Destination")]
#' )
#' base[, c("ln_dist_dm", "ln_euros_dm")] <- X_demean
#'
#' est <- feols(ln_euros_dm ~ ln_dist_dm, base)
#' est_fe <- feols(ln_euros ~ ln_dist | Origin + Destination, base)
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
#' base <- iris
#' names(base) <- c("y", "x1", "x2", "x3", "species")
#'
#' #
#' # We center y and x1 on species and x2 * species
#'
#' # using a formula
#' base_dm <- demean(y + x1 ~ species[x2], data = base)
#'
#' # using vectors
#' base_dm_bis <- demean(
#'   X = base[, c("y", "x1")], f = base$species,
#'   slope.vars = base$x2, slope.flag = 1
#' )
#'
#' # Let's look at the equivalences
#' res_vs_1 <- feols(y ~ x1 + species + x2:species, base)
#' res_vs_2 <- feols(y ~ x1, base_dm)
#' res_vs_3 <- feols(y ~ x1, base_dm_bis)
#'
#' # only the small sample adj. differ in the SEs
#' etable(res_vs_1, res_vs_2, res_vs_3, keep = "x1")
#'
#' #
#' # center on x2 * species and on another FE
#'
#' base$fe <- rep(1:5, 10)
#'
#' # using a formula => double square brackets!
#' base_dm <- demean(y + x1 ~ fe + species[[x2]], data = base)
#'
#' # using vectors => note slope.flag!
#' base_dm_bis <- demean(
#'   X = base[, c("y", "x1")], f = base[, c("fe", "species")],
#'   slope.vars = base$x2, slope.flag = c(0, -1)
#' )
#'
#' # Explanations slope.flag = c(0, -1):
#' # - the first 0: the first factor (fe) is associated to no variable
#' # - the "-1":
#' #    * |-1| = 1: the second factor (species) is associated to ONE variable
#' #    *   -1 < 0: the second factor should not be included as such
#'
#' # Let's look at the equivalences
#' res_vs_1 <- feols(y ~ x1 + i(fe) + x2:species, base)
#' res_vs_2 <- feols(y ~ x1, base_dm)
#' res_vs_3 <- feols(y ~ x1, base_dm_bis)
#'
#' # only the small sample adj. differ in the SEs
#' etable(res_vs_1, res_vs_2, res_vs_3, keep = "x1")
#'
demean <- function(X, f, slope.vars, slope.flag, data, weights,
                   nthreads = getFixest_nthreads(), notes = getFixest_notes(),
                   iter = 2000, tol = 1e-6, fixef.reorder = TRUE, na.rm = TRUE,
                   as.matrix = is.atomic(X),
                   im_confident = FALSE, ...) {
  ANY_NA <- FALSE
  # SK: to reassign class if X is data.frame, data.table or tibble. Optimally you would preserve all attributes,
  # but using attributes<- is slow on data.frames. What I did in collapse is export the SET_ATTRIB and DUPLICATE_ATTRIB
  # functions from the C-API to use them internally in R -> copy attributes without checks at 0 cost, even for large data.frames.

  # LB: next line is needed if input data is matrix and as.matrix is set to FALSE
  clx <- NULL
  is_fixest <- inherits(X, "fixest")
  if (lX <- is.list(X) && !is_fixest) {
    clx <- oldClass(X)
    # SK: oldClass is faster and returns NULL when the list is plain. class returns the implicit class "list".
    # This is the key to fast R code -> all data.frame methods are super slow and should duly be avoided in internal code.
    oldClass(X) <- NULL
  }
  # LB: The reassignment of attributes to data.frames is actually a very good idea, thanks Seb!

  # LB: To avoid delayed evaluation problems (due to new default is.atomic(X))
  as_matrix <- as.matrix

  # internal argument
  fe_info <- FALSE

  # Step 1: formatting the input
  if (!im_confident) {
    check_arg(X, "numeric vmatrix | list | formula | class(fixest) mbt")
    check_arg(iter, "integer scalar GE{1}")
    check_arg(tol, "numeric scalar GT{0}")
    check_arg(notes, "logical scalar")

    validate_dots(valid_args = "fe_info", stop = TRUE)
    dots <- list(...)
    fe_info <- isTRUE(dots$fe_info)

    data_mbt <- TRUE
    if (inherits(X, "fixest")) {
      data_mbt <- FALSE
      if (!identical(X$method, "feols")) {
        stop("This function only works for 'feols' estimations (not for ", X$method, ").")
      }

      if (!fe_info && !is.null(X$y_demeaned)) {
        return(cbind(X$y_demeaned, X$X_demeaned))
      }

      if (!"fixef_vars" %in% names(X)) {
        stop("To apply demeaning, the estimation must contain fixed-effects, this is currently not the case.")
      }

      # Otherwise => we create data and X: formula
      data <- fetch_data(X)
      fml_linear <- X$fml_all$linear
      fml_fixef <- X$fml_all$fixef

      fml_vars <- .xpd(~ ..lhs + ..rhs,
        ..lhs = fml_linear[[2]],
        ..rhs = fml_linear[[3]]
      )

      if (!is.null(X$fml_all$iv)) {
        fml_iv <- X$fml_all$iv
        fml_vars <- .xpd(~ ..vars + ..iv_lhs + ..iv_rhs,
          ..vars = fml_vars,
          ..iv_lhs = fml_iv[[2]],
          ..iv_rhs = fml_iv[[3]]
        )
      }

      fml <- .xpd(lhs = fml_vars, rhs = fml_fixef)

      mc <- match.call()
      if (!"tol" %in% names(mc)) tol <- X$fixef.tol
      if (!"iter" %in% names(mc)) iter <- X$fixef.iter

      weights <- X[["weights"]]

      X <- fml

      as_matrix <- TRUE
    }

    #
    # X
    #

    fe_done <- FALSE
    if (is.call(X)) {
      if (data_mbt) check_arg(data, "data.frame mbt")
      check_arg(X, "ts formula var(data)", .data = data)

      # Extracting the information
      terms_fixef <- fixef_terms(.xpd(rhs = X[[3L]]))
      # We add the intercept only for only slope models, otherwise this would be void since the result would be always 0
      X <- fixest_model_matrix(.xpd(lhs = quote(y), rhs = X[[2L]]), data, fake_intercept = any(terms_fixef$slope_flag >= 0))
      var_names <- dimnames(X)[[2]]

      lX <- FALSE # Needed for the rest of the code to work

      # FE
      fe_done <- TRUE
      f <- unclass(error_sender(
        prepare_df(terms_fixef$fe_vars, data),
        "Error when evaluating the fixed-effects variables: "
      ))

      isSlope <- any(terms_fixef$slope_flag != 0)

      if (isSlope) {
        slope.vars <- unclass(error_sender(
          prepare_df(terms_fixef$slope_vars, data),
          "Error when evaluating the variable with varying slopes: "
        ))

        if (anyDuplicated(terms_fixef$slope_vars)) {
          # we need it!!!!
          slope.vars <- slope.vars[terms_fixef$slope_vars]
        }

        slope.flag <- terms_fixef$slope_flag

        is_numeric <- vapply(`attributes<-`(slope.vars, NULL), is.numeric, TRUE)
        if (!all(is_numeric)) {
          stop("In the fixed-effects part of the formula, variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope.vars)[!is_numeric], "s.is.quote"), " not.")
        }
      } else {
        slope.vars <- list(0)
        slope.flag <- rep(0L, length(f))
      }
    } else if (lX) {
      var_names <- names(X)
      if (!any(clx == "data.frame")) {
        n_all <- lengths(X)
        if (!all(n_all == n_all[1L])) {
          n_unik <- unique(n_all)
          stop("In argument 'X' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
        }
      }

      # SK: If your using Rcpp and the input is NumericVector or NumericMatrix, you should get a C++ error for wrongly typed data, so I don't see the need for this
      # LB: I disagree, such error messages are poor and hard to debug for end users.

      is_numeric <- vapply(`attributes<-`(X, NULL), is.numeric, TRUE)
      if (!all(is_numeric)) {
        # Faster than any(!is_numeric) => LB: yet a few nano seconds saved... ;-)
        n_non_num <- sum(!is_numeric)
        stop("All values in 'X' must be numeric. This is not the case for ", n_non_num, " variable", plural(n_non_num), ".")
      }
    } else {
      # Here: numeric matrix or vector
      var_names <- dimnames(X)[[2L]]
    }

    # nobs: useful for checks
    nobs <- if (is.list(X)) length(.subset2(X, 1L)) else NROW(X)

    #
    # f + slope.vars
    #

    if (!fe_done) {
      check_arg(f, "vmatrix | list mbt")
      check_arg(slope.vars, "numeric vmatrix | list")
      check_arg(slope.flag, "integer vector no na")

      if (is.list(f)) {
        if (!is.data.frame(f)) {
          n_all <- lengths(f)
          if (!all(n_all == n_all[1L])) {
            n_unik <- unique(n_all)
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
        f <- if (!is.array(f)) list(f) else unclass(as.data.frame(f))
      }

      is_pblm <- vapply(`attributes<-`(f, NULL), function(x) !(is.numeric(x) || is.character(x)), TRUE)
      if (any(is_pblm)) {
        for (i in which(is_pblm)) f[[i]] <- as.character(f[[i]])
      }

      if (length(f[[1L]]) != nobs) {
        stop("The number of observations in 'X' and in 'f' don't match (", if (lX) length(X[[1L]]) else NROW(X), " vs ", length(f[[1L]]), ").")
      }

      # Now the slopes
      isSlope <- FALSE
      if (!missnull(slope.vars) || !missnull(slope.flag)) {
        isSlope <- TRUE

        if (missnull(slope.vars)) stop("If argument 'slope.flag' is provided, you must also provide the argument 'slope.vars'.")
        if (missnull(slope.flag)) stop("If argument 'slope.vars' is provided, you must also provide the argument 'slope.flag'.")

        if (is.list(slope.vars)) {
          if (!is.data.frame(slope.vars)) {
            n_all <- lengths(slope.vars)
            if (!all(n_all == n_all[1L])) {
              n_unik <- unique(n_all)
              stop("In argument 'slope.vars' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
            }
          } else {
            oldClass(slope.vars) <- NULL
          }

          is_numeric <- vapply(`attributes<-`(slope.vars, NULL), is.numeric, TRUE)
          # SK: Much faster than !sapply(slope.vars, is.numeric)
          # LB: 10us diff, only impedes readability IMO
          if (!all(is_numeric)) stop("In the argument 'slope.vars', the variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope.vars)[!is_numeric], "s.is.quote"), " not.")
        } else {
          # SK: Same as above.. Necessary to convert ??
          # LB: The new C++ code accepts lists of vector/matrices or vector/matrices directly // So if we're here (meaning vector or matrix), that's fine
        }

        if (length(slope.flag) != length(f)) stop("The argument 'slope.flag' must be a vector representing, for each fixed-effect, the number of variables with varying slopes associated to it. Problem: the lengths of slope.flag and the fixed-effect differ: ", length(slope.flag), " vs ", length(f), ".")


        # SK: if(sum(abs(slope.flag)) != length(slope.vars))  : This means that slope flag can only be 0, 1 or -1. In the documentation you only talk about positive and negative values. The change I made reflects this.
        # LB: Cmon Sebastian why change sum(abs(slope.flag)) != length(slope.vars) into sum(slope.flag != 0L) != length(slope.vars)?  The time gains lie in a handful of nano second and you've just introduced a bug!

        n_sv <- if (is.list(slope.vars)) length(slope.vars) else NCOL(slope.vars)
        if (sum(abs(slope.flag)) != n_sv) stop("The argument 'slope.flag' must be a vector representing, for each fixed-effect, the number of variables with varying slopes associated to it. Problem: currently the number of variables in 'slope.flag' differ from the number of variables in 'slope.vars': ", sum(abs(slope.flag)), " vs ", n_sv, ".")
      } else {
        slope.vars <- list(0)
        slope.flag <- rep(0L, length(f))
      }
    }


    ## weights
    check_arg_plus(weights, "NULL numeric conv vector len(value) GE{0}", .value = nobs)
    is_weight <- !missing(weights) && !is.null(weights)

    ## nthreads (Also seems unnecessarily costly: 250 microseconds)
    # LB: I know... but it has to be checked
    if (!missing(nthreads)) nthreads <- check_set_nthreads(nthreads)

    #
    # FORMATTING
    #

    # NAs
    is_NA <- FALSE
    info_X <- cpp_which_na_inf(X, nthreads)

    if (info_X$any_na_inf) {
      is_NA <- info_X$is_na_inf
      n_na_x <- 1L
    } else {
      n_na_x <- 0L
    }


    if (anyNA(f)) {
      is_na_fe <- !complete.cases(f)
      n_na_fe <- 1L
      is_NA <- is_NA | is_na_fe
    } else {
      n_na_fe <- 0L
    }


    is_na_slope <- 0L
    if (isSlope) {
      info_slopes <- cpp_which_na_inf(slope.vars, nthreads)

      if (info_slopes$any_na_inf) {
        is_na_slope <- info_slopes$is_na_inf
        is_NA <- is_NA | is_na_slope
      }
    }


    is_na_weights <- 0L
    if (is_weight) {
      info_weights <- cpp_which_na_inf(weights, nthreads)

      if (info_weights$any_na_inf) {
        is_na_weights <- info_weights$is_na_inf
        is_NA <- is_NA | is_na_weights
      }
    }

    n_na <- sum(is_NA)
    if (n_na && notes) {
      # Here is all the error message stuff now: Only evaluated if needed
      if (n_na_x) n_na_x <- sum(info_X$is_na_inf)
      if (n_na_fe) n_na_fe <- sum(is_na_fe)
      slopes_msg <- if (isSlope) paste0(", slope vars: ", sum(is_na_slope)) else ""
      weight_msg <- if (is_weight) paste0(", weights: ", sum(is_na_weights)) else ""

      if (n_na == length(is_NA)) stop("All observations contain NA values (Breakup: X: ", n_na_x, ", f: ", n_na_fe, slopes_msg, weight_msg, ").")
      message("NOTE: ", fsignif(n_na), " observation", plural(n_na), " removed because of NA values (Breakup: X: ", n_na_x, ", f: ", n_na_fe, slopes_msg, weight_msg, ").")
    }

    if (n_na) {
      ANY_NA <- TRUE
      # SK: subsetting integer is faster than logical !!
      # LB: Indeed, just checked, quite a diff. Good to know!
      cc <- which(!is_NA)
      X <- select_obs(X, cc)

      f <- lapply(f, `[`, cc)
      if (isSlope) slope.vars <- lapply(slope.vars, `[`, cc)
      if (is_weight) weights <- weights[cc]
    }

    if (!is_weight) weights <- 1
  } else {
    # Need this here, otherwise error:
    var_names <- if (lX) names(X) else dimnames(X)[[2L]]
    if (missing(weights) || is.null(weights)) weights <- 1
    if (missnull(slope.vars) || missnull(slope.flag)) {
      slope.vars <- list(0)
      slope.flag <- rep(0L, length(unclass(f)))
      # SK: unclass gives extra speed if f is data.frame, and no cost if f is list.
    }
  }

  #
  # Unclassing fes
  #

  quf_info_all <- quf_table_sum(
    x = f, y = 0, do_sum_y = FALSE, rm_0 = FALSE,
    rm_1 = FALSE, rm_single = FALSE, do_refactor = FALSE,
    r_x_sizes = 0L, obs2keep = 0L, only_slope = slope.flag < 0L,
    nthreads = as.integer(nthreads)
  )

  # table/sum_y/sizes
  fixef_table <- quf_info_all$table
  fixef_sizes <- lengths(fixef_table)

  if (fixef.reorder) {
    new_order <- order(fixef_sizes, decreasing = TRUE)
    if (is.unsorted(new_order)) {
      fixef_table <- fixef_table[new_order]
      fixef_sizes <- fixef_sizes[new_order]
      quf_info_all$quf <- quf_info_all$quf[new_order]

      # We reorder slope.vars only if needed (because it is a pain)
      if (sum(slope.flag != 0) >= 2) {
        # here slope.vars have 2 or more variables:
        # either matrix or data.frame/list
        if (!is.list(slope.vars)) {
          slope.vars <- as.data.frame(slope.vars)
        }
        # we ensure it's a plain list
        slope.vars <- unclass(slope.vars)

        n_fixef <- length(slope.flag)
        slope_vars_list <- vector("list", n_fixef)
        id_current <- 0
        for (i in 1:n_fixef) {
          n_slopes <- abs(slope.flag[i])

          if (n_slopes == 0) {
            slope_vars_list[[i]] <- 0
          } else {
            slope_vars_list[[i]] <- slope.vars[1:n_slopes + id_current]
            id_current <- id_current + n_slopes
          }
        }

        # Reordering => currently it's quite clumsy (but I HATE VSs!)
        slope.vars <- vector("list", sum(abs(slope.flag)))
        slope_vars_list_reordered <- slope_vars_list[new_order]
        id_current <- 0
        for (i in 1:n_fixef) {
          vs <- slope_vars_list_reordered[[i]]
          if (is.list(vs)) {
            n_vs <- length(vs)
            for (j in 1:n_vs) {
              slope.vars[[id_current + j]] <- vs[[j]]
            }
            id_current <- id_current + n_vs
          }
        }

        # slope.flag must be reordered afterwards
        slope.flag <- slope.flag[new_order]
      }
    }
  }

  fixef_table_vector <- unlist(fixef_table)
  if (!is.integer(slope.flag)) slope.flag <- as.integer(slope.flag)

  #
  # The algorithm
  #

  # y => list, X => matrix
  if (as_matrix) {
    y <- 0
  } else {
    # Quirk => y returns a list only if NCOL(y) > 1 or is.list(y)
    y <- if (lX) X else list(X)
    # SK: This does the same, a bit more frugal
    # LB: I tend to prefer late evaluations instead of lX which requires bookkeeping
    X <- 0
  }

  vars_demean <- cpp_demean(y, X, weights,
    iterMax = iter,
    diffMax = tol, r_nb_id_Q = fixef_sizes,
    fe_id_list = quf_info_all$quf, table_id_I = fixef_table_vector,
    slope_flag_Q = slope.flag, slope_vars_list = slope.vars,
    r_init = 0, nthreads = nthreads, save_fixef = FALSE
  )

  # Internal call
  if (fe_info) {
    res <- list(
      y = vars_demean$y_demean, X = vars_demean$X_demean,
      weights = weights, fixef_id_list = quf_info_all$quf,
      slope_vars = slope.vars, slope_flag = slope.flag, varnames = colnames(X)
    )

    return(res)
  }


  if (as_matrix) {
    if (ANY_NA && na.rm == FALSE) {
      K <- NCOL(vars_demean$X_demean)
      res <- matrix(NA_real_, length(is_NA), K)
      res[cc, ] <- vars_demean$X_demean
    } else {
      res <- vars_demean$X_demean
    }

    dimnames(res)[[2L]] <- var_names
  } else {
    if (ANY_NA && na.rm == FALSE) {
      n <- length(is_NA)
      # SK: One-line efficient solution to the task:
      res <- lapply(vars_demean$y_demean, function(x) `[<-`(rep(NA_real_, n), cc, value = x))
      # LB: nice compactness
    } else {
      res <- vars_demean$y_demean
    }

    names(res) <- var_names
    # SK: This makes your function class-agnostic to lists, data.frame's, data.table's and tibbles (with the exception for any custom data.frame row.names which are not preserved, but the as.data.frame solution also did not do that)
    if (!lX || length(clx)) attr(res, "row.names") <- .set_row_names(length(res[[1L]]))
    # SK: Here !lX || length(clx) means add row.names if either X was a matrix or is classed that is a data.frame
    oldClass(res) <- if (lX) clx else "data.frame"
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
#' base <- iris
#' names(base) <- c("y", "x1", "x2", "x3", "species")
#' base$y[1:5] <- NA
#'
#' # Split sample estimations
#' est_split <- feols(y ~ x1, base, split = ~species)
#' (obs_setosa <- obs(est_split$setosa))
#' (obs_versi <- obs(est_split$versicolor))
#'
#' est_versi <- feols(y ~ x1, base, subset = obs_versi)
#'
#' etable(est_split, est_versi)
#'
obs <- function(x) {
  check_arg(x, "class(fixest)")

  if (isTRUE(x$lean)) {
    stop("obs() does not work with models estimated with 'lean = TRUE'.")
  }

  id <- 1:x$nobs_origin

  for (i in seq_along(x$obs_selection)) {
    id <- id[x$obs_selection[[i]]]
  }

  return(id)
}




#' Check the fixed-effects convergence of a \code{feols} estimation
#'
#' Checks the convergence of a \code{feols} estimation by computing the first-order conditions of all fixed-effects (all should be close to 0)
#'
#' @param x A \code{\link[fixest2]{feols}} estimation that should contain fixed-effects.
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
#' base <- setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' base$FE <- rep(1:30, 5)
#'
#' # one estimation with fixed-effects + varying slopes
#' est <- feols(y ~ x1 | species[x2] + FE[x3], base)
#'
#' # Checking the convergence
#' conv <- check_conv_feols(est)
#'
#' # We can check that al values are close to 0
#' summary(conv)
#'
#' summary(conv, "detail")
#'
check_conv_feols <- function(x) {
  check_arg(x, "class(fixest) mbt")
  if (!identical(x$method, "feols")) {
    stop("This function only works for 'feols' estimations (not for ", x$method, ").")
  }

  if (!"fixef_vars" %in% names(x)) {
    message("This function only works with fixed-effects (which are currently not present).")
    return(NULL)
  }

  # fixef names for information
  new_order <- x$fe.reorder
  fixef_vars <- x$fixef_vars[new_order]
  if (!is.null(x$fixef_terms)) {
    slope_variables <- x$slope_variables_reordered
    slope_flag <- x$slope_flag_reordered

    # We reconstruct the terms
    fixef_names <- c()
    start <- c(0, cumsum(abs(slope_flag)))
    for (i in seq_along(slope_flag)) {
      sf <- slope_flag[i]
      if (sf >= 0) {
        fixef_names <- c(fixef_names, fixef_vars[i])
      }

      if (abs(sf) > 0) {
        fixef_names <- c(fixef_names, paste0(fixef_vars[i], "[[", names(slope_variables)[start[i] + 1:abs(sf)], "]]"))
      }
    }
  } else {
    fixef_names <- fixef_vars
  }


  info <- demean(x, fe_info = TRUE)

  res <- check_conv(
    y = info$y, X = info$X, fixef_id_list = info$fixef_id_list,
    slope_flag = info$slope_flag, slope_vars = info$slope_vars,
    weights = info$weights, full = TRUE, fixef_names = fixef_names
  )

  names(res) <- info$varnames

  class(res) <- "fixest_check_conv"

  res
}

#' Replicates \code{fixest} objects
#'
#' Simple function that replicates \code{fixest} objects while (optionally) computing different standard-errors. Useful mostly in combination with \code{\link[fixest2]{etable}} or \code{\link[fixest2]{coefplot}}.
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
#' est <- feols(Ozone ~ Solar.R + Wind + Temp, data = airquality)
#'
#' my_vcov <- list(~Month, ~Day, ~ Day + Month)
#'
#' etable(rep(est, vcov = my_vcov))
#'
#' coefplot(rep(est, vcov = my_vcov), drop = "Int")
#'
#' #
#' # To rep multiple objects, you need to use .l()
#' #
#'
#' est_bis <- feols(Ozone ~ Solar.R + Wind + Temp | Month, airquality)
#'
#' etable(rep(.l(est, est_bis), vcov = my_vcov))
#'
#' # using each
#' etable(rep(.l(est, est_bis), each = 3, vcov = my_vcov))
#'
rep.fixest <- function(x, times = 1, each = 1, vcov, ...) {
  # each is applied first, then times
  # x can be either a list of fixest objects, either a fixest object

  check_arg(x, "class(fixest, fixest_list) mbt")
  check_arg(times, "integer scalar GE{1} | integer vector no na GE{0}")
  check_arg(each, "integer scalar GE{1} | logical scalar")
  check_arg(vcov, "class(list)")

  if (is_user_level_call()) {
    validate_dots(suggest_args = c("times", "each"), stop = TRUE)
  }

  # Checking the arguments
  IS_LIST <- FALSE
  if ("fixest_list" %in% class(x)) {
    IS_LIST <- TRUE
    class(x) <- "list"

    n <- length(x)
  } else {
    n <- 1
  }

  if (is.logical(each) && each == FALSE) {
    stop("Argument 'each' cannot be equal to FALSE.")
  }

  IS_MULTI_VCOV <- !missing(vcov)
  if (IS_MULTI_VCOV) {
    n_vcov <- length(vcov)

    if (times == 1 && each == 1) {
      if (isTRUE(each)) {
        each <- n_vcov
      } else {
        times <- n_vcov
      }
    }
  }

  res_int <- rep(1:n, times = times, each = each)
  n_res <- length(res_int)

  if (IS_MULTI_VCOV) {
    # Checking and expanding

    vcov_mapping <- 1:n_res
    if (times == 1) {
      if (n_vcov != each && n_vcov != n_res) {
        stop("In rep, the argument 'vcov' (currently of length ", n_vcov, ") must be a list either of length ", each, " or of length ", n_res, ".")
      }

      if (n_vcov == each) vcov_mapping <- rep(1:each, times = n)
    } else if (each == 1) {
      if (n_vcov != times && n_vcov != n_res) {
        stop("In rep, the argument 'vcov' (currently of length ", n_vcov, ") must be a list either of length ", times, " or of length ", n_res, ".")
      }

      if (n_vcov == times) vcov_mapping <- rep(1:n_vcov, each = n)
    } else {
      if (n_vcov != n_res) {
        stop("In rep, the argument 'vcov' (currently of length ", n_vcov, ") must be a list either of length ", n_res, ".")
      }
    }
  }

  res <- vector("list", length(res_int))

  if (IS_MULTI_VCOV) {
    for (i in 1:n_res) {
      if (IS_LIST) {
        res[[i]] <- summary(x[[res_int[i]]], vcov = vcov[[vcov_mapping[i]]])
      } else {
        res[[i]] <- summary(x, vcov = vcov[[vcov_mapping[i]]])
      }
    }
  } else {
    for (i in unique(res_int)) {
      if (IS_LIST) {
        res[res_int == i] <- x[[i]]
      } else {
        res[res_int == i] <- list(x)
      }
    }
  }

  for (i in 1:n_res) {
    res[[i]]$model_id <- res_int[i]
  }

  res
}

#' @rdname rep.fixest
rep.fixest_list <- function(x, times = 1, each = 1, vcov, ...) {
  rep.fixest(x, times = times, each = each, vcov = vcov, ...)
}

#' @rdname rep.fixest
.l <- function(...) {
  check_arg(..., "mbt class(fixest) | list")

  dots <- list(...)
  if (all(sapply(dots, function(x) "fixest" %in% class(x)))) {
    class(dots) <- "fixest_list"

    return(dots)
  }

  if (length(dots) == 1) {
    if ("fixest_multi" %in% class(dots[[1]])) {
      res <- attr(dots[[1]], "data")
      class(res) <- "fixest_list"
      return(res)
    }

    if (all(sapply(dots[[1]], function(x) "fixest" %in% class(x)))) {
      res <- dots[[1]]
      class(res) <- "fixest_list"
      return(res)
    }
  }

  res <- list()
  for (i in seq_along(dots)) {
    if ("fixest" %in% class(dots[[i]])) {
      res[[length(res) + 1]] <- dots[[i]]
    } else {
      obj <- dots[[i]]

      if (class(obj) %in% "fixest_multi") {
        data <- attr(obj, "data")
        for (j in seq_along(data)) {
          res[[length(res) + 1]] <- data[[j]]
        }
      } else {
        for (j in seq_along(obj)) {
          if (!"fixest" %in% class(obj[[j]])) {
            stop("In .l(...), each argument must be either a fixest object, or a list of fixest objects. Problem: The ", n_th(j), " element of the ", n_th(i), " argument (the latter being a list) is not a fixest object.")
          }

          res[[length(res) + 1]] <- obj[[j]]
        }
      }
    }
  }

  class(res) <- "fixest_list"
  res
}
