#' A print facility for \code{fixest} objects.
#'
#' This function is very similar to usual \code{summary} functions as it provides the table of coefficients along with other information on the fit of the estimation. The type of output can be customized by the user (using function \code{setFixest_print}).
#'
#' @method print fixest
#'
#' @param x A \code{fixest} object. Obtained using the methods \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param n Integer, number of coefficients to display. By default, only the first 8 coefficients are displayed if \code{x} does not come from \code{\link[fixest]{summary.fixest}}.
#' @param type Either \code{"table"} (default) to display the coefficients table or \code{"coef"} to display only the coefficients.
#' @param fitstat A formula or a character vector representing which fit statistic to display. The types must be valid types of the function \code{\link[fixest]{fitstat}}. The default fit statistics depend on the type of estimation (OLS, GLM, IV, with/without fixed-effect). Providing the argument \code{fitstat} overrides the default fit statistics, you can however use the point "." to summon them back. Ex 1: \code{fitstat = ~ . + ll} adds the log-likelihood to the default values. Ex 2: \code{fitstat = ~ ll + pr2} only displays the log-likelihood and the pseudo-R2.
#' @param ... Other arguments to be passed to \code{\link[fixest]{vcov.fixest}}.
#'
#' @details
#'  It is possible to set the default values for the arguments \code{type} and \code{fitstat} by using the function \code{setFixest_print}.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade
#' #   => we account for 3 fixed-effects (FEs)
#' est_pois <- fepois(Euros ~ log(dist_km) | Origin + Destination + Product, trade)
#'
#' # displaying the results
#' #  (by default SEs are clustered if FEs are used)
#' print(est_pois)
#'
#' # By default the coefficient table is displayed.
#' #  If the user wished to display only the coefficents, use option type:
#' print(est_pois, type = "coef")
#'
#' # To permanently display coef. only, use setFixest_print:
#' setFixest_print(type = "coef")
#' est_pois
#' # back to default:
#' setFixest_print(type = "table")
#'
#' #
#' # fitstat
#' #
#'
#' # We modify which fit statistic to display
#' print(est_pois, fitstat = ~ . + lr)
#'
#' # We add the LR test to the default (represented by the ".")
#'
#' # to show only the LR stat:
#' print(est_pois, fitstat = ~ . + lr.stat)
#'
#' # To modify the defaults:
#' setFixest_print(fitstat = ~ . + lr.stat + rmse)
#' est_pois
#'
#' # Back to default (NULL == default)
#' setFixest_print(fitstat = NULL)
#'
print.fixest <- function(x, n, type = "table", fitstat = NULL, ...) {
  # checking the arguments
  if (is_user_level_call()) {
    validate_dots(
      suggest_args = c("n", "type", "vcov"),
      valid_args = c("vcov", "se", "cluster", "ssc", "forceCovariance", "keepBounded")
    )
  }

  # The objects from the estimation and the summary are identical, except regarding the vcov
  fromSummary <- isTRUE(x$summary)

  if (!missnull(fitstat)) {
    fitstat <- fitstat_validate(fitstat, TRUE)
  }

  # User options
  set_defaults("fixest_print")

  # if NOT from summary, we consider the argument 'type'
  if (!fromSummary) {
    # checking argument type
    check_arg_plus(type, "match(coef, table)")

    if (type == "coef") {
      print(coef(x))
      return(invisible())
    }
  }

  isNegbin <- x$method == "fenegbin" || (x$method %in% c("femlm", "feNmlm") && x$family == "negbin")

  x <- summary(x, fromPrint = TRUE, ...)

  check_arg(n, "integer scalar GE{1}")

  msgRemaining <- ""
  nb_coef <- length(coef(x)) - isNegbin
  if (missing(n) && is.null(x$n_print)) {
    if (fromSummary && !isTRUE(x$summary_from_fit)) {
      n <- Inf
    } else {
      if (nb_coef <= 10) {
        n <- 10
      } else {
        n <- 8
        msgRemaining <- paste0("... ", nb_coef - n, " coefficients remaining (display them with summary() or use argument n)\n")
      }
    }
  } else {
    if (!is.null(x$n_print)) n <- x$n_print

    if (n < nb_coef) {
      msgRemaining <- paste0("... ", nb_coef - n, " coefficients remaining\n")
    }
  }

  # We also add the collinearity message
  collinearity_msg <- ""
  if (!is.null(x$collin.var)) {
    n_collin <- length(x$collin.var)
    collinearity_msg <- paste0("... ", n_collin, " variable", plural(n_collin, "s.was"), " removed because of collinearity (", enumerate_items(x$collin.var, nmax = 3), ifelse(n_collin > 3, " [full set in $collin.var]", ""), ")\n")
    if (isTRUE(x$iv) && any(grepl("^fit_", x$collin.var))) {
      if (!any(grepl("^fit_", names(x$coefficients)))) {
        iv_msg <- "NOTE: all endogenous regressors were removed.\n"
      } else {
        n_rm <- sum(grepl("^fit_", x$collin.var))
        iv_msg <- paste0("Important note: ", n_letter(n_rm), " endogenous regressor", plural(n_rm, "s.was"), " removed => IV estimation not valid.\n")
      }

      collinearity_msg <- paste0(collinearity_msg, iv_msg)
    }
  }

  if (isFALSE(x$convStatus)) {
    last_warn <- getOption("fixest_last_warning")
    if (is.null(last_warn) || (proc.time() - last_warn)[3] > 1) {
      if (x$method %in% c("femlm", "feNmlm", "fenegbin")) {
        warning("The optimization algorithm did not converge, the results are not reliable. (", x$message, ")", call. = FALSE)
      } else if (x$method_type == "feols") {
        warning("The demeaning algorithm did not converge, the results are not reliable. (", x$message, ")", call. = FALSE)
      } else {
        warning("The GLM algorithm did not converge, the results are not reliable. (", x$message, ")", call. = FALSE)
      }
    }
  }

  coeftable <- x$coeftable

  # The type of SE
  se.type <- attr(coeftable, "type")
  if (is.null(se.type)) se.type <- "Custom"

  if (x$method_type %in% c("femlm", "feNmlm")) {
    family_format <- c(poisson = "Poisson", negbin = "Negative Binomial", logit = "Logit", gaussian = "Gaussian")
    msg <- ifelse(is.null(x$call$NL.fml), "", "Non-linear ")
    half_line <- paste0(msg, "ML estimation, family = ", family_format[x$family])
  } else if (x$method %in% c("feglm", "feglm.fit")) {
    fam_call <- x$call$family
    if (is.null(names(fam_call))) {
      half_line <- paste0("GLM estimation, family = ", x$family$family)
    } else {
      half_line <- paste0("GLM estimation, family = ", deparse_long(fam_call))
    }
  } else if (x$method == "fepois") {
    half_line <- "Poisson estimation"
  } else if (x$method == "fenegbin") {
    half_line <- "Negative Binomial ML estimation"
  } else {
    half_line <- "OLS estimation"
  }

  if (isTRUE(x$iv)) {
    glue <- function(...) paste(..., collapse = ", ")
    first_line <- paste0("TSLS estimation, Dep. Var.: ", as.character(x$fml)[[2]], ", Endo.: ", glue(get_vars(x$iv_endo_fml)), ", Instr.: ", glue(x$iv_inst_names), "\n")
    second_line <- paste0(ifunit(x$iv_stage, "First", "Second"), " stage: Dep. Var.: ", as.character(x$fml)[[2]], "\n")
    cat(first_line, second_line, sep = "")
  } else {
    cat(half_line, ", Dep. Var.: ", as.character(x$fml)[[2]], "\n", sep = "")
  }


  cat("Observations:", addCommas(x$nobs), "\n")
  if (!is.null(x$fixef_terms)) {
    terms_full <- extract_fe_slope(x$fixef_terms)
    fixef_vars <- terms_full$fixef_vars

    if (length(fixef_vars) > 0) {
      cat("Fixed-effects: ", paste0(fixef_vars, ": ", addCommas(x$fixef_sizes[fixef_vars]), collapse = ",  "), "\n", sep = "")
    }

    cat("Varying slopes: ", paste0(terms_full$slope_vars, " (", terms_full$slope_fe, ": ", addCommas(x$fixef_sizes[terms_full$slope_fe]), ")", collapse = ",  "), "\n", sep = "")
  } else {
    if (!is.null(x$fixef_sizes)) cat("Fixed-effects: ", paste0(x$fixef_vars, ": ", addCommas(x$fixef_sizes), collapse = ",  "), "\n", sep = "")
  }


  if (is.null(x$onlyFixef)) {
    cat("Standard-errors:", se.type, "\n")

    last_line <- paste0(msgRemaining, collinearity_msg)

    # The matrix of coefficients
    if (isNegbin) {
      if (nrow(coeftable) == 2) {
        new_table <- coeftable[1, , drop = FALSE]
      } else {
        new_table <- coeftable[-nrow(coeftable), ]
      }

      print_coeftable(head(new_table, n), lastLine = last_line)

      theta <- coeftable[".theta", 1]
      noDispInfo <- ifelse(theta > 1000, "(theta >> 0, no sign of overdispersion, you may consider a Poisson model)", "")
      cat("Over-dispersion parameter: theta =", theta, noDispInfo, "\n")
    } else {
      print_coeftable(head(coeftable, n), lastLine = last_line)
    }
  }

  if (isTRUE(x$NA_model)) {
    return(invisible())
  }

  if (!is.null(fitstat) && identical(fitstat, NA)) {
    # No fitstat
  } else {
    if (is.null(fitstat) || "." %in% fitstat) {
      if (x$method_type == "feols") {
        default_fit <- c("rmse", "ar2")

        if (!is.null(x$fixef_sizes) && is.null(x$onlyFixef)) {
          default_fit <- c(default_fit, "wr2")
        }

        if (isTRUE(x$iv)) {
          default_fit <- c(default_fit, "ivf1", "wh", "sargan")
        }
      } else {
        default_fit <- c("ll", "apr2", "bic", "cor2")
      }

      if ("." %in% fitstat) {
        fitstat <- setdiff(c(default_fit, fitstat), ".")
      } else {
        fitstat <- default_fit
      }
    }

    print(fitstat(x, fitstat), na.rm = TRUE, group.solo = TRUE)
  }

  if (isFALSE(x$convStatus)) {
    iter_format <- x$iterations
    if (length(iter_format) == 1) {
      iter_format <- paste0("lhs: ", iter_format)
    } else {
      n_iter <- length(iter_format)
      iter_format <- paste0("lhs: ", iter_format[n_iter], ", rhs: ", paste0(head(iter_format, min(n_iter - 1, n)), collapse = ", "))
    }
    cat("# Evaluations:", iter_format, "--", x$message, "\n")
  }
}

#' @rdname n_unik
print.vec_n_unik <- function(x, ...) {
  hash <- "## "

  x_names <- sfill(names(x))
  na.info <- attr(x, "na.info")
  na.info_format <- sfill(fsignif(na.info))
  x_value <- sfill(fsignif(x))

  n <- length(x)
  na_col <- paste0("(# NAs: ", na.info_format, ")")
  na_col[na.info == 0] <- ""

  res <- paste0(hash, x_names, ": ", x_value, " ", na_col)
  cat(res, sep = "\n")
}

#' @rdname n_unik
print.list_n_unik <- function(x, ...) {
  # I can't use #> anymore!!! The auto complete by the new version
  # of Rstudio drives me nuts!!!
  hash <- "## "

  x_all <- x
  n_x <- length(x_all)

  if (n_x == 1) {
    return(x_all[[1]])
  }

  info_pairs <- attr(x, "info_pairs")
  IS_PAIR <- !is.null(info_pairs)

  # First, the variable names

  all_names <- names(x_all[[1]])
  qui_0 <- nchar(all_names) == 0
  i <- 2
  while (any(qui_0) && i <= n_x) {
    names_i <- names(x_all[[i]])
    all_names[qui_0] <- names_i[qui_0]
    qui_0 <- nchar(all_names) == 0
    i <- i + 1
  }

  # id_vars: used in info_pairs
  id_vars <- seq_along(all_names)

  if (all(qui_0)) {
    stop("Not any valid information to display: please check that all your variables exist in all data sets.")
  } else if (any(qui_0)) {
    warning("Some variables could not be evaluated in any data set: please check that all your variables exist in all data sets.")

    all_names <- all_names[!qui_0]
    id_vars <- id_vars[!qui_0]

    for (i in 1:n_x) {
      x <- x_all[[i]]
      na.info <- attr(x, "na.info")

      x <- x[!qui_0]
      na.info <- na.info[!qui_0]
      attr(x, "na.info") <- na.info

      x_all[[i]] <- x
    }
  }


  x_names <- sfill(all_names)

  # If ambiguous pairs: we add data set suffixes
  if (n_x > 2 && IS_PAIR) {
    add_suffix <- function(x, i) paste0(x, " (", letters[i], ")")
    var_width <- max(nchar("Exclusive "), nchar(x_names[1]))
  } else {
    add_suffix <- function(x, i) x
    var_width <- nchar(x_names[1])
  }

  var_col <- paste0(hash, x_names, ":")
  na_intro <- paste0(hash, sfill(" ", var_width), "|NAs:")
  var_col <- insert_in_between(var_col, na_intro)
  # we add the first row of the data sets names + format
  var_col <- sfill(c(hash, var_col), right = TRUE)

  print_mat <- var_col

  KEEP_NA_ROW <- rep(FALSE, length(x_names))

  for (i in 1:n_x) {
    data_name <- add_suffix(names(x_all)[i], i)

    x <- x_all[[i]]

    na.info <- attr(x, "na.info")
    KEEP_NA_ROW <- KEEP_NA_ROW | na.info > 0
    na.info_format <- fsignif(na.info)
    x_value <- fsignif(x)
    x_value[is.na(x)] <- "--"
    na.info_format[is.na(x)] <- "--"


    width <- max(nchar(na.info_format), nchar(x_value))

    na.info_format <- sfill(na.info_format, width)
    x_value <- sfill(x_value, width)

    x_col <- insert_in_between(x_value, na.info_format)
    x_col <- sfill(c(data_name, x_col))

    print_mat <- cbind(print_mat, x_col)
  }

  if (!any(KEEP_NA_ROW)) {
    print_mat[, 1] <- substr(print_mat[, 1], 1, nchar(print_mat[1, 1]) - 4)
  }

  keep <- c(TRUE, insert_in_between(TRUE, KEEP_NA_ROW))

  print_mat <- print_mat[keep, ]

  # PAIR information
  if (IS_PAIR) {
    # identifiers used for insertion
    id_vars_all <- c(0, insert_in_between(id_vars, 0))
    id_vars_all <- id_vars_all[keep]

    # we add two columns: id and common
    print_mat <- cbind(print_mat, id_vars_all, "")

    insert_row_after_id <- function(mat, row, col_id) {
      for (i in 1:nrow(row)) {
        i_id_mat <- max(which(mat[, col_id] == row[i, col_id]))
        tmp_before <- mat[1:i_id_mat, ]

        if (i_id_mat < nrow(mat)) {
          tmp_after <- mat[(i_id_mat + 1):nrow(mat), ]
        } else {
          tmp_after <- NULL
        }

        mat <- rbind(tmp_before, row[i, ], tmp_after)
      }

      mat
    }

    for (pair in info_pairs) {
      data_id <- attr(pair, "data_id")
      pair <- as.data.frame(pair)

      if (n_x == 2) {
        # data_col = paste0(hash, sfill(" ", var_width, right = TRUE), "|Excl:")
        data_col <- paste0(hash, sfill("# Exclusive ", var_width, right = TRUE), "|")
      } else {
        data_col <- paste0(
          hash, sfill("# Exclusive ", var_width, right = TRUE),
          "|", paste0(letters[data_id], collapse = ":")
        )
      }

      # formatting the NAs
      for (i in 1:(ncol(pair) - 2)) {
        pair[[i]] <- fsignif(pair[[i]])
      }
      pair[is.na(pair)] <- " "
      pair <- as.matrix(pair)

      # adding the data col
      pair <- cbind(data_col, pair)

      print_mat <- insert_row_after_id(print_mat, pair, n_x + 2)
    }

    # Formatting
    print_mat <- print_mat[, -(n_x + 2)]
    print_mat[, 1] <- sfill(print_mat[, 1], right = TRUE)
    for (j in 2:(n_x + 1)) {
      print_mat[, j] <- sfill(print_mat[, j])
    }

    # Last column
    common <- print_mat[, n_x + 2]
    is_num <- grepl("\\d", common)
    common[is_num] <- sfill(fsignif(as.numeric(common[is_num])))
    common[is_num] <- paste0("# Common: ", common[is_num])
    print_mat[, n_x + 2] <- common
  }

  vec2print <- apply(print_mat, 1, paste, collapse = " ")

  cat(vec2print, sep = "\n")
}

#' @rdname osize
print.osize <- function(x, ...) {
  cat(x, "\n")
}

#' @rdname print.fixest
setFixest_print <- function(type = "table", fitstat = NULL) {
  check_arg_plus(type, "match(coef, table)")

  if (!missnull(fitstat)) {
    fitstat <- fitstat_validate(fitstat, TRUE)
  }

  # Getting the existing defaults
  opts <- getOption("fixest_print")

  if (is.null(opts)) {
    opts <- list()
  } else if (!is.list(opts)) {
    warning("Wrong formatting of option 'fixest_print', all options are reset.")
    opts <- list()
  }

  # Saving the default values
  mc <- match.call()
  args_default <- names(mc)[-1]

  # NOTA: we don't allow delayed evaluation => all arguments must have hard values
  for (v in args_default) {
    opts[[v]] <- eval(as.name(v))
  }

  options(fixest_print = opts)
}


#' @rdname print.fixest
getFixest_print <- function() {
  x <- getOption("fixest_print")
  if (!(is.null(x) || is.list(x))) {
    stop("The value of getOption(\"fixest_print\") is currently not legal. Please use function setFixest_print to set it to an appropriate value. ")
  }

  x
}
