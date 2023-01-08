#' Summary of a \code{fixest} object. Computes different types of standard errors.
#'
#' This function is similar to \code{print.fixest}. It provides the table of coefficients along with other information on the fit of the estimation. It can compute different types of standard errors. The new variance covariance matrix is an object returned.
#'
#' @inheritParams feNmlm
#' @inheritParams aggregate.fixest
#'
#' @method summary fixest
#' @param vcov Versatile argument to specify the VCOV. In general, it is either a character scalar equal to a VCOV type, either a formula of the form: \code{vcov_type ~ variables}. The VCOV types implemented are: "iid", "hetero" (or "HC1"), "cluster", "twoway", "NW" (or "newey_west"), "DK" (or "driscoll_kraay"), and "conley". It also accepts object from \code{\link[fixest2]{vcov_cluster}}, \code{\link[fixest:vcov_hac]{vcov_NW}}, \code{\link[fixest:vcov_hac]{NW}}, \code{\link[fixest:vcov_hac]{vcov_DK}}, \code{\link[fixest:vcov_hac]{DK}}, \code{\link[fixest2]{vcov_conley}} and \code{\link[fixest:vcov_conley]{conley}}. It also accepts covariance matrices computed externally. Finally it accepts functions to compute the covariances. See the `vcov` documentation in the \href{https://lrberge.github.io/fixest/articles/fixest_walkthrough.html#the-vcov-argument-1}{vignette}.
#' @param se Character scalar. Which kind of standard error should be computed: \dQuote{standard}, \dQuote{hetero}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}? By default if there are clusters in the estimation: \code{se = "cluster"}, otherwise \code{se = "iid"}. Note that this argument is deprecated, you should use \code{vcov} instead.
#' @param cluster Tells how to cluster the standard-errors (if clustering is requested). Can be either a list of vectors, a character vector of variable names, a formula or an integer vector. Assume we want to perform 2-way clustering over \code{var1} and \code{var2} contained in the data.frame \code{base} used for the estimation. All the following \code{cluster} arguments are valid and do the same thing: \code{cluster = base[, c("var1", "var2")]}, \code{cluster = c("var1", "var2")}, \code{cluster = ~var1+var2}. If the two variables were used as fixed-effects in the estimation, you can leave it blank with \code{vcov = "twoway"} (assuming \code{var1} [resp. \code{var2}] was the 1st [res. 2nd] fixed-effect). You can interact two variables using \code{^} with the following syntax: \code{cluster = ~var1^var2} or \code{cluster = "var1^var2"}.
#' @param stage Can be equal to \code{2} (default), \code{1}, \code{1:2} or \code{2:1}. Only used if the object is an IV estimation: defines the stage to which \code{summary} should be applied. If \code{stage = 1} and there are multiple endogenous regressors or if \code{stage} is of length 2, then an object of class \code{fixest_multi} is returned.
#' @param object A \code{fixest} object. Obtained using the functions \code{\link[fixest2]{femlm}}, \code{\link[fixest2]{feols}} or \code{\link[fixest2]{feglm}}.
#' @param ssc An object of class \code{ssc.type} obtained with the function \code{\link[fixest2]{ssc}}. Represents how the degree of freedom correction should be done.You must use the function \code{\link[fixest2]{ssc}} for this argument. The arguments and defaults of the function \code{\link[fixest2]{ssc}} are: \code{adj = TRUE}, \code{fixef.K="nested"}, \code{cluster.adj = TRUE}, \code{cluster.df = "min"}, \code{t.df = "min"}, \code{fixef.force_exact=FALSE)}. See the help of the function \code{\link[fixest2]{ssc}} for details.
#' @param .vcov A user provided covariance matrix or a function computing this matrix. If a matrix, it must be a square matrix of the same number of rows as the number of variables estimated. If a function, it must return the previously mentioned matrix.
#' @param lean Logical, default is \code{FALSE}. Used to reduce the (memory) size of the summary object. If \code{TRUE}, then all objects of length N (the number of observations) are removed from the result. Note that some \code{fixest} methods may consequently not work when applied to the summary.
#' @param forceCovariance (Advanced users.) Logical, default is \code{FALSE}. In the peculiar case where the obtained Hessian is not invertible (usually because of collinearity of some variables), use this option to force the covariance matrix, by using a generalized inverse of the Hessian. This can be useful to spot where possible problems come from.
#' @param keepBounded (Advanced users -- \code{feNmlm} with non-linear part and bounded coefficients only.) Logical, default is \code{FALSE}. If \code{TRUE}, then the bounded coefficients (if any) are treated as unrestricted coefficients and their S.E. is computed (otherwise it is not).
#' @param n Integer, default is 1000. Number of coefficients to display when the print method is used.
#' @param ... Only used if the argument \code{.vocv} is provided and is a function: extra arguments to be passed to that function.
#'
#' @section Compatibility with \pkg{sandwich} package:
#' The VCOVs from \code{sandwich} can be used with \code{feols}, \code{feglm} and \code{fepois} estimations. If you want to have a \code{sandwich} VCOV when using \code{summary.fixest}, you can use the argument \code{vcov} to specify the VCOV function to use (see examples).
#' Note that if you do so and you use a formula in the \code{cluster} argument, an innocuous warning can pop up if you used several non-numeric fixed-effects in the estimation (this is due to the function \code{\link[stats]{expand.model.frame}} used in \code{sandwich}).
#'
#' @return
#' It returns a \code{fixest} object with:
#' \item{cov.scaled}{The new variance-covariance matrix (computed according to the argument \code{se}).}
#' \item{se}{The new standard-errors (computed according to the argument \code{se}).}
#' \item{coeftable}{The table of coefficients with the new standard errors.}
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest2]{femlm}}, \code{\link[fixest2]{feols}} or \code{\link[fixest2]{feglm}}. Use \code{\link[fixest2]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest2]{etable}} to visualize the results of multiple estimations.
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
#' est_pois <- fepois(Euros ~ log(dist_km) | Origin + Destination + Product, trade)
#'
#' # Comparing different types of standard errors
#' sum_standard <- summary(est_pois, vcov = "iid")
#' sum_hetero <- summary(est_pois, vcov = "hetero")
#' sum_oneway <- summary(est_pois, vcov = "cluster")
#' sum_twoway <- summary(est_pois, vcov = "twoway")
#'
#' etable(sum_standard, sum_hetero, sum_oneway, sum_twoway)
#'
#' # Alternative ways to cluster the SE:
#' summary(est_pois, vcov = cluster ~ Product + Origin)
#' summary(est_pois, vcov = ~ Product + Origin)
#' summary(est_pois, cluster = ~ Product + Origin)
#'
#' # You can interact the clustering variables "live" using the var1 ^ var2 syntax.#'
#' summary(est_pois, vcov = ~ Destination^Product)
#'
#' #
#' # Newey-West and Driscoll-Kraay SEs
#' #
#'
#' data(base_did)
#' # Simple estimation on a panel
#' est <- feols(y ~ x1, base_did)
#'
#' # --
#' # Newey-West
#' # Use the syntax NW ~ unit + time
#' summary(est, NW ~ id + period)
#'
#' # Now take a lag of 3:
#' summary(est, NW(3) ~ id + period)
#'
#' # --
#' # Driscoll-Kraay
#' # Use the syntax DK ~ time
#' summary(est, DK ~ period)
#'
#' # Now take a lag of 3:
#' summary(est, DK(3) ~ period)
#'
#' #--
#' # Implicit deductions
#' # When the estimation is done with a panel.id, you don't need to
#' # specify these values.
#'
#' est_panel <- feols(y ~ x1, base_did, panel.id = ~ id + period)
#'
#' # Both methods, NM and DK, now work automatically
#' summary(est_panel, "NW")
#' summary(est_panel, "DK")
#'
#' #
#' # VCOVs robust to spatial correlation
#' #
#'
#' data(quakes)
#' est_geo <- feols(depth ~ mag, quakes)
#'
#' # --
#' # Conley
#' # Use the syntax: conley(cutoff) ~ lat + lon
#' # with lat/lon the latitude/longitude variable names in the data set
#' summary(est_geo, conley(100) ~ lat + long)
#'
#' # Change the cutoff, and how the distance is computed
#' summary(est_geo, conley(200, distance = "spherical") ~ lat + long)
#'
#' # --
#' # Implicit deduction
#' # By default the latitude and longitude are directly fetched in the data based
#' # on pattern matching. So you don't have to specify them.
#' # Further an automatic cutoff is computed by default.
#'
#' # The following works
#' summary(est_geo, "conley")
#'
#'
#'
#' #
#' # Compatibility with sandwich
#' #
#'
#' # You can use the VOCVs from sandwich by using the argument .vcov:
#' library(sandwich)
#' summary(est_pois, .vcov = vcovCL, cluster = trade[, c("Destination", "Product")])
#'
summary.fixest <- function(object, vcov = NULL, cluster = NULL, ssc = NULL, .vcov = NULL,
                           stage = NULL, lean = FALSE, agg = NULL, forceCovariance = FALSE,
                           se = NULL, keepBounded = FALSE, n = 1000,
                           nthreads = getFixest_nthreads(), ...) {
  # computes the clustered SEs and returns the modified vcov and coeftable
  # NOTA: if the object is already a summary

  if (isTRUE(object$onlyFixef) || isTRUE(object$NA_model)) {
    # means that the estimation is done without variables
    return(object)
  }

  mc <- match.call()

  dots <- list(...)

  check_arg(n, "integer scalar GE{1}")
  if (!missing(n)) {
    object$n_print <- n
  }

  # we need this to save the summary flags
  # All three arguments se+cluster+.vcov are formatted into a valid vcov arg.
  vcov_in <- vcov <- oldargs_to_vcov(se, cluster, vcov, .vcov)

  check_arg(lean, "logical scalar")
  check_arg(stage, "NULL integer vector no na len(,2) GE{1} LE{2}")

  skip_vcov <- FALSE
  if (isTRUE(object$summary)) {
    if ("fromPrint" %in% names(dots)) {
      # From print
      return(object)
    } else if (is.null(vcov) && is.null(ssc)) {
      # We return directly the object ONLY if not any other argument has been passed

      skip_vcov <- TRUE
      if (is.null(agg) && is.null(stage)) {
        if (!lean || (lean && isTRUE(object$lean))) {
          # No modification required
          object$summary_from_fit <- FALSE
          return(object)
        }
      }
    }

    assign_flags(object$summary_flags, vcov = vcov, ssc = ssc, agg = agg)
  }

  # Checking arguments in ...
  if (is_user_level_call()) {
    if (!is_function_in_it(vcov)) {
      validate_dots(suggest_args = c("se", "cluster", "ssc"), valid_args = "dof")
    }
  }

  if (is.null(stage)) stage <- 2


  # IV
  if (isTRUE(object$iv) && !isTRUE(dots$iv)) {
    stage <- unique(stage)
    res <- list()

    # if lean, we still compute the summary for the first stage,
    #  then we will inject it in the iv_first_stage object of the 2nd stage
    # => this is critical to get the right Wald stat (of the 1st stage),
    #  otherwise it won't be possible to get it.
    remove_stage_1 <- FALSE
    if (lean && !1 %in% stage) {
      remove_stage_1 <- TRUE
      stage <- 1:2
    }

    stage_names <- c()

    for (s in seq_along(stage)) {
      if (stage[s] == 1) {
        for (i in seq_along(object$iv_first_stage)) {
          res[[length(res) + 1]] <- summary(object$iv_first_stage[[i]],
            vcov = vcov, ssc = ssc, lean = lean,
            forceCovariance = forceCovariance,
            n = n, nthreads = nthreads, iv = TRUE
          )

          stage_names[length(stage_names) + 1] <- paste0("First stage: ", names(object$iv_first_stage)[i])
        }
      } else {
        # We keep the information on clustering => matters for wald tests of 1st stage
        my_res <- summary(object,
          vcov = vcov, ssc = ssc, lean = lean,
          forceCovariance = forceCovariance, n = n, nthreads = nthreads, iv = TRUE
        )

        res[[length(res) + 1]] <- my_res
        stage_names[length(stage_names) + 1] <- "Second stage"
      }
    }

    if (lean && 2 %in% stage) {
      # we inject the summary of the first stage into the iv_first_stage
      qui_1st <- which(grepl("^First", stage_names))
      qui_2nd <- which(stage_names == "Second stage")

      tmp_1st <- res[qui_1st]
      names(tmp_1st) <- names(object$iv_first_stage)

      res[[qui_2nd]][["iv_first_stage"]] <- tmp_1st
    }

    if (remove_stage_1) {
      qui_2nd <- which(stage_names == "Second stage")
      return(res[[qui_2nd]])
    }

    if (length(res) == 1) {
      return(res[[1]])
    }

    index <- list("iv" = length(res))
    all_names <- list("iv" = stage_names)
    res_multi <- setup_multi(index, all_names, res)
    attr(res_multi, "print_request") <- "long"

    return(res_multi)
  }


  # The new VCOV
  if (skip_vcov) {
    vcov <- object$cov.scaled
  } else {
    vcov <- vcov.fixest(object,
      vcov = vcov, ssc = ssc, forceCovariance = forceCovariance,
      keepBounded = keepBounded, nthreads = nthreads,
      attr = TRUE, se = se, cluster = cluster, ...
    )
  }

  # NOTA:
  # I need to add se and cluster even if they're not needed only to ensure it
  # works fine when vcov is a function and cluster/se are arguments

  sd2 <- diag(vcov)
  sd2[sd2 < 0] <- NA
  se <- sqrt(sd2)

  # used to handle the case of bounded parameters
  params <- names(object$coefficients)
  if (length(se) != length(params)) {
    se <- se[params]
  }
  names(se) <- params

  # The coeftable is modified accordingly
  coeftable <- object$coeftable

  # th z & p values
  zvalue <- object$coefficients / se
  pvalue <- fixest_pvalue(object, zvalue, vcov)

  # update of se if bounded
  se_format <- se
  isBounded <- object$isBounded
  if (!is.null(isBounded) && any(isBounded)) {
    if (!keepBounded) {
      se_format[!isBounded] <- decimalFormat(se_format[!isBounded])
      se_format[isBounded] <- attr(isBounded, "type")
    }
  }

  # modifs of the table
  coeftable <- cbind(
    "Estimate" = object$coefficients, "Std. Error" = se_format,
    "t value" = zvalue, "Pr(>|t|)" = pvalue
  )

  attr(coeftable, "type") <- attr(se, "type") <- attr(vcov, "type")

  object$cov.scaled <- vcov
  object$coeftable <- coeftable
  object$se <- se

  if (lean) {
    var2clean <- c(
      "fixef_id", "residuals", "fitted.values", "scores", "sumFE",
      "slope_variables_reordered", "y", "weights", "irls_weights",
      "obs_selection", "iv_residuals", "fitted.values_demean",
      "working_residuals", "linear.predictors"
    )

    object[var2clean] <- NULL

    object$lean <- TRUE
  }

  object$summary <- TRUE

  # We save the arguments used to construct the summary
  if ("summary_flags" %in% names(dots)) {
    # If here => this is a call from an estimation (=fit)
    object$summary_flags <- dots$summary_flags
    object$summary_from_fit <- TRUE
  } else {
    # build_flags does not accept missing arguments
    if (missing(ssc)) ssc <- NULL

    if (lean) {
      size_KB <- as.numeric(object.size(vcov)) / 8 / 1000

      if (size_KB > 100) {
        # Here => means the user has manually provided a cluster => will be of size N at least
        # To respect lean = TRUE we keep no memory of this choice
        vcov_in <- NULL
      }
    }

    object$summary_flags <- build_flags(mc, vcov = vcov_in, ssc = ssc)
    object$summary_from_fit <- NULL
  }

  # agg
  if (!missnull(agg)) {
    agg_result <- aggregate(object, agg, full = TRUE, from_summary = TRUE)
    object$coeftable <- agg_result$coeftable
    object$model_matrix_info <- agg_result$model_matrix_info
    object$is_agg <- TRUE
  }

  return(object)
}


#' @rdname summary.fixest
summary.fixest_list <- function(object, se, cluster, ssc = getFixest_ssc(), .vcov, stage = 2, lean = FALSE, n, ...) {
  dots <- list(...)

  res <- list()
  for (i in seq_along(object)) {
    my_res <- summary(object[[i]], se = se, cluster = cluster, ssc = ssc, .vcov = .vcov, stage = stage, lean = lean, n = n)

    # we unroll in case of IV
    if ("fixest_multi" %in% class(my_res)) {
      data <- attr(my_res, "data")
      for (j in seq_along(data)) {
        res[[length(res) + 1]] <- data[[j]]
      }
    } else {
      res[[length(res) + 1]] <- my_res
    }
  }

  # We return a simple list
  class(res) <- NULL

  res
}

#' Summary method for fixed-effects coefficients
#'
#' This function summarizes the main characteristics of the fixed-effects coefficients. It shows the number of fixed-effects that have been set as references and the first elements of the fixed-effects.
#'
#' @method summary fixest.fixef
#'
#' @param object An object returned by the function \code{\link[fixest2]{fixef.fixest}}.
#' @param n Positive integer, defaults to 5. The \code{n} first fixed-effects for each fixed-effect dimension are reported.
#' @param ... Not currently used.
#'
#' @return
#' It prints the number of fixed-effect coefficients per fixed-effect dimension, as well as the number of fixed-effects used as references for each dimension, and the mean and variance of the fixed-effect coefficients. Finally, it reports the first 5 (arg. \code{n}) elements of each fixed-effect.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest2]{femlm}}, \code{\link[fixest2]{fixef.fixest}}, \code{\link[fixest2]{plot.fixest.fixef}}.
#'
#' @examples
#'
#' data(trade)
#'
#' # We estimate the effect of distance on trade
#' # => we account for 3 fixed-effects effects
#' est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination + Product, trade)
#'
#' # obtaining the fixed-effects coefficients
#' fe_trade <- fixef(est_pois)
#'
#' # printing some summary information on the fixed-effects coefficients:
#' summary(fe_trade)
#'
summary.fixest.fixef <- function(object, n = 5, ...) {
  # This function shows some generic information on the fixed-effect coefficients

  # checking arguments in dots
  if (is_user_level_call()) {
    validate_dots(suggest_args = "n")
  }

  Q <- length(object)
  fixef_names <- names(object)
  slope_flag <- grepl("\\[", fixef_names)
  fe <- gsub("\\[.+", "", fixef_names)
  slope <- gsub(".+\\[|\\].+", "", fixef_names)

  isSlope <- any(slope_flag)
  isFE <- any(!slope_flag)
  info <- as.character(10 * isFE + isSlope)

  # we rework the names
  fixef_names[slope_flag] <- paste0(slope[slope_flag], " (slopes: ", fe[slope_flag], ")")

  isRegular <- TRUE
  if (Q > 1) {
    nb_ref <- attr(object, "references")
    nb_per_cluster <- sapply(object, length)
    mean_per_cluster <- sd_per_cluster <- c()
    for (i in 1:Q) {
      mean_per_cluster[i] <- as.character(signif(mean(object[[i]]), 3))
      sd_per_cluster[i] <- as.character(signif(sd(object[[i]]), 3))
    }
    res <- as.data.frame(rbind(nb_per_cluster, nb_ref, mean_per_cluster, sd_per_cluster))

    row_1 <- paste0("Number of ", switch(info,
      "11" = "fixed-effects/slopes",
      "10" = "fixed-effects",
      "1" = "slopes"
    ))

    rownames(res) <- c(row_1, "Number of references", "Mean", "Standard-deviation")

    colnames(res) <- fixef_names

    if (sum(nb_ref) > Q - 1) {
      isRegular <- FALSE
    }
  }

  # The message

  my_title <- paste0(switch(info,
    "11" = "Fixed-effects/Slope",
    "10" = "Fixed_effects",
    "1" = "Slope"
  ), " coefficients\n")
  cat(my_title)
  if (Q == 1) {
    x1 <- object[[1]]
    if (slope_flag) {
      cat("Number of slope coefficients for variable ", slope, " (slope: ", fe, ") is ", length(x1), ".\n", sep = "")
    } else {
      cat("Number of fixed-effects for variable ", fixef_names, " is ", length(x1), ".\n", sep = "")
    }

    cat("\tMean = ", signif(mean(x1), 3), "\tVariance = ", signif(var(x1), 3), "\n", sep = "")
  } else {
    print(res)
  }

  # We print the first 5 elements of each fixed-effect
  cat("\nCOEFFICIENTS:\n")
  for (i in 1:Q) {
    m <- head(object[[i]], n)

    m_char <- as.data.frame(t(as.data.frame(c("", as.character(signif(m, 4))))))
    names(m_char) <- c(paste0(fixef_names[i], ":"), names(m))
    rownames(m_char) <- " "

    n_cluster <- length(object[[i]])
    if (n_cluster > n) {
      m_char[["   "]] <- paste0("... ", addCommas(n_cluster - n), " remaining")
    }

    print(m_char)
    if (i != Q) cat("-----\n")
  }
}

#' @rdname check_conv_feols
summary.fixest_check_conv <- function(object, type = "short", ...) {
  check_arg_plus(type, "match(short, detail)")
  if (is_user_level_call()) {
    validate_dots(suggest_args = "type")
  }

  if (type == "short") {
    info_max_abs <- lapply(object, function(x) sapply(x, function(y) max(abs(y))))

    info_max_abs <- do.call("rbind", info_max_abs)

    cat("Maximum absolute value of the first-order conditions:\n\n")
    print(info_max_abs)
  } else {
    extract_and_format <- function(x) {
      # x: list

      x_sorted <- lapply(x, sort)

      is_short <- sapply(x, function(y) length(y) <= 4)

      x_small <- list()
      for (i in seq_along(x)) {
        xi <- x_sorted[[i]]
        if (is_short[i]) {
          x_small[[i]] <- c(xi, rep(NA, 4 - length(xi)))
        } else {
          x_small[[i]] <- c(head(xi, 2), tail(xi, 2))
        }
      }

      x_mat <- format(do.call("rbind", x_small), digits = 3)

      x_fmt <- c()
      for (i in seq_along(x)) {
        xi <- x_mat[i, ]
        if (is_short[i]) {
          x_fmt[[i]] <- gsub(", +NA", "", paste0(xi, collapse = ", "))
        } else {
          x_fmt[[i]] <- paste0(xi[1], ", ", xi[2], ", ..., ", xi[3], ", ", xi[4])
        }
      }

      x_names <- sfill(names(x))

      res <- paste0(x_names, ": ", x_fmt)

      res
    }

    info <- lapply(object, extract_and_format)

    cat("Smallest and largest values of the first-order conditions:\n\n")
    var_names <- sfill(names(object), right = TRUE)
    intro <- sfill("|", n = nchar(var_names[1]) + 1)
    n_var <- length(object)
    for (i in 1:n_var) {
      cat(var_names[i], "\n")
      value <- paste0(intro, info[[i]])
      cat(value, sep = "\n")
      cat(gsub(" ", "_", intro), "\n")

      if (i < n_var) cat("\n")
    }
  }
}
