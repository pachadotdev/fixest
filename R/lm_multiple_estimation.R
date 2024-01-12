# Multiple estimation tools ----

multi_split <- function(env, fun) {
  split <- get("split", env)
  split.items <- get("split.items", env)
  split.id <- get("split.id", env)
  split.name <- get("split.name", env)

  assign("do_split", FALSE, env)

  n_split <- length(split.id)
  res_all <- vector("list", n_split)
  I <- 1
  for (i in split.id) {
    if (i == 0) {
      my_env <- reshape_env(env)
      my_res <- fun(env = my_env)
    } else {
      my_res <- fun(env = reshape_env(env, obs2keep = which(split == i)))
    }

    res_all[[I]] <- my_res
    I <- I + 1
  }

  values <- list(sample = split.items)

  # result
  res_multi <- setup_multi(res_all, values, var = split.name)

  return(res_multi)
}



multi_LHS_RHS <- function(env, fun) {
  do_multi_lhs <- get("do_multi_lhs", env)
  do_multi_rhs <- get("do_multi_rhs", env)

  assign("do_multi_lhs", FALSE, env)
  assign("do_multi_rhs", FALSE, env)

  nthreads <- get("nthreads", env)

  # IMPORTANT NOTE:
  # contrary to feols, the preprocessing is only a small fraction of the
  # computing time in ML models
  # Therefore we don't need to optimize processing as hard as in FEOLS
  # because the gains are only marginal

  fml <- get("fml", env)

  # LHS
  lhs_names <- get("lhs_names", env)
  lhs <- get("lhs", env)

  if (do_multi_lhs == FALSE) {
    lhs <- list(lhs)
  }

  # RHS
  if (do_multi_rhs) {
    rhs_info_stepwise <- get("rhs_info_stepwise", env)
    multi_rhs_fml_full <- rhs_info_stepwise$fml_all_full
    multi_rhs_fml_sw <- rhs_info_stepwise$fml_all_sw
    multi_rhs_cumul <- rhs_info_stepwise$is_cumul

    linear_core <- get("linear_core", env)
    rhs_sw <- get("rhs_sw", env)
  } else {
    multi_rhs_fml_full <- list(.xpd(rhs = fml[[3]]))
    multi_rhs_cumul <- FALSE
    linear.mat <- get("linear.mat", env)
    linear_core <- list(left = linear.mat, right = 1)
    rhs_sw <- list(1)
  }

  isLinear_left <- length(linear_core$left) > 1
  isLinear_right <- length(linear_core$right) > 1

  n_lhs <- length(lhs)
  n_rhs <- length(rhs_sw)
  res <- vector("list", n_lhs * n_rhs)

  rhs_names <- sapply(multi_rhs_fml_full, function(x) as.character(x)[[2]])

  for (i in seq_along(lhs)) {
    for (j in seq_along(rhs_sw)) {
      # reshaping the env => taking care of the NAs

      # Forming the RHS
      my_rhs <- linear_core[1]

      if (multi_rhs_cumul) {
        if (j > 1) {
          # if j == 1, the RHS is already in the data
          my_rhs[1 + (2:j - 1)] <- rhs_sw[2:j]
        }
      } else {
        my_rhs[2] <- rhs_sw[j]
      }

      if (isLinear_right) {
        my_rhs[[length(my_rhs) + 1]] <- linear_core$right
      }

      n_all <- lengths(my_rhs)
      if (any(n_all == 1)) {
        my_rhs <- my_rhs[n_all > 1]
      }

      if (length(my_rhs) == 0) {
        my_rhs <- 1
      } else {
        my_rhs <- do.call("cbind", my_rhs)
      }

      if (length(my_rhs) == 1) {
        is_na_current <- !is.finite(lhs[[i]])
      } else {
        is_na_current <- !is.finite(lhs[[i]]) | cpppar_which_na_inf_mat(my_rhs, nthreads)$is_na_inf
      }

      my_fml <- .xpd(lhs = lhs_names[i], rhs = multi_rhs_fml_full[[j]])

      if (any(is_na_current)) {
        my_env <- reshape_env(env, which(!is_na_current), lhs = lhs[[i]], rhs = my_rhs, fml_linear = my_fml)
      } else {
        # We still need to check the RHS (only 0/1)
        my_env <- reshape_env(env, lhs = lhs[[i]], rhs = my_rhs, fml_linear = my_fml, check_lhs = TRUE)
      }

      my_res <- fun(env = my_env)

      res[[index_2D_to_1D(i, j, n_rhs)]] <- my_res
    }
  }

  # Meta information for fixest_multi
  values <- list(
    lhs = rep(lhs_names, each = n_rhs),
    rhs = rep(rhs_names, n_lhs)
  )

  # result
  res_multi <- setup_multi(res, values)

  return(res_multi)
}


multi_fixef <- function(env, estfun) {
  # Honestly had I known it was so painful, I wouldn't have done it...
  assign("do_multi_fixef", FALSE, env)

  multi_fixef_fml_full <- get("multi_fixef_fml_full", env)
  combine.quick <- get("combine.quick", env)
  fixef.rm <- get("fixef.rm", env)
  family <- get("family", env)
  origin_type <- get("origin_type", env)
  nthreads <- get("nthreads", env)

  data <- get("data", env)

  n_fixef <- length(multi_fixef_fml_full)

  data_results <- vector("list", n_fixef)
  for (i in 1:n_fixef) {
    fml_fixef <- multi_fixef_fml_full[[i]]

    if (length(all.vars(fml_fixef)) > 0) {
      #
      # Evaluation of the fixed-effects
      #

      fixef_terms_full <- fixef_terms(fml_fixef)

      # fixef_terms_full computed in the formula section
      fixef_terms <- fixef_terms_full$fml_terms

      # FEs
      fixef_df <- error_sender(
        prepare_df(fixef_terms_full$fe_vars, data, combine.quick),
        "Problem evaluating the fixed-effects part of the formula:\n"
      )

      fixef_vars <- names(fixef_df)

      # Slopes
      isSlope <- any(fixef_terms_full$slope_flag != 0)
      slope_vars_list <- list(0)
      if (isSlope) {
        slope_df <- error_sender(
          prepare_df(fixef_terms_full$slope_vars, data),
          "Problem evaluating the variables with varying slopes in the fixed-effects part of the formula:\n"
        )

        slope_flag <- fixef_terms_full$slope_flag
        slope_vars <- fixef_terms_full$slope_vars
        slope_vars_list <- fixef_terms_full$slope_vars_list

        # Further controls
        not_numeric <- !sapply(slope_df, is.numeric)
        if (any(not_numeric)) {
          stop("In the fixed-effects part of the formula (i.e. in ", as.character(fml_fixef[2]), "), variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope_df)[not_numeric], "s.is.quote"), " not.")
        }

        # slope_flag: 0: no Varying slope // > 0: varying slope AND fixed-effect // < 0: varying slope WITHOUT fixed-effect
        onlySlope <- all(slope_flag < 0)
      }

      # fml update
      fml_fixef <- .xpd(rhs = fixef_terms)

      #
      # NA
      #

      for (j in seq_along(fixef_df)) {
        if (!is.numeric(fixef_df[[j]]) && !is.character(fixef_df[[j]])) {
          fixef_df[[j]] <- as.character(fixef_df[[j]])
        }
      }

      is_NA <- !complete.cases(fixef_df)

      if (isSlope) {
        # Convert to double
        who_not_double <- which(sapply(slope_df, is.integer))
        for (j in who_not_double) {
          slope_df[[j]] <- as.numeric(slope_df[[j]])
        }

        info <- cpppar_which_na_inf_df(slope_df, nthreads)
        if (info$any_na_inf) {
          is_NA <- is_NA | info$is_na_inf
        }
      }

      if (any(is_NA)) {
        # Remember that isFixef is FALSE so far => so we only change the reg vars
        my_env <- reshape_env(env = env, obs2keep = which(!is_NA))

        # NA removal in fixef
        fixef_df <- fixef_df[!is_NA, , drop = FALSE]

        if (isSlope) {
          slope_df <- slope_df[!is_NA, , drop = FALSE]
        }
      } else {
        my_env <- new.env(parent = env)
      }

      # We remove the linear part if needed

      if (get("do_multi_rhs", env)) {
        linear_core <- get("linear_core", my_env)
        if ("(Intercept)" %in% colnames(linear_core$left)) {
          int_col <- which("(Intercept)" %in% colnames(linear_core$left))
          if (ncol(linear_core$left) == 1) {
            linear_core$left <- 1
          } else {
            linear_core$left <- linear_core$left[, -int_col, drop = FALSE]
          }
          assign("linear_core", linear_core, my_env)
        }
      } else {
        linear.mat <- get("linear.mat", my_env)
        if ("(Intercept)" %in% colnames(linear.mat)) {
          int_col <- which("(Intercept)" %in% colnames(linear.mat))
          if (ncol(linear.mat) == 1) {
            assign("linear.mat", 1, my_env)
          } else {
            assign("linear.mat", linear.mat[, -int_col, drop = FALSE], my_env)
          }
        }
      }

      # We assign the fixed-effects
      lhs <- get("lhs", my_env)

      # We delay the computation by using isSplit = TRUE and split.full = FALSE
      # Real QUF will be done in the last reshape env
      info_fe <- setup_fixef(fixef_df = fixef_df, lhs = lhs, fixef_vars = fixef_vars, fixef.rm = fixef.rm, family = family, isSplit = TRUE, split.full = FALSE, origin_type = origin_type, isSlope = isSlope, slope_flag = slope_flag, slope_df = slope_df, slope_vars_list = slope_vars_list, nthreads = nthreads)

      fixef_id <- info_fe$fixef_id
      fixef_names <- info_fe$fixef_names
      fixef_sizes <- info_fe$fixef_sizes
      fixef_table <- info_fe$fixef_table
      sum_y_all <- info_fe$sum_y_all
      lhs <- info_fe$lhs

      obs2remove <- info_fe$obs2remove
      fixef_removed <- info_fe$fixef_removed
      message_fixef <- info_fe$message_fixef

      slope_variables <- info_fe$slope_variables
      slope_flag <- info_fe$slope_flag

      fixef_id_res <- info_fe$fixef_id_res
      fixef_sizes_res <- info_fe$fixef_sizes_res
      new_order <- info_fe$new_order

      assign("isFixef", TRUE, my_env)
      assign("new_order_original", new_order, my_env)
      assign("fixef_names", fixef_names, my_env)
      assign("fixef_vars", fixef_vars, my_env)

      assign_fixef_env(env, family, origin_type, fixef_id, fixef_sizes, fixef_table, sum_y_all, slope_flag, slope_variables, slope_vars_list)

      #
      # Formatting the fixef stuff from res
      #

      # fml & fixef_vars => other stuff will be taken care of in reshape
      res <- get("res", my_env)
      res$fml_all$fixef <- fml_fixef
      res$fixef_vars <- fixef_vars
      if (isSlope) {
        res$fixef_terms <- fixef_terms
      }
      assign("res", res, my_env)

      #
      # Last reshape
      #

      my_env_est <- reshape_env(my_env, assign_fixef = TRUE)
    } else {
      # No fixed-effect // new.env is indispensable => otherwise multi RHS/LHS not possible
      my_env_est <- reshape_env(env)
    }

    data_results[[i]] <- estfun(env = my_env_est)
  }

  index <- list(fixef = n_fixef)
  fixef_names <- sapply(multi_fixef_fml_full, function(x) as.character(x)[[2]])
  values <- list(fixef = fixef_names)

  res_multi <- setup_multi(data_results, values)

  if ("lhs" %in% names(attr(res_multi, "meta")$index)) {
    res_multi <- res_multi[lhs = TRUE]
  }

  return(res_multi)
}
