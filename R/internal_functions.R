# Internal functions ----

as.character.formula <- function(x, ...) as.character.default(x, ...)

parse_macros <- function(..., reset = FALSE, from_xpd = FALSE, check = TRUE, frame = NULL) {
  set_up(1)

  if (check) {
    check_arg(..., "dotnames os formula | character vector no na | numeric scalar | class(call, name)", .message = paste0("Each element of '...' must be a one-sided formula, and the name of each argument must start with two dots (ex: ", ifelse(from_xpd, "xpd(fml, ..ctrl = ~ x5 + x6)", "setFixest_fml(..ctrl = ~ x5 + x6)"), ").\nAlternatively it can be a character vector of variable names, or a numeric scalar."))
  }

  # Original macros
  fml_macro <- getOption("fixest_fml_macro")
  if (reset || is.null(fml_macro)) {
    fml_macro <- list()
  } else if (!is.list(fml_macro)) {
    warn_up("The value of getOption(\"fixest_fml_macro\") wasn't legal, it has been reset.")
    fml_macro <- list()
    options(fixest_fml_macro = list())
  }

  # We check the names
  if (...length() == 0) {
    return(fml_macro)
  }

  dots <- list(...)

  # Setting the new macros
  if (check) {
    qui_pblm <- !grepl("^\\.\\.", names(dots))
    if (any(qui_pblm)) {
      arg_name <- names(dots)[qui_pblm][1]
      correct_name <- gsub("^\\.*", "\\.\\.", arg_name)
      stop_up("Each argument name must start with two dots. Use '", correct_name, "' instead of '", arg_name, "'.")
    }
  }

  for (v in names(dots)) {
    fml_macro[[v]] <- value2stringCall(dots[[v]], check = check, frame = frame)
  }

  fml_macro
}

value2stringCall <- function(value_raw, call = FALSE, check = FALSE, frame = NULL) {
  if (any(c("call", "name") %in% class(value_raw))) {
    res <- if (call) value_raw else deparse_long(value_raw)
  } else if (inherits(value_raw, "formula")) {
    res <- if (call) value_raw[[2]] else as.character(value_raw)[[2]]
  } else {
    if (check) {
      value_raw <- grep("[[:alnum:]]", value_raw, value = TRUE)
      if (length(value_raw)) {
        # We need to check that it leads to a valid formula => otherwise problems later

        for (i in seq_along(value_raw)) {
          value_raw[i] <- dot_square_bracket(value_raw[i], frame)
        }

        value_raw <- paste(value_raw, collapse = " + ")

        my_call <- error_sender(str2lang(value_raw), "The value '", value_raw, "' does not lead to a valid formula: ", up = 2)
        res <- if (call) my_call else value_raw
      } else {
        res <- if (call) 1 else "1"
      }
    } else {
      value_raw <- value_raw[nzchar(value_raw)]
      if (length(value_raw)) {
        value_raw <- paste(value_raw, collapse = "+")
        res <- if (call) str2lang(value_raw) else value_raw
      } else {
        res <- if (call) 1 else "1"
      }
    }
  }

  if (check && "[" %in% all.vars(res, functions = TRUE)) {
    res_txt <- dot_square_bracket(deparse_long(res), frame)
    res <- as.formula(res_txt)
  }

  res
}

dot_square_bracket <- function(x, frame = .GlobalEnv, regex = FALSE, text = FALSE,
                               forced_merge = FALSE, up = 0) {
  # transforms "x.[i]" into x1 if i==1
  # z = "XX" ; x = ".[z] + x.[1:5] + y.[1:2]_t"
  # x = "x.[a] | fe[b] + m.[o]"
  # .[,stuff] => means aggregation is comma based
  # if text, we allow:
  # - before.["text", stuff]after AND .[stuff, "text"]
  #   * .["text", stuff] => collapses stuff with "text"
  #     ie paste0(before, paste0(stuff, collapse = "text"), after)
  #   * .["text", stuff] => collapses stuff WITH the previous string with "text",
  #     ie paste0(paste0(before, stuff, collapse = "text), after)
  # regex = FALSE ; frame = .GlobalEnv ; text = FALSE
  # x = 'c(.[,y]) ~ sw(x.[,1:3])'
  # name = c("Jon", "Juliet") ; x = "bonjour .[' and ', name]" ; text = TRUE
  # name = c("Jon", "Juliet") ; x = "bonjour .[name, ' and ']" ; text = TRUE

  if (!grepl(".[", x, fixed = TRUE)) {
    return(x)
  }

  x_split_open <- strsplit(x, ".[", fixed = TRUE)[[1]]

  # comma aggregation
  is_comma <- grepl("^,", x_split_open[-1])
  if (any(is_comma) && !text) {
    x_split_open[-1] <- gsub("^, *", "", x_split_open[-1])
  }

  # nesting: a ~ .["x.[1:2]_sq"] + b
  any_nested <- any(grepl("^\"", x_split_open[-1])) && !text
  if (any_nested) {
    i_open_quote <- setdiff(which(grepl("^\"", x_split_open)), 1)
    i_close_quote <- which(grepl("\"\\]", x_split_open))

    x_split_new <- x_split_open[1]
    for (i in i_open_quote) {
      j <- i_close_quote[i_close_quote >= i]

      xi_new_all <- x_split_open[i:j]
      xi_new_all <- gsub("\\]", "__close__", xi_new_all)

      xi_new <- paste0(xi_new_all, collapse = "__open__")
      xi_new <- gsub("\"__close__", "\"]", xi_new)

      x_split_new[[length(x_split_new) + 1]] <- xi_new
    }

    n <- length(x_split_open)
    if (j < n) {
      for (i in (j + 1):n) {
        x_split_new[[length(x_split_new) + 1]] <- x_split_open[i]
      }
    }

    x_split_open <- x_split_new
  }

  # Finding the elements left/right of the closing bracket
  # we take care of indexing within brackets (e.g. x1 + .[m[1]])
  b_open <- gregexpr("[", x_split_open[-1], fixed = TRUE)
  b_close <- gregexpr("]", x_split_open[-1], fixed = TRUE)

  x_split <- character(length(x_split_open) * 2 - 1)
  n <- length(x_split)
  x_split[1] <- x_split_open[1]

  for (i in seq_along(b_open)) {
    if (b_open[[i]][1] != -1) {
      # means there was an open [
      n_extend <- length(b_close[[i]]) - length(b_open[[i]])
      index_closing <- b_close[[i]][which.max(b_close[[i]] < c(b_open[[i]], rep(Inf, n_extend)))]
    } else {
      index_closing <- b_close[[i]][1]
    }

    x_SO <- x_split_open[i + 1]
    x_split_close_left <- substr(x_SO, 1, index_closing - 1)
    x_split_close_right <- substr(x_SO, index_closing + 1, nchar(x_SO))

    j <- 2 * i
    x_split[[j]] <- x_split_close_left
    x_split[[j + 1]] <- x_split_close_right
  }

  if (any_nested) {
    x_split <- gsub("__open__", ".[", x_split)
    x_split <- gsub("__close__", "]", x_split)
  }

  if (text) {
    # NA means no operation
    operator <- rep(NA_character_, n)
    do_split <- rep(FALSE, n)
  }

  res <- as.list(x_split)
  for (i in (1:n)[(1:n) %% 2 == 0]) {
    # catching the operation
    do_lang <- TRUE
    if (text) {
      x_txt <- trimws(x_split[i])

      op <- NULL
      is_agg <- is_split <- FALSE
      if (grepl("^,", x_txt)) {
        # default value of the aggregation
        op <- ""
        x_txt <- str_trim(x_txt, 1)
        is_agg <- TRUE
      } else if (grepl("^/", x_txt)) {
        # default value of the split
        op <- "@, *"
        x_txt <- str_trim(x_txt, 1)
        is_split <- TRUE
      }

      if (!is_split && !is_agg) {
        info_op <- regexpr("^(('[^']*')|(\"[^\"]*\"))[,/]", x_txt)
        if (info_op != -1) {
          arg_pos <- attr(info_op, "match.length")
          op <- substr(x_txt, 2, arg_pos - 2)
          arg <- substr(x_txt, arg_pos, arg_pos)
          x_txt <- str_trim(x_txt, arg_pos)

          is_agg <- arg == ","
          is_split <- !is_agg
        }
      }

      if (is_split) {
        do_lang <- FALSE

        operator[i] <- op
        do_split[i] <- TRUE
        my_call <- x_txt
      } else if (is_agg) {
        do_lang <- FALSE

        x_txt <- paste0("dsb_check_set_agg(", x_txt, ")")
        my_list <- list(dsb_check_set_agg = dsb_check_set_agg)
        my_call <- error_sender(eval(str2lang(x_txt), my_list, frame),
          "Dot square bracket operator: Evaluation of '.[",
          x_split[i], "]' led to an error:",
          up = up + 1
        )

        operator[i] <- op
      }
    }


    if (do_lang) {
      my_call <- error_sender(str2lang(x_split[i]),
        "Dot square bracket operator: Evaluation of '.[",
        x_split[i], "]' led to an error:",
        up = up + 1
      )
    }

    if (is.character(my_call) && grepl(".[", my_call, fixed = TRUE)) {
      # Nested call
      value <- .dsb(my_call, frame = frame)
    } else {
      # Informative error message
      value <- error_sender(as.character(eval(my_call, frame)),
        "Dot square bracket operator: Evaluation of '.[",
        x_split[i], "]' led to an error:",
        up = up + 1
      )

      if (length(value) == 2 && value[1] == "~") {
        value <- char_to_vars(value[2])
      }
    }


    res[[i]] <- value
  }

  if (any(is_comma)) {
    is_comma <- insert_in_between(FALSE, is_comma)
  } else if (!text) {
    is_comma <- rep(FALSE, length(res))
  }

  if (max(lengths(res)) == 1) {
    if (text && any(do_split)) {
      res_txt <- res[[1]]
      for (i in 2:length(res)) {
        if (do_split[i]) {
          res_txt <- paste0(res_txt, str_split(res[[i]], operator[i])[[1]])
        } else {
          res_txt <- paste0(res_txt, res[[i]])
        }
      }

      res <- res_txt
    } else {
      res <- paste(res, collapse = "")
    }
  } else {
    # first value is NEVER a vector ("" is added automatically in split)
    res_txt <- res[[1]]
    i <- 2
    while (i <= n) {
      if (length(res[[i]]) == 1) {
        if (text && do_split[i]) {
          res_txt <- paste0(res_txt, str_split(res[[i]], operator[i])[[1]])
        } else {
          res_txt <- paste0(res_txt, res[[i]])
        }

        i <- i + 1
      } else if (text) {
        if (is.na(operator[i])) {
          if (forced_merge) {
            res_txt <- paste0(res_txt, paste0(res[[i]], collapse = "_MERGE_"))
          } else {
            res_txt <- paste0(res_txt, res[[i]])
          }
        } else {
          # If we're here => must be an aggregation requested
          # (if it was a split, it would be of length 1)
          res_txt <- paste0(res_txt, paste0(res[[i]], collapse = operator[i]))
        }

        i <- i + 1
      } else if (regex) {
        after <- if (i != n) res[[i + 1]] else ""
        res_txt <- paste0(res_txt, res[[i]], after, collapse = "|")
        i <- i + 2
      } else {
        before_no_var <- gsub("[[:alnum:]_\\.]+$", "", res_txt)
        var_before <- substr(res_txt, nchar(before_no_var) + 1, nchar(res_txt))

        coll <- if (is_comma[i]) ", " else " + "

        if (i != n && grepl("[[:alnum:]_\\.]", substr(res[[i + 1]], 1, 1))) {
          after <- res[[i + 1]]
          after_no_var <- gsub("^[[:alnum:]_\\.]+", "", after)
          var_after <- substr(after, 1, nchar(after) - nchar(after_no_var))

          res_txt <- paste0(before_no_var, paste0(var_before, res[[i]], var_after, collapse = coll), after_no_var)
          i <- i + 2
        } else {
          res_txt <- paste0(before_no_var, paste0(var_before, res[[i]], collapse = coll))
          i <- i + 1
        }
      }
    }

    res <- res_txt
  }


  # Recursivity prevention: no recursivity if text = TRUE
  DSB_RECURSIVE <- TRUE
  is_rec <- exists("DSB_RECURSIVE", parent.frame(), inherits = FALSE)

  if (!text && !is_rec && grepl(".[", res, fixed = TRUE)) {
    res <- dot_square_bracket(res, frame, regex)
  }

  res
}


dsb_check_set_agg <- function(...) {
  # Here we should have only one element, everything has been cleaned up before

  if (...length() > 1) {
    if (...length() == 2) {
      stop_up("The operator .[] accepts only up to two elements, if so the first one MUST be a string literal (eg: .['text', y]). Problem: it is not a string literal.")
    }
    stop_up("You cannot have more than two elements in between .[]. Currently there are ", ...length() + 1, ".")
  }

  mc <- match.call()

  return(mc[[2]])
}

print_coeftable <- function(coeftable, lastLine = "", show_signif = TRUE) {
  # Simple function that does as the function coeftable but handles special cases
  # => to take care of the case when the coefficient is bounded

  if (!is.data.frame(coeftable)) {
    class(coeftable) <- NULL
    ct <- as.data.frame(coeftable)
  } else {
    ct <- coeftable
  }

  signifCode <- c("***" = 0.001, "** " = 0.01, "*  " = 0.05, ".  " = 0.1)

  pvalues <- ct[, 4]

  stars <- cut(pvalues, breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
  stars[is.na(stars)] <- ""

  whoIsLow <- !is.na(pvalues) & pvalues < 2.2e-16

  # Note that it's a bit different than format => I don't like xxe-yy numbers, very hard to read: you can't see large/small nbers at first sight
  for (i in 1:3) {
    ct[, i] <- decimalFormat(ct[, i])
  }

  ct[!whoIsLow, 4] <- format(ct[!whoIsLow, 4], digits = 5)

  ct[whoIsLow, 4] <- "< 2.2e-16"
  ct[is.na(ct[, 4]), 4] <- "NA"

  ct[, 5] <- stars
  names(ct)[5] <- ""

  print(ct)

  cat(lastLine)

  if (show_signif) {
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
}

plot_single_cluster <- function(x, n = 5, addExp = FALSE, fe_name, ...) {
  # It plots the n first and last most notable FEs

  isSlope <- grepl("\\[", fe_name)

  # we compare with the average of the coefficients
  if (isSlope == FALSE) x <- sort(x) - mean(x)

  x_name <- names(x)

  k <- length(x)
  nb_show <- min(k, 2 * n + 1)
  mid_value <- floor(nb_show / 2)
  xlim <- c(1, nb_show)
  xlim <- xlim + c(-1, 1) * diff(xlim) / 30

  if (k <= 2 * n + 1) {
    # no need of space to split the data
    x_values <- 1:k
    y_values <- x
    isSplit <- FALSE
  } else {
    nb_show <- nb_show - 1 # because we don't want to show the middle point
    x_values <- c(1:n, (n + 2):(2 * n + 1))
    y_values <- c(head(x, n), tail(x, n))
    isSplit <- TRUE
  }

  # very specific case where the axis confonds with the boxing
  ylim <- range(y_values)
  if (pretty(range(y_values))[1] < min(y_values) || tail(pretty(range(y_values)), 1) > max(y_values)) {
    ylim <- range(y_values) + c(-1, 1) * diff(range(y_values)) / 30
  }

  plot(x_values, y_values, ann = FALSE, axes = FALSE, xlim = xlim, ylim = ylim, col = 0)

  # display
  box()
  y <- axis(2)
  abline(h = y, lty = 4, col = "lightgrey")

  # name information & points
  points(x_values, y_values)
  text(head(x_values, mid_value), head(y_values, mid_value), head(x_name, mid_value), pos = 3)
  text(tail(x_values, nb_show - mid_value), tail(y_values, nb_show - mid_value), tail(x_name, nb_show - mid_value), pos = 1)

  if (isSlope) {
    fe <- gsub("\\[.+", "", fe_name)
    slope <- gsub(".+\\[|\\].+", "", fe_name)
    title(xlab = paste0("Slope Coefficients (", fe, ")"), line = 1)
    title(main = slope)
  } else {
    title(xlab = "Centered Fixed-Effects", line = 1)
    title(main = fe_name)
  }


  if (isSplit) {
    abline(v = c(n + 0.75, n + 1.25), lty = 2)

    axis(1, at = c(n + 0.75, n + 1.25), col = "grey", labels = NA, lwd.ticks = 0)
    axis(3, at = c(n + 0.75, n + 1.25), col = "grey", labels = NA, lwd.ticks = 0)
  }

  coord <- par("usr")
  # axis(1, at = coord[1:2], labels = c("coef", "exp(coef)"), lwd = 0, lwd.ticks = 1)

  if (addExp && !isSlope) {
    axis(1, at = coord[1], labels = "coef", lwd = 0, lwd.ticks = 1, line = -1)
    axis(1, at = coord[2], labels = "exp(coef)", lwd = 0, lwd.ticks = 1, line = -1)
    axis(4, at = y, labels = signif(exp(y), 2))
  }
}

prepare_matrix <- function(fml, base, fake_intercept = FALSE) {
  # This function is way faster than model.matrix but does not accept factors
  # The argument fml **MUST** not have factors!

  rhs <- if (length(fml) == 3) fml[c(1, 3)] else fml

  t <- terms(rhs, data = base)

  all_var_names <- attr(t, "term.labels")

  # We take care of interactions: references can be multiple, then ':' is legal
  all_vars <- cpp_colon_to_star(all_var_names)

  # Forming the call
  if (attr(t, "intercept") == 1 && !fake_intercept) {
    n <- nrow(base)
    all_vars_call <- str2lang(paste0("list('(Intercept)' = rep(1, ", n, "), ", paste0(all_vars, collapse = ", "), ")"))
    all_var_names <- c("(Intercept)", all_var_names)
  } else {
    all_vars_call <- str2lang(paste0("list(", paste0(all_vars, collapse = ", "), ")"))
  }

  # evaluation
  data_list <- eval(all_vars_call, base)

  # Handling the multi columns case (ex: bs(x1), splines)
  # NOTA: I need to add a check for i() because of 1 value interactions
  #       not caught by the != nber of obs
  qui_i <- grepl("\\bi\\(", all_var_names)
  if (any(lengths(data_list) != nrow(base)) || any(qui_i)) {
    all_n <- as.vector(lengths(data_list) / nrow(base))

    qui_pblm <- which(all_n %% 1 != 0)
    if (length(qui_pblm) > 0) {
      what <- data_list[[qui_pblm]]
      reason <- ifelse(is.null(nrow(what)), paste0("of length ", length(what)), paste0("with ", nrow(what), " rows"))

      stop("Evaluation of ", all_var_names[qui_pblm], " returns an object ", reason, " while the data set has ", nrow(base), " rows.", call. = FALSE)
    }

    all_n_vector <- rep(all_n, all_n)

    new_names <- as.list(all_var_names)
    for (i in which(all_n > 1 | qui_i)) {
      my_names <- colnames(data_list[[i]])
      if (is.null(my_names)) {
        my_names <- 1:all_n[i]
      }
      new_names[[i]] <- paste0(all_var_names[i], my_names)
    }

    all_var_names <- unlist(new_names)
  }

  res <- do.call("cbind", data_list)

  colnames(res) <- all_var_names

  res
}


fixest_model_matrix <- function(fml, data, fake_intercept = FALSE, i_noref = FALSE, mf = NULL) {
  # This functions takes in the formula of the linear part and the
  # data
  # It reformulates the formula (ie with lags and interactions)
  # then either apply a model.matrix
  # either applies an evaluation (which can be faster)
  #
  # fake_intercept => whether to add the intercept, only to make sure
  #  the factors are well created

  # fml = ~a*b+c+i(x1)+Temp:i(x2)+i(x3)/Wind

  if (length(fml) == 3) fml <- fml[c(1, 3)]

  # We need to check a^b otherwise error is thrown in terms()
  rhs_txt <- deparse_long(fml[[2]])
  if (grepl("\\^[[:alpha:]]", rhs_txt)) {
    stop("The special operator '^' can only be used in the fixed-effects part of the formula. Please use ':' instead.")
  }

  #
  # Evaluation
  #

  t_fml <- terms(fml)
  tl <- attr(t_fml, "term.labels")

  if (length(tl) == 0) {
    if (fake_intercept) {
      return(1)
    }

    res <- matrix(1, nrow = nrow(data), ncol = 1, dimnames = list(NULL, "(Intercept)"))
    return(res)
  }

  # We check for calls to i()
  qui_i <- grepl("(^|[^[:alnum:]_\\.])i\\(", tl)
  IS_I <- any(qui_i)
  if (IS_I) {
    # OMG... why do I always have to reinvent the wheel???
    is_intercept <- fake_intercept || (attr(t_fml, "intercept") == 1)
    i_naked <- which(is_naked_fun(tl[qui_i], "i"))

    if (i_noref) {
      for (i in seq_along(i_naked)) {
        j <- i_naked[i]
        txt <- gsub("(^|(?<=[^[:alnum:]\\._]))i\\(", "i_noref(", tl[qui_i][j], perl = TRUE)
        tl[qui_i][j] <- eval(str2lang(txt))
      }
    } else {
      for (i in seq_along(i_naked)) {
        if (!is_intercept && i == 1) next

        j <- i_naked[i]
        txt <- gsub("(^|(?<=[^[:alnum:]\\._]))i\\(", "i_ref(", tl[qui_i][j], perl = TRUE)
        tl[qui_i][j] <- eval(str2lang(txt))
      }
    }

    fml_no_inter <- .xpd(rhs = tl[!qui_i])

    if (!is_intercept) tl <- c("-1", tl)
    fml <- .xpd(rhs = tl)
  }

  # Are there factors NOT in i()? If so => model.matrix is used
  dataNames <- names(data)

  if (!is.null(mf)) {
    useModel.matrix <- TRUE
  } else {
    useModel.matrix <- FALSE
    if (IS_I) {
      linear.varnames <- all.vars(fml_no_inter[[2]])
      is_num <- sapply(data[, dataNames %in% linear.varnames, FALSE], is.numeric)
      if (length(is_num) > 0 && (any(!is_num) || grepl("factor", deparse_long(fml_no_inter)))) {
        useModel.matrix <- TRUE
      }
    } else {
      linear.varnames <- all.vars(fml[[2]])
      is_num <- sapply(data[, dataNames %in% linear.varnames, FALSE], is.numeric)
      if (length(is_num) == 0 || any(!is_num) || grepl("factor", deparse_long(fml))) {
        useModel.matrix <- TRUE
      }
    }
  }

  if (useModel.matrix) {
    # to catch the NAs, model.frame needs to be used....
    if (is.null(mf)) {
      mf <- stats::model.frame(fml, data, na.action = na.pass)
    }

    linear.mat <- stats::model.matrix(fml, mf)

    if (fake_intercept) {
      who_int <- which("(Intercept)" %in% colnames(linear.mat))
      if (length(who_int) > 0) {
        linear.mat <- linear.mat[, -who_int, drop = FALSE]
      }
    }
  } else {
    linear.mat <- prepare_matrix(fml, data, fake_intercept)
  }

  if (any(grepl("__CLEAN__", colnames(linear.mat), fixed = TRUE))) {
    new_names <- clean_interact_names(colnames(linear.mat))

    colnames(linear.mat) <- new_names
  }

  if (is.integer(linear.mat)) {
    linear.mat <- 1 * linear.mat
  }

  linear.mat
}


fixest_model_matrix_extra <- function(object, newdata, original_data, fml, fake_intercept = FALSE, i_noref = FALSE, subset = FALSE) {
  # Only used within model.matrix and predict
  # Overlay of fixest_model_matrix to take care of special things, eg:
  # - poly
  # - subset
  # - ?

  #
  # subset
  #

  # subset => we allow the extraction of only some variables
  all_vars <- all_vars_with_i_prefix(fml[[3]])
  if (isFALSE(subset)) {
    if (!original_data && any(!all_vars %in% names(newdata))) {
      pblm <- setdiff(all_vars, names(newdata))
      stop("In 'model.matrix', the variable", enumerate_items(pblm, "is.s.quote"), " in the formula but not in the argument 'data'. Use 'subset = TRUE' to enable the creation of partial data.")
    }
  } else {
    vars_keep <- names(newdata)

    if (is.character(subset)) {
      # ex: subset = c("x1$", "x2$")

      vars_keep <- keep_apply(vars_keep, subset)
      if (length(vars_keep) == 0) {
        stop("The variables in 'subset' do not match any variable in the 'data'.")
      }

      if (isFALSE(keep_apply("(Intercept)", subset, logical = TRUE))) {
        fake_intercept <- TRUE
      }
    } else {
      # intercept always removed if subset = TRUE!!!
      # that's the point of subset.

      fake_intercept <- TRUE
    }

    if (!all(all_vars %in% vars_keep)) {
      terms_all <- attr(terms(fml), "term.labels")
      # We first check pure variables (because non pure variables are slower to check)
      is_var <- grepl("^[\\.[:alpha:]][[:alnum:]\\._]*$", terms_all)
      terms_drop <- is_var & !terms_all %in% vars_keep

      for (i in which(!is_var)) {
        if (any(!all_vars_with_i_prefix(str2lang(terms_all[i])) %in% vars_keep)) {
          terms_drop[i] <- TRUE
        }
      }

      if (all(terms_drop)) {
        stop("Due to the use of the argument 'subset', not a single variable is left.")
      }

      fml <- .xpd(rhs = terms_all[!terms_drop])
    }
  }

  #
  # Extra functions that need raw data-evaluation + single valued factors
  #

  mf <- NULL
  if (!original_data) {
    # What I don't like here is:
    # to take care of only two particular cases, namely
    #  + functions using the original data
    #  + single value factors
    # we incur a lot of computing time to EVERYONE
    # That's really not good.
    # The problem is that it's not so easy to catch those cases 100% without error.
    # I shall find a solution at some point.

    # if lean = TRUE, we should be avoiding that
    # => I don't know of a solution yet...

    if (length(fml) == 3) fml <- fml[c(1, 3)]

    # We apply model.frame to the original data
    data <- fetch_data(object, "To apply 'model.matrix.fixest', ")

    panel__meta__info <- set_panel_meta_info(object, data)

    mf <- model.frame(fml, data, na.action = na.pass)

    rm(panel__meta__info) # needed, so that when the data is recreated for real

    t_mf <- terms(mf)
    xlev <- .getXlevels(t_mf, mf)

    if (!identical(attr(t_mf, "variables"), attr(t_mf, "predvars")) || length(xlev) > 0) {
      mf <- model.frame(t_mf, newdata, xlev = xlev, na.action = na.pass)
    } else {
      mf <- NULL
    }
  }

  GLOBAL_fixest_mm_info <- list()

  new_matrix <- fixest_model_matrix(fml, newdata, fake_intercept, i_noref, mf = mf)

  if (length(GLOBAL_fixest_mm_info) > 0) {
    attr(new_matrix, "model_matrix_info") <- GLOBAL_fixest_mm_info
  }

  new_matrix
}





fixef_terms <- function(fml, stepwise = FALSE, origin_type = "feols") {
  # separate all terms of fml into fixed effects and varying slopes
  # fml: one sided formula
  # can also be a vector of terms
  # Je fais tout ce tralala a cause de ce putain de terms() qui ne marche pas avec a^b!!!! fait chier!

  # fixef_terms(~dum_1[[x1]] + dum_2[x2])
  # fixef_terms(~dum_1[[x1]] + dum_2 + dum_3[x2, x3])

  if (!is.vector(fml)) {
    # we need to take care of the ^ used to combine variables
    fml_char <- as.character(fml[2])
    if (grepl("^", fml_char, fixed = TRUE)) {
      fml_char_new <- gsub("^", "%^%", fml_char, fixed = TRUE)
      fml <- as.formula(paste0("~", fml_char_new))
    }

    if (!origin_type %in% c("feols", "feglm") && grepl("[", fml_char, fixed = TRUE)) {
      stop("The use of varying slopes is available only for the functions feols, feglm or fepois.")
    }

    t <- terms(fml)

    if (stepwise) {
      sw_info <- extract_stepwise(tms = t)
      return(sw_info)
    }

    my_vars <- attr(t, "term.labels")
  } else {
    my_vars <- fml
  }

  # Further error checking (misuse of varying slopes)
  if (any(grepl("]", my_vars, fixed = TRUE))) {
    var2check <- my_vars[grepl("]", my_vars, fixed = TRUE)]
    var2check_double <- var2check[grepl("]]", var2check, fixed = TRUE)]
    var2check_single <- var2check[!grepl("]]", var2check, fixed = TRUE)]
    if (length(var2check_double) > 0) {
      qui <- !grepl("\\]\\]$", var2check_double) | !grepl("\\[\\[", var2check_double) | lengths(strsplit(var2check_double, "\\[")) != 3
      if (any(qui)) {
        item_pblm <- var2check_double[qui]
        msg <- paste0("Square bracket are special characters used **only** to designate varying slopes (see help). They are currenlty misused (it concerns ", enumerate_items(item_pblm), ").")
        class(msg) <- "try-error"
        return(msg)
      }
    }

    if (length(var2check_single) > 0) {
      qui <- !grepl("\\]$", var2check_single) | lengths(strsplit(var2check_single, "\\[")) != 2
      if (any(qui)) {
        item_pblm <- var2check_single[qui]
        msg <- paste0("Square bracket are special characters used **only** to designate varying slopes (see help). They are currenlty misused (it concerns ", enumerate_items(item_pblm), ").")
        class(msg) <- "try-error"
        return(msg)
      }
    }
  }

  # And yet again some error checking => i() should NOT be used
  if (any(grepl("^i\\(", my_vars))) {
    # We create an error instead of simply correcting the syntax => this is because the function i is
    # very different and should not be confused

    var_pblm <- my_vars[grepl("^i\\(", my_vars)][1]

    get_new_var <- function(var, f, f2, ...) match.call()

    what <- eval(str2lang(gsub("^i", "get_new_var", var_pblm)))
    n_var <- sum(c("var", "f", "f2") %in% names(what))
    msg <- if (n_var == 1) "Using i() to create fixed-effects is not possible, use directly the variable." else paste0("To interact fixed-effects, use the syntax fe1^fe2 (in your case ", deparse(what[[2]]), "^", deparse(what[[3]]), ").")

    stop("The function i() should not be used in the fixed-effects part of the formula. ", msg)
  }

  # Internal function
  reformulate_varslope <- function(dum, ..., add_dum = TRUE) {
    mc <- match.call()
    dum <- deparse_long(mc[["dum"]])

    if (add_dum) {
      res <- dum
    } else {
      res <- c()
    }

    for (i in 3:(length(mc) - !add_dum)) {
      value <- deparse_long(mc[[i]])
      res <- c(res, paste0(dum, "[[", value, "]]"))
    }
    res
  }

  qui_slope <- grepl("[", my_vars, fixed = TRUE)
  if (any(qui_slope)) {
    new_terms <- c()

    vars_slope <- my_vars[qui_slope]
    for (i in seq_along(my_vars)) {
      if (qui_slope[i]) {
        v <- my_vars[i]
        v_mod <- paste0("reformulate_varslope(", gsub("\\[+", ", ", v))
        add_dum <- ifelse(grepl("\\[\\[", v), ", add_dum = FALSE", "")
        v_mod <- gsub("\\]+", paste0(add_dum, ")"), v_mod)

        new_v <- eval(str2lang(v_mod))
        new_terms <- c(new_terms, new_v)
      } else {
        new_terms <- c(new_terms, my_vars[i])
      }
    }

    my_vars <- new_terms
  }

  # we return fml_terms, fe_vars, slope_vars, slope_flag

  # we need to unique it
  my_vars <- unique(my_vars)

  # we put the ^ back
  my_vars <- gsub(" %^% ", "^", my_vars, fixed = TRUE)

  res <- list(fml_terms = my_vars)
  # OLD version
  fe_vars_all <- gsub("\\[.+", "", my_vars)
  is_slope_all <- grepl("\\[", my_vars)
  slope_vars_all <- rep(NA_character_, length(my_vars))
  slope_vars_all[is_slope_all] <- gsub(".+\\[|\\]", "", my_vars[is_slope_all])

  # New version following the reworking of demeaning
  fe_vars <- unique(fe_vars_all)
  # fe_vars: unique of the fixed-effects variables
  slope_flag <- rep(0, length(fe_vars))
  # slope_flag: 0: no Varying slope // > 0: varying slope AND fixed-effect // < 0: varying slope WITHOUT fixed-effect

  slope_vars_list <- vector("list", length(fe_vars))
  names(slope_vars_list) <- fe_vars
  for (i in seq_along(fe_vars)) {
    qui <- which(fe_vars_all == fe_vars[i])
    if (any(is_slope_all[qui])) {
      nb <- sum(is_slope_all[qui])
      slope_flag[i] <- nb * (1 - 2 * (length(qui) == nb))
      slope_vars_list[[i]] <- slope_vars_all[qui][is_slope_all[qui]]
    }
  }

  res$fe_vars <- fe_vars
  res$slope_flag <- as.integer(slope_flag)
  res$slope_vars_list <- slope_vars_list
  res$slope_vars <- unlist(slope_vars_list, use.names = FALSE)

  res
}

prepare_df <- function(vars, base, fastCombine = NA) {
  # vars: vector of variables to evaluate

  # we drop NAs and make it unique
  vars <- unique(vars[!is.na(vars)])
  all_var_names <- vars

  do_combine <- !is.na(fastCombine)

  changeNames <- FALSE
  if (do_combine && any(grepl("^", vars, fixed = TRUE))) {
    # special indicator to combine factors
    # ^ is a special character: only used to combine variables!!!

    vars <- fml_combine(vars, fastCombine, vars = TRUE)

    changeNames <- TRUE
  }

  all_vars <- cpp_colon_to_star(vars)

  if (all(all_vars %in% names(base))) {
    res <- base[, all_vars, drop = FALSE]
  } else {
    all_vars_call <- str2lang(paste0("list(", paste0(all_vars, collapse = ", "), ")"))
    data_list <- try(eval(all_vars_call, base))

    # if error: we send it back to the main function
    if ("try-error" %in% class(data_list)) {
      return(data_list)
    }

    names(data_list) <- all_var_names
    data_list$stringsAsFactors <- FALSE

    res <- do.call("data.frame", data_list)

    if (changeNames) {
      all_var_names <- rename_hat(all_var_names)
    }

    names(res) <- all_var_names
  }

  res
}


# fml_char = "x + y + u^factor(v1, v2) + x5"
fml_combine <- function(fml_char, fastCombine, vars = FALSE) {
  # function that transforms "hat" interactions into a proper function call:
  # Origin^Destination^Product + Year becomes ~combine_clusters(Origin, Destination, Product) + Year

  fun2combine <- ifelse(fastCombine, "combine_clusters_fast", "combine_clusters")

  # we need to change ^ into %^% otherwise terms sends error
  labels <- attr(terms(.xpd(rhs = gsub("\\^(?=[^0-9])", "%^%", fml_char, perl = TRUE))), "term.labels")

  # now we work this out
  for (i in seq_along(labels)) {
    lab <- labels[i]
    if (grepl("^", lab, fixed = TRUE)) {
      lab_split <- trimws(strsplit(lab, "%^%", fixed = TRUE)[[1]])
      if (grepl("(", lab, fixed = TRUE)) {
        # we add some error control -- imperfect, but... it's enough
        lab_collapsed <- gsub("\\([^\\)]+\\)", "", lab)
        if (length(lab_split) != length(strsplit(lab_collapsed, "%^%", fixed = TRUE)[[1]])) {
          msg <- "Wrong formatting of the fixed-effects interactions. The '^' operator should not be within parentheses."
          stop(msg)
        }
      }
      labels[i] <- paste0(fun2combine, "(", paste0(lab_split, collapse = ", "), ")")
    }
  }

  if (vars) {
    return(labels)
  }

  fml <- .xpd(rhs = labels)

  fml
}

# x = c('combine_clusters(bin(fe1, "!bin::2"), fe2)', 'fe3')
rename_hat <- function(x) {
  qui <- grepl("combine_clusters", x, fixed = TRUE)
  if (!any(qui)) {
    return(x)
  }

  for (i in which(qui)) {
    xi_new <- gsub("combine_clusters(_fast)?", "sw", x[i])
    sw_only <- extract_fun(xi_new, "sw")
    sw_eval <- eval(str2lang(sw_only$fun))
    new_var <- paste0(sw_only$before, paste0(sw_eval, collapse = "^"), sw_only$after)
    x[i] <- new_var
  }

  x
}

prepare_cluster_mat <- function(fml, base, fastCombine) {
  # prepares the data.frame of the cluster variables

  fml_char <- as.character(fml[2])
  changeNames <- FALSE
  if (grepl("^", fml_char, fixed = TRUE)) {
    # special indicator to combine factors
    fml <- fml_combine(fml_char, fastCombine)
    changeNames <- TRUE
  }

  t <- terms(fml, data = base)

  all_var_names <- attr(t, "term.labels")
  all_vars <- gsub(":", "*", all_var_names)

  if (all(all_vars %in% names(base))) {
    res <- base[, all_vars, drop = FALSE]
  } else {
    all_vars_call <- str2lang(paste0("list(", paste0(all_vars, collapse = ", "), ")"))
    data_list <- eval(all_vars_call, base)
    names(data_list) <- all_var_names
    data_list$stringsAsFactors <- FALSE

    res <- do.call("data.frame", data_list)

    if (changeNames) {
      all_var_names <- rename_hat(all_var_names)
    }

    names(res) <- all_var_names
  }

  res
}

combine_clusters_fast <- function(...) {
  # This functions creates a new cluster from several clusters
  # basically: paste(cluster1, cluster2, ... etc, sep = "__")
  # but it tries to be faster than that because paste is very slow on large datasets

  cluster <- list(...)
  Q <- length(cluster)

  # The problem is that the clusters can be of ANY type...
  # So too much of a pain to take care of them in c++
  # => when I have time I'll do it, but pfff...

  # Not super efficient, but that's life
  ANY_NA <- FALSE
  if (any(who_NA <- sapply(cluster, anyNA))) {
    ANY_NA <- TRUE
    who_NA <- which(who_NA)
    IS_NA <- is.na(cluster[[who_NA[1]]])
    for (i in who_NA[-1]) {
      IS_NA <- IS_NA | is.na(cluster[[i]])
    }

    # Nice error message comes later
    if (all(IS_NA)) {
      return(rep(NA, length(IS_NA)))
    }

    # we recreate the clusters
    for (i in 1:Q) {
      cluster[[i]] <- cluster[[i]][!IS_NA]
    }
  }

  # First we unclass
  for (i in 1:Q) {
    cluster[[i]] <- quickUnclassFactor(cluster[[i]])
  }

  # Then we combine
  power <- floor(1 + log10(sapply(cluster, max)))

  if (sum(power) > 14) {
    order_index <- do.call(order, cluster)
    index <- cpp_combine_clusters(cluster, order_index)
  } else {
    # quicker, but limited by the precision of doubles
    index <- cluster[[1]]
    for (q in 2:Q) {
      index <- index + cluster[[q]] * 10**sum(power[1:(q - 1)])
    }
  }

  if (ANY_NA) {
    # we recreate the return vector with appropriate NAs
    res <- rep(NA_real_, length(IS_NA))
    res[!IS_NA] <- index
  } else {
    res <- index
  }

  return(res)
}

combine_clusters <- function(...) {
  # This functions creates a new cluster from several clusters
  # basically: paste(cluster1, cluster2, ... etc, sep = "_")
  # No, it's mow much more fatser than that!!!! Thanks for my new algo
  # => the paste only applies to the unique number of items


  clusters <- to_integer(..., add_items = TRUE, items.list = TRUE, internal = TRUE)

  res <- clusters$items[clusters$x]

  return(res)
}

listDefault <- function(x, variable, value) {
  # This function puts 'value' into the element 'variable' of list 'x'
  # IF it does not already exists in 'x'

  x_name <- deparse(substitute(x))

  if (is.null(x[[variable]])) {
    x[[variable]] <- value

    assign(x_name, x, envir = parent.frame(n = 1))
  }
}

hgrid <- function(lty = 3, col = "darkgray", ymin = -Inf, ymax = Inf, ...) {
  # simple function that draws an horizontal grid

  # Finding the coordinates
  y <- axis(2, lwd = 0, labels = NA)

  y <- y[y > ymin & y < ymax]

  # now drawing the lines
  if (length(y) > 0) {
    abline(h = y, col = col, lty = lty, ...)
  }
}

vgrid <- function(lty = 3, col = "darkgray", xmin = -Inf, xmax = Inf, ...) {
  # simple function that draws an horizontal grid

  # Finding the coordinates
  x <- axis(1, lwd = 0, labels = NA)

  x <- x[x > xmin & x < xmax]

  # now drawing the lines
  if (length(x) > 0) {
    abline(v = x, col = col, lty = lty, ...)
  }
}

shade_area <- function(y1, y2, x, xmin, xmax, col = "grey", ...) {
  # fonction plus pratique que polygon
  # elle permet de griser une partie delimitee par
  # y1 et y2 pour chacune des valeurs de x
  # on doit avoir la meme longueur de y1,y2 et x
  # exemple:
  # a=curve(x**2,-5,5)
  # shade_area(a$y+1,a$y-1,a$x)
  # qqes parametres graphiques:
  # lwd / border (couleur du bord, peut etre NA) / lty

  n <- length(x)
  stopifnot(length(y1) == n | length(y1) == 1)
  stopifnot(length(y2) == n | length(y2) == 1)

  if (length(y1) == 1) y1 <- rep(y1, n)
  if (length(y2) == 1) y2 <- rep(y2, n)

  if (missing(xmin)) xmin <- min(x)
  if (missing(xmax)) xmax <- max(x)

  ind <- which(x >= xmin & x <= xmax)
  x1 <- x[ind]
  x2 <- x[rev(ind)]
  polygon(c(x1, x2), c(y1[ind], y2[rev(ind)]), col = col, ...)
}

clean_interact_names <- function(x) {
  # GOD, why is it so complicated? Why?
  # Why do I have to spend so much time on that crap?
  #
  # x = c("i(Month, keep = 4:7)__CLEAN__Month::5", "i(Month):Wind__CLEAN__Month::5", "Temp:i(bonjour(test), drop = 5:12):Wind__CLEAN__bonjour::5", "i(a, pi(b))__CLEAN__a::o")
  #
  # Speed consideration:
  # for 2000 names: 10ms simple case
  #                 28ms complex case
  # price to pay is OK

  res <- x

  who2clean <- grepl("__CLEAN__", x, fixed = TRUE)

  x2clean <- x[who2clean]

  x_split <- strsplit(x2clean, "(^|(?<=[^[:alnum:]\\._]))(i|sunab(_att)?)\\(|__CLEAN__", perl = TRUE)

  x_left <- sapply(x_split, function(v) v[1])
  x_right <- sapply(x_split, function(v) v[2])
  x_alias <- sapply(x_split, function(v) v[3])

  x_right <- gsub("^[^\\(\\)]+(\\(|\\))", "\\1", x_right)

  qui_ok <- substr(x_right, 1, 1) == ")"

  x_right[qui_ok] <- substr(x_right[qui_ok], 2, 500)

  if (any(!qui_ok)) {
    clean_x_right_hard <- function(letter_vec) {
      open <- 1 + cumsum(letter_vec == "(")
      close <- cumsum(letter_vec == ")")
      letter_vec <- letter_vec[-(1:which.max(close - open == 0))]
      paste(letter_vec, collapse = "")
    }

    x_right[!qui_ok] <- sapply(strsplit(x_right[!qui_ok], ""), clean_x_right_hard)
  }

  res[who2clean] <- paste0(x_left, x_alias, x_right)

  return(res)
}

is_naked_fun <- function(x, fun_pattern) {
  # Why is it always so complicated... There must be an easier way
  # x = c("i(x1)", "i(I(x3))", "i(x3, x4, TRUE, drop = c(1, 3:5))", "Temp:i(x2)", "i(x3):Wind")

  x_split <- strsplit(x, paste0("(^|(?<=[^[:alnum:]\\._]))", fun_pattern, "\\("), perl = TRUE)

  left <- sapply(x_split, function(x) x[1])
  right <- sapply(x_split, function(x) x[2])

  right_naked <- function(r) {
    if (!grepl("(", r, fixed = TRUE) && grepl("\\)$", r)) {
      return(TRUE)
    }

    letter_vec <- strsplit(r, "")[[1]]
    open <- 1 + cumsum(letter_vec == "(")
    close <- cumsum(letter_vec == ")")
    which.max(close - open == 0) == length(letter_vec)
  }

  left_ok <- nchar(left) == 0
  right_ok <- rep(FALSE, length(right))
  if (any(left_ok)) {
    right_ok[left_ok] <- sapply(right[left_ok], right_naked)
  }

  left_ok & right_ok
}

set_defaults <- function(opts_name) {
  opts <- getOption(opts_name)
  if (is.null(opts) || length(opts) == 0) {
    return(NULL)
  }

  sysOrigin <- sys.parent()
  mc <- match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin), expand.dots = FALSE)
  args_in <- names(mc)

  args2set <- setdiff(names(opts), args_in)

  for (v in args2set) {
    assign(v, opts[[v]], parent.frame())
  }
}

fetch_data <- function(x, prefix = "", suffix = "") {
  # x: fixest estimation
  # We try different strategies:
  # 1) using the environment where the estimation was done
  # 2) the "parent.frame()" defined as the frame on top of ALL fixest functions
  # 3) the global environment, if it wasn't in 1)

  # Maybe I should keep only 1) => is there a reason to add the others?

  # 1) safest
  # 2) less safe but OK => note ???
  # 3) kind of dangerous => warning() ???

  if (is.null(x$call$data)) {
    return(NULL)
  }

  # 1) Environment of the call

  data <- NULL
  try(data <- eval(x$call$data, x$call_env), silent = TRUE)

  if (!is.null(data)) {
    return(data)
  }

  # 2) First non fixest frame

  fixest_funs <- ls(getNamespace("fixest2"))

  i <- 2
  sysOrigin <- sys.parent(i)
  while (sysOrigin != 0 && as.character(sys.call(sysOrigin)[[1]]) %in% fixest_funs) {
    i <- i + 1
    sysOrigin <- sys.parent(i)
  }

  if (sysOrigin != 0) {
    # We try again...
    try(data <- eval(x$call$data, parent.frame(sysOrigin)), silent = TRUE)

    if (!is.null(data)) {
      return(data)
    }
  }

  # 3) Global environment

  if (!identical(parent.env(x$call_env), .GlobalEnv)) {
    # ...and again
    try(data <- eval(x$call$data, .GlobalEnv), silent = TRUE)

    if (!is.null(data)) {
      return(data)
    }
  }

  # => Error message

  if (nchar(prefix) == 0) {
    msg <- "W"
  } else {
    s <- ifelse(grepl(" $", prefix), "", " ")
    if (grepl("\\. *$", prefix)) {
      msg <- paste0(prefix, s, "W")
    } else {
      msg <- paste0(prefix, s, "w")
    }
  }

  if (nchar(prefix) == 0) {
    msg <- "W"
  } else if (grepl("\\. *$", prefix)) {
    msg <- paste0(gsub(" +$", "", prefix), " W")
  } else {
    msg <- paste0(gsub(prefix, " +$", ""), " w")
  }

  if (nchar(suffix) > 0) {
    suffix <- gsub("^ +", "", suffix)
  }

  stop_up(msg, "e fetch the data in the enviroment where the estimation was made, but the data does not seem to be there any more (btw it was ", charShorten(deparse(x$call$data)[1], 15), "). ", suffix)
}

is_large_object <- function(x) {
  if (is.list(x)) {
    if (length(x[[1]]) > 10000) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if (length(x)[1] > 10000) {
    return(TRUE)
  }

  FALSE
}

build_flags <- function(mc, ..., call_env) {
  # Returns a list of arguments
  # If the arguments are too large => we use calls instead and save the calling environment
  # Note that the arguments that will be passed in ... must NOT be missing // they must be initialized to NULL
  # All that stuff just to avoid too large objects....

  dots <- list(...)

  args <- names(dots)

  res <- list()
  for (i in seq_along(dots)) {
    x <- dots[[i]]
    if (is.null(x)) {
      # nothing // NULL arguments are still in 'dots', so we want to avoid them
    } else if (is_large_object(x)) {
      res[[args[i]]] <- mc[[args[i]]]
      if (is.null(res$call_env)) {
        if (missing(call_env)) {
          res$call_env <- new.env(parent = parent.frame(2))
        } else {
          res$call_env <- call_env
        }
      }
    } else {
      res[[args[i]]] <- x
    }
  }

  res
}

assign_flags <- function(flags, ...) {
  # flags == arguments
  # We get the arguments used when the function was originally called
  # we assign these arguments to the previous frame

  dots <- list(...)

  res <- list()

  for (arg in setdiff(names(flags), "call_env")) {
    if (is.null(dots[[arg]])) {
      # we only assign values that were not provided by the user
      x <- flags[[arg]]
      if (is.name(x) || is.call(x)) {
        x <- eval(x, flags$call_env)
      }
      assign(arg, x, parent.frame())
    }
  }
}

items_to_drop <- function(items, x, varname, keep = FALSE, argname,
                          keep_first = FALSE, up = 1, no_error = FALSE,
                          valid_ref = FALSE) {
  # selection of items
  # the selection depends on the type of x
  # always returns the IDs of the items to drop

  set_up(up)

  if (missing(argname)) {
    argname <- deparse(substitute(x))
  }

  ref <- argname == "ref" && !valid_ref

  if (keep_first && (keep || ref)) {
    stop("Internal error: keep_first should not be used with 'keep' or 'ref'.")
  }

  if (is.factor(items)) {
    items <- as.character(items)
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }


  if (inherits(x, "formula")) {
    if (!identical(all.vars(x), "x")) {
      stop_up("In argument '", argname, "', the formula must contain a single variable name: 'x'. So far '", deparse_long(x), "' is not valid.")
    }

    if (length(x) > 2) {
      stop_up("In argument '", argname, "', if a formula, it must be one-sided. Problem: '", deparse_long(x), "' is two-sided.")
    }

    is_here <- error_sender(
      eval(x[[2]], list(x = items)),
      "In argument '", argname, "', the evaluation of the formula led to an error:"
    )
    if (length(is_here) != length(items)) {
      stop_up("In argument '", argname, "', the evaluation of the formula must return a logical vector of the same length as 'x'. Problem: '", deparse_long(x), "' returns a vector of length ", length(is_here), " (expected: ", length(items), ").")
    }

    if (!is.logical(is_here)) {
      stop_up("In argument '", argname, "', the evaluation of the formula must return a logical vector. Problem: '", deparse_long(x), "' is not logical (instead it is of class ", enumerate_items(class(is_here)), ").")
    }

    is_here <- !is.na(is_here) & is_here

    if (!no_error && !any(is_here)) {
      stop_up("In argument '", argname, "', the evaluation of the formula must match at least one value. Problem: '", deparse_long(x), "' does not match any.")
    }

    if (keep) {
      id_drop <- which(!is_here)
    } else {
      id_drop <- which(is_here)
    }
  } else if (is.character(x)) {
    all_x <- c()
    for (i in seq_along(x)) {
      my_x <- x[i]
      if (grepl("^@", my_x)) {
        # A) regex
        pattern <- substr(my_x, 2, nchar(my_x))
        new_x <- grep(pattern, items, value = TRUE)
        if (length(new_x) == 0 && !no_error) {
          # strong checking!
          stop_up("In argument '", argname, "', the regular expression '", pattern, "' does not match any value of '", varname, "'.")
        }
        all_x <- c(all_x, new_x)
      } else {
        # B) partial matching
        if (no_error) {
          qui_ok <- pmatch(tolower(my_x), tolower(items), duplicates.ok = TRUE)
          my_x <- items[qui_ok]
          my_x <- my_x[!is.na(my_x)]
        } else {
          check_value_plus(my_x, "match", .choices = unique(items), .message = paste0("The argument '", argname, "' should contain values of the variable '", varname, "'."))
        }

        all_x <- c(all_x, my_x)
      }

      if (keep_first && i == 1) {
        all_x_first <- all_x
      }
    }

    if (ref) {
      if (length(all_x) == 1) {
        id_drop <- which(items %in% all_x)
      } else {
        id_drop <- c(which(items %in% all_x[1]), which(items %in% all_x[-1]))
      }
    } else if (keep) {
      id_drop <- which(!items %in% all_x)
    } else {
      id_drop <- which(items %in% all_x)

      if (keep_first) {
        id_drop <- unique(c(which(items %in% all_x_first), id_drop))
      }
    }
  } else {
    # exact matching
    if (keep) {
      id_drop <- which(!items %in% x)
    } else {
      if (ref) {
        if (length(x) == 1) {
          id_drop <- which(items %in% x)
        } else {
          id_drop <- c(which(items %in% x[1]), which(items %in% x[-1]))
        }
      } else {
        id_drop <- which(items %in% x)

        if (keep_first) {
          if (!x[1] %in% items && !no_error) {
            stop_up("In argument '", argname, "', the value ", x[1], " does not match any value of '", varname, "'.")
          }

          id_drop <- unique(c(which(items %in% x[1]), id_drop))
        }
      }

      if (length(id_drop) == 0 && !no_error) {
        stop_up("In argument '", argname, "', the value", plural_len(x, "s.don't"), " match any value of '", varname, "'.")
      }
    }
  }

  id_drop
}



bin_factor <- function(bin, x, varname, no_error = FALSE) {
  # x = base_did$period
  # bin = list("9" = c(9, 2, 3), "7" = "@7|8", 5:6)
  # varname = "period" ; argname = "bin"

  set_up(1)

  argname <- deparse(substitute(bin))
  check_arg(bin, "list | vector | os formula")

  #
  # DSB expansion
  #

  bin_dp <- deparse_long(bin)
  if (grepl(".[", bin_dp, fixed = TRUE)) {
    bin_names <- names(bin)
    if (!is.null(bin_names)) {
      qui <- which(grepl(".[", bin_names, fixed = TRUE))
      for (i in qui) {
        bin_names[i] <- error_sender(
          dsb(bin_names[i],
            frame = parent.frame(2), nest = FALSE,
            vectorize = TRUE, collapse = ""
          ),
          dsb(
            "Error when binning: the name (.[bin_names[i]]) expanded",
            " with '.[]' led to an error:"
          )
        )
      }
    }
    names(bin) <- bin_names

    qui <- which(sapply(bin, function(x) is.character(x) && any(grepl(".[", x, fixed = TRUE))))
    for (i in qui) {
      value <- as.list(bin[[i]])
      qui_bis <- which(grepl(".[", value, fixed = TRUE))
      for (j in qui_bis) {
        value[[j]] <- error_sender(
          .dsb(value[[j]], frame = parent.frame(2), nest = FALSE),
          dsb(
            "Error when binning: the name (.[value[[j]]]) expanded",
            " with '.[]' led to an error:"
          )
        )
      }
      bin[[i]] <- unlist(value)
    }
  }

  #
  # cut:: => special treatment, we short circuit
  #

  if (!is.list(bin) && (is.character(bin) && any(grepl("^cut::", bin)))) {
    if (!grepl("^cut::", bin[1])) {
      stop_up("To use the special binning 'cut::values', it must be in the first position of the vector, possibly followed by custom names. Currently 'cut::values' is in position ", which(grepl("^cut::", bin))[1], ".")
    }

    n_names <- length(bin) - 1
    if (n_names > 1) {
      my_cut <- trimws(sub("^cut::", "", bin[1]))
      if (grepl("^[[:digit:]]+$", my_cut)) {
        n_cuts <- as.numeric(my_cut)
      } else {
        n_cuts <- length(strsplit(my_cut, "\\[|\\]")[[1]]) + 1
      }

      if (n_names > n_cuts) {
        stop_up("In the special binning 'cut::values', the number of custom names is greater than the number of bins (", n_names, " vs ", n_cuts, ").")
      }
    }

    if (!is.numeric(x)) {
      stop_up("To use the special binning 'cut::values', the variable '", varname, "' must be numeric. Currently this is not the case (it is of class ", enumerate_items(class(x)), " instead).")
    }

    return(cut_vector(x, bin))
  }

  x_int <- to_integer(x, add_items = TRUE, items.list = TRUE, sorted = TRUE)
  x_items <- x_int$items

  do_factor <- FALSE
  if (is.factor(x_items)) {
    do_factor <- TRUE
    x_items <- as.character(x_items)
  }

  if (!is.list(bin)) {
    if (is.character(bin) && any(grepl("^!?!?bin::", bin))) {
      if (length(bin) > 1) {
        stop_up("To use the special binning 'bin::digit', the argument '", argname, "' must be of length 1. Currently it is of length ", length(bin), ".")
      }

      d <- gsub("^!?!?bin::", "", bin)
      if (any(grepl("[^[:digit:]]", d))) {
        bin_type <- gsub("^(!?!?bin).*", "\\1", bin)
        stop_up("In the argument bin, the special binning must be of the form '", bin_type, "::digit'. Currently this is not the case for '", bin_type, "::", d, "'.")
      }
      d <- as.numeric(d)

      consecutive <- grepl("^!", bin)
      from_last <- grepl("^!!", bin)

      if (!consecutive) {
        if (!is.numeric(x_items)) {
          stop_up("To use the special binning 'bin::digit', the variable '", varname, "' must be numeric. Currently this is not the case (it is of class ", enumerate_items(class(x_items)), " instead).")
        }

        new_x <- (x_items %/% d) * d
        new_x_unik <- unique(new_x)
      } else {
        n <- length(x_items)
        x_seq <- if (from_last) n:1 else 1:n
        new_x <- ((x_seq - 1) %/% d) * d + 1
        new_x_unik <- unique(new_x)
      }

      bin <- list()
      for (i in seq_along(new_x_unik)) {
        bin[[i]] <- x_items[new_x == new_x_unik[i]]
      }
    } else {
      bin <- list(bin)
    }
  }

  # we catch "@" as first item
  if (identical(bin[[1]], "@") && length(bin) > 1) {
    # this means that the user wants to specify
    # the binned values as first elements of the new factor
    bin[[1]] <- NULL
    if (is.null(names(bin))) {
      names(bin) <- paste0("@", seq_along(bin), " ", names(bin))
    } else {
      qui_ok <- !grepl("^@\\d+", names(bin))
      names(bin)[qui_ok] <- paste0("@", seq(sum(qui_ok)), " ", names(bin)[qui_ok])
    }
  }

  x_map <- x_range <- seq_along(x_int$items)
  id_bin <- list()
  for (i in seq_along(bin)) {
    id_bin[[i]] <- items_to_drop(x_items, bin[[i]], varname,
      up = 2, argname = argname,
      keep_first = TRUE, no_error = no_error, valid_ref = TRUE
    )
  }

  # sanity check
  if (any(table(unlist(id_bin)) > 1)) {
    t_id <- table(unlist(id_bin))
    n_max <- max(t_id)
    pblm <- x_items[as.numeric(names(t_id)[t_id == n_max])]
    stop_up("In 'bin', some values are binned in different bins, it's of course not allowed. The value '", pblm, "' is in ", n_max, " bins.")
  }

  # recreating the factor

  # case no_error = TRUE
  if (any(lengths(id_bin) == 0)) {
    qui_ok <- lengths(id_bin) > 0
    id_bin <- id_bin[qui_ok]
    if (length(id_bin) == 0) {
      return(x)
    }
    bin <- bin[qui_ok]
  }

  id_not_core <- FALSE
  for (i in seq_along(bin)) {
    id_bin_i <- id_bin[[i]]
    id_core <- id_bin_i[1]

    id_change <- x_range %in% setdiff(id_bin_i, id_core)
    id_not_core <- id_not_core | id_change

    x_map <- x_map - cumsum(id_change)
  }

  # final pass
  for (i in seq_along(bin)) {
    id_bin_i <- id_bin[[i]]
    id_core <- id_bin_i[1]
    id_change <- x_range %in% id_bin_i

    x_map[id_change] <- x_map[id_core]
  }

  # changing the item values
  x_items_new <- x_items

  bin_names <- names(bin)
  custom_loc <- FALSE # custom location
  if (is.null(bin_names)) {
    bin_names <- character(length(bin))
  } else if (any(grepl("^@\\d", bin_names))) {
    custom_loc <- TRUE
    do_factor <- TRUE
    x_items_new <- as.character(x_items_new)

    qui <- grepl("^@\\d", bin_names)
    bin_location <- sub("^@(\\d+).*", "\\1", bin_names)
    bin_location[!qui] <- NA
    bin_location <- as.numeric(bin_location)

    bin_names <- sub("^@\\d+ *", "", bin_names)
  }

  for (i in seq_along(bin)) {
    b_name <- bin_names[[i]]
    if (nchar(b_name) == 0) {
      b_name <- as.character(x_items_new[id_bin[[i]][1]])
      bin_names[[i]] <- b_name # useful for new_loc
    }

    if (is.numeric(x_items_new)) {
      b_name_num <- tryCatch(as.numeric(b_name), warning = function(x) "not numeric")
      if (identical(b_name_num, "not numeric")) {
        # we convert to character
        do_factor <- TRUE
        x_items_new <- as.character(x_items_new)
        b_name_num <- b_name
      }

      x_items_new[id_bin[[i]]] <- b_name_num
    } else {
      x_items_new[id_bin[[i]]] <- b_name
    }
  }

  x_items_new <- x_items_new[!id_not_core]

  if (custom_loc) {
    n_bins <- length(bin_names)
    new_loc <- rep(NA_real_, length(x_items_new))
    for (i in 1:n_bins) {
      b_name <- bin_names[i]
      new_loc[x_items_new == b_name] <- bin_location[i]
    }

    loc_no_na <- new_loc[!is.na(new_loc)]
    while (anyDuplicated(loc_no_na)) {
      # regularization
      is_dup <- duplicated(loc_no_na)
      value <- min(loc_no_na[is_dup])
      i_first <- which.max(loc_no_na == value)
      loc_no_na[loc_no_na >= value] <- loc_no_na[loc_no_na >= value] + 1
      loc_no_na[i_first] <- value
    }

    n_max <- length(x_items_new)
    if (any(loc_no_na > n_max)) {
      # regularization

      n_lnona <- length(loc_no_na)
      if (n_lnona == 1) {
        loc_no_na <- n_max
      } else {
        order_loc_no_na <- order(loc_no_na)

        loc_sorted <- loc_no_na[order_loc_no_na]
        loc_sorted[n_lnona] <- n_max
        for (i in n_lnona:2) {
          if (loc_sorted[i] <= loc_sorted[i - 1]) {
            loc_sorted[i - 1] <- loc_sorted[i] - 1
          }
        }

        loc_no_na <- loc_sorted[order(order_loc_no_na)]
      }
    }

    new_loc[!is.na(new_loc)] <- loc_no_na

    for (i in seq_along(new_loc)) {
      if (!i %in% new_loc) {
        j <- which(is.na(new_loc))[1]
        new_loc[j] <- i
      }
    }

    x_map <- new_loc[x_map]
    x_items_new <- x_items_new[order(new_loc)]
  }

  if (do_factor) {
    x_items_new <- factor(x_items_new, levels = x_items_new)
  }

  x_new <- x_items_new[x_map[x_int$x]]

  return(x_new)
}


cut_vector <- function(x, bin) {
  # We cut a vector into pieces
  # **only numeric vectors**
  # cut::a]b[c] cuts the vectors into four slices: [-Inf, a]; ]a, b[; [b, c]; ]c, Inf]
  # a, b, c should be numeric values
  # they can be replaced with: pXX or qX, percentiles and quartiles
  # cut::n splits the data into n equal sizes

  set_up(2)

  bin_names_usr <- bin[-1]
  bin <- bin[1]

  # checking bin
  if (!is.character(bin) || !grepl("^cut::", bin)) {
    stop("Internal bug: Argument 'bin' should be equal to 'cut::stg' over here -- this is not the case.")
  }

  # cleaning x
  n <- length(x)
  ANY_NA <- anyNA(x)
  IS_NA <- FALSE
  if (ANY_NA) {
    IS_NA <- is.na(x)

    if (all(IS_NA)) {
      return(x)
    }

    x <- x[!IS_NA]
  }

  # we sort x (needed later)
  x_order <- order(x)
  x_sorted <- x[x_order]

  my_cut <- gsub("^cut::", "", bin)

  if (is_numeric_in_char(my_cut)) {
    my_cut <- as.numeric(my_cut)
    pct <- cumsum(rep(1 / my_cut, max(my_cut - 1, 1)))
    cut_points <- quantile(x_sorted, pct)

    bounds <- rep("]", length(cut_points))
  } else {
    bounds <- regmatches(my_cut, gregexpr("\\[|\\]", my_cut))[[1]]
    values <- strsplit(my_cut, "\\[|\\]")[[1]]

    n_cuts <- length(values)

    if (n_cuts != length(bounds)) {
      stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles/percentiles. Problem: each number must be followed by an open (or closed) square bracket, this is currentlly not the case.")
    }

    cut_points <- numeric(n_cuts)
    for (i in seq_along(values)) {
      v <- trimws(values[i])

      if (is_numeric_in_char(v)) {
        cut_points[i] <- as.numeric(v)
      } else if (grepl("^p|P|q|Q", v)) {
        p <- gsub("^.", "", v)
        if (!is_numeric_in_char(p)) {
          stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value '", v, "', in '", bin, "', is incorrect.")
        }

        p <- as.numeric(p)


        if (grepl("^q|Q", v)) {
          # we transform the quartile into a percentile
          if (!p %in% c(0:4)) {
            stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value '", v, "', in '", bin, "', is incorrect. \n  The quartile must be an integer between 0 and 4.")
          }

          p <- c(0, 25, 50, 75, 100)[p + 1]
        }

        if (!p %in% 0:100) {
          stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value '", v, "', in '", bin, "', is incorrect. \n  The percentile must be an integer between 0 and 100.")
        }

        cut_points[i] <- quantile(x, p / 100)
      } else {
        stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value '", v, "', in '", bin, "', is incorrect. This is not a percentile nor a number.")
      }
    }

    if (n_cuts > 1 && any(cut_points[1:(n_cuts - 1)] > cut_points[1 + 1:(n_cuts - 1)])) {
      i_pblm <- which(cut_points[1:(n_cuts - 1)] > cut_points[1 + 1:(n_cuts - 1)])[1]
      stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The values 'a', 'b', etc should be increasing, but the ", n_th(i_pblm), " value (", cut_points[i_pblm], ") is larger than the ", n_th(i_pblm + 1), " (", cut_points[i_pblm + 1], ").")
    }
  }

  bin_names <- character(length(cut_points) + 1)
  for (i in seq_along(bin_names_usr)) {
    bi <- bin_names_usr[i]
    if (!is.na(bi)) {
      bin_names[i] <- bi
    }
  }

  is_included <- 1 * (bounds == "]")
  x_cut <- cpp_cut(x_sorted, cut_points, is_included)

  x_int <- x_cut$x_int
  isnt_empty <- x_cut$isnt_empty
  value_min <- x_cut$value_min
  value_max <- x_cut$value_max
  is_int <- x_cut$is_int

  n_bins <- sum(isnt_empty)

  if (any(isnt_empty == 0)) {
    i_empty <- which(isnt_empty == 0)
    if (getFixest_notes()) {
      message("When binning: in '", bin, "', the ", enumerate_items(n_th(i_empty)), " bin", plural_len(i_empty, "s.is"), " empty.")
    }

    x_int <- cumsum(isnt_empty)[x_int]
    value_min <- value_min[-i_empty]
    value_max <- value_max[-i_empty]
    bin_names <- bin_names[-i_empty]
  }

  # creating the labels
  labels <- character(n_bins)
  if (is_int) {
    # format A-B

    # we don't want to add a comma to the years! But to large numbers: yes
    fmt_fun <- if (value_max[n_bins] > 9999) dreamerr::fsignif else function(x) as.character(x)

    for (i in 1:n_bins) {
      if (value_min[i] == value_max[i]) {
        labels[i] <- fmt_fun(value_min[i])
      } else {
        labels[i] <- paste0(fmt_fun(value_min[i]), "-", fmt_fun(value_max[i]))
      }
    }
  } else {
    # format [A, B]
    # we always write in inclusion, too hard to do the exclusion mentally
    d_min <- NULL
    for (i in 1:n_bins) {
      vmin <- value_min[i]
      vmax <- value_max[i]

      l10_diff <- if (vmax == vmin) 0 else log10(vmax - vmin)
      # d: nber of digits
      d <- if (l10_diff >= 0) 1 else ceiling(abs(l10_diff))
      if (!is.null(d_min) && d < d_min) d <- d_min

      # we check consistency problems
      d_min <- NULL
      vmax_fmt <- sprintf("%.*f", d, vmax)
      if (i != n_bins) {
        vmin_next_fmt <- sprintf("%.*f", d, value_min[i + 1])
        if (vmax_fmt == vmin_next_fmt) {
          l10_diff <- log10(value_min[i + 1] - vmax)
          d <- ceiling(abs(l10_diff))
          d_min <- d
        }
        vmax_fmt <- cpp_add_commas(vmax, d, FALSE)
      }
      vmin_fmt <- cpp_add_commas(vmin, d, FALSE)

      labels[i] <- paste0("[", vmin_fmt, "; ", vmax_fmt, "]")
    }
  }

  for (i in seq_along(bin_names)) {
    bi <- bin_names[i]
    if (nchar(bi) > 0) {
      labels[i] <- bi
    }
  }


  # Te result
  if (ANY_NA) {
    res <- rep(NA_real_, n)
    res[!IS_NA] <- x_int[order(x_order)]
  } else {
    res <- x_int[order(x_order)]
  }

  res <- factor(res, labels = labels)

  return(res)
}


fixest_pvalue <- function(x, zvalue, vcov) {
  # compute the pvalue for a fixest estimation

  if (use_t_distr(x)) {
    if (missing(vcov)) {
      stop("Internal error (=bug): the argument 'vcov' should not be missing in fixest_pvalue().")
    } else if (is.null(attr(vcov, "dof.K"))) {
      stop("Internal error (=bug): the attribute 'dof.K' from 'vcov' should not be NULL.")
    }

    # df.t is always an attribute of the vcov
    df.t <- attr(vcov, "df.t")
    if (is.null(df.t)) {
      df.t <- max(nobs(x) - attr(vcov, "dof.K"), 1)
    }

    pvalue <- 2 * pt(-abs(zvalue), df.t)
  } else {
    pvalue <- 2 * pnorm(-abs(zvalue))
  }

  pvalue
}

fixest_CI_factor <- function(x, level, vcov) {
  val <- (1 - level) / 2
  val <- c(val, 1 - val)

  if (use_t_distr(x)) {
    if (missing(vcov)) {
      stop("Internal error (=bug): the argument 'vcov' should not be missing in fixest_CI_factor().")
    } else if (is.null(attr(vcov, "dof.K"))) {
      stop("Internal error (=bug): the attribute 'dof.K' from 'vcov' should not be NULL.")
    }

    # df.t is always an attribute of the vcov
    df.t <- attr(vcov, "df.t")
    if (is.null(df.t)) {
      df.t <- max(nobs(x) - attr(vcov, "dof.K"), 1)
    }

    fact <- qt(val, df.t)
  } else {
    fact <- qnorm(val)
  }

  fact
}

#  x = ~a + (b + c) + d*k
extract_vars_simple <- function(x) {
  # We don't apply terms!
  # => m * (a + b) is not expanded
  # a bit slow

  check_arg(x, "formula")

  res <- c()
  x <- x[[length(x)]]

  while (is_operator(x, "+")) {
    res <- c(res, deparse_long(x[[3]]))
    x <- x[[2]]
  }

  res <- c(res, deparse_long(x))

  rev(res)
}

# x = "a + (b + c) + d*k"
char_to_vars <- function(x) {
  check_arg(x, "character scalar")

  if (!grepl("(", x, fixed = TRUE)) {
    res <- trimws(strsplit(x, "+", fixed = TRUE)[[1]])
    return(res)
  }

  res <- c()
  x <- str2lang(x)

  while (is_operator(x, "+")) {
    res <- c(res, deparse_long(x[[3]]))
    x <- x[[2]]
  }

  res <- c(res, deparse_long(x))

  rev(res)
}
