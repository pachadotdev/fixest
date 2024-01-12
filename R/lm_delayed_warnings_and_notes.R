# Delayed warnings and notes ----

setup_multi_notes <- function() {
  # We reset all the notes
  options("fixest_multi_notes" = NULL)
}

stack_multi_notes <- function(note) {
  all_notes <- getOption("fixest_multi_notes")
  options("fixest_multi_notes" = c(all_notes, note))
}

release_multi_notes <- function() {
  notes <- getOption("fixest_multi_notes")

  if (length(notes) > 0) {
    dict <- c(
      "\\$collin\\.var" = "Some variables have been removed because of collinearity (see $collin.var).",
      "collinear with the fixed-effects" = "All the variables are collinear with the fixed-effects (the estimation is void).",
      "virtually constant and equal to 0" = "All the variables are virtually constant and equal to 0 (the estimation is void).",
      "Very high value of theta" = "Very high value of theta, there is no sign of overdispersion, you may consider a Poisson model."
    )

    for (i in seq_along(dict)) {
      d_value <- dict[i]
      d_name <- names(dict)[i]

      x_i <- grepl(d_name, notes)
      x <- notes[x_i]
      # we prefer specific messages
      if (length(unique(x)) == 1) {
        next
      } else {
        notes[x_i] <- d_value
      }
    }

    tn <- table(notes)
    new_notes <- paste0("[x ", tn, "] ", names(tn))

    message("Notes from the estimations:\n", dsb("'\n'c?new_notes"))
  }

  options("fixest_multi_notes" = NULL)
}


warn_fixef_iter <- function(env, stack_multi = FALSE) {
  # Show warnings related to the nber of times the maximum of iterations was reached

  # For fixed-effect
  fixef.iter <- get("fixef.iter", env)
  fixef.iter.limit_reached <- get("fixef.iter.limit_reached", env)
  origin <- get("origin", env)
  warn <- get("warn", env)

  if (!warn) {
    return(invisible(NULL))
  }

  goWarning <- FALSE
  warning_msg <- ""
  if (fixef.iter.limit_reached > 0) {
    goWarning <- TRUE
    warning_msg <- paste0(origin, ": [Getting the fixed-effects] iteration limit reached (", fixef.iter, ").", ifelse(fixef.iter.limit_reached > 1, paste0(" (", fixef.iter.limit_reached, " times.)"), " (Once.)"))
  }

  # For the fixed-effect derivatives
  deriv.iter <- get("deriv.iter", env)
  deriv.iter.limit_reached <- get("deriv.iter.limit_reached", env)

  if (deriv.iter.limit_reached > 0) {
    prefix <- ifelse(goWarning, paste0("\n", sprintf("% *s", nchar(origin) + 2, " ")), paste0(origin, ": "))
    warning_msg <- paste0(warning_msg, prefix, "[Getting fixed-effects derivatives] iteration limit reached (", deriv.iter, ").", ifelse(deriv.iter.limit_reached > 1, paste0(" (", deriv.iter.limit_reached, " times.)"), " (Once.)"))
    goWarning <- TRUE
  }

  if (goWarning) {
    if (stack_multi) {
      stack_multi_notes(warning_msg)
    } else {
      warning(warning_msg, call. = FALSE, immediate. = TRUE)
    }
  }
}


warn_step_halving <- function(env, stack_multi = FALSE) {
  nb_sh <- get("nb_sh", env)
  warn <- get("warn", env)

  if (!warn) {
    return(invisible(NULL))
  }

  if (nb_sh > 0) {
    msg <- paste0("feglm: Step halving due to non-finite deviance (", ifelse(nb_sh > 1, paste0(nb_sh, " times"), "once"), ").")
    if (stack_multi) {
      stack_multi_notes(msg)
    } else {
      warning(msg, call. = FALSE, immediate. = TRUE)
    }
  }
}


format_error_msg <- function(x, origin) {
  # Simple formatting of the error msg

  # LATER:
  # - for object not found: provide a better error msg by calling the name of the missing
  #   argument => likely I'll need a match.call argument

  x <- gsub("\n+$", "", x)

  if (grepl("^Error (in|:|: in) (fe|fixest|fun|fml_split)[^\n]+\n", x)) {
    res <- gsub("^Error (in|:|: in) (fe|fixest|fun|fml_split)[^\n]+\n *(.+)", "\\3", x)
  } else if (grepl("[Oo]bject '.+' not found", x) || grepl("memory|cannot allocate", x)) {
    res <- x
  } else {
    res <- paste0(x, "\nThis error was unforeseen by the author of the function ", origin, ". If you think your call to the function is legitimate, could you report?")
  }
  res
}
