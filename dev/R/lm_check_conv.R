# TODO: port this to C++
check_conv <- function(y, X, fixef_id_list, slope_flag, slope_vars, weights, full = FALSE, fixef_names = NULL) {
  # y, X => variables that were demeaned

  # For each variable: we compute the optimal FE coefficient
  # it should be 0 if the algorithm converged

  Q <- length(slope_flag)

  use_y <- TRUE
  if (length(X) == 1) {
    K <- 1
    nobs <- length(y)
  } else {
    nobs <- NROW(X)

    if (length(y) <= 1) {
      use_y <- FALSE
    }

    K <- NCOL(X) + use_y
  }

  info_max <- list()
  info_full <- vector("list", K)

  for (k in 1:K) {
    if (k == 1 && use_y) {
      x <- y
    } else if (use_y) {
      x <- X[, k - use_y]
    } else {
      x <- X[, k]
    }

    info_max_tmp <- c()
    info_full_tmp <- list()
    index_slope <- 1
    for (q in 1:Q) {
      fixef_id <- fixef_id_list[[q]]

      if (slope_flag[q] >= 0) {
        x_means <- tapply(weights * x, fixef_id, mean)
        info_max_tmp <- c(info_max_tmp, max(abs(x_means)))

        if (full) {
          info_full_tmp[[length(info_full_tmp) + 1]] <- x_means
        }
      }

      n_slopes <- abs(slope_flag[q])
      if (n_slopes > 0) {
        for (i in 1:n_slopes) {
          var <- slope_vars[[index_slope]]

          num <- tapply(weights * x * var, fixef_id, sum)
          denom <- tapply(weights * var^2, fixef_id, sum)

          x_means <- num / denom
          info_max_tmp <- c(info_max_tmp, max(abs(x_means)))

          if (full) {
            info_full_tmp[[length(info_full_tmp) + 1]] <- x_means
          }

          index_slope <- index_slope + 1
        }
      }
    }

    if (!is.null(fixef_names)) {
      names(info_full_tmp) <- fixef_names
    }

    info_max[[k]] <- info_max_tmp

    if (full) {
      info_full[[k]] <- info_full_tmp
    }
  }

  info_max <- do.call("rbind", info_max)

  if (full) {
    res <- info_full
  } else {
    res <- info_max
  }

  res
}
