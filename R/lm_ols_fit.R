ols_fit <- function(y, X, w, correct_0w = FALSE, collin.tol, nthreads, xwx = NULL,
                    xwy = NULL, only.coef = FALSE) {
  # No control here -- done before

  if (is.null(xwx)) {
    info_products <- cpp_sparse_products(X, w, y, correct_0w, nthreads)
    xwx <- info_products$XtX
    xwy <- info_products$Xty
  }

  multicol <- FALSE
  info_inv <- cpp_cholesky(xwx, collin.tol, nthreads)

  if (!is.null(info_inv$all_removed)) {
    # Means all variables are collinear! => can happen when using FEs
    return(list(all_removed = TRUE))
  }

  xwx_inv <- info_inv$XtX_inv
  is_excluded <- info_inv$id_excl

  multicol <- any(is_excluded)

  if (only.coef) {
    if (multicol) {
      beta <- rep(NA_real_, length(is_excluded))
      beta[!is_excluded] <- as.vector(xwx_inv %*% xwy[!is_excluded])
    } else {
      beta <- as.vector(xwx_inv %*% xwy)
    }

    return(beta)
  }

  if (multicol) {
    beta <- as.vector(xwx_inv %*% xwy[!is_excluded])
    fitted.values <- cpppar_xbeta(X[, !is_excluded, drop = FALSE], beta, nthreads)
  } else {
    # avoids copies
    beta <- as.vector(xwx_inv %*% xwy)
    fitted.values <- cpppar_xbeta(X, beta, nthreads)
  }

  residuals <- y - fitted.values

  res <- list(
    xwx = xwx, coefficients = beta, fitted.values = fitted.values,
    xwx_inv = xwx_inv, multicol = multicol, residuals = residuals,
    is_excluded = is_excluded, collin.min_norm = info_inv$min_norm
  )

  res
}
