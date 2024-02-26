#include "05_0_misc.hpp"

simple_vec_double::simple_vec_double(SEXP x) {
  n = Rf_length(x);

  if (TYPEOF(x) == REALSXP) {
    px_double = REAL(x);
    is_real = true;
  } else if (TYPEOF(x) == INTSXP) {
    px_int = INTEGER(x);
    is_real = false;
  } else {
    stop("Error: Wrong argument type in cpp_factor_matrix.");
  }
}

double simple_vec_double::operator[](int i) {
  if (i >= n) {
    return 1;
  } else if (is_real) {
    return px_double[i];
  } else {
    return static_cast<double>(px_int[i]);
  }
}

[[cpp11::register]] doubles_matrix<> cpp_factor_matrix_(integers fact,
                                                        logicals is_na_all,
                                                        integers who_is_dropped,
                                                        SEXP var,
                                                        strings col_names) {
  // fact: integer vector from 1 (!) to K, can contain NAs
  // Checking Na is cheap as opposed to populating the matrix, but having an
  // argument avoids creating a new object

  int n = fact.size();
  int K = 0;

  // Finding out the TOTAL nber of cols (before removal)
  for (int i = 0; i < n; ++i) {
    if (!is_na_all[i] && K < fact[i]) {
      K = fact[i];
    }
  }

  // Are there values to remove? If so, we create a mapping
  int n_drop = who_is_dropped[0] == -1 ? 0 : Rf_length(who_is_dropped);
  bool IS_REMOVAL = n_drop > 0;
  std::vector<int> mapping;

  if (IS_REMOVAL) {
    mapping.resize(K);

    for (int k = 0; k < K; ++k) {
      mapping[k] = k;
    }

    // In mapping: values equal to -1 mean dropped value

    int k_drop = 0;
    for (int k = 0; k < K; ++k) {
      if (k_drop < n_drop && k + 1 == who_is_dropped[k_drop]) {
        ++k_drop;
        mapping[k] = -1;
      } else {
        mapping[k] -= k_drop;
      }
    }

    K -= k_drop;
  }

  writable::doubles_matrix<> res(n, K);

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < K; ++k) {
      res(i, k) = 0;
    }
  }

  // The interacted var
  simple_vec_double my_var(var);

  // Filling the matrix
  for (int i = 0; i < n; ++i) {
    if (is_na_all[i]) {
      // we fill the row
      for (int k = 0; k < K; ++k) {
        res(i, k) += NA_REAL;
      }
    } else if (IS_REMOVAL) {
      if (mapping[fact[i] - 1] != -1) {
        res(i, mapping[fact[i] - 1]) = my_var[i];
      }
    } else {
      res(i, fact[i] - 1) = my_var[i];
    }
  }

  res.attr("dimnames") = writable::list({R_NilValue, col_names});

  return res;
}

[[cpp11::register]] doubles cpp_diag_XUtX_(doubles_matrix<> X,
                                           doubles_matrix<> U) {
  // computes the diagonal of X %*% U %*% t(X)

  int n = X.nrow();
  int K = X.ncol();

  writable::doubles res(n);

  for (int i = 0; i < n; ++i) {
    res[i] = 0;
  }

  for (int i = 0; i < n; ++i) {
    double res_i = 0;
    for (int k = 0; k < K; ++k) {
      double xk = 0;
      for (int k2 = 0; k2 < K; ++k2) {
        xk += X(i, k2) * U(k, k2);
      }

      res_i += xk * X(i, k);
    }

    res[i] = res_i;
  }

  return res;
}
