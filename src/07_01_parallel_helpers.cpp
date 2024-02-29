#include "07_0_parallel.hpp"

std::vector<int> set_parallel_scheme_bis(int N, int nthreads) {
  // => this concerns only the parallel application on a 1-Dimensional matrix
  // takes in the nber of observations of the vector and the nber of threads
  // gives back a vector of the length the nber of threads + 1 giving the
  // start/stop of each threads

  std::vector<int> res(nthreads + 1, 0);
  double N_rest = N;

  for (int i = 0; i < nthreads; ++i) {
    res[i + 1] = ceil(N_rest / (nthreads - i));
    N_rest -= res[i + 1];
    res[i + 1] += res[i];
  }

  return res;
}

[[cpp11::register]] int cpp_get_nb_threads_() { return omp_get_max_threads(); }

[[cpp11::register]] list cpp_which_na_inf_vec_(SEXP x, int nthreads) {
  /*
    This function takes a vector and looks at whether it contains NA or infinite
    values return: flag for na/inf + logical vector of obs that are na/inf x is
    ALWAYS a numeric vector std::isnan, std::isinf are OK since cpp11 required
    do_any_na_inf: if high suspicion of NA present: we go directly constructing
    the vector is_na_inf in the "best" case (default expected), we need not
    construct is_na_inf
  */

  int nobs = Rf_length(x);
  double *px = REAL(x);
  bool anyNAInf = false;
  bool any_na = false;   // return value
  bool any_inf = false;  // return value

  /*
    we make parallel the anyNAInf loop
    why? because we want that when there's no NA (default) it works as fast as
    possible if there are NAs, single threaded mode is faster, but then we
    circumvent with the do_any_na_inf flag
  */

  // no need to care about the race condition
  // "trick" to make a break in a multi-threaded section

  std::vector<int> bounds = set_parallel_scheme_bis(nobs, nthreads);
#pragma omp parallel for num_threads(nthreads)
  for (int t = 0; t < nthreads; ++t) {
    for (int i = bounds[t]; i < bounds[t + 1] && !anyNAInf; ++i) {
      if (std::isnan(px[i]) || std::isinf(px[i])) {
        anyNAInf = true;
      }
    }
  }

  // object to return: is_na_inf
  writable::logicals is_na_inf(nobs);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < nobs; ++i) {
    is_na_inf[i] = false;
  }

  if (anyNAInf) {
// again: no need to care about race conditions
#pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < nobs; ++i) {
      double x_tmp = px[i];
      if (std::isnan(x_tmp)) {
        is_na_inf[i] = true;
        any_na = true;
      } else if (std::isinf(x_tmp)) {
        is_na_inf[i] = true;
        any_inf = true;
      }
    }
  }

  // Return
  return writable::list({"any_na"_nm = any_na, "any_inf"_nm = any_inf,
                         "any_na_inf"_nm = any_na || any_inf,
                         "is_na_inf"_nm = is_na_inf});
}

[[cpp11::register]] list cpp_which_na_inf_mat_(doubles_matrix<> mat,
                                               int nthreads) {
  // almost identical to cpp_which_na_inf_vec but for R matrices. Changes:
  // - main argument becomes doubles_matrix<>
  // - k-for loop within the i-for loop
  /*
     This function takes a matrix and looks at whether it contains NA or
     infinite values return: flag for na/inf + logical vector of obs that are
     Na/inf std::isnan, std::isinf are OK since cpp11 required do_any_na_inf: if
     high suspicion of NA present: we go directly constructing the vector
     is_na_inf in the "best" case (default expected), we need not construct
     is_na_inf
  */

  int nobs = mat.nrow();
  int K = mat.ncol();
  bool anyNAInf = false;
  bool any_na = false;   // return value
  bool any_inf = false;  // return value

  /*
    we make parallel the anyNAInf loop
    why? because we want that when there's no NA (default) it works as fast as
    possible if there are NAs, single threaded mode is faster, but then we
    circumvent with the do_any_na_inf flag
  */

  // no need to care about the race condition
  // "trick" to make a break in a multi-threaded section

  std::vector<int> bounds = set_parallel_scheme_bis(nobs, nthreads);

#pragma omp parallel for num_threads(nthreads)
  for (int t = 0; t < nthreads; ++t) {
    for (int k = 0; k < K; ++k) {
      for (int i = bounds[t]; i < bounds[t + 1] && !anyNAInf; ++i) {
        if (std::isnan(mat(i, k)) || std::isinf(mat(i, k))) {
          anyNAInf = true;
        }
      }
    }
  }

  // object to return: is_na_inf
  writable::logicals is_na_inf(nobs);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < nobs; ++i) {
    is_na_inf[i] = false;
  }

  if (anyNAInf) {
#pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < nobs; ++i) {
      double x_tmp = 0;
      for (int k = 0; k < K; ++k) {
        x_tmp = mat(i, k);
        if (std::isnan(x_tmp)) {
          is_na_inf[i] = true;
          any_na = true;
          break;
        } else if (std::isinf(x_tmp)) {
          is_na_inf[i] = true;
          any_inf = true;
          break;
        }
      }
    }
  }

  // Return
  return writable::list({"any_na"_nm = any_na, "any_inf"_nm = any_inf,
                         "any_na_inf"_nm = any_na || any_inf,
                         "is_na_inf"_nm = is_na_inf});
}

[[cpp11::register]] list cpp_which_na_inf_df_(SEXP df, int nthreads) {
  // almost identical to cpp_which_na_inf_vec but for R **numeric** data frames.
  // Changes:
  // - main argument becomes SEXP
  // - k-for loop within the i-for loop
  /*
   This function takes a df and looks at whether it contains NA or infinite
   values return: flag for na/inf + logical vector of obs that are Na/inf
   std::isnan, std::isinf are OK since cpp11 required
   in the "best" case (default expected), we need not construct is_na_inf
   */

  int K = Rf_length(df);
  int nobs = Rf_length(VECTOR_ELT(df, 0));
  bool anyNAInf = false;
  bool any_na = false;   // return value
  bool any_inf = false;  // return value

  // The Mapping of the data
  std::vector<double *> df_data(K);
  for (int k = 0; k < K; ++k) {
    df_data[k] = REAL(VECTOR_ELT(df, k));
  }

  /*
   we make parallel the anyNAInf loop
   why? because we want that when there's no NA (default) it works as fast as
   possible if there are NAs, single threaded mode is faster, but then we
   circumvent with the do_any_na_inf flag
   */

  // no need to care about the race condition
  // "trick" to make a break in a multi-threaded section

  std::vector<int> bounds = set_parallel_scheme_bis(nobs, nthreads);

#pragma omp parallel for num_threads(nthreads)
  for (int t = 0; t < nthreads; ++t) {
    for (int k = 0; k < K; ++k) {
      for (int i = bounds[t]; i < bounds[t + 1] && !anyNAInf; ++i) {
        if (std::isnan(df_data[k][i]) || std::isinf(df_data[k][i])) {
          anyNAInf = true;
        }
      }
    }
  }

  // object to return: is_na_inf
  writable::logicals is_na_inf(nobs);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < nobs; ++i) {
    is_na_inf[i] = false;
  }

  if (anyNAInf) {
#pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < nobs; ++i) {
      double x_tmp = 0;
      for (int k = 0; k < K; ++k) {
        x_tmp = df_data[k][i];
        if (std::isnan(x_tmp)) {
          is_na_inf[i] = true;
          any_na = true;
          break;
        } else if (std::isinf(x_tmp)) {
          is_na_inf[i] = true;
          any_inf = true;
          break;
        }
      }
    }
  }

  // Return
  return writable::list({"any_na"_nm = any_na, "any_inf"_nm = any_inf,
                         "any_na_inf"_nm = any_na || any_inf,
                         "is_na_inf"_nm = is_na_inf});
}

[[cpp11::register]] integers cpp_check_only_0_(doubles_matrix<> x_mat,
                                               int nthreads) {
  // returns a 0/1 vectors => 1 means only 0

  int n = x_mat.nrow();
  int K = x_mat.ncol();

  writable::integers res(K);

#pragma omp parallel for num_threads(nthreads)
  for (int k = 0; k < K; ++k) {
    bool is_zero = true;
    for (int i = 0; i < n; ++i) {
      if (x_mat(i, k) != 0) {
        is_zero = false;
        break;
      }
    }

    res[k] = is_zero;
  }

  return res;
}
