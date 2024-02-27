/*******************************************************************
 * ________________________                                        *
 * || Parallel functions ||                                        *
 * ------------------------                                        *
 *                                                                 *
 * Author: Laurent R. Berge                                        *
 *                                                                 *
 * Group of functions doing simple things... but in parallel.      *
 *                                                                 *
 * The functions don't do much more than what their names suggest. *
 *                                                                 *
 ******************************************************************/

#include <math.h>

#include <cfloat>
#include <cpp11.hpp>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#include <pthread.h>
#else
#define omp_get_max_threads() 0
#endif

#include <Rmath.h>
#include <stdio.h>

#include <cmath>

using namespace cpp11;

// This file contains misc fixest functions parallelized with the omp library

// Regarding fork detection => I don't know enough yet
// The code below doesn't seem to work. And I can't check since I'm on
// Windows.... Safer not to include it for now.

// static bool fixest_in_fork = false;
//
// // Trick taken from data.table to detect forking
// // Actually I don't know if it works, and I can't check...
// void when_fork() {
//   fixest_in_fork = true;
// }
//
// void after_fork() {
//   fixest_in_fork = false;
// }
//
// [[cpp11::register]]
// void cpp_setup_fork_presence() {
//  // Called only once at startup
//  #ifdef _OPENMP
//     pthread_atfork(&when_fork, &after_fork, NULL);
//  #endif
// }
//
// [[cpp11::register]]
// bool cpp_is_in_fork(){
//     return fixest_in_fork;
// }

[[cpp11::register]] int cpp_get_nb_threads_() { return omp_get_max_threads(); }

// The following function is already defined in lm_related (I know...)
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

[[cpp11::register]] doubles cpp_exp_(doubles x, int nthreads) {
  // parallel exponentiation using omp

  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    res[i] = exp(x[i]);
  }

  return (res);
}

[[cpp11::register]] doubles cpp_log_(doubles x, int nthreads) {
  // parallel exponentiation using omp

  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    res[i] = log(x[i]);
  }

  return (res);
}

[[cpp11::register]] doubles cpp_log_a_exp_(double a, doubles mu, doubles exp_mu,
                                           int nthreads) {
  // faster this way

  int n = mu.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    if (mu[i] < 200) {
      res[i] = log(a + exp_mu[i]);
    } else {
      res[i] = mu[i];
    }
  }

  return (res);
}

[[cpp11::register]] doubles cpp_lgamma_(doubles x, int nthreads) {
  // parallel lgamma using omp

  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    res[i] = lgamma(x[i]);
  }

  return (res);
}

[[cpp11::register]] doubles cpp_digamma_(doubles x, int nthreads) {
  // parallel digamma using omp

  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    res[i] = Rf_digamma(x[i]);
  }

  return (res);
}

[[cpp11::register]] doubles cpp_trigamma_(doubles x, int nthreads) {
  // parallel trigamma using omp

  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    res[i] = Rf_trigamma(x[i]);
  }

  return (res);
}

inline double poisson_linkinv(double x) {
  return x < -36 ? DBL_EPSILON : exp(x);
}

[[cpp11::register]] doubles cpp_poisson_linkinv_(doubles x, int nthreads) {
  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    res[i] = poisson_linkinv(x[i]);
  }

  return (res);
}

[[cpp11::register]] bool cpp_poisson_validmu_(SEXP x, int nthreads) {
  int n = Rf_length(x);
  double *px = REAL(x);
  bool res = true;

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    double x_tmp = px[i];
    if (std::isinf(x_tmp) || x_tmp <= 0) {
      res = false;
    }
  }

  return res;
}

[[cpp11::register]] doubles cpp_logit_linkfun_(doubles x, int nthreads) {
  // parallel trigamma using omp

  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    double x_tmp = x[i];
    res[i] = log(x_tmp) - log(1 - x_tmp);
  }

  return (res);
}

inline double logit_linkinv(double x) {
  return x < -30    ? DBL_EPSILON
         : (x > 30) ? 1 - DBL_EPSILON
                    : 1 / (1 + 1 / exp(x));
}

[[cpp11::register]] doubles cpp_logit_linkinv_(doubles x, int nthreads) {
  // parallel trigamma using omp

  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    // res[i] = 1 / (1 + 1 / exp(x[i]));
    res[i] = logit_linkinv(x[i]);
  }

  return (res);
}

inline double logit_mueta(double x) {
  if (fabs(x) > 30) {
    return DBL_EPSILON;
  } else {
    double exp_x = exp(x);
    return (1 / ((1 + 1 / exp_x) * (1 + exp_x)));
  }
}

[[cpp11::register]] doubles cpp_logit_mueta_(doubles x, int nthreads) {
  // parallel trigamma using omp

  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    // double exp_x = exp(x[i]);
    // res[i] = 1 / ((1 + 1 / exp_x) * (1 + exp_x));
    res[i] = logit_mueta(x[i]);
  }

  return (res);
}

[[cpp11::register]] doubles cpp_logit_devresids_(doubles y, doubles mu,
                                                 doubles wt, int nthreads) {
  int n = mu.size();
  writable::doubles res(n);
  bool isWeight = wt.size() != 1;

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    if (y[i] == 1) {
      res[i] = -2 * log(mu[i]);
    } else if (y[i] == 0) {
      res[i] = -2 * log(1 - mu[i]);
    } else {
      res[i] = 2 * (y[i] * log(y[i] / mu[i]) +
                    (1 - y[i]) * log((1 - y[i]) / (1 - mu[i])));
    }

    if (isWeight)
      res[i] *= wt[i];
  }

  return (res);
}

[[cpp11::register]] doubles cpp_xwy_(doubles_matrix<> X, doubles y, doubles w,
                                     int nthreads) {
  int N = X.nrow();
  int K = X.ncol();

  bool isWeight = false;
  if (w.size() > 1) {
    isWeight = true;
  }

  writable::doubles res(K);

  // computation
#pragma omp parallel for num_threads(nthreads)
  for (int k = 0; k < K; ++k) {
    double val = 0;
    if (isWeight) {
      for (int i = 0; i < N; ++i) {
        val += X(i, k) * w[i] * y[i];
      }
    } else {
      for (int i = 0; i < N; ++i) {
        val += X(i, k) * y[i];
      }
    }
    res[k] = val;
  }

  return (res);
}

[[cpp11::register]] doubles cpp_xbeta_(doubles_matrix<> X, doubles beta,
                                       int nthreads) {
  int N = X.nrow();
  int K = X.ncol();

  writable::doubles res(N);

  // computation
#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < N; ++i) {
    double val = 0;
    for (int k = 0; k < K; ++k) {
      val += X(i, k) * beta[k];
    }
    res[i] = val;
  }

  return (res);
}

[[cpp11::register]] doubles_matrix<>
cpp_matprod_(doubles_matrix<> x, doubles_matrix<> y, int nthreads) {
  // => simply x %*% y

  int N = x.nrow();
  int K = x.ncol();

  writable::doubles_matrix<> xy(N, K);

  // computing xy
#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < K; ++k) {
      double value = 0;
      for (int l = 0; l < K; ++l) {
        value += x(i, l) * y(l, k);
      }
      xy(i, k) = value;
    }
  }

  return (xy);
}

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
  bool any_na = false;  // return value
  bool any_inf = false; // return value

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
  writable::logicals is_na_inf(anyNAInf ? nobs : 1);

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
  bool any_na = false;  // return value
  bool any_inf = false; // return value

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
  writable::logicals is_na_inf(anyNAInf ? nobs : 1);

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
  bool any_na = false;  // return value
  bool any_inf = false; // return value

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
  writable::logicals is_na_inf(anyNAInf ? nobs : 1);

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

[[cpp11::register]] list cpp_cond_means_(doubles_matrix<> mat_vars,
                                         integers treat, int nthreads) {
  // conditional means: function did_means

  int N = mat_vars.nrow();
  int K = mat_vars.ncol();

  // objects to return:
  writable::integers na_vect(K);
  writable::doubles_matrix<> mean_mat(K, 2);
  writable::doubles_matrix<> sd_mat(K, 2);
  writable::integers_matrix<> n_mat(K, 2);
  writable::integers n_01(2);

  // computation
#pragma omp parallel for num_threads(nthreads)
  for (int k = 0; k < K; ++k) {
    double sum_0 = 0, sum_1 = 0;
    double sum2_0 = 0, sum2_1 = 0;
    int n_0 = 0, n_1 = 0, n_na = 0;
    double x_tmp = 0;

    for (int i = 0; i < N; ++i) {
      x_tmp = mat_vars(i, k);

      if (std::isnan(x_tmp) || std::isinf(x_tmp)) {
        ++n_na;
      } else {
        if (treat[i] == 0) {
          sum_0 += x_tmp;
          sum2_0 += x_tmp * x_tmp;
          ++n_0;
        } else {
          sum_1 += x_tmp;
          sum2_1 += x_tmp * x_tmp;
          ++n_1;
        }
      }
    }

    // saving
    double m_0 = sum_0 / n_0;
    double m_1 = sum_1 / n_1;
    mean_mat(k, 0) = m_0;
    mean_mat(k, 1) = m_1;

    sd_mat(k, 0) = sqrt(sum2_0 / (n_0 - 1) - m_0 * sum_0 / (n_0 - 1));
    sd_mat(k, 1) = sqrt(sum2_1 / (n_1 - 1) - m_1 * sum_1 / (n_1 - 1));

    n_mat(k, 0) = n_0;
    n_mat(k, 1) = n_1;

    na_vect[k] = n_na;
  }

  // number of obs per treat case
  for (int i = 0; i < N; ++i) {
    if (treat[i] == 0) {
      ++n_01[0];
    } else {
      ++n_01[1];
    }
  }

  return writable::list({"means"_nm = mean_mat, "sd"_nm = sd_mat,
                         "n"_nm = n_mat, "n_01"_nm = n_01, "na"_nm = na_vect});
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
