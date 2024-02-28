#include "07_0_parallel.hpp"

[[cpp11::register]] doubles cpppar_exp_(doubles x, int nthreads = 1) {
  // parallel exponentiation using omp

  // int n = x.length();
  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    res[i] = exp(x[i]);
  }

  return res;
}

[[cpp11::register]] doubles cpppar_log_(doubles x, int nthreads = 1) {
  // parallel exponentiation using omp

  int n = x.size();
  writable::doubles res(n);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < n; ++i) {
    res[i] = log(x[i]);
  }

  return res;
}

[[cpp11::register]] doubles cpppar_log_a_exp_(int nthreads, double a,
                                              doubles mu, doubles exp_mu) {
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

  return res;
}

[[cpp11::register]] list cpppar_cond_means_(doubles_matrix<> mat_vars,
                                            integers treat, int nthreads = 1) {
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

      if (isnan(x_tmp) || isinf(x_tmp)) {
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
