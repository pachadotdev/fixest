#include "07_0_parallel.hpp"

[[cpp11::register]] doubles cpppar_xwy_(doubles_matrix<> X, doubles y,
                                        doubles w, int nthreads = 1) {
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

  return res;
}

[[cpp11::register]] doubles cpppar_xbeta_(doubles_matrix<> X, doubles beta,
                                          int nthreads = 1) {
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

  return res;
}

[[cpp11::register]] doubles_matrix<> cpppar_matprod_(doubles_matrix<> x,
                                                     doubles_matrix<> y,
                                                     int nthreads = 1) {
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
