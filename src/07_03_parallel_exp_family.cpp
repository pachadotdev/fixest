#include "07_0_parallel.hpp"

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

  return res;
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

    if (isWeight) res[i] *= wt[i];
  }

  return (res);
}
