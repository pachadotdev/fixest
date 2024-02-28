#include "05_0_misc.hpp"

// similar to table but faster
[[cpp11::register]] doubles cpp_table_(int Q, integers dum) {
  // Q: nber of classes
  // dum: the N vector of clusters

  int N = dum.size();

  writable::doubles res(Q);
  int i, q;

  for (i = 0; i < N; i++) {
    q = dum[i] - 1;  // we take 1 off => different indexation in C
    res[q]++;
  }

  return res;
}

[[cpp11::register]] double cpp_ssr_null_(doubles y, doubles w) {
  // simple fun to compute the ssr of the null ols model
  // 2/3 times faster than pure r

  bool is_weight = w.size() > 1;

  int n = y.size();
  double denom = 0;

  // the "mean"
  double y_mean = 0;
  for (int i = 0; i < n; ++i) {
    if (is_weight) {
      y_mean += y[i] * w[i];
      denom += w[i];
    } else {
      y_mean += y[i];
    }
  }

  if (is_weight) {
    y_mean = y_mean / denom;
  } else {
    y_mean = y_mean / n;
  }

  double res = 0, value = 0;
  for (int i = 0; i < n; ++i) {
    value = y[i] - y_mean;
    if (is_weight) {
      res += value * value * w[i];
    } else {
      res += value * value;
    }
  }

  return res;
}

[[cpp11::register]] double cpp_ssq_(doubles x, doubles w) {
  // simple fun to compute the sum of the square of the elt of a vector
  // 30% faster than pure r (twice faster with weights)

  bool is_weight = w.size() > 1;

  int n = x.size();

  // the mean
  double res = 0;
  for (int i = 0; i < n; i++) {
    if (is_weight) {
      res += x[i] * x[i] * w[i];
    } else {
      res += x[i] * x[i];
    }
  }

  return res;
}

[[cpp11::register]] int cpp_constant_dum_(int k, doubles x, integers dum,
                                          bool only_0 = false) {
  // number of values of dum for which x is constant

  int n_obs = dum.size();

  double ref = x[0];
  int dum_current = dum[0];

  bool found_different = only_0 ? ref != 0 : false;
  int nb_constant = 0;

  for (int i = 1; i < n_obs; ++i) {
    if (dum[i] != dum_current) {
      // new guy
      dum_current = dum[i];
      if (found_different == false) {
        ++nb_constant;
      }

      ref = x[i];
      found_different = only_0 ? ref != 0 : false;
    } else if (!found_different) {
      if (x[i] != ref) {
        found_different = true;
      }
    }
  }

  if (found_different == false) {
    ++nb_constant;
  }

  return nb_constant;
}
