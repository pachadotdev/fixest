#include "04_0_linear_model.hpp"

void invert_tri(writable::doubles_matrix<> &R, int K, int nthreads = 1) {
  // Strategy: we invert by bands (b) => better for parallelization

  // initialization of R prime
  for (int i = 0; i < K; ++i) {
    for (int j = i + 1; j < K; ++j) {
      // std::swap(R(i, j), R(j, i));
      double temp = R(i, j);
      R(j, i) = temp;
    }
  }

  // b0
  for (int i = 0; i < K; ++i) {
    R(i, i) = 1.0 / R(i, i);
  }

  // Check for interrupts
  // number of computations is (K - b) * (b + 1) => max is (K + 1)**2 / 2
  const double flop = (K + 1) * (K + 1) / 2.0;
  const int iterSecond =
      ceil(2000000000 / flop / 5); // nber iter per 1/5 second

  for (int b = 1; b < K; ++b) {
    if (b % iterSecond == 0) {
      check_user_interrupt();
    }

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
    for (int i = 0; i < K - b; ++i) {
      double numerator = 0;
      const int col = i + b;
      for (int k = i + 1; k <= col; ++k) {
        numerator -= R(k, i) * R(k, col);
      }
      R(i, col) = numerator * R(i, i);
    }
  }
}

void tproduct_tri(writable::doubles_matrix<> &RRt,
                  writable::doubles_matrix<> &R, int nthreads = 1) {
  const int K = RRt.ncol();

  // initialization of R prime
  for (int i = 0; i < K; ++i) {
    for (int j = i + 1; j < K; ++j) {
      double temp = R(i, j);
      R(j, i) = temp;
    }
  }

  // Check for interrupts
  // we do the same as for the invert_tri
  double flop = (K + 1) * (K + 1) / 2.0;
  int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second
  int n_iter_main = 0;

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
  for (int i = 0; i < K; ++i) {
    if (omp_get_thread_num() == 0 && n_iter_main % iterSecond == 0) {
      check_user_interrupt();
      ++n_iter_main;
    }

    for (int j = i; j < K; ++j) {
      double value = 0;
      int k_start = i < j ? j : i;
      for (int k = k_start; k < K; ++k) {
        value += R(k, j) * R(k, i);
      }
      RRt(i, j) = value;
      RRt(j, i) = value;
    }
  }
}

[[cpp11::register]] list cpp_cholesky_(doubles_matrix<> X, double tol,
                                       int nthreads) {
  // X is symmetric and semi-definite positive
  // rank-revealing on-the-fly

  writable::list res;

  int K = X.ncol();

  writable::doubles_matrix<> R(K, K);
  for (int i = 0; i < K; ++i) {
    for (int j = 0; j < K; ++j) {
      R(i, j) = 0;
    }
  }

  writable::logicals id_excl(K);
  for (int k = 0; k < K; ++k) {
    id_excl[k] = false;
  }

  int n_excl = 0;

  // we check for interrupt every 1s when it's the most computationally
  // intensive at each iteration we have K * (j+1) - j**2 - 2*j - 1
  // multiplications max => K**2/4
  double flop = K * K / 4.0;
  int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second
  double min_norm = X(0, 0);

  for (int j = 0; j < K; ++j) {
    if (j % iterSecond == 0) {
      check_user_interrupt();
    }

    // implicit pivoting (it's like 0-rank variables are stacked in the end)

    double R_jj = X(j, j);
    for (int k = 0; k < j; ++k) {
      if (id_excl[k] == true) {
        continue;
      }
      R_jj -= std::pow(R(k, j), 2);
    }

    if (R_jj < tol) {
      n_excl++;
      id_excl[j] = true;

      // Corner case, may happen:
      if (n_excl == K) {
        writable::list res;
        res.push_back({"all_removed"_nm = true});
        return res;
      }

      continue;
    }

    if (min_norm > R_jj) {
      min_norm = R_jj;
    }

    R_jj = sqrt(R_jj);

    R(j, j) = R_jj;

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
    for (int i = j + 1; i < K; ++i) {
      double value = X(i, j);
      for (int k = 0; k < j; ++k) {
        if (id_excl[k] == true) {
          continue;
        }
        value -= R(k, i) * R(k, j);
      }
      R(j, i) = value / R_jj;
    }
  }

  // we reconstruct the R matrix, otherwise it's just a mess to handle excluded
  // variables on the fly changes are in place

  if (n_excl > 0) {
    int n_j_excl = 0;

    // The first chunk of the matrix is identical
    // we find when it starts being different
    int j_start = 0;
    for (; id_excl[j_start] == false; ++j_start) {
      for (int j = j_start; j < K; ++j) {
        if (id_excl[j] == true) {
          ++n_j_excl;
          continue;
        }

        int n_i_excl = 0;
        for (int i = 0; i <= j; ++i) {
          if (id_excl[i] == true) {
            ++n_i_excl;
            continue;
          }
          R(i - n_i_excl, j - n_j_excl) += R(i, j);
        }
      }
    }

    K -= n_excl;
  }

  // Inversion of R, in place
  invert_tri(R, K, nthreads);

  writable::doubles_matrix<> XtX_inv(K, K);
  tproduct_tri(XtX_inv, R, nthreads);

  return writable::list({"XtX_inv"_nm = XtX_inv, "id_excl"_nm = id_excl,
                         "min_norm"_nm = min_norm});
}
