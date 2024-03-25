#include "04_0_linear_model.hpp"

bool sparse_check(const doubles_matrix<> &X) {
  // Super cheap sparsity test
  // works even for small K large N

  int N = X.nrow();
  int K = X.ncol();

  if (K < 5) return false;
  if (N < 1000 && K < 100) return false;
  if (N < 100)
    return false;  // avoids corner cases where you have more vars than obs

  // Let's look at some rows
  // If there are factors, there should be 0s
  //
  // Let's put it differently: if there are no 0s, then it's not factors for
  // sure

  int n0_first_row = 0;
  int n0_middle_row = 0;
  int n0_last_row = 0;

  int mid = N / 2;

  for (int k = 0; k < K; ++k) {
    n0_first_row += X(0, k) == 0;
    n0_middle_row += X(mid, k) == 0;
    n0_last_row += X(N - 1, k) == 0;
  }

  // Rcout << "1st: " << n0_first_row << ", mid: " << n0_middle_row << ", last:
  // " << n0_last_row << "\n";

  if (n0_first_row > K / 2 && n0_middle_row > K / 2 && n0_last_row > K / 2)
    return true;

  return false;
}

void set_sparse(vector<int> &n_j, vector<int> &start_j, vector<int> &all_i,
                vector<double> &x, const doubles_matrix<> &X,
                const doubles &w) {
  int N = X.nrow();
  int K = X.ncol();

  bool isWeight = w.size() > 1;

  int n_obs = 0;
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < N; ++i) {
      if (X(i, k) != 0) {
        ++n_j[k];
        all_i.push_back(i);

        if (isWeight) {
          // BEWARE: here there is no square root for w!
          // (as opposed to the non-sparse case)
          x.push_back(X(i, k) * w[i]);
        } else {
          x.push_back(X(i, k));
        }
      }
    }

    n_obs += n_j[k];
    start_j[k + 1] = n_obs;
  }
}