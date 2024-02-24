#include "05_0_misc.hpp"

[[cpp11::register]] doubles
cpp_partialDerivative_other_(int iterMax, int Q, int N, double epsDeriv,
                             doubles ll_d2, doubles dx_dother, doubles init,
                             integers_matrix<> dumMat, integers nbCluster) {
  // takes in:
  // dumMat: the matrix of dummies (n X c) each obs => cluster // must be in cpp
  // index!!! init: the initialisation of the sum of derivatives vector ll_d2:
  // the second derivative dx_dother: the vector of dx_dother

  int i, q, c;
  int index;
  int sum_cases = 0;
  double new_value;

  // OLD: int start[Q], end[Q];
  // NEW: using ISO C++
  std::vector<int> start, end;
  start.reserve(Q);
  end.reserve(Q);

  for (q = 0; q < Q; q++) {
    // the total number of clusters (eg if man/woman and 10 countries: total of
    // 12 cases)
    sum_cases += nbCluster[q];
    if (q == 0) {
      start[q] = 0;
      end[q] = nbCluster[q];
    } else {
      start[q] = start[q - 1] + nbCluster[q - 1];
      end[q] = end[q - 1] + nbCluster[q];
    }
  }

  // the derivatives per cluster, stacked in one vector
  writable::doubles clusterDeriv(sum_cases);
  for (int i = 0; i < sum_cases; i++) {
    clusterDeriv[i] = 0;
  }

  writable::doubles sum_lld2(sum_cases);
  for (int i = 0; i < sum_cases; i++) {
    sum_lld2[i] = 0;
  }

  // Creation of the sum_lld2
  for (i = 0; i < N; i++) {
    for (q = 0; q < Q; q++) {
      index = start[q] + dumMat(i, q);
      sum_lld2[index] += ll_d2[i];
    }
  }

  // the result to return
  writable::doubles S(N);
  for (i = 0; i < N; i++) {
    S[i] = init[i];
  }

  bool ok = true;
  int iter = 0;
  while (ok & (iter < iterMax)) {
    iter++;
    ok = false;

    for (q = 0; q < Q; q++) {
      check_user_interrupt();

      for (i = 0; i < N; i++) {
        index = start[q] + dumMat(i, q);
        clusterDeriv[index] += dx_dother[i] + S[i] * ll_d2[i];
      }

      // finish calculating clusterDeriv + control
      for (c = start[q]; c < end[q]; c++) {
        new_value = (-1.0) * clusterDeriv[c] / sum_lld2[c];
        clusterDeriv[c] = new_value;
        if (fabs(new_value) > epsDeriv) {
          ok = true;
        }
      }

      // add the derivative to S
      for (i = 0; i < N; i++) {
        index = start[q] + dumMat(i, q);
        S[i] += clusterDeriv[index];
      }
    }
  }

  // Rprintf("other, nb iter=%i\n", iter);
  if (iter == iterMax) {
    Rprintf("[Getting cluster deriv. other] Max iterations reached (%i)\n",
            iterMax);
  }

  return (S);
}
[[cpp11::register]] doubles
cpp_partialDerivative_other_(int iterMax, int Q, int N, double epsDeriv,
                             doubles ll_d2, doubles dx_dother, doubles init,
                             integers_matrix<> dumMat, integers nbCluster) {
  // takes in:
  // dumMat: the matrix of dummies (n X c) each obs => cluster // must be in cpp
  // index!!! init: the initialisation of the sum of derivatives vector ll_d2:
  // the second derivative dx_dother: the vector of dx_dother

  int iter;

  int i, q, c;
  int index;
  int sum_cases = 0;
  bool ok;
  double new_value;

  std::vector<int> start, end;
  start.reserve(Q);
  end.reserve(Q);

  for (q = 0; q < Q; q++) {
    // the total number of clusters (eg if man/woman and 10 countries: total of
    // 12 cases)
    sum_cases += nbCluster[q];
    if (q == 0) {
      start[q] = 0;
      end[q] = nbCluster[q];
    } else {
      start[q] = start[q - 1] + nbCluster[q - 1];
      end[q] = end[q - 1] + nbCluster[q];
    }
  }

  // the derivatives per cluster, stacked in one vector
  writable::doubles clusterDeriv(sum_cases);
  for (int i = 0; i < sum_cases; i++) {
    clusterDeriv[i] = 0;
  }

  writable::doubles sum_lld2(sum_cases);
  for (int i = 0; i < sum_cases; i++) {
    sum_lld2[i] = 0;
  }

  // Creation of the sum_lld2
  for (i = 0; i < N; i++) {
    for (q = 0; q < Q; q++) {
      index = start[q] + dumMat(i, q);
      sum_lld2[index] += ll_d2[i];
    }
  }

  // the result to return
  writable::doubles S(N);
  for (i = 0; i < N; i++) {
    S[i] = init[i];
  }

  ok = true;
  iter = 0;
  while (ok & (iter < iterMax)) {
    iter++;
    ok = false;

    for (q = 0; q < Q; q++) {
      check_user_interrupt();

      // init of the vector
      for (c = start[q]; c < end[q]; c++) {
        clusterDeriv[c] = 0;
      }

      for (i = 0; i < N; i++) {
        index = start[q] + dumMat(i, q);
        clusterDeriv[index] += dx_dother[i] + S[i] * ll_d2[i];
      }

      // finish calculating clusterDeriv + control
      for (c = start[q]; c < end[q]; c++) {
        new_value = -1.0 * clusterDeriv[c] / sum_lld2[c];
        clusterDeriv[c] = new_value;
        if (fabs(new_value) > epsDeriv) {
          ok = true;
        }
      }

      // on ajoute la derivee a S:
      for (i = 0; i < N; i++) {
        index = start[q] + dumMat(i, q);
        S[i] += clusterDeriv[index];
      }
    }
  }

  // Rprintf("other, nb iter=%i\n", iter);
  if (iter == iterMax) {
    Rprintf("[Getting cluster deriv. other] Max iterations reached (%i)\n",
            iterMax);
  }

  return (S);
}
