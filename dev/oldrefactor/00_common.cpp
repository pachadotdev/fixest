#include "00_common.hpp"

#include "01_0_convergence.hpp"
#include "02_0_demeaning.hpp"
#include "03_0_dot_square_bracket.hpp"
#include "04_0_linear_model.hpp"
#include "07_0_parallel.hpp"
#include "08_0_strings.hpp"
#include "09_0_vcov.hpp"

bool continue_criterion(double a, double b, double diffMax) {
  // continuing criterion of the algorithm
  double diff = fabs(a - b);
  return ((diff > diffMax) && (diff / (0.1 + fabs(a)) > diffMax));
}

bool stopping_criterion(double a, double b, double diffMax) {
  // stopping criterion of the algorithm
  double diff = fabs(a - b);
  return ((diff < diffMax) || (diff / (0.1 + fabs(a)) < diffMax));
}

// IT update + returns numerical convergence indicator
// nb_coef_no_KQ = K or Q
bool update_X_IronsTuck(int nb_coef_no_KQ, vector<double> &X,
                        const vector<double> &GX, const vector<double> &GGX,
                        vector<double> &delta_GX, vector<double> &delta2_X) {
  double vprod = 0, ssq = 0;

  for (int i = 0; i < nb_coef_no_KQ; ++i) {
    delta_GX[i] = GGX[i] - GX[i];
    delta2_X[i] = delta_GX[i] - GX[i] + X[i];

    vprod += delta_GX[i] * delta2_X[i];
    ssq += delta2_X[i] * delta2_X[i];
  }

  if (ssq == 0) {
    return true;
  }

  double coef = vprod / ssq;

  // update of X:
  for (int i = 0; i < nb_coef_no_KQ; ++i) {
    X[i] = GGX[i] - coef * delta_GX[i];
  }

  return false;
}

// => this concerns only the parallel application on a 1-Dimensional matrix
// takes in the nber of observations of the vector and the nber of threads
// gives back a vector of the length the nber of threads + 1 giving the
// start/stop of each threads
vector<int> set_parallel_scheme(int N, int nthreads = 1) {
  vector<int> res(nthreads + 1, 0);
  double N_rest = N;

  for (int i = 0; i < nthreads; ++i) {
    int division = ceil(N_rest / (nthreads - i));
    N_rest -= division;
    res[i + 1] = division + res[i];
  }

  return res;
}
