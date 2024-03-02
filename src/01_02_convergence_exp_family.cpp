// list of objects that will be used to
// lighten the writting of the functions
// Contains only parameters that are fixed throguhout ALL the process

#include "01_0_convergence.hpp"

void CCC_poisson(int n_obs, int nb_cluster, double *cluster_coef,
                 double *exp_mu, double *sum_y, int *dum) {
  // compute cluster coef, poisson

  // initialize cluster coef
  for (int m = 0; m < nb_cluster; ++m) {
    cluster_coef[m] = 0;
  }

  // looping sequentially over exp_mu
  for (int i = 0; i < n_obs; ++i) {
    cluster_coef[dum[i]] += exp_mu[i];
  }

  // calculating cluster coef
  for (int m = 0; m < nb_cluster; ++m) {
    cluster_coef[m] = sum_y[m] / cluster_coef[m];
  }

  // "output" is the update of my_cluster_coef
}

void CCC_poisson_2(const std::vector<double> &pcluster_origin,
                   std::vector<double> &pcluster_destination, int n_i, int n_j,
                   int n_cells, const std::vector<int> &mat_row,
                   std::vector<int> &mat_col, std::vector<double> &mat_value,
                   const std::vector<double> &ca, const std::vector<double> &cb,
                   std::vector<double> &alpha) {
  double *beta = pcluster_destination.data() + n_i;

  std::fill(alpha.begin(), alpha.begin() + n_i, 0.0);
  std::fill(beta, beta + n_j, 0.0);

  for (int obs = 0; obs < n_cells; ++obs) {
    beta[mat_col[obs]] += mat_value[obs] * pcluster_origin[mat_row[obs]];
  }

  for (int j = 0; j < n_j; ++j) {
    beta[j] = cb[j] / beta[j];
  }

  for (int obs = 0; obs < n_cells; ++obs) {
    alpha[mat_row[obs]] += mat_value[obs] * beta[mat_col[obs]];
  }

  for (int i = 0; i < n_i; ++i) {
    pcluster_destination[i] = ca[i] / alpha[i];
  }
}

void CCC_poisson_log(int n_obs, int nb_cluster, double *cluster_coef,
                     double *mu, double *sum_y, int *dum) {
  // compute cluster coef, poisson
  // This is log poisson => we are there because classic poisson did not work
  // Thus high chance there are very high values of the cluster coefs (in abs
  // value) we need to take extra care in computing it => we apply trick of
  // substracting the max in the exp

  std::vector<double> mu_max(nb_cluster, 0.0);
  std::vector<bool> doInit(nb_cluster, true);

  // finding the max mu for each cluster
  int d;
  for (int i = 0; i < n_obs; ++i) {
    d = dum[i];
    if (doInit[d]) {
      mu_max[d] = mu[i];
      doInit[d] = false;
    } else if (mu[i] > mu_max[d]) {
      mu_max[d] = mu[i];
    }
  }

  // looping sequentially over exp_mu
  for (int i = 0; i < n_obs; ++i) {
    d = dum[i];
    cluster_coef[d] += exp(mu[i] - mu_max[d]);
  }

  // calculating cluster coef
  for (int m = 0; m < nb_cluster; ++m) {
    cluster_coef[m] = log(sum_y[m]) - log(cluster_coef[m]) - mu_max[m];
  }

  // "output" is the update of my_cluster_coef
}

// GAUSSIAN

void CCC_gaussian(int n_obs, int nb_cluster, double *cluster_coef, double *mu,
                  double *sum_y, int *dum, int *table) {
  // compute cluster coef, gaussian

  // initialize cluster coef
  std::fill_n(cluster_coef, nb_cluster, 0.0);

  // looping sequentially over mu
  for (int i = 0; i < n_obs; ++i) {
    cluster_coef[dum[i]] += mu[i];
  }

  // calculating cluster coef
  for (int m = 0; m < nb_cluster; ++m) {
    cluster_coef[m] = (sum_y[m] - cluster_coef[m]) / table[m];
  }

  // "output" is the update of my_cluster_coef
}

void CCC_gaussian_2(const std::vector<double> &pcluster_origin,
                    std::vector<double> &pcluster_destination, int n_i, int n_j,
                    int n_cells, int *mat_row, int *mat_col,
                    double *mat_value_Ab, double *mat_value_Ba,
                    const std::vector<double> &a_tilde,
                    std::vector<double> &beta) {
  std::copy(a_tilde.begin(), a_tilde.begin() + n_i,
            pcluster_destination.begin());
  std::fill(beta.begin(), beta.begin() + n_j, 0);

  for (int obs = 0; obs < n_cells; ++obs) {
    beta[mat_col[obs]] += mat_value_Ba[obs] * pcluster_origin[mat_row[obs]];
  }

  for (int obs = 0; obs < n_cells; ++obs) {
    pcluster_destination[mat_row[obs]] +=
        mat_value_Ab[obs] * beta[mat_col[obs]];
  }
}

// NEGATIVE BINOMIAL and LOGIT (helpers + negbin + logit)

// Helpers

void computeBounds(bool negbin, std::vector<double> &lower_bound,
                   std::vector<double> &upper_bound, const double *mu,
                   const int *cumtable, const int *obsCluster,
                   const double *sum_y, const int *table, int nb_cluster) {
  int u0;
  double value, mu_min, mu_max;

  for (int m = 0; m < nb_cluster; ++m) {
    u0 = (m == 0 ? 0 : cumtable[m - 1]);
    mu_min = mu[obsCluster[u0]];
    mu_max = mu[obsCluster[u0]];
    for (int u = 1 + u0; u < cumtable[m]; ++u) {
      value = mu[obsCluster[u]];
      if (value < mu_min) {
        mu_min = value;
      } else if (value > mu_max) {
        mu_max = value;
      }
    }

    if (negbin) {
      lower_bound[m] = log(sum_y[m]) - log(table[m]) - mu_max;
    } else {
      lower_bound[m] = log(sum_y[m]) - log(table[m] - sum_y[m]) - mu_max;
    }

    upper_bound[m] = lower_bound[m] + (mu_max - mu_min);
  }
}

double computeValue(bool negbin, int m, int u0, const int *cumtable,
                    const int *obsCluster, const double *mu, double x1,
                    const double *sum_y, double theta = 0.0,
                    const double *lhs = nullptr) {
  double value = sum_y[m];
  int i;
  if (negbin) {
    for (int u = u0; u < cumtable[m]; ++u) {
      i = obsCluster[u];
      value -= (theta + lhs[i]) / (1 + theta * exp(-x1 - mu[i]));
    }
  } else {
    for (int u = u0; u < cumtable[m]; ++u) {
      i = obsCluster[u];
      value -= 1 / (1 + exp(-x1 - mu[i]));
    }
  }
  return value;
}

double computeDerivative(bool negbin, int m, int u0, const int *cumtable,
                         const int *obsCluster, const double *mu, double x1,
                         double theta = 0, const double *lhs = nullptr) {
  double derivative = 0, exp_mu;
  int i;
  for (int u = u0; u < cumtable[m]; ++u) {
    if (negbin) {
      i = obsCluster[u];
      exp_mu = exp(x1 + mu[i]);
      derivative -=
          theta * (theta + lhs[i]) / ((theta / exp_mu + 1) * (theta + exp_mu));
    } else {
      exp_mu = exp(x1 + mu[obsCluster[u]]);
      derivative -= 1 / ((1 / exp_mu + 1) * (1 + exp_mu));
    }
  }
  return derivative;
}

// Negbin and Logit

void CCC_negbin(int nthreads, int nb_cluster, double theta, double diffMax_NR,
                double *cluster_coef, double *mu, double *lhs, double *sum_y,
                int *obsCluster, int *table, int *cumtable) {
  int iterMax = 100, iterFullDicho = 10;
  std::vector<double> lower_bound(nb_cluster), upper_bound(nb_cluster);

  // Use computeBounds to replace the initial bounds computation
  computeBounds(true, lower_bound, upper_bound, mu, cumtable, obsCluster, sum_y,
                table, nb_cluster);

#pragma omp parallel for num_threads(nthreads)
  for (int m = 0; m < nb_cluster; ++m) {
    double x1 = 0, lower = lower_bound[m], upper = upper_bound[m];
    bool keepGoing = true;
    int iter = 0, u0 = (m == 0 ? 0 : cumtable[m - 1]);
    double value, x0;

    if (x1 >= upper || x1 <= lower) {
      x1 = (lower + upper) / 2;
    }

    while (keepGoing) {
      ++iter;
      value = computeValue(true, m, u0, cumtable, obsCluster, mu, x1, sum_y,
                           theta, lhs);

      if (value > 0) {
        lower = x1;
      } else {
        upper = x1;
      }

      x0 = x1;
      if (value == 0 || iter <= iterFullDicho) {
        double derivative = computeDerivative(true, m, u0, cumtable, obsCluster,
                                              mu, x1, theta, lhs);
        x1 = x0 - value / derivative;
        if (x1 >= upper || x1 <= lower) {
          x1 = (lower + upper) / 2;
        }
      } else {
        x1 = (lower + upper) / 2;
      }

      if (iter == iterMax || stopping_criterion(x0, x1, diffMax_NR)) {
        keepGoing = false;
      }
    }
    cluster_coef[m] = x1;
  }
}

void CCC_logit(int nthreads, int nb_cluster, double diffMax_NR,
               double *cluster_coef, double *mu, double *sum_y, int *obsCluster,
               int *table, int *cumtable) {
  int iterMax = 100, iterFullDicho = 10;
  std::vector<double> lower_bound(nb_cluster), upper_bound(nb_cluster);

  // Use computeBounds to replace the initial bounds computation
  computeBounds(false, lower_bound, upper_bound, mu, cumtable, obsCluster,
                sum_y, table, nb_cluster);

#pragma omp parallel for num_threads(nthreads)
  for (int m = 0; m < nb_cluster; ++m) {
    double x1 = 0, lower = lower_bound[m], upper = upper_bound[m];
    bool keepGoing = true;
    int iter = 0, u0 = (m == 0 ? 0 : cumtable[m - 1]);
    double value, x0;

    if (x1 >= upper || x1 <= lower) {
      x1 = (lower + upper) / 2;
    }

    while (keepGoing) {
      ++iter;
      value = computeValue(false, m, u0, cumtable, obsCluster, mu, x1, sum_y);

      if (value > 0) {
        lower = x1;
      } else {
        upper = x1;
      }

      x0 = x1;
      if (value == 0 || iter <= iterFullDicho) {
        double derivative =
            computeDerivative(false, m, u0, cumtable, obsCluster, mu, x1);
        x1 = x0 - value / derivative;
        if (x1 >= upper || x1 <= lower) {
          x1 = (lower + upper) / 2;
        }
      } else {
        x1 = (lower + upper) / 2;
      }

      if (iter == iterMax) {
        keepGoing = false;
        Rprintf("[Getting cluster coefficients nber %i] max iterations reached "
                "(%i).\n",
                m, iterMax);
        Rprintf("Value Sum Deriv (NR) = %f. Difference = %f.\n", value,
                fabs(x0 - x1));
      }

      if (stopping_criterion(x0, x1, diffMax_NR)) {
        keepGoing = false;
      }
    }

    cluster_coef[m] = x1;
  }
}
