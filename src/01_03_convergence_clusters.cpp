#include "01_0_convergence.hpp"

void computeClusterCoef_single(int family, int n_obs, int nb_cluster,
                               double theta, double diffMax_NR,
                               double *cluster_coef, double *mu, double *lhs,
                               double *sum_y, int *dum, int *obsCluster,
                               int *table, int *cumtable, int nthreads) {
  // Leads to the appropriate function
  // we update the cluster "in place" (ie using pointers)
  // switch is peculiar... we need break too

  // comment on mutlithreading:
  // I know I could use only the CCC_par_negbin/CCC_par_logit functions but I'm
  // scared that in some instance openmp creates some bug, thus when single
  // core, there is no reference to openmp at all

  switch (family) {
  case 1:
    CCC_poisson(n_obs, nb_cluster, cluster_coef, mu, sum_y, dum);
    break;
  case 2: // Negbin
    CCC_negbin(nthreads, nb_cluster, theta, diffMax_NR, cluster_coef, mu, lhs,
               sum_y, obsCluster, table, cumtable);
    break;
  case 3: // logit
    CCC_logit(nthreads, nb_cluster, diffMax_NR, cluster_coef, mu, sum_y,
              obsCluster, table, cumtable);
    break;
  case 4: // Gaussian
    CCC_gaussian(n_obs, nb_cluster, cluster_coef, mu, sum_y, dum, table);
    break;
  case 5: // log poisson
    CCC_poisson_log(n_obs, nb_cluster, cluster_coef, mu, sum_y, dum);
    break;
  }
}

void updateMuWithCoef(int family, int n_obs, double *mu_with_coef,
                      double *my_cluster_coef, int *my_dum) {
  auto operation = (family == 1) ? [](double &mu, double coef) { mu *= coef; }
                                 : [](double &mu, double coef) { mu += coef; };

  for (int i = 0; i < n_obs; ++i) {
    operation(mu_with_coef[i], my_cluster_coef[my_dum[i]]);
  }
}

void computeClusterCoef(std::vector<double *> &pcluster_origin,
                        std::vector<double *> &pcluster_destination,
                        PARAM_CCC *args) {
  // update of the cluster coefficients
  // first we update mu, then we update the cluster coefficicents

  //
  // Loading the variables
  //

  int family = args->family;
  int n_obs = args->n_obs;
  int K = args->K;
  int nthreads = args->nthreads;
  double theta = args->theta;
  double diffMax_NR = args->diffMax_NR;

  int *pcluster = args->pcluster;
  double *lhs = args->lhs;
  double *mu_init = args->mu_init;

  std::vector<int *> &pdum = args->pdum;
  std::vector<int *> &ptable = args->ptable;
  std::vector<double *> &psum_y = args->psum_y;
  std::vector<int *> &pobsCluster = args->pobsCluster;
  std::vector<int *> &pcumtable = args->pcumtable;

  // value that will be modified
  double *mu_with_coef = args->mu_with_coef;

  // We update each cluster coefficient, starting from K

  // we first set the value of mu_with_coef
  memcpy(mu_with_coef, mu_init, n_obs * sizeof(double));

  for (int k = 0; k < (K - 1); ++k) {
    updateMuWithCoef(family, n_obs, mu_with_coef, pcluster_origin[k],
                     args->pdum[k]);
  }

  for (int k = K - 1; k >= 0; k--) {
    check_user_interrupt();

    // computing the optimal cluster coef -- given mu_with_coef
    double *my_cluster_coef = pcluster_destination[k];
    int *my_table = ptable[k];
    double *my_sum_y = psum_y[k];
    int *my_dum = pdum[k];
    int *my_cumtable = pcumtable[k];
    int *my_obsCluster = pobsCluster[k];
    int nb_cluster = pcluster[k];

    // update of the cluster coefficients
    computeClusterCoef_single(family, n_obs, nb_cluster, theta, diffMax_NR,
                              my_cluster_coef, mu_with_coef, lhs, my_sum_y,
                              my_dum, my_obsCluster, my_table, my_cumtable,
                              nthreads);

    // updating the value of mu_with_coef (only if necessary)
    if (k != 0) {
      memcpy(mu_with_coef, mu_init, n_obs * sizeof(double));

      for (int h = 0; h < K; h++) {
        if (h == k - 1) {
          continue;
        }

        double *my_cluster_coef =
            (h < k - 1) ? pcluster_origin[h] : pcluster_destination[h];
        updateMuWithCoef(family, n_obs, mu_with_coef, my_cluster_coef,
                         args->pdum[h]);
      }
    }
  }

  // In the end, the array pcluster_coef is fully updated, starting from K to 1
}

[[cpp11::register]] SEXP
update_mu_single_cluster_(int family, int nb_cluster, double theta,
                          double diffMax_NR, SEXP mu_in, SEXP lhs, SEXP sum_y,
                          SEXP dum, SEXP obsCluster, SEXP table, SEXP cumtable,
                          int nthreads) {
  // Function used to compute the cluster coefficient for ONE fixed-effect only
  // and then add it to the existing mu

  int n_obs = Rf_length(mu_in);

  // Getting the pointers
  int *pdum = INTEGER(dum);
  int *pobsCluster = INTEGER(obsCluster);
  int *ptable = INTEGER(table);
  int *pcumtable = INTEGER(cumtable);
  double *plhs = REAL(lhs);
  double *psum_y = REAL(sum_y);
  double *pmu_in = REAL(mu_in);

  // The cluster coeffficients
  std::vector<double> cluster_coef(nb_cluster);

  computeClusterCoef_single(family, n_obs, nb_cluster, theta, diffMax_NR,
                            cluster_coef.data(), pmu_in, plhs, psum_y, pdum,
                            pobsCluster, ptable, pcumtable, nthreads);

  // the result
  SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
  double *pmu = REAL(mu);

  if (family == 1) {
    // Poisson, exp form
    for (int i = 0; i < n_obs; ++i) {
      pmu[i] = pmu_in[i] * cluster_coef[pdum[i]];
    }
  } else {
    for (int i = 0; i < n_obs; ++i) {
      pmu[i] = pmu_in[i] + cluster_coef[pdum[i]];
    }
  }

  UNPROTECT(1);

  return mu;
}

[[cpp11::register]] int get_n_cells_(integers index_i, integers index_j) {
  int n = index_i.size();

  // we count the nber of different elements
  int index_current = 0;

  for (int i = 1; i < n; ++i) {
    if (index_j[i] != index_j[i - 1] || index_i[i] != index_i[i - 1]) {
      // new row
      index_current++;
    }
  }

  index_current++;

  return index_current;
}
