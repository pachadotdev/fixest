#include "01_0_convergence.hpp"

void computeClusterCoef_single(int family, int n_obs, int nb_cluster,
                               double theta, double diffMax_NR,
                               double *cluster_coef, double *mu, double *lhs,
                               double *sum_y, int *dum, int *obsCluster,
                               int *table, int *cumtable, int nthreads = 1) {
  // Leads to the appropriate function
  // we update the cluster "in place" (ie using pointers)

  // comment on mutlithreading:
  // I know I could use only the CCC_par_negbin/CCC_par_logit functions but I'm
  // scared that in some instance openmp creates some bug, thus when single
  // core, there is no reference to openmp at all

  switch (family) {
    case 1:
      CCC_poisson(n_obs, nb_cluster, cluster_coef, mu, sum_y, dum);
      break;
    case 2:  // Negbin
      CCC_negbin(nthreads, nb_cluster, theta, diffMax_NR, cluster_coef, mu, lhs,
                 sum_y, obsCluster, table, cumtable);
      break;
    case 3:  // logit
      CCC_logit(nthreads, nb_cluster, diffMax_NR, cluster_coef, mu, sum_y,
                obsCluster, table, cumtable);
      break;
    case 4:  // Gaussian
      CCC_gaussian(n_obs, nb_cluster, cluster_coef, mu, sum_y, dum, table);
      break;
    case 5:  // log poisson
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

void computeClusterCoef(vector<double *> &pcluster_origin,
                        vector<double *> &pcluster_destination,
                        PARAM_CCC *args) {
  // Loading the variables
  int family = args->family;
  int n_obs = args->n_obs;
  int K = args->K;
  double *mu_with_coef = args->mu_with_coef;
  double *mu_init = args->mu_init;

  // we first set the value of mu_with_coef
  memcpy(mu_with_coef, mu_init, n_obs * sizeof(double));

  for (int k = 0; k < (K - 1); ++k) {
    updateMuWithCoef(family, n_obs, mu_with_coef, pcluster_origin[k],
                     args->pdum[k]);
  }

  for (int k = K - 1; k >= 0; k--) {
    R_CheckUserInterrupt();

    // computing the optimal cluster coef -- given mu_with_coef
    computeClusterCoef_single(
        family, n_obs, args->pcluster[k], args->theta, args->diffMax_NR,
        pcluster_destination[k], mu_with_coef, args->lhs, args->psum_y[k],
        args->pdum[k], args->pobsCluster[k], args->ptable[k],
        args->pcumtable[k], args->nthreads);

    // updating the value of mu_with_coef (only if necessary)
    if (k != 0) {
      memcpy(mu_with_coef, mu_init, n_obs * sizeof(double));

      for (int h = 0; h < K; h++) {
        if (h == k - 1) continue;

        double *my_cluster_coef =
            (h < k - 1) ? pcluster_origin[h] : pcluster_destination[h];
        updateMuWithCoef(family, n_obs, mu_with_coef, my_cluster_coef,
                         args->pdum[h]);
      }
    }
  }
}

struct ClusterData {
  double *mu;
  double *lhs;
  double *sum_y;
  int *dum;
  int *obsCluster;
  int *table;
  int *cumtable;
};

[[cpp11::register]] SEXP compute_cluster_coef_r_(
    int family, int nb_coef, double theta, double diffMax_NR, SEXP r_mu,
    SEXP r_lhs, SEXP r_sum_y, SEXP r_dum, SEXP r_obsCluster, SEXP r_table,
    SEXP r_cumtable, int nthreads = 1) {
  int n_obs = Rf_length(r_mu);

  ClusterData data = {
      REAL(r_mu),         REAL(r_lhs),           REAL(r_sum_y),
      INTEGER(r_dum),     INTEGER(r_obsCluster), INTEGER(r_table),
      INTEGER(r_cumtable)};

  SEXP res = PROTECT(Rf_allocVector(REALSXP, nb_coef));
  double *pcoef = REAL(res);

  computeClusterCoef_single(family, n_obs, nb_coef, theta, diffMax_NR, pcoef,
                            data.mu, data.lhs, data.sum_y, data.dum,
                            data.obsCluster, data.table, data.cumtable,
                            nthreads);

  UNPROTECT(1);

  return res;
}

[[cpp11::register]] SEXP update_mu_single_cluster_(
    int family, int nb_cluster, double theta, double diffMax_NR, SEXP mu_in,
    SEXP lhs, SEXP sum_y, SEXP dum, SEXP obsCluster, SEXP table, SEXP cumtable,
    int nthreads = 1) {
  // Function used to compute the cluster coefficient for ONE fixed-effect only
  // and then add it to the existing mu

  int n_obs = Rf_length(mu_in);

  ClusterData data = {REAL(mu_in),      REAL(lhs),           REAL(sum_y),
                      INTEGER(dum),     INTEGER(obsCluster), INTEGER(table),
                      INTEGER(cumtable)};

  // The cluster coeffficients
  vector<double> cluster_coef(nb_cluster);

  computeClusterCoef_single(family, n_obs, nb_cluster, theta, diffMax_NR,
                            cluster_coef.data(), data.mu, data.lhs, data.sum_y,
                            data.dum, data.obsCluster, data.table,
                            data.cumtable, nthreads);

  // the result
  SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
  double *pmu = REAL(mu);

  if (family == 1) {
    // Poisson, exp form
    for (int i = 0; i < n_obs; ++i) {
      pmu[i] = data.mu[i] * cluster_coef[data.dum[i]];
    }
  } else {
    for (int i = 0; i < n_obs; ++i) {
      pmu[i] = data.mu[i] + cluster_coef[data.dum[i]];
    }
  }

  UNPROTECT(1);

  return mu;
}

[[cpp11::register]] int get_n_cells_(integers index_i, integers index_j) {
  int n = index_i.size();

  // we count the number of different elements
  int index_current = 0;

  for (int i = 1; i < n; ++i) {
    if (index_j[i] != index_j[i - 1] || index_i[i] != index_i[i - 1]) {
      // new row
      index_current++;
    }
  }

  return (index_current + 1);
}
