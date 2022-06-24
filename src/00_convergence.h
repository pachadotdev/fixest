/***************************************************************************
 * _________________                                                       *
 * || Convergence ||                                                       *
 * -----------------                                                       *
 *                                                                         *
 * Author: Laurent R. Berge                                                *
 *                                                                         *
 * Compute the optimal set of cluster coefficients based on Berge (2018),  *
 * in a ML framework.                                                      *
 *                                                                         *
 * The code concerns both the cluster coefficients (i.e. fixed-effects)    *
 * and the coefficients of the derivatives of the cluster coefficients.    *
 *                                                                         *
 * 2-FEs trick:                                                            *
 * For both the Poisson and the Gaussian methods I use a trick to hasten   *
 * computation. Instead of looping on the number of observations, the      *
 * loop is on the number of unique cases. For balanced panels, this is     *
 * not useful, but for strongly unbalanced panels, this makes a big        *
 * difference.                                                             *
 *                                                                         *
 * Logit/Negbin:                                                           *
 * To obtain the cluster coefficients for these two likelihoods, I use     *
 * a Newton-Raphson + dichotomy algorithm. This ensures convergence        *
 * even when the NR algo would go astray (which may be the case in some    *
 * situations).                                                            *
 *                                                                         *
 * Drawback:                                                               *
 * The big drawback of the ML methods is that, as opposed to GLM methods,  *
 * parallel computing cannot be leveraged.                                 *
 *                                                                         *
 **************************************************************************/

#pragma once

#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

using namespace cpp11;
using std::vector;
using std::fabs;

// HELPERS

bool continue_criterion(double a, double b, double diffMax);

bool stopping_criterion(double a, double b, double diffMax);

bool update_X_IronsTuck(int nb_coef_no_K, vector<double> &X,
						const vector<double> &GX, const vector<double> &GGX,
						vector<double> &delta_GX, vector<double> &delta2_X);

// STAT FAMILIES

struct PARAM_CCC
{
	int family;
	int n_obs;
	int K;
	double theta;
	double diffMax_NR;
	int nthreads;

	// vectors from R
	double *mu_init;
	int *pcluster;
	double *lhs;

	// vectors of pointers
	vector<int *> pdum;
	vector<int *> ptable;
	vector<double *> psum_y;
	vector<int *> pobsCluster;
	vector<int *> pcumtable;

	// value that will vary
	double *mu_with_coef;
};

void CCC_poisson(int n_obs, int nb_cluster,
				 double *cluster_coef, double *exp_mu,
				 double *sum_y, int *dum);

void CCC_poisson_2(const vector<double> &pcluster_origin, vector<double> &pcluster_destination,
				   int n_i, int n_j, int n_cells,
				   const vector<int> &mat_row, vector<int> &mat_col, vector<double> &mat_value,
				   const vector<double> &ca, const vector<double> &cb,
				   vector<double> &alpha);

void CCC_negbin(int nthreads, int nb_cluster, double theta, double diffMax_NR,
				double *cluster_coef, double *mu,
				double *lhs, double *sum_y, int *obsCluster, int *table, int *cumtable);

void CCC_logit(int nthreads, int nb_cluster, double diffMax_NR,
			   double *cluster_coef, double *mu,
			   double *sum_y, int *obsCluster, int *table, int *cumtable);

void CCC_gaussian(int n_obs, int nb_cluster,
				  double *cluster_coef, double *mu,
				  double *sum_y, int *dum, int *table);

void CCC_gaussian_2(const vector<double> &pcluster_origin, vector<double> &pcluster_destination,
					int n_i, int n_j, int n_cells,
					int *mat_row, int *mat_col,
					double *mat_value_Ab, double *mat_value_Ba,
					const vector<double> &a_tilde, vector<double> &beta);

void CCC_poisson_log(int n_obs, int nb_cluster,
					 double *cluster_coef, double *mu,
					 double *sum_y, int *dum);

// CLUSTERS

void computeClusterCoef_single(int family, int n_obs, int nb_cluster, double theta, double diffMax_NR,
							   double *cluster_coef, double *mu,
							   double *lhs, double *sum_y,
							   int *dum, int *obsCluster, int *table, int *cumtable, int nthreads);

void computeClusterCoef(vector<double *> &pcluster_origin, vector<double *> &pcluster_destination,
						PARAM_CCC *args);
                        