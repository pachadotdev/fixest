// list of objects that will be used to
// lighten the writting of the functions
// Contains only parameters that are fixed throguhout ALL the process

#include "01_0_convergence.h"

void CCC_poisson(int n_obs, int nb_cluster,
				 double *cluster_coef, double *exp_mu,
				 double *sum_y, int *dum)
{
	// compute cluster coef, poisson
	// Rprintf("in gaussian\n");

	// initialize cluster coef
	for (int m = 0; m < nb_cluster; ++m)
	{
		cluster_coef[m] = 0;
	}

	// looping sequentially over exp_mu
	for (int i = 0; i < n_obs; ++i)
	{
		cluster_coef[dum[i]] += exp_mu[i];
	}

	// calculating cluster coef
	for (int m = 0; m < nb_cluster; ++m)
	{
		cluster_coef[m] = sum_y[m] / cluster_coef[m];
	}

	// "output" is the update of my_cluster_coef
}

void CCC_poisson_2(const vector<double> &pcluster_origin, vector<double> &pcluster_destination,
				   int n_i, int n_j, int n_cells,
				   const vector<int> &mat_row, vector<int> &mat_col, vector<double> &mat_value,
				   const vector<double> &ca, const vector<double> &cb,
				   vector<double> &alpha)
{

	// alpha = ca / (Ab %m% (cb / (Ab %tm% alpha)))

	double *beta = pcluster_destination.data() + n_i;

	for (int i = 0; i < n_i; ++i)
	{
		alpha[i] = 0;
	}

	for (int j = 0; j < n_j; ++j)
	{
		beta[j] = 0;
	}

	for (int obs = 0; obs < n_cells; ++obs)
	{
		beta[mat_col[obs]] += mat_value[obs] * pcluster_origin[mat_row[obs]];
	}

	for (int j = 0; j < n_j; ++j)
	{
		beta[j] = cb[j] / beta[j];
	}

	for (int obs = 0; obs < n_cells; ++obs)
	{
		alpha[mat_row[obs]] += mat_value[obs] * beta[mat_col[obs]];
	}

	for (int i = 0; i < n_i; ++i)
	{
		pcluster_destination[i] = ca[i] / alpha[i];
	}
}

void CCC_poisson_log(int n_obs, int nb_cluster,
					 double *cluster_coef, double *mu,
					 double *sum_y, int *dum)
{
	// compute cluster coef, poisson
	// This is log poisson => we are there because classic poisson did not work
	// Thus high chance there are very high values of the cluster coefs (in abs value)
	// we need to take extra care in computing it => we apply trick of substracting the max in the exp

	vector<double> mu_max(nb_cluster);
	vector<bool> doInit(nb_cluster);

	// initialize cluster coef
	for (int m = 0; m < nb_cluster; ++m)
	{
		cluster_coef[m] = 0;
		doInit[m] = true;
	}

	// finding the max mu for each cluster
	int d;
	for (int i = 0; i < n_obs; ++i)
	{
		d = dum[i];
		if (doInit[d])
		{
			mu_max[d] = mu[i];
			doInit[d] = false;
		}
		else if (mu[i] > mu_max[d])
		{
			mu_max[d] = mu[i];
		}
	}

	// looping sequentially over exp_mu
	for (int i = 0; i < n_obs; ++i)
	{
		d = dum[i];
		cluster_coef[d] += exp(mu[i] - mu_max[d]);
	}

	// calculating cluster coef
	for (int m = 0; m < nb_cluster; ++m)
	{
		cluster_coef[m] = log(sum_y[m]) - log(cluster_coef[m]) - mu_max[m];
	}

	// "output" is the update of my_cluster_coef
}

// GAUSSIAN

void CCC_gaussian(int n_obs, int nb_cluster,
				  double *cluster_coef, double *mu,
				  double *sum_y, int *dum, int *table)
{
	// compute cluster coef, gaussian

	// initialize cluster coef
	for (int m = 0; m < nb_cluster; ++m)
	{
		cluster_coef[m] = 0;
	}

	// looping sequentially over mu
	for (int i = 0; i < n_obs; ++i)
	{
		cluster_coef[dum[i]] += mu[i];
	}

	// calculating cluster coef
	for (int m = 0; m < nb_cluster; ++m)
	{
		cluster_coef[m] = (sum_y[m] - cluster_coef[m]) / table[m];
	}

	// "output" is the update of my_cluster_coef
}

void CCC_gaussian_2(const vector<double> &pcluster_origin, vector<double> &pcluster_destination,
					int n_i, int n_j, int n_cells,
					int *mat_row, int *mat_col,
					double *mat_value_Ab, double *mat_value_Ba,
					const vector<double> &a_tilde, vector<double> &beta)
{

	// alpha = a_tilde + (Ab %m% (Ba %m% alpha))

	// for(int i=0 ; i<n_i ; ++i){
	// 	alpha[i] = 0;
	// }
	//
	// for(int j=0 ; j<n_j ; ++j){
	// 	beta[j] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_cells ; ++obs){
	// 	beta[mat_col[obs]] += mat_value_Ba[obs]*pcluster_origin[mat_row[obs]];
	// }
	//
	// for(int obs=0 ; obs<n_cells ; ++obs){
	// 	alpha[mat_row[obs]] += mat_value_Ab[obs]*beta[mat_col[obs]];
	// }
	//
	// for(int i=0 ; i<n_i ; ++i){
	// 	pcluster_destination[i] = a_tilde[i] + alpha[i];
	// }

	for (int i = 0; i < n_i; ++i)
	{
		pcluster_destination[i] = a_tilde[i];
	}

	for (int j = 0; j < n_j; ++j)
	{
		beta[j] = 0;
	}

	for (int obs = 0; obs < n_cells; ++obs)
	{
		beta[mat_col[obs]] += mat_value_Ba[obs] * pcluster_origin[mat_row[obs]];
	}

	for (int obs = 0; obs < n_cells; ++obs)
	{
		pcluster_destination[mat_row[obs]] += mat_value_Ab[obs] * beta[mat_col[obs]];
	}
}

// NEGATIVE BINOMIAL

void CCC_negbin(int nthreads, int nb_cluster, double theta, double diffMax_NR,
				double *cluster_coef, double *mu,
				double *lhs, double *sum_y, int *obsCluster, int *table, int *cumtable)
{
	// compute cluster coefficients negbin
	// This is not direct: needs to apply dichotomy+NR algorithm

	// first we find the min max for each cluster to get the bounds
	int iterMax = 100, iterFullDicho = 10;

	// finding the max/min values of mu for each cluster
	vector<double> borne_inf(nb_cluster);
	vector<double> borne_sup(nb_cluster);
	// attention borne_inf => quand mu est maximal
	int u0;
	double value, mu_min, mu_max;

	for (int m = 0; m < nb_cluster; ++m)
	{
		// the min/max of mu
		u0 = (m == 0 ? 0 : cumtable[m - 1]);
		mu_min = mu[obsCluster[u0]];
		mu_max = mu[obsCluster[u0]];
		for (int u = 1 + u0; u < cumtable[m]; ++u)
		{
			value = mu[obsCluster[u]];
			if (value < mu_min)
			{
				mu_min = value;
			}
			else if (value > mu_max)
			{
				mu_max = value;
			}
		}

		// computing the bounds
		borne_inf[m] = log(sum_y[m]) - log(static_cast<double>(table[m])) - mu_max;
		borne_sup[m] = borne_inf[m] + (mu_max - mu_min);
	}

	// Rprintf("inf: %f -- sup: %f -- middle: %f\n", borne_inf[0], borne_sup[0], (borne_inf[0] + borne_sup[0])/2);

    // TODO: OMP functions
    #pragma omp parallel for num_threads(nthreads)
	for (int m = 0; m < nb_cluster; ++m)
	{
		// we loop over each cluster

		// we initialise the cluster coefficient at 0 (it should converge to 0 at some point)
		double x1 = 0;
		bool keepGoing = true;
		int iter = 0, i;
		int u0 = (m == 0 ? 0 : cumtable[m - 1]);

		double value, x0, derivee = 0, exp_mu;

		// the bounds
		double lower_bound = borne_inf[m];
		double upper_bound = borne_sup[m];

		// Update of the value if it goes out of the boundaries
		// because we dont know ex ante if 0 is within the bounds
		if (x1 >= upper_bound || x1 <= lower_bound)
		{
			x1 = (lower_bound + upper_bound) / 2;
		}

		while (keepGoing)
		{
			++iter;

			// 1st step: initialisation des bornes

			// computing the value of f(x)
			value = sum_y[m];
			for (int u = u0; u < cumtable[m]; ++u)
			{
				i = obsCluster[u];
				value -= (theta + lhs[i]) / (1 + theta * exp(-x1 - mu[i]));
			}

			// update of the bounds.
			if (value > 0)
			{
				lower_bound = x1;
			}
			else
			{
				upper_bound = x1;
			}

			// 2nd step: NR iteration or Dichotomy
			x0 = x1;
			if (value == 0)
			{
				keepGoing = false;
			}
			else if (iter <= iterFullDicho)
			{
				// computing the derivative
				derivee = 0;
				for (int u = u0; u < cumtable[m]; ++u)
				{
					i = obsCluster[u];
					exp_mu = exp(x1 + mu[i]);
					derivee -= theta * (theta + lhs[i]) / ((theta / exp_mu + 1) * (theta + exp_mu));
				}

				x1 = x0 - value / derivee;
				// Rprintf("x1: %5.2f\n", x1);

				// 3rd step: dichotomy (if necessary)
				// Update of the value if it goes out of the boundaries
				if (x1 >= upper_bound || x1 <= lower_bound)
				{
					x1 = (lower_bound + upper_bound) / 2;
				}
			}
			else
			{
				x1 = (lower_bound + upper_bound) / 2;
			}

			// the stopping criterion
			if (iter == iterMax)
			{
				keepGoing = false;
				if (omp_get_thread_num() == 0)
				{
					Rprintf("[Getting cluster coefficients nber %i] max iterations reached (%i).\n", m, iterMax);
					Rprintf("Value Sum Deriv (NR) = %f. Difference = %f.\n", value, fabs(x0 - x1));
				}
			}

			// if(fabs(x0-x1) / (0.1 + fabs(x1)) < diffMax_NR){
			if (stopping_criterion(x0, x1, diffMax_NR))
			{
				keepGoing = false;
			}
		}

		// if(m == 0) Rprintf("diffMax = %f, value = %f, test = %f\n", fabs(x0-x1), value, fabs(5.1 - 5));
		// if(m == 0) Rprintf("total iter = %i, final x1 = %f\n", iter, x1);

		// after convegence: only update of cluster coef
		cluster_coef[m] = x1;
	}
}

void CCC_logit(int nthreads, int nb_cluster, double diffMax_NR,
			   double *cluster_coef, double *mu,
			   double *sum_y, int *obsCluster, int *table, int *cumtable)
{
	// compute cluster coefficients negbin
	// This is not direct: needs to apply dichotomy+NR algorithm

	// first we find the min max for each cluster to get the bounds
	int iterMax = 100, iterFullDicho = 10;

	// finding the max/min values of mu for each cluster
	vector<double> borne_inf(nb_cluster);
	vector<double> borne_sup(nb_cluster);
	// attention borne_inf => quand mu est maximal

	int u0;
	double value, mu_min, mu_max;
	for (int m = 0; m < nb_cluster; ++m)
	{
		// the min/max of mu
		u0 = (m == 0 ? 0 : cumtable[m - 1]);
		mu_min = mu[obsCluster[u0]];
		mu_max = mu[obsCluster[u0]];
		for (int u = 1 + u0; u < cumtable[m]; ++u)
		{
			value = mu[obsCluster[u]];
			if (value < mu_min)
			{
				mu_min = value;
			}
			else if (value > mu_max)
			{
				mu_max = value;
			}
		}

		// computing the "bornes"
		borne_inf[m] = log(sum_y[m]) - log(table[m] - sum_y[m]) - mu_max;
		borne_sup[m] = borne_inf[m] + (mu_max - mu_min);
	}

	//
	// Parallel loop
	//

    // TODO: OMP functions
    #pragma omp parallel for num_threads(nthreads)
	for (int m = 0; m < nb_cluster; ++m)
	{
		// we loop over each cluster

		// we initialise the cluster coefficient at 0 (it should converge to 0 at some point)
		double x1 = 0;
		bool keepGoing = true;
		int iter = 0;
		int u0 = (m == 0 ? 0 : cumtable[m - 1]);

		double value, x0, derivee = 0, exp_mu;

		// the bounds
		double lower_bound = borne_inf[m];
		double upper_bound = borne_sup[m];

		// Update of the value if it goes out of the boundaries
		// because we dont know ex ante if 0 is within the bounds
		if (x1 >= upper_bound || x1 <= lower_bound)
		{
			x1 = (lower_bound + upper_bound) / 2;
		}

		while (keepGoing)
		{
			++iter;

			// 1st step: initialisation des bornes

			// computing the value of f(x)
			value = sum_y[m];
			for (int u = u0; u < cumtable[m]; ++u)
			{
				value -= 1 / (1 + exp(-x1 - mu[obsCluster[u]]));
			}

			// update of the bounds.
			if (value > 0)
			{
				lower_bound = x1;
			}
			else
			{
				upper_bound = x1;
			}

			// 2nd step: NR iteration or Dichotomy
			x0 = x1;
			if (value == 0)
			{
				keepGoing = false;
			}
			else if (iter <= iterFullDicho)
			{
				// computing the derivative
				derivee = 0;
				for (int u = u0; u < cumtable[m]; ++u)
				{
					exp_mu = exp(x1 + mu[obsCluster[u]]);
					derivee -= 1 / ((1 / exp_mu + 1) * (1 + exp_mu));
				}

				x1 = x0 - value / derivee;
				// Rprintf("x1: %5.2f\n", x1);

				// 3rd step: dichotomy (if necessary)
				// Update of the value if it goes out of the boundaries
				if (x1 >= upper_bound || x1 <= lower_bound)
				{
					x1 = (lower_bound + upper_bound) / 2;
				}
			}
			else
			{
				x1 = (lower_bound + upper_bound) / 2;
			}

			// the stopping criteria
			if (iter == iterMax)
			{
				keepGoing = false;
				Rprintf("[Getting cluster coefficients nber %i] max iterations reached (%i).\n", m, iterMax);
				Rprintf("Value Sum Deriv (NR) = %f. Difference = %f.\n", value, fabs(x0 - x1));
			}

			// if(fabs(x0-x1) / (0.1 + fabs(x1)) < diffMax_NR){
			if (stopping_criterion(x0, x1, diffMax_NR))
			{
				keepGoing = false;
			}
		}
		// Rprintf("iter=%i.\n", iter);
		// res[k] = iter;
		// res[k] = value;

		// after convegence: only update of cluster coef
		cluster_coef[m] = x1;
	}
}
