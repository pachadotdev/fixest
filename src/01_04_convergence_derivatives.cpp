//
// Maintenant la convergence des derivees
//

#include "01_0_convergence.h"

[[cpp11::register]] list cpp_derivconv_seq_gnl(int iterMax, double diffMax, int n_vars, SEXP nb_cluster_all, SEXP ll_d2,
											   SEXP jacob_vector, SEXP deriv_init_vector, SEXP dum_vector)
{

	int n_obs = Rf_length(ll_d2);
	int K = Rf_length(nb_cluster_all);

	int *pcluster = INTEGER(nb_cluster_all);
	double *pll_d2 = REAL(ll_d2);

	int nb_coef = 0;
	for (int k = 0; k < K; ++k)
	{
		nb_coef += pcluster[k];
	}

	// Setting up the vectors on variables
	vector<double *> pjac(n_vars);
	pjac[0] = REAL(jacob_vector);
	for (int v = 1; v < n_vars; ++v)
	{
		pjac[v] = pjac[v - 1] + n_obs;
	}

	// setting up the vectors on clusters
	vector<int *> pdum(K);
	pdum[0] = INTEGER(dum_vector);

	for (int k = 1; k < K; ++k)
	{
		pdum[k] = pdum[k - 1] + n_obs;
	}

	// Setting up the target => deriv
	vector<double> deriv(n_obs * n_vars);
	double *my_init = REAL(deriv_init_vector);
	for (int i = 0; i < n_obs * n_vars; ++i)
	{
		deriv[i] = my_init[i];
	}
	// pointers to deriv
	vector<double *> pderiv(n_vars);
	pderiv[0] = deriv.data();
	for (int v = 1; v < n_vars; ++v)
	{
		pderiv[v] = pderiv[v - 1] + n_obs;
	}

	// the deriv coefficients & sum_ll_d2
	vector<double> deriv_coef(nb_coef);
	vector<double> sum_ll_d2(nb_coef);

	vector<double *> pderiv_coef(K);
	vector<double *> psum_ll_d2(K);

	pderiv_coef[0] = deriv_coef.data();
	psum_ll_d2[0] = sum_ll_d2.data();
	for (int k = 1; k < K; ++k)
	{
		pderiv_coef[k] = pderiv_coef[k - 1] + pcluster[k - 1];
		psum_ll_d2[k] = psum_ll_d2[k - 1] + pcluster[k - 1];
	}

	// computing the sum_ll_d2
	for (int m = 0; m < nb_coef; ++m)
	{
		// init
		sum_ll_d2[m] = 0;
	}

	for (int k = 0; k < K; ++k)
	{
		double *my_sum_ll_d2 = psum_ll_d2[k];
		int *my_dum = pdum[k];
		for (int i = 0; i < n_obs; ++i)
		{
			my_sum_ll_d2[my_dum[i]] += pll_d2[i];
		}
	}

	//
	// Loop
	//

	bool keepGoing = true;
	int iter = 0, iter_all_max = 0;

	// We first loop on each variable
	for (int v = 0; v < n_vars; ++v)
	{

		// loading the required data
		double *my_deriv = pderiv[v];
		double *my_jac = pjac[v];

		keepGoing = true;
		iter = 0;

		while (keepGoing && iter < iterMax)
		{
			++iter;
			keepGoing = false;

			// we update the clusters sequentially
			// we loop over all clusters => from K to 1
			for (int k = (K - 1); k >= 0; k--)
			{
				// Rprintf("k=%i ", k);
				R_CheckUserInterrupt();

				// loading the required info
				double *my_deriv_coef = pderiv_coef[k];
				int *my_dum = pdum[k];
				double *my_sum_ll_d2 = psum_ll_d2[k];
				int nb_cluster = pcluster[k];

				// init the deriv coef
				for (int m = 0; m < nb_cluster; ++m)
				{
					my_deriv_coef[m] = 0;
				}

				// sum the jac and deriv
				for (int i = 0; i < n_obs; ++i)
				{
					my_deriv_coef[my_dum[i]] += (my_jac[i] + my_deriv[i]) * pll_d2[i];
				}

				// divide by the LL sum
				for (int m = 0; m < nb_cluster; ++m)
				{
					my_deriv_coef[m] /= -my_sum_ll_d2[m];
				}

				// we update deriv
				for (int i = 0; i < n_obs; ++i)
				{
					my_deriv[i] += my_deriv_coef[my_dum[i]];
				}

				// the stopping criterion
				if (keepGoing == false)
				{
					for (int m = 0; m < nb_cluster; ++m)
					{
						if (fabs(my_deriv_coef[m]) > diffMax)
						{
							keepGoing = true;
							break;
						}
					}
				}
			}
		}

		if (iter > iter_all_max)
		{
			iter_all_max = iter;
		}
	}

	//
	// The result
	//

	writable::doubles_matrix<> dxi_dbeta(n_obs, n_vars);

	for (int v = 0; v < n_vars; ++v)
	{
		double *my_deriv = pderiv[v];
		for (int i = 0; i < n_obs; ++i)
		{
			dxi_dbeta(i, v) = my_deriv[i];
		}
	}

	writable::list res;
	res["dxi_dbeta"] = dxi_dbeta;
	res.push_back({"iter"_nm = iter_all_max});

	return (res);
}

// easier to handle in the acceleration
struct PARAM_DERIV_COEF
{
	int n_obs;
	int K;

	// vectors of pointers (length K)
	vector<int *> pdum;
	vector<double *> psum_ll_d2;
	vector<double *> psum_jac_lld2;

	// simple vectors using pointers
	int *pcluster;
	double *ll_d2;

	// only value that will vary
	double *deriv_with_coef; // => vector length n_obs
};

void computeDerivCoef(vector<double *> &pcoef_origin, vector<double *> &pcoef_destination,
					  double *my_deriv_init, PARAM_DERIV_COEF *args)
{

	//
	// Loading the arguments
	//

	int n_obs = args->n_obs;
	int K = args->K;
	vector<int *> &pdum = args->pdum;
	vector<double *> &psum_ll_d2 = args->psum_ll_d2;
	vector<double *> &psum_jac_lld2 = args->psum_jac_lld2;
	int *pcluster = args->pcluster;
	double *ll_d2 = args->ll_d2;
	double *deriv_with_coef = args->deriv_with_coef;

	// compute all the first two coefficients
	// thus you need to start by the last item

	for (int i = 0; i < n_obs; ++i)
	{
		deriv_with_coef[i] = my_deriv_init[i];
	}

	for (int k = 0; k < (K - 1); ++k)
	{
		int *my_dum = pdum[k];
		double *my_deriv_coef = pcoef_origin[k];

		for (int i = 0; i < n_obs; ++i)
		{
			deriv_with_coef[i] += my_deriv_coef[my_dum[i]];
		}
	}

	for (int k = K - 1; k >= 0; k--)
	{
		R_CheckUserInterrupt();

		// computing the optimal cluster coef -- given mu_with_coef
		double *my_deriv_coef = pcoef_destination[k];
		double *my_sum_ll_d2 = psum_ll_d2[k];
		double *my_sum_jac_lld2 = psum_jac_lld2[k];
		int *my_dum = pdum[k];
		int nb_cluster = pcluster[k];

		//
		// update of the deriv coefficients
		//

		// init the deriv coef
		for (int m = 0; m < nb_cluster; ++m)
		{
			my_deriv_coef[m] = my_sum_jac_lld2[m];
		}

		// sum the jac and deriv
		for (int i = 0; i < n_obs; ++i)
		{
			my_deriv_coef[my_dum[i]] += deriv_with_coef[i] * ll_d2[i];
		}

		// divide by the LL sum
		for (int m = 0; m < nb_cluster; ++m)
		{
			my_deriv_coef[m] /= -my_sum_ll_d2[m];
		}

		// updating the value of deriv_with_coef (only if necessary)
		if (k != 0)
		{

			for (int i = 0; i < n_obs; ++i)
			{
				deriv_with_coef[i] = my_deriv_init[i];
			}

			int *my_dum;
			double *my_deriv_coef;
			for (int h = 0; h < K; h++)
			{
				if (h == k - 1)
					continue;

				my_dum = pdum[h];

				if (h < k - 1)
				{
					my_deriv_coef = pcoef_origin[h];
				}
				else
				{
					my_deriv_coef = pcoef_destination[h];
				}

				for (int i = 0; i < n_obs; ++i)
				{
					deriv_with_coef[i] += my_deriv_coef[my_dum[i]];
				}
			}
		}
	}
}

[[cpp11::register]] list cpp_derivconv_acc_gnl(int iterMax, double diffMax, int n_vars, SEXP nb_cluster_all, SEXP ll_d2,
											   SEXP jacob_vector, SEXP deriv_init_vector, SEXP dum_vector)
{

	int n_obs = Rf_length(ll_d2);
	int K = Rf_length(nb_cluster_all);

	int *pcluster = INTEGER(nb_cluster_all);
	double *pll_d2 = REAL(ll_d2);

	int nb_coef = 0;
	for (int k = 0; k < K; ++k)
	{
		nb_coef += pcluster[k];
	}

	// Setting up the vectors on variables
	vector<double *> pjac(n_vars);
	pjac[0] = REAL(jacob_vector);
	for (int v = 1; v < n_vars; ++v)
	{
		pjac[v] = pjac[v - 1] + n_obs;
	}

	// setting up the vectors on clusters
	vector<int *> pdum(K);
	pdum[0] = INTEGER(dum_vector);
	for (int k = 1; k < K; ++k)
	{
		pdum[k] = pdum[k - 1] + n_obs;
	}

	// pointers to deriv_init_vector (length n_vars*n_obs)
	vector<double *> pderiv_init(n_vars);
	pderiv_init[0] = REAL(deriv_init_vector);
	for (int v = 1; v < n_vars; ++v)
	{
		pderiv_init[v] = pderiv_init[v - 1] + n_obs;
	}

	// the deriv coefficients & sum_ll_d2 & sum_jac_lld2
	vector<double> deriv_coef(nb_coef);
	vector<double> sum_ll_d2(nb_coef);
	vector<double> sum_jac_lld2(nb_coef);
	vector<double *> pderiv_coef(K);
	vector<double *> psum_ll_d2(K);
	vector<double *> psum_jac_lld2(K);
	pderiv_coef[0] = deriv_coef.data();
	psum_ll_d2[0] = sum_ll_d2.data();
	psum_jac_lld2[0] = sum_jac_lld2.data();
	for (int k = 1; k < K; ++k)
	{
		pderiv_coef[k] = pderiv_coef[k - 1] + pcluster[k - 1];
		psum_ll_d2[k] = psum_ll_d2[k - 1] + pcluster[k - 1];
		psum_jac_lld2[k] = psum_jac_lld2[k - 1] + pcluster[k - 1];
	}

	// computing the sum_ll_d2
	for (int m = 0; m < nb_coef; ++m)
	{
		// init
		sum_ll_d2[m] = 0;
	}

	for (int k = 0; k < K; ++k)
	{
		double *my_sum_ll_d2 = psum_ll_d2[k];
		int *my_dum = pdum[k];
		for (int i = 0; i < n_obs; ++i)
		{
			my_sum_ll_d2[my_dum[i]] += pll_d2[i];
		}
	}

	vector<double> deriv_with_coef(n_obs);

	//
	// Sending the information
	//

	PARAM_DERIV_COEF args;
	args.n_obs = n_obs;
	args.K = K;
	args.pdum = pdum;
	args.psum_ll_d2 = psum_ll_d2;
	args.psum_jac_lld2 = psum_jac_lld2;
	args.pcluster = pcluster;
	args.ll_d2 = pll_d2;
	args.deriv_with_coef = deriv_with_coef.data();

	//
	// IT iteration (preparation)
	//

	// variables on 1:K
	vector<double> X(nb_coef);
	vector<double> GX(nb_coef);
	vector<double> GGX(nb_coef);
	// pointers:
	vector<double *> pX(K);
	vector<double *> pGX(K);
	vector<double *> pGGX(K);

	pX[0] = X.data();
	pGX[0] = GX.data();
	pGGX[0] = GGX.data();

	for (int k = 1; k < K; ++k)
	{
		pX[k] = pX[k - 1] + pcluster[k - 1];
		pGX[k] = pGX[k - 1] + pcluster[k - 1];
		pGGX[k] = pGGX[k - 1] + pcluster[k - 1];
	}

	// variables on 1:(K-1)
	int nb_coef_no_K = 0;
	for (int k = 0; k < (K - 1); ++k)
	{
		nb_coef_no_K += pcluster[k];
	}
	vector<double> delta_GX(nb_coef_no_K);
	vector<double> delta2_X(nb_coef_no_K);

	//
	// Loop
	//

	writable::doubles_matrix<> dxi_dbeta(n_obs, n_vars);

	bool keepGoing = true;
	int iter = 0, iter_all_max = 0;

	// We first loop on each variable
	for (int v = 0; v < n_vars; ++v)
	{

		// loading the required data
		double *my_deriv_init = pderiv_init[v];
		double *my_jac = pjac[v];

		//
		// we update the values of sum_jac_lld2
		//

		for (int k = 0; k < K; ++k)
		{
			int *my_dum = pdum[k];
			double *my_sum_jac_lld2 = psum_jac_lld2[k];

			for (int m = 0; m < pcluster[k]; ++m)
			{
				my_sum_jac_lld2[m] = 0;
			}

			// sum the jac and deriv
			for (int i = 0; i < n_obs; ++i)
			{
				my_sum_jac_lld2[my_dum[i]] += my_jac[i] * pll_d2[i];
			}
		}

		//
		// The IT loop
		//

		// we initialize the deriv coefficients
		for (int m = 0; m < nb_coef; ++m)
		{
			X[m] = 0;
		}

		computeDerivCoef(pX, pGX, my_deriv_init, &args);

		keepGoing = true;
		iter = 0;
		bool numconv = false;

		while (keepGoing && iter < iterMax)
		{
			++iter;

			// origin: GX, destination: GGX
			computeDerivCoef(pGX, pGGX, my_deriv_init, &args);

			// X ; update of the cluster coefficient
			numconv = update_X_IronsTuck(nb_coef_no_K, X, GX, GGX, delta_GX, delta2_X);
			if (numconv)
				break;

			// origin: X, destination: GX
			computeDerivCoef(pX, pGX, my_deriv_init, &args);

			keepGoing = false;
			// the stopping criterion
			for (int m = 0; m < nb_coef_no_K; ++m)
			{
				// if(fabs(X[m] - GX[m]) / (0.1 + fabs(GX[m])) > diffMax){
				if (continue_criterion(X[m], GX[m], diffMax))
				{
					// Rprintf("diffMax = %f\n", fabs(X[m] - GX[m]));
					keepGoing = true;
					break;
				}
			}
		}

		// save iteration info
		if (iter > iter_all_max)
		{
			iter_all_max = iter;
		}

		//
		// Save the deriv into res
		//

		// we compute the deriv based on the cluster coefs
		for (int i = 0; i < n_obs; ++i)
		{
			deriv_with_coef[i] = my_deriv_init[i];
		}

		for (int k = 0; k < K; ++k)
		{
			int *my_dum = pdum[k];
			double *my_deriv_coef = pGX[k];
			for (int i = 0; i < n_obs; ++i)
			{
				deriv_with_coef[i] += my_deriv_coef[my_dum[i]];
			}
		}

		// save
		for (int i = 0; i < n_obs; ++i)
		{
			dxi_dbeta(i, v) = deriv_with_coef[i];
		}
	}

	writable::list res;
	res["dxi_dbeta"] = dxi_dbeta;
	res.push_back({"iter"_nm = iter_all_max});

	return (res);
}

void computeDerivCoef_2(vector<double> &alpha_origin, vector<double> &alpha_destination,
						int n_i, int n_j, int n_cells,
						const vector<double> &a_tilde,
						const vector<int> &mat_row, const vector<int> &mat_col,
						const vector<double> &mat_value_Ab, const vector<double> &mat_value_Ba,
						vector<double> &beta)
{

	// a_tile + Ab * Ba * alpha

	for (int m = 0; m < n_i; ++m)
	{
		alpha_destination[m] = a_tilde[m];
	}

	for (int m = 0; m < n_j; ++m)
	{
		beta[m] = 0;
	}

	for (int obs = 0; obs < n_cells; ++obs)
	{
		beta[mat_col[obs]] += mat_value_Ba[obs] * alpha_origin[mat_row[obs]];
	}

	for (int obs = 0; obs < n_cells; ++obs)
	{
		alpha_destination[mat_row[obs]] += mat_value_Ab[obs] * beta[mat_col[obs]];
	}
}

[[cpp11::register]] list cpp_derivconv_acc_2(int iterMax, double diffMax, int n_vars, SEXP nb_cluster_all,
											 int n_cells, SEXP index_i, SEXP index_j, SEXP ll_d2, SEXP order,
											 SEXP jacob_vector, SEXP deriv_init_vector, SEXP dum_vector)
{

	int n_obs = Rf_length(ll_d2);

	int *pcluster = INTEGER(nb_cluster_all);
	double *pll_d2 = REAL(ll_d2);
	int n_i = pcluster[0], n_j = pcluster[1];

	// Setting up the vectors on variables
	vector<double *> pjac(n_vars);
	pjac[0] = REAL(jacob_vector);
	for (int v = 1; v < n_vars; ++v)
	{
		pjac[v] = pjac[v - 1] + n_obs;
	}

	// setting up the vectors on clusters
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;

	// pointers to deriv_init_vector (length n_vars*n_obs)
	vector<double *> pderiv_init(n_vars);
	pderiv_init[0] = REAL(deriv_init_vector);
	for (int v = 1; v < n_vars; ++v)
	{
		pderiv_init[v] = pderiv_init[v - 1] + n_obs;
	}

	// the sum_ll_d2
	vector<double> sum_ll_d2_i(n_i);
	vector<double> sum_ll_d2_j(n_j);

	// computing the sum_ll_d2
	for (int m = 0; m < n_i; ++m)
	{
		sum_ll_d2_i[m] = 0;
	}

	for (int m = 0; m < n_j; ++m)
	{
		sum_ll_d2_j[m] = 0;
	}

	for (int obs = 0; obs < n_obs; ++obs)
	{
		sum_ll_d2_i[dum_i[obs]] += pll_d2[obs];
		sum_ll_d2_j[dum_j[obs]] += pll_d2[obs];
	}

	//
	// Setting up the matrices A and B, and vector a_tilde
	//

	vector<int> mat_row(n_cells);
	vector<int> mat_col(n_cells);
	vector<double> mat_value_Ab(n_cells);
	vector<double> mat_value_Ba(n_cells);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);
	int *new_order = INTEGER(order);

	int index_current = 0;

	double value_Ab = pll_d2[new_order[0]];
	double value_Ba = pll_d2[new_order[0]];

	for (int obs = 1; obs < n_obs; ++obs)
	{
		if (pindex_j[obs] != pindex_j[obs - 1] || pindex_i[obs] != pindex_i[obs - 1])
		{
			// save the value if we change index
			mat_row[index_current] = pindex_i[obs - 1];
			mat_col[index_current] = pindex_j[obs - 1];
			mat_value_Ab[index_current] = value_Ab / -sum_ll_d2_i[pindex_i[obs - 1]];
			mat_value_Ba[index_current] = value_Ba / -sum_ll_d2_j[pindex_j[obs - 1]];

			// new row
			index_current++;

			// the new value
			value_Ab = pll_d2[new_order[obs]];
			value_Ba = pll_d2[new_order[obs]];
		}
		else
		{
			value_Ab += pll_d2[new_order[obs]];
			value_Ba += pll_d2[new_order[obs]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs - 1];
	mat_col[index_current] = pindex_j[n_obs - 1];
	mat_value_Ab[index_current] = value_Ab / -sum_ll_d2_i[pindex_i[n_obs - 1]];
	;
	mat_value_Ba[index_current] = value_Ba / -sum_ll_d2_j[pindex_j[n_obs - 1]];
	;

	//
	// IT iteration (preparation)
	//

	vector<double> X(n_i);
	vector<double> GX(n_i);
	vector<double> GGX(n_i);
	vector<double> delta_GX(n_i);
	vector<double> delta2_X(n_i);
	vector<double> beta(n_j);
	vector<double> alpha_final(n_i);
	vector<double> beta_final(n_j);

	//
	// Loop
	//

	writable::doubles_matrix<> dxi_dbeta(n_obs, n_vars);

	bool keepGoing = true;
	int iter = 0, iter_all_max = 0;

	// We first loop on each variable
	for (int v = 0; v < n_vars; ++v)
	{

		// loading the required data
		double *my_deriv_init = pderiv_init[v];
		double *my_jac = pjac[v];

		//
		// we update the constants and then alpha tilde
		//

		// we compute the constants
		vector<double> a(n_i);
		vector<double> b(n_j);

		for (int m = 0; m < n_i; ++m)
		{
			a[m] = 0;
		}
		for (int m = 0; m < n_j; ++m)
		{
			b[m] = 0;
		}
		for (int obs = 0; obs < n_obs; ++obs)
		{
			a[dum_i[obs]] += (my_jac[obs] + my_deriv_init[obs]) * pll_d2[obs];
			b[dum_j[obs]] += (my_jac[obs] + my_deriv_init[obs]) * pll_d2[obs];
		}
		for (int m = 0; m < n_i; ++m)
		{
			a[m] /= -sum_ll_d2_i[m];
		}
		for (int m = 0; m < n_j; ++m)
		{
			b[m] /= -sum_ll_d2_j[m];
		}

		// a_tilde: a_tilde = a + (Ab %m% b)
		vector<double> a_tilde(n_i);
		for (int m = 0; m < n_i; ++m)
		{
			a_tilde[m] = a[m];
		}
		for (int i = 0; i < n_cells; ++i)
		{
			a_tilde[mat_row[i]] += mat_value_Ab[i] * b[mat_col[i]];
		}

		//
		// The IT loop
		//

		// we initialize the deriv coefficients
		for (int m = 0; m < n_i; ++m)
		{
			X[m] = 0;
		}

		computeDerivCoef_2(X, GX, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);

		keepGoing = true;
		iter = 0;
		bool numconv = false;

		while (keepGoing && iter < iterMax)
		{
			++iter;

			// origin: GX, destination: GGX
			computeDerivCoef_2(GX, GGX, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);

			// X ; update of the cluster coefficient
			numconv = update_X_IronsTuck(n_i, X, GX, GGX, delta_GX, delta2_X);
			if (numconv)
				break;

			// origin: X, destination: GX
			computeDerivCoef_2(X, GX, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);

			keepGoing = false;
			// the stopping criterion
			for (int m = 0; m < n_i; ++m)
			{
				// if(fabs(X[m] - GX[m]) / (0.1 + fabs(GX[m])) > diffMax){
				if (continue_criterion(X[m], GX[m], diffMax))
				{
					keepGoing = true;
					break;
				}
			}
		}

		// save iteration info
		if (iter > iter_all_max)
		{
			iter_all_max = iter;
		}

		//
		// Save the deriv into res
		//

		// we compute the last alpha and beta
		for (int m = 0; m < n_i; ++m)
		{
			alpha_final[m] = a[m];
		}

		for (int m = 0; m < n_j; ++m)
		{
			beta_final[m] = b[m];
		}

		for (int obs = 0; obs < n_cells; ++obs)
		{
			beta_final[mat_col[obs]] += mat_value_Ba[obs] * GX[mat_row[obs]];
		}

		for (int obs = 0; obs < n_cells; ++obs)
		{
			alpha_final[mat_row[obs]] += mat_value_Ab[obs] * beta_final[mat_col[obs]];
		}

		// save
		for (int obs = 0; obs < n_obs; ++obs)
		{
			dxi_dbeta(obs, v) = my_deriv_init[obs] + alpha_final[dum_i[obs]] + beta_final[dum_j[obs]];
		}
	}

	writable::list res;
	res["dxi_dbeta"] = dxi_dbeta;
	res.push_back({"iter"_nm = iter_all_max});

	return (res);
}

[[cpp11::register]] list cpp_derivconv_seq_2(int iterMax, double diffMax, int n_vars, SEXP nb_cluster_all,
											 int n_cells, SEXP index_i, SEXP index_j, SEXP order, SEXP ll_d2,
											 SEXP jacob_vector, SEXP deriv_init_vector, SEXP dum_vector)
{

	int n_obs = Rf_length(ll_d2);

	int *pcluster = INTEGER(nb_cluster_all);
	double *pll_d2 = REAL(ll_d2);
	int n_i = pcluster[0], n_j = pcluster[1];

	// Rprintf("n_i = %i ; n_j = %i ; n_obs = %i ; n_cells = %i\n", n_i, n_j, n_obs, n_cells);

	// Setting up the vectors on variables
	vector<double *> pjac(n_vars);
	pjac[0] = REAL(jacob_vector);
	for (int v = 1; v < n_vars; ++v)
	{
		pjac[v] = pjac[v - 1] + n_obs;
	}

	// setting up the vectors on clusters
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;

	// pointers to deriv_init_vector (length n_vars*n_obs)
	vector<double *> pderiv_init(n_vars);
	pderiv_init[0] = REAL(deriv_init_vector);
	for (int v = 1; v < n_vars; ++v)
	{
		pderiv_init[v] = pderiv_init[v - 1] + n_obs;
	}

	// the sum_ll_d2
	vector<double> sum_ll_d2_i(n_i);
	vector<double> sum_ll_d2_j(n_j);

	// computing the sum_ll_d2
	for (int m = 0; m < n_i; ++m)
	{
		sum_ll_d2_i[m] = 0;
	}

	for (int m = 0; m < n_j; ++m)
	{
		sum_ll_d2_j[m] = 0;
	}

	for (int obs = 0; obs < n_obs; ++obs)
	{
		sum_ll_d2_i[dum_i[obs]] += pll_d2[obs];
		sum_ll_d2_j[dum_j[obs]] += pll_d2[obs];
	}

	//
	// Setting up the matrices A and B, and vector a_tilde
	//

	vector<int> mat_row(n_cells);
	vector<int> mat_col(n_cells);
	vector<double> mat_value_Ab(n_cells);
	vector<double> mat_value_Ba(n_cells);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);
	int *new_order = INTEGER(order);

	int index_current = 0;

	double value_current = pll_d2[new_order[0]];

	for (int obs = 1; obs < n_obs; ++obs)
	{
		if (pindex_j[obs] != pindex_j[obs - 1] || pindex_i[obs] != pindex_i[obs - 1])
		{
			// save the value if we change index
			mat_row[index_current] = pindex_i[obs - 1];
			mat_col[index_current] = pindex_j[obs - 1];
			mat_value_Ab[index_current] = value_current / -sum_ll_d2_i[pindex_i[obs - 1]];
			mat_value_Ba[index_current] = value_current / -sum_ll_d2_j[pindex_j[obs - 1]];

			// new row
			index_current++;

			// the new value
			value_current = pll_d2[new_order[obs]];
		}
		else
		{
			value_current += pll_d2[new_order[obs]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs - 1];
	mat_col[index_current] = pindex_j[n_obs - 1];
	mat_value_Ab[index_current] = value_current / -sum_ll_d2_i[pindex_i[n_obs - 1]];
	mat_value_Ba[index_current] = value_current / -sum_ll_d2_j[pindex_j[n_obs - 1]];

	//
	// IT iteration (preparation)
	//

	vector<double> X(n_i);
	vector<double> X_new(n_i);
	double *X_final;
	vector<double> beta(n_j);
	vector<double> alpha_final(n_i);
	vector<double> beta_final(n_j);

	//
	// Loop
	//

	writable::doubles_matrix<> dxi_dbeta(n_obs, n_vars);

	bool keepGoing = true;
	int iter = 0, iter_all_max = 0;

	// We first loop on each variable
	for (int v = 0; v < n_vars; ++v)
	{

		// loading the required data
		double *my_deriv_init = pderiv_init[v];
		double *my_jac = pjac[v];

		//
		// we update the constants and then alpha tilde
		//

		// we compute the constants
		vector<double> a(n_i);
		vector<double> b(n_j);
		for (int m = 0; m < n_i; ++m)
		{
			a[m] = 0;
		}
		for (int m = 0; m < n_j; ++m)
		{
			b[m] = 0;
		}
		for (int obs = 0; obs < n_obs; ++obs)
		{
			a[dum_i[obs]] += (my_jac[obs] + my_deriv_init[obs]) * pll_d2[obs];
			b[dum_j[obs]] += (my_jac[obs] + my_deriv_init[obs]) * pll_d2[obs];
		}
		for (int m = 0; m < n_i; ++m)
		{
			a[m] /= -sum_ll_d2_i[m];
		}
		for (int m = 0; m < n_j; ++m)
		{
			b[m] /= -sum_ll_d2_j[m];
		}

		// a_tilde: a_tilde = a + (Ab %m% b)
		vector<double> a_tilde(n_i);
		for (int m = 0; m < n_i; ++m)
		{
			a_tilde[m] = a[m];
		}
		for (int i = 0; i < n_cells; ++i)
		{
			a_tilde[mat_row[i]] += mat_value_Ab[i] * b[mat_col[i]];
		}

		// init of X (alpha)
		for (int m = 0; m < n_i; ++m)
		{
			X[m] = 0;
		}

		keepGoing = true;
		iter = 0;

		while (keepGoing && iter < iterMax)
		{
			++iter;

			if (iter % 2 == 1)
			{
				computeDerivCoef_2(X, X_new, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);
			}
			else
			{
				computeDerivCoef_2(X_new, X, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);
			}

			keepGoing = false;
			// the stopping criterion
			for (int m = 0; m < n_i; ++m)
			{
				// if(fabs(X[m] - X_new[m]) / (0.1 + fabs(X_new[m])) > diffMax){
				if (continue_criterion(X[m], X_new[m], diffMax))
				{
					keepGoing = true;
					break;
				}
			}
		}

		// save iteration info
		if (iter > iter_all_max)
		{
			iter_all_max = iter;
		}

		//
		// Save the deriv into res
		//

		X_final = (iter % 2 == 1) ? X_new.data() : X.data();

		// we compute the last alpha and beta
		for (int m = 0; m < n_i; ++m)
		{
			alpha_final[m] = a[m];
		}

		for (int m = 0; m < n_j; ++m)
		{
			beta_final[m] = b[m];
		}

		for (int obs = 0; obs < n_cells; ++obs)
		{
			beta_final[mat_col[obs]] += mat_value_Ba[obs] * X_final[mat_row[obs]];
		}

		for (int obs = 0; obs < n_cells; ++obs)
		{
			alpha_final[mat_row[obs]] += mat_value_Ab[obs] * beta_final[mat_col[obs]];
		}

		// save
		for (int obs = 0; obs < n_obs; ++obs)
		{
			dxi_dbeta(obs, v) = my_deriv_init[obs] + alpha_final[dum_i[obs]] + beta_final[dum_j[obs]];
		}
	}

	writable::list res;
	res["dxi_dbeta"] = dxi_dbeta;
	res.push_back({"iter"_nm = iter_all_max});

	return (res);
}

[[cpp11::register]] doubles_matrix<> update_deriv_single(int n_vars, int nb_coef,
														 SEXP r_ll_d2, SEXP r_jacob_vector, SEXP r_dum_vector)
{

	int n_obs = Rf_length(r_ll_d2);

	// loading variables
	double *ll_d2 = REAL(r_ll_d2);
	int *dum = INTEGER(r_dum_vector);

	vector<double *> pjac(n_vars);
	pjac[0] = REAL(r_jacob_vector);
	for (int v = 1; v < n_vars; ++v)
	{
		pjac[v] = pjac[v - 1] + n_obs;
	}

	// vector sum_ll_d2
	vector<double> sum_ll_d2(nb_coef, 0);
	for (int obs = 0; obs < n_obs; ++obs)
	{
		sum_ll_d2[dum[obs]] += ll_d2[obs];
	}

	// the vector of coefficients
	vector<double> coef_deriv(nb_coef);

	// the result
	writable::doubles_matrix<> res(n_obs, n_vars); // init at 0

	for (int v = 0; v < n_vars; ++v)
	{
		double *my_jac = pjac[v];

		//
		// 1) we compute the coef
		//

		for (int m = 0; m < nb_coef; ++m)
		{
			coef_deriv[m] = 0;
		}

		for (int obs = 0; obs < n_obs; ++obs)
		{
			coef_deriv[dum[obs]] += my_jac[obs] * ll_d2[obs];
		}

		// divide is costly, we do it only nb_coef times
		for (int m = 0; m < nb_coef; ++m)
		{
			coef_deriv[m] /= -sum_ll_d2[m];
		}

		//
		// 2) We save the value in res
		//

		for (int obs = 0; obs < n_obs; ++obs)
		{
			res(obs, v) = coef_deriv[dum[obs]];
		}
	}

	return (res);
}
