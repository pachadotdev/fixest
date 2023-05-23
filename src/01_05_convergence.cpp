#include "01_0_convergence.hpp"

[[cpp11::register]] list cpp_conv_acc_gnl_(int family, int iterMax, double diffMax, double diffMax_NR, double theta, SEXP nb_cluster_all,
										  SEXP lhs, SEXP mu_init, SEXP dum_vector, SEXP tableCluster_vector,
										  SEXP sum_y_vector, SEXP cumtable_vector, SEXP obsCluster_vector, int nthreads = 1)
{

	// initial variables
	int K = Rf_length(nb_cluster_all);
	int *pcluster = INTEGER(nb_cluster_all);
	int n_obs = Rf_length(mu_init);
	double *pmu_init = REAL(mu_init);

	int nb_coef = 0;
	for (int k = 0; k < K; ++k)
	{
		nb_coef += pcluster[k];
	}

	// setting the pointers
	// table, coef, sum_y: length nb_coef
	// dum_all: length K*n_obs
	vector<int *> ptable(K);
	vector<double *> psum_y(K);
	ptable[0] = INTEGER(tableCluster_vector);
	psum_y[0] = REAL(sum_y_vector);

	for (int k = 1; k < K; ++k)
	{
		ptable[k] = ptable[k - 1] + pcluster[k - 1];
		psum_y[k] = psum_y[k - 1] + pcluster[k - 1];
	}

	// cumtable (points to nothing for non negbin/logit families)
	vector<int *> pcumtable(K);
	if (family == 2 || family == 3)
	{
		pcumtable[0] = INTEGER(cumtable_vector);
		for (int k = 1; k < K; ++k)
		{
			pcumtable[k] = pcumtable[k - 1] + pcluster[k - 1];
		}
	}

	// obsCluster (points to nothing for non negbin/logit families)
	vector<int *> pobsCluster(K);
	if (family == 2 || family == 3)
	{
		pobsCluster[0] = INTEGER(obsCluster_vector);
		for (int k = 1; k < K; ++k)
		{
			pobsCluster[k] = pobsCluster[k - 1] + n_obs;
		}
	}

	// cluster of each observation
	vector<int *> pdum(K);
	pdum[0] = INTEGER(dum_vector);
	for (int k = 1; k < K; ++k)
	{
		pdum[k] = pdum[k - 1] + n_obs;
	}

	// lhs (only negbin will use it)
	double *plhs = REAL(lhs);

	//
	// Sending variables to envir
	//

	PARAM_CCC args;

	args.family = family;
	args.n_obs = n_obs;
	args.K = K;
	args.nthreads = nthreads;
	args.theta = (family == 2 ? theta : 1); // theta won't be used if family not negbin
	args.diffMax_NR = diffMax_NR;
	args.pdum = pdum;
	args.mu_init = pmu_init;
	args.ptable = ptable;
	args.psum_y = psum_y;
	args.pcluster = pcluster;
	args.pcumtable = pcumtable;
	args.pobsCluster = pobsCluster;
	args.lhs = plhs;

	// value that will be modified
	vector<double> mu_with_coef(n_obs);
	args.mu_with_coef = mu_with_coef.data();

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
	// the main loop
	//

	// initialisation of X and then pGX
	if (family == 1)
	{
		for (int i = 0; i < nb_coef; ++i)
		{
			X[i] = 1;
		}
	}
	else
	{
		for (int i = 0; i < nb_coef; ++i)
		{
			X[i] = 0;
		}
	}

	// first iteration
	computeClusterCoef(pX, pGX, &args);

	// Rprintf("init: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", GX[i]);
	// }
	// Rprintf("\n");

	// stop("kkk");

	// flag for problem with poisson
	bool any_negative_poisson = false;

	// check whether we should go into the loop
	bool keepGoing = false;
	for (int i = 0; i < nb_coef; ++i)
	{
		// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
		if (continue_criterion(X[i], GX[i], diffMax))
		{
			keepGoing = true;
			break;
		}
	}

	int iter = 0;
	bool numconv = false;
	while (keepGoing && iter < iterMax)
	{
		++iter;

		// GGX -- origin: GX, destination: GGX
		computeClusterCoef(pGX, pGGX, &args);

		// Rprintf("iter: %i ; ", iter);
		//
		// if(iter >= iterMax - 3){
		// 	Rprintf("GGX: ");
		// 	for(int i=0 ; i<8 ; ++i){
		// 		Rprintf("%2.7f ", GGX[i]);
		// 	}
		// 	Rprintf("\n");
		// }

		// X ; update of the cluster coefficient
		numconv = update_X_IronsTuck(nb_coef_no_K, X, GX, GGX, delta_GX, delta2_X);
		if (numconv)
			break;

		// if(iter >= iterMax - 3){
		// 	Rprintf(" IT: ");
		// 	for(int i=0 ; i<8 ; ++i){
		// 		Rprintf("%2.7f ", X[i]);
		// 	}
		// 	Rprintf("\n");
		// }

		if (family == 1)
		{
			// We control for possible problems with poisson
			for (int i = 0; i < nb_coef_no_K; ++i)
			{
				if (X[i] <= 0)
				{
					any_negative_poisson = true;
					break;
				}
			}

			if (any_negative_poisson)
			{
				break; // we quit the loop
					   // update of mu is OK, it's only IT iteration that leads to negative values
			}
		}

		// GX -- origin: X, destination: GX
		computeClusterCoef(pX, pGX, &args);

		keepGoing = false;
		for (int i = 0; i < nb_coef_no_K; ++i)
		{
			// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
			if (continue_criterion(X[i], GX[i], diffMax))
			{
				keepGoing = true;
				break;
			}
		}
	}

	//
	// We update mu => result
	//

	SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(mu);
	for (int i = 0; i < n_obs; ++i)
	{
		pmu[i] = pmu_init[i];
	}

	// // regular code (using last iteration)
	// for(int k=0 ; k<K ; ++k){
	// 	int *my_dum = pdum[k];
	// 	double *my_cluster_coef = pGX[k]; // GX is last updated
	// 	if(family == 1){ // "quick" Poisson case
	// 		for(int i=0 ; i<n_obs ; ++i){
	// 			pmu[i] *= my_cluster_coef[my_dum[i]];
	// 		}
	// 	} else {
	// 		for(int i=0 ; i<n_obs ; ++i){
	// 			pmu[i] += my_cluster_coef[my_dum[i]];
	// 		}
	// 	}
	// }

	// To obtain stg identical to the R code:
	computeClusterCoef(pGX, pGGX, &args);
	for (int k = 0; k < K; ++k)
	{
		int *my_dum = pdum[k];
		double *my_cluster_coef = pGGX[k];

		if (family == 1)
		{ // "quick" Poisson case
			for (int i = 0; i < n_obs; ++i)
			{
				pmu[i] *= my_cluster_coef[my_dum[i]];
			}
		}
		else
		{
			for (int i = 0; i < n_obs; ++i)
			{
				pmu[i] += my_cluster_coef[my_dum[i]];
			}
		}
	}

	UNPROTECT(1);

	writable::list res;
	res.push_back({"mu_new"_nm = mu});
	res.push_back({"iter"_nm = iter});
	res.push_back({"any_negative_poisson"_nm = any_negative_poisson});

	return (res);
}

[[cpp11::register]] list cpp_conv_seq_gnl_(int family, int iterMax, double diffMax, double diffMax_NR, double theta, SEXP nb_cluster_all,
										  SEXP lhs, SEXP mu_init, SEXP dum_vector, SEXP tableCluster_vector,
										  SEXP sum_y_vector, SEXP cumtable_vector, SEXP obsCluster_vector, int nthreads = 1)
{

	// initial variables
	int K = Rf_length(nb_cluster_all);
	int *pcluster = INTEGER(nb_cluster_all);
	int n_obs = Rf_length(mu_init);
	double *pmu_init = REAL(mu_init);

	int nb_coef = 0;
	for (int k = 0; k < K; ++k)
	{
		nb_coef += pcluster[k];
	}

	// setting the pointers
	// table, coef, sum_y: length nb_coef
	// dum_all: length K*n_obs
	vector<int *> ptable(K);
	vector<double *> psum_y(K);
	ptable[0] = INTEGER(tableCluster_vector);
	psum_y[0] = REAL(sum_y_vector);

	for (int k = 1; k < K; ++k)
	{
		ptable[k] = ptable[k - 1] + pcluster[k - 1];
		psum_y[k] = psum_y[k - 1] + pcluster[k - 1];
	}

	// cumtable (points to nothing for non negbin/logit families)
	vector<int *> pcumtable(K);
	if (family == 2 || family == 3)
	{
		pcumtable[0] = INTEGER(cumtable_vector);
		for (int k = 1; k < K; ++k)
		{
			pcumtable[k] = pcumtable[k - 1] + pcluster[k - 1];
		}
	}

	// obsCluster (points to nothing for non negbin/logit families)
	vector<int *> pobsCluster(K);
	if (family == 2 || family == 3)
	{
		pobsCluster[0] = INTEGER(obsCluster_vector);
		for (int k = 1; k < K; ++k)
		{
			pobsCluster[k] = pobsCluster[k - 1] + n_obs;
		}
	}

	// cluster of each observation
	vector<int *> pdum(K);
	pdum[0] = INTEGER(dum_vector);
	for (int k = 1; k < K; ++k)
	{
		pdum[k] = pdum[k - 1] + n_obs;
	}

	// lhs (only negbin will use it)
	double *plhs = REAL(lhs);

	// the variable with the value of mu
	vector<double> mu_with_coef(n_obs);
	// initialization of mu_with_coef
	for (int i = 0; i < n_obs; ++i)
	{
		mu_with_coef[i] = pmu_init[i];
	}

	// The cluster coefficients
	// variables on 1:K
	vector<double> cluster_coef(nb_coef);
	// pointers:
	vector<double *> pcluster_coef(K);
	pcluster_coef[0] = cluster_coef.data();
	for (int k = 1; k < K; ++k)
	{
		pcluster_coef[k] = pcluster_coef[k - 1] + pcluster[k - 1];
	}

	//
	// the main loop
	//

	// initialisation of the cluster coefficients
	if (family == 1)
	{ // Poisson
		for (int i = 0; i < nb_coef; ++i)
		{
			cluster_coef[i] = 1;
		}
	}
	else
	{
		for (int i = 0; i < nb_coef; ++i)
		{
			cluster_coef[i] = 0;
		}
	}

	// the main loop
	bool keepGoing = true;
	int iter = 1;
	while (keepGoing && iter <= iterMax)
	{
		// Rprintf("\n%i.", iter);
		++iter;
		keepGoing = false;

		/// we loop over all clusters => from K to 1
		for (int k = (K - 1); k >= 0; k--)
		{
			// Rprintf("k=%i ", k);
			R_CheckUserInterrupt();

			//
			// 1) computing the cluster coefficient
			//

			// computing the optimal cluster coef -- given mu_with_coef
			double *my_cluster_coef = pcluster_coef[k];
			int *my_table = ptable[k];
			double *my_sum_y = psum_y[k];
			int *my_dum = pdum[k];
			int *my_cumtable = pcumtable[k];
			int *my_obsCluster = pobsCluster[k];
			int nb_cluster = pcluster[k];

			// update of the cluster coefficients
			computeClusterCoef_single(family, n_obs, nb_cluster, theta, diffMax_NR,
									  my_cluster_coef, mu_with_coef.data(), plhs, my_sum_y,
									  my_dum, my_obsCluster, my_table, my_cumtable, nthreads);

			//
			// 2) Updating the value of mu
			//

			// we add the computed value
			if (family == 1)
			{
				for (int i = 0; i < n_obs; ++i)
				{
					mu_with_coef[i] *= my_cluster_coef[my_dum[i]];
				}
			}
			else
			{
				for (int i = 0; i < n_obs; ++i)
				{
					mu_with_coef[i] += my_cluster_coef[my_dum[i]];
				}
			}

			// Stopping criterion
			if (keepGoing == false)
			{

				// Rprintf("nb: %i ", nb_cluster);
				// Rprintf("value: %f", my_cluster_coef[nb_cluster]);
				// stop("hhhh");

				if (family == 1)
				{
					for (int m = 0; m < nb_cluster; ++m)
					{
						if (fabs(my_cluster_coef[m] - 1) > diffMax)
						{
							keepGoing = true;
							break;
						}
					}
				}
				else
				{
					for (int m = 0; m < nb_cluster; ++m)
					{
						if (fabs(my_cluster_coef[m]) > diffMax)
						{
							keepGoing = true;
							break;
						}
					}
				}
			}
		}
	}

	//
	// We update mu => result
	//

	SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(mu);
	for (int i = 0; i < n_obs; ++i)
	{
		pmu[i] = mu_with_coef[i];
	}

	UNPROTECT(1);

	writable::list res;
	res.push_back({"mu_new"_nm = mu});
	res.push_back({"iter"_nm = iter});

	return (res);
}

[[cpp11::register]] list cpp_conv_acc_poi_2_(int n_i, int n_j, int n_cells, SEXP index_i, SEXP index_j,
											SEXP dum_vector, SEXP sum_y_vector,
											int iterMax, double diffMax, SEXP exp_mu_in, SEXP order)
{

	// values that will be used later
	vector<double> alpha(n_i);

	// We compute the matrix Ab
	int index_current = 0;
	vector<int> mat_row(n_cells);
	vector<int> mat_col(n_cells);
	vector<double> mat_value(n_cells);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);

	int n_obs = Rf_length(exp_mu_in);

	int *new_order = INTEGER(order);
	double *pexp_mu_in = REAL(exp_mu_in);
	// double value = pexp_mu_in[0];
	double value = pexp_mu_in[new_order[0]];

	for (int i = 1; i < n_obs; ++i)
	{
		if (pindex_j[i] != pindex_j[i - 1] || pindex_i[i] != pindex_i[i - 1])
		{
			// save the value if we change index
			mat_row[index_current] = pindex_i[i - 1];
			mat_col[index_current] = pindex_j[i - 1];
			mat_value[index_current] = value;

			// new row
			index_current++;

			// the new value
			value = pexp_mu_in[new_order[i]];
		}
		else
		{
			value += pexp_mu_in[new_order[i]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs - 1];
	mat_col[index_current] = pindex_j[n_obs - 1];
	mat_value[index_current] = value;

	//
	// IT iteration (preparation)
	//

	vector<double> X(n_i + n_j);
	vector<double> GX(n_i + n_j);
	vector<double> GGX(n_i + n_j);
	vector<double> delta_GX(n_i);
	vector<double> delta2_X(n_i);

	//
	// the main loop
	//

	// initialisation of X and then GX
	for (int i = 0; i < n_i; ++i)
	{
		X[i] = 1;
	}

	// ca and cb
	double *psum_y = REAL(sum_y_vector);

	vector<double> ca(n_i);
	vector<double> cb(n_j);
	for (int i = 0; i < n_i; ++i)
	{
		ca[i] = psum_y[i];
	}

	double *psum_y_j = psum_y + n_i;
	for (int j = 0; j < n_j; ++j)
	{
		cb[j] = psum_y_j[j];
	}

	// first iteration
	CCC_poisson_2(X, GX, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);

	// Rprintf("init: ");
	// for(i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", GX[i]);
	// }
	// Rprintf("\n");

	// stop("kkk");

	// flag for problem with poisson
	bool any_negative_poisson = false;

	// check whether we should go into the loop => always 1 iter
	bool keepGoing = true;
	for (int i = 0; i < n_i; ++i)
	{
		// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
		if (continue_criterion(X[i], GX[i], diffMax))
		{
			keepGoing = true;
			break;
		}
	}

	bool numconv = false;
	int iter = 0;
	while (keepGoing && iter < iterMax)
	{
		// Rprintf("%i.", iter);
		++iter;

		// GGX -- origin: GX, destination: GGX
		CCC_poisson_2(GX, GGX, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);

		// Rprintf("ggx: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f -- ", GGX[i]);
		// }
		// Rprintf("\n");

		// X ; update of the cluster coefficient
		numconv = update_X_IronsTuck(n_i, X, GX, GGX, delta_GX, delta2_X);
		if (numconv)
			break;

		// Control for negative values
		for (int i = 0; i < n_i; ++i)
		{
			if (X[i] <= 0)
			{
				any_negative_poisson = true;
				break;
			}
		}

		if (any_negative_poisson)
		{
			break; // we quit the loop
				   // update of mu is OK, it's only IT iteration that leads to negative values
		}

		// Rprintf("  x: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f -- ", X[i]);
		// }
		// Rprintf("\n");

		// GX -- origin: X, destination: GX
		CCC_poisson_2(X, GX, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);

		// Rprintf("   gx: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f\t", GX[i]);
		// }
		// Rprintf("\n");

		keepGoing = false;
		for (int i = 0; i < n_i; ++i)
		{
			// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
			if (continue_criterion(X[i], GX[i], diffMax))
			{
				keepGoing = true;
				break;
			}
		}
	}

	SEXP exp_mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(exp_mu);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;

	// double *beta = GX + n_i;
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	pmu[obs] = pexp_mu_in[obs] * GX[dum_i[obs]] * beta[dum_j[obs]];
	// }

	// pour avoir identique a acc_pois
	CCC_poisson_2(GX, X, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);
	double *beta = X.data() + n_i;
	for (int obs = 0; obs < n_obs; ++obs)
	{
		pmu[obs] = pexp_mu_in[obs] * X[dum_i[obs]] * beta[dum_j[obs]];
	}

	// Rprintf("alpha: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f\t", X[i]);
	// }
	// Rprintf("\n");

	UNPROTECT(1);

	writable::list res;
	res.push_back({"mu_new"_nm = exp_mu});
	res.push_back({"iter"_nm = iter});
	res.push_back({"any_negative_poisson"_nm = any_negative_poisson});

	return (res);
}

[[cpp11::register]] list cpp_conv_seq_poi_2_(int n_i, int n_j, int n_cells, SEXP index_i, SEXP index_j,
											SEXP dum_vector, SEXP sum_y_vector,
											int iterMax, double diffMax, SEXP exp_mu_in, SEXP order)
{

	// values that will be used later
	vector<double> alpha(n_i);

	// We compute the matrix Ab
	int index_current = 0;
	vector<int> mat_row(n_cells);
	vector<int> mat_col(n_cells);
	vector<double> mat_value(n_cells);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);

	int n_obs = Rf_length(exp_mu_in);

	int *new_order = INTEGER(order);
	double *pexp_mu_in = REAL(exp_mu_in);
	// double value = pexp_mu_in[0];
	double value = pexp_mu_in[new_order[0]];

	for (int i = 1; i < n_obs; ++i)
	{
		if (pindex_j[i] != pindex_j[i - 1] || pindex_i[i] != pindex_i[i - 1])
		{
			// save the value if we change index
			mat_row[index_current] = pindex_i[i - 1];
			mat_col[index_current] = pindex_j[i - 1];
			mat_value[index_current] = value;

			// new row
			index_current++;

			// the new value
			value = pexp_mu_in[new_order[i]];
		}
		else
		{
			value += pexp_mu_in[new_order[i]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs - 1];
	mat_col[index_current] = pindex_j[n_obs - 1];
	mat_value[index_current] = value;

	// X, X_new => the vector of coefficients
	vector<double> X_new(n_i + n_j);
	vector<double> X(n_i + n_j);

	//
	// the main loop
	//

	// initialisation of X and then GX
	for (int i = 0; i < n_i; ++i)
	{
		X[i] = 1;
	}

	// ca and cb
	double *psum_y = REAL(sum_y_vector);

	vector<double> ca(n_i);
	vector<double> cb(n_j);
	for (int i = 0; i < n_i; ++i)
	{
		ca[i] = psum_y[i];
	}

	double *psum_y_j = psum_y + n_i;
	for (int j = 0; j < n_j; ++j)
	{
		cb[j] = psum_y_j[j];
	}

	bool keepGoing = true;
	int iter = 0;
	while (keepGoing && iter < iterMax)
	{
		// Rprintf("%i.", iter);
		++iter;

		// This way I don't need to update the values of X with a loop
		if (iter % 2 == 1)
		{
			CCC_poisson_2(X, X_new, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);
		}
		else
		{
			CCC_poisson_2(X_new, X, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);
		}

		// double *X_current = (iter % 2 == 1 ? X_new : X);
		// Rprintf("alpha: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f -- ", X_current[i]);
		// }
		// Rprintf("\n");

		keepGoing = false;
		for (int i = 0; i < n_i; ++i)
		{
			// if(fabs(X[i] - X_new[i]) / (0.1 + fabs(X_new[i])) > diffMax){
			if (continue_criterion(X[i], X_new[i], diffMax))
			{
				keepGoing = true;
				break;
			}
		}
	}

	double *X_final = (iter % 2 == 1 ? X_new.data() : X.data());

	SEXP exp_mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(exp_mu);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;

	double *beta = X_final + n_i;
	for (int obs = 0; obs < n_obs; ++obs)
	{
		pmu[obs] = pexp_mu_in[obs] * X_final[dum_i[obs]] * beta[dum_j[obs]];
	}

	// Rprintf(" beta: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", beta[i]);
	// }
	// Rprintf("\n");
	//
	// Rprintf("alpha: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", X_final[i]);
	// }
	// Rprintf("\n");

	UNPROTECT(1);

	writable::list res;
	res.push_back({"mu_new"_nm = exp_mu});
	res.push_back({"iter"_nm = iter});

	return (res);
}

[[cpp11::register]] list cpp_conv_acc_gau_2_(int n_i, int n_j, int n_cells,
											SEXP r_mat_row, SEXP r_mat_col, SEXP r_mat_value_Ab, SEXP r_mat_value_Ba,
											SEXP dum_vector, SEXP lhs, SEXP invTableCluster_vector,
											int iterMax, double diffMax, SEXP mu_in)
{

	//
	// Setting up
	//

	int n_obs = Rf_length(mu_in);

	int *mat_row = INTEGER(r_mat_row);
	int *mat_col = INTEGER(r_mat_col);
	double *mat_value_Ab = REAL(r_mat_value_Ab);
	double *mat_value_Ba = REAL(r_mat_value_Ba);

	vector<double> resid(n_obs);
	double *plhs = REAL(lhs), *pmu_in = REAL(mu_in);
	for (int obs = 0; obs < n_obs; ++obs)
	{
		resid[obs] = plhs[obs] - pmu_in[obs];
	}

	//
	// const_a and const_b
	vector<double> const_a(n_i, 0);
	vector<double> const_b(n_j, 0);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;
	double *invTable_i = REAL(invTableCluster_vector);
	double *invTable_j = invTable_i + n_i;

	// for(int i=0 ; i<n_i ; ++i){
	// 	// init
	// 	const_a[i] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	const_a[dum_i[obs]] += resid[obs] * invTable_i[dum_i[obs]];
	// }
	//
	// for(int j=0 ; j<n_j ; ++j){
	// 	// init
	// 	const_b[j] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	const_b[dum_j[obs]] += resid[obs] * invTable_j[dum_j[obs]];
	// }

	for (int obs = 0; obs < n_obs; ++obs)
	{
		double resid_tmp = resid[obs];
		int d_i = dum_i[obs];
		int d_j = dum_j[obs];

		const_a[d_i] += resid_tmp * invTable_i[d_i];

		const_b[d_j] += resid_tmp * invTable_j[d_j];
	}

	// Rprintf("const_a: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", const_a[i]);
	// }
	// Rprintf("\n");
	//
	//
	// Rprintf("const_b: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", const_b[i]);
	// }
	// Rprintf("\n");

	// values that will be used later
	vector<double> beta(n_j);

	//
	// some tests
	//

	// double test_Ab[n_i], test_Ba[n_j];
	// for(int i=0 ; i<n_i ; ++i){
	// 	test_Ab[i] = 0;
	// }
	// for(int j=0 ; j<n_j ; ++j){
	// 	test_Ba[j] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_cells ; ++obs){
	// 	test_Ab[mat_row[obs]] += mat_value_Ab[obs];
	// 	test_Ba[mat_col[obs]] += mat_value_Ba[obs];
	// }

	// Rprintf("Ab x 1: \n");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f \t", test_Ab[i]);
	// }
	// Rprintf("\n");
	//
	// Rprintf("Ba x 1: \n");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f \t", test_Ba[i]);
	// }
	// Rprintf("\n");

	//
	// alpha_tilde:
	// a_tilde = const_a - (Ab %m% const_b)
	// vector<double> a_tilde(n_i);
	//
	// for(int i=0 ; i<n_i ; ++i){
	// 	a_tilde[i] = const_a[i];
	// }

	vector<double> a_tilde(const_a); // init at const_a

	for (int obs = 0; obs < n_cells; ++obs)
	{
		a_tilde[mat_row[obs]] -= mat_value_Ab[obs] * const_b[mat_col[obs]];
	}

	//
	// IT iteration (preparation)
	//

	vector<double> X(n_i);
	vector<double> GX(n_i);
	vector<double> GGX(n_i);
	vector<double> delta_GX(n_i);
	vector<double> delta2_X(n_i);

	//
	// the main loop
	//

	// initialisation of X and then GX
	for (int i = 0; i < n_i; ++i)
	{
		// X[i] = a_tilde[i];
		X[i] = 0;
	}

	// first iteration
	CCC_gaussian_2(X, GX, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);

	// Rprintf("  X: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f\t", X[i]);
	// }
	// Rprintf("\n");
	//
	//
	// Rprintf(" GX: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f\t", GX[i]);
	// }
	// Rprintf("\n");

	// stop("kkk");

	bool numconv = false;
	bool keepGoing = true;
	int iter = 0;
	while (keepGoing && iter < iterMax)
	{
		// Rprintf("%i.", iter);
		++iter;

		// GGX -- origin: GX, destination: GGX
		CCC_gaussian_2(GX, GGX, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);

		// Rprintf("ggx: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f \t", GGX[i]);
		// }
		// Rprintf("\n");

		// X ; update of the cluster coefficient
		numconv = update_X_IronsTuck(n_i, X, GX, GGX, delta_GX, delta2_X);
		if (numconv)
			break;

		// Rprintf("  x: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f \t", X[i]);
		// }
		// Rprintf("\n");

		// GX -- origin: X, destination: GX
		CCC_gaussian_2(X, GX, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);

		// Rprintf("   gx: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f\t", GX[i]);
		// }
		// Rprintf("\n");

		keepGoing = false;
		for (int i = 0; i < n_i; ++i)
		{
			// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
			if (continue_criterion(X[i], GX[i], diffMax))
			{
				keepGoing = true;
				break;
			}
		}
	}

	SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(mu);

	// we need to compute beta, and then alpha

	//	beta = const_b - (Ba %m% alpha)
	// vector<double> beta_final(n_j);
	// for(int j=0 ; j<n_j ; ++j){
	// 	beta_final[j] = const_b[j];
	// }
	vector<double> beta_final(const_b);

	for (int obs = 0; obs < n_cells; ++obs)
	{
		beta_final[mat_col[obs]] -= mat_value_Ba[obs] * GX[mat_row[obs]];
	}

	// alpha = const_a - (Ab %m% beta)
	// vector<double> alpha_final(n_i);
	// for(int i=0 ; i<n_i ; ++i){
	// 	alpha_final[i] = const_a[i];
	// }
	vector<double> alpha_final(const_a);

	for (int obs = 0; obs < n_cells; ++obs)
	{
		alpha_final[mat_row[obs]] -= mat_value_Ab[obs] * beta_final[mat_col[obs]];
	}

	// mu final
	for (int obs = 0; obs < n_obs; ++obs)
	{
		pmu[obs] = pmu_in[obs] + alpha_final[dum_i[obs]] + beta_final[dum_j[obs]];
	}

	// Rprintf("alpha: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f\t", alpha_final[i]);
	// }
	// Rprintf("\n");

	UNPROTECT(1);

	writable::list res;
	res.push_back({"mu_new"_nm = mu});
	res.push_back({"iter"_nm = iter});

	return (res);
}

[[cpp11::register]] list cpp_conv_seq_gau_2_(int n_i, int n_j, int n_cells,
											SEXP r_mat_row, SEXP r_mat_col, SEXP r_mat_value_Ab, SEXP r_mat_value_Ba,
											SEXP dum_vector, SEXP lhs, SEXP invTableCluster_vector,
											int iterMax, double diffMax, SEXP mu_in)
{

	//
	// Setting up
	//

	int n_obs = Rf_length(mu_in);

	int *mat_row = INTEGER(r_mat_row);
	int *mat_col = INTEGER(r_mat_col);
	double *mat_value_Ab = REAL(r_mat_value_Ab);
	double *mat_value_Ba = REAL(r_mat_value_Ba);

	vector<double> resid(n_obs);
	double *plhs = REAL(lhs), *pmu_in = REAL(mu_in);
	for (int obs = 0; obs < n_obs; ++obs)
	{
		resid[obs] = plhs[obs] - pmu_in[obs];
	}

	//
	// const_a and const_b
	vector<double> const_a(n_i, 0);
	vector<double> const_b(n_j, 0);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;
	double *invTable_i = REAL(invTableCluster_vector);
	double *invTable_j = invTable_i + n_i;

	// for(int i=0 ; i<n_i ; ++i){
	// 	// init
	// 	const_a[i] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	const_a[dum_i[obs]] += resid[obs] * invTable_i[dum_i[obs]];
	// }
	//
	// for(int j=0 ; j<n_j ; ++j){
	// 	// init
	// 	const_b[j] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	const_b[dum_j[obs]] += resid[obs] * invTable_j[dum_j[obs]];
	// }

	for (int obs = 0; obs < n_obs; ++obs)
	{
		double resid_tmp = resid[obs];
		int d_i = dum_i[obs];
		int d_j = dum_j[obs];

		const_a[d_i] += resid_tmp * invTable_i[d_i];

		const_b[d_j] += resid_tmp * invTable_j[d_j];
	}

	// values that will be used later
	vector<double> beta(n_j);

	//
	// alpha_tilde:
	// a_tilde = const_a - (Ab %m% const_b)
	// vector<double> a_tilde(n_i);
	//
	// for(int i=0 ; i<n_i ; ++i){
	// 	a_tilde[i] = const_a[i];
	// }

	vector<double> a_tilde(const_a);

	for (int obs = 0; obs < n_cells; ++obs)
	{
		a_tilde[mat_row[obs]] -= mat_value_Ab[obs] * const_b[mat_col[obs]];
	}

	//
	// the main loop
	//

	// X, X_new => the vector of coefficients
	vector<double> X(n_i);
	vector<double> X_new(n_i);

	// initialisation of X
	for (int i = 0; i < n_i; ++i)
	{
		X[i] = a_tilde[i];
	}

	bool keepGoing = true;
	int iter = 0;
	while (keepGoing && iter < iterMax)
	{
		++iter;

		// This way I don't need to update the values of X with a loop
		if (iter % 2 == 1)
		{
			CCC_gaussian_2(X, X_new, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);
		}
		else
		{
			CCC_gaussian_2(X_new, X, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);
		}

		keepGoing = false;
		for (int i = 0; i < n_i; ++i)
		{
			// if(fabs(X[i] - X_new[i]) / (0.1 + fabs(X_new[i])) > diffMax){
			if (continue_criterion(X[i], X_new[i], diffMax))
			{
				keepGoing = true;
				break;
			}
		}
	}

	double *X_final = (iter % 2 == 1 ? X_new.data() : X.data());

	SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(mu);

	// we need to compute beta, and then alpha

	//	beta = const_b - (Ba %m% alpha)
	// vector<double> beta_final(n_j);
	// for(int j=0 ; j<n_j ; ++j){
	// 	beta_final[j] = const_b[j];
	// }
	vector<double> beta_final(const_b);

	for (int obs = 0; obs < n_cells; ++obs)
	{
		beta_final[mat_col[obs]] -= mat_value_Ba[obs] * X_final[mat_row[obs]];
	}

	// alpha = const_a - (Ab %m% beta)
	// vector<double> alpha_final(n_i);
	// for(int i=0 ; i<n_i ; ++i){
	// 	alpha_final[i] = const_a[i];
	// }
	vector<double> alpha_final(const_a);

	for (int obs = 0; obs < n_cells; ++obs)
	{
		alpha_final[mat_row[obs]] -= mat_value_Ab[obs] * beta_final[mat_col[obs]];
	}

	// mu final
	for (int obs = 0; obs < n_obs; ++obs)
	{
		pmu[obs] = pmu_in[obs] + alpha_final[dum_i[obs]] + beta_final[dum_j[obs]];
	}

	UNPROTECT(1);

	writable::list res;
	res.push_back({"mu_new"_nm = mu});
	res.push_back({"iter"_nm = iter});

	return (res);
}
