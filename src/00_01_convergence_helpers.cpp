// Stopping / continuing criteria
// Functions used inside all loops

#include "00_convergence.h"

inline bool continue_criterion(double a, double b, double diffMax)
{
	// continuing criterion of the algorithm
	double diff = fabs(a - b);
	return ((diff > diffMax) && (diff / (0.1 + fabs(a)) > diffMax));
}

inline bool stopping_criterion(double a, double b, double diffMax)
{
	// stopping criterion of the algorithm
	double diff = fabs(a - b);
	return ((diff < diffMax) || (diff / (0.1 + fabs(a)) < diffMax));
}

// IT update + returns numerical convergence indicator
bool update_X_IronsTuck(int nb_coef_no_K, vector<double> &X,
						const vector<double> &GX, const vector<double> &GGX,
						vector<double> &delta_GX, vector<double> &delta2_X)
{

	for (int i = 0; i < nb_coef_no_K; ++i)
	{
		double GX_tmp = GX[i];
		delta_GX[i] = GGX[i] - GX_tmp;
		delta2_X[i] = delta_GX[i] - GX_tmp + X[i];
		// delta_GX[i] = GGX[i] - GX[i];
		// delta2_X[i] = delta_GX[i] - GX[i] + X[i];
	}

	// delta_GX %*% delta2_X and crossprod(delta2_X)
	double vprod = 0, ssq = 0;
	for (int i = 0; i < nb_coef_no_K; ++i)
	{
		double delta2_X_tmp = delta2_X[i];
		vprod += delta_GX[i] * delta2_X_tmp;
		ssq += delta2_X_tmp * delta2_X_tmp;
		// vprod += delta_GX[i] * delta2_X[i];
		// ssq += delta2_X[i] * delta2_X[i];
	}

	bool res = false;

	if (ssq == 0)
	{
		res = true;
	}
	else
	{
		double coef = vprod / ssq;

		// update of X:
		for (int i = 0; i < nb_coef_no_K; ++i)
		{
			X[i] = GGX[i] - coef * delta_GX[i];
		}
	}

	return (res);
}

[[cpp11::register]] list cpp_fixed_cost_gaussian(int n_i, int n_cells, SEXP index_i, SEXP index_j, SEXP order,
												 SEXP invTableCluster_vector, SEXP dum_vector)
{

	// conversion of R objects
	int n_obs = Rf_length(index_i);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;
	double *invTable_i = REAL(invTableCluster_vector);
	double *invTable_j = invTable_i + n_i;

	// We compute the matrix Ab and Ba
	int index_current = 0;

	// the objects returned
	SEXP r_mat_row = PROTECT(Rf_allocVector(INTSXP, n_cells));
	SEXP r_mat_col = PROTECT(Rf_allocVector(INTSXP, n_cells));
	SEXP r_mat_value_Ab = PROTECT(Rf_allocVector(REALSXP, n_cells));
	SEXP r_mat_value_Ba = PROTECT(Rf_allocVector(REALSXP, n_cells));
	int *mat_row = INTEGER(r_mat_row);
	int *mat_col = INTEGER(r_mat_col);
	double *mat_value_Ab = REAL(r_mat_value_Ab);
	double *mat_value_Ba = REAL(r_mat_value_Ba);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);
	int *new_order = INTEGER(order);

	// double value_Ab = invTable_i[dum_i[0]];
	// double value_Ba = invTable_j[dum_j[0]];
	double value_Ab = invTable_i[dum_i[new_order[0]]];
	double value_Ba = invTable_j[dum_j[new_order[0]]];

	for (int obs = 1; obs < n_obs; ++obs)
	{
		if (pindex_j[obs] != pindex_j[obs - 1] || pindex_i[obs] != pindex_i[obs - 1])
		{
			// save the value if we change index
			mat_row[index_current] = pindex_i[obs - 1];
			mat_col[index_current] = pindex_j[obs - 1];
			mat_value_Ab[index_current] = value_Ab;
			mat_value_Ba[index_current] = value_Ba;

			// new row
			index_current++;

			// the new value
			value_Ab = invTable_i[dum_i[new_order[obs]]];
			value_Ba = invTable_j[dum_j[new_order[obs]]];
		}
		else
		{
			value_Ab += invTable_i[dum_i[new_order[obs]]];
			value_Ba += invTable_j[dum_j[new_order[obs]]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs - 1];
	mat_col[index_current] = pindex_j[n_obs - 1];
	mat_value_Ab[index_current] = value_Ab;
	mat_value_Ba[index_current] = value_Ba;

	// the object returned
	writable::list res;
	res["mat_row"] = r_mat_row;
	res["mat_col"] = r_mat_col;
	res["mat_value_Ab"] = r_mat_value_Ab;
	res["mat_value_Ba"] = r_mat_value_Ba;

	UNPROTECT(4);

	return (res);
}
