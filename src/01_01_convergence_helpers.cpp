#include "01_0_convergence.h"

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
	res.push_back({"mat_row"_nm = r_mat_row});
	res.push_back({"mat_col"_nm = r_mat_col});
	res.push_back({"mat_value_Ab"_nm = r_mat_value_Ab});
	res.push_back({"mat_value_Ba"_nm = r_mat_value_Ba});

	UNPROTECT(4);

	return (res);
}
