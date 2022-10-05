// Now we start a big chunk => computing the varying slopes coefficients
// That's a big job. To simplify it, I created the class FEClass that takes care of it.

#include "02_0_demeaning.h"

FEClass::FEClass(int n_obs, int Q, SEXP r_weights, SEXP fe_id_list, SEXP r_nb_id_Q, SEXP table_id_I, SEXP slope_flag_Q, SEXP slope_vars_list)
{
    // The constructor does the job of creating the pre-solved system of equations

    // Information on p_slope_flag_Q:
    // - vector of length Q
    // - say sf = p_slope_flag_Q[q]
    //  * abs(sf) number of VARIABLES with varying slope
    //  * sf > 0: we add the fixed-effect
    //  * sf < 0: the FE should NOT be included

    this->n_obs = n_obs;
    this->Q = Q;

    //
    // Step 0: General information
    //

    // fixed-effect id for each observation
    p_fe_id.resize(Q);
    // New version => dum_vector (a vector) is replaced by fe_id_list (a list)
    for (int q = 0; q < Q; ++q)
    {
        p_fe_id[q] = INTEGER(VECTOR_ELT(fe_id_list, q));
    }

    nb_id_Q = INTEGER(r_nb_id_Q);

    //
    // Step 1: we check if slopes are needed
    //

    // nb_slopes: number of variables with varying slopes (the FE does not count!)

    int nb_slopes = 0;
    int sf = 0; // slope flag
    int *p_slope_flag_Q = INTEGER(slope_flag_Q);
    vector<bool> is_slope_Q(Q, false);
    vector<bool> is_slope_fe_Q(Q, false);
    vector<int> nb_vs_Q(Q, 0);
    vector<int> nb_vs_noFE_Q(Q, 0);
    vector<int> nb_coef_Q(Q);
    int nb_coef_T = 0;

    for (int q = 0; q < Q; ++q)
    {
        //   0: no slope
        // < 0: slope but no fixed-effect
        // > 0: slope WITH fixed-effect
        // here we count the number of slopes only, we exclude the FEs (that's why there's the substraction)

        sf = p_slope_flag_Q[q];
        if (sf != 0)
        {
            nb_slopes += abs(sf);
            is_slope_Q[q] = true;
            nb_vs_Q[q] = abs(sf);
            nb_vs_noFE_Q[q] = abs(sf);

            if (sf > 0)
            {
                ++nb_vs_Q[q];
                is_slope_fe_Q[q] = true;
            }

            // There is n_vs coefficients (and not n_vs squared => this is only in the systems of eq)
            nb_coef_Q[q] = nb_vs_Q[q] * nb_id_Q[q];
        }
        else
        {
            nb_coef_Q[q] = nb_id_Q[q];
        }

        nb_coef_T += nb_coef_Q[q];
    }

    // where to start the coefficients
    vector<int> coef_start_Q(Q, 0);
    for (int q = 1; q < Q; ++q)
        coef_start_Q[q] = coef_start_Q[q - 1] + nb_coef_Q[q - 1];

    // Copying (tiny objects)
    this->is_slope_Q = is_slope_Q;
    this->is_slope_fe_Q = is_slope_fe_Q;
    this->nb_vs_Q = nb_vs_Q;
    this->nb_vs_noFE_Q = nb_vs_noFE_Q;
    this->nb_coef_Q = nb_coef_Q;
    this->nb_coef_T = nb_coef_T;
    this->coef_start_Q = coef_start_Q;

    //
    // Step 2: precomputed stuff for non slopes and slopes
    //

    // Weights
    bool is_weight = Rf_length(r_weights) != 1;
    this->is_weight = is_weight;
    p_weights = REAL(r_weights);

    is_slope = nb_slopes > 0;

    // First the non slope coefficients
    // We create the sum of weights

    int nb_coef_noVS_T = 0;
    for (int q = 0; q < Q; ++q)
    {
        if (is_slope_Q[q] == false)
        {
            nb_coef_noVS_T += nb_id_Q[q];
        }
    }

    sum_weights_noVS_C.resize(nb_coef_noVS_T > 0 ? nb_coef_noVS_T : 1);
    std::fill(sum_weights_noVS_C.begin(), sum_weights_noVS_C.end(), 0);

    p_sum_weights_noVS_C.resize(Q);
    p_sum_weights_noVS_C[0] = sum_weights_noVS_C.data();
    for (int q = 1; q < Q; ++q)
    {
        p_sum_weights_noVS_C[q] = p_sum_weights_noVS_C[q - 1] + (is_slope_Q[q - 1] == true ? 0 : nb_id_Q[q - 1]);
    }

    // table is already computed
    vector<int *> p_table_id_I(Q);
    p_table_id_I[0] = INTEGER(table_id_I);
    for (int q = 1; q < Q; ++q)
    {
        p_table_id_I[q] = p_table_id_I[q - 1] + nb_id_Q[q - 1];
    }

    for (int q = 0; q < Q; ++q)
    {
        if (is_slope_Q[q] == true)
        {
            continue;
        }

        double *my_SW = p_sum_weights_noVS_C[q];

        if (is_weight)
        {
            int *my_fe = p_fe_id[q];
            for (int obs = 0; obs < n_obs; ++obs)
            {
                my_SW[my_fe[obs] - 1] += p_weights[obs];
            }
        }
        else
        {
            int nb_coef = nb_id_Q[q];
            int *my_table = p_table_id_I[q];
            for (int i = 0; i < nb_coef; ++i)
            {
                my_SW[i] = my_table[i];
            }
        }
    }

    if (is_weight && nb_coef_noVS_T > 0)
    {
        // Checking the 0-weights => we set them to 1 to wavoid division by 0
        for (int c = 0; c < nb_coef_noVS_T; ++c)
        {
            if (sum_weights_noVS_C[c] == 0)
            {
                sum_weights_noVS_C[c] = 1;
            }
        }
    }

    // Then the slopes
    if (nb_slopes > 0)
    {

        // A) Meta variables => the ones containing the main information

        // slope_vars_list: R list
        // p_vs_vars.resize(nb_slopes);
        // for(int v=0 ; v<nb_slopes ; ++v){
        //     p_vs_vars[v] = REAL(VECTOR_ELT(slope_vars_list, v));
        // }
        sMat m_slopes(slope_vars_list);
        p_vs_vars.resize(nb_slopes);
        for (int v = 0; v < nb_slopes; ++v)
        {
            p_vs_vars[v] = m_slopes[v];
        }

        // B) Computing the coefficients of the systems of equations

        int nb_vs_coef_T = 0;
        for (int q = 0; q < Q; ++q)
        {
            nb_vs_coef_T += nb_vs_Q[q] * nb_vs_Q[q] * nb_id_Q[q];
        }

        // The sys of eqs, all coefs to 0
        eq_systems_VS_C.resize(nb_vs_coef_T);
        std::fill(eq_systems_VS_C.begin(), eq_systems_VS_C.end(), 0);

        p_eq_systems_VS_C.resize(Q);
        p_eq_systems_VS_C[0] = eq_systems_VS_C.data();
        for (int q = 1; q < Q; ++q)
        {
            p_eq_systems_VS_C[q] = p_eq_systems_VS_C[q - 1] + nb_vs_Q[q - 1] * nb_vs_Q[q - 1] * nb_id_Q[q - 1];
        }

        for (int q = 0; q < Q; ++q)
        {
            if (is_slope_Q[q] == false)
                continue;

            simple_mat_of_vs_vars VS_mat(this, q);
            simple_mat_with_id my_system(p_eq_systems_VS_C[q], nb_vs_Q[q]);

            int V = nb_vs_Q[q];
            int *my_fe = p_fe_id[q];
            int nb_coef = nb_id_Q[q];

            for (int i = 0; i < n_obs; ++i)
            {
                for (int v1 = 0; v1 < V; ++v1)
                {
                    for (int v2 = 0; v2 <= v1; ++v2)
                    {
                        if (is_weight)
                        {
                            my_system(my_fe[i] - 1, v1, v2) += VS_mat(i, v1) * VS_mat(i, v2) * p_weights[i];
                        }
                        else
                        {
                            my_system(my_fe[i] - 1, v1, v2) += VS_mat(i, v1) * VS_mat(i, v2);
                        }
                    }
                }
            }

            // Finishing the computation of the system (symmetry)
            for (int c = 0; c < nb_coef; ++c)
            {
                for (int v1 = 0; v1 < V; ++v1)
                {
                    for (int v2 = 0; v2 < v1; ++v2)
                    {
                        my_system(c, v2, v1) = my_system(c, v1, v2);
                    }
                }
            }

            // Precomputing the solver coefficients
            double my_row_coef = 0;
            for (int c = 0; c < nb_coef; ++c)
            {
                for (int v = 0; v < V; ++v)
                {
                    if (my_system(c, v, v) == 0)
                    {
                        // The pivot is equal to 0 => overidentified system
                        for (int i = v + 1; i < V; ++i)
                        {
                            my_system(c, i, v) = 0;
                        }
                    }
                    else
                    {
                        for (int i = v + 1; i < V; ++i)
                        {
                            my_row_coef = my_system(c, i, v) / my_system(c, v, v);
                            my_system(c, i, v) = my_row_coef;
                            for (int j = v + 1; j < V; ++j)
                            {
                                my_system(c, i, j) -= my_row_coef * my_system(c, v, j);
                            }
                        }
                    }
                }
            }

            // We end up with all the (pre-solved) systems of equations
        }
    }
}

// Overloaded versions
void FEClass::compute_fe_coef(double *fe_coef_C, sVec &mu_in_N)
{
    // mu: length n_obs, vector giving sum_in_out
    // fe_coef: vector receiving the cluster coefficients

    // single FE version
    compute_fe_coef_internal(0, fe_coef_C, true, mu_in_N, nullptr, nullptr);
}

void FEClass::compute_fe_coef(int q, double *fe_coef_C, double *sum_other_coef_N, double *in_out_C)
{
    // mu: length n_obs, vector giving sum_in_out
    // fe_coef: vector receiving the cluster coefficients

    // multiple FE version
    compute_fe_coef_internal(q, fe_coef_C, false, nullptr, sum_other_coef_N, in_out_C);
}

void FEClass::compute_fe_coef_internal(int q, double *fe_coef_C, bool is_single, sVec mu_in_N, double *sum_other_coef_N, double *in_out_C)
{
    // mu: length n_obs, vector giving sum_in_out
    // fe_coef: vector receiving the cluster coefficients

    int V = nb_vs_Q[q];
    int *my_fe = p_fe_id[q];
    int nb_coef = nb_coef_Q[q];
    int nb_id = nb_id_Q[q];

    double *my_fe_coef = fe_coef_C + coef_start_Q[q];

    if (is_slope_Q[q] == false)
    {

        double *my_SW = p_sum_weights_noVS_C[q];

        if (is_single)
        {

            for (int obs = 0; obs < n_obs; ++obs)
            {
                if (is_weight)
                {
                    my_fe_coef[my_fe[obs] - 1] += p_weights[obs] * mu_in_N[obs];
                }
                else
                {
                    my_fe_coef[my_fe[obs] - 1] += mu_in_N[obs];
                }
            }
        }
        else
        {

            double *sum_in_out = in_out_C + coef_start_Q[q];

            // initialize cluster coef
            for (int m = 0; m < nb_coef; ++m)
            {
                my_fe_coef[m] = sum_in_out[m];
            }

            // looping sequentially over the sum of other coefficients
            for (int i = 0; i < n_obs; ++i)
            {
                my_fe_coef[my_fe[i] - 1] -= sum_other_coef_N[i];
            }
        }

        // Rcout << "is_weight:" << is_weight << ", is_single: " << is_single << "\n";

        // Rcout << "Coefs:\n ";
        for (int m = 0; m < nb_coef; ++m)
        {
            // Rcout << my_fe_coef[m] << " ("  << my_SW[m] << "), ";

            my_fe_coef[m] /= my_SW[m];
        }
        // Rcout << "\n\n";
    }
    else
    {
        // Setting up => computing the raw coefficient of the last column
        simple_mat_of_vs_vars VS_mat(this, q);
        simple_mat_with_id my_system(p_eq_systems_VS_C[q], nb_vs_Q[q]);

        // We use the vector of FE coefficients
        simple_mat_with_id rhs(my_fe_coef, nb_vs_Q[q], 1);

        if (is_single)
        {

            for (int i = 0; i < n_obs; ++i)
            {
                for (int v1 = 0; v1 < V; ++v1)
                {
                    if (is_weight)
                    {
                        rhs(my_fe[i] - 1, v1) += VS_mat(i, v1) * mu_in_N[i] * p_weights[i];
                    }
                    else
                    {
                        rhs(my_fe[i] - 1, v1) += VS_mat(i, v1) * mu_in_N[i];
                    }
                }
            }
        }
        else
        {

            // initialization
            double *sum_in_out = in_out_C + coef_start_Q[q];
            for (int m = 0; m < nb_coef; ++m)
            {
                my_fe_coef[m] = sum_in_out[m];
            }

            for (int i = 0; i < n_obs; ++i)
            {
                for (int v = 0; v < V; ++v)
                {
                    rhs(my_fe[i] - 1, v) -= VS_mat(i, v) * sum_other_coef_N[i];
                }
            }
        }

        // "Correcting" the last column
        for (int c = 0; c < nb_id; ++c)
        {
            for (int v = 0; v < V; ++v)
            {
                for (int v1 = v + 1; v1 < V; ++v1)
                {
                    // -= row_coef * value from the v row
                    rhs(c, v1) -= my_system(c, v1, v) * rhs(c, v);
                }
            }
        }

        // Solving => we store the solution in rhs
        double val = 0;
        for (int c = 0; c < nb_id; ++c)
        {
            // We backward solve
            for (int v = V - 1; v >= 0; --v)
            {
                if (my_system(c, v, v) == 0)
                {
                    rhs(c, v) = 0;
                }
                else
                {
                    val = rhs(c, v);
                    for (int v_done = v + 1; v_done < V; ++v_done)
                    {
                        val -= rhs(c, v_done) * my_system(c, v, v_done);
                    }
                    rhs(c, v) = val / my_system(c, v, v);
                }
            }
        }
    }
}

void FEClass::compute_fe_coef_2_internal(double *fe_coef_in_out_C, double *fe_coef_tmp, double *in_out_C, bool step_2 = false)
{
    // The aim here is to update the value of b given the value of a
    // fe_coef_tmp is always of size the 2nd FE
    // fe_in/out are always the size the 1st FE
    // a: origin
    // b: destination

    int index_a = 0;
    int index_b = 1;

    double *my_fe_coef_a;
    double *my_fe_coef_b;

    if (step_2)
    {
        index_a = 1;
        index_b = 0;

        my_fe_coef_a = fe_coef_tmp;
        my_fe_coef_b = fe_coef_in_out_C;
    }
    else
    {

        my_fe_coef_a = fe_coef_in_out_C;
        my_fe_coef_b = fe_coef_tmp;
    }

    int V_a = nb_vs_Q[index_a];
    int V_b = nb_vs_Q[index_b];

    int *my_fe_a = p_fe_id[index_a];
    int *my_fe_b = p_fe_id[index_b];

    // int nb_coef_a = nb_coef_Q[index_a];
    int nb_coef_b = nb_coef_Q[index_b];

    // int nb_id_a = nb_id_Q[index_a];
    int nb_id_b = nb_id_Q[index_b];

    // double *my_in_out_a = in_out_C + coef_start_Q[index_a];
    double *my_in_out_b = in_out_C + coef_start_Q[index_b];

    bool is_slope_a = is_slope_Q[index_a];
    bool is_slope_b = is_slope_Q[index_b];

    // VSlope utilities
    simple_mat_of_vs_vars VS_mat_a(this, index_a);
    simple_mat_with_id my_vs_coef_a(my_fe_coef_a, V_a, 1);

    simple_mat_of_vs_vars VS_mat_b(this, index_b);
    simple_mat_with_id my_vs_coef_b(my_fe_coef_b, V_b, 1);

    // initialize cluster coef
    for (int m = 0; m < nb_coef_b; ++m)
    {
        my_fe_coef_b[m] = my_in_out_b[m];
    }

    if (is_slope_b == false)
    {

        // Adding the values of the first FE
        for (int i = 0; i < n_obs; ++i)
        {
            if (is_slope_a)
            {
                for (int v = 0; v < V_a; ++v)
                {
                    if (is_weight)
                    {
                        my_fe_coef_b[my_fe_b[i] - 1] -= my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v) * p_weights[i];
                    }
                    else
                    {
                        my_fe_coef_b[my_fe_b[i] - 1] -= my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v);
                    }
                }
            }
            else if (is_weight)
            {
                my_fe_coef_b[my_fe_b[i] - 1] -= my_fe_coef_a[my_fe_a[i] - 1] * p_weights[i];
            }
            else
            {
                my_fe_coef_b[my_fe_b[i] - 1] -= my_fe_coef_a[my_fe_a[i] - 1];
            }
        }

        double *my_SW = p_sum_weights_noVS_C[index_b];
        for (int m = 0; m < nb_coef_b; ++m)
        {
            my_fe_coef_b[m] /= my_SW[m];
        }
    }
    else
    {
        simple_mat_with_id my_system_b(p_eq_systems_VS_C[index_b], V_b);
        simple_mat_with_id rhs_b(my_fe_coef_b, V_b, 1);

        double tmp = 0;

        for (int i = 0; i < n_obs; ++i)
        {

            // if(is_slope_a){
            //     tmp = 0;
            //     for(int v=0 ; v<V_a ; ++v){
            //         if(is_weight){
            //             tmp += my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v) * p_weights[i];
            //         } else {
            //             tmp += my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v);
            //         }
            //     }
            //
            // } else if(is_weight){
            //     tmp = my_fe_coef_a[my_fe_a[i] - 1] * p_weights[i];
            // } else {
            //     tmp = my_fe_coef_a[my_fe_a[i] - 1];
            // }
            if (is_slope_a)
            {
                tmp = 0;
                for (int v = 0; v < V_a; ++v)
                {
                    tmp += my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v);
                }
            }
            else
            {
                tmp = my_fe_coef_a[my_fe_a[i] - 1];
            }

            for (int v = 0; v < V_b; ++v)
            {
                if (is_weight)
                {
                    rhs_b(my_fe_b[i] - 1, v) -= VS_mat_b(i, v) * tmp * p_weights[i];
                }
                else
                {
                    rhs_b(my_fe_b[i] - 1, v) -= VS_mat_b(i, v) * tmp;
                }
            }
        }

        // "Correcting" the last column
        for (int c = 0; c < nb_id_b; ++c)
        {
            for (int v = 0; v < V_b; ++v)
            {
                for (int v1 = v + 1; v1 < V_b; ++v1)
                {
                    // -= row_coef * value from the v row
                    rhs_b(c, v1) -= my_system_b(c, v1, v) * rhs_b(c, v);
                }
            }
        }

        // Solving => we store the solution in rhs
        double val = 0;
        for (int c = 0; c < nb_id_b; ++c)
        {
            // We backward solve
            for (int v = V_b - 1; v >= 0; --v)
            {
                if (my_system_b(c, v, v) == 0)
                {
                    rhs_b(c, v) = 0;
                }
                else
                {
                    val = rhs_b(c, v);
                    for (int v_done = v + 1; v_done < V_b; ++v_done)
                    {
                        val -= rhs_b(c, v_done) * my_system_b(c, v, v_done);
                    }
                    rhs_b(c, v) = val / my_system_b(c, v, v);
                }
            }
        }
    }
}

void FEClass::compute_fe_coef_2(double *fe_coef_in_C, double *fe_coef_out_C, double *fe_coef_tmp, double *in_out_C)
{
    // Specific to the 2-FEs case
    // This way we avoid creating and using a temp object of length N

    //
    // Step 1: Updating b
    //

    compute_fe_coef_2_internal(fe_coef_in_C, fe_coef_tmp, in_out_C);

    // Rcout << "Coefs IN:\nFirst dim: ";
    // for(int i=0 ; i<nb_coef_Q[0] ; ++i){
    //     Rcout << fe_coef_in_C[i] << ", ";
    // }
    // Rcout << "\n";

    //
    // Step 2: Updating a
    //

    compute_fe_coef_2_internal(fe_coef_out_C, fe_coef_tmp, in_out_C, true);

    // Rcout << "Coefs OUT:\nFirst dim: ";
    // for(int i=0 ; i<nb_coef_Q[0] ; ++i){
    //     Rcout << fe_coef_out_C[i] << ", ";
    // }
    // Rcout << "\n";
}

void FEClass::add_wfe_coef_to_mu(int q, double *fe_coef_C, double *out_N)
{
    // We add the weighted FE coefficient to each observation

    add_wfe_coef_to_mu_internal(q, fe_coef_C, out_N, true);
}

void FEClass::add_fe_coef_to_mu(int q, double *fe_coef_C, double *out_N)
{
    // We add the FE coefficient to each observation -- NO WEIGHTS!!!

    // Single FE version
    add_wfe_coef_to_mu_internal(q, fe_coef_C, out_N, false);
}

void FEClass::add_wfe_coef_to_mu_internal(int q, double *fe_coef_C, double *out_N, bool add_weights)
{
    // We just add the weighted (or not) FE coefficient to each observation
    // We need to be careful for VS

    int V = nb_vs_Q[q];
    int *my_fe = p_fe_id[q];

    double *my_fe_coef = fe_coef_C + coef_start_Q[q];

    bool use_weights = add_weights && is_weight;

    if (is_slope_Q[q] == false)
    {

        for (int i = 0; i < n_obs; ++i)
        {
            if (use_weights)
            {
                out_N[i] += my_fe_coef[my_fe[i] - 1] * p_weights[i];
            }
            else
            {
                out_N[i] += my_fe_coef[my_fe[i] - 1];
            }
        }
    }
    else
    {
        simple_mat_of_vs_vars VS_mat(this, q);
        simple_mat_with_id my_vs_coef(my_fe_coef, nb_vs_Q[q], 1);

        for (int i = 0; i < n_obs; ++i)
        {
            for (int v = 0; v < V; ++v)
            {
                if (use_weights)
                {
                    out_N[i] += my_vs_coef(my_fe[i] - 1, v) * VS_mat(i, v) * p_weights[i];
                }
                else
                {
                    out_N[i] += my_vs_coef(my_fe[i] - 1, v) * VS_mat(i, v);
                }
            }
        }
    }
}

void FEClass::compute_in_out(int q, double *in_out_C, sVec &in_N, double *out_N)
{
    // output: vector of length the number of coefficients

    int V = nb_vs_Q[q];
    int *my_fe = p_fe_id[q];

    double *sum_in_out = in_out_C + coef_start_Q[q];

    if (is_slope_Q[q] == false)
    {

        for (int i = 0; i < n_obs; ++i)
        {
            if (is_weight)
            {
                sum_in_out[my_fe[i] - 1] += (in_N[i] - out_N[i]) * p_weights[i];
            }
            else
            {
                sum_in_out[my_fe[i] - 1] += (in_N[i] - out_N[i]);
            }
        }
    }
    else
    {
        simple_mat_of_vs_vars VS_mat(this, q);
        simple_mat_with_id my_vs_sum_in_out(sum_in_out, nb_vs_Q[q], 1);

        for (int i = 0; i < n_obs; ++i)
        {
            for (int v = 0; v < V; ++v)
            {
                if (is_weight)
                {
                    my_vs_sum_in_out(my_fe[i] - 1, v) += (in_N[i] - out_N[i]) * VS_mat(i, v) * p_weights[i];
                }
                else
                {
                    my_vs_sum_in_out(my_fe[i] - 1, v) += (in_N[i] - out_N[i]) * VS_mat(i, v);
                }
            }
        }
    }
}

FEClass::simple_mat_of_vs_vars::simple_mat_of_vs_vars(const FEClass *FE_info, int q)
{
    // We set up the matrix
    int start = 0;
    for (int l = 0; l < q; ++l)
        start += FE_info->nb_vs_noFE_Q[l];

    int K = FE_info->nb_vs_noFE_Q[q];
    pvars.resize(K);
    for (int k = 0; k < K; ++k)
    {
        pvars[k] = FE_info->p_vs_vars[start + k];
    }

    if (FE_info->is_slope_fe_Q[q])
    {
        K_fe = K;
    }
    else
    {
        K_fe = -1;
    }
}

double FEClass::simple_mat_of_vs_vars::operator()(int i, int k)
{
    if (k == K_fe)
    {
        return 1;
    }

    return pvars[k][i];
}

void compute_fe_gnl(double *p_fe_coef_origin, double *p_fe_coef_destination,
                    double *p_sum_other_means, double *p_sum_in_out, PARAM_DEMEAN *args)
{
    // update of the cluster coefficients
    // first we update mu, then we update the cluster coefficicents

    //
    // Loading the variables
    //

    int n_obs = args->n_obs;
    int Q = args->Q;

    // fe_info
    FEClass &FE_info = *(args->p_FE_info);

    // We update each cluster coefficient, starting from Q (the smallest one)

    std::fill_n(p_sum_other_means, n_obs, 0);

    // we start with Q-1
    for (int q = 0; q < (Q - 1); ++q)
    {
        FE_info.add_wfe_coef_to_mu(q, p_fe_coef_origin, p_sum_other_means);
    }

    // Rcout << "Head of sum_other_means: ";
    // for(int i=0 ; i<10 ; ++i) Rcout << p_sum_other_means[i] << ", ";
    // Rcout << "\n";

    for (int q = Q - 1; q >= 0; q--)
    {

        FE_info.compute_fe_coef(q, p_fe_coef_destination, p_sum_other_means, p_sum_in_out);

        // if(int q == 0){
        // 	Rprintf("p_fe_coef_destination: %.3f, %.3f, %.3f, %.3f\n", my_fe_coef[0], my_fe_coef[1], my_fe_coef[2], my_fe_coef[3]);
        // }

        // updating the value of p_sum_other_means (only if necessary)
        if (q != 0)
        {

            // We recompute it from scratch (only way -- otherwise precision problems arise)

            std::fill_n(p_sum_other_means, n_obs, 0);

            double *my_fe_coef;
            for (int h = 0; h < Q; h++)
            {
                if (h == q - 1)
                    continue;

                if (h < q - 1)
                {
                    my_fe_coef = p_fe_coef_origin;
                }
                else
                {
                    my_fe_coef = p_fe_coef_destination;
                }

                FE_info.add_wfe_coef_to_mu(h, my_fe_coef, p_sum_other_means);
            }
        }
    }

    // In the end, the array p_fe_coef_destination is fully updated, starting from Q to 1
}
