//
// Demeans each variable in input
// The method is based on obtaining the optimal cluster coefficients
//

// Works but only the master thread can call that function!
// What happens if the master thread has finished its job but the lower thread is in an "infinite" loop?
// this is tricky, as such we cannot stop it.
// Solution: create a function keeping the threads idle waiting for the complete job to be done
// BUT I need to add static allocation of threads => performance cost

// List of objects, used to
// lighten the writting of the functions

#include "02_demeaning.h"

void FEClass::add_2_fe_coef_to_mu(double *fe_coef_a, double *fe_coef_b, double *in_out_C, double *out_N, bool update_beta = true)
{
    // We add the value of the FE coefficients to each observation

    //
    // Step 1: we update the coefficients of b
    //

    if (update_beta)
    {
        compute_fe_coef_2_internal(fe_coef_a, fe_coef_b, in_out_C, out_N);
    }

    //
    // Step 2: we add the value of each coef
    //

    for (int q = 0; q < 2; ++q)
    {

        double *my_fe_coef = q == 0 ? fe_coef_a : fe_coef_b;
        int *my_fe = p_fe_id[q];
        bool is_slope = is_slope_Q[q];
        int V = nb_vs_Q[q];

        simple_mat_of_vs_vars VS_mat(this, q);
        simple_mat_with_id my_vs_coef(my_fe_coef, nb_vs_Q[q], 1);

        for (int i = 0; i < n_obs; ++i)
        {
            if (is_slope)
            {
                for (int v = 0; v < V; ++v)
                {
                    out_N[i] += my_vs_coef(my_fe[i] - 1, v) * VS_mat(i, v);
                }
            }
            else
            {
                out_N[i] += my_fe_coef[my_fe[i] - 1];
            }
        }
    }
}

void demean_single_1(int v, PARAM_DEMEAN *args)
{
    // v: variable identifier to demean

    // Q == 1: nothing to say, just compute the closed form

    // loading the data
    int nb_coef_T = args->nb_coef_T;

    vector<sVec> &p_input = args->p_input;
    vector<double *> &p_output = args->p_output;

    // fe_info
    FEClass &FE_info = *(args->p_FE_info);

    // vector of fixed-effects coefficients initialized at 0
    vector<double> fe_coef(nb_coef_T, 0);
    double *p_fe_coef = fe_coef.data();

    // interruption handling

    // TODO: pending_interrupt is not supported in cpp11

    // define boolean isMaster which is true if omp_get_thread_num() == 0
    // and false otherwise
    // bool isMaster = omp_get_thread_num() == 0;

    // bool *pStopNow = args->stopnow;

    // if(isMaster){
    //     if(pending_interrupt()){
    //         *pStopNow = true;
    //     }
    // }

    // the input & output
    sVec &input = p_input[v];
    double *output = p_output[v];

    // We compute the FEs
    FE_info.compute_fe_coef(p_fe_coef, input);

    // Output:
    FE_info.add_fe_coef_to_mu(0, p_fe_coef, output);

    // saving the fixef coefs
    double *fixef_values = args->fixef_values;
    if (args->save_fixef)
    {
        for (int m = 0; m < nb_coef_T; ++m)
        {
            fixef_values[m] = fe_coef[m];
        }
    }
}

void demean_acc_2(int v, int iterMax, PARAM_DEMEAN *args)
{

    //
    // Setting up
    //

    // loading data
    int n_obs = args->n_obs;
    double diffMax = args->diffMax;

    // input output
    vector<sVec> &p_input = args->p_input;
    vector<double *> &p_output = args->p_output;
    sVec &input = p_input[v];
    double *output = p_output[v];

    int *iterations_all = args->p_iterations_all;

    // FE info
    FEClass &FE_info = *(args->p_FE_info);

    // vector of FE coefs => needs not be equal to nb_coef_T (because Q can be higher than 2)
    int n_a = FE_info.nb_coef_Q[0];
    int n_b = FE_info.nb_coef_Q[1];
    int nb_coef_2 = n_a + n_b;

    // conditional sum of input minus output
    vector<double> sum_input_output(nb_coef_2, 0);
    double *p_sum_in_out = sum_input_output.data();

    for (int q = 0; q < 2; ++q)
    {
        FE_info.compute_in_out(q, p_sum_in_out, input, output);
    }

    // Rcout << "Value of in_out:\n";
    // for(int c=0 ; c<nb_coef_2 ; ++c){
    //     if(c == n_a) Rcout << "\n";
    //     Rcout << p_sum_in_out[c] << ", ";
    //
    // }
    // Rcout << "\n";

    // Rprintf("a_tilde: %.3f, %.3f, %.3f, %.3f\n", a_tilde[0], a_tilde[1], a_tilde[2], a_tilde[3]);

    // TODO: interruption handling
    // bool isMaster = omp_get_thread_num() == 0;
    // bool *pStopNow = args->stopnow;
    // double flop = 20.0 * static_cast<double>(n_obs); // rough estimate nber operation per iter

    // int iterSecond = ceil(2000000000 / flop / 5);
    // nber iter per 1/5 second

    //
    // IT iteration (preparation)
    //

    vector<double> X(n_a, 0);
    vector<double> GX(n_a);
    vector<double> GGX(n_a);

    double *p_X = X.data();
    double *p_GX = GX.data();
    double *p_GGX = GGX.data();

    vector<double> delta_GX(n_a);
    vector<double> delta2_X(n_a);

    vector<double> coef_beta(n_b);
    double *p_coef_beta = coef_beta.data();

    //
    // the main loop
    //

    // first iteration => update GX
    FE_info.compute_fe_coef_2(p_X, p_GX, p_coef_beta, p_sum_in_out);

    // Rcout << "Coefs at first iteration:\na: ";
    // for(int i=0 ; i<n_a ; ++i){
    //     Rcout << p_GX[i] << ", ";
    // }
    // Rcout << "\nb: ";
    // for(int i=0 ; i<n_b ; ++i){
    //     Rcout << p_coef_beta[i] << ", ";
    // }
    // Rcout << "\n";

    // For the stopping criterion on total addition
    // vector<double> mu_last(n_obs, 0);
    // double input_mean = 0;
    double ssr = 0;

    bool numconv = false;
    bool keepGoing = true;
    int iter = 1;

    // TODO: re-enable pStopNow
    // while(!*pStopNow && keepGoing && iter<=iterMax){
    while (keepGoing && iter <= iterMax)
    {

        // if(isMaster && iter % iterSecond == 0){
        //     if(pending_interrupt()){
        //         *pStopNow = true;
        //         break;
        //     }
        // }

        ++iter;

        // GGX -- origin: GX, destination: GGX
        FE_info.compute_fe_coef_2(p_GX, p_GGX, p_coef_beta, p_sum_in_out);

        // X ; update of the fixed-effects coefficients
        numconv = update_X_IronsTuck(n_a, X, GX, GGX, delta_GX, delta2_X);
        if (numconv)
            break;

        // GX -- origin: X, destination: GX
        FE_info.compute_fe_coef_2(p_X, p_GX, p_coef_beta, p_sum_in_out);

        keepGoing = false;
        for (int i = 0; i < n_a; ++i)
        {
            if (continue_criterion(X[i], GX[i], diffMax))
            {
                keepGoing = true;
                break;
            }
        }

        // Other stopping criterion: change to SSR very small
        if (iter % 50 == 0)
        {

            vector<double> mu_current(n_obs, 0);
            double *p_mu = mu_current.data();

            FE_info.add_2_fe_coef_to_mu(p_GX, p_coef_beta, p_sum_in_out, p_mu);

            // init ssr if iter == 50 / otherwise, comparison
            if (iter == 50)
            {
                ssr = 0;
                double resid;
                for (int i = 0; i < n_obs; ++i)
                {
                    resid = input[i] - mu_current[i];
                    ssr += resid * resid;
                }
            }
            else
            {
                double ssr_old = ssr;

                // we compute the new SSR
                ssr = 0;
                double resid;
                for (int i = 0; i < n_obs; ++i)
                {
                    resid = input[i] - mu_current[i];
                    ssr += resid * resid;
                }

                // if(isMaster) Rprintf("iter %i -- SSR = %.0f (diff = %.0f)\n", iter, ssr, ssr_old - ssr);

                if (stopping_criterion(ssr_old, ssr, diffMax))
                {
                    break;
                }
            }
        }
    }

    //
    // we update the result (output)
    //

    // we end with a last iteration
    double *p_alpha_final = p_GX;
    double *p_beta_final = p_coef_beta;

    FE_info.compute_fe_coef_2(p_alpha_final, p_alpha_final, p_beta_final, p_sum_in_out);

    // Rcout << "Output before:\n";
    // for(int i=0 ; i<6 ; ++i){
    //     Rcout << output[i] << ",";
    // }

    // we update the output
    FE_info.add_2_fe_coef_to_mu(p_alpha_final, p_beta_final, p_sum_in_out, output, false);

    // Rcout << "\nOutput after:\n";
    // for(int i=0 ; i<6 ; ++i){
    //     Rcout << output[i] << ",";
    // }
    // Rcout << "\n";
    //
    // Rcout << "Final coefs:\na: ";
    // for(int i=0 ; i<n_a ; ++i){
    //     Rcout << p_alpha_final[i] << ", ";
    // }
    // Rcout << "\nb: ";
    // for(int i=0 ; i<n_b ; ++i){
    //     Rcout << p_beta_final[i] << ", ";
    // }
    // Rcout << "\n\n";

    // keeping track of iterations
    iterations_all[v] += iter;

    // saving the fixef coefs
    double *fixef_values = args->fixef_values;
    if (args->save_fixef)
    {

        for (int i = 0; i < n_a; ++i)
        {
            fixef_values[i] += p_alpha_final[i];
        }

        for (int j = 0; j < n_b; ++j)
        {
            fixef_values[n_a + j] += p_beta_final[j];
        }
    }
}

bool demean_acc_gnl(int v, int iterMax, PARAM_DEMEAN *args)
{

    //
    // data
    //

    int n_obs = args->n_obs;
    int nb_coef_T = args->nb_coef_T;
    int Q = args->Q;
    double diffMax = args->diffMax;

    // fe_info
    FEClass &FE_info = *(args->p_FE_info);

    // input output
    vector<sVec> &p_input = args->p_input;
    vector<double *> &p_output = args->p_output;
    sVec &input = p_input[v];
    double *output = p_output[v];

    // temp var:
    vector<double> sum_other_means(n_obs);
    double *p_sum_other_means = sum_other_means.data();

    // conditional sum of input minus output
    vector<double> sum_input_output(nb_coef_T, 0);
    double *p_sum_in_out = sum_input_output.data();

    for (int q = 0; q < Q; ++q)
    {
        FE_info.compute_in_out(q, p_sum_in_out, input, output);
    }

    // Rcout << "Sum in out:\n";
    // int c_next = 0;
    // int q_current = 0;
    // for(int c=0 ; c<nb_coef_T ; ++c){
    //     if(c == c_next){
    //         if(c != 0) Rcout << "\n";
    //
    //         Rcout << "q = " << q_current << ": ";
    //         c_next += FE_info.nb_coef_Q[q_current];
    //         q_current++;
    //     }
    //
    //     Rcout << p_sum_in_out[c] << ", ";
    // }
    // Rcout << "\n";

    // TODO: interruption handling
    // bool isMaster = omp_get_thread_num() == 0;
    // bool *pStopNow = args->stopnow;
    // I overcast to remember the lesson
    // double flop = 4.0*(5 + 12*(Q-1) + 4*(Q-1)*(Q-1))*static_cast<double>(n_obs);
    // rough estimate nber operation per iter
    // int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second

    //
    // IT iteration (preparation)
    //

    // variables on 1:K
    vector<double> X(nb_coef_T, 0);
    vector<double> GX(nb_coef_T);
    vector<double> GGX(nb_coef_T);
    // pointers:
    double *p_X = X.data();
    double *p_GX = GX.data();
    double *p_GGX = GGX.data();

    // variables on 1:(Q-1)
    int nb_coef_no_Q = 0;
    for (int q = 0; q < (Q - 1); ++q)
    {
        nb_coef_no_Q += FE_info.nb_coef_Q[q];
    }
    vector<double> delta_GX(nb_coef_no_Q);
    vector<double> delta2_X(nb_coef_no_Q);

    //
    // the main loop
    //

    // first iteration
    compute_fe_gnl(p_X, p_GX, p_sum_other_means, p_sum_in_out, args);

    // Rcout << "Coef first iter:\n";
    // c_next = 0;
    // q_current = 0;
    // for(int c=0 ; c<nb_coef_T ; ++c){
    //     if(c == c_next){
    //         if(c != 0) Rcout << "\n";
    //
    //         Rcout << "q = " << q_current << ": ";
    //         c_next += FE_info.nb_coef_Q[q_current];
    //         q_current++;
    //     }
    //
    //     Rcout << p_GX[c] << ", ";
    // }
    // Rcout << "\n\n\n";

    // check whether we should go into the loop
    bool keepGoing = false;
    for (int i = 0; i < nb_coef_T; ++i)
    {
        if (continue_criterion(X[i], GX[i], diffMax))
        {
            keepGoing = true;
            break;
        }
    }

    // For the stopping criterion on total addition
    // vector<double> mu_last(n_obs, 0);
    // double input_mean = 0;
    double ssr = 0;

    int iter = 0;
    bool numconv = false;
    // TODO: re-enable pstopnow
    // while(!*pStopNow && keepGoing && iter<iterMax){
    while (keepGoing && iter < iterMax)
    {

        // if(isMaster && iter % iterSecond == 0){
        //     if(pending_interrupt()){
        //         *pStopNow = true;
        //         break;
        //     }
        // }

        iter++;

        // GGX -- origin: GX, destination: GGX
        compute_fe_gnl(p_GX, p_GGX, p_sum_other_means, p_sum_in_out, args);

        // X ; update of the cluster coefficient
        numconv = update_X_IronsTuck(nb_coef_no_Q, X, GX, GGX, delta_GX, delta2_X);
        if (numconv)
            break;

        // GX -- origin: X, destination: GX
        compute_fe_gnl(p_X, p_GX, p_sum_other_means, p_sum_in_out, args);

        keepGoing = false;
        for (int i = 0; i < nb_coef_no_Q; ++i)
        {
            if (continue_criterion(X[i], GX[i], diffMax))
            {
                keepGoing = true;
                break;
            }
        }

        // Other stopping criterion: change to SSR very small
        if (iter % 50 == 0)
        {

            // mu_current is the vector of means
            vector<double> mu_current(n_obs, 0);
            double *p_mu = mu_current.data();
            for (int q = 0; q < Q; ++q)
            {

                FE_info.add_fe_coef_to_mu(q, p_GX, p_mu);
            }

            // init ssr if iter == 50 / otherwise, comparison
            if (iter == 50)
            {

                ssr = 0;
                double resid;
                for (int i = 0; i < n_obs; ++i)
                {
                    resid = input[i] - mu_current[i];
                    ssr += resid * resid;
                }
            }
            else
            {
                double ssr_old = ssr;

                // we compute the new SSR
                ssr = 0;
                double resid;
                for (int i = 0; i < n_obs; ++i)
                {
                    resid = input[i] - mu_current[i];
                    ssr += resid * resid;
                }

                // if(isMaster) Rprintf("iter %i -- SSR = %.0f (diff = %.0f)\n", iter, ssr, ssr_old - ssr);

                if (stopping_criterion(ssr_old, ssr, diffMax))
                {
                    break;
                }
            }
        }
    }

    //
    // Updating the output
    //

    for (int q = 0; q < Q; ++q)
    {
        FE_info.add_fe_coef_to_mu(q, p_GX, output);
    }

    // keeping track of iterations
    int *iterations_all = args->p_iterations_all;
    iterations_all[v] += iter;

    // saving the fixef coefs
    double *fixef_values = args->fixef_values;
    if (args->save_fixef)
    {
        for (int m = 0; m < nb_coef_T; ++m)
        {
            fixef_values[m] += GX[m];
        }
    }

    bool conv = iter == iterMax ? false : true;

    return (conv);
}

void demean_single_gnl(int v, PARAM_DEMEAN *args)
{
    // v: variable identifier to demean

    // Algorithm to quickly get the means
    // Q >= 3 => acceleration for 15 iter
    // if no convergence: conv 2 clusters
    // then acceleration again

    // data
    int iterMax = args->iterMax;
    int Q = args->Q;

    if (Q == 2)
    {
        demean_acc_2(v, iterMax, args);
        // demean_acc_gnl(v, iterMax, args);
    }
    else
    {
        // 15 iterations
        bool conv = demean_acc_gnl(v, 15, args);

        if (conv == false && iterMax > 15)
        {
            // 2 convergence
            demean_acc_2(v, iterMax / 2 - 15, args);

            if (Q > 2)
            {
                // re-acceleration
                demean_acc_gnl(v, iterMax / 2, args);
            }
        }
    }

    // if(args->p_iterations_all[v] > 50) balance_FEs(v, args);

    int *jobdone = args->jobdone;
    jobdone[v] = 1;
}

// Loop over demean_single
[[cpp11::register]] list cpp_demean(SEXP y, SEXP X_raw, SEXP r_weights, int iterMax, double diffMax, SEXP r_nb_id_Q,
                                    SEXP fe_id_list, SEXP table_id_I, SEXP slope_flag_Q, SEXP slope_vars_list,
                                    SEXP r_init, int nthreads, bool save_fixef = false)
{
    // main fun that calls demean_single
    // preformats all the information needed on the fixed-effects
    // y: the dependent variable
    // X_raw: the matrix of the explanatory variables -- can be "empty"

    // when including weights: recreate table values
    // export weights and is_weight bool

    // slope_flag: whether a FE is a varying slope
    // slope_var: the associated variables with varying slopes

    // initial variables
    int Q = Rf_length(r_nb_id_Q);

    // info on y
    sMat m_y(y);
    int n_vars_y = m_y.ncol();
    bool useY = n_vars_y > 0;
    bool is_y_list = n_vars_y > 1 || TYPEOF(y) == VECSXP;
    int n_obs = m_y.nrow();

    // info on X
    sMat m_X(X_raw);
    int n_vars_X = m_X.ncol();
    if (useY == false)
    {
        n_obs = m_X.nrow();
    }
    bool useX = n_vars_X > 0;

    int n_vars = n_vars_y + n_vars_X;

    // initialisation if needed (we never initialize when only one FE, except if asked explicitly)
    bool isInit = Rf_xlength(r_init) != 1 && Q > 1;
    double *init = REAL(r_init);
    bool saveInit = ((isInit || init[0] != 0) && Q > 1) || init[0] == 666;

    // Creating the object containing all information on the FEs
    FEClass FE_info(n_obs, Q, r_weights, fe_id_list, r_nb_id_Q, table_id_I, slope_flag_Q, slope_vars_list);
    int nb_coef_T = FE_info.nb_coef_T;

    // output vector: (Note that if the means are provided, we use that vector and will modify it in place)
    int64_t n_total = static_cast<int64_t>(n_obs) * n_vars;
    SEXP output_values = PROTECT(Rf_allocVector(REALSXP, isInit ? 1 : n_total));

    double *p_output_origin = isInit ? init : REAL(output_values);
    if (isInit == false)
    {
        std::fill_n(p_output_origin, n_total, 0);
    }

    //
    // vector of pointers: input/output
    //

    vector<double *> p_output(n_vars);
    p_output[0] = p_output_origin;
    for (int v = 1; v < n_vars; v++)
    {
        p_output[v] = p_output[v - 1] + n_obs;
    }

    vector<sVec> p_input(n_vars);

    for (int k = 0; k < n_vars_X; ++k)
    {
        p_input[k] = m_X[k];
    }

    for (int k = 0; k < n_vars_y; ++k)
    {
        p_input[n_vars_X + k] = m_y[k];
    }

    // keeping track of iterations
    vector<int> iterations_all(n_vars, 0);
    int *p_iterations_all = iterations_all.data();

    // save fixef option
    if (useX && save_fixef)
    {
        stop("save_fixef can be used only when there is no Xs.");
    }

    vector<double> fixef_values(save_fixef ? nb_coef_T : 1, 0);
    double *p_fixef_values = fixef_values.data();

    //
    // Sending variables to envir
    //

    PARAM_DEMEAN args;

    args.n_obs = n_obs;
    args.iterMax = iterMax;
    args.diffMax = diffMax;
    args.Q = Q;
    args.nb_coef_T = nb_coef_T;
    args.p_input = p_input;
    args.p_output = p_output;
    args.p_iterations_all = p_iterations_all;

    // save FE_info
    args.p_FE_info = &FE_info;

    // save fixef:
    args.save_fixef = save_fixef;
    args.fixef_values = p_fixef_values;

    // stopping flag + indicator that job is finished
    bool stopnow = false;
    args.stopnow = &stopnow;
    vector<int> jobdone(n_vars, 0);
    int *pjobdone = jobdone.data();
    args.jobdone = pjobdone;

    int counter = 0;
    int *pcounter = &counter;

    //
    // the main loop
    //

    int nthreads_current = nthreads > n_vars ? n_vars : nthreads;

// enlever les rprintf dans les nthreads jobs
#pragma omp parallel for num_threads(nthreads_current) schedule(static, 1)
    for (int v = 0; v < (n_vars + nthreads_current); ++v)
    {
        // demean_single is the workhorse
        // you get the "mean"

        if (!*(args.stopnow))
        {
            if (v < n_vars)
            {
                if (Q == 1)
                {
                    demean_single_1(v, &args);
                }
                else
                {
                    demean_single_gnl(v, &args);
                }
            }
            else if (true && Q != 1)
            {
                stayIdleCheckingInterrupt(&stopnow, jobdone, n_vars, pcounter);
            }
        }
    }

    if (*(args.stopnow))
    {
        UNPROTECT(1);
        stop("cpp_demean: User interrupt.");
    }

    // Rprintf("Master checking: %i\n", *pcounter);

    //
    // save
    //

    writable::list res; // a vector and a matrix

    int nrow = useX ? n_obs : 1;
    int ncol = useX ? n_vars_X : 1;
    writable::doubles_matrix<> X_demean(nrow, ncol);

    sVec p_input_tmp;
    double *p_output_tmp;
    for (int k = 0; k < n_vars_X; ++k)
    {
        p_input_tmp = p_input[k];
        p_output_tmp = p_output[k];

        for (int i = 0; i < n_obs; ++i)
        {
            X_demean(i, k) = p_input_tmp[i] - p_output_tmp[i];
        }
    }

    res["X_demean"] = X_demean;

    if (is_y_list && useY)
    {
        writable::list y_demean(n_vars_y);

        for (int v = 0; v < n_vars_y; ++v)
        {
            p_input_tmp = p_input[n_vars_X + v];
            p_output_tmp = p_output[n_vars_X + v];

            writable::doubles y_demean_tmp(n_obs);
            for (int i = 0; i < n_obs; ++i)
            {
                y_demean_tmp[i] = p_input_tmp[i] - p_output_tmp[i];
            }

            y_demean[v] = y_demean_tmp;
        }

        res.push_back({"y_demean"_nm = y_demean});
    }
    else
    {
        writable::doubles y_demean(useY ? n_obs : 1);
        if (useY)
        {
            // y is always the last variable
            p_input_tmp = p_input[n_vars - 1];
            p_output_tmp = p_output[n_vars - 1];
            for (int i = 0; i < n_obs; ++i)
            {
                y_demean[i] = p_input_tmp[i] - p_output_tmp[i];
            }
        }
        res.push_back({"y_demean"_nm = y_demean});
    }

    // iterations
    writable::integers iter_final(n_vars);
    for (int v = 0; v < n_vars; ++v)
    {
        iter_final[v] = p_iterations_all[v];
    }

    res.push_back({"iterations"_nm = iter_final});

    // if save is requested
    if (saveInit)
    {
        if (isInit)
        {
            res["means"] = r_init;
        }
        else
        {
            res["means"] = output_values;
        }
    }
    else
    {
        res.push_back({"means"_nm = 0.0});
    }

    // save fixef coef
    if (save_fixef)
    {
        res.push_back({"fixef_coef"_nm = fixef_values});
    }

    UNPROTECT(1);

    return (res);
}
