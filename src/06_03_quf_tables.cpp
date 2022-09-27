#include "06_0_quf.h"

void quf_table_sum_single(void *px_in, string &x_type, int n, int q, int *x_quf,
                          vector<double> &x_unik, vector<int> &x_table, double *py,
                          vector<double> &sum_y, bool do_sum_y, bool rm_0, bool rm_1,
                          bool rm_single, vector<int> &any_pblm, vector<bool> &id_pblm,
                          bool check_pblm, bool do_refactor, int x_size, integers &obs2keep)
{

    // check_pblm => FALSE only if only_slope = TRUE

    // Rcout << "q = " << q << ", n = " << n;

    // UFing
    if (do_refactor)
    {
        // refactoring of previous qufing => in one pass, we create table on the way
        int *px_in_int = (int *)px_in;

        // return;

        quf_refactor(px_in_int, x_size, obs2keep, n, x_quf, x_unik, x_table);

        if (obs2keep[0] != 0)
        {
            n = obs2keep.size();
        }
        // Rcout << "new n = " << n << "\n";
        // return;
    }
    else
    {
        quf_single(px_in, x_type, n, x_quf, x_unik);
    }

    // table + sum_y

    // You need sum_y to remove FEs (even if you don't need sum_y in the end)
    bool compute_sum_y = do_sum_y || rm_0;

    int D = x_unik.size();

    if (!do_refactor)
    {
        x_table.resize(D);
    }

    // Rcout << "D = " << D << "\n";
    // Rcout << "length(x_table) = " << x_table.size() << "\n";
    // Rcout << "compute_sum_y = " << compute_sum_y << "\n";
    // return;

    sum_y.resize(compute_sum_y > 0 ? D : 1);
    fill(sum_y.begin(), sum_y.end(), 0);

    // Rcout << "size sum_y = " << sum_y.size() << "\n";

    int obs;
    if (compute_sum_y || !do_refactor)
    {
        for (int i = 0; i < n; ++i)
        {
            obs = x_quf[i] - 1;
            if (!do_refactor)
                ++x_table[obs];
            if (compute_sum_y)
                sum_y[obs] += py[i];
        }
    }

    if ((rm_0 || rm_single) && check_pblm)
    {

        // 1) We check whether there is any problem
        int d_end = 0;
        for (int d = 0; d < D; ++d)
        {
            if ((rm_0 && sum_y[d] == 0) || (rm_1 && sum_y[d] == x_table[d]) || (rm_single && x_table[d] == 1))
            {
                any_pblm[q] = 1;
                d_end = d;
                break;
            }
        }

        // Rcout << "problem = " << any_pblm[q] << "\n";
        // return;

        // 2) If so: we find them
        if (any_pblm[q])
        {
            id_pblm.resize(D);
            fill(id_pblm.begin(), id_pblm.end(), false);

            for (int d = d_end; d < D; ++d)
            {
                // Rcout << "d = " << d << "\n";
                if ((rm_0 && sum_y[d] == 0) || (rm_1 && sum_y[d] == x_table[d]) || (rm_single && x_table[d] == 1))
                {
                    id_pblm[d] = true;
                }
            }

            // int sum_problems = accumulate(id_pblm.begin(), id_pblm.end(), 0);
            // Rcout << ", sum_problems = " << sum_problems;
        }
    }

    // Rcout << "\n";
}

[[cpp11::register]] writable::list cpppar_quf_table_sum(SEXP x, SEXP y, bool do_sum_y, bool rm_0, bool rm_1,
                          bool rm_single, writable::integers only_slope, int nthreads,
                          bool do_refactor, SEXP r_x_sizes, writable::integers obs2keep){

    // x: List of vectors of IDs (type int/num or char only)
    // y: dependent variable
    // rm_0: remove FEs where dep var is only 0
    // rm_1: remove FEs where dep var is only 0 or 1
    // rm_single: remove FEs with only one observation
    // do_sum_y: should we compute the sum_y?

    // When the data is refactored (ie x is a fixef_id_list, each element ranging from 1 to n_items):
    // - do_refactor
    // - r_x_sizes => a vector of length Q, the number of items for each FE
    // - obs2keep => vector of observations to keep. The neutral is 0.

    int Q = Rf_length(x);
    SEXP xq = VECTOR_ELT(x, 0);
    int n = Rf_length(xq);

    vector<bool> check_pblm(Q, true);
    if(only_slope.size() == Q){
        for(int q=0 ; q<Q ; ++q){
            // we check for pblm only if NOT only slope
            check_pblm[q] = only_slope[q] == false;
        }
    }

    // Rcout << "Q = " << Q << "\n";

    double *py = nullptr;
    if(TYPEOF(y) == REALSXP){
        py = REAL(y);
    } else {
        // => there will be no use of y, so nullptr is OK
        // but I must ensure that beforehand: do_sum_y = rm_0 = rm_1 = false
        if(do_sum_y || rm_0){
            stop("y should not be a list when its values are assessed.");
        }
    }


    // I create a fake vector to avoid conditional calls later on
    vector<int> x_sizes_fake(Q, 0);
    int *px_sizes = do_refactor ? INTEGER(r_x_sizes) : x_sizes_fake.data();

    int n_keep = obs2keep.size();
    bool identical_x = do_refactor && obs2keep[0] == 0;

    // the vectors of qufed
    writable::list res_x_quf_all(Q);
    vector<int*> p_x_quf_all(Q);
    for(int q=0 ; q<Q ; ++q){
        if(identical_x){
            SEXP x_val = VECTOR_ELT(x, q);
            res_x_quf_all[q] = x_val;
            p_x_quf_all[q]   = INTEGER(x_val);
        } else {
            res_x_quf_all[q] = PROTECT(Rf_allocVector(INTSXP, do_refactor ? n_keep : n));
            p_x_quf_all[q]   = INTEGER(res_x_quf_all[q]);
        }
    }

    vector< vector<int> > x_table_all(Q);
    vector< vector<double> > x_unik_all(Q);
    vector<int> any_pblm(Q, 0);
    vector< vector<bool> > id_pblm_all(Q);
    // The following may not be needed:
    vector< vector<double> > sum_y_all(Q);
    vector<bool> obs_removed;
    vector< vector<double> > x_removed_all(Q);

    // New code => we avoid passing SEXP within the parallel loop
    // We get the types of the R vectors
    // For strings => we convert into numeric without parallel (sigh)

    vector<void *> px_all(Q);
    vector<std::string> x_type_all(Q);
    // vector to store modified strings
    vector< vector<unsigned long long> > x_ull_all(Q);

    for(int q=0 ; q<Q ; ++q){

        xq = VECTOR_ELT(x, q);

        if(TYPEOF(xq) == INTSXP){
            x_type_all[q] = "int";
            px_all[q] = INTEGER(xq);

        } else if(TYPEOF(xq) == REALSXP){
            x_type_all[q] = "double";
            px_all[q] = REAL(xq);

        } else if(TYPEOF(xq) == STRSXP){
            // We make the conversion to unsigned long long
            x_type_all[q] = "string";

            std::uintptr_t xi_uintptr;

            for(int i=0 ; i<n ; ++i){
                const char *pxi = CHAR(STRING_ELT(xq, i));
                xi_uintptr = reinterpret_cast<std::uintptr_t>(pxi);

                x_ull_all[q].push_back(static_cast<unsigned long long>(xi_uintptr));

                // Rcout << xi_uintptr << "  ----  " << xi_ull << "\n";
            }

            px_all[q] = x_ull_all[q].data();
        } else {
            // We NEVER end here
            stop("Error: wrong type.");
        }
    }

    #pragma omp parallel for num_threads(nthreads)
    for(int q=0 ; q<Q ; ++q){
        quf_table_sum_single(px_all[q], x_type_all[q], n, q, p_x_quf_all[q],
                             x_unik_all[q], x_table_all[q],
                             py, sum_y_all[q], do_sum_y, rm_0, rm_1,
                             rm_single, any_pblm, id_pblm_all[q], check_pblm[q],
                             do_refactor, px_sizes[q], obs2keep);
    }


    bool is_pblm = false;
    for(int q=0 ; q<Q ; ++q){
        if(any_pblm[q]){
            is_pblm = true;
            break;
        }
    }

    // Rcout << "Any problem: " << is_pblm << "\n";


    if(obs2keep[0] != 0){
        n = n_keep;
    }

    if(is_pblm){

        //
        // finding which observation to remove based on the removed fixed-effects
        //

        // New scheme:
        // - if parallel, we create an int vector gathering the obs removed
        //   => this will allows to avoid a critical section which makes the
        //      parallel useless
        //   we then copy the values into the bool vector
        // - if not parallel, we fill the bool vector directly
        //

        int n_problems = std::accumulate(any_pblm.begin(), any_pblm.end(), 0);
        bool is_parallel = n_problems > 1 && nthreads > 1;

        // creating the obs2remove vector
        obs_removed.resize(n);
        std::fill(obs_removed.begin(), obs_removed.end(), false);

        if(is_parallel){

            // using an int vector removes the race condition issue
            // => we can't use a boolean vector directly bc of the way
            // boolean vectors are stored (allocation of new values depends
            // on the position of the existing values => creating a race condition
            // problem)
            //
            // using a omp critical section is a no go, renders parallel useless
            vector<int> obs_removed_int(n, 0);

            // TODO: OMP functions #pragma omp parallel for num_threads(nthreads)
            for(int q=0 ; q<Q ; ++q){
                if(any_pblm[q]){
                    vector<bool> &id_pblm = id_pblm_all[q];
                    int *pquf = p_x_quf_all[q];
                    for(int i=0 ; i<n ; ++i){
                        if(id_pblm[pquf[i] - 1]){
                            obs_removed_int[i] = 1;
                        }
                    }
                }
            }

            // we copy the values into the boolean vector
            for(int i=0 ; i<n ; ++i){
                if(obs_removed_int[i]){
                    obs_removed[i] = true;
                }
            }

        } else {

            // we fill the boolean vector directly
            for(int q=0 ; q<Q ; ++q){
                if(any_pblm[q]){
                    vector<bool> &id_pblm = id_pblm_all[q];
                    int *pquf = p_x_quf_all[q];
                    for(int i=0 ; i<n ; ++i){
                        if(id_pblm[pquf[i] - 1]){
                            obs_removed[i] = true;
                        }
                    }
                }
            }
        }


        // refactoring, recomputing all the stuff
        // ie x_quf, x_unik, x_table, sum_y (if needed)
        // but also: x_removed

        // The values of x_table and sum_y are changed in place
        // I create the addition new_quf and new_unik because we need the
        // them old items for the construction of the new.

        int n_removed = std::accumulate(obs_removed.begin(), obs_removed.end(), 0);
        int n_new = n - n_removed;

        writable::list res_x_new_quf_all(Q);
        vector<int*> p_x_new_quf_all(Q);
        for(int q=0 ; q<Q ; ++q){
            res_x_new_quf_all[q] = PROTECT(Rf_allocVector(INTSXP, n_new));
            p_x_new_quf_all[q] = INTEGER(res_x_new_quf_all[q]);
        }

        vector< vector<double> > x_new_unik_all(Q);
        bool stop_now = false;
        bool *pstop_now = &stop_now;

        // TODO: OMP functions #pragma omp parallel for num_threads(nthreads)
        for(int q=0 ; q<Q ; ++q){
            quf_refactor_table_sum_single(n, p_x_quf_all[q], p_x_new_quf_all[q], obs_removed,
                                          x_unik_all[q], x_new_unik_all[q], x_removed_all[q],
                                          x_table_all[q], py, sum_y_all[q], do_sum_y,
                                          rm_1, rm_single, id_pblm_all[q], check_pblm[q],
                                          pstop_now, Q);
        }

        if(*pstop_now){
            UNPROTECT((Q * (!identical_x)) + Q);
            stop("The dependent variable is fully explained by the fixed-effects.");
        }

        // Saving the values in the appropriate locations
        res_x_quf_all = res_x_new_quf_all;
        x_unik_all = x_new_unik_all;

    }

    // The object to be returned

    writable::list res;

    // x: UF (the largest object => no copy)
    res.push_back({"quf"_nm = res_x_quf_all});

    // x: Unik
    writable::list res_tmp(Q);
    res.push_back({"items"_nm = x_unik_all});

    // table
    res.push_back({"table"_nm = x_table_all});

    // sum y

    // Rcpp
    // for(int q=0 ; q<Q ; ++q){
    //     if(do_sum_y){
    //         res_tmp[q] = sum_y_all[q];
    //     } else {
    //         res_tmp[q] = 0;
    //     }
    // }
    // res["sum_y"] = clone(res_tmp);

    // cpp11
    writable::doubles copy_sum_y_all = sum_y_all;
    for (int q = 0; q < Q; ++q){
        if (do_sum_y) {
            copy_sum_y_all[q] = sum_y_all[q];
        } else {
            copy_sum_y_all[q] = 0.0;
        }
    }
    res.push_back({"sum_y"_nm = copy_sum_y_all});

    //
    // IF PROBLEM ONLY
    //

    if (is_pblm)
    {
        // The removed observations
        res.push_back({"obs_removed"_nm = obs_removed});

        // The removed clusters
        res.push_back({"fe_removed"_nm = copy_sum_y_all});
    }

    UNPROTECT((Q * (!identical_x)) + (Q * is_pblm));

    return res;
}
