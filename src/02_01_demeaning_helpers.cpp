#include "02_0_demeaning.h"

void stayIdleCheckingInterrupt(bool *stopnow, vector<int> &jobdone, int n_vars, int *counterInside)
{
    // function that keeps the master thread busy until everything is done

    int nbDone = 0, iter = 1;
    bool isMaster = omp_get_thread_num() == 0;

    while (isMaster && nbDone < n_vars && !(*stopnow))
    {
        ++iter;

        if (iter % 500000000 == 0)
        {
            // if(pending_interrupt()){
            //     (*counterInside)++;
            //     *stopnow = true;
            //     break;
            // } else {
            //     // to avoid int overflow:
            //     iter = 0;
            // }
            iter = 0;
        }

        if (iter % 1000000 == 0)
        {
            nbDone = 0;
            for (int v = 0; v < n_vars; v++)
            {
                nbDone += jobdone[v];
            }
        }
    }
}

std::vector<int> set_parallel_scheme_ter(int N, int nthreads)
{
    // => this concerns only the parallel application on a 1-Dimensional matrix
    // takes in the nber of observations of the vector and the nber of threads
    // gives back a vector of the length the nber of threads + 1 giving the start/stop of each threads

    std::vector<int> res(nthreads + 1, 0);
    double N_rest = N;

    for (int i = 0; i < nthreads; ++i)
    {
        res[i + 1] = ceil(N_rest / (nthreads - i));
        N_rest -= res[i + 1];
        res[i + 1] += res[i];
    }

    return res;
}

[[cpp11::register]] list cpp_which_na_inf_(SEXP x, int nthreads)
{
    // x: vector, matrix, data.frame // double or integer

    /*
     This function takes a matrix and looks at whether it contains NA or infinite values
     return: flag for na/inf + logical vector of obs that are Na/inf
     std::isnan, std::isinf are OK since cpp11 required
     do_any_na_inf: if high suspicion of NA present: we go directly constructing the vector is_na_inf
     in the "best" case (default expected), we need not construct is_na_inf
     */

    sMat mat(x);

    int nobs = mat.nrow();
    int K = mat.ncol();
    bool anyNAInf = false;
    bool any_na = false;  // return value
    bool any_inf = false; // return value

    /*
     we make parallel the anyNAInf loop
     why? because we want that when there's no NA (default) it works as fast as possible
     if there are NAs, single threaded mode is faster, but then we circumvent with the do_any_na_inf flag
     */

    // no need to care about the race condition
    // "trick" to make a break in a multi-threaded section

    std::vector<int> bounds = set_parallel_scheme_ter(nobs, nthreads);

    // TODO: OMP functions
    #pragma omp parallel for num_threads(nthreads)
    for (int t = 0; t < nthreads; ++t)
    {
        for (int k = 0; k < K; ++k)
        {
            for (int i = bounds[t]; i < bounds[t + 1] && !anyNAInf; ++i)
            {

                if (mat[k].is_int)
                {
                    if (mat(i, k) == -2147483648.0)
                    {
                        anyNAInf = true;
                    }
                }
                else if (std::isnan(mat(i, k)) || std::isinf(mat(i, k)))
                {
                    anyNAInf = true;
                }
            }
        }
    }

    // object to return: is_na_inf
    writable::logicals is_na_inf(anyNAInf ? nobs : 1);

    // fill is_na_inf with false
    if (anyNAInf)
    {
        std::fill(is_na_inf.begin(), is_na_inf.end(), false);
    }

    if (anyNAInf)
    {
    // TODO: OMP functions
    #pragma omp parallel for num_threads(nthreads)
        for (int i = 0; i < nobs; ++i)
        {
            double x_tmp = 0;
            for (int k = 0; k < K; ++k)
            {
                x_tmp = mat(i, k);
                if (mat[k].is_int)
                {
                    if (mat(i, k) == -2147483648.0)
                    {
                        is_na_inf[i] = true;
                        any_na = true;
                        break;
                    }
                }
                else if (std::isnan(x_tmp))
                {
                    is_na_inf[i] = true;
                    any_na = true;
                    break;
                }
                else if (std::isinf(x_tmp))
                {
                    is_na_inf[i] = true;
                    any_inf = true;
                    break;
                }
            }
        }
    }

    // Return
    writable::list res;
    res.push_back({"any_na"_nm = any_na});
    res.push_back({"any_inf"_nm = any_inf});
    res.push_back({"any_na_inf"_nm = any_na || any_inf});
    res.push_back({"is_na_inf"_nm = is_na_inf});

    return res;
}
